package context;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import augmentedTree.IntervalTree;
import tools.BufferedRandomAccessFile;
import tools.FileSorter7;
import tools.UnixSort;
import main.Context;
import main.IOCoordinator;
import main.InitialRead;
import main.InitialReadLocation;
import main.Microbe;
import main.Read;
import main.ReadLocation;
import main.UnsynchronizedBufferedWriter;

public class ContextProcessingCoordinator extends Thread {
	
	private MappingProcessor mappingProcessor;
	private IOCoordinator ioCoordinator;
	//private BufferedRandomAccessFile rmapReader;
	private File chrFile;
	private Context context;
	private String outputDirPath;
	
	private int maxMissmatchDifference;
	private int readLength;
	private int maxMissmatches;
	private int minStartPositionOverlaps;
	private int maxContextSize;
	private int maxGapSize;
	private boolean preferExtensionsWithKnownSpliceSignal;
	private boolean skipDenovoJunctions;
	private boolean skipNonCanonicalJunctions;
	
	private ArrayList<Integer> windowSizes = new ArrayList<Integer>();
	private IntervalTree<Microbe> intervalTree;
	private boolean updateQueue;
	private boolean pairedEnd;
	private int updateInterval;
	private boolean verbose;
	private boolean developer;
	private boolean isLargeContext;
	
	
	
	
	private ArrayList<ActionListener> listeners;
	
	public ContextProcessingCoordinator(MappingProcessor mappingProcessor,IOCoordinator sortCoordinator,Context context,File chrFile, String outputDirPath,IntervalTree<Microbe> intervalTree, ArrayList<Integer> windowSizes,int maxMissmatchDifference,int maxMissmatches, int minStartPositionOverlaps, int readLength,int maxContextSize, int maxGapSize, boolean preferExtensionsWithKnownSpliceSignal, boolean skipDenovoJunctions, boolean skipNonCanonicalJunctions, boolean updateQueue, boolean pairedEnd, int updateInterval, boolean verbose,boolean developer, boolean isLargeContext) {
		super();
		try {
			
			this.listeners = new ArrayList<ActionListener>();
			this.mappingProcessor = mappingProcessor;
			this.ioCoordinator = sortCoordinator;
			this.chrFile = chrFile;
			this.context = context;
			this.outputDirPath = outputDirPath;
			this.maxMissmatchDifference = maxMissmatchDifference;
			this.readLength = readLength;
			this.maxMissmatches = maxMissmatches;
			this.minStartPositionOverlaps = minStartPositionOverlaps;
			this.maxContextSize = maxContextSize;
			this.maxGapSize = maxGapSize;
			this.windowSizes = windowSizes;
			this.preferExtensionsWithKnownSpliceSignal = preferExtensionsWithKnownSpliceSignal;
			this.skipDenovoJunctions = skipDenovoJunctions;
			this.skipNonCanonicalJunctions = skipNonCanonicalJunctions;
			this.updateQueue = updateQueue;
			this.pairedEnd = pairedEnd;
			this.updateInterval = updateInterval;
			this.verbose = verbose;
			this.developer = developer;
			this.isLargeContext = isLargeContext;
			this.intervalTree = intervalTree;
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public void run() {
		try {
			
			File outputDir = new File(this.outputDirPath);
			if(!outputDir.isDirectory())
				outputDir.mkdirs();
			
			//build local reference
			StringBuilder localReference = new StringBuilder();
			int localStart = Math.max(1,this.context.getStart() - this.maxGapSize);
			int localEnd = Math.min(this.mappingProcessor.getCurrentChromosome().length(), this.context.getEnd() + this.maxGapSize);
			localReference.append(this.mappingProcessor.getCurrentChromosome().substring(localStart - 1,localEnd).toUpperCase());
			fireAction(new ActionEvent(this,0,"chromosome_unlocked"));
			
			
			
			BufferedRandomAccessFile rmapReader = new BufferedRandomAccessFile(this.chrFile, "r",10240);
			
			if(this.verbose) System.out.println(String.format("%s\tprocessing context:\t%s (containing %s candidates)", this.toString(),this.context.getId(),this.context.getContainedReads()));
			long prevTimePoint = System.currentTimeMillis();
			long currentTimePoint;
			double usedTime;
			double overallUsedTime = 0.0;
			FileSorter7 fileSorter;
			String alignmentsFilePath = null;

			
			//if context is too large, we sort the alignments of the current context by read id and process the reads in seperate blocks
			if(this.isLargeContext) {
				writeAlignments(rmapReader,this.outputDirPath + "/alignments.rmap");
				fileSorter = new FileSorter7(this.outputDirPath + "/alignments.rmap", this.outputDirPath + "/alignments.rmap.sorted",null,new int[]{0},250, "\t",false);
				fileSorter.sortFile(this.outputDirPath + "/alignments.rmap", this.outputDirPath + "/alignments.rmap.sorted",null);
				fileSorter = null;
				alignmentsFilePath = this.outputDirPath + "/alignments.rmap.sorted";
			}
			
			//get mapping
			if(this.verbose) System.out.print(this.toString() + "\tparsing mapping...");
			this.mappingProcessor.parseMapping(this.context,rmapReader, this.pairedEnd,localReference, localStart,alignmentsFilePath);
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			//extend mapping
			String multiMappingOutputPath = this.outputDirPath + "/multi_mapping.txt";
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tdetermining split alignments...");
			this.mappingProcessor.extendMapping(this.context,localReference,localStart,multiMappingOutputPath,this.outputDirPath + "/multi_mapping_with_multi_splits.txt",this.minStartPositionOverlaps, alignmentsFilePath);
			this.context.resetReads();
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			
			//extracting best matching mappings
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\textracting best matching mappings...");
			new MappingProcessor(multiMappingOutputPath, outputDirPath + "/multi_mapping.best.matchings.txt", maxMissmatchDifference,preferExtensionsWithKnownSpliceSignal,skipDenovoJunctions,skipNonCanonicalJunctions).extractBestMatchingAlignmentsInLocalResolution(context,true);
			context.clearPartialAndFullReads2filePointer();
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			
			//sorting best matching mappings by start positions
			if(this.verbose) System.out.print(this.toString() + "\tsorting best mappings by start positions...");
			prevTimePoint = System.currentTimeMillis();
			

			while(!this.ioCoordinator.hasFreeSlot()) {
				Thread.sleep(500);
			}
			this.ioCoordinator.inputOutputOperationStarted();
			fileSorter = new FileSorter7(this.outputDirPath + "/multi_mapping.best.matchings.txt", this.outputDirPath + "/multi_mapping.best.matchings.txt.startposition.sorted",null,new int[]{4},(!this.isLargeContext)?45:250, "\t",true);
			fileSorter.sortFile(this.outputDirPath + "/multi_mapping.best.matchings.txt", this.outputDirPath + "/multi_mapping.best.matchings.txt.startposition.sorted",null);
			fileSorter = null;
			this.ioCoordinator.inputOutputOperationEnded();
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			
			//resolving pairwise overlapping splice sites
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tresolving overlapping splice sites...");
			this.mappingProcessor.resolveOverlappingSpliceSites(this.outputDirPath + "/multi_mapping.best.matchings.txt.startposition.sorted",this.readLength, this.maxMissmatches,localReference,localStart, this.outputDirPath + "/multi_mapping.resolved.splicesites.txt");
			currentTimePoint = System.currentTimeMillis();
			usedTime = (currentTimePoint - prevTimePoint)/1000;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			
			
			/**
			 * NEW: Multi-Split Detection Block
			 */
			
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\textending full and split candidates...");
			this.mappingProcessor.extendSplitAndFullAlignments(this.outputDirPath + "/multi_mapping.resolved.splicesites.txt", this.context,localReference,localStart, this.outputDirPath + "/multi_mapping_with_multi_splits.txt", this.minStartPositionOverlaps);
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			while(!this.ioCoordinator.hasFreeSlot()) {
				Thread.sleep(500);
			}
			
			this.ioCoordinator.inputOutputOperationStarted();
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tsorting candidates by read id...");
			fileSorter = new FileSorter7(this.outputDirPath + "/multi_mapping.resolved.splicesites.txt", this.outputDirPath + "/multi_mapping.resolved.splicesites.txt.sorted",null,new int[]{1},(!this.isLargeContext)?45:250, "\t",false);
			fileSorter.sortFile(this.outputDirPath + "/multi_mapping.resolved.splicesites.txt", this.outputDirPath + "/multi_mapping.resolved.splicesites.txt.sorted",null);
			fileSorter = null;
			this.ioCoordinator.inputOutputOperationEnded();
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tdetermining multi split candidates...");
			this.mappingProcessor.buildAllMultiSplitCombinations(this.outputDirPath + "/multi_mapping.resolved.splicesites.txt.sorted", this.outputDirPath + "/multi_mapping_with_multi_splits.txt", pairedEnd,localReference,localStart);
			//unlock read sequences
			fireAction(new ActionEvent(this,0,"reads_unlocked"));
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tsorting candidates by read id...");
			fileSorter = new FileSorter7(this.outputDirPath + "/multi_mapping_with_multi_splits.txt", this.outputDirPath + "/multi_mapping_with_multi_splits.txt.sorted",null,new int[]{1},(!this.isLargeContext)?45:250, "\t",false);
			fileSorter.sortFile(this.outputDirPath + "/multi_mapping_with_multi_splits.txt", this.outputDirPath + "/multi_mapping_with_multi_splits.txt.sorted",null);
			fileSorter = null;
			this.ioCoordinator.inputOutputOperationEnded();
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\textracting best matching alignments...");
			new MappingProcessor(this.outputDirPath + "/multi_mapping_with_multi_splits.txt.sorted",this.outputDirPath + "/multi_mapping_with_multi_splits.txt.bestmatches",this.maxMissmatchDifference,this.preferExtensionsWithKnownSpliceSignal,this.skipDenovoJunctions,this.skipNonCanonicalJunctions).extractBestMatchingAlignments(false);
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			
			//resolving context
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tresolving context...");
			ContextResolver contextResolver;
			if(!this.pairedEnd)
				contextResolver = new LocalContextResolverSingleEnd(this.outputDirPath + "/multi_mapping_with_multi_splits.txt.bestmatches", String.format("%s/resolved_%s.txt",this.outputDirPath,context.getId()),this.context.getUpstreamCoverage(),this.context.getDownstreamCoverage(),this.intervalTree, windowSizes,readLength,this.maxContextSize,updateQueue,updateInterval,verbose,developer,isLargeContext);
			else
				contextResolver = new LocalContextResolverPairedEnd(this.outputDirPath + "/multi_mapping_with_multi_splits.txt.bestmatches", String.format("%s/resolved_%s.txt",this.outputDirPath,context.getId()),this.context.getUpstreamCoverage(),this.context.getDownstreamCoverage(),this.intervalTree, windowSizes,readLength,this.maxContextSize,updateQueue,updateInterval,verbose,developer,isLargeContext);
			
			contextResolver.resolve();
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			
			
			//sorting result
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tsorting resolved context...");
			new UnixSort().sort(String.format("%s/resolved_%s.txt",this.outputDirPath,context.getId()),String.format("%s/resolved_%s.txt.sorted",this.outputDirPath,context.getId()),this.outputDirPath + "/tmp","\t",2,100,false,false,this.verbose);
			new File(String.format("%s/resolved_%s.txt",this.outputDirPath,context.getId())).delete();
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			
			
			//move result to tmp folder
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tmoving result to output folder...");
			File fileToMove = new File(String.format("%s/resolved_%s.txt.sorted",this.outputDirPath,context.getId()));
			String newPath = this.outputDirPath.substring(0,this.outputDirPath.lastIndexOf(System.getProperty("file.separator")) + 1) + "resolved_local_contexts/" + fileToMove.getName();
			fileToMove.renameTo(new File(newPath));
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));			
			
			//now clear the content of the current outputDir
			prevTimePoint = System.currentTimeMillis();
			if(this.verbose) System.out.print(this.toString() + "\tdeleting temporary files...");
			File[] fileContent = new File(this.outputDirPath).listFiles();
			for(File f : fileContent) {
				f.delete();
			}
			new File(this.outputDirPath).delete();
		
		
			rmapReader.close();
			this.listeners.clear();
			this.context.clearUpAndDownstreamCoverages();
			
			currentTimePoint = System.currentTimeMillis();
			usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
			overallUsedTime += usedTime;
			if(this.verbose) System.out.println(String.format("done (%s sec).",usedTime));
			if(this.verbose) System.out.println(String.format(this.toString() + "\tfinished. took %s seconds.",overallUsedTime));
			if(this.verbose) System.out.println();
			
		}
		catch(Exception e) {
			e.printStackTrace();
			System.err.println("Error occured in context: " + this.context.getId());
			fireAction(new ActionEvent(this,0,"unlocked"));
		}
	}
	
	
	private void writeAlignments(BufferedRandomAccessFile rmapReader, String outputFilePath) throws Exception {
		PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
		String currentLine;
		rmapReader.seek(this.context.getPointerToFirstRead());
		while((currentLine = rmapReader.getNextLine()) != null) {
			pw.println(currentLine);
			
			if(rmapReader.getFilePointer() == this.context.getPointerToLastRead())
				break;
		}
		pw.close();
	}
	
	
	public void addListener(ActionListener listener) {
		this.listeners.add(listener);
	}
	
	public void fireAction(ActionEvent e) {
		for(ActionListener listener : this.listeners) {
			listener.actionPerformed(e);
		}
	}
}
