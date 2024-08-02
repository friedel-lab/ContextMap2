package main;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import org.mapdb.BTreeKeySerializer;
import org.mapdb.BTreeMap;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Fun;
import org.mapdb.HTreeMap;
import org.mapdb.Pump;
import org.mapdb.Serializer;
import org.mapdb.DB.BTreeMapMaker;

import alignment.AlignmentCoordinator;
import augmentedTree.IntervalTree;
import context.ContextExtractor;
import context.ContextProcessingCoordinator;
import context.ContextResolver;
import context.GlobalContextResolverPairedEnd;
import context.GlobalContextResolverSingleEnd;
import context.MappingProcessor;
import tools.BufferedRandomAccessFile;
import tools.FileSorter7;
import tools.ReadAligner;
import tools.RmapProcessor;
import tools.SamProcessor;
import tools.String2Bitset;
import tools.UnixSort;


public class ContextMap implements ActionListener {

	/*
	 * global stuff
	 */
	private String contextmapVersion;
	private String readFilePath;
	private String samFilePath;
	private String outputDirPath;
	private String tmpOutputDirPath;
	private String readFormat;
	private int readLength;
	private final int maxIOThreads = Integer.MAX_VALUE;
	
	/*
	 * alignment stuff
	 */
	private String alignerName;
	private String alignerBinPath;
	private String alignerIndexerBinPath;
	private String alignerTmpDirPath;
	private ArrayList<String> genomeIndexBasePaths;
	private AlignmentCoordinator slidingContext;
	
	
	//optional fields (default settings are defined in the constructor)
	private int seedLength;
	private int seedMissmatches;
	private int maxMismatches;
	private int[] splitSeedSizes;
	private int[] splitSeedMismatches;
	private int maxIntronLength;
	private boolean autosetAlignmentOptions;

	private int maxHits;
	private String skippedReadsPath;
	private String unalignedReadsPath;
	private int threads;
	private boolean skipPass1;
	private ArrayList<Boolean> skipSplitDetection;
	private ArrayList<Boolean> skipMultiSplitDetection;
	private boolean skipRealignment;
	
	/*
	 * context processing stuff (default settings are defined in the constructor)
	 */
	private int maxMismatchDifference;
	private String referencesDir;
	private String annotationFilePath;
	private String gtfFilePath;
	private String indexDirPath;
	private int minDistanceBetweenContext;
	private int maxContextSize;
	
	private int maxGapSize;
	private int minGapSize;
	private int maxDelSize;
	
	//polyA prediction parameters
	private int minPolyALength;
	private int minPolyAReadCount;
	private int maxConsideredClippingLength;
	private double upperPolyACutoff;
	private double lowerPolyACutoff;
	
	
	private int minNumberOfReadsInContext;
	//at least #minStartPositionOverlaps have to overlap with a full/partial alignment before it will be made any effort to extend it to a split candidate. 
	private int minStartPositionOverlaps;
	private ArrayList<Integer> windowSizes;
	private int completeWindowSize;
	private boolean strandSpecific;
	private boolean pairedEnd;
	private boolean preferExtensionsWithKnownSpliceSignals;
	private boolean skipDenovoJunctions;
	private boolean skipNonCanonicalJunctions;
	private boolean updateQueue;
	private boolean clipping;
	private boolean polyA;
	private boolean strandedPolyA;
	private boolean printMultiMappings;
	private int updateInterval;
	
	private boolean printSecondBestChr;
	private boolean verbose;
	private boolean developer;
	private boolean keepTmp;
	
	private boolean writeSequenceDB;
	private int databaseSortBatchSize;
	private final int read2sequenceSize = 20000000;
	
	
	private AtomicInteger contextsNeedingBufferedChromosomeSequence;
	private AtomicInteger contextsNeedingBufferedReadSequences;
	
	private AtomicInteger threadsBufferingReadSequences;
	private AtomicInteger threadsBufferingChromosomeSequences;
	
	public ContextMap(String contextmapVersion, String readFilePath, String samFilePath, String alignerName, String alignerBinPath, String alignerTmpDirPath, String alignerIndexerBinPath, ArrayList<String> genomeIndexBasePaths, String referencesDir, String annotationFilePath, String gtfFilePath, String indexDirPath, String outputDirPath, String readFormat, int readLength,
			int seedLength, int seedMissmatches, int maxMismatches, int maxMismatchDifference, int[] splitSeedSizes, int[] splitSeedMismatches, int maxIntronLength, int maxIntronCount, int maxHits, int threads,
			boolean skipPass1,ArrayList<Boolean> skipSplitDetection, ArrayList<Boolean> skipMultiSplitDetection, boolean skipRealignment, int minDistanceBetweenContext, int maxContextSize, int maxGapSize, int minGapSize, int maxDelSize, int minNumberOfReadsInContext,boolean strandSpecific, boolean pairedEnd, boolean preferExtensionsWithKnownSpliceSignals, boolean skipDenovoJunctions, boolean skipNonCanonicalJunctions, boolean updateQueue, int updateInterval, boolean printMultiMappings, boolean printSecondBestChr, boolean lowCoverage, boolean clipping, boolean polyA, boolean strandedPolyA, int minPolyALength, int minPolyAReadCount, double upperPolyACutoff, double lowerPolyACutoff, int maxConsideredClippingLength, boolean verbose, boolean developer, boolean keepTmp, boolean writeSequenceDB, boolean autosetAlignmentOptions, int databaseSortBatchSize) {
		
		//init the global settings
		this.contextmapVersion = contextmapVersion;
		this.readFilePath = readFilePath;
		this.samFilePath = samFilePath;
		this.outputDirPath = outputDirPath;
		this.tmpOutputDirPath = outputDirPath + "/tmp";
		this.alignerName = alignerName;
		this.alignerBinPath = alignerBinPath;
		this.alignerIndexerBinPath = alignerIndexerBinPath;
		this.genomeIndexBasePaths = genomeIndexBasePaths;
		this.referencesDir = referencesDir;
		this.annotationFilePath = annotationFilePath;
		this.gtfFilePath = gtfFilePath;
		this.indexDirPath = indexDirPath;
		this.readFormat = readFormat;
		this.readLength = readLength;
		
		//init default alignment settings
		if(alignerTmpDirPath == null)
			this.alignerTmpDirPath = this.tmpOutputDirPath + "/aligner_tmp";
		else
			this.alignerTmpDirPath = alignerTmpDirPath;
		
		this.seedLength = seedLength;
		this.seedMissmatches = seedMissmatches;
		this.maxMismatches = maxMismatches;
		this.maxMismatchDifference = maxMismatchDifference;
		this.splitSeedSizes = splitSeedSizes;
		this.splitSeedMismatches = splitSeedMismatches;
		this.maxIntronLength = maxIntronLength;
		this.maxHits = maxHits;
		this.threads = threads;
		this.autosetAlignmentOptions = autosetAlignmentOptions;
		
		this.skipPass1 = skipPass1;
		this.skipSplitDetection = skipSplitDetection;
		this.skipMultiSplitDetection = skipMultiSplitDetection;
		this.skipRealignment = skipRealignment;
	
		
		//init context processing stuff
		this.minDistanceBetweenContext = minDistanceBetweenContext;
		this.maxContextSize = maxContextSize;
		
		this.maxGapSize = maxGapSize;
		this.minGapSize = minGapSize;
		this.maxDelSize = maxDelSize;
		
		
		this.minPolyALength= minPolyALength;
		this.minPolyAReadCount = minPolyAReadCount;
		this.upperPolyACutoff = upperPolyACutoff;
		this.lowerPolyACutoff = lowerPolyACutoff;
		this.maxConsideredClippingLength = maxConsideredClippingLength;
		
		this.minNumberOfReadsInContext = minNumberOfReadsInContext;
		this.strandSpecific = strandSpecific;
		this.pairedEnd = pairedEnd;
		this.preferExtensionsWithKnownSpliceSignals = preferExtensionsWithKnownSpliceSignals;
		this.skipDenovoJunctions = skipDenovoJunctions;
		this.skipNonCanonicalJunctions = skipNonCanonicalJunctions;
		this.clipping = clipping;
		this.polyA = polyA;
		this.strandedPolyA = strandedPolyA;
		this.windowSizes = new ArrayList<Integer>();
		this.windowSizes.add(25);
		this.windowSizes.add(75);
		this.completeWindowSize = 100;
		this.updateQueue = updateQueue;
		this.updateInterval = updateInterval;
		this.printMultiMappings = printMultiMappings;
		this.printSecondBestChr = printSecondBestChr;
		if(lowCoverage)
			this.minStartPositionOverlaps = 1;
		else
			this.minStartPositionOverlaps = 4;
		
		this.contextsNeedingBufferedChromosomeSequence = new AtomicInteger(0);
		this.contextsNeedingBufferedReadSequences = new AtomicInteger(0);
		this.threadsBufferingReadSequences = new AtomicInteger(0);
		this.threadsBufferingChromosomeSequences = new AtomicInteger(0);
		this.verbose = verbose;
		this.developer = developer;
		this.keepTmp = keepTmp;
		this.writeSequenceDB = writeSequenceDB;
		this.databaseSortBatchSize = databaseSortBatchSize;
		
		
		
		this.slidingContext = new AlignmentCoordinator(this.alignerName, this.alignerBinPath, this.alignerIndexerBinPath, this.alignerTmpDirPath, this.referencesDir,
				   this.maxMismatches, this.seedLength, this.seedMissmatches, this.maxHits, this.splitSeedSizes, this.splitSeedMismatches, this.maxContextSize, this.maxIntronLength, this.maxGapSize, this.minGapSize, this.maxDelSize, this.threads,
				   this.pairedEnd, this.skipPass1, this.autosetAlignmentOptions, this.verbose);
		
		
		
		this.seedLength =this.slidingContext.getSeedSize();
		this.seedMissmatches = this.slidingContext.getSeedMismatches();
		this.maxHits = this.slidingContext.getMaxHits();
		this.splitSeedSizes = this.slidingContext.getSplitSeedSizes();
		this.splitSeedMismatches = this.slidingContext.getSplitSeedMismatches();
		
		
	}
	
	public void start() {
		try {
			Date date;
			//first check if output destiniations exists
			File tmpFile = new File(this.tmpOutputDirPath);
			if(!tmpFile.exists()) tmpFile.mkdirs();
			tmpFile = new File(this.alignerTmpDirPath);
			if(!tmpFile.exists()) tmpFile.mkdirs();
			tmpFile = new File(this.tmpOutputDirPath + "/splitted_rmap");
			if(!tmpFile.exists()) tmpFile.mkdir();
			tmpFile = new File(this.tmpOutputDirPath + "/splitted_rmap_sorted");
			if(!tmpFile.exists()) tmpFile.mkdir();
			tmpFile = new File(this.tmpOutputDirPath + "/resolved_local_contexts");
			if(!tmpFile.exists()) tmpFile.mkdir();
			
			ReadAligner readAligner = new ReadAligner();
			RmapProcessor rmapProcessor = new RmapProcessor();
			String genomeIndexBasePath;
			String tmpAlignmentOutputPath;
			String tmpMultiSplitSeqsOutputPath;
			ArrayList<String> rmaps = new ArrayList<String>();
			ArrayList<String> multiSplitCandidates = new ArrayList<String>();
			File[] fileContent;
			int numberOfThreads = Integer.valueOf(this.threads);
			int maxSortThreads;
			//ContextMap on SAM version
			if(this.samFilePath != null) {
				SamProcessor samFileProcessor = new SamProcessor();
				if(!this.skipRealignment) {
						date = new Date();
						System.out.println(String.format("[%s]\tExtracting unmapped and multi mapped reads from sam file.",date.toLocaleString()));
						//first extract unmapped and multi mapped reads from sam file
						samFileProcessor.extractUnmappedAndMultimappedReads(this.samFilePath, this.readFilePath, this.readFormat, this.readLength, this.tmpOutputDirPath);
						
						samFileProcessor.convertSamToRmap(this.tmpOutputDirPath + "/uniquely_mapped_reads.sam", this.tmpOutputDirPath + "/uniquely_mapped_reads.rmap");
						rmaps.add(this.tmpOutputDirPath + "/uniquely_mapped_reads.rmap");
						//now align those reads
						
						for(int i = 0; i < this.genomeIndexBasePaths.size(); i++) {
							genomeIndexBasePath = this.genomeIndexBasePaths.get(i);
							tmpAlignmentOutputPath = this.tmpOutputDirPath + String.format("/unmapped_reads_%s.rmap",i);
							tmpMultiSplitSeqsOutputPath = this.tmpOutputDirPath + String.format("/multi_split_candidates_%s.seq",i);
							slidingContext.start(this.tmpOutputDirPath + "/unmapped_reads.fa", tmpAlignmentOutputPath,tmpMultiSplitSeqsOutputPath, genomeIndexBasePath, this.skipSplitDetection.get(i),this.skipMultiSplitDetection.get(i));
							
							rmaps.add(tmpAlignmentOutputPath);
							multiSplitCandidates.add(tmpMultiSplitSeqsOutputPath);
						}
						
						
						//concatenate uniquely mapped reads from sam file with the newly aligned reads
						date = new Date();
						System.out.println(String.format("[%s]\tConcatenating uniquely mappeds reads from input mapping with realigned reads.",date.toLocaleString()));
						rmapProcessor.concatenateFilesWithNIO(rmaps, this.tmpOutputDirPath + "/all_reads.rmap");
						rmapProcessor.concatenateFilesWithNIO(multiSplitCandidates, this.tmpOutputDirPath + "/multi_split_candidates.seq");
						new File(this.alignerTmpDirPath).delete();
				}
				
				else {
					samFileProcessor.convertSamToRmap(this.samFilePath, this.tmpOutputDirPath + "/all_reads.rmap");
				}
			}
			
			//ContextMap standalone version
			else {
				date = new Date();
				
				for(int i = 0; i < this.genomeIndexBasePaths.size(); i++) {
					genomeIndexBasePath = this.genomeIndexBasePaths.get(i);
					tmpAlignmentOutputPath = this.tmpOutputDirPath + String.format("/unmapped_reads_%s.rmap",i);
					tmpMultiSplitSeqsOutputPath = this.tmpOutputDirPath + String.format("/multi_split_candidates_%s.seq",i);
					slidingContext.start(this.readFilePath, tmpAlignmentOutputPath,tmpMultiSplitSeqsOutputPath, genomeIndexBasePath, this.skipSplitDetection.get(i),this.skipMultiSplitDetection.get(i));
					
					if(this.genomeIndexBasePaths.size() == 1) {
						new File(tmpAlignmentOutputPath).renameTo(new File(this.tmpOutputDirPath + "/all_reads.rmap"));
						new File(tmpMultiSplitSeqsOutputPath).renameTo(new File(this.tmpOutputDirPath + "/multi_split_candidates.seq"));
					}
					else {	
						rmaps.add(tmpAlignmentOutputPath);
						multiSplitCandidates.add(tmpMultiSplitSeqsOutputPath);
					}
				}
				new File(this.alignerTmpDirPath).delete();
				if(rmaps.size() > 0) {
					rmapProcessor.concatenateFilesWithNIO(rmaps, this.tmpOutputDirPath + "/all_reads.rmap");
					rmapProcessor.concatenateFilesWithNIO(multiSplitCandidates, this.tmpOutputDirPath + "/multi_split_candidates.seq");
				}
				
			}
			
			
			for(int i = 0; i < rmaps.size(); i++) {
				if(new File(rmaps.get(i)).exists())
					new File(rmaps.get(i)).delete();
				if(new File(multiSplitCandidates.get(i)).exists())
					new File(multiSplitCandidates.get(i)).delete();
			}
			rmaps.clear();
			multiSplitCandidates.clear();
			
			
			
			
			 // independent of running in standalone or in sam mode, all important stuff is now in tmpOutputPath + "/all_reads.rmap"
			// In case we did not find any alignment we can abort the run here.
			 
			BufferedReader checkReader = new BufferedReader(new FileReader(new File(this.tmpOutputDirPath + "/all_reads.rmap")));
			if(checkReader.readLine() == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tNo alignments found.",date.toLocaleString()));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				checkReader.close();
				System.exit(1);
			}
			checkReader.close();
			
			//sorting rmap file by read ids in order to pre filter split candidates.
			File[] splittedFiles;
			int threadIndex;
			File[] filesToDelete;
			maxSortThreads = numberOfThreads;
			
			ExecutorService sortExecutor = Executors.newFixedThreadPool(maxSortThreads);
			ArrayList<Future> futures = new ArrayList<Future>();
		
			date = new Date();
			System.out.println(String.format("[%s]\tSplitting and sorting candidate mappings by read id.",date.toLocaleString()));
			rmapProcessor.splitByChromosome(this.tmpOutputDirPath + "/all_reads.rmap",this.tmpOutputDirPath + "/splitted_rmap",numberOfThreads);
		
			splittedFiles = new File(this.tmpOutputDirPath + "/splitted_rmap").listFiles();
			threadIndex = 0; 
			for(File rmapFile : splittedFiles) {
				//futures.add(sortExecutor.submit(new FileSorter7(rmapFile.getAbsolutePath(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + rmapFile.getName(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + threadIndex,new int[]{0},100,"\t",false)));
				futures.add(sortExecutor.submit(new UnixSort(rmapFile.getAbsolutePath(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + rmapFile.getName(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + threadIndex,"\t",1,500,false,true,this.verbose)));
				threadIndex++;
			}
			sortExecutor.shutdown();
			//check if there is still a thread running
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			//delete unsorted rmap files
			filesToDelete = new File(this.tmpOutputDirPath + "/splitted_rmap").listFiles();
			for(File toDel : filesToDelete) {
				toDel.delete();
			}
			
			//delete tmp directories
			filesToDelete = new File(this.tmpOutputDirPath + "/splitted_rmap_sorted/").listFiles();
			for(File folderToDel : filesToDelete) {
				if(folderToDel.isDirectory())
					deleteFolderWithContent(folderToDel);
			}
			
			//now merge sorted files
			new FileSorter7(new File(this.tmpOutputDirPath + "/splitted_rmap_sorted"),new int[]{0},100,"\t",false).mergeFiles(new File(this.tmpOutputDirPath + "/splitted_rmap_sorted"), this.tmpOutputDirPath + "/all_reads.rmap.sorted",100);
			filesToDelete = new File(this.tmpOutputDirPath + "/splitted_rmap_sorted").listFiles();
			for(File toDel : filesToDelete) {
				toDel.delete();
			}
			
			
			
			
			
			//pre filter split candidates with too many mismatches
			String tmpAlignmentPath = this.tmpOutputDirPath + "/all_reads.rmap.sorted";
			
			date = new Date();
			System.out.println(String.format("[%s]\tRemoving split candidates with mismatch differences larger than %s (mmdiff) mismatches to the best alignment.",date.toLocaleString(),this.maxMismatchDifference));
			//new MappingProcessor().filterSplitCandidates(this.tmpOutputDirPath + "/all_reads.rmap.sorted", this.tmpOutputDirPath + "/all_reads.rmap.filtered",this.maxMismatchDifference,this.readLength);
			String filteredCandidatesDir = this.tmpOutputDirPath + "/filtered_rmaps";
			new File(filteredCandidatesDir).mkdir();
			new MappingProcessor().filterSplitCandidates(this.tmpOutputDirPath + "/all_reads.rmap.sorted", filteredCandidatesDir,this.tmpOutputDirPath + "/all_reads.rmap.filtered",this.maxMismatchDifference, this.readLength, this.threads,this.pairedEnd);
			
			new File(this.tmpOutputDirPath + "/all_reads.rmap.sorted").delete();
			tmpAlignmentPath = this.tmpOutputDirPath + "/all_reads.rmap.filtered";
			
			//in paired end mode remove read alignments which do not fit to a valid pair (only for mates with at least one valid pair)
			if(this.pairedEnd) {
				date = new Date();
				System.out.println(String.format("[%s]\tFiltering alignments of mate pairs.",date.toLocaleString()));
				new MappingProcessor().filterPairedEndCandidates(tmpAlignmentPath, filteredCandidatesDir, this.tmpOutputDirPath + "/all_reads.rmap.filtered.pairedend", this.readLength, this.maxContextSize, this.threads);
				//new MappingProcessor().filterPairedEndCandidates(tmpAlignmentPath, this.tmpOutputDirPath + "/all_reads.rmap.filtered.pairedend",this.readLength,this.maxContextSize);
				
				new File(this.tmpOutputDirPath + "/all_reads.rmap.filtered").delete();
			}
			


			//sorting rmap file by start positions
			date = new Date();
			System.out.println(String.format("[%s]\tSorting (filtered) candidate mappings by start positions.",date.toLocaleString()));
			//sorting splitted rmap files by start positions
			String allReadsPath = this.tmpOutputDirPath + "/all_reads.rmap.filtered";
			
			if(this.pairedEnd)
				allReadsPath = this.tmpOutputDirPath + "/all_reads.rmap.filtered.pairedend";
			
			rmapProcessor.splitByChromosome(allReadsPath,this.tmpOutputDirPath + "/splitted_rmap",numberOfThreads);			
		
			splittedFiles = new File(this.tmpOutputDirPath + "/splitted_rmap").listFiles();
			sortExecutor = Executors.newFixedThreadPool(maxSortThreads);
			futures.clear();
			threadIndex = 0; 
			for(File rmapFile : splittedFiles) {
				//futures.add(sortExecutor.submit(new FileSorter7(rmapFile.getAbsolutePath(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + rmapFile.getName(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + threadIndex,new int[]{3},100,"\t",true)));
				futures.add(sortExecutor.submit(new UnixSort(rmapFile.getAbsolutePath(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + rmapFile.getName(),this.tmpOutputDirPath + "/splitted_rmap_sorted/" + threadIndex,"\t",4,500,true,true,this.verbose)));
				threadIndex++;
			}
			sortExecutor.shutdown();
			//check if there is still a thread running
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			
			deleteFolderWithContent(new File(this.tmpOutputDirPath + "/splitted_rmap"));
			//delete tmp directories
			filesToDelete = new File(this.tmpOutputDirPath + "/splitted_rmap_sorted/").listFiles();
			for(File folderToDel : filesToDelete) {
				if(folderToDel.isDirectory())
					deleteFolderWithContent(folderToDel);
			}
			
			
			
			
			
			
			HashMap<String,File> chr2rmap = new HashMap<String,File>();
			String tmpChr;
			for(File rmapFile : new File(this.tmpOutputDirPath + "/splitted_rmap_sorted").listFiles()) {
				tmpChr = rmapFile.getName().substring(0,rmapFile.getName().lastIndexOf("."));
				chr2rmap.put(tmpChr, rmapFile);
				
			}

			
			//user is able to store read_id -> sequence mapping onto disk (recommended if ContextMap needs too much memory for the underlying data)
			ConcurrentMap<String,long[]> read2sequence = null;
			DB db = null;
			if(this.writeSequenceDB) {
				date = new Date();
				System.out.println(String.format("[%s]\tWriting readId -> sequence mapping to disk.",date.toLocaleString()));
				
				
				//init DB
				db = DBMaker.newFileDB(new File(this.tmpOutputDirPath + "/sequences.mapdb")).transactionDisable().closeOnJvmShutdown().mmapFileEnableIfSupported().asyncWriteEnable().cacheSize(500000).make();
				new File(this.tmpOutputDirPath + "/db_tmp").mkdir();
				System.setProperty("java.io.tmpdir", this.tmpOutputDirPath + "/db_tmp");
				

				Iterator<Fun.Tuple2<String,long[]>> readIdIt = new KeyIterator(this.readFilePath,false);
				Iterator<Fun.Tuple2<String,long[]>> mscIdIt = new KeyIterator(this.tmpOutputDirPath + "/multi_split_candidates.seq",true);
				Iterator<Fun.Tuple2<String,long[]>> mergedIt = Pump.merge(readIdIt,mscIdIt);
				
				mergedIt = Pump.sort(mergedIt,true, this.databaseSortBatchSize,Collections.reverseOrder(new Fun.Tuple2Comparator<String, long[]>(Fun.COMPARATOR,Fun.LONG_ARRAY_COMPARATOR)), db.getDefaultSerializer());				
				BTreeMapMaker mapMaker = db.createTreeMap("read2sequence").valueSerializer(Serializer.LONG_ARRAY).keySerializer(BTreeKeySerializer.STRING).nodeSize(12);
				mapMaker = mapMaker.pumpSource(mergedIt);
				read2sequence = mapMaker.make();
				
				db.commit();				
			}
			
			
			
			//splitting mscs by chr
			if(new File(this.tmpOutputDirPath + "/multi_split_candidates.seq").exists()) {
				File splittedMscsDir = new File(this.tmpOutputDirPath + "/splitted_mscs");
				splittedMscsDir.mkdir();
				splitMSCs(this.tmpOutputDirPath + "/multi_split_candidates.seq",this.tmpOutputDirPath + "/splitted_mscs");
				new File(this.tmpOutputDirPath + "/multi_split_candidates.seq").delete();
			}
			
			
			
			//extracting local contexts
			date = new Date();
			if(!this.strandSpecific) System.out.println(String.format("[%s]\tExtracting local contexts.",date.toLocaleString()));
			else System.out.println(String.format("[%s]\tExtracting strand specific local contexts.",date.toLocaleString()));
			ContextExtractor contextExtractor = new ContextExtractor(this.tmpOutputDirPath + "/splitted_rmap_sorted", this.indexDirPath, this.minDistanceBetweenContext,this.maxContextSize,this.minNumberOfReadsInContext,this.readLength,this.strandSpecific,this.completeWindowSize);
			HashMap<String,ArrayList<Context>> contexts;
			contexts = contextExtractor.getLocalContexts(numberOfThreads);
			int foundContexts = 0;
			
			for(String c : contexts.keySet()) {
				foundContexts += contexts.get(c).size();
				Collections.sort(contexts.get(c),new LocalContextComparator());
			}
			
			date = new Date();
			System.out.println(String.format("[%s]\tFound %s contexts.",date.toLocaleString(),foundContexts));
			if(foundContexts == 0) {
				System.out.println(String.format("[%s]\tTry to increase the minimum number of reads a genomic region has to contain for being regarded as a local context (-minsize option)",date.toLocaleString()));
				System.out.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				if(!this.keepTmp)
					deleteFolderWithContent(new File(this.tmpOutputDirPath));
				System.exit(0);
				
			}
			
			
			date = new Date();
			System.out.println(String.format("[%s]\tProcessing local contexts.",date.toLocaleString()));
			
			MappingProcessor mappingProcessor = new MappingProcessor(this.tmpOutputDirPath + "/splitted_rmap_sorted",this.referencesDir,this.annotationFilePath,this.gtfFilePath,this.readFilePath,this.readFormat,this.readLength,Integer.valueOf(this.maxMismatches),Integer.valueOf(this.maxMismatchDifference),Integer.valueOf(this.maxIntronLength), this.maxDelSize,Integer.valueOf(this.seedLength),Integer.valueOf(this.seedMissmatches),this.preferExtensionsWithKnownSpliceSignals,this.skipDenovoJunctions,this.skipNonCanonicalJunctions,this.strandSpecific,this.clipping, this.polyA, this.strandedPolyA, this.minPolyALength,this.minPolyAReadCount,this.upperPolyACutoff,this.lowerPolyACutoff,this.maxConsideredClippingLength);
			IOCoordinator ioCoordinator = new IOCoordinator(this.maxIOThreads);
			boolean chrSequenceAvailable;
			
			ExecutorService executor = Executors.newFixedThreadPool(numberOfThreads);
			futures.clear();
			Future tmpFuture;
			
			HashMap<String,IntervalTree<Microbe>> chr2intervalTree = null;
			IntervalTree<Microbe> currentIntervalTree;
			if(this.indexDirPath != null) {
				chr2intervalTree = buildIntervalTrees(this.indexDirPath);
				filterIntervalTrees(chr2intervalTree,contexts);
				
			}
			
			ArrayList<String> chrNames = new ArrayList<String>();
			chrNames.addAll(contexts.keySet());
			Collections.sort(chrNames,new GlobalContextComparator(contexts));
			
			
			String chrName;
			
			
			File[] chromosomeFiles = new File(this.referencesDir).listFiles();
			HashMap<String,File> chrName2refFile = new HashMap<String,File>();
			for(File chrFile : chromosomeFiles) {
				if(chrFile.isFile() && chrFile.getName().contains(".")) {
					chrName = chrFile.getName().substring(0,chrFile.getName().lastIndexOf("."));
					chrName2refFile.put(chrName,chrFile);
				}
			}
			
			HashMap<String,String> chrName2rmapPath = new HashMap<String,String>();
			File[] rmapFiles = new File(this.tmpOutputDirPath + "/splitted_rmap_sorted").listFiles();
			for(File f : rmapFiles) {
				chrName = f.getName().substring(0,f.getName().lastIndexOf("."));
				chrName2rmapPath.put(chrName, f.getAbsolutePath());
			}
			
			HashMap<String,String> chr2Name2MscPath = new HashMap<String,String>();
			File[] mscFiles = new File(this.tmpOutputDirPath + "/splitted_mscs").listFiles();
			if(mscFiles != null) {
				for(File f : mscFiles) {
					chrName = f.getName().substring(0,f.getName().lastIndexOf("."));
					chr2Name2MscPath.put(chrName, f.getAbsolutePath());
				}
			}
			
			//mapping chrName -> sequence
			StringBuilder sb = new StringBuilder();
			ChromosomeSequenceParser csp;
			if(chrName2refFile.containsKey(chrNames.get(0))) {
				csp = new ChromosomeSequenceParser(sb,chrName2refFile.get(chrNames.get(0)));
				csp.addListener(this);
				executor.submit(csp);
				this.threadsBufferingChromosomeSequences.incrementAndGet();
			}
			boolean chrPreBuffered;
			
			
			//mapping readId -> sequence
			ArrayList<Pair<Long,Long>> readFps = null;
			Set<String> readIds = null;
			int alreadyProcessedContexts = 0;
			int newlyBufferedContexts = contexts.get(chrNames.get(0)).size();
			ReadSequenceParser rsp = null;
			MscSequenceParser mscp = null;
			//readIds = new HashSet<String>();
			DB readIdDb = DBMaker.newMemoryDB().transactionDisable().closeOnJvmShutdown().cacheSize(1000000).make();
			readIds = readIdDb.createHashSet("readIds").counterEnable().serializer(Serializer.STRING).makeOrGet();
			
			if(!this.writeSequenceDB) {
				readFps = getReadFilePointers(this.readFilePath, this.threads);
				
				//db = DBMaker.newMemoryDB().transactionDisable().make();
				//read2sequence = db.createHashMap("read2sequence").keySerializer(Serializer.STRING).valueSerializer(Serializer.STRING).makeOrGet();
				read2sequence = new ConcurrentHashMap<String,long[]>();
				
			
				alreadyProcessedContexts = 0;
				newlyBufferedContexts = getNextReadIds(contexts.get(chrNames.get(0)),readIds, chrName2rmapPath.get(chrNames.get(0)), alreadyProcessedContexts,this.read2sequenceSize);
				for(Pair<Long,Long> fpPair : readFps) {
					rsp = new ReadSequenceParser(readIds, read2sequence, fpPair.getFirst(), fpPair.getSecond(), this.readFilePath);
					rsp.addListener(this);
					executor.submit(rsp);
					this.threadsBufferingReadSequences.incrementAndGet();
				}
				mscp = new MscSequenceParser(readIds, read2sequence, chr2Name2MscPath.get(chrNames.get(0)));
				mscp.addListener(this);
				executor.submit(mscp);
				this.threadsBufferingReadSequences.incrementAndGet();
				
				while(this.threadsBufferingReadSequences.get() != 0) {
					Thread.sleep(500);
				}
				readIds.clear();
			}
			
			//mapping chrName -> length (for SAM header)
			HashMap<String,Integer> chr2length = new HashMap<String,Integer>();
			boolean foundLargeContext = false;
			
			for(int k = 0; k < chrNames.size(); k++) {
				
				while(this.threadsBufferingChromosomeSequences.get() != 0) {
					Thread.sleep(500);
				}
				chrPreBuffered = false;
				mappingProcessor.setCurrentChromosome(sb);
				
				chrName = chrNames.get(k);
				chr2length.put(chrName,sb.length());
				
				if(chr2intervalTree != null && chr2intervalTree.containsKey(chrName))
					currentIntervalTree = chr2intervalTree.get(chrName);
				else
					currentIntervalTree = null;
				
				ArrayList<Context> currentContexts = contexts.get(chrName);
				
				if(sb.length() == 0) {
					date = new Date();
					System.err.println(String.format("[%s]\tWARNING, could not find sequence for chr: %s",date.toLocaleString(),chrName));
					if(k < chrNames.size() - 1) {
						
						sb = new StringBuilder();
						if(chrName2refFile.containsKey(chrNames.get(k+1))) {
							csp = new ChromosomeSequenceParser(sb,chrName2refFile.get(chrNames.get(k+1)));
							csp.addListener(this);
							executor.submit(csp);
							this.threadsBufferingChromosomeSequences.incrementAndGet();
						}
						
						
						if(!this.writeSequenceDB) {
							read2sequence = new ConcurrentHashMap<String,long[]>();
							//db = DBMaker.newMemoryDB().transactionDisable().make();
							//read2sequence = db.createHashMap("read2sequence").keySerializer(Serializer.STRING).valueSerializer(Serializer.STRING).makeOrGet();
							
							newlyBufferedContexts = getNextReadIds(contexts.get(chrNames.get(k+1)),readIds, chrName2rmapPath.get(chrNames.get(k+1)), 0,read2sequenceSize);
							for(Pair<Long,Long> fpPair : readFps) {
								rsp = new ReadSequenceParser(readIds, read2sequence, fpPair.getFirst(), fpPair.getSecond(), this.readFilePath);
								rsp.addListener(this);
								executor.submit(rsp);
								this.threadsBufferingReadSequences.incrementAndGet();
							}
							mscp = new MscSequenceParser(readIds, read2sequence, chr2Name2MscPath.get(chrNames.get(k+1)));
							mscp.addListener(this);
							executor.submit(mscp);
							this.threadsBufferingReadSequences.incrementAndGet();
						}
						else {
							newlyBufferedContexts = contexts.get(chrNames.get(k+1)).size();
						}
						
					}
					continue;
				}
				
				alreadyProcessedContexts = 0;
				while(alreadyProcessedContexts < currentContexts.size()) {
					
					while(this.threadsBufferingReadSequences.get() != 0) {
						Thread.sleep(500);
					}
					readIds.clear();
					mappingProcessor.setCurrentReadSequences(read2sequence);
					
					foundLargeContext = false;
					for(int i = alreadyProcessedContexts; i < alreadyProcessedContexts + newlyBufferedContexts; i++) {
						ContextProcessingCoordinator coordinator = new ContextProcessingCoordinator(mappingProcessor,ioCoordinator,currentContexts.get(i),chr2rmap.get(chrName), String.format("%s/%s_thread_%s",this.tmpOutputDirPath,chrName,i),currentIntervalTree,windowSizes,maxMismatchDifference,Integer.valueOf(maxMismatches),minStartPositionOverlaps,readLength,this.maxContextSize,this.maxGapSize,this.preferExtensionsWithKnownSpliceSignals,skipDenovoJunctions,this.skipNonCanonicalJunctions,this.updateQueue, this.pairedEnd, this.updateInterval,this.verbose,this.developer,(currentContexts.get(i).getContainedReads() > this.read2sequenceSize));
						coordinator.addListener(this);
						futures.add(executor.submit(coordinator));
						this.contextsNeedingBufferedChromosomeSequence.incrementAndGet();
						if(!this.writeSequenceDB)
							this.contextsNeedingBufferedReadSequences.incrementAndGet();
						
						if(currentContexts.get(i).getContainedReads() > this.read2sequenceSize)
							foundLargeContext = true;
					}
					
					if(foundLargeContext) {
						while(this.contextsNeedingBufferedReadSequences.get() > 0) {
							Thread.sleep(500);
						}
						
						mappingProcessor.getCurrentReadSequence().clear();
						for(int i = 0; i < futures.size(); i++) {
							tmpFuture = futures.get(i);
							tmpFuture.get();
						}
					}
					
					
					
					if(!chrPreBuffered && k < chrNames.size() - 1) {
						sb = new StringBuilder();
						if(chrName2refFile.containsKey(chrNames.get(k+1))) {
							csp = new ChromosomeSequenceParser(sb,chrName2refFile.get(chrNames.get(k+1)));
							csp.addListener(this);
							executor.submit(csp);
							this.threadsBufferingChromosomeSequences.incrementAndGet();
						}
						chrPreBuffered = true;
					}
					
					
					alreadyProcessedContexts += newlyBufferedContexts;
					
					if(alreadyProcessedContexts < currentContexts.size()) {
						if(!writeSequenceDB) {
							read2sequence = new ConcurrentHashMap<String,long[]>();
							newlyBufferedContexts = getNextReadIds(contexts.get(chrName),readIds, chrName2rmapPath.get(chrName), alreadyProcessedContexts,read2sequenceSize);
							for(Pair<Long,Long> fpPair : readFps) {
								rsp = new ReadSequenceParser(readIds, read2sequence, fpPair.getFirst(), fpPair.getSecond(), this.readFilePath);
								rsp.addListener(this);
								executor.submit(rsp);
								this.threadsBufferingReadSequences.incrementAndGet();
							}
							mscp = new MscSequenceParser(readIds, read2sequence, chr2Name2MscPath.get(chrName));
							mscp.addListener(this);
							executor.submit(mscp);
							this.threadsBufferingReadSequences.incrementAndGet();
						}
						else {
							newlyBufferedContexts = contexts.get(chrName).size();
						}
					}
					
					else if(k < chrNames.size() - 1) {
						if(!writeSequenceDB) {
							read2sequence = new ConcurrentHashMap<String,long[]>();
							//db = DBMaker.newMemoryDB().transactionDisable().make();
							//read2sequence = db.createHashMap("read2sequence").keySerializer(Serializer.STRING).valueSerializer(Serializer.STRING).makeOrGet();
							
							newlyBufferedContexts = getNextReadIds(contexts.get(chrNames.get(k+1)),readIds, chrName2rmapPath.get(chrNames.get(k+1)), 0,read2sequenceSize);
							for(Pair<Long,Long> fpPair : readFps) {
								rsp = new ReadSequenceParser(readIds, read2sequence, fpPair.getFirst(), fpPair.getSecond(), this.readFilePath);
								rsp.addListener(this);
								executor.submit(rsp);
								this.threadsBufferingReadSequences.incrementAndGet();
							}
							mscp = new MscSequenceParser(readIds, read2sequence, chr2Name2MscPath.get(chrNames.get(k+1)));
							mscp.addListener(this);
							executor.submit(mscp);
							this.threadsBufferingReadSequences.incrementAndGet();
						}
						else {
							newlyBufferedContexts = contexts.get(chrNames.get(k+1)).size();
						}
						
					}
					
					//wait until the last thread unlocks buffered read sequences
					while(this.contextsNeedingBufferedReadSequences.get() > 0) {
						Thread.sleep(500);
						
						//free mem
						if(this.threadsBufferingReadSequences.get() == 0) 
							readIds.clear();
					}
					
				}

				//wait until the last thread unlocks the actual buffered chromosome sequence
				while(this.contextsNeedingBufferedChromosomeSequence.get() > 0) {
					Thread.sleep(500);
					
					//free mem
					if(this.threadsBufferingReadSequences.get() == 0) 
						readIds.clear();
				}
			}
			executor.shutdown();
			for(int i = 0; i < futures.size(); i++) {
				tmpFuture = futures.get(i);
				tmpFuture.get();
			}
			futures.clear();
			contexts.clear();
			contexts = null;
			mappingProcessor.emptyReadSequences();
			mappingProcessor.emptyChrSequence();
			mappingProcessor = null;
			System.gc();
			
			
			if(new File(this.tmpOutputDirPath + "/resolved_local_contexts").listFiles().length == 0) {
				date = new Date();
				System.out.println(String.format("[%s]\tDid not find any alignments. Aborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			//now resolve the global context and output the most likely read mapping
			date = new Date();
			System.out.println(String.format("[%s]\tResolving global contexts.",date.toLocaleString())); 
			String allResolvedLocalContexts = this.tmpOutputDirPath + "/all_resolved_local_contexts.txt";
			
			//concatenateResolvedContexts(this.tmpOutputDirPath + "/resolved_local_contexts",allResolvedLocalContexts);
			//FileSorter7 fileSorter = new FileSorter7(allResolvedLocalContexts,allResolvedLocalContexts + ".sorted",null,new int[]{1},200, "\t",false);
			//fileSorter.sortFile(allResolvedLocalContexts, allResolvedLocalContexts + ".sorted",null);
			
			new UnixSort().merge(this.tmpOutputDirPath + "/resolved_local_contexts", allResolvedLocalContexts, this.tmpOutputDirPath +"/sort_tmp", "\t", 2, 500, false, this.verbose);
			
			
			
			if(!this.pairedEnd)
				new MappingProcessor(allResolvedLocalContexts, allResolvedLocalContexts + ".sorted.best.matchings", this.maxMismatchDifference,this.preferExtensionsWithKnownSpliceSignals,skipDenovoJunctions,skipNonCanonicalJunctions).extractBestMatchingAlignments(true);
			else
				new MappingProcessor(allResolvedLocalContexts, allResolvedLocalContexts + ".sorted.best.matchings", this.maxMismatchDifference,this.preferExtensionsWithKnownSpliceSignals,skipDenovoJunctions,skipNonCanonicalJunctions).extractBestMatchingPairedEndAlignments(this.maxContextSize, this.polyA);
			
			
			ContextResolver contextResolver;
			if(!this.pairedEnd)
				contextResolver = new GlobalContextResolverSingleEnd(allResolvedLocalContexts + ".sorted.best.matchings", this.outputDirPath + "/mapping.sam",chr2length,chr2intervalTree,this.windowSizes,this.readLength,this.maxContextSize,this.maxDelSize,false,this.updateQueue,this.updateInterval,numberOfThreads,this.printMultiMappings,printSecondBestChr,verbose,this.strandSpecific);
			else
				contextResolver = new GlobalContextResolverPairedEnd(allResolvedLocalContexts + ".sorted.best.matchings", this.outputDirPath + "/mapping.sam",chr2length,chr2intervalTree,this.windowSizes,this.readLength,this.maxContextSize,this.maxDelSize,false,this.updateQueue,this.updateInterval,numberOfThreads,this.printMultiMappings,printSecondBestChr,verbose,this.strandSpecific);
			
			contextResolver.resolve();
			
			
			
			if(this.writeSequenceDB) {
				db.close();
			}
			
			//remove tmp folder with content
			if(!this.keepTmp)
				deleteFolderWithContent(new File(this.tmpOutputDirPath));
			if(new File(this.outputDirPath + "/tmp_reads.fa").exists())
				new File(this.outputDirPath + "/tmp_reads.fa").delete();
			
			
			
			
			
			//if polyA-tail search is enabled we extract the found polyA reads here and write the result into a separate file
			if(this.polyA) {
				date = new Date();
				System.out.println(String.format("[%s]\tExtracting coordinates of found polyA tails.",date.toLocaleString()));
				extractPolyAtails(this.outputDirPath + "/mapping.sam",this.outputDirPath + "/polyA_tails.bed");
			}
			
			date = new Date();
			System.out.println(String.format("[%s]\tContextMap run finished.",date.toLocaleString()));
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private class Sequence2DBWriter extends Thread {
		
		private HTreeMap<String,String> read2sequence; 
		private String readFilePath;
		private long startFilePointer; 
		private long stopFilePointer;
		private boolean isMscFile;
		private final int maxTmpElements = 1000000;
		
		public Sequence2DBWriter(HTreeMap<String,String> read2sequence, String readFilePath, long startFilePointer, long stopFilePointer, boolean isMscFile) {
			this.read2sequence = read2sequence;
			this.readFilePath = readFilePath;
			this.startFilePointer = startFilePointer;
			this.stopFilePointer = stopFilePointer;
			this.isMscFile = isMscFile;
		}
		
		public void run() {
			try {
				BufferedRandomAccessFile br = new BufferedRandomAccessFile(new File(this.readFilePath),"r",10240);
				br.seek(this.startFilePointer);
				
				String currentLine;
				Pattern tabPattern = Pattern.compile("\t");
				String[] splittedLine;
				while((currentLine = br.getNextLine()) != null) {
					if(!this.isMscFile) {
						this.read2sequence.put(currentLine.substring(1), br.getNextLine());
					}
					else {
						splittedLine = tabPattern.split(currentLine);
						this.read2sequence.put(splittedLine[0], splittedLine[1]);
					}
					
					
					if(br.getFilePointer() == this.stopFilePointer)
						break;
				}
				br.close();
				
			}
			
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private void writeSequences2DB(HTreeMap<String,String> read2sequence, String readFilePath, boolean isMscFile, int threads) {
		try {
			

			BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(readFilePath),"r",1024);
			long fileSize = braf.getChannel().size();
			long incRate = fileSize/threads;
			long prevPosition = 0;
			long currentPosition = 0;
			
			
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			String currentLine;
			
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				braf.seek(currentPosition);
				braf.getNextLine();
				currentLine = braf.getNextLine();
				if(currentLine != null && currentLine.charAt(0) == '>')
					braf.getNextLine();
				
				currentPosition = braf.getFilePointer();
				
				futures.add(executor.submit(new Sequence2DBWriter(read2sequence,readFilePath, prevPosition, currentPosition, isMscFile)));
				
				prevPosition = currentPosition;
				
			}
			braf.close();
			
			executor.shutdown();
			for(Future f : futures) {
				f.get();
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	
private class KeyIterator implements Iterator<Fun.Tuple2<String, long[]>> {
		
		private BufferedReader br;
		private boolean isMSC;
		private final Pattern tabPattern = Pattern.compile("\t");
		private String2Bitset string2bitset;
		
		public KeyIterator(String readFilePath, boolean isMSC) {
			try {
				this.br = new BufferedReader(new FileReader(new File(readFilePath)));
				this.isMSC = isMSC;
				this.string2bitset = new String2Bitset();
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		

		@Override
		public boolean hasNext() {
			try {
				return this.br.ready();
			}
			catch(Exception e) {
				e.printStackTrace();
				return false;
			}
		}

		@Override
		public Fun.Tuple2<String, long[]> next() {
			try {
				if(!this.isMSC) {
					String id = br.readLine().substring(1);
					String sequence = br.readLine();
					Fun.Tuple2<String, long[]> tmpTuple2 = new Fun.Tuple2<String,long[]>(id,this.string2bitset.compress(sequence));
					return tmpTuple2;
				}
				else {
					String[] splittedLine = tabPattern.split(br.readLine());
					Fun.Tuple2<String, long[]> tmpTuple2 = new Fun.Tuple2<String,long[]>(splittedLine[0],this.string2bitset.compress(splittedLine[1]));
					return tmpTuple2;
				}
			}
			catch(Exception e) {
				e.printStackTrace();
				return null;
			}
		}
		
	}
	
	//output in bed format: chr	start	end	name	read_count	transcript_strand	start	end	255,0,0	1	1	0
	private void extractPolyAtails(String mappingFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(mappingFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern doublePointPattern = Pattern.compile(":");
			String[] splittedPolyAinfo;
			int alignmentStart;
			int polyAStart;
			String strand = null;
			String key;
			int flag;
			boolean isPaired = false;
			boolean isFirstMate = false;
			
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> chr2site2StartPositionsFwd = new HashMap<String,TreeMap<Integer,HashSet<Integer>>>();
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> chr2site2StartPositionsRev = new HashMap<String,TreeMap<Integer,HashSet<Integer>>>();
			HashSet<String> chrs = new HashSet<String>();
			HashMap<String,HashSet<Integer>> site2StartPositions = new HashMap<String,HashSet<Integer>>();
			HashMap<String,MutableInt> site2count = new HashMap<String,MutableInt>();
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				
				polyAStart = -1;
				for(int i = 10; i < splittedLine.length; i++) {
					if(splittedLine[i].contains("PT:i:")) {
						splittedPolyAinfo = doublePointPattern.split(splittedLine[i]);
						polyAStart = Integer.valueOf(splittedPolyAinfo[splittedPolyAinfo.length - 1]);
						break;
					}
				}
				
				if(polyAStart != -1) {
					
					isPaired = false;
					flag = Integer.valueOf(splittedLine[1]);
					if((flag & (1L << 6)) != 0) {
						isPaired = true;
						isFirstMate = true;
					}
					else if((flag & (1L << 7)) != 0) {
						isPaired = true;
						isFirstMate = false;
					}
					
					chrs.add(splittedLine[2]);
					alignmentStart = Integer.valueOf(splittedLine[3]);
					//clipped at the start -> transcript on the reverse strand
					if(polyAStart + 1 == alignmentStart) {
						alignmentStart -= Integer.valueOf(splittedLine[5].split("S")[0]);
						if(strandedPolyA && isPaired && isFirstMate) {
							strand = "+";
							if(!chr2site2StartPositionsFwd.containsKey(splittedLine[2]))
								chr2site2StartPositionsFwd.put(splittedLine[2], new TreeMap<Integer,HashSet<Integer>>());
							
							if(!chr2site2StartPositionsFwd.get(splittedLine[2]).containsKey(polyAStart))
								chr2site2StartPositionsFwd.get(splittedLine[2]).put(polyAStart, new HashSet<Integer>());
							
							chr2site2StartPositionsFwd.get(splittedLine[2]).get(polyAStart).add(alignmentStart);
						}
						
						else {
						strand = "-";
						if(!chr2site2StartPositionsRev.containsKey(splittedLine[2]))
							chr2site2StartPositionsRev.put(splittedLine[2], new TreeMap<Integer,HashSet<Integer>>(new IntComparator()));
						
						if(!chr2site2StartPositionsRev.get(splittedLine[2]).containsKey(polyAStart))
							chr2site2StartPositionsRev.get(splittedLine[2]).put(polyAStart, new HashSet<Integer>());
						
						
						chr2site2StartPositionsRev.get(splittedLine[2]).get(polyAStart).add(alignmentStart);
						}
					}
					
					//clipped at the end -> transcript on the forward strand (not in stranded mode if the first mate is found
					else {
						if(strandedPolyA && isPaired && isFirstMate) {
							strand = "-";
							if(!chr2site2StartPositionsRev.containsKey(splittedLine[2]))
								chr2site2StartPositionsRev.put(splittedLine[2], new TreeMap<Integer,HashSet<Integer>>(new IntComparator()));
							
							if(!chr2site2StartPositionsRev.get(splittedLine[2]).containsKey(polyAStart))
								chr2site2StartPositionsRev.get(splittedLine[2]).put(polyAStart, new HashSet<Integer>());
							
							
							chr2site2StartPositionsRev.get(splittedLine[2]).get(polyAStart).add(alignmentStart);
						}
						
						else {
							strand = "+";
							
							if(!chr2site2StartPositionsFwd.containsKey(splittedLine[2]))
								chr2site2StartPositionsFwd.put(splittedLine[2], new TreeMap<Integer,HashSet<Integer>>());
							
							if(!chr2site2StartPositionsFwd.get(splittedLine[2]).containsKey(polyAStart))
								chr2site2StartPositionsFwd.get(splittedLine[2]).put(polyAStart, new HashSet<Integer>());
							
							chr2site2StartPositionsFwd.get(splittedLine[2]).get(polyAStart).add(alignmentStart);
						}
					}
					
					key = String.format("%s::%s::%s",splittedLine[2],strand,polyAStart);
					if(!site2count.containsKey(key)) {
						site2StartPositions.put(key,new HashSet<Integer>());
						site2count.put(key,new MutableInt(0));
					}
					
					
					
					
					site2StartPositions.get(key).add(alignmentStart);
					site2count.get(key).increment();
					
				}
			}
			br.close();
			
			
			ArrayList<Integer> currentStarts = new ArrayList<Integer>();
			ArrayList<Integer> startsToDelete = new ArrayList<Integer>();
			int clusterSize = 25;
			int currentStartPositions = 0;
			int clusterStart;
			int clusterStop;
			
			for(String chr : chrs) {
				if(chr2site2StartPositionsFwd.containsKey(chr)) {
					startsToDelete.clear();
					currentStarts.clear();
					currentStartPositions = 0;
					clusterStart = chr2site2StartPositionsFwd.get(chr).firstKey();
					clusterStop = clusterStart + clusterSize - 1;
					for(int start : chr2site2StartPositionsFwd.get(chr).keySet()) {
						if(start >= clusterStart && start <= clusterStop) {
							currentStarts.add(start);
							currentStartPositions += chr2site2StartPositionsFwd.get(chr).get(start).size();
						}
						
						//open a new cluster
						else if(start > clusterStop) {
							if(currentStartPositions < this.minPolyAReadCount)
								startsToDelete.addAll(currentStarts);
							
							clusterStart = start;
							clusterStop = clusterStart + clusterSize - 1;
							currentStarts.clear();
							currentStarts.add(start);
							currentStartPositions = chr2site2StartPositionsFwd.get(chr).get(start).size();
						}
					}
					
					if(!currentStarts.isEmpty()) {
						if(currentStartPositions < this.minPolyAReadCount)
							startsToDelete.addAll(currentStarts);
					}
					
					for(int start : startsToDelete) {
						chr2site2StartPositionsFwd.get(chr).remove(start);
					}
				}
				
				if(chr2site2StartPositionsRev.containsKey(chr)) {
					startsToDelete.clear();
					currentStarts.clear();
					currentStartPositions = 0;
					clusterStart = chr2site2StartPositionsRev.get(chr).firstKey() - clusterSize + 1;
					clusterStop = chr2site2StartPositionsRev.get(chr).firstKey();
					
					for(int start : chr2site2StartPositionsRev.get(chr).keySet()) {
						if(start >= clusterStart && start <= clusterStop) {
							currentStarts.add(start);
							currentStartPositions += chr2site2StartPositionsRev.get(chr).get(start).size();
						}
						
						//open a new cluster
						else if(start < clusterStop) {
							if(currentStartPositions < this.minPolyAReadCount)
								startsToDelete.addAll(currentStarts);
							
							clusterStart = start - clusterSize + 1;
							clusterStop = start;
							
							currentStarts.clear();
							currentStarts.add(start);
							currentStartPositions = chr2site2StartPositionsRev.get(chr).get(start).size();
						}
					}
					
					if(!currentStarts.isEmpty()) {
						if(currentStartPositions < this.minPolyAReadCount)
							startsToDelete.addAll(currentStarts);
					}
					
					for(int start : startsToDelete) {
						chr2site2StartPositionsRev.get(chr).remove(start);
					}
				}
				
			}
			
			//output in bed format: chr	start	end	name	read_count	transcript_strand	start	end	255,0,0	1	1	0
			int polyACounter = 0;
			for(String chr : chrs) {
				if(chr2site2StartPositionsFwd.containsKey(chr)) {
					for(int start : chr2site2StartPositionsFwd.get(chr).keySet()) {
						//if(chr2site2StartPositionsFwd.get(chr).get(start).size() < this.minPolyAReadCount)
						//	continue;
						
						pw.println(chr + "\t" + (start - 1) + "\t" + start + "\tpolyA_site_" + (++polyACounter) + "\t" + site2count.get(String.format("%s::%s::%s",chr,"+",start)).intValue() + "\t" +   "+" + "\t" +(start - 1) + "\t" + start + "\t255,0,0	1	1	0");
					}
				}
				
				if(chr2site2StartPositionsRev.containsKey(chr)) {
					for(int start : chr2site2StartPositionsRev.get(chr).keySet()) {
						//if(chr2site2StartPositionsRev.get(chr).get(start).size() < this.minPolyAReadCount)
						//	continue;
						
						pw.println(chr + "\t" + (start - 1) + "\t" + start + "\tpolyA_site_" + (++polyACounter) + "\t" + site2count.get(String.format("%s::%s::%s",chr,"-",start)).intValue() + "\t" +   "-" + "\t" +(start - 1) + "\t" + start + "\t255,0,0	1	1	0");
					}
				}
			}
			
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private class IntComparator implements Comparator {

		public IntComparator() {
			
		}
		
		@Override
		public int compare(Object o1, Object o2) {
			Integer i1 = (Integer)o1;
			Integer i2 = (Integer)o2;
			
			return i2.compareTo(i1);
		}
		
	}
	
	private void splitMSCs(String inputFilePath, String outputDir) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			String currentLine;
			String [] splittedLine;
			Pattern doublePointPattern = Pattern.compile("::");
			
			HashMap<String,PrintWriter> chr2writer = new HashMap<String,PrintWriter>();
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = doublePointPattern.split(currentLine);
				if(!chr2writer.containsKey(splittedLine[2]))
					chr2writer.put(splittedLine[2], new PrintWriter(new FileWriter(new File(outputDir + "/" + splittedLine[2] + ".msc"))));
				
				chr2writer.get(splittedLine[2]).println(currentLine);
			}
			br.close();
			for(PrintWriter pw : chr2writer.values())
				pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private class ChromosomeSequenceParser extends Thread {
		
		private StringBuilder sb;
		private File chrFile;
		
		private ArrayList<ActionListener> listeners;
		
		
		public ChromosomeSequenceParser(StringBuilder sb, File chrFile) {
			this.sb = sb;
			this.chrFile = chrFile;
			this.listeners = new ArrayList<ActionListener>();
		}
		
		public void run() {
			try {
				BufferedReader br = new BufferedReader(new FileReader(this.chrFile));
				this.sb.setLength(0);
				br.readLine();
				while(br.ready()) {
					sb.append(br.readLine());
				}
				br.close();
				fireAction(new ActionEvent(this,0,"chromosome_parsed"));
			}
			
			catch(Exception e) {
				e.printStackTrace();
			}
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
	
	
	private int getNextReadIds(ArrayList<Context> contexts,Set<String> readIds, String rmapFilePath, int alreadyProcessedContexts, int maxBufferSize) {
		try {
			BufferedRandomAccessFile contextReader = new BufferedRandomAccessFile(new File(rmapFilePath),"r",1000000);
			int newlyBufferedContexts = 0;
			
			readIds.clear();
			StringTokenizer st;
			String currentLine;
			String readId;
			String[] splittedId;
			Context context;
			boolean isMsc;
			for(int i = alreadyProcessedContexts; i < contexts.size(); i++) {
				context = contexts.get(i);
				contextReader.seek(context.getPointerToFirstRead());
				while((currentLine = contextReader.getNextLine()) != null) {
					st = new StringTokenizer(currentLine,"\t");
					readId = st.nextToken();
					isMsc = readId.contains("::MSC");
					if(!readIds.contains(readId)) {
						readIds.add(readId);
					}
					
					if(isMsc) {
						splittedId = readId.split("::");
						readId = splittedId[0];
						
						if(!readIds.contains(readId)) {
							readIds.add(readId);
						}
					}
					
					if(contextReader.getFilePointer() == context.getPointerToLastRead()) {
						break;
					}
				}
				newlyBufferedContexts++;
				
				if(readIds.size() >= maxBufferSize)
					break;
			}
			contextReader.close();
			return newlyBufferedContexts;
		}
		catch(Exception e) {
			e.printStackTrace();
			return 0;
		}
	}
	
	
	private ArrayList<Pair<Long,Long>> getReadFilePointers(String readFilePath, int threads) {
		try {
			ArrayList<Pair<Long,Long>> pointers = new ArrayList<Pair<Long,Long>>();
			BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(readFilePath),"r",1024);
			long fileSize = braf.getChannel().size();
			long incRate = fileSize/threads;
			long prevPosition = 0;
			long currentPosition = 0;
			String currentLine;
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				braf.seek(currentPosition);
				braf.getNextLine();
				currentLine = braf.getNextLine();
				
				if(currentLine != null && currentLine.charAt(0) == '>')
					braf.getNextLine();
				
				currentPosition = braf.getFilePointer();
				pointers.add(new Pair(prevPosition,currentPosition));
				
				prevPosition = currentPosition;
			}
			braf.close();
			return pointers;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	
	private class MscSequenceParser extends Thread {
		private Set<String> readIds;
		private ConcurrentMap<String,long[]> read2sequence;
		private String mscFilePath;
		private String2Bitset string2bitset;
		
		private ArrayList<ActionListener> listeners;
		
		public MscSequenceParser(Set<String> readIds, ConcurrentMap<String,long[]> read2sequence, String mscFilePath) {
			this.readIds = readIds;
			this.read2sequence = read2sequence;
			this.mscFilePath = mscFilePath;
			this.listeners = new ArrayList<ActionListener>();
			this.string2bitset = new String2Bitset();
		}
		
		public void run() {
			try {
				if(this.mscFilePath != null) {
					BufferedReader br = new BufferedReader(new FileReader(new File(this.mscFilePath)));
					String currentLine;
					String[] splittedLine;
					Pattern tabPattern = Pattern.compile("\t");
					while((currentLine = br.readLine()) != null) {
						splittedLine = tabPattern.split(currentLine);
						
						if(this.readIds.contains(splittedLine[0])) {
							this.read2sequence.put(splittedLine[0],this.string2bitset.compress(splittedLine[1]));
							
						}
						
					}
					br.close();
				}
				fireAction(new ActionEvent(this,0,"sequences_parsed"));
			}
			
			catch(Exception e) {
				e.printStackTrace();
			}
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
	
	
	private class ReadSequenceParser extends Thread {
		
		private Set<String> readIds;
		private ConcurrentMap<String,long[]> read2sequence;
		private long startPointer;
		private long stopPointer;
		private String readFilePath;
		
		private ArrayList<ActionListener> listeners;
		private String2Bitset string2bitset;
		
		public ReadSequenceParser(Set<String> readIds, ConcurrentMap<String,long[]> read2sequence, long startPointer, long stopPointer, String readFilePath) {
			this.readIds = readIds;
			this.read2sequence = read2sequence;
			this.startPointer  = startPointer;
			this.stopPointer = stopPointer;
			this.readFilePath = readFilePath;
			this.string2bitset = new String2Bitset();
			this.listeners = new ArrayList<ActionListener>();
		}
		
		public void run() {
			try {
				
				BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(readFilePath),"r",10000);
				String currentLine;
				braf.seek(this.startPointer);
				while((currentLine = braf.getNextLine()) != null) {
					if(currentLine.charAt(0) == '>') {
						if(this.readIds.contains(currentLine.substring(1))) {
							this.read2sequence.put(currentLine.substring(1), string2bitset.compress(braf.getNextLine()));
						}
					}
					
					if(braf.getFilePointer() == this.stopPointer)
						break;
				}
				braf.close();
				fireAction(new ActionEvent(this,0,"sequences_parsed"));
			}
			catch(Exception e) {
				e.printStackTrace();
			}
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
	
	
	
	private int getFreeSlot(Thread[] threads) {
		for(int i = 0; i < threads.length; i++) {
			Thread thread = threads[i];
			if(thread == null)
				return i;
			else if(!thread.isAlive())
				return i;
		}
		return -1;
	}
	
	private ArrayList<Integer> getAllFreeSlots(Thread[] threads) {
		ArrayList<Integer> freeSlots = new ArrayList<Integer>();
		for(int i = 0; i < threads.length; i++) {
			Thread thread = threads[i];
			if(thread == null) 
				freeSlots.add(i);
			
			else if(!thread.isAlive())
				freeSlots.add(i);
		}
		return freeSlots;
	}
	
	private void concatenateResolvedContexts(String inputDirPath, String outputFilePath) throws Exception {
		File outputFile = new File(outputFilePath);
		if(outputFile.isFile())
			outputFile.delete();
		
		FileOutputStream fos = new FileOutputStream(outputFilePath,true);
		FileChannel writeChannel = fos.getChannel();
		RandomAccessFile rf;
		FileChannel readChannel;
		File[] tmpFiles = new File(inputDirPath).listFiles();
		long currentChannelSize;
		long transferedBytes;
		for(File tmpFile : tmpFiles) {
			rf = new RandomAccessFile(tmpFile,"r");
			readChannel = rf.getChannel();
			currentChannelSize = readChannel.size();
			transferedBytes = 0;
			while(transferedBytes < currentChannelSize) {
				transferedBytes += readChannel.transferTo(transferedBytes, readChannel.size(), writeChannel);
			}
			rf.close();
		}
		fos.close();
	}
	
	private void deleteFolderWithContent(File file) {
		if(file.isDirectory()) {
			File[] content = file.listFiles();
			for(File f : content)
				deleteFolderWithContent(f);
		}
		file.delete();
	}
	
	private class LocalContextComparator implements Comparator {
		@Override
		public int compare(Object o1, Object o2) {
			Context c1 = (Context)o1;
			Context c2 = (Context)o2;
			return Integer.valueOf(c2.getContainedReads()).compareTo(c1.getContainedReads());
		}
		
	}
	
	private class GlobalContextComparator implements Comparator {
		private HashMap<String,ArrayList<Context>> chrName2contexts;
		
		public GlobalContextComparator(HashMap<String,ArrayList<Context>> chrName2contexts) {
			this.chrName2contexts = chrName2contexts;
		}
		
		public int compare(Object o1, Object o2) {
			String chrNameA = (String)o1;
			String chrNameB = (String)o2;
			return Integer.valueOf(this.chrName2contexts.get(chrNameB).get(0).getContainedReads()).compareTo(this.chrName2contexts.get(chrNameA).get(0).getContainedReads());
		}
	}

	@Override
	public void actionPerformed(ActionEvent event) {
		
		if(event.getActionCommand().equals("sequences_parsed"))
			this.threadsBufferingReadSequences.decrementAndGet();
		
		else if(event.getActionCommand().equals("chromosome_parsed"))
			this.threadsBufferingChromosomeSequences.decrementAndGet();
		
		else if(event.getActionCommand().equals("chromosome_unlocked"))
			this.contextsNeedingBufferedChromosomeSequence.decrementAndGet();
		
		else if(event.getActionCommand().equals("reads_unlocked"))
			this.contextsNeedingBufferedReadSequences.decrementAndGet();
	}
	
	
	private HashMap<String,IntervalTree<Microbe>> buildIntervalTrees(String indexDirPath) {
		try {
			HashMap<String,IntervalTree<Microbe>> index2genomeTree = new HashMap<String,IntervalTree<Microbe>>();
			File[] indexFiles = new File(indexDirPath).listFiles();
			
			BufferedReader br;
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					IntervalTree<Microbe> genomeTree = new IntervalTree<Microbe>();
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						Microbe tmpMicrobe = new Microbe(splittedLine[0],Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]));
						genomeTree.add(tmpMicrobe);
					}
					br.close();
					index2genomeTree.put(indexFile.getName().substring(0, indexFile.getName().lastIndexOf('.')), genomeTree);
				}
			}
			return index2genomeTree;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	
	private void filterIntervalTrees(HashMap<String,IntervalTree<Microbe>> index2genomeTree, HashMap<String,ArrayList<Context>> contexts) {
		try {
			int offset;
			IntervalTree<Microbe> currentTree;
			IntervalTree<Microbe> filteredTree;
			HashSet<Microbe> spanningMicrobes = new HashSet<Microbe>();
			boolean foundMicrobe;
			ArrayList<Microbe> tmpHits = new ArrayList<Microbe>();
			for(String chr : contexts.keySet()) {
				if(index2genomeTree.containsKey(chr)) {
					currentTree = index2genomeTree.get(chr);
					filteredTree = new IntervalTree<Microbe>();
					for(Context context : contexts.get(chr)) {
						offset = 0;
						foundMicrobe = false;
						while(context.getStart() + offset <= context.getEnd()) {
							tmpHits.clear();
							tmpHits = currentTree.getIntervalsSpanning(context.getStart() + offset,tmpHits);
							if(!tmpHits.isEmpty()) {
								spanningMicrobes.add(tmpHits.get(0));
								foundMicrobe = true;
							}
							
							offset++;
						}
						
						if(!foundMicrobe) {
							throw new Exception(String.format("Could not find microbial/viral genome for context: [%s,%s-%s]. Aborting ContextMap run.",chr,context.getStart(),context.getEnd()));
						}
					}
					
					for(Microbe m : spanningMicrobes)
						filteredTree.add(m);
					
					index2genomeTree.put(chr,filteredTree);
				}
			}
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	
	
}
