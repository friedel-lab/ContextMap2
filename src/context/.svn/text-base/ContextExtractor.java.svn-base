package context;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.SortedMap;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;

import main.Context;
import main.MutableDouble;
import tools.BufferedRandomAccessFile;

/**
 * reads from sorted rmap files (sorted by start positions) and defines for every chromosome the
 * contained context. For each context we save the start and end position as well as the lines in the
 * rmap file, where the first and last reads of the each context occurs.
 * 
 * input: folder to rmap files, for each chromosome a single file (folder should not contain anything else)
 * output: context on every chr
 * @author bonfert
 *
 */
public class ContextExtractor extends Thread {

	private final String rmapDirPath;
	private final String indexDirPath;
	private final int minDistanceBetweenContexts;
	private final int maxContextSize;
	private final int minNumberOfReadsInContext;
	private final int readLength;
	private final int completeWindowSize;
	
	private File chrFile;
	private HashMap<String,ArrayList<Context>> contexts;
	private IntervalTree<Microbe> genomeTree;
	
	private boolean strandSpecific;
	

	public ContextExtractor(String rmapDirPath, String indexDirPath, int minDistanceBetweenContexts, int maxContextSize, int minNumberOfReadsInContext, int readLength,boolean strandSpecific, int completeWindowSize) {
		this.rmapDirPath = rmapDirPath;
		this.indexDirPath = indexDirPath;
		this.minDistanceBetweenContexts = minDistanceBetweenContexts;
		this.maxContextSize = maxContextSize;
		this.minNumberOfReadsInContext = minNumberOfReadsInContext;
		this.readLength = readLength;
		this.completeWindowSize = completeWindowSize;
		this.strandSpecific = strandSpecific;
	}
	
	
	public void setChrFile(File f) {
		this.chrFile = f;
	}
		
	public void setGenomeTree(IntervalTree<Microbe> genomeTree) {
		this.genomeTree = genomeTree;
	}
	
	public void setContexts(HashMap<String,ArrayList<Context>> contexts) {
		this.contexts = contexts;
	}
	
	
	public void run() {
		if(!this.strandSpecific) {
			getLocalContextForStrandUnspecificReads(this.contexts,this.chrFile);
		}
		else {
			getLocalContextForStrandSpecificReads(this.contexts,this.chrFile);
		}
			
	}
	
	public HashMap<String,ArrayList<Context>> getLocalContexts(int threads) {
		try {
			
			/* in case we have given an index directory, we want to restrict the size of a local context
			 * by the genome sizes given in the index files
			 * index files have to be named like chr1.idx or microbes_0.idx
			 * further it is important that the chromosome names in the bowtie indices are equivalent.
			 * 
		    */
			HashMap<String,IntervalTree<Microbe>> chr2genomeTree = new HashMap<String,IntervalTree<Microbe>>();
			String chr;
			if(this.indexDirPath != null) {
				chr2genomeTree = buildIntervalTrees(indexDirPath);
			}
			
			HashMap<String,ArrayList<Context>> contexts = new HashMap<String,ArrayList<Context>>();
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			for(File f : new File(this.rmapDirPath).listFiles()) {
				ContextExtractor tmpExtractor = new ContextExtractor(this.rmapDirPath,this.indexDirPath,this.minDistanceBetweenContexts,this.maxContextSize,this.minNumberOfReadsInContext,this.readLength,this.strandSpecific,this.completeWindowSize);
				tmpExtractor.setContexts(contexts);
				tmpExtractor.setChrFile(f);
				chr = f.getName().substring(0,f.getName().lastIndexOf("."));
				
				/*
				 * for the default contextmap runs we don't have an index dir from the input and
				 * therefore we always set the index file to null here (chr2genomeTree is empty per default).
				 */
				if(chr2genomeTree.containsKey(chr)) {
					tmpExtractor.setGenomeTree(chr2genomeTree.get(chr));
				}
				
				else {
					tmpExtractor.setGenomeTree(null);
				}
				
				
				futures.add(executor.submit(tmpExtractor));
			}
			executor.shutdown();
			for(Future future : futures)
				future.get();
			
			return contexts;
		}
		
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	
	private void getLocalContextForStrandUnspecificReads(HashMap<String,ArrayList<Context>> contexts, File rmapFile) {
		try {
			BufferedRandomAccessFile rmapReader;
			String fileName;
			String chr;
			int containedReads;
			long pointerToLineOfFirstRead;
			long pointerToLineOfLastRead;
			int readStart;
			int readEnd;
			String readEndAsString;
			int readLength = this.readLength;
			int overallMappingCount;
			double coverageWeight;
			String currentLine;
			String mappingType;
			StringTokenizer st;
			
			int currentContextStart = -1;
			int currentContextEnd = -1;
			boolean foundFullRead = false;
			
			String currentSpecies = null;
			String prevSpecies = null;
			ArrayList<Microbe> microbes;
			TreeMap<Integer,MutableDouble> currentCoverage = new TreeMap<Integer,MutableDouble>();
			TreeMap<Integer,MutableDouble> prevCoverage = new TreeMap<Integer,MutableDouble>();
			boolean checkDownstreamCoverage;
			
			fileName = rmapFile.getName();
			chr = fileName.substring(0,fileName.lastIndexOf("."));
			rmapReader = new BufferedRandomAccessFile(rmapFile, "r",1024 * 1024);
			//set the inital context coordinate
			pointerToLineOfFirstRead = rmapReader.getFilePointer();
			currentLine = rmapReader.getNextLine();
			pointerToLineOfLastRead = rmapReader.getFilePointer();
			st = new StringTokenizer(currentLine,"\t");
			st.nextToken();
			mappingType = st.nextToken();
			st.nextToken();
			readStart = Integer.valueOf(st.nextToken());
			readEndAsString = st.nextToken();
			//skip strand and mismatch count
			st.nextToken();
			st.nextToken();
			readLength = Integer.valueOf(st.nextToken());
			overallMappingCount = Integer.valueOf(st.nextToken());
			coverageWeight = 1.0/(double)overallMappingCount;
			
			
			
			if(mappingType.equals("F") || mappingType.equals("P")) {
				readEnd = readStart + readLength - 1;
				
				//coverage updates only for full and split reads (partial reads are probably discarded later)
				if(mappingType.equals("F")) {
					foundFullRead = true;
					
					for(int i = readStart; i <= readEnd; i++) {
						if(currentCoverage.containsKey(i))
							currentCoverage.get(i).add(coverageWeight);
						else
							currentCoverage.put(i, new MutableDouble(coverageWeight));
					}
				}
			}
			else {
				readEnd = Integer.valueOf(readEndAsString);
				
				//here we have a split. since we don't know the exact split position we approximate the coverage here
				for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
					if(currentCoverage.containsKey(i))
						currentCoverage.get(i).add(coverageWeight);
					else
						currentCoverage.put(i, new MutableDouble(coverageWeight));
				}
				for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
					if(currentCoverage.containsKey(i))
						currentCoverage.get(i).add(coverageWeight);
					else
						currentCoverage.put(i, new MutableDouble(coverageWeight));
				}
			}
			
			currentContextStart = readStart;
			currentContextEnd = readEnd;
			containedReads = 1;
			
			//genome tree is null if we either haven't an index dir in the input or an index file for the actual chr. 
			if(this.genomeTree != null) {
				microbes = this.genomeTree.getIntervalsSpanning(readEnd, new ArrayList<Microbe>());
				if(!microbes.isEmpty())
					currentSpecies = microbes.get(0).getId();
				prevSpecies = currentSpecies;
			}
			
			//now go through the whole file and define the contexts
			while((currentLine = rmapReader.getNextLine()) != null) {
				st = new StringTokenizer(currentLine,"\t");
				//read_id	mapping_type	chr	start	end
				st.nextToken();
				mappingType = st.nextToken();
				st.nextToken();
				readStart = Integer.valueOf(st.nextToken());
				readEndAsString = st.nextToken();
				//skip strand and mismatch count
				st.nextToken();
				st.nextToken();
				readLength = Integer.valueOf(st.nextToken());
				overallMappingCount = Integer.valueOf(st.nextToken());
				coverageWeight = 1.0/(double)overallMappingCount;
				
				
				if(mappingType.equals("F") || mappingType.equals("P"))
					readEnd = readStart + readLength - 1;
				else
					readEnd = Integer.valueOf(readEndAsString);
				
				if(this.genomeTree != null) {
					microbes = this.genomeTree.getIntervalsSpanning(readEnd, new ArrayList<Microbe>());
					if(!microbes.isEmpty()) {
						currentSpecies = microbes.get(0).getId();
					}
				}
				
				//current read extends the actual context
				if(readStart - currentContextEnd <= this.minDistanceBetweenContexts &&
				(readEnd - currentContextStart + 1) <= this.maxContextSize && 
				(this.genomeTree == null || currentSpecies.equals(prevSpecies))) {
					if(readEnd > currentContextEnd)
						currentContextEnd = readEnd;
					containedReads++;
					if(mappingType.equals("F")) {
						foundFullRead = true;
						
						//update coverage
						for(int i = readStart; i <= readEnd; i++) {
							if(currentCoverage.containsKey(i))
								currentCoverage.get(i).add(coverageWeight);
							else
								currentCoverage.put(i, new MutableDouble(coverageWeight));
						}
					}
					else if(!mappingType.equals("P")) {
						for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
							if(currentCoverage.containsKey(i))
								currentCoverage.get(i).add(coverageWeight);
							else
								currentCoverage.put(i, new MutableDouble(coverageWeight));
						}
						for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
							if(currentCoverage.containsKey(i))
								currentCoverage.get(i).add(coverageWeight);
							else
								currentCoverage.put(i, new MutableDouble(coverageWeight));
						}
					}
				}
					
				//read is too far away from actual context end or the context would get too long, spawn a new context
				else {
					//check if enough reads are in the actual context and if at least one full read has been mapped to this context
					//if(containedReads >= this.minNumberOfReadsInContext && foundFullRead) {
					if(containedReads >= this.minNumberOfReadsInContext) {
						Context context = new Context(chr,currentContextStart,currentContextEnd, "both");
						context.setContainedReads(containedReads);
						context.setPointerToFirstRead(pointerToLineOfFirstRead);
						context.setPointerToLastRead(pointerToLineOfLastRead);
						
						checkDownstreamCoverage = true;
						if(this.genomeTree != null && !currentSpecies.equals(prevSpecies))
							checkDownstreamCoverage = false;
						
						setUpAndDownstreamCoverages(context,rmapReader,prevCoverage,this.completeWindowSize,checkDownstreamCoverage,prevSpecies,this.genomeTree,"both");
						
						if(contexts.containsKey(chr))
							contexts.get(chr).add(context);
						else {
							ArrayList<Context> tmpList = new ArrayList<Context>();
							tmpList.add(context);
							contexts.put(chr, tmpList);
						}
					}
					
					pointerToLineOfFirstRead = pointerToLineOfLastRead;
					currentContextStart = readStart;
					currentContextEnd = readEnd;
					containedReads = 1;
					prevCoverage.clear();
					prevCoverage.putAll(currentCoverage);
					currentCoverage.clear();
					
					if(mappingType.equals("F")) {
						foundFullRead = true;
						for(int i = readStart; i <= readEnd; i++) {
							if(currentCoverage.containsKey(i))
								currentCoverage.get(i).add(coverageWeight);
							else
								currentCoverage.put(i, new MutableDouble(coverageWeight));
						}
					}
					else {
						foundFullRead = false;
						if(!mappingType.equals("P")) {
							for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
								if(currentCoverage.containsKey(i))
									currentCoverage.get(i).add(coverageWeight);
								else
									currentCoverage.put(i, new MutableDouble(coverageWeight));
							}
							for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
								if(currentCoverage.containsKey(i))
									currentCoverage.get(i).add(coverageWeight);
								else
									currentCoverage.put(i, new MutableDouble(coverageWeight));
							}
						}
					}
					
					if(this.genomeTree != null) {
						//microbes = this.genomeTree.getIntervalsSpanning(readEnd, new ArrayList<Microbe>());
						//if(!microbes.isEmpty())
							//currentSpecies = microbes.get(0).getId();
						
						if(!currentSpecies.equals(prevSpecies))
							prevCoverage.clear();
						
						prevSpecies = currentSpecies;
					}
				}
				
				pointerToLineOfLastRead = rmapReader.getFilePointer();
			}
			
			//add the last context
			//if(containedReads >= this.minNumberOfReadsInContext && foundFullRead) {
			if(containedReads >= this.minNumberOfReadsInContext) {
				Context context = new Context(chr,currentContextStart,currentContextEnd, "both");
				context.setContainedReads(containedReads);
				context.setPointerToFirstRead(pointerToLineOfFirstRead);
				context.setPointerToLastRead(pointerToLineOfLastRead);
				setUpAndDownstreamCoverages(context,rmapReader,prevCoverage,this.completeWindowSize,false,null,null,"both");
				if(contexts.containsKey(chr))
					contexts.get(chr).add(context);
				else {
					ArrayList<Context> tmpList = new ArrayList<Context>();
					tmpList.add(context);
					contexts.put(chr, tmpList);
				}
			}
			
			rmapReader.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void setUpAndDownstreamCoverages(Context context, BufferedRandomAccessFile br, TreeMap<Integer,MutableDouble> prevCoverage,int completeWindowSize,boolean checkDownstreamCoverage,String prevSpecies, IntervalTree<Microbe> genomeTree, String strand) {
		try {
			//setting upstream coverages
			context.setUpstreamCoverage(prevCoverage.subMap(context.getStart() - completeWindowSize,context.getStart()));
			
			//setting downstream coverage
			if(checkDownstreamCoverage) {
				long initialPointer = br.getFilePointer();
				String currentLine;
				StringTokenizer st;
				String mappingType;
				int readStart;
				int readEnd;
				String readEndAsString;
				int readLength = this.readLength;
				int overallMappingCount;
				double coverageWeight;
				String currentStrand;
				br.seek(context.getPointerToLastRead());
				TreeMap<Integer,MutableDouble> downstreamCoverage = new TreeMap<Integer,MutableDouble>();
				ArrayList<Microbe> microbes;
				String currentSpecies;
				while((currentLine = br.getNextLine()) != null) {
					st = new StringTokenizer(currentLine,"\t");
					//read_id	mapping_type	chr	start	end	strand	mismatches	(readLength)
					st.nextToken();
					mappingType = st.nextToken();
					st.nextToken();
					readStart = Integer.valueOf(st.nextToken());
					readEndAsString = st.nextToken();
					currentStrand = st.nextToken();
					//skip mismatch count
					st.nextToken();
					readLength = Integer.valueOf(st.nextToken());
					overallMappingCount = Integer.valueOf(st.nextToken());
					coverageWeight = 1.0/(double)overallMappingCount;
					
					if(readStart < context.getEnd())
						continue;
					
					if(readStart > (context.getEnd() + completeWindowSize)) {
						break;
					}
					
					if(this.genomeTree != null) {
						microbes = this.genomeTree.getIntervalsSpanning(readStart, new ArrayList<Microbe>());
						if(!microbes.isEmpty()) {
							if(!microbes.get(0).getId().equals(prevSpecies))
								break;
						}
					}
					
					
					if(mappingType.equals("F")) {
						readEnd = readStart + readLength - 1;
						

						if(strand.equals("both") || currentStrand.equals(strand)) {
							for(int i = readStart; i <= readEnd; i++) {
								if(downstreamCoverage.containsKey(i))
									downstreamCoverage.get(i).add(coverageWeight);
								else
									downstreamCoverage.put(i, new MutableDouble(coverageWeight));
							}
						}
					}
					else if(!mappingType.equals("P")) {
						readEnd = Integer.valueOf(readEndAsString);
						//here we have a split. since we don't know the exact split position we approximate the coverage here
						if(strand.equals("both") || currentStrand.equals(strand)) {
							for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
								if(downstreamCoverage.containsKey(i))
									downstreamCoverage.get(i).add(coverageWeight);
								else
									downstreamCoverage.put(i, new MutableDouble(coverageWeight));
							}
							
							
							for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
								if(i > (context.getEnd() + completeWindowSize)) {
									break;
								}
								if(downstreamCoverage.containsKey(i))
									downstreamCoverage.get(i).add(coverageWeight);
								else
									downstreamCoverage.put(i, new MutableDouble(coverageWeight));
							}
						}
					}
				}
				br.seek(initialPointer);
				context.setDownstreamCoverage(downstreamCoverage.subMap(context.getEnd() + 1, context.getEnd() + completeWindowSize));
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	

	
	
	//TODO re-implement this function to remove duplicated code.
	//TODO prevent contexts from overlapping between different genomes (microbial indices) -> DONE. CHECK THIS!
	//TODO set up and downstream coverages -> DONE. CHECK THIS!
	public void getLocalContextForStrandSpecificReads(HashMap<String,ArrayList<Context>> contexts, File rmapFile) {
		try {
			BufferedRandomAccessFile br;
			String fileName;
			String chr;
			String strand;
			int forwardContainedReads = 0;
			int reverseContainedReads = 0;
			long currentPointer;
			long pointerToLineOfFirstForwardRead = Long.MIN_VALUE;
			long pointerToLineOfLastForwardRead = Long.MIN_VALUE;
			long pointerToLineOfFirstReverseRead = Long.MIN_VALUE;
			long pointerToLineOfLastReverseRead = Long.MIN_VALUE;
			
			TreeMap<Integer,MutableDouble> currentCoverageForwardStrand = new TreeMap<Integer,MutableDouble>();
			TreeMap<Integer,MutableDouble> currentCoverageReverseStrand = new TreeMap<Integer,MutableDouble>();
			TreeMap<Integer,MutableDouble> prevCoverageForwardStrand = new TreeMap<Integer,MutableDouble>();
			TreeMap<Integer,MutableDouble> prevCoverageReverseStrand = new TreeMap<Integer,MutableDouble>();
			
			int readStart;
			int readEnd;
			int readLength = this.readLength;
			int overallMappingCount;
			double coverageWeight;
			String readEndAsString;
			
			String currentLine = "";
			String mappingType;
			StringTokenizer st;
			
			int currentContextStart = -1;
			int currentContextEnd = -1;
			
			int currentForwardContextStart = -1;
			int currentForwardContextEnd = -1;
			int currentReverseContextStart = -1;
			int currentReverseContextEnd = -1;
			
			boolean foundForwardFullRead = false;
			boolean foundReverseFullRead = false;
			boolean checkDownstreamCoverage;
			
			String currentSpecies = null;
			String prevSpecies = null;
			String prevSpeciesForwardStrand = null;
			String prevSpeciesReverseStrand = null;
			ArrayList<Microbe> microbes;

			fileName = rmapFile.getName();
			chr = fileName.substring(0,fileName.lastIndexOf("."));
			
			//set the inital context coordinates
			br = new BufferedRandomAccessFile(rmapFile, "r",1024 * 1024);
			while(pointerToLineOfFirstForwardRead == Long.MIN_VALUE || pointerToLineOfFirstReverseRead == Long.MIN_VALUE) {
				currentPointer = br.getFilePointer();
				currentLine = br.getNextLine();
				
				if(currentLine == null)
					break;
				
				st = new StringTokenizer(currentLine,"\t");
				//read_id	type	chr	start	end	strand
				st.nextToken();
				mappingType = st.nextToken();
				st.nextToken();
				readStart = Integer.valueOf(st.nextToken());
				readEndAsString = st.nextToken();
				strand = st.nextToken();
				st.nextToken();
				readLength = Integer.valueOf(st.nextToken());
				overallMappingCount = Integer.valueOf(st.nextToken());
				coverageWeight = 1.0/(double)overallMappingCount;
				
				if(mappingType.equals("F") || mappingType.equals("P")) {
					readEnd = readStart + readLength - 1;
				}
				else
					readEnd = Integer.valueOf(readEndAsString);
				
				if(strand.equals("+")) {
					pointerToLineOfFirstForwardRead = currentPointer;
					pointerToLineOfLastForwardRead = br.getFilePointer();
					currentForwardContextStart = readStart;
					currentForwardContextEnd = readEnd;
					
					//genome tree is null if we either haven't an index dir in the input or an index file for the actual chr. 
					if(this.genomeTree != null) {
						microbes = this.genomeTree.getIntervalsSpanning(readEnd, new ArrayList<Microbe>());
						if(!microbes.isEmpty())
							prevSpeciesForwardStrand = microbes.get(0).getId();
					}
				}
				else {
					pointerToLineOfFirstReverseRead = currentPointer;
					pointerToLineOfLastReverseRead = br.getFilePointer();
					currentReverseContextStart = readStart;
					currentReverseContextEnd = readEnd;
					
					//genome tree is null if we either haven't an index dir in the input or an index file for the actual chr. 
					if(this.genomeTree != null) {
						microbes = this.genomeTree.getIntervalsSpanning(readEnd, new ArrayList<Microbe>());
						if(!microbes.isEmpty())
							prevSpeciesReverseStrand = microbes.get(0).getId();
					}
				}
			}
			br.close();
			
			br = new BufferedRandomAccessFile(rmapFile, "r",1024 * 1024);
			//now go through the whole file and define the contexts
			while((currentLine = br.getNextLine()) != null) {
				st = new StringTokenizer(currentLine,"\t");
				//read_id	mapping_type	chr	start	end
				st.nextToken();
				mappingType = st.nextToken();
				st.nextToken();
				readStart = Integer.valueOf(st.nextToken());
				readEndAsString = st.nextToken();
				strand = st.nextToken();
				st.nextToken();
				readLength = Integer.valueOf(st.nextToken());
				overallMappingCount = Integer.valueOf(st.nextToken());
				coverageWeight = 1.0/(double)overallMappingCount;
				
				
				if(mappingType.equals("F") || mappingType.equals("P")) {
					readEnd = readStart + readLength - 1;
				}
			
				else
					readEnd = Integer.valueOf(readEndAsString);
				
				if(strand.equals("+")) {
					currentContextStart = currentForwardContextStart;
					currentContextEnd = currentForwardContextEnd;
					
					if(this.genomeTree != null) {
						microbes = this.genomeTree.getIntervalsSpanning(readEnd, new ArrayList<Microbe>());
						if(!microbes.isEmpty()) {
							currentSpecies = microbes.get(0).getId();
						}
						prevSpecies = prevSpeciesForwardStrand;
					}
				}
				else {
					currentContextStart = currentReverseContextStart;
					currentContextEnd = currentReverseContextEnd;
					
					if(this.genomeTree != null) {
						microbes = this.genomeTree.getIntervalsSpanning(readEnd, new ArrayList<Microbe>());
						if(!microbes.isEmpty()) {
							currentSpecies = microbes.get(0).getId();
						}
						prevSpecies = prevSpeciesReverseStrand;
					}
				}
				
				//current read extends the actual context
				if(readStart - currentContextEnd <= this.minDistanceBetweenContexts &&
				(readEnd - currentContextStart + 1) <= this.maxContextSize && 
				(this.genomeTree == null || currentSpecies.equals(prevSpecies))) {
					if(strand.equals("+")) {
						if(readEnd > currentContextEnd)
							currentForwardContextEnd = readEnd;
						forwardContainedReads++;
						if(mappingType.equals("F")) {
							foundForwardFullRead = true;
							
							//coverage updates only for full and split reads (partial reads are probably discarded later)
							for(int i = readStart; i <= readEnd; i++) {
								if(currentCoverageForwardStrand.containsKey(i))
									currentCoverageForwardStrand.get(i).add(coverageWeight);
								else
									currentCoverageForwardStrand.put(i, new MutableDouble(coverageWeight));
							}
						}
						
						else if(!mappingType.equals("P")) {
							for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
								if(currentCoverageForwardStrand.containsKey(i))
									currentCoverageForwardStrand.get(i).add(coverageWeight);
								else
									currentCoverageForwardStrand.put(i, new MutableDouble(coverageWeight));
							}
							for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
								if(currentCoverageForwardStrand.containsKey(i))
									currentCoverageForwardStrand.get(i).add(coverageWeight);
								else
									currentCoverageForwardStrand.put(i, new MutableDouble(coverageWeight));
							}
						}
						
						pointerToLineOfLastForwardRead = br.getFilePointer();
					}
					else {
						if(readEnd > currentContextEnd)
							currentReverseContextEnd = readEnd;
						reverseContainedReads++;
						if(mappingType.equals("F")) {
							foundReverseFullRead = true;
							
							//coverage updates only for full and split reads (partial reads are probably discarded later)
							for(int i = readStart; i <= readEnd; i++) {
								if(currentCoverageReverseStrand.containsKey(i))
									currentCoverageReverseStrand.get(i).add(coverageWeight);
								else
									currentCoverageReverseStrand.put(i, new MutableDouble(coverageWeight));
							}
						}
						
						else if(!mappingType.equals("P")) {
							for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
								if(currentCoverageReverseStrand.containsKey(i))
									currentCoverageReverseStrand.get(i).add(coverageWeight);
								else
									currentCoverageReverseStrand.put(i, new MutableDouble(coverageWeight));
							}
							for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
								if(currentCoverageReverseStrand.containsKey(i))
									currentCoverageReverseStrand.get(i).add(coverageWeight);
								else
									currentCoverageReverseStrand.put(i, new MutableDouble(coverageWeight));
							}
						}
						
						pointerToLineOfLastReverseRead = br.getFilePointer();
					}
				}
					
				//read is too far away from actual context end or the context would get too long, spawn a new context
				else {
					if(strand.equals("+")) {
						//check if enough reads are in the actual context and if at least one full read has been mapped to this context
						//if(forwardContainedReads >= this.minNumberOfReadsInContext && foundForwardFullRead) {						
						if(forwardContainedReads >= this.minNumberOfReadsInContext) {
							Context context = new Context(chr,currentForwardContextStart,currentForwardContextEnd, "+");
							context.setContainedReads(forwardContainedReads);
							context.setPointerToFirstRead(pointerToLineOfFirstForwardRead);
							context.setPointerToLastRead(pointerToLineOfLastForwardRead);
							
							checkDownstreamCoverage = true;
							if(this.genomeTree != null && !currentSpecies.equals(prevSpecies))
								checkDownstreamCoverage = false;
							
							setUpAndDownstreamCoverages(context,br,prevCoverageForwardStrand,this.completeWindowSize,checkDownstreamCoverage,prevSpecies,this.genomeTree,"+");
							
							
							if(contexts.containsKey(chr))
								contexts.get(chr).add(context);
							else {
								ArrayList<Context> tmpList = new ArrayList<Context>();
								tmpList.add(context);
								contexts.put(chr, tmpList);
							}
						}
						
						pointerToLineOfFirstForwardRead = pointerToLineOfLastForwardRead;
						pointerToLineOfLastForwardRead = br.getFilePointer();
						currentForwardContextStart = readStart;
						currentForwardContextEnd = readEnd;
						forwardContainedReads = 1;
						prevCoverageForwardStrand.clear();
						prevCoverageForwardStrand.putAll(currentCoverageForwardStrand);
						currentCoverageForwardStrand.clear();
						if(mappingType.equals("F")) {
							foundForwardFullRead = true;
							for(int i = readStart; i <= readEnd; i++) {
								if(currentCoverageForwardStrand.containsKey(i))
									currentCoverageForwardStrand.get(i).add(coverageWeight);
								else
									currentCoverageForwardStrand.put(i, new MutableDouble(coverageWeight));
							}
						}
						else {
							foundForwardFullRead = false;
							if(!mappingType.equals("P")) {
								for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
									if(currentCoverageForwardStrand.containsKey(i))
										currentCoverageForwardStrand.get(i).add(coverageWeight);
									else
										currentCoverageForwardStrand.put(i, new MutableDouble(coverageWeight));
								}
								for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
									if(currentCoverageForwardStrand.containsKey(i))
										currentCoverageForwardStrand.get(i).add(coverageWeight);
									else
										currentCoverageForwardStrand.put(i, new MutableDouble(coverageWeight));
								}
							}
						}
						if(this.genomeTree != null) {							
							if(!currentSpecies.equals(prevSpecies))
								prevCoverageForwardStrand.clear();
							
							prevSpeciesForwardStrand = currentSpecies;
						}
					}
					else {
						//check if enough reads are in the actual context and if at least one full read has been mapped to this context
						//if(reverseContainedReads >= this.minNumberOfReadsInContext && foundReverseFullRead) {
						if(reverseContainedReads >= this.minNumberOfReadsInContext) {
							Context context = new Context(chr,currentReverseContextStart,currentReverseContextEnd, "-");
							context.setContainedReads(reverseContainedReads);
							context.setPointerToFirstRead(pointerToLineOfFirstReverseRead);
							context.setPointerToLastRead(pointerToLineOfLastReverseRead);
							checkDownstreamCoverage = true;
							if(this.genomeTree != null && !currentSpecies.equals(prevSpecies))
								checkDownstreamCoverage = false;
							
							setUpAndDownstreamCoverages(context,br,prevCoverageReverseStrand,this.completeWindowSize,checkDownstreamCoverage,prevSpecies,this.genomeTree,"-");
							
							if(contexts.containsKey(chr))
								contexts.get(chr).add(context);
							else {
								ArrayList<Context> tmpList = new ArrayList<Context>();
								tmpList.add(context);
								contexts.put(chr, tmpList);
							}
						}
						
						pointerToLineOfFirstReverseRead = pointerToLineOfLastReverseRead;
						pointerToLineOfLastReverseRead = br.getFilePointer();
						currentReverseContextStart = readStart;
						currentReverseContextEnd = readEnd;
						reverseContainedReads = 1;
						prevCoverageReverseStrand.clear();
						prevCoverageReverseStrand.putAll(currentCoverageForwardStrand);
						currentCoverageReverseStrand.clear();
						if(mappingType.equals("F")) {
							foundReverseFullRead = true;
							
							for(int i = readStart; i <= readEnd; i++) {
								if(currentCoverageReverseStrand.containsKey(i))
									currentCoverageReverseStrand.get(i).add(coverageWeight);
								else
									currentCoverageReverseStrand.put(i, new MutableDouble(coverageWeight));
							}
						}
						else {
							foundReverseFullRead = false;
							
							if(!mappingType.equals("P")) {
								for(int i = readStart; i <= (readStart + (readLength/2) - 1);i++) {
									if(currentCoverageReverseStrand.containsKey(i))
										currentCoverageReverseStrand.get(i).add(coverageWeight);
									else
										currentCoverageReverseStrand.put(i, new MutableDouble(coverageWeight));
								}
								for(int i = (readEnd - (readLength/2) + 1); i <= readEnd;i++) {
									if(currentCoverageReverseStrand.containsKey(i))
										currentCoverageReverseStrand.get(i).add(coverageWeight);
									else
										currentCoverageReverseStrand.put(i, new MutableDouble(coverageWeight));
								}
							}
						}
						if(this.genomeTree != null) {							
							if(!currentSpecies.equals(prevSpecies))
								prevCoverageReverseStrand.clear();
							
							prevSpeciesReverseStrand = currentSpecies;
						}
					}
				}
			}
			
			//add the last contexts
			
			//if(forwardContainedReads >= this.minNumberOfReadsInContext && foundForwardFullRead) {
			if(forwardContainedReads >= this.minNumberOfReadsInContext) {
				Context context = new Context(chr,currentForwardContextStart,currentForwardContextEnd, "+");
				context.setContainedReads(forwardContainedReads);
				context.setPointerToFirstRead(pointerToLineOfFirstForwardRead);
				context.setPointerToLastRead(pointerToLineOfLastForwardRead);
				setUpAndDownstreamCoverages(context,br,prevCoverageForwardStrand,this.completeWindowSize,false,null,null,"+");
				
				if(contexts.containsKey(chr))
					contexts.get(chr).add(context);
				else {
					ArrayList<Context> tmpList = new ArrayList<Context>();
					tmpList.add(context);
					contexts.put(chr, tmpList);
				}
			}
			
			//if(reverseContainedReads >= this.minNumberOfReadsInContext && foundReverseFullRead) {
			if(reverseContainedReads >= this.minNumberOfReadsInContext) {
				Context context = new Context(chr,currentReverseContextStart,currentReverseContextEnd, "-");
				context.setContainedReads(reverseContainedReads);
				context.setPointerToFirstRead(pointerToLineOfFirstReverseRead);
				context.setPointerToLastRead(pointerToLineOfLastReverseRead);
				setUpAndDownstreamCoverages(context,br,prevCoverageReverseStrand,this.completeWindowSize,false,null,null,"-");
				if(contexts.containsKey(chr))
					contexts.get(chr).add(context);
				else {
					ArrayList<Context> tmpList = new ArrayList<Context>();
					tmpList.add(context);
					contexts.put(chr, tmpList);
				}
			}
			
			br.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private HashMap<String,IntervalTree<Microbe>> buildIntervalTrees(String indexDirPath) throws Exception {
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
	
	public void printContexts(HashMap<String,ArrayList<Context>> contexts) {
		for(String chr : contexts.keySet()) {
			for(Context context : contexts.get(chr)) {
				System.out.println(String.format("%s\t%s\t%s\t%s\t%s\t%s",context.getChr(),context.getStrand(),context.getStart(),context.getEnd(),context.getPointerToFirstRead(),context.getPointerToLastRead()));
			}
		}
	}
	
	
	private void printContextInformation(Context context, String currentSpecies, String prevSpecies) {
		int upstreamCovStart = -1;
		int upstreamCovEnd = -1;
		if(!context.getUpstreamCoverage().isEmpty()) {
			upstreamCovStart = context.getUpstreamCoverage().firstKey();
			upstreamCovEnd = context.getUpstreamCoverage().lastKey();
		}
		
		int downstreamCovStart = -1;
		int downstreamCovEnd = -1;
		if(!context.getDownstreamCoverage().isEmpty()) {
			downstreamCovStart = context.getDownstreamCoverage().firstKey();
			downstreamCovEnd = context.getDownstreamCoverage().lastKey();
		}
		
		System.out.println("##---------------------------------------------------------------------##");
		System.out.println(String.format("context: %s\tstart: %s\tend: %s\tchr: %s\tcontext_species: %s\tnext_species: %s\tupstream_cov_start: %s\tupstream_cov_end: %s\tupstream_positions: %s\tdownstream_cov_start: %s\tdownstream_cov_end: %s\tdownstream_positions: %s",context.getId(),context.getStart(),context.getEnd(),context.getChr(),prevSpecies,currentSpecies,upstreamCovStart,upstreamCovEnd,context.getUpstreamCoverage().size(),downstreamCovStart,downstreamCovEnd,context.getDownstreamCoverage().size()));
		System.out.println();
		System.out.println("##---------------------------------------------------------------------##");
		System.out.println();
	}
	
private static class Microbe implements Interval {
		
		private String id;
		private int start;
		private int stop;
		
		public Microbe(String id, int start, int end) {
			this.id = id;
			this.start = start;
			this.stop = end;
			
		}
		
		public String getId() {
			return this.id;
		}
		
		public int getStart() {
			return this.start;
		}
		
		public int getStop() {
			return this.stop;
		}
		
		public int getLength() {
			return this.stop - this.start;
		}
		
	}
	
}
