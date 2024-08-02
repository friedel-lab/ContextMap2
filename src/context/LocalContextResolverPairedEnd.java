package context;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.NavigableMap;
import java.util.SortedMap;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.regex.Pattern;

import augmentedTree.IntervalTree;
import main.Microbe;
import main.MutableDouble;
import main.Pair;
import main.ReadLocation;
import main.ReadPair;
import main.SparseRead;
import main.MultiSparseReadLocation;
import main.UnsynchronizedBufferedWriter;
import tools.BufferedRandomAccessFile;
import tools.MaxPriorityQueue;

public class LocalContextResolverPairedEnd implements ContextResolver {
	
	private String multiMappingFilePath;
	private String outputPath;
	
	private SortedMap<Integer,MutableDouble> upstreamCoverage;
	private SortedMap<Integer,MutableDouble> downstreamCoverage;
	private IntervalTree<Microbe> intervalTree;
	private ArrayList<Integer> windowSizes;
	private int completeWindowIntervall;
	private int maxContextSize;
	private PairScoreComparator pairScoreComparator;
	private SparseLocationScoreComparator sparseReadLocationComparator;
	private final int maxAllowedValidPairs = 1000;
	private boolean updateQueue;
	private boolean developer;
	private int updateInterval;
	private boolean verbose;
	
	private boolean isLargeContext;
	
	private ArrayList<ReadPair<SparseRead,SparseRead>> uniquelyMappingReads;
	private ArrayList<ReadPair<SparseRead,SparseRead>> multiMappingReadsUpdateNeeded;
	private ArrayList<ReadPair<SparseRead,SparseRead>> multiMappingReads;

	public LocalContextResolverPairedEnd(String multiMappingFilePath, String outputPath,TreeMap<Integer,MutableDouble> upstreamCoverage, TreeMap<Integer,MutableDouble> downstreamCoverage, IntervalTree<Microbe> intervalTree, ArrayList<Integer> windowSizes, int readLength, int maxContextSize,boolean updateQueue, int updateInterval, boolean verbose, boolean developer, boolean isLargeContext) {
		this.multiMappingFilePath = multiMappingFilePath;
		this.windowSizes = windowSizes;
		this.completeWindowIntervall = 0;
		for(int i = 0; i < windowSizes.size(); i++) {
			this.completeWindowIntervall += windowSizes.get(i);
		}
		this.outputPath = outputPath;
		this.maxContextSize = maxContextSize;
		
		
		this.uniquelyMappingReads = new ArrayList<ReadPair<SparseRead,SparseRead>>();
		this.multiMappingReads = new ArrayList<ReadPair<SparseRead,SparseRead>>();
		this.multiMappingReadsUpdateNeeded = new ArrayList<ReadPair<SparseRead,SparseRead>>();
		this.pairScoreComparator = new PairScoreComparator();
		this.sparseReadLocationComparator = new SparseLocationScoreComparator();
		
		this.updateQueue = updateQueue;
		this.updateInterval = updateInterval;
		
		
		this.verbose = verbose;
		this.developer = developer;
		
		this.isLargeContext = isLargeContext;
		
		this.upstreamCoverage = upstreamCoverage;
		this.downstreamCoverage = downstreamCoverage;
		this.intervalTree = intervalTree;
	}
	
	public void resolve() {
		try {
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(this.outputPath)),102400);
			BufferedRandomAccessFile bufferedReader = new BufferedRandomAccessFile(new File(this.multiMappingFilePath),"r",10240);
			Pair<Integer,Integer> contextRange = getContextRange(this.multiMappingFilePath);
			int coverageOffsetIndex = contextRange.getFirst() - this.completeWindowIntervall;
			this.upstreamCoverage = this.upstreamCoverage.subMap(contextRange.getFirst() - this.completeWindowIntervall, contextRange.getFirst());
			this.downstreamCoverage = this.downstreamCoverage.subMap(contextRange.getSecond() + 1, contextRange.getSecond() + 1 +  this.completeWindowIntervall);
			
			//init the coverage array
			MutableDouble[] coverage = new MutableDouble[(contextRange.getSecond() - contextRange.getFirst() + 1) + (2* this.completeWindowIntervall)];
			for(int i = 0; i < coverage.length; i++) {
				coverage[i] = new MutableDouble();
			}
			
			//getting the coverage
			getCoverage(coverage,coverageOffsetIndex,this.uniquelyMappingReads,this.multiMappingReads,this.multiMappingReadsUpdateNeeded,bufferedReader,pw);

			
			//we determine for every read the best location (in terms of support score) and update the coverage values by discarding all other alignments of that read
			MutableDouble[] coverageUpdated = new MutableDouble[coverage.length];
			for(int i = 0; i < coverage.length; i++) {
				coverageUpdated[i] = new MutableDouble(coverage[i].doubleValue());
			}
			//determines a score for every location, a score for the read and updates the coverage values
			determineBestMapping(coverage,coverageUpdated,coverageOffsetIndex,this.multiMappingReads,bufferedReader,pw,false,false);
			
			
			//now the updated coverage values are contained in 'coverageUpdated'. we re-determine scores for global resolution and print the alignments
			determineBestMapping(coverageUpdated,null,coverageOffsetIndex,this.uniquelyMappingReads,bufferedReader,pw,true,true);
			determineBestMapping(coverageUpdated,null,coverageOffsetIndex,this.multiMappingReads,bufferedReader,pw,false,true);
			
			
			bufferedReader.close();
			pw.flush();
			pw.close();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private void processReadPairItem(SparseRead read, MultiSparseReadLocation topLocation, MultiSparseReadLocation secondBestLocation, MutableDouble[] coverage,int coverageOffsetIndex) {
		double coverageWeight;
		for(MultiSparseReadLocation tmpLocation : read.getLocations()) {
			coverageWeight = 1.0/(double)tmpLocation.getAmbiguityFactor();
			
			tmpLocation.setAmbiguityFactor(0);
			
			
			//update coverages
			for(Pair<Integer,Integer> coordinate : tmpLocation.getCoordinates()) {
				if(coordinate.getFirst() <= 0)
					break;
				
				for(int i = coordinate.getFirst(); i <= coordinate.getSecond(); i++) {
					coverage[i - coverageOffsetIndex].subtract(coverageWeight);
					
					if(coverage[i - coverageOffsetIndex].toDouble() < 0)
						coverage[i - coverageOffsetIndex].setValue(0.0);
				}
			}
			
		}
		
		
		//topLocation.setAmbiguityFactor(1 + (topLocation.getAmbiguityFactor() - read.getLocations().size()));
		topLocation.setAmbiguityFactor(read.getOverallMappingCount());
		coverageWeight = 1.0/(double)topLocation.getAmbiguityFactor();
				
		for(Pair<Integer,Integer> coordinate : topLocation.getCoordinates()) {
			if(coordinate.getFirst() <= 0)
				break;
			
			for(int i = coordinate.getFirst(); i <= coordinate.getSecond(); i++) {
				coverage[i - coverageOffsetIndex].add(coverageWeight);
			}
		}
	}
	
	private void printLocationToInternalFormat(long filePointer, double readScore,double locationScore, BufferedRandomAccessFile br, StringBuilder tmpLineBuilder,UnsynchronizedBufferedWriter pw) {
		try {
			//print the best hit
			br.seek(filePointer);
			String locationInfo = br.getNextLine();
			StringTokenizer st = new StringTokenizer(locationInfo,"\t");
			//context id
			st.nextToken();
			//read id
			String readId = st.nextToken();
			String chr = st.nextToken();
			String strand = st.nextToken();
		
			ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
			Pattern commaPattern = Pattern.compile(",");
			String[] splittedCoordinates = commaPattern.split(st.nextToken());
			for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
				coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
			}
			
			String mismatches = st.nextToken();
			String hasSpliceSignal = st.nextToken();
			String overlapsKnownJunction = st.nextToken();
			
			
			int overallMappingCount = Integer.valueOf(st.nextToken());
			int overallValidPairsCount = Integer.valueOf(st.nextToken());
			char overlapsPolyA = st.nextToken().charAt(0);
			
			
			tmpLineBuilder.setLength(0);
			tmpLineBuilder.append("global\t").append(readId).append("\t").append(chr).append("\t").append(strand);
			tmpLineBuilder.append("\t").append(coordinates.get(0).getFirst()).append(",").append(coordinates.get(0).getSecond());
			for(int i = 1; i < coordinates.size(); i++) {
				tmpLineBuilder.append(",").append(coordinates.get(i).getFirst()).append(",").append(coordinates.get(i).getSecond());
			}
			tmpLineBuilder.append("\t").append(mismatches).append("\t").append(hasSpliceSignal).append("\t").append(overlapsKnownJunction).append("\t").append(trimDouble(readScore)).append("\t").append(trimDouble(locationScore)).append("\t").append(overallMappingCount).append("\t").append(overallValidPairsCount).append("\t").append(overlapsPolyA);
			
			
			pw.write(tmpLineBuilder.toString());
			pw.newLine();
			
			
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private class LocationContainer {
		
		private long filePointerToBestLocation;
		private long filePointerToSecondBestLocation;
		private double readScore;
		private double locationScore;
		
		public LocationContainer(long filePointer, long filePointerSecondBestLocation, double readScore, double locationScore) {
			this.filePointerToBestLocation = filePointer;
			this.filePointerToSecondBestLocation = filePointerSecondBestLocation;
			this.readScore = readScore;
			this.locationScore = locationScore;
		}
		
		public long getFilePointerToBestLocation() {
			return this.filePointerToBestLocation;
		}
		
		public long getFilePointerToSecondBestLocation() {
			return this.filePointerToSecondBestLocation;
		}
		
		public double getReadScore() {
			return this.readScore;
		}
		
		public double getLocationScore() {
			return this.locationScore;
		}
	}
	
	private class LocationContainerComparator implements Comparator<LocationContainer> {

		@Override
		public int compare(LocationContainer o1, LocationContainer o2) {
			return Long.valueOf(o1.getFilePointerToBestLocation()).compareTo(o2.getFilePointerToBestLocation());
		}
		
	}
	
	
	private String trimDouble(double inValue){
		DecimalFormat twoDec = new DecimalFormat("0.00000",new DecimalFormatSymbols(Locale.US));
		twoDec.setGroupingUsed(false);
		return twoDec.format(inValue);
		}
	
	
	public void setReadScore(ReadPair<SparseRead,SparseRead> readPair) {
		double currentDifference;
		ArrayList<Pair<Integer,Integer>> validPairs;
		double scoreA;
		double scoreB;
		if(readPair.getSecond() != null && !readPair.getValidPairs().isEmpty()) {
			validPairs = readPair.getValidPairs();
			for(Pair<Integer,Integer> validPair : validPairs) {
				scoreA = readPair.getFirst().getLocations().get(validPair.getFirst()).getScore();
				scoreB = 0.0;
				if(readPair.getSecond() != null)
					scoreB = readPair.getSecond().getLocations().get(validPair.getSecond()).getScore();
				
				validPair.setScore(scoreA + scoreB);
			}
			
			Collections.sort(readPair.getValidPairs(),this.pairScoreComparator);
			readPair.setTopScoringPair(readPair.getValidPairs().get(readPair.getValidPairs().size() - 1));
			currentDifference = readPair.getTopScoringPair().getScore();
			if(readPair.getValidPairs().size() > 1)
				currentDifference -= readPair.getValidPairs().get(readPair.getValidPairs().size() -2).getScore();
			readPair.setScore(Math.abs(currentDifference));
		}
		
		else {
			Collections.sort(readPair.getFirst().getLocations(),this.sparseReadLocationComparator);
			scoreA = readPair.getFirst().getLocations().get(readPair.getFirst().getLocations().size() -1).getScore();
			scoreB = 0;
			if(readPair.getFirst().getLocations().size() > 1)
				scoreB = readPair.getFirst().getLocations().get(readPair.getFirst().getLocations().size() -2).getScore();
			
			currentDifference = scoreA - scoreB;
			readPair.setScore(currentDifference);
			
			//resort the pointers
			if(this.isLargeContext) {
				readPair.getFirst().getLocationPointers().clear();
				for(MultiSparseReadLocation l : readPair.getFirst().getLocations()) {
					readPair.getFirst().getLocationPointers().add(l.getFilePointer());
				}
			}
		}
	}
	
	
	private void updateCoverage(MutableDouble[] coverage, ReadPair<SparseRead,SparseRead> readPair, int coverageOffsetIndex) {
		MultiSparseReadLocation topLocation;
		MultiSparseReadLocation secondBestLocation;
		
		/*
		 * process first item
		 * topScoringPair == null <-> no valid pair found
		 * if no valid pair was found, we generate for each mate a ReadPair object containing only one SparseRead (first)
		 * and treat them as single end reads
		 */
		if(readPair.getTopScoringPair() != null) {
			topLocation = readPair.getFirst().getLocations().get(readPair.getTopScoringPair().getFirst());
		}
		else {
			topLocation = readPair.getFirst().getLocations().get(readPair.getFirst().getLocations().size() - 1);
		}
		secondBestLocation = null;
		if(readPair.getFirst().getLocations().size() > 1) {
			if(readPair.getTopScoringPair() != null && readPair.getValidPairs().size() > 1) {
				secondBestLocation = readPair.getFirst().getLocations().get(readPair.getValidPairs().get(readPair.getValidPairs().size() - 2).getFirst());
			}
			else {
				secondBestLocation = readPair.getFirst().getLocations().get(readPair.getFirst().getLocations().size() - 1);
			}
		}
		
		processReadPairItem(readPair.getFirst(),topLocation,secondBestLocation,coverage,coverageOffsetIndex);
		
		//process second item
		if(readPair.getSecond() != null) {
			
			//TODO check if we need this case here. actually if the read pair has a second element here, we also should have a topscoring pair
			if(readPair.getTopScoringPair() != null) {
				topLocation = readPair.getSecond().getLocations().get(readPair.getTopScoringPair().getSecond());
			}
			else {
				topLocation = readPair.getSecond().getLocations().get(readPair.getSecond().getLocations().size() - 1);
			}
			secondBestLocation = null;
			
			if(readPair.getSecond().getLocations().size() > 1) {
				if(readPair.getTopScoringPair() != null && readPair.getValidPairs().size() > 1) {
					secondBestLocation = readPair.getSecond().getLocations().get(readPair.getValidPairs().get(readPair.getValidPairs().size() - 2).getSecond());
				}
				else {
					secondBestLocation = readPair.getSecond().getLocations().get(readPair.getSecond().getLocations().size() - 1);
				}
			}
			processReadPairItem(readPair.getSecond(),topLocation,secondBestLocation,coverage,coverageOffsetIndex);
		}
	}

	private void determineBestMapping(MutableDouble[] coverage,MutableDouble[] coverageUpdated, int coverageOffsetIndex, ArrayList<ReadPair<SparseRead,SparseRead>> reads, BufferedRandomAccessFile bufferedReader,UnsynchronizedBufferedWriter pw, boolean isUniquelyMapped, boolean printReads) throws Exception {
	
		int currentWindowStart;
		int currentWindowEnd;
		int processedWindowSizes;
		int coordinatesIndexBeforeClipping;
		
		double currentMaxCover;
		double tmpCoverage = Double.MIN_VALUE;
		double currentMaxCoverA = Double.MIN_VALUE;
		double currentMaxCoverB = Double.MIN_VALUE;
		
		Microbe currentMicrobe;
		int microbeGenomeStart = -1;
		int microbeGenomeEnd = -1;
		boolean outOfGenome;
		
		ArrayList<Double> windowCoverages = new ArrayList<Double>();
		ArrayList<SparseRead> currentReads = new ArrayList<SparseRead>();
		HashMap<String,Double> window2coverage = new HashMap<String,Double>();
		String windowKey;
		StringBuilder keyBuilder = new StringBuilder();
		StringBuilder tmpLineBuilder = new StringBuilder();
		
		for(ReadPair<SparseRead,SparseRead> readPair : reads) {
			if(readPair.getFirst() == null && readPair.getSecond() == null)
				continue;
			
			currentReads.clear();
			currentReads.add(readPair.getFirst());
			if(readPair.getSecond() != null)
				currentReads.add(readPair.getSecond());
			
			for(SparseRead read : currentReads) {
				
				if(this.isLargeContext) {
					parseReadLocations(read,bufferedReader);
				}
				
				
				if(this.intervalTree != null) {
					if(!this.intervalTree.getIntervalsSpanning(read.getLocations().get(0).getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).isEmpty())
						currentMicrobe = this.intervalTree.getIntervalsSpanning(read.getLocations().get(0).getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).get(0);
					else
						currentMicrobe = this.intervalTree.getIntervalsSpanning(read.getLocations().get(0).getCoordinates().get(0).getSecond(),new ArrayList<Microbe>()).get(0);
						
					microbeGenomeStart = currentMicrobe.getStart();
					microbeGenomeEnd = currentMicrobe.getStop();
				}
				
				for(MultiSparseReadLocation location : read.getLocations()) {
					if(location.getAmbiguityFactor() == 0)
						continue;
					
					if(window2coverage.size() > 2000000) {
						window2coverage.clear();
						
					}
					
					windowCoverages.clear();
					currentMaxCover = Double.MIN_VALUE;
					coordinatesIndexBeforeClipping = location.getCoordinates().size() - 1;
					for(int i = 0; i < location.getCoordinates().size(); i++) {
						Pair<Integer,Integer> segment = location.getCoordinates().get(i);
						
						if(segment.getFirst() <= 0) {
							coordinatesIndexBeforeClipping = i - 1;
							break;
						}
						
						keyBuilder.setLength(0);
						windowKey = keyBuilder.append(segment.getFirst()).append("_").append(segment.getSecond()).toString();
						if(window2coverage.containsKey(windowKey)) {
							tmpCoverage = window2coverage.get(windowKey);
						}
						else {
							tmpCoverage = getMaximumReadCoverage(coverage,coverageOffsetIndex,segment.getFirst(),segment.getSecond());
							window2coverage.put(windowKey, tmpCoverage);
						}
						
						if(tmpCoverage > currentMaxCover) {
							currentMaxCover = tmpCoverage;
						}
					}
					windowCoverages.add(currentMaxCover);
					
					//now go through the windows
					processedWindowSizes = 0;
					for(int j = 0; j < windowSizes.size(); j++) {
						//check upstream coverage
						 currentWindowStart = location.getCoordinates().get(0).getFirst() - processedWindowSizes - windowSizes.get(j);
						 currentWindowEnd = location.getCoordinates().get(0).getFirst() - 1 - processedWindowSizes;
						 
						 
						 outOfGenome = false;
						 if(this.intervalTree != null) {
							 
							 //check upstream genome start
							 //in case the current window is out of the genome, we try to use the mirrored window from the downstream part of the read
							 //if this window is not valid, we use the previous valid window
							 if(currentWindowStart < microbeGenomeStart && currentWindowEnd <= microbeGenomeStart) {
								 //set the window positions to the mirrored version of the current window
								
								currentWindowStart = location.getCoordinates().get(coordinatesIndexBeforeClipping).getSecond() + 1 + processedWindowSizes;
								currentWindowEnd = location.getCoordinates().get(coordinatesIndexBeforeClipping).getSecond() + processedWindowSizes + windowSizes.get(j);
							
								 //mirrored window not valid use previous window coverage
								 if(currentWindowEnd > microbeGenomeEnd && currentWindowStart >= microbeGenomeEnd) {
									if(windowCoverages.size() == 1)
										currentMaxCoverA = windowCoverages.get(0);
									else
										currentMaxCoverA = windowCoverages.get(windowCoverages.size() - 2);
									
									outOfGenome = true;
								 }
								 //mirrored window overlaps with genome end
								 else if(currentWindowEnd > microbeGenomeEnd) {
									 currentWindowEnd = microbeGenomeEnd;
								 }
							 }
							 
							 //in case the window overlaps with the genome start we set the window start to genome start...
							 else if(currentWindowStart < microbeGenomeStart) {
								 currentWindowStart = microbeGenomeStart;
							 }
						 }
						 
						 if(!outOfGenome) {
							 keyBuilder.setLength(0);
							 windowKey = keyBuilder.append(currentWindowStart).append("_").append(currentWindowEnd).toString();
							 if(window2coverage.containsKey(windowKey)) {
								 currentMaxCoverA = window2coverage.get(windowKey);
							 }
							 else {
								 currentMaxCoverA = getMaximumReadCoverage(coverage,coverageOffsetIndex,currentWindowStart,currentWindowEnd);
								 //currentMaxCoverA = getMeanReadCoverage(coverage,currentWindowStart,currentWindowEnd);
								 window2coverage.put(windowKey, currentMaxCoverA);
							 }
						 }
						
						//check downstream coverage
						currentWindowStart = location.getCoordinates().get(coordinatesIndexBeforeClipping).getSecond()  + 1 + processedWindowSizes;
						currentWindowEnd = location.getCoordinates().get(coordinatesIndexBeforeClipping).getSecond() + processedWindowSizes + windowSizes.get(j);
						
						outOfGenome = false;
						if(this.intervalTree != null) {
							//check downstream genome end
							if(currentWindowEnd > microbeGenomeEnd && currentWindowStart >= microbeGenomeEnd) {
								//set the window positions to the mirrored version of the current window
								currentWindowStart = location.getCoordinates().get(0).getFirst() - processedWindowSizes - windowSizes.get(j);
								currentWindowEnd = location.getCoordinates().get(0).getFirst() - 1 - processedWindowSizes;
								if(currentWindowStart < microbeGenomeStart && currentWindowEnd <= microbeGenomeStart) {
									currentMaxCoverB = windowCoverages.get(windowCoverages.size() - 1);
									outOfGenome = true;
								}
								else if(currentWindowStart < microbeGenomeStart) {
									currentWindowStart = microbeGenomeStart;
								}
							 }
							
							 else if(currentWindowEnd > microbeGenomeEnd) {
								 currentWindowEnd = microbeGenomeEnd;
							 }
						 }
						
						if(!outOfGenome) {
							keyBuilder.setLength(0);
							windowKey = keyBuilder.append(currentWindowStart).append("_").append(currentWindowEnd).toString();
							if(window2coverage.containsKey(windowKey)) {
								currentMaxCoverB = window2coverage.get(windowKey);
							}
							else {
								currentMaxCoverB = getMaximumReadCoverage(coverage,coverageOffsetIndex,currentWindowStart,currentWindowEnd);
								//currentMaxCoverB = getMeanReadCoverage(coverage,currentWindowStart,currentWindowEnd);
								window2coverage.put(windowKey, currentMaxCoverB);
							}
						}
						
						currentMaxCover = Math.max(currentMaxCoverA, currentMaxCoverB);
						//windowCoverages.add(currentMaxCover);
						windowCoverages.add(currentMaxCoverA);
						windowCoverages.add(currentMaxCoverB);
						processedWindowSizes += windowSizes.get(j);
					}
					
					location.setScore(calculateLocationScore(windowCoverages));
					
				}
				
			}
			
			
			if(!printReads) {
				setReadScore(readPair);
				updateCoverage(coverageUpdated,readPair,coverageOffsetIndex);
			}
			
			else {
				MultiSparseReadLocation tmpLocation;
				
				if(isUniquelyMapped) {
					if(readPair.getFirst() != null) {
						if(!readPair.getValidPairs().isEmpty()) {
							tmpLocation = readPair.getFirst().getLocations().get(readPair.getValidPairs().get(0).getFirst());
						}
						else {
							tmpLocation = readPair.getFirst().getLocations().get(0);
						}
						printLocationToInternalFormat(tmpLocation.getFilePointer(),0.0,tmpLocation.getScore(),bufferedReader,tmpLineBuilder,pw);
					}
					
					
					if(readPair.getSecond() != null) {
						if(!readPair.getValidPairs().isEmpty()) {
							tmpLocation = readPair.getSecond().getLocations().get(readPair.getValidPairs().get(0).getSecond());
						}
						else {
							
							tmpLocation = readPair.getSecond().getLocations().get(0);
						}
						
						printLocationToInternalFormat(tmpLocation.getFilePointer(),0.0,tmpLocation.getScore(),bufferedReader,tmpLineBuilder,pw);
					}
				}
				
				else {
					if(readPair.getTopScoringPair() != null) {
						tmpLocation = readPair.getFirst().getLocations().get(readPair.getTopScoringPair().getFirst());
					}
					else {
						tmpLocation = readPair.getFirst().getLocations().get(readPair.getFirst().getLocations().size() - 1);
					}
					printLocationToInternalFormat(tmpLocation.getFilePointer(),0.0,tmpLocation.getScore(),bufferedReader,tmpLineBuilder,pw);
					
					
					if(readPair.getSecond() != null) {
						if(readPair.getTopScoringPair() != null) {
							tmpLocation = readPair.getSecond().getLocations().get(readPair.getTopScoringPair().getSecond());
						}
						else {
							tmpLocation = readPair.getSecond().getLocations().get(readPair.getSecond().getLocations().size() - 1);
						}
						printLocationToInternalFormat(tmpLocation.getFilePointer(),0.0,tmpLocation.getScore(),bufferedReader,tmpLineBuilder,pw);
					}
				}
			}
			
			
			if(this.isLargeContext) {
				for(SparseRead read : currentReads) {
					read.getLocations().clear();
					read.setLocations(null);
				}
			}
		}
	}
	
	
	private void parseReadLocations(SparseRead read, BufferedRandomAccessFile braf) throws Exception {
		String currentLine;
		StringBuilder tmpLineBuilder = new StringBuilder();
		StringTokenizer st;
		char strand;
		char mappingType;
		Pattern commaPattern = Pattern.compile(",");
		String[] splittedCoordinates;
		ArrayList<Pair<Integer,Integer>> coordinates;
		String currentSplitKey;
		HashSet<String> currentSplitKeys = new HashSet<String>();
		MultiSparseReadLocation currentSparseLocation;
		ArrayList<MultiSparseReadLocation> locations = new ArrayList<MultiSparseReadLocation>();
		for(long pointer : read.getLocationPointers()) {
			braf.seek(pointer);
			currentLine = braf.getNextLine();
			st = new StringTokenizer(currentLine,"\t");
			//contextId	readId chr
			st.nextToken();
			st.nextToken();
			st.nextToken();
			strand = st.nextToken().toCharArray()[0];
			splittedCoordinates = commaPattern.split(st.nextToken());
			coordinates = new ArrayList<Pair<Integer,Integer>>(2);
			currentSplitKey = "";
			for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
				coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
				currentSplitKey += String.format("%s_%s_",splittedCoordinates[i],splittedCoordinates[i+1]);
			}
			
			if(coordinates.size() == 1 || coordinates.get(1).getFirst() <= 0)
				mappingType = 'F';
			else
				mappingType = 'S';
			
			
			st.nextToken();
			st.nextToken();
			st.nextToken();
			int overallMappingCount = Integer.valueOf(st.nextToken());
			
			if(mappingType != 'S' || !currentSplitKeys.contains(currentSplitKey)) {
				currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
				currentSparseLocation.setFilePointer(pointer);
				currentSparseLocation.setAmbiguityFactor(Math.max(overallMappingCount,read.getLocationPointers().size()));
				locations.add(currentSparseLocation);
				
				if(mappingType == 'S')
					currentSplitKeys.add(currentSplitKey);
				
			}
		}
		
		read.setLocations(locations);
	}
	
	/*
	 * expects for the read region one coverage score and for every other region two scores (up- and downstream)
	 */
	private double calculateLocationScore(ArrayList<Double> coverages) {
		double score = 0.0;
		int windowCount = 1 + this.windowSizes.size();
		if(coverages.get(0) > 0)
			score += Math.pow(2, windowCount) * Math.log(coverages.get(0) + 1);
		windowCount--;
		for(int i = 1; i < coverages.size() - 1; i+=2) {
			if(coverages.get(i) > 0)
				score += Math.pow(2, windowCount) * Math.log(coverages.get(i) + 1);
			if(coverages.get(i+1) > 0)
				score += Math.pow(2, windowCount) * Math.log(coverages.get(i+1) + 1);
			windowCount--;
		}
		return score;
	}

	
/*	private double calculateLocationScore(ArrayList<Double> coverages) {
		double score = 0.0;
		int windowCount = 1 + this.windowSizes.size();
		for(int i = 0; i < coverages.size(); i++) {
			if(coverages.get(i) > 0)
				score += Math.pow(2, windowCount) * Math.log(coverages.get(i) + 1);
			windowCount--;
		}
		
		return score;
	}
*/
	

	private double getMaximumReadCoverage(MutableDouble[] coverages, int coverageIndexOffset, int windowStart, int windowEnd) {
		double maxCover = 0.0;
		for(int i = windowStart - coverageIndexOffset; i <= windowEnd - coverageIndexOffset; i++) {
			if(coverages[i].doubleValue() > maxCover)
				maxCover = coverages[i].doubleValue();
		}
		return maxCover;
		
	}
	
	
	private double getMeanReadCoverage(TreeMap<Integer,MutableDouble> coverages, int windowStart, int windowEnd) {
		double meanCoverage = 0.0;
		//http://www-01.ibm.com/support/docview.wss?uid=swg1IZ47493 (java problem with submap iterator)
		try {
		NavigableMap<Integer,MutableDouble> subMap = coverages.subMap(windowStart,true, windowEnd,true);
		double sum = 0.0;
		for(MutableDouble coverage : subMap.values()) {
				sum += coverage.doubleValue();
		}
		
		return sum/(double)subMap.size();
		}
		catch(Exception e) {
			//e.printStackTrace();
			return meanCoverage;
		}
	}
	
	
	private Pair<Integer,Integer> getContextRange(String mappingFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(mappingFilePath)));
			String currentLine;
			String[] splittedLine;
			String[] splittedCoordinates;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern commaPattern = Pattern.compile(",");
			int start;
			int end;
			int tmpEnd;
			
			int minStart = Integer.MAX_VALUE;
			int maxEnd = Integer.MIN_VALUE; 
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				splittedCoordinates = commaPattern.split(splittedLine[4]);
				
				start = Integer.valueOf(splittedCoordinates[0]);
				end = Integer.valueOf(splittedCoordinates[1]);
				for(int i = 3; i < splittedCoordinates.length; i+=2) {
					tmpEnd = Integer.valueOf(splittedCoordinates[i]);
					if(tmpEnd > 0) {
						end = tmpEnd;
					}
				}
				
				if(start < minStart)
					minStart = start;
				if(end > maxEnd)
					maxEnd = end;
				
			}
			br.close();
			
			Pair<Integer,Integer> result = new Pair<Integer,Integer>();
			result.setFirst(minStart);
			result.setSecond(maxEnd);
			return result;
		}
		
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
		
	}
	
	
	
	private void getCoverage(MutableDouble[] coverage,int coverageIndexOffset,ArrayList<ReadPair<SparseRead,SparseRead>> uniquelyMappingReads,ArrayList<ReadPair<SparseRead,SparseRead>> multiMappingReads, ArrayList<ReadPair<SparseRead,SparseRead>> multiMappingReadsUpdateNeeded, BufferedRandomAccessFile braf, UnsynchronizedBufferedWriter pw) {
		try {
			multiMappingReads.clear();
			BufferedReader br = new BufferedReader(new FileReader(new File(this.multiMappingFilePath)));
			
			long filePointer = 0;
			long tmpPointer;
			
			String currentLine = br.readLine();
			if(currentLine == null)
				return;
				
			filePointer += currentLine.length() + 1;
			
			int directlyPrintedReads = 0;
			
			
			StringTokenizer st = new StringTokenizer(currentLine,"\t");
			Pattern commaPattern = Pattern.compile(",");
			String[] splittedCoordinates;
			String currentContextId = st.nextToken();
			String prevContextId = currentContextId;
			String readId = st.nextToken();
			String readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
			String prevReadId = readId;
			String prevReadIdPrefix = readIdPrefix;
			String chr = st.nextToken();
			char strand = st.nextToken().toCharArray()[0];
			splittedCoordinates = commaPattern.split(st.nextToken());
			ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>(2);
			HashSet<String> currentSplitKeys = new HashSet<String>();
			String currentSplitKey = "";
			for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
				coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
				currentSplitKey += String.format("%s_%s_",splittedCoordinates[i],splittedCoordinates[i+1]);
			}
			
			char mappingType;
			if(coordinates.size() == 1 ||coordinates.get(1).getFirst() <= 0)
				mappingType = 'F';
			else {
				mappingType = 'S';
				currentSplitKeys.add(currentSplitKey);
			}
			
			int mismatches = Integer.valueOf(st.nextToken());
			char strandOfSpliceSignal = st.nextToken().charAt(0);
			boolean overlapsKnownJunction = Boolean.valueOf(st.nextToken());
			int overallMappingCount = Integer.valueOf(st.nextToken());
			int overallMappingCountA = overallMappingCount;
			int overallMappingCountB = Integer.MIN_VALUE;
			
			int overallValidPairsCount = Integer.valueOf(st.nextToken());
			int prevOverallValidPairCounts = overallValidPairsCount;
			
			
			char overlapsPolyA = st.nextToken().charAt(0);
			char overlapsPolyAFirstMate = overlapsPolyA;
			char overlapsPolyASecondMate = '0';
			
			
			ArrayList<MultiSparseReadLocation> currentSparseReadLocations = new ArrayList<MultiSparseReadLocation>();
			ArrayList<ReadLocation> uniquePair = new ArrayList<ReadLocation>();
			MultiSparseReadLocation tmpSparseLocation;
			MultiSparseReadLocation currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
			currentSparseLocation.setFilePointer(0);
			
			
			currentSparseReadLocations.add(currentSparseLocation);
			
			ReadLocation readLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
			readLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
			readLocation.setOverlapsKnownJunction(overlapsKnownJunction);
			uniquePair.add(readLocation);
			
			SparseRead sparseRead = new SparseRead();
			Pair<String,String> currentReadIds = new Pair<String,String>();
			currentReadIds.setFirst(readId);
			
			if(overallValidPairsCount > 0) 
				sparseRead.setOverallMappingCount(overallValidPairsCount);
			else
				sparseRead.setOverallMappingCount(overallMappingCount);
			
			ReadPair<SparseRead,SparseRead> currentPair = new ReadPair<SparseRead,SparseRead>(readIdPrefix);
			ReadPair<SparseRead,SparseRead> tmpPair;
			currentPair.setFirst(sparseRead);
			
			StringBuilder tmpLineBuilder = new StringBuilder();
			while((currentLine = br.readLine()) != null) {
				
				st = new StringTokenizer(currentLine,"\t");
				currentContextId = st.nextToken();
				readId = st.nextToken();
				readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
				
				chr = st.nextToken();
				strand = st.nextToken().toCharArray()[0];
				splittedCoordinates = commaPattern.split(st.nextToken());
				coordinates = new ArrayList<Pair<Integer,Integer>>(2);
				currentSplitKey = "";
				for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
					coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
					currentSplitKey += String.format("%s_%s_",splittedCoordinates[i],splittedCoordinates[i+1]);
				}
				
				if(coordinates.size() == 1 || coordinates.get(1).getFirst() <= 0)
					mappingType = 'F';
				else
					mappingType = 'S';
				
				//mismatches
				mismatches = Integer.valueOf(st.nextToken());
				//known splice signal
				strandOfSpliceSignal = st.nextToken().charAt(0);
				//overlaps known junction
				overlapsKnownJunction = Boolean.valueOf(st.nextToken());
				overallMappingCount = Integer.valueOf(st.nextToken());
				overallValidPairsCount = Integer.valueOf(st.nextToken());
				overlapsPolyA = st.nextToken().charAt(0);

				if(readId.equals(prevReadId)) {
					if(mappingType != 'S' || !currentSplitKeys.contains(currentSplitKey)) {
						currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
						currentSparseLocation.setFilePointer(filePointer);
						currentSparseReadLocations.add(currentSparseLocation);
						
						if(mappingType == 'S')
							currentSplitKeys.add(currentSplitKey);
					}
				}
				
				//second read from the actual fragment
				else if(readIdPrefix.equals(prevReadIdPrefix)) {
					
					overallMappingCountB = overallMappingCount;
					overlapsPolyASecondMate = overlapsPolyA;
					
					//the actual sparse read is readA of the current pair
					for(int i = 0; i < currentSparseReadLocations.size(); i++) {
						//currentSparseReadLocations.get(i).setAmbiguityFactor(currentSparseReadLocations.size());
						currentSparseReadLocations.get(i).setAmbiguityFactor(Math.max(overallMappingCountA,currentSparseReadLocations.size()));
						sparseRead.addLocation(currentSparseReadLocations.get(i));
					}
					updateCoverages(coverage,coverageIndexOffset,currentSparseReadLocations);
					
					//now process the second mate
					sparseRead = new SparseRead();
					currentReadIds.setSecond(readId);
					
					if(overallValidPairsCount > 0) 
						sparseRead.setOverallMappingCount(overallValidPairsCount);
					else
						sparseRead.setOverallMappingCount(overallMappingCount);
					
					currentPair.setSecond(sparseRead);
					currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
					currentSparseLocation.setFilePointer(filePointer);
					currentSparseReadLocations.clear();
					currentSparseReadLocations.add(currentSparseLocation);
					
					readLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
					readLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
					readLocation.setOverlapsKnownJunction(overlapsKnownJunction);
					uniquePair.add(readLocation);
					
					prevReadId = readId;
					currentSplitKeys.clear();
					if(mappingType == 'S')
						currentSplitKeys.add(currentSplitKey);
				}
				
				else {
					
					for(int i = 0; i < currentSparseReadLocations.size(); i++) {
						//currentSparseReadLocations.get(i).setAmbiguityFactor(currentSparseReadLocations.size());
						if(currentPair.getSecond() != null)
							currentSparseReadLocations.get(i).setAmbiguityFactor(Math.max(overallMappingCountB,currentSparseReadLocations.size()));
						else
							currentSparseReadLocations.get(i).setAmbiguityFactor(Math.max(overallMappingCountA,currentSparseReadLocations.size()));
						
						sparseRead.addLocation(currentSparseReadLocations.get(i));
					}
					//update coverages
					updateCoverages(coverage,coverageIndexOffset,currentSparseReadLocations);
					
					//both mates found. check for valid pairs
					if(currentPair.getSecond() != null) {
						getValidReadPairs(currentPair);
						
						//no valid pair found (possibly due to overlapping contexts). -> treat mates as single end reads (if this happens again when considering the global context, we will discard both mates)
						if(currentPair.getValidPairs().isEmpty() || currentPair.getValidPairs().size() > this.maxAllowedValidPairs) {
													
							currentPair.clearValidPairs();
							tmpPair = new ReadPair<SparseRead,SparseRead>(currentReadIds.getSecond());
							tmpPair.setFirst(currentPair.getSecond());
							currentPair.setSecond(null);
							
							if(currentPair.getFirst().getLocations().size() == 1) {
								
								//here we know that this read also appears in another context or it is part of a valid pair (but not in this context). Thus, we have to calculate the final support score later
								if(overallMappingCountA > 1 || prevOverallValidPairCounts > 1) {
									if(this.isLargeContext)
										switchToPointers(currentPair.getFirst());
									this.uniquelyMappingReads.add(currentPair);
								}
								
								else {
									
									directlyPrintedReads++;
									
									tmpLineBuilder.setLength(0);
									//tmpLineBuilder.append("global\t").append(currentPair.getFirst().getId()).append("\t").append(uniquePair.get(0).getChr()).append("\t").append(uniquePair.get(0).getStrand()); 
									tmpLineBuilder.append("global\t").append(currentReadIds.getFirst()).append("\t").append(uniquePair.get(0).getChr()).append("\t").append(uniquePair.get(0).getStrand());
									tmpLineBuilder.append("\t").append(uniquePair.get(0).getCoordinates().get(0).getFirst()).append(",").append(uniquePair.get(0).getCoordinates().get(0).getSecond());
									for(int i = 1; i < uniquePair.get(0).getCoordinates().size(); i++) {
										tmpLineBuilder.append(",").append(uniquePair.get(0).getCoordinates().get(i).getFirst()).append(",").append(uniquePair.get(0).getCoordinates().get(i).getSecond());
									}
									tmpLineBuilder.append("\t").append(uniquePair.get(0).getMismatches()).append("\t").append(uniquePair.get(0).getStrandOfSpliceSignal()).append("\t").append(uniquePair.get(0).overlapsKnownJunction()).append("\t0.0\t0.0\t1\t1\t").append(overlapsPolyAFirstMate);
									pw.write(tmpLineBuilder.toString());
									pw.newLine();
								}
							}
							else {
								if(this.isLargeContext)
									switchToPointers(currentPair.getFirst());
								this.multiMappingReads.add(currentPair);
							}
							
							if(tmpPair.getFirst().getLocations().size() == 1) { 
								
								if(overallMappingCountB > 1 || prevOverallValidPairCounts > 1) {
									if(this.isLargeContext)
										switchToPointers(tmpPair.getFirst());
									this.uniquelyMappingReads.add(tmpPair);
								}
								
								else {
									
									directlyPrintedReads++;
									
									tmpLineBuilder.setLength(0);
									//tmpLineBuilder.append("global\t").append(tmpPair.getFirst().getId()).append("\t").append(uniquePair.get(1).getChr()).append("\t").append(uniquePair.get(1).getStrand());
									tmpLineBuilder.append("global\t").append(currentReadIds.getSecond()).append("\t").append(uniquePair.get(1).getChr()).append("\t").append(uniquePair.get(1).getStrand());
									tmpLineBuilder.append("\t").append(uniquePair.get(1).getCoordinates().get(0).getFirst()).append(",").append(uniquePair.get(1).getCoordinates().get(0).getSecond());
									for(int i = 1; i < uniquePair.get(1).getCoordinates().size(); i++) {
										tmpLineBuilder.append(",").append(uniquePair.get(1).getCoordinates().get(i).getFirst()).append(",").append(uniquePair.get(1).getCoordinates().get(i).getSecond());
									}
									tmpLineBuilder.append("\t").append(uniquePair.get(1).getMismatches()).append("\t").append(uniquePair.get(1).getStrandOfSpliceSignal()).append("\t").append(uniquePair.get(1).overlapsKnownJunction()).append("\t0.0\t0.0\t1\t1\t").append(overlapsPolyASecondMate);
									pw.write(tmpLineBuilder.toString());
									pw.newLine();
									
								}
							}
							else {
								if(this.isLargeContext)
									switchToPointers(tmpPair.getFirst());
								this.multiMappingReads.add(tmpPair);

							}
							
						}
						
						
						//only found one valid pair. no need for further calculations, print mates.
						else if(currentPair.getValidPairs().size() == 1) {
							if(prevOverallValidPairCounts > 1) {
								if(this.isLargeContext) {
									switchToPointers(currentPair.getFirst());
									switchToPointers(currentPair.getSecond());
								}
								this.uniquelyMappingReads.add(currentPair);
							}
							
							
							else {
								directlyPrintedReads += 2;
								printLocationToInternalFormat(currentPair.getFirst().getLocations().get(currentPair.getValidPairs().get(0).getFirst()).getFilePointer(), -1,-1, braf, tmpLineBuilder, pw);
								printLocationToInternalFormat(currentPair.getSecond().getLocations().get(currentPair.getValidPairs().get(0).getSecond()).getFilePointer(), -1,-1, braf, tmpLineBuilder, pw);
							}
							
							
							
						}
						
						//at least two valid pairs found. add mates to the multi mapping fraction 
						else {
							if(this.isLargeContext) {
								switchToPointers(currentPair.getFirst());
								switchToPointers(currentPair.getSecond());
							}
							this.multiMappingReads.add(currentPair);							
						}
					}

					//found only one mate, but with several mapping positions.
					else if(currentPair.getFirst().getLocations().size() > 1) {
						if(this.isLargeContext)
							switchToPointers(currentPair.getFirst());
						this.multiMappingReads.add(currentPair);

					}
					
					else {
						if(overallMappingCountA > 1 || prevOverallValidPairCounts >= 1) {
							if(this.isLargeContext)
								switchToPointers(currentPair.getFirst());
							this.uniquelyMappingReads.add(currentPair);
						}
						else {
							directlyPrintedReads++;
							printLocationToInternalFormat(currentPair.getFirst().getLocations().get(0).getFilePointer(), -1,-1, braf, tmpLineBuilder, pw);
						}
					}
					
					
					overallMappingCountA = overallMappingCount;
					overallMappingCountB = Integer.MIN_VALUE;
					prevOverallValidPairCounts = overallValidPairsCount;
					
					overlapsPolyAFirstMate = overlapsPolyA;
					overlapsPolyASecondMate = '0';
					
					prevReadId = readId;
					prevReadIdPrefix = readIdPrefix;
					sparseRead = new SparseRead();
					currentReadIds.setFirst(readId);
					
					if(overallValidPairsCount > 0) 
						sparseRead.setOverallMappingCount(overallValidPairsCount);
					else
						sparseRead.setOverallMappingCount(overallMappingCount);
					currentPair = new ReadPair<SparseRead,SparseRead>(readIdPrefix);
					currentPair.setFirst(sparseRead);
					currentSparseReadLocations.clear();
					currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
					currentSparseLocation.setFilePointer(filePointer);
					currentSparseReadLocations.add(currentSparseLocation);
					
					readLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
					readLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
					readLocation.setOverlapsKnownJunction(overlapsKnownJunction);
					
					uniquePair.clear();
					uniquePair.add(readLocation);
					
					currentSplitKeys.clear();
					if(mappingType == 'S')
						currentSplitKeys.add(currentSplitKey);
					
				}
				filePointer += currentLine.length() + 1;
				
			}
			
			
			
			//process last bunch...
			for(int i = 0; i < currentSparseReadLocations.size(); i++) {
				//currentSparseReadLocations.get(i).setAmbiguityFactor(currentSparseReadLocations.size());
				if(currentPair.getSecond() != null)
					currentSparseReadLocations.get(i).setAmbiguityFactor(Math.max(overallMappingCountB,currentSparseReadLocations.size()));
				else
					currentSparseReadLocations.get(i).setAmbiguityFactor(Math.max(overallMappingCountA,currentSparseReadLocations.size()));
				
				
				sparseRead.addLocation(currentSparseReadLocations.get(i));
			}
			//update coverages
			updateCoverages(coverage,coverageIndexOffset,currentSparseReadLocations);
			
			
			//both mates found. check for valid pairs
			if(currentPair.getSecond() != null) {

				getValidReadPairs(currentPair);
				
				
				//no valid pair found. -> treat mates as single end reads
				if(currentPair.getValidPairs().isEmpty()) {
					tmpPair = new ReadPair<SparseRead,SparseRead>(currentReadIds.getSecond());
					tmpPair.setFirst(currentPair.getSecond());
					currentPair.setSecond(null);
					
					if(currentPair.getFirst().getLocations().size() == 1) {
						if(overallMappingCountA > 1 || prevOverallValidPairCounts > 1) {
							if(this.isLargeContext)
								switchToPointers(currentPair.getFirst());
							this.uniquelyMappingReads.add(currentPair);
						}
						
						else {
							
							directlyPrintedReads++;
							
							tmpLineBuilder.setLength(0);
							//tmpLineBuilder.append("global\t").append(currentPair.getFirst().getId()).append("\t").append(uniquePair.get(0).getChr()).append("\t").append(uniquePair.get(0).getStrand());
							tmpLineBuilder.append("global\t").append(currentReadIds.getFirst()).append("\t").append(uniquePair.get(0).getChr()).append("\t").append(uniquePair.get(0).getStrand());
							tmpLineBuilder.append("\t").append(coordinates.get(0).getFirst()).append(",").append(coordinates.get(0).getSecond());
							for(int i = 1; i < coordinates.size(); i++) {
								tmpLineBuilder.append(",").append(coordinates.get(i).getFirst()).append(",").append(coordinates.get(i).getSecond());
							}
							tmpLineBuilder.append("\t").append(uniquePair.get(0).getMismatches()).append("\t").append(uniquePair.get(0).getStrandOfSpliceSignal()).append("\t").append(uniquePair.get(0).overlapsKnownJunction()).append("\t0.0\t0.0\t1\t1\t").append(overlapsPolyAFirstMate);
							pw.write(tmpLineBuilder.toString());
							pw.newLine();
						}
						
					}
					else {
						if(this.isLargeContext)
							switchToPointers(currentPair.getFirst());
						this.multiMappingReads.add(currentPair);
					}
					
					if(tmpPair.getFirst().getLocations().size() == 1) {
						if(overallMappingCountB > 1 || prevOverallValidPairCounts > 1) {
							if(this.isLargeContext)
								switchToPointers(tmpPair.getFirst());
							this.uniquelyMappingReads.add(tmpPair);
						}
						
						else{
							
							directlyPrintedReads++;
							
							tmpLineBuilder.setLength(0);
							//tmpLineBuilder.append("global\t").append(tmpPair.getFirst().getId()).append("\t").append(uniquePair.get(1).getChr()).append("\t").append(uniquePair.get(1).getStrand());
							tmpLineBuilder.append("global\t").append(currentReadIds.getSecond()).append("\t").append(uniquePair.get(1).getChr()).append("\t").append(uniquePair.get(1).getStrand());
							tmpLineBuilder.append("\t").append(coordinates.get(0).getFirst()).append(",").append(coordinates.get(0).getSecond());
							for(int i = 1; i < coordinates.size(); i++) {
								tmpLineBuilder.append(",").append(coordinates.get(i).getFirst()).append(",").append(coordinates.get(i).getSecond());
							}
							tmpLineBuilder.append("\t").append(uniquePair.get(1).getMismatches()).append("\t").append(uniquePair.get(1).getStrandOfSpliceSignal()).append("\t").append(uniquePair.get(1).overlapsKnownJunction()).append("\t0.0\t0.0\t1\t1\t").append(overlapsPolyASecondMate);
							pw.write(tmpLineBuilder.toString());
							pw.newLine();
						}
					}
					else {
						if(this.isLargeContext)
							switchToPointers(tmpPair.getFirst());
						this.multiMappingReads.add(tmpPair);

					}
					
				}
				
				//only found one valid pair. no need for further calculations, print mates.
				else if(currentPair.getValidPairs().size() == 1) {
					if(prevOverallValidPairCounts > 1) {
						if(this.isLargeContext) {
							switchToPointers(currentPair.getFirst());
							switchToPointers(currentPair.getSecond());
						}
						this.uniquelyMappingReads.add(currentPair);
					}
					
					else {
						
						directlyPrintedReads += 2;
						printLocationToInternalFormat(currentPair.getFirst().getLocations().get(currentPair.getValidPairs().get(0).getFirst()).getFilePointer(), -1,-1, braf, tmpLineBuilder, pw);
						printLocationToInternalFormat(currentPair.getSecond().getLocations().get(currentPair.getValidPairs().get(0).getSecond()).getFilePointer(), -1,-1, braf, tmpLineBuilder, pw);
					}
					
				}
				
				//at least two valid pairs found. add mates to the multi mapping fraction 
				else {
					if(this.isLargeContext) {
						switchToPointers(currentPair.getFirst());
						switchToPointers(currentPair.getSecond());
					}
					this.multiMappingReads.add(currentPair);
				}
				
			}

			//found only one mate, but with several mapping positions.
			else if(currentPair.getFirst().getLocations().size() > 1) {
				if(this.isLargeContext)
					switchToPointers(currentPair.getFirst());
				this.multiMappingReads.add(currentPair);
			}
			
			else {
				if(overallMappingCountA > 1 || prevOverallValidPairCounts >= 1) {
					if(this.isLargeContext)
						switchToPointers(currentPair.getFirst());
					this.uniquelyMappingReads.add(currentPair);
				}
				else {
					directlyPrintedReads++;
					printLocationToInternalFormat(currentPair.getFirst().getLocations().get(0).getFilePointer(), -1,-1, braf, tmpLineBuilder, pw);
				}
			}
			
			
			
			//finally add coverages up and downstream of the context start and end, respectively.
			for(int i : this.upstreamCoverage.keySet()) {
				coverage[i - coverageIndexOffset].setValue(this.upstreamCoverage.get(i));
			}
			for(int i : this.downstreamCoverage.keySet()) {
				coverage[i - coverageIndexOffset].setValue(this.downstreamCoverage.get(i));
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	private void switchToPointers(SparseRead read) {
		ArrayList<Long> locationPointers = new ArrayList<Long>();
		for(MultiSparseReadLocation location : read.getLocations()) {
			locationPointers.add(location.getFilePointer());
		}
		
		read.getLocations().clear();
		read.setLocations(null);
		read.setLocationPointers(locationPointers);
	}
	
	
	/**
	 * currently this function only uses the following constraints:
	 * 
	 * - mates have to be aligned to opposite strands (first: + -> second: - or first: - -> second: +)
	 * - the distance between two mates is not allowed to be larger than the maximum context size (see variable maxContextSize).
	 *
	 * 
	 * TODO: calculate a more accurate fragment length distribution from fully aligned reads to further restrict the distance criterion.
	 * TODO: Think about a maximum allowed context-score difference between two mates. It makes no sense to considers cases where one mate
	 * has almost no support by other reads whereas the other one is seems to originate from a highly expressed region.
	 * 
	 * 
	 *
	 * 
	 **/
	
	public void getValidReadPairs(ReadPair<SparseRead,SparseRead> pair) {
		try {
			ArrayList<MultiSparseReadLocation> locationsA;
			ArrayList<MultiSparseReadLocation> locationsB;
			locationsA = pair.getFirst().getLocations();
			int currentStart;
			int currentEnd;
			int currentDistance;
			
			int secondLocationEnd;
			
			int coordinatesIndexBeforeClippingA;
			int coordinatesIndexBeforeClippingB;
			
			if(pair.getSecond() != null) {
				locationsB = pair.getSecond().getLocations();
				
				for(int i = 0; i < locationsA.size(); i++) {
					
					coordinatesIndexBeforeClippingA = locationsA.get(i).getCoordinates().size() - 1;
					for(int k = 0; k < locationsA.get(i).getCoordinates().size(); k++) {
						Pair<Integer,Integer> segment = locationsA.get(i).getCoordinates().get(k);
						
						if(segment.getFirst() <= 0) {
							coordinatesIndexBeforeClippingA = k - 1;
							break;
						}
					}
					
					for(int j = 0; j < locationsB.size(); j++) {
						
						coordinatesIndexBeforeClippingB = locationsB.get(j).getCoordinates().size() - 1;
						for(int k = 0; k < locationsB.get(j).getCoordinates().size(); k++) {
							Pair<Integer,Integer> segment = locationsB.get(j).getCoordinates().get(k);
							
							if(segment.getFirst() <= 0) {
								coordinatesIndexBeforeClippingB = k - 1;
								break;
							}
						}
						
						//strand criterion fulfilled
						if(locationsA.get(i).getStrand() != locationsB.get(j).getStrand()) {
							
							//A starts upstream of B
							if(locationsA.get(i).getCoordinates().get(0).getFirst() <= locationsB.get(j).getCoordinates().get(0).getFirst()) {
								
								
								currentStart = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getFirst();
								currentEnd = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getSecond();
								
/*								if(locationsA.get(i).getCoordinates().size() > 1 && locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size()-1).getSecond() > 0) {
									currentStart = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size()-1).getFirst();
									currentEnd = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size()-1).getSecond();
								}
*/								
								currentDistance = locationsB.get(j).getCoordinates().get(0).getFirst() - currentEnd;
								
								//distance criterion fulfilled
								if(currentDistance >= 0 && currentDistance <= this.maxContextSize) {
									pair.addValidPair(new Pair<Integer,Integer>(i,j));
								}
								else if(currentDistance < 0) {
									secondLocationEnd = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getSecond();
/*									if(locationsB.get(j).getCoordinates().size() > 1 && locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size() - 1).getSecond() > 0)
										secondLocationEnd = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size() - 1).getSecond();
*/									
									if(currentEnd < secondLocationEnd && secondLocationEnd - currentEnd <= this.maxContextSize) {
										pair.addValidPair(new Pair<Integer,Integer>(i,j));
									}
								}
							}
							
							//B starts upstream of A
							else if(locationsA.get(i).getCoordinates().get(0).getFirst() > locationsB.get(j).getCoordinates().get(0).getFirst()) {
								
								currentStart = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getFirst();
								currentEnd = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getSecond();
/*								if(locationsB.get(j).getCoordinates().size() > 1 && locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size()-1).getSecond() > 0) {
									currentStart = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size()-1).getFirst();
									currentEnd = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size()-1).getSecond();
								}
*/								
								currentDistance = locationsA.get(i).getCoordinates().get(0).getFirst() - currentEnd;
								
								if(currentDistance >= 0 && currentDistance <= this.maxContextSize) {
									pair.addValidPair(new Pair<Integer,Integer>(i,j));
								}
								else if(currentDistance < 0) {
									secondLocationEnd = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getSecond();
/*									if(locationsA.get(i).getCoordinates().size() > 1 && locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size() - 1).getSecond() > 0)
										secondLocationEnd = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size() - 1).getSecond();
*/									
									if(currentEnd < secondLocationEnd && secondLocationEnd - currentEnd <= this.maxContextSize) {
										pair.addValidPair(new Pair<Integer,Integer>(i,j));
									}
								}
							}
						}
					}
					
					if(pair.getValidPairs().size() > this.maxAllowedValidPairs)
						return;
				}
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	
	private void updateCoverages(MutableDouble[] coverage,int coverageIndexOffset,ArrayList<MultiSparseReadLocation> currentSparseReadLocations) {
		double coverageWeight = 1.0/(double)currentSparseReadLocations.get(0).getAmbiguityFactor();
		
		for(MultiSparseReadLocation location : currentSparseReadLocations) {
			for(Pair<Integer,Integer> segment : location.getCoordinates()) {
				if(segment.getFirst() <= 0)
					break;
				
				for(int i = segment.getFirst(); i <= segment.getSecond(); i++) {
					coverage[i - coverageIndexOffset].add(coverageWeight);
				}
			}
		}
	}

	

	
	private class SparseLocationScoreComparator implements Comparator<MultiSparseReadLocation> {
		public int compare(MultiSparseReadLocation l1, MultiSparseReadLocation l2) {
			double score1 = l1.getScore();
			double score2 = l2.getScore();
			return ((Double)score1).compareTo(score2);
		}
		
	}
	
	private class PairScoreComparator implements Comparator<Pair> {
		public int compare(Pair r1, Pair r2) {
			double score1 = r1.getScore();
			double score2 = r2.getScore();
			return ((Double)score1).compareTo(score2);
		}
		
	}
	
}
