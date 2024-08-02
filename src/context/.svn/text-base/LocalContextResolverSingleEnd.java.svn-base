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

public class LocalContextResolverSingleEnd implements ContextResolver {
	
	private String multiMappingFilePath;
	private String outputPath;
	
	private SortedMap<Integer,MutableDouble> upstreamCoverage;
	private SortedMap<Integer,MutableDouble> downstreamCoverage;
	private IntervalTree<Microbe> intervalTree;
	private ArrayList<Integer> windowSizes;
	private int completeWindowIntervall;
	private int maxContextSize;
	private SparseLocationScoreComparator sparseLocationScoreComparator;
	
	private boolean updateQueue;
	private boolean developer;
	private int updateInterval;
	private boolean verbose;

	private boolean isLargeContext;
	
	public LocalContextResolverSingleEnd(String multiMappingFilePath, String outputPath,TreeMap<Integer,MutableDouble> upstreamCoverage, TreeMap<Integer,MutableDouble> downstreamCoverage, IntervalTree<Microbe> intervalTree, ArrayList<Integer> windowSizes, int readLength, int maxContextSize,boolean updateQueue, int updateInterval, boolean verbose, boolean developer, boolean isLargeContext) {
		this.multiMappingFilePath = multiMappingFilePath;
		this.windowSizes = windowSizes;
		this.completeWindowIntervall = 0;
		for(int i = 0; i < windowSizes.size(); i++) {
			this.completeWindowIntervall += windowSizes.get(i);
		}
		this.outputPath = outputPath;
		this.maxContextSize = maxContextSize;
		
		this.sparseLocationScoreComparator = new SparseLocationScoreComparator();
		
		
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
			
			
			ArrayList<SparseRead> uniquelyMappingReads = new ArrayList<SparseRead>();
			ArrayList<SparseRead> multiMappingReads = new ArrayList<SparseRead>();
			ArrayList<SparseRead> multiMappingReadsUpdateNeeded = new ArrayList<SparseRead>();
			
			
			Pair<Integer,Integer> contextRange = getContextRange(this.multiMappingFilePath);
			int coverageOffsetIndex = contextRange.getFirst() - this.completeWindowIntervall;
			this.upstreamCoverage = this.upstreamCoverage.subMap(contextRange.getFirst() - this.completeWindowIntervall, contextRange.getFirst());
			this.downstreamCoverage = this.downstreamCoverage.subMap(contextRange.getSecond() + 1, contextRange.getSecond() + 1 +  this.completeWindowIntervall);
			
			MutableDouble[] coverage = new MutableDouble[(contextRange.getSecond() - contextRange.getFirst() + 1) + (2* this.completeWindowIntervall)];
			for(int i = 0; i < coverage.length; i++) {
				coverage[i] = new MutableDouble();
			}
			getCoverage(coverage,coverageOffsetIndex,uniquelyMappingReads,multiMappingReads,multiMappingReadsUpdateNeeded,bufferedReader,pw);

			
			//we determine for every read the best location (in terms of support score) and update the coverage values by discarding all other alignments of that read
			MutableDouble[] coverageUpdated = new MutableDouble[coverage.length];
			for(int i = 0; i < coverage.length; i++) {
				coverageUpdated[i] = new MutableDouble(coverage[i].doubleValue());
			}
			getLocationScores(coverage,coverageUpdated,coverageOffsetIndex,multiMappingReads,bufferedReader,pw,false,false);
			getLocationScores(coverage,null,coverageOffsetIndex,uniquelyMappingReads,bufferedReader,pw,true,true);
			getLocationScores(coverage,null,coverageOffsetIndex,multiMappingReads,bufferedReader,pw,false,true);
			
			bufferedReader.close();
			pw.flush();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
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
	
	
	private void setTopScoringLocation(SparseRead read) {
		Collections.sort(read.getLocations(),this.sparseLocationScoreComparator);
		read.setTopScoringLocation(read.getLocations().get(read.getLocations().size() - 1));
		
		//re-sort file pointers
		if(this.isLargeContext) {
			read.getLocationPointers().clear();
			for(MultiSparseReadLocation l : read.getLocations()) {
				read.getLocationPointers().add(l.getFilePointer());
			}
		}
		
	}
	
	private void getLocationScores(MutableDouble[] coverage,MutableDouble[] coverageUpdated, int coverageOffsetIndex, ArrayList<SparseRead> reads, BufferedRandomAccessFile bufferedReader,UnsynchronizedBufferedWriter pw, boolean isUniquelyMapped, boolean printReads) throws Exception {
		
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
		HashMap<String,Double> window2coverage = new HashMap<String,Double>();
		String windowKey;
		StringBuilder keyBuilder = new StringBuilder();
		StringBuilder tmpLineBuilder = new StringBuilder();
		for(SparseRead read : reads) {
			if(read == null)
				continue;
			
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
				
				if(window2coverage.size() > 1000000)
					window2coverage.clear();
				
				windowCoverages.clear();
				currentMaxCover = Double.MIN_VALUE;			
				coordinatesIndexBeforeClipping = location.getCoordinates().size() - 1;
				//the first window consists of the read coordinates itself
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
							 window2coverage.put(windowKey, currentMaxCoverA);
						 }
					 }
					
					//check downstream coverage
					currentWindowStart = location.getCoordinates().get(coordinatesIndexBeforeClipping).getSecond() + 1 + processedWindowSizes;
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
			
			if(!printReads) {
				setTopScoringLocation(read);
				updateCoverage(coverageUpdated,read,coverageOffsetIndex);
			}
			
			else {
				MultiSparseReadLocation tmpLocation;
				if(isUniquelyMapped) {
					tmpLocation = read.getLocations().get(0);
				}
				else {
					tmpLocation = read.getTopScoringLocation();
				}
				printLocationToInternalFormat(tmpLocation.getFilePointer(),0.0,tmpLocation.getScore(),bufferedReader,tmpLineBuilder,pw);
			}
			
			if(this.isLargeContext) {
				read.getLocations().clear();
				read.setLocations(null);
			}
		}
	}
	
	
	private void updateCoverage(MutableDouble[] coverage, SparseRead highScoringRead, int coverageOffsetIndex) {
		double coverageWeight;
		for(MultiSparseReadLocation tmpLocation : highScoringRead.getLocations()) {
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
	
		highScoringRead.getTopScoringLocation().setAmbiguityFactor(highScoringRead.getOverallMappingCount());
		coverageWeight = 1.0/(double)highScoringRead.getTopScoringLocation().getAmbiguityFactor();
		for(Pair<Integer,Integer> coordinate : highScoringRead.getTopScoringLocation().getCoordinates()) {
			if(coordinate.getFirst() <= 0)
				break;
			
			for(int i = coordinate.getFirst(); i <= coordinate.getSecond(); i++) {
				coverage[i - coverageOffsetIndex].add(coverageWeight);
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
			
			
			if(mappingType != 'S' || !currentSplitKeys.contains(currentSplitKey)) {
				currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
				currentSparseLocation.setFilePointer(pointer);
				currentSparseLocation.setAmbiguityFactor(Math.max(read.getOverallMappingCount(),read.getLocationPointers().size()));
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



	
	private class SparseLocationScoreComparator implements Comparator<MultiSparseReadLocation> {
		public int compare(MultiSparseReadLocation l1, MultiSparseReadLocation l2) {
			double score1 = l1.getScore();
			double score2 = l2.getScore();
			return ((Double)score1).compareTo(score2);
		}
		
	}
	
	
	private void getCoverage(MutableDouble[] coverage,int coverageOffsetIndex,ArrayList<SparseRead> uniquelyMappingReads, ArrayList<SparseRead> multiMappingReads,ArrayList<SparseRead> multiMappingReadsUpdateNeeded, BufferedRandomAccessFile br, UnsynchronizedBufferedWriter pw) throws Exception {
		multiMappingReads.clear();
		
		
		long prevFilePointer = br.getFilePointer();
		String currentLine = br.getNextLine();
		if(currentLine == null)
			return;
			
		double coverageWeight;
		long filePointer = br.getFilePointer();
		StringTokenizer st = new StringTokenizer(currentLine,"\t");
		String currentContextId = st.nextToken();
		String prevContextId = currentContextId;
		String readId = st.nextToken();
		String prevReadId = readId;
		String chr = st.nextToken();
		char strand = st.nextToken().toCharArray()[0];
		Pattern commaPattern = Pattern.compile(",");
		String[] splittedCoordinates = commaPattern.split(st.nextToken());
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
		int prevOverallMappingCount = overallMappingCount;
		//skip overall valid pair count
		st.nextToken();
		char overlapsPolyA = st.nextToken().charAt(0);
		char prevOverlapsPolyA = overlapsPolyA;
		
		
		ArrayList<MultiSparseReadLocation> currentSparseReadLocations = new ArrayList<MultiSparseReadLocation>();
		MultiSparseReadLocation tmpSparseLocation;
		MultiSparseReadLocation currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
		currentSparseLocation.setFilePointer(prevFilePointer);
		currentSparseReadLocations.add(currentSparseLocation);
		
		ReadLocation readLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
		readLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
		readLocation.setOverlapsKnownJunction(overlapsKnownJunction);
		
		StringBuilder tmpLineBuilder = new StringBuilder();
		while((currentLine = br.getNextLine()) != null) {
			
			
			String[] splittedLine = currentLine.split("\t");
			if(splittedLine.length < 11) {
				System.err.println("Missing field(s) in line: ");
				System.err.println(currentLine);
				System.exit(1);
			}
			
			
			st = new StringTokenizer(currentLine,"\t");
			currentContextId = st.nextToken();
			readId = st.nextToken();
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
			//skip overall valid pair count
			st.nextToken();
			overlapsPolyA = st.nextToken().charAt(0);
			
			if(readId.equals(prevReadId)) {
				if(mappingType != 'S' || !currentSplitKeys.contains(currentSplitKey)) {
					currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
					currentSparseReadLocations.add(currentSparseLocation);
					currentSparseLocation.setFilePointer(filePointer);
				}
				if(mappingType == 'S')
					currentSplitKeys.add(currentSplitKey);
			}
			
			else {
				
				//update coverages
				coverageWeight = 1.0/(double)Math.max(prevOverallMappingCount,currentSparseReadLocations.size());
				//coverageWeight = 1.0/(double)currentSparseReadLocations.size();
				for(MultiSparseReadLocation location : currentSparseReadLocations) {
					for(Pair<Integer,Integer> segment : location.getCoordinates()) {
						if(segment.getFirst() <= 0)
							break;
						
						for(int i = segment.getFirst(); i <= segment.getSecond(); i++) {
							coverage[i - coverageOffsetIndex].add(coverageWeight);
						}
					}
				}
				
				
				//in case we have more than one read location for the current read we add this read to the actual context
				if(currentSparseReadLocations.size() > 1) {
					SparseRead currentRead = new SparseRead();
					currentRead.setOverallMappingCount(prevOverallMappingCount);
					for(int i = 0; i < currentSparseReadLocations.size(); i++) {
						tmpSparseLocation = currentSparseReadLocations.get(i);
						tmpSparseLocation.setAmbiguityFactor(Math.max(prevOverallMappingCount,currentSparseReadLocations.size()));
						currentRead.addLocation(tmpSparseLocation);
						
					}
					
					if(this.isLargeContext)
						switchToPointers(currentRead);
					
					multiMappingReads.add(currentRead);
				}
				
				else {
					
					if(prevOverallMappingCount > 1) {
						SparseRead currentRead = new SparseRead();
						currentRead.setOverallMappingCount(prevOverallMappingCount);
						currentSparseReadLocations.get(0).setAmbiguityFactor(Math.max(prevOverallMappingCount,currentSparseReadLocations.size()));
						currentRead.addLocation(currentSparseReadLocations.get(0));
						
						if(this.isLargeContext)
							switchToPointers(currentRead);
						uniquelyMappingReads.add(currentRead);
						
					}
					else {
						tmpLineBuilder.setLength(0);
						tmpLineBuilder.append("global\t").append(prevReadId).append("\t").append(readLocation.getChr()).append("\t").append(readLocation.getStrand());
						tmpLineBuilder.append("\t").append(readLocation.getCoordinates().get(0).getFirst()).append(",").append(readLocation.getCoordinates().get(0).getSecond());
						for(int i = 1; i < readLocation.getCoordinates().size(); i++) {
							tmpLineBuilder.append(",").append(readLocation.getCoordinates().get(i).getFirst()).append(",").append(readLocation.getCoordinates().get(i).getSecond());
						}
						
						tmpLineBuilder.append("\t").append(readLocation.getMismatches()).append("\t").append(readLocation.getStrandOfSpliceSignal()).append("\t").append(readLocation.overlapsKnownJunction()).append("\t0.0\t0.0\t1\t1\t").append(prevOverlapsPolyA);
						pw.write(tmpLineBuilder.toString());
						pw.newLine();
					}

				}
				
				prevReadId = readId;
				prevOverallMappingCount = overallMappingCount;
				currentSparseReadLocations.clear();
				currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
				currentSparseLocation.setFilePointer(filePointer);
				currentSparseReadLocations.add(currentSparseLocation);
				
				readLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
				readLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
				readLocation.setOverlapsKnownJunction(overlapsKnownJunction);
				prevOverlapsPolyA = overlapsPolyA;
				
				currentSplitKeys.clear();
				if(mappingType == 'S')
					currentSplitKeys.add(currentSplitKey);
			}
			filePointer = br.getFilePointer();
		}
		
		
		//process last bunch...
		
		//update coverages
		coverageWeight = 1.0/(double)Math.max(prevOverallMappingCount,currentSparseReadLocations.size());
		//coverageWeight = 1.0/(double)currentSparseReadLocations.size();
		for(MultiSparseReadLocation location : currentSparseReadLocations) {
			for(Pair<Integer,Integer> segment : location.getCoordinates()) {

				if(segment.getFirst() <= 0)
					break;
				
				for(int i = segment.getFirst(); i <= segment.getSecond(); i++) {
					coverage[i - coverageOffsetIndex].add(coverageWeight);
				}
			}
			
		}
		
		
		if(currentSparseReadLocations.size() == 1) {
			
			if(prevOverallMappingCount > 1) {
				SparseRead currentRead = new SparseRead();
				currentRead.setOverallMappingCount(prevOverallMappingCount);
				currentSparseReadLocations.get(0).setAmbiguityFactor(Math.max(prevOverallMappingCount,currentSparseReadLocations.size()));
				currentRead.addLocation(currentSparseReadLocations.get(0));
				
				if(this.isLargeContext)
					switchToPointers(currentRead);
				
				uniquelyMappingReads.add(currentRead);
				
			}
			
			else {
				tmpLineBuilder.setLength(0);
				tmpLineBuilder.append("global\t").append(prevReadId).append("\t").append(readLocation.getChr()).append("\t").append(readLocation.getStrand());
				tmpLineBuilder.append("\t").append(readLocation.getCoordinates().get(0).getFirst()).append(",").append(readLocation.getCoordinates().get(0).getSecond());
				for(int i = 1; i < readLocation.getCoordinates().size(); i++) {
					tmpLineBuilder.append(",").append(readLocation.getCoordinates().get(i).getFirst()).append(",").append(readLocation.getCoordinates().get(i).getSecond());
				}
				
				tmpLineBuilder.append("\t").append(readLocation.getMismatches()).append("\t").append(readLocation.getStrandOfSpliceSignal()).append("\t").append(readLocation.overlapsKnownJunction()).append("\t0.0\t0.0\t1\t1\t").append(overlapsPolyA);
				pw.write(tmpLineBuilder.toString());
				pw.newLine();
			}

		}
		
		else {
			SparseRead currentRead = new SparseRead();
			currentRead.setOverallMappingCount(prevOverallMappingCount);
			for(int i = 0; i < currentSparseReadLocations.size(); i++) {
				tmpSparseLocation = currentSparseReadLocations.get(i);
				tmpSparseLocation.setAmbiguityFactor(Math.max(prevOverallMappingCount,currentSparseReadLocations.size()));
				currentRead.addLocation(tmpSparseLocation);
			}
			
			if(this.isLargeContext)
				switchToPointers(currentRead);
			
			multiMappingReads.add(currentRead);
		}
		
		
		//finally add coverages up and downstream of the context start and end, respectively.
		for(int i : this.upstreamCoverage.keySet()) {
			coverage[i - coverageOffsetIndex].setValue(this.upstreamCoverage.get(i));
		}
		for(int i : this.downstreamCoverage.keySet()) {
			coverage[i - coverageOffsetIndex].setValue(this.downstreamCoverage.get(i));
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
	
	
	private class SparseReadComparator implements Comparator {
		public SparseReadComparator() {
			
		}
		@Override
		public int compare(Object o1, Object o2) {
			SparseRead l1 = (SparseRead)o1;
			SparseRead l2 = (SparseRead)o2;
			return(Integer.valueOf(l1.getLocations().get(0).getCoordinates().get(0).getFirst()).compareTo(l2.getLocations().get(0).getCoordinates().get(0).getFirst()));
		}
	}
	
}
