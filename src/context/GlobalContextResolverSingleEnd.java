package context;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
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
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import main.Context;
import main.Microbe;
import main.MutableDouble;
import main.Pair;
import main.Read;
import main.ReadLocation;
import main.SparseRead;
import main.MultiSparseReadLocation;
import main.SparseReadLocation;
import tools.BufferedRandomAccessFile;
import tools.MaxPriorityQueue;

public class GlobalContextResolverSingleEnd implements ContextResolver {
	
	private String multiMappingFilePath;
	private String outputPath;
	private ArrayList<Integer> windowSizes;
	private int completeWindowIntervall;
	private int maxContextSize;
	private int maxDelSize;
	private int numberOfThreads;
	private SparseLocationScoreComparator sparseLocationScoreComparator;
	
	private boolean updateQueue;
	private boolean printMultiMappings;
	private boolean printSecondBestChr;
	private int updateInterval;
	private boolean verbose;
	private boolean strandSpecific;
	
	HashMap<String,Integer> chr2length;
	HashMap<String,IntervalTree<Microbe>> chr2intervalTree;
	
	private double scoreDiffCutoff;
	

	public GlobalContextResolverSingleEnd(String multiMappingFilePath, String outputPath,HashMap<String,Integer> chr2length, HashMap<String,IntervalTree<Microbe>> chr2intervalTree, ArrayList<Integer> windowSizes, int readLength, int maxContextSize, int maxDelSize, boolean localContexts,boolean updateQueue, int updateInterval, int numberOfThreads, boolean printMultiMappings, boolean printSecondBestChr, boolean verbose, boolean strandSpecific) {
		this.multiMappingFilePath = multiMappingFilePath;
		this.windowSizes = windowSizes;
		this.completeWindowIntervall = 0;
		for(int i = 0; i < windowSizes.size() && i < 1; i++) {
			this.completeWindowIntervall += windowSizes.get(i);
		}
		
		this.chr2length = chr2length;
		this.chr2intervalTree = chr2intervalTree;
		
		this.outputPath = outputPath;
		this.maxContextSize = maxContextSize;
		this.maxDelSize = maxDelSize;
		this.sparseLocationScoreComparator = new SparseLocationScoreComparator();
		this.updateQueue = updateQueue;
		this.updateInterval = updateInterval;
		this.numberOfThreads = numberOfThreads;
		this.printMultiMappings = printMultiMappings;
		this.printSecondBestChr = printSecondBestChr;
		this.verbose = verbose;
		this.strandSpecific = strandSpecific;
		
		//this.scoreDiffCutoff = new CoverageCalculator(this.windowSizes).getScoreCutoff();
		this.scoreDiffCutoff = 0;
		
	}
	
	public void resolve() {
		try {			
			HashMap<String,BufferedWriter> chr2writer = new HashMap<String,BufferedWriter>();
			processGlobalContext(chr2writer);
			
			for(BufferedWriter bw : chr2writer.values()) {
				bw.close();
			}
			
			
			//finally write the header and concatenate the written files
			ArrayList<String> chrNames = new ArrayList<String>();
			File tmpFile;
			
			for(String chrName : this.chr2length.keySet()) {
				tmpFile = new File(this.outputPath + "." + chrName);
				if(tmpFile.exists()) {
					if(tmpFile.length() > 0) {
						chrNames.add(chrName);
					}
					else
						tmpFile.delete();
				}
			}
		
			
			Collections.sort(chrNames);
			ArrayList<String> tmpFilePaths = new ArrayList<String>();
			for(String chrName : chrNames) {
				tmpFilePaths.add(this.outputPath + "." + chrName);
			}
			

			writeSamHeader(chrNames,this.chr2length,this.outputPath);
			concatenateFilesWithNIO(tmpFilePaths, this.outputPath);
			
			for(String tmpPath : tmpFilePaths) {
				new File(tmpPath).delete();
			}

		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private void writeSamHeader(ArrayList<String> chrNames, HashMap<String,Integer> chr2length, String outputPath) throws Exception {
		PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
		for(String chrName : chrNames) {
			pw.println(String.format("@SQ\tSN:%s\tLN:%s",chrName,chr2length.get(chrName)));
		}
		pw.close();
		
	}
		
	private void concatenateFilesWithNIO(ArrayList<String> filePaths, String outputPath) {
		try {
			FileOutputStream fos = new FileOutputStream(outputPath,true);
			FileChannel writeChannel = fos.getChannel();
			RandomAccessFile rf;
			FileChannel readChannel;
			long currentChannelSize;
			long transferedBytes;
			for(String filePath : filePaths) {
				if(!new File(filePath).exists())
					continue;
				
				rf = new RandomAccessFile(filePath,"r");
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
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void printLocationToSamFormatFromBufferedReader(long filePointerToBestLocation,long filePointerToSecondBestLocation,double scoreOfBestLocation, double scoreOfSecondBestLocation, BufferedRandomAccessFile br, boolean isMultiMapping, HashMap<String,BufferedWriter> chr2writer) throws Exception {
		
		if(filePointerToSecondBestLocation != -1 && !this.printSecondBestChr && (scoreOfBestLocation - scoreOfSecondBestLocation) <= this.scoreDiffCutoff) {
			return;
		}
		
		//print the best hit
		long currentFilePointer = br.getFilePointer();
		br.seek(filePointerToBestLocation);
		String locationInfo = br.getNextLine();
		StringTokenizer st = new StringTokenizer(locationInfo,"\t");
		//context id
		st.nextToken();
		//read id
		String readId = st.nextToken();
		String chr = st.nextToken();
		String strand = st.nextToken();
		int flag = 0;
		if(strand.equals("-"))
			flag += 16;
		
		int startA;
		ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
		Pattern commaPattern = Pattern.compile(",");
		String[] splittedCoordinates = commaPattern.split(st.nextToken());
		for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
			coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
		}
		int editDistance = Integer.valueOf(st.nextToken());
		
		//tmpLineBuilder.append("\t").append(mismatches).append("\t").append(hasSpliceSignal).append("\t").append(overlapsKnownJunction).append("\t").append(trimDouble(readScore)).append("\t").append(trimDouble(locationScore)).append("\t").append(overallMappingCount).append("\t").append(overallValidPairsCount);
		//splice signal
		char strandOfSpliceSignal = st.nextToken().charAt(0);
		
		//known junction
		st.nextToken();
		//read score
		st.nextToken();
		//location score
		st.nextToken();
		
		int overallMappingCount = Integer.valueOf(st.nextToken());
		
		//overall valid pair count
		st.nextToken();
		
		char overlapsPolyAtail = st.nextToken().charAt(0);
		
		int mappingQuality = 255;
		String cigar = "";
		int readLength = 0;
		int upstreamClippingIndex = 0;
		int downstreamClippingIndex = 0;
		int downstreamClippingLength = 0;
		
		boolean clippedAtStart = false;
		boolean clippedAtEnd = false;
		
		//check if alignment is clipped at the start
		if(coordinates.get(coordinates.size() - 1).getFirst() < 0) {
			upstreamClippingIndex = coordinates.get(coordinates.size() - 1).getFirst() * -1;
			readLength += upstreamClippingIndex;
			cigar = String.format("%sS",upstreamClippingIndex);
			clippedAtStart = true;
		}
		
		//check if alignment is clipped at the end
		if(coordinates.get(coordinates.size() - 1).getSecond() < 0) {
			downstreamClippingIndex = coordinates.get(coordinates.size() - 1).getSecond() * -1;
			clippedAtEnd = true;
		}
		
		
		//process matching parts and indels
		char currentCigSymbol;
		int segmentLength = coordinates.get(0).getSecond() - coordinates.get(0).getFirst() + 1;
		if(clippedAtStart)
			segmentLength -= upstreamClippingIndex;
		
		if(clippedAtEnd && coordinates.size() == 2) {
			//final read length already reached with the first matching block
			downstreamClippingLength = (readLength + segmentLength) - (downstreamClippingIndex  + 1);
			segmentLength -= downstreamClippingLength;
		}
		
		int insertionLength = 0;
		
		if(segmentLength > 0) {
			cigar += String.format("%sM", segmentLength);
			readLength += segmentLength;
		}
		
		for(int i = 1; i < coordinates.size(); i++) {
			if(coordinates.get(i).getFirst() <= 0)
				break;
		
			insertionLength = 0;
			if(coordinates.get(i).getFirst() - coordinates.get(i-1).getSecond() <= 0) {
				segmentLength = (-1 *(coordinates.get(i).getFirst() - coordinates.get(i-1).getSecond() - 1));
				insertionLength = segmentLength;
				editDistance += segmentLength;
				currentCigSymbol = 'I';
				readLength += segmentLength;
			}
			else {
				segmentLength = (coordinates.get(i).getFirst() - coordinates.get(i-1).getSecond() - 1);
				if(segmentLength > maxDelSize)
					currentCigSymbol = 'N';
				else {
					currentCigSymbol = 'D';
					editDistance += segmentLength;
				}
			}
			cigar += String.format("%s%s", segmentLength, currentCigSymbol);
			
			
			segmentLength = coordinates.get(i).getSecond() - coordinates.get(i).getFirst() + 1 - insertionLength;
			
			if(clippedAtEnd && i == coordinates.size() - 2) {
				//final read length reached with this last matching block
				downstreamClippingLength = (readLength + segmentLength) - (downstreamClippingIndex  + 1);
				segmentLength -= downstreamClippingLength;
			}
			
			if(segmentLength > 0) {
				cigar += String.format("%sM", segmentLength);
				readLength += segmentLength;
			}
			
		}
		
		//add clipping to alignment end
		if(clippedAtEnd && downstreamClippingLength > 0) {
			cigar += String.format("%sS",downstreamClippingLength);
			readLength += downstreamClippingLength;
		}
		
		//determined downstream clipping length == 0
		else if(clippedAtEnd)
			clippedAtEnd = false;
		

		if(!chr2writer.containsKey(chr)) {
			chr2writer.put(chr, new BufferedWriter(new FileWriter(new File(this.outputPath + "." + chr))));
		}
		BufferedWriter pw = chr2writer.get(chr);
		
		pw.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNM:i:%s\tNH:i:%s",readId,flag,chr,coordinates.get(0).getFirst() + upstreamClippingIndex,mappingQuality,cigar,"*",0,0,"*","*",editDistance,overallMappingCount));
		if(strandOfSpliceSignal != '0') {
			pw.write(String.format("\tXS:A:%s",strandOfSpliceSignal));
		}
		if(overlapsPolyAtail == '1') {
			if(clippedAtEnd) {
				int alignmentEnd = getAlignmentEnd(coordinates.get(0).getFirst() + upstreamClippingIndex,cigar);
				pw.write(String.format("\tPT:i:%s",alignmentEnd + 1));
			}
			else
				pw.write(String.format("\tPT:i:%s", coordinates.get(0).getFirst() + upstreamClippingIndex - 1));
		}
		
		
		if(filePointerToSecondBestLocation != -1) {
		//print second best hit
		br.seek(filePointerToSecondBestLocation);
		locationInfo = br.getNextLine();
		st = new StringTokenizer(locationInfo,"\t");
		//context id
		st.nextToken();
		//read id
		st.nextToken();
		chr = st.nextToken();
		//strand
		st.nextToken();
		//startA
		startA = Integer.valueOf(commaPattern.split(st.nextToken())[0]);
		
		pw.write(String.format("\tCC:Z:%s\tCP:i:%s\tS1:f:%s\tS2:f:%s", chr,startA,trimDouble(scoreOfBestLocation),trimDouble(scoreOfSecondBestLocation)));
		}
		
		pw.newLine();
		//set reader back to previous position
		if(!isMultiMapping)
			br.seek(currentFilePointer);
	}
	
	
	private int getAlignmentEnd(int alignmentStart, String cigar) {
		Pattern matchingPattern = Pattern.compile("[0-9]+[M|X|=|N|D|I]");
    	ArrayList<String> matches = new ArrayList<String>();
		Matcher matcher;
		
		
		matcher = matchingPattern.matcher(cigar);
		int currentPosition = alignmentStart;
		String currentMatch;
		int currentLength;
		String currentOp;
		while(matcher.find()) {
			matches.add(matcher.group());
		}
		
		
		for(int i = 0; i < matches.size(); i++) {
			currentMatch = matches.get(i);
			currentLength = Integer.valueOf(currentMatch.substring(0,currentMatch.length() - 1));
			currentOp = currentMatch.substring(currentMatch.length() - 1);
			
			if(currentOp.equals("M") || currentOp.equals("X") || currentOp.equals("="))
				currentPosition += currentLength - 1;
			
			else if(currentOp.equals("N") || currentOp.equals("D")) {
				currentPosition += currentLength + 1;
			}
			
			else if(currentOp.equals("I") && i < matches.size() - 1 && !matches.get(i+1).substring(matches.get(i+1).length() - 1).equals("N") && !matches.get(i+1).substring(matches.get(i+1).length() - 1).equals("D"))
				currentPosition++;
			
		}
		return currentPosition;
	}
	
	private class LocationContainer {
		
		private long pointerToBestLocation;
		private long pointerToSecondBestLocation;
		private double scoreOfBestLocation;
		private double scoreOfSecondBestLocation;
		
		public LocationContainer(long pointerToBestLocation, long pointerToSecondBestLocation, double scoreOfBestLocation, double scoreOfSecondBestLocation) {
			this.pointerToBestLocation = pointerToBestLocation;
			this.pointerToSecondBestLocation = pointerToSecondBestLocation;
			this.scoreOfBestLocation = scoreOfBestLocation;
			this.scoreOfSecondBestLocation = scoreOfSecondBestLocation;
		}
		
		public long getPointerToBestLocation() {
			return this.pointerToBestLocation;
		}
		
		public long getPointerToSecondBestLocation() {
			return this.pointerToSecondBestLocation;
		}
		
		public double getScoreOfBestLocation() {
			return this.scoreOfBestLocation;
		}
		
		public double getScoreOfSecondBestLocation() {
			return this.scoreOfSecondBestLocation;
		}
		
		
	}
	
	private class LocationContainerComparator implements Comparator<LocationContainer> {

		@Override
		public int compare(LocationContainer o1, LocationContainer o2) {
			return Long.valueOf(o1.getPointerToBestLocation()).compareTo(o2.getPointerToBestLocation());
		}
		
	}
	
	private String trimDouble(double inValue){
		DecimalFormat twoDec = new DecimalFormat("0.00",new DecimalFormatSymbols(Locale.US));
		twoDec.setGroupingUsed(false);
		return twoDec.format(inValue);
		}
	
	
	
	private void setReadScore(SparseRead read) {
		Collections.sort(read.getLocations(),this.sparseLocationScoreComparator);
		read.setTopScoringLocation(read.getLocations().get(read.getLocations().size() - 1));
	}
	
	
	private void getWindowCoverages(HashMap<String,ArrayList<SparseReadLocation>> chr2locationsStartSorted, HashMap<String,ArrayList<SparseReadLocation>> chr2locationsEndSorted,HashMap<String,ArrayList<SparseReadLocation>> chr2ReverselocationsStartSorted, HashMap<String,ArrayList<SparseReadLocation>> chr2ReverselocationsEndSorted) {
		try {
			int maxThreads = this.numberOfThreads;
			//if(maxThreads > 3) maxThreads = 3;
			ExecutorService executor = Executors.newFixedThreadPool(maxThreads);
			ArrayList<Future> futures = new ArrayList<Future>();
			IntervalTree intervalTree;
			//process forward strand
			for(String chr : chr2locationsStartSorted.keySet()) {
				intervalTree = null;
				if(this.chr2intervalTree != null && this.chr2intervalTree.containsKey(chr))
					intervalTree = this.chr2intervalTree.get(chr);
				
				futures.add(executor.submit(new CoverageCalculator(this.windowSizes, this.completeWindowIntervall, this.maxContextSize, chr2locationsStartSorted.get(chr), chr2locationsEndSorted.get(chr),intervalTree)));
			}
			
			//process reverse complemented strand (in case we don't have strand specific reads, we have an empty hashmap here).
			for(String chr : chr2ReverselocationsStartSorted.keySet()) {
				intervalTree = null;
				if(this.chr2intervalTree != null && this.chr2intervalTree.containsKey(chr))
					intervalTree = this.chr2intervalTree.get(chr);
				futures.add(executor.submit(new CoverageCalculator(this.windowSizes, this.completeWindowIntervall, this.maxContextSize, chr2ReverselocationsStartSorted.get(chr), chr2ReverselocationsEndSorted.get(chr),intervalTree)));
			}
			
			executor.shutdown();
			for(Future future : futures)
				future.get();
			
			futures.clear();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private class SparseLocationStartComparator implements Comparator<SparseReadLocation> {
		public int compare(SparseReadLocation l1, SparseReadLocation l2) {
			int startA = l1.getCoordinates().get(0).getFirst();
			int startB = l2.getCoordinates().get(0).getFirst();
			return ((Integer)startA).compareTo(startB);
		}
		
	}
	
	private class SparseLocationEndComparator implements Comparator<SparseReadLocation> {
		public int compare(SparseReadLocation l1, SparseReadLocation l2) {
			int end1 = l1.getCoordinates().get(0).getSecond();
			if(l1.getCoordinates().size() > 1 && l1.getCoordinates().get(1).getFirst() != -1)
				end1 = l1.getCoordinates().get(l1.getCoordinates().size()-1).getSecond();
			
			int end2 = l2.getCoordinates().get(0).getSecond();
			if(l2.getCoordinates().size() > 1 && l2.getCoordinates().get(1).getFirst() != -1)
				end2 = l2.getCoordinates().get(l2.getCoordinates().size()-1).getSecond();
			
			return ((Integer)end1).compareTo(end2);
		}
		
	}
	
	private class SparseLocationScoreComparator implements Comparator<SparseReadLocation> {
		public int compare(SparseReadLocation l1, SparseReadLocation l2) {
			double score1 = l1.getScore();
			double score2 = l2.getScore();
			return ((Double)score1).compareTo(score2);
		}
		
	}
	
	
	private void processGlobalContext(HashMap<String,BufferedWriter> chr2writer) throws Exception {
		BufferedRandomAccessFile br = new BufferedRandomAccessFile(new File(this.multiMappingFilePath),"r",10240);
		BufferedRandomAccessFile bufferedReaderForMultiMappings = new BufferedRandomAccessFile(new File(this.multiMappingFilePath),"r",10240);
		
		long prevFilePointer = br.getFilePointer();
		String currentLine = br.getNextLine();
		if(currentLine == null) {
			br.close();
			bufferedReaderForMultiMappings.close();
			return;
		}
			
		long filePointer = br.getFilePointer();
		int overallParsedReads = 0;
		int overallAddedLocations = 0;
		int overallAddedReads = 0;
		StringTokenizer st = new StringTokenizer(currentLine,"\t");
		Pattern commaPattern = Pattern.compile(",");
		String[] splittedCoordinates;
		String currentContextId = st.nextToken();
		String prevContextId = currentContextId;
		String readId = st.nextToken();
		String prevReadId = readId;
		String chr = st.nextToken();
		char strand = st.nextToken().toCharArray()[0];
		splittedCoordinates = commaPattern.split(st.nextToken());
		ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
		for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
			coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
		}
		char mappingType;
		if(coordinates.size() == 1 || coordinates.get(1).getFirst() == -1)
			mappingType = 'F';
		else
			mappingType = 'S';
		
		
		//parse location score
				st.nextToken();
				st.nextToken();
				st.nextToken();
				st.nextToken();
				
		double locationScore = Double.valueOf(st.nextToken());
		
		
		ArrayList<MultiSparseReadLocation> currentSparseReadLocations = new ArrayList<MultiSparseReadLocation>();
		ArrayList<String> currentChromosomes = new ArrayList<String>();
		MultiSparseReadLocation tmpSparseLocation;
		MultiSparseReadLocation currentSparseLocation = new MultiSparseReadLocation(coordinates, strand);
		SparseReadLocation tmpUniqueSparseReadLocation;
		currentSparseLocation.setFilePointer(prevFilePointer);
		currentSparseLocation.setScore(locationScore);
		currentSparseReadLocations.add(currentSparseLocation);
		currentChromosomes.add(chr);
		
		ArrayList<LocationContainer> locations2write = new ArrayList<LocationContainer>();
		ArrayList<LocationContainer> multiLocations2write = new ArrayList<LocationContainer>();
		LocationContainerComparator locationContainerComparator = new LocationContainerComparator();
		
		while((currentLine = br.getNextLine()) != null) {
			st = new StringTokenizer(currentLine,"\t");
			currentContextId = st.nextToken();
			readId = st.nextToken();
			chr = st.nextToken();
			strand = st.nextToken().toCharArray()[0];
			splittedCoordinates = commaPattern.split(st.nextToken());
			coordinates = new ArrayList<Pair<Integer,Integer>>(5);
			for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
				coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
			}
			if(coordinates.size() == 1 || coordinates.get(1).getFirst() == -1)
				mappingType = 'F';
			else
				mappingType = 'S';
			//mismatches
			st.nextToken();
			//known splice signal
			st.nextToken();
			//overlaps known junction
			st.nextToken();
			//read score
			st.nextToken();
			
			locationScore = Double.valueOf(st.nextToken());
			
			
			if(readId.equals(prevReadId)) {
				currentSparseLocation = new MultiSparseReadLocation(coordinates, strand);
				currentSparseReadLocations.add(currentSparseLocation);
				currentSparseLocation.setFilePointer(filePointer);
				currentSparseLocation.setScore(locationScore);
				currentChromosomes.add(chr);
			}
			
			else {
				overallParsedReads++;
				//in case we have more than one read location for the current read we add this read to the actual context
				if(currentSparseReadLocations.size() > 1) {
					SparseRead currentRead = new SparseRead();
					for(int i = 0; i < currentSparseReadLocations.size(); i++) {
						tmpSparseLocation = currentSparseReadLocations.get(i);
						tmpSparseLocation.setAmbiguityFactor(currentSparseReadLocations.size());
						currentRead.addLocation(tmpSparseLocation);
						overallAddedLocations++;
					}
					overallAddedReads++;
					//reads.add(currentRead);
					
					setReadScore(currentRead);
					processMultiMappedRead(currentRead, bufferedReaderForMultiMappings,locations2write, multiLocations2write, chr2writer,locationContainerComparator);
					
				}
				
				else {
					printLocationToSamFormatFromBufferedReader(currentSparseReadLocations.get(0).getFilePointer(),-1,-1,-1,br,false,chr2writer);
				}
				
				prevReadId = readId;
				currentSparseReadLocations.clear();
				currentChromosomes.clear();
				
				
				currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
				currentSparseLocation.setFilePointer(filePointer);
				currentSparseLocation.setScore(locationScore);
				currentSparseReadLocations.add(currentSparseLocation);
				currentChromosomes.add(chr);
			}
			filePointer = br.getFilePointer();
		}
		
		//process last bunch...
		if(currentSparseReadLocations.size() == 1) {
			printLocationToSamFormatFromBufferedReader(currentSparseReadLocations.get(0).getFilePointer(),-1,-1,-1,br,false,chr2writer);
		}
		
		else {
			SparseRead currentRead = new SparseRead();
			for(int i = 0; i < currentSparseReadLocations.size(); i++) {
				tmpSparseLocation = currentSparseReadLocations.get(i);
				tmpSparseLocation.setAmbiguityFactor(currentSparseReadLocations.size());
				currentRead.addLocation(tmpSparseLocation);
				overallAddedLocations++;
			}
			overallAddedReads++;
			//reads.add(currentRead);
			
			setReadScore(currentRead);
			processMultiMappedRead(currentRead, bufferedReaderForMultiMappings,locations2write, multiLocations2write, chr2writer,locationContainerComparator);
		}
		
		//write last bunch of reads
		Collections.sort(locations2write, locationContainerComparator);
		for(LocationContainer container : locations2write) {
			printLocationToSamFormatFromBufferedReader(container.getPointerToBestLocation(),container.getPointerToSecondBestLocation(),container.getScoreOfBestLocation(),container.getScoreOfSecondBestLocation(),br,true,chr2writer);
		}
		locations2write.clear();
		br.close();
		bufferedReaderForMultiMappings.close();
		
	}
	
	
	private void processMultiMappedRead(SparseRead read, BufferedRandomAccessFile bufferedReader,ArrayList<LocationContainer> locations2write, ArrayList<LocationContainer> multiLocations2write, HashMap<String,BufferedWriter> chr2writer,LocationContainerComparator locationContainerComparator) throws Exception {
				
		IntervalTree<Microbe> genomeTree;
		String chr;
		String nextHitChr;
		String bestHitMicrobeId;
		String nextHitMicrobeId;
		MultiSparseReadLocation nextHitMicrobe;
		ArrayList<MultiSparseReadLocation> locations;
		LocationContainer tmpContainer = new LocationContainer(read.getTopScoringLocation().getFilePointer(),read.getLocations().get(read.getLocations().size()-2).getFilePointer(),read.getTopScoringLocation().getScore(),read.getLocations().get(read.getLocations().size()-2).getScore());
		
		
		if(this.printSecondBestChr && this.chr2intervalTree != null) {
			//check if there is a hit outside the current chr (species)
			//currentFilePointer = bufferedReader.getFilePointer();
			bufferedReader.seek(read.getTopScoringLocation().getFilePointer());
			chr = bufferedReader.getNextLine().split("\t")[2];
			if(this.chr2intervalTree.containsKey(chr)) {
				genomeTree = this.chr2intervalTree.get(chr);
				locations = read.getLocations();
				if(!genomeTree.getIntervalsSpanning(read.getTopScoringLocation().getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).isEmpty())
					bestHitMicrobeId = genomeTree.getIntervalsSpanning(read.getTopScoringLocation().getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).get(0).getId();
				else
					bestHitMicrobeId = genomeTree.getIntervalsSpanning(read.getTopScoringLocation().getCoordinates().get(0).getSecond(),new ArrayList<Microbe>()).get(0).getId();
				
				for(int i = locations.size() - 2; i >= 0; i--) {
					nextHitMicrobe = locations.get(i);
					bufferedReader.seek(nextHitMicrobe.getFilePointer());
					nextHitChr = bufferedReader.getNextLine().split("\t")[2];
					if(nextHitChr.equals(chr)) {
						if(!genomeTree.getIntervalsSpanning(nextHitMicrobe.getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).isEmpty())
							nextHitMicrobeId = genomeTree.getIntervalsSpanning(nextHitMicrobe.getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).get(0).getId();
						else
							nextHitMicrobeId = genomeTree.getIntervalsSpanning(nextHitMicrobe.getCoordinates().get(0).getSecond(),new ArrayList<Microbe>()).get(0).getId();
							
						if(!nextHitMicrobeId.equals(bestHitMicrobeId)) {
							tmpContainer = new LocationContainer(read.getTopScoringLocation().getFilePointer(),nextHitMicrobe.getFilePointer(),read.getTopScoringLocation().getScore(),nextHitMicrobe.getScore());
							break;
						}
					}
					else {
						tmpContainer = new LocationContainer(read.getTopScoringLocation().getFilePointer(),nextHitMicrobe.getFilePointer(),read.getTopScoringLocation().getScore(),nextHitMicrobe.getScore());
						break;
					}
				}
			}
			//bufferedReader.seek(currentFilePointer);
		}
		
		else if(this.printSecondBestChr) {
			//currentFilePointer = bufferedReader.getFilePointer();
			bufferedReader.seek(read.getTopScoringLocation().getFilePointer());
			chr = bufferedReader.getNextLine().split("\t")[2];
			locations = read.getLocations();
			for(int i = locations.size() - 2; i >= 0; i--) {
				nextHitMicrobe = locations.get(i);
				bufferedReader.seek(nextHitMicrobe.getFilePointer());
				nextHitChr = bufferedReader.getNextLine().split("\t")[2];
				if(!chr.equals(nextHitChr)) {
					tmpContainer = new LocationContainer(read.getTopScoringLocation().getFilePointer(),nextHitMicrobe.getFilePointer(),read.getTopScoringLocation().getScore(),nextHitMicrobe.getScore());
					break;
				}
			}
			//bufferedReader.seek(currentFilePointer);
		}
		
		
		
		locations2write.add(tmpContainer);
		
		if(locations2write.size() >= 2000000) {
			//sort by file pointer, this way we can use the buffer for seek operations
			Collections.sort(locations2write, locationContainerComparator);
			for(LocationContainer container : locations2write) {
				printLocationToSamFormatFromBufferedReader(container.getPointerToBestLocation(),container.getPointerToSecondBestLocation(),container.getScoreOfBestLocation(),container.getScoreOfSecondBestLocation(),bufferedReader,true,chr2writer);
			}
			locations2write.clear();
		}
		
		//output multi mapping here
		if(this.printMultiMappings) {
			for(MultiSparseReadLocation tmpLocation : read.getLocations()) {
				if(!tmpLocation.equals(read.getTopScoringLocation())) {
					tmpContainer = new LocationContainer(tmpLocation.getFilePointer(),-1,tmpLocation.getScore(),Double.MIN_VALUE);
					multiLocations2write.add(tmpContainer);
				}
			}
			if(multiLocations2write.size() >= 2000000) {
				//sort by file pointer, this way we can use the buffer for seek operations
				Collections.sort(multiLocations2write, locationContainerComparator);
				for(LocationContainer container : multiLocations2write) {
					printLocationToSamFormatFromBufferedReader(container.getPointerToBestLocation(),container.getPointerToSecondBestLocation(),container.getScoreOfBestLocation(),container.getScoreOfSecondBestLocation(),bufferedReader,true,chr2writer);
				}
				multiLocations2write.clear();
			}
		}
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
