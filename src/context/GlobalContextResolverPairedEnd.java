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
import main.Read;
import main.ReadLocation;
import main.ReadPair;
import main.SparseRead;
import main.MultiSparseReadLocation;
import main.Pair;
import main.SparseReadLocation;
import tools.BufferedRandomAccessFile;
import tools.MaxPriorityQueue;

public class GlobalContextResolverPairedEnd implements ContextResolver {
	
	private String multiMappingFilePath;
	private String outputPath;
	private ArrayList<Integer> windowSizes;
	private int completeWindowIntervall;
	private int maxContextSize;
	private int maxDelSize;
	private int numberOfThreads;
	private int printedReads;
	private int overallProcessedReads;
	private PairScoreComparator pairScoreComparator;
	private SparseLocationScoreComparator sparseReadLocationComparator;
	
	private boolean updateQueue;
	private boolean printMultiMappings;
	private boolean printSecondBestChr;
	private int updateInterval;
	private boolean verbose;
	private boolean strandSpecific;
	
	private ArrayList<ReadPair<SparseRead,SparseRead>> reads;
	private MaxPriorityQueue queue;
	
	private HashMap<String,Integer> chr2length;
	private HashMap<String,IntervalTree<Microbe>> chr2intervalTree;
		
	private LocationContainerComparator locationContainerComparator;
	
	
	private double scoreDiffCutoff;
	private DecimalFormat twoDec;
	

	public GlobalContextResolverPairedEnd(String multiMappingFilePath, String outputPath,HashMap<String,Integer> chr2length, HashMap<String,IntervalTree<Microbe>> chr2intervalTree, ArrayList<Integer> windowSizes, int readLength, int maxContextSize, int maxDelSize, boolean localContexts,boolean updateQueue, int updateInterval, int numberOfThreads, boolean printMultiMappings, boolean printSecondBestChr, boolean verbose, boolean strandSpecific) {
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
		
		this.pairScoreComparator = new PairScoreComparator();
		this.sparseReadLocationComparator = new SparseLocationScoreComparator();
		
		this.updateQueue = updateQueue;
		this.updateInterval = updateInterval;
		this.numberOfThreads = numberOfThreads;
		this.printMultiMappings = printMultiMappings;
		this.printSecondBestChr = printSecondBestChr;
		this.verbose = verbose;
		this.strandSpecific = strandSpecific;
		this.printedReads = 0;
		this.overallProcessedReads = 0;
		
		this.reads = new ArrayList<ReadPair<SparseRead,SparseRead>>();
		this.queue = null;
		
		this.locationContainerComparator = new LocationContainerComparator();
		
		//this.scoreDiffCutoff = new CoverageCalculator(this.windowSizes).getScoreCutoff();
		this.scoreDiffCutoff = 0;
		this.twoDec = new DecimalFormat("0.00",new DecimalFormatSymbols(Locale.US));
		
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
	
	
	private void processReadPair(SparseRead firstMate, SparseRead secondMate, ArrayList<Pair<Integer,Integer>> validPairs, MultiSparseReadLocation topLocationFirstMate, MultiSparseReadLocation secondBestLocationFirstMate, 
			MultiSparseReadLocation topLocationSecondMate, MultiSparseReadLocation secondBestLocationSecondMate, BufferedRandomAccessFile bufferedReader,ArrayList<LocationContainer> locations2write,
									  ArrayList<LocationContainer> multiLocations2write,
									  HashMap<String,BufferedWriter> chr2writer, boolean nextSegmentUnmapped) throws Exception {
		
		LocationContainer tmpContainer = null;
		
		//get all relevant file pointers and score for the first mate
		long bestLocationFirstMateFilePointer = topLocationFirstMate.getFilePointer();
		double scoreBestLocationFirstMate = topLocationFirstMate.getScore();
		if(!validPairs.isEmpty())
			scoreBestLocationFirstMate = validPairs.get(validPairs.size()-1).getScore();
		
		double scoreSecondBestLocationFirstMate = Double.MIN_VALUE;
		long secondBestLocationFirstMateFilePointer = -1;
		if(secondBestLocationFirstMate != null) {
			scoreSecondBestLocationFirstMate = secondBestLocationFirstMate.getScore();
			if(!validPairs.isEmpty())
				scoreSecondBestLocationFirstMate = validPairs.get(validPairs.size()-2).getScore();
			
			secondBestLocationFirstMateFilePointer = secondBestLocationFirstMate.getFilePointer();
		}
		
		//get all relevant file pointers and score for the second mate
		long bestLocationSecondMateFilePointer = -1;
		double scoreBestLocationSecondMate = Double.MIN_VALUE;
		double scoreSecondBestLocationSecondMate = Double.MIN_VALUE;
		long secondBestLocationSecondMateFilePointer = -1;
		if(secondMate != null) {
			bestLocationSecondMateFilePointer = topLocationSecondMate.getFilePointer();
			scoreBestLocationSecondMate = topLocationSecondMate.getScore();
			if(!validPairs.isEmpty())
				scoreBestLocationSecondMate = validPairs.get(validPairs.size()-1).getScore();
			
			if(secondBestLocationSecondMate != null) {
				scoreSecondBestLocationSecondMate = secondBestLocationSecondMate.getScore();
				if(!validPairs.isEmpty())
					scoreSecondBestLocationSecondMate = validPairs.get(validPairs.size()-2).getScore();
				
				secondBestLocationSecondMateFilePointer = secondBestLocationSecondMate.getFilePointer();
			}
		}
		
		
		
		tmpContainer = new LocationContainer(bestLocationFirstMateFilePointer,secondBestLocationFirstMateFilePointer,bestLocationSecondMateFilePointer,secondBestLocationSecondMateFilePointer,scoreBestLocationFirstMate,scoreSecondBestLocationFirstMate,scoreBestLocationSecondMate,scoreSecondBestLocationSecondMate,nextSegmentUnmapped);		
		locations2write.add(tmpContainer);
		this.printedReads++;
		
		if(locations2write.size() >= 2000000) {
			//sort by file pointer, this way we can use the buffer for seek operations
			Collections.sort(locations2write, this.locationContainerComparator);
			for(LocationContainer container : locations2write) {
				if(container.getPointerToSecondBestLocationFirstMate() != -1)
					printLocationToSamFormatFromBufferedReader(container.getPointerToBestLocationFirstMate(),container.getPointerToSecondBestLocationFirstMate(),container.getPointerToBestLocationSecondMate(),container.getPointerToSecondBestLocationSecondMate(), 
															   container.getScoreOfBestLocationFirstMate(),container.getScoreOfSecondBestLocationFirstMate(),container.getScoreOfBestLocationSecondMate(),container.getScoreOfSecondBestLocationSecondMate(),
															   bufferedReader,true,chr2writer,container.isNextSegmentUnmapped());
				else
					printLocationToSamFormatFromBufferedReader(container.getPointerToBestLocationFirstMate(),container.getPointerToSecondBestLocationFirstMate(),container.getPointerToBestLocationSecondMate(),container.getPointerToSecondBestLocationSecondMate(), 
							   container.getScoreOfBestLocationFirstMate(),container.getScoreOfSecondBestLocationFirstMate(),container.getScoreOfBestLocationSecondMate(),container.getScoreOfSecondBestLocationSecondMate(),
							   bufferedReader,false,chr2writer,container.isNextSegmentUnmapped());
			}
			locations2write.clear();
		}
	}
	

	
	
	private void printLocationToSamFormatFromBufferedReader(long filePointerToBestLocationFirstMate,long filePointerToSecondBestLocationFirstMate, long filePointerToBestLocationSecondMate,long filePointerToSecondBestLocationSecondMate, double scoreOfBestLocationFirstMate, double scoreOfSecondBestLocationFirstMate,
			                                                  double scoreOfBestLocationSecondMate, double scoreOfSecondBestLocationSecondMate, BufferedRandomAccessFile br, 
			                                                  boolean isMultiMapping, HashMap<String,BufferedWriter> chr2writer, boolean nextSegmentUnmapped) throws Exception {
		
		if(filePointerToSecondBestLocationFirstMate != -1 && !this.printSecondBestChr && (scoreOfBestLocationFirstMate - scoreOfSecondBestLocationFirstMate) <= this.scoreDiffCutoff) {
			return;
			
		}
		

		long currentFilePointer = br.getFilePointer();
		
		br.seek(filePointerToBestLocationFirstMate);
		String locationInfoBestHitFirstMate = br.getNextLine();
		SamRecord samRecordBestHitFirstMate = new SamRecord(locationInfoBestHitFirstMate,this.maxDelSize,nextSegmentUnmapped);
		
		
		String locationInfoBestHitSecondMate = null;
		SamRecord samRecordBestHitSecondMate = null;
		if(filePointerToBestLocationSecondMate != -1) {
			br.seek(filePointerToBestLocationSecondMate);
			locationInfoBestHitSecondMate = br.getNextLine();
			samRecordBestHitSecondMate = new SamRecord(locationInfoBestHitSecondMate,this.maxDelSize,nextSegmentUnmapped);
		}
		
		
		String locationInfoSecondBestHitFirstMate = null;
		SamRecord samRecordSecondBestHitFirstMate = null;
		if(filePointerToSecondBestLocationFirstMate != -1) {
			br.seek(filePointerToSecondBestLocationFirstMate);
			locationInfoSecondBestHitFirstMate = br.getNextLine();
			samRecordSecondBestHitFirstMate = new SamRecord(locationInfoSecondBestHitFirstMate,this.maxDelSize,nextSegmentUnmapped);
		}
		String locationInfoSecondBestHitSecondMate = null;
		SamRecord samRecordSecondBestHitSecondMate = null;
		if(filePointerToSecondBestLocationSecondMate != -1) {
			br.seek(filePointerToSecondBestLocationSecondMate);
			locationInfoSecondBestHitSecondMate = br.getNextLine();
			samRecordSecondBestHitSecondMate = new SamRecord(locationInfoSecondBestHitSecondMate,this.maxDelSize,nextSegmentUnmapped);
		}
		
		
		String rNext = "*";
		int pNextFirstMate = 0;
		int pNextSecondMate = 0;
		int tLenFirstMate = 0;
		int tLenSecondMate = 0;
		int rightmostBase;
		
		
		if(samRecordBestHitSecondMate != null) {
			if(samRecordBestHitFirstMate.getReferenceName().equals(samRecordBestHitSecondMate.getReferenceName())) {
				rNext = "=";
				pNextFirstMate = samRecordBestHitSecondMate.getAlignmentStart();
				pNextSecondMate = samRecordBestHitFirstMate.getAlignmentStart();
				
				if(samRecordBestHitFirstMate.getAlignmentStart() <= samRecordBestHitSecondMate.getAlignmentStart()) {
					//leftmost base is startA. It could happen that the end of A is larger than of B (i.e. by a split alignment)
					rightmostBase = (samRecordBestHitSecondMate.getAlignmentEnd() >= samRecordBestHitFirstMate.getAlignmentEnd())?samRecordBestHitSecondMate.getAlignmentEnd():samRecordBestHitFirstMate.getAlignmentEnd();
					tLenFirstMate = rightmostBase - samRecordBestHitFirstMate.getAlignmentStart() + 1;
					tLenSecondMate = -1 * tLenFirstMate;
				}
				else {
					//leftmost base is startB.
					rightmostBase = (samRecordBestHitFirstMate.getAlignmentEnd() >= samRecordBestHitSecondMate.getAlignmentEnd())?samRecordBestHitFirstMate.getAlignmentEnd():samRecordBestHitSecondMate.getAlignmentEnd();
					tLenSecondMate = rightmostBase - samRecordBestHitSecondMate.getAlignmentStart() + 1;
					tLenFirstMate = -1 * tLenSecondMate;
				}
			}
			else {
				rNext = samRecordBestHitSecondMate.getReferenceName();
			}
			
			
		}
		
		
		
		//write the first mate
		
		if(!chr2writer.containsKey(samRecordBestHitFirstMate.getReferenceName())) {
			chr2writer.put(samRecordBestHitFirstMate.getReferenceName(), new BufferedWriter(new FileWriter(new File(this.outputPath + "." + samRecordBestHitFirstMate.getReferenceName())),10240));
		}
		BufferedWriter pw = chr2writer.get(samRecordBestHitFirstMate.getReferenceName());
		
		
		pw.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNM:i:%s\tNH:i:%s",samRecordBestHitFirstMate.getId().substring(0,samRecordBestHitFirstMate.getId().lastIndexOf('/')),samRecordBestHitFirstMate.getFlag(),samRecordBestHitFirstMate.getReferenceName(),samRecordBestHitFirstMate.getAlignmentStart(),samRecordBestHitFirstMate.getMappingQuality(),samRecordBestHitFirstMate.getCigar(),rNext,pNextFirstMate,tLenFirstMate,samRecordBestHitFirstMate.getSeq(),samRecordBestHitFirstMate.getQual(),samRecordBestHitFirstMate.getEditDistance(),samRecordBestHitFirstMate.getNumberOfHits()));
		
		if(samRecordBestHitFirstMate.getXsTag() != '0')
			pw.write(String.format("\tXS:A:%s",samRecordBestHitFirstMate.getXsTag()));
		
		if(samRecordBestHitFirstMate.overlapsPolyAtail() == '1') {
			if(samRecordBestHitFirstMate.isClippedAtEnd())
				pw.write(String.format("\tPT:i:%s",samRecordBestHitFirstMate.getAlignmentEnd() + 1));
			else
				pw.write(String.format("\tPT:i:%s",samRecordBestHitFirstMate.getAlignmentStart() - 1));
		}
		
		if(samRecordSecondBestHitFirstMate != null)
			pw.write(String.format("\tCC:Z:%s\tCP:i:%s\tS1:f:%s\tS2:f:%s", samRecordSecondBestHitFirstMate.getReferenceName(),samRecordSecondBestHitFirstMate.getAlignmentStart(),trimDouble(scoreOfBestLocationFirstMate),trimDouble(scoreOfSecondBestLocationFirstMate)));
		
		pw.newLine();
		
		
		//write the second mate
		
		
		if(samRecordBestHitSecondMate != null) {
			
			//the current definition of a valid read pair does not allow this scenario... but we never know....
			if(!samRecordBestHitFirstMate.getReferenceName().equals(samRecordBestHitSecondMate)) {
				if(!chr2writer.containsKey(samRecordBestHitSecondMate.getReferenceName())) {
					chr2writer.put(samRecordBestHitSecondMate.getReferenceName(), new BufferedWriter(new FileWriter(new File(this.outputPath + "." + samRecordBestHitSecondMate.getReferenceName())),10240));
				}
				pw = chr2writer.get(samRecordBestHitSecondMate.getReferenceName());
			}
			
			pw.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNM:i:%s\tNH:i:%s",samRecordBestHitSecondMate.getId().substring(0,samRecordBestHitSecondMate.getId().lastIndexOf('/')),samRecordBestHitSecondMate.getFlag(),samRecordBestHitSecondMate.getReferenceName(),samRecordBestHitSecondMate.getAlignmentStart(),samRecordBestHitSecondMate.getMappingQuality(),samRecordBestHitSecondMate.getCigar(),rNext,pNextSecondMate,tLenSecondMate,samRecordBestHitSecondMate.getSeq(),samRecordBestHitSecondMate.getQual(),samRecordBestHitSecondMate.getEditDistance(),samRecordBestHitSecondMate.getNumberOfHits()));

			if(samRecordBestHitSecondMate.getXsTag() != '0')
				pw.write(String.format("\tXS:A:%s",samRecordBestHitSecondMate.getXsTag()));
			
			if(samRecordBestHitSecondMate.overlapsPolyAtail() == '1') {
				if(samRecordBestHitSecondMate.isClippedAtEnd())
					pw.write(String.format("\tPT:i:%s",samRecordBestHitSecondMate.getAlignmentEnd() + 1));
				else
					pw.write(String.format("\tPT:i:%s",samRecordBestHitSecondMate.getAlignmentStart() - 1));
			}
			
			if(samRecordSecondBestHitSecondMate != null)
				pw.write(String.format("\tCC:Z:%s\tCP:i:%s\tS1:f:%s\tS2:f:%s", samRecordSecondBestHitSecondMate.getReferenceName(),samRecordSecondBestHitSecondMate.getAlignmentStart(),trimDouble(scoreOfBestLocationSecondMate),trimDouble(scoreOfSecondBestLocationSecondMate)));
			
			
			pw.newLine();
		}
		
		//set reader back to previous position
		if(!isMultiMapping)
			br.seek(currentFilePointer);
	}
	
	
	private class SamRecord {
		
		private String id;
		private int flag;
		private String referenceName;
		private int alignmentStart;
		private int alignmentEnd;
		private int mappingQuality;
		private String cigar;
		private String rNext;
		private int pnext;
		private int tlen;
		private String seq;
		private String qual;
		private int editDistance;
		private int numberOfHits;
		private char xsTag;
		private char overlapsPolyAtail;
		
		private boolean clippedAtStart;
		private boolean clippedAtEnd;
		
		
		
		public SamRecord(String locationInfoContextMap, int maxDelSize, boolean nextSegmentUnmapped) {
			setSamTags(locationInfoContextMap, maxDelSize,nextSegmentUnmapped);
		}
		


		private void setSamTags(String locationInfoContextMap, int maxDelSize, boolean nextSegmentUnmapped) {
			StringTokenizer st = new StringTokenizer(locationInfoContextMap,"\t");
			//context id
			st.nextToken();
			//read id
			String readId = st.nextToken();
			boolean isFirstSegment;
			if(readId.substring(readId.lastIndexOf('/') + 1).equals("1"))
				isFirstSegment = true;
			else
				isFirstSegment = false;
			
			String chr = st.nextToken();
			String strand = st.nextToken();
			String rNext = "*";
			//template having multiple segments in sequencing
			int flag = 1;
			if(strand.equals("-"))
				//SEQ being reverse complemented
				flag += 16;
			
			if(isFirstSegment) {
				//the first segment in the template
				flag += 64;
			}
			
			else {
				//the last segment in the template
				flag += 128;
			}
			
			
			if(nextSegmentUnmapped) {
				//next segment in the template unmapped
				flag += 8;
			}
			else {
				//each segment properly aligned according to the aligner
				flag += 2;
				rNext = "=";
				if(strand.equals("+")) {
					//SEQ of the next segment in the template being reversed
					flag += 32;
				}
			}
			
			
			ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
			Pattern commaPattern = Pattern.compile(",");
			String[] splittedCoordinates = commaPattern.split(st.nextToken());
			for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
				coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
			}
			int editDistance = Integer.valueOf(st.nextToken());
			
			//splice signal
			char strandOfSpliceSignal = st.nextToken().charAt(0);
			//known junction
			st.nextToken();
			//read score
			st.nextToken();
			//location score
			st.nextToken();
			
			int overallMappingCount = Integer.valueOf(st.nextToken());
			int overallValidPairCount = Integer.valueOf(st.nextToken());
			char overlapsPolyAtail = st.nextToken().charAt(0);
			
			int numberOfHits = (overallValidPairCount > 0 ?overallValidPairCount:overallMappingCount);
			
			
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
			
			int start = coordinates.get(0).getFirst() + upstreamClippingIndex;
			
			//set the fields
			this.id = readId;
			this.flag = flag; 
			this.referenceName = chr;
			this.alignmentStart = start;
			this.alignmentEnd = getAlignmentEnd(start, cigar);
			this.mappingQuality = mappingQuality;
			this.cigar = cigar;
			this.rNext = rNext;
			this.pnext = 0;
			this.seq = "*";
			this.qual = "*";
			this.editDistance = editDistance;
			this.numberOfHits = numberOfHits;
			this.xsTag = strandOfSpliceSignal;
			this.overlapsPolyAtail = overlapsPolyAtail;
			this.clippedAtStart = clippedAtStart;
			this.clippedAtEnd = clippedAtEnd;
			
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
		
		
		public String getId() {
			return id;
		}

		public int getFlag() {
			return flag;
		}

		public String getReferenceName() {
			return referenceName;
		}

		public int getAlignmentStart() {
			return alignmentStart;
		}

		public int getAlignmentEnd() {
			return alignmentEnd;
		}

		public int getMappingQuality() {
			return mappingQuality;
		}

		public String getCigar() {
			return cigar;
		}

		public String getrNext() {
			return rNext;
		}

		public int getPnext() {
			return pnext;
		}

		public int getTlen() {
			return tlen;
		}

		public String getSeq() {
			return seq;
		}

		public String getQual() {
			return qual;
		}

		public int getEditDistance() {
			return editDistance;
		}

		public int getNumberOfHits() {
			return numberOfHits;
		}
		
		public char getXsTag() {
			return this.xsTag;
		}
		
		public char overlapsPolyAtail() {
			return this.overlapsPolyAtail;
		}
		
		public boolean isClippedAtStart() {
			return this.clippedAtStart;
		}
		
		public boolean isClippedAtEnd() {
			return this.clippedAtEnd;
		}
		
	}
	
	private class LocationContainer {
		
		private long pointerToBestLocationFirstMate;
		private long pointerToSecondBestLocationFirstMate;
		
		private long pointerToBestLocationSecondMate;
		private long pointerToSecondBestLocationSecondMate;
		
		private double scoreOfBestLocationFirstMate;
		private double scoreOfSecondBestLocationFirstMate;
		
		private double scoreOfBestLocationSecondMate;
		private double scoreOfSecondBestLocationSecondMate;
		
		private boolean nextSegmentUnmapped;
		
		public LocationContainer(long pointerToBestLocationFirstMate, long pointerToSecondBestLocationFirstMate, long pointerToBestLocationSecondMate, long pointerToSecondBestLocationSecondMate, double scoreOfBestLocationFirstMate, double scoreOfSecondBestLocationFirstMate, double scoreOfBestLocationSecondMate, double scoreOfSecondBestLocationSecondMate, boolean nextSegmentUnmapped) {
			this.pointerToBestLocationFirstMate = pointerToBestLocationFirstMate;
			this.pointerToSecondBestLocationFirstMate = pointerToSecondBestLocationFirstMate;
			this.scoreOfBestLocationFirstMate = scoreOfBestLocationFirstMate;
			this.scoreOfSecondBestLocationFirstMate = scoreOfSecondBestLocationFirstMate;
			
			
			this.pointerToBestLocationSecondMate = pointerToBestLocationSecondMate;
			this.pointerToSecondBestLocationSecondMate = pointerToSecondBestLocationSecondMate;
			this.scoreOfBestLocationSecondMate = scoreOfBestLocationSecondMate;
			this.scoreOfSecondBestLocationSecondMate = scoreOfSecondBestLocationSecondMate;
			
			
			this.nextSegmentUnmapped = nextSegmentUnmapped;
		}
		
		public boolean isNextSegmentUnmapped() {
			return this.nextSegmentUnmapped;
		}
		
		public long getPointerToBestLocationFirstMate() {
			return this.pointerToBestLocationFirstMate;
		}
		
		public long getPointerToBestLocationSecondMate() {
			return this.pointerToBestLocationSecondMate;
		}
		
		public long getPointerToSecondBestLocationFirstMate() {
			return this.pointerToSecondBestLocationFirstMate;
		}
		
		public long getPointerToSecondBestLocationSecondMate() {
			return this.pointerToSecondBestLocationSecondMate;
		}
		
		public double getScoreOfBestLocationFirstMate() {
			return this.scoreOfBestLocationFirstMate;
		}
		
		public double getScoreOfBestLocationSecondMate() {
			return this.scoreOfBestLocationSecondMate;
		}
		
		public double getScoreOfSecondBestLocationFirstMate() {
			return this.scoreOfSecondBestLocationFirstMate;
		}
		
		public double getScoreOfSecondBestLocationSecondMate() {
			return this.scoreOfSecondBestLocationSecondMate;
		}
		
		
	}
	
	private class LocationContainerComparator implements Comparator<LocationContainer> {

		@Override
		public int compare(LocationContainer o1, LocationContainer o2) {
			return Long.valueOf(o1.getPointerToBestLocationFirstMate()).compareTo(o2.getPointerToBestLocationFirstMate());
		}
		
	}
	
	private String trimDouble(double inValue){
		
		twoDec.setGroupingUsed(false);
		return twoDec.format(inValue);
		}
	
	
	public void setReadScores(ArrayList<ReadPair<SparseRead,SparseRead>> reads) {
		double currentDifference;
		ArrayList<Pair<Integer,Integer>> validPairs;
		double scoreA;
		double scoreB;
		for(ReadPair<SparseRead,SparseRead> readPair : reads) {
			
			if(readPair.getSecond() != null && !readPair.getValidPairs().isEmpty()) {
				validPairs = readPair.getValidPairs();
				for(Pair<Integer,Integer> validPair : validPairs) {
					scoreA = readPair.getFirst().getLocations().get(validPair.getFirst()).getScore();
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
			}
		}
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
		}
	}
	
	
	private void updateReadScores(ArrayList<ReadPair<SparseRead,SparseRead>> reads, MaxPriorityQueue queue) {
		double currentDifference;
		double oldScore;
		double scoreA;
		double scoreB;
		ArrayList<ReadPair<SparseRead,SparseRead>> readsWithIncreasedKeys = new ArrayList<ReadPair<SparseRead,SparseRead>>();
		ArrayList<ReadPair<SparseRead,SparseRead>> readsWithDecreasedKeys = new ArrayList<ReadPair<SparseRead,SparseRead>>();
		for(ReadPair<SparseRead,SparseRead> readPair : reads) {
			if(readPair == null) {
				continue;
			}
			
			oldScore = readPair.getScore();
			
			if(readPair.getSecond() != null && !readPair.getValidPairs().isEmpty()) {
				for(Pair<Integer,Integer> validPair : readPair.getValidPairs()) {
					scoreA = readPair.getFirst().getLocations().get(validPair.getFirst()).getScore();
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
			}
			
			if(currentDifference > oldScore) {
				readsWithIncreasedKeys.add(readPair);
			}
			
			else if(currentDifference < oldScore) {
				readsWithDecreasedKeys.add(readPair);
			}
		}
		
		//now update
		for(ReadPair<SparseRead,SparseRead> readPair : readsWithIncreasedKeys) {
			queue.increaseKey(readPair);
		}
		for(ReadPair<SparseRead,SparseRead> readPair : readsWithDecreasedKeys) {
			queue.decreaseKey(readPair);
		}
		readsWithIncreasedKeys.clear();
		readsWithDecreasedKeys.clear();
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
	
	private class PairScoreComparator implements Comparator<Pair> {
		public int compare(Pair r1, Pair r2) {
			double score1 = r1.getScore();
			double score2 = r2.getScore();
			return ((Double)score1).compareTo(score2);
		}
		
	}
	
	
	private void processGlobalContext(HashMap<String,BufferedWriter> chr2writer) throws Exception {
		
		BufferedRandomAccessFile br = new BufferedRandomAccessFile(new File(this.multiMappingFilePath),"r",10240);
		BufferedRandomAccessFile brForMultiMappings = new BufferedRandomAccessFile(new File(this.multiMappingFilePath),"r",10240);
		
		ArrayList<LocationContainer> locations2write = new ArrayList<LocationContainer>();
		ArrayList<LocationContainer> multiLocations2write = new ArrayList<LocationContainer>();
		
		long prevFilePointer = br.getFilePointer();
		String currentLine = br.getNextLine();
		if(currentLine == null) {
			br.close();
			brForMultiMappings.close();
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
		String readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
		String prevReadId = readId;
		String prevReadIdPrefix = readIdPrefix;
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
		MultiSparseReadLocation tmpSparseLocation;
		MultiSparseReadLocation currentSparseLocation = new MultiSparseReadLocation(coordinates, strand);
		currentSparseLocation.setFilePointer(prevFilePointer);
		currentSparseLocation.setScore(locationScore);
		currentSparseReadLocations.add(currentSparseLocation);
		
		SparseRead sparseRead = new SparseRead();
		ReadPair<SparseRead,SparseRead> currentPair = new ReadPair<SparseRead,SparseRead>(readIdPrefix);
		ReadPair<SparseRead,SparseRead> tmpPair = null;
		currentPair.setFirst(sparseRead);
		
		ArrayList<String> chrA = new ArrayList<String>();
		ArrayList<String> chrB = new ArrayList<String>();
		Pair<ArrayList<String>,ArrayList<String>> chrPair = new Pair<ArrayList<String>,ArrayList<String>>();
		chrPair.setFirst(chrA);
		chrPair.setSecond(chrB);
		ArrayList<String> currentChrList = chrA;
		currentChrList.add(chr);
		
		while((currentLine = br.getNextLine()) != null) {
			st = new StringTokenizer(currentLine,"\t");
			currentContextId = st.nextToken();
			
			readId = st.nextToken();
			readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
			
			chr = st.nextToken();
			strand = st.nextToken().toCharArray()[0];
			splittedCoordinates = commaPattern.split(st.nextToken());
			coordinates = new ArrayList<Pair<Integer,Integer>>();
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
			
			//parse location score
			locationScore = Double.valueOf(st.nextToken());
			
			
			
			if(readId.equals(prevReadId)) {
				currentSparseLocation = new MultiSparseReadLocation(coordinates, strand);
				currentSparseReadLocations.add(currentSparseLocation);
				currentSparseLocation.setFilePointer(filePointer);
				currentSparseLocation.setScore(locationScore);
				currentChrList.add(chr);
				
				
			}
			
			//second read from the actual fragment
			else if(readIdPrefix.equals(prevReadIdPrefix)) {
				//the actual sparse read is readA of the current pair
				for(int i = 0; i < currentSparseReadLocations.size(); i++) {
					currentSparseReadLocations.get(i).setAmbiguityFactor(currentSparseReadLocations.size());
					sparseRead.addLocation(currentSparseReadLocations.get(i));
					overallAddedLocations++;
				}
				
				sparseRead = new SparseRead();
				currentPair.setSecond(sparseRead);
				currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
				currentSparseLocation.setFilePointer(filePointer);
				currentSparseLocation.setScore(locationScore);
				currentSparseReadLocations.clear();
				currentSparseReadLocations.add(currentSparseLocation);
				currentChrList = chrB;
				currentChrList.add(chr);
				
				prevReadId = readId;
				
				
			}
			
			//here we have a new fragment
			else {
				
				//first process the old fragment/read
				for(int i = 0; i < currentSparseReadLocations.size(); i++) {
					currentSparseReadLocations.get(i).setAmbiguityFactor(currentSparseReadLocations.size());
					sparseRead.addLocation(currentSparseReadLocations.get(i));
					overallAddedLocations++;
				}
				
				
				//both mates found. check for valid pairs
				if(currentPair.getSecond() != null) {
					
					getValidReadPairs(currentPair,chrPair);
					
					//only found one valid pair. no need for further calculations, print mates.
					if(currentPair.getValidPairs().size() == 1) {
						printLocationToSamFormatFromBufferedReader(currentPair.getFirst().getLocations().get(currentPair.getValidPairs().get(0).getFirst()).getFilePointer(),-1,currentPair.getSecond().getLocations().get(currentPair.getValidPairs().get(0).getSecond()).getFilePointer(),-1,-1,-1,-1,-1,br,false,chr2writer,false);
					}
					
					//at least two valid pairs found. add mates to the multi mapping fraction 
					else if(currentPair.getValidPairs().size() > 1) {
						//this.reads.add(currentPair);
						setReadScore(currentPair);
						processMultiMappedReads(currentPair,locations2write,multiLocations2write,brForMultiMappings,chr2writer);
						
					}
				}

				//found only one mate, but with several mapping positions.
				else if(currentPair.getFirst().getLocations().size() > 1) {
					//reads.add(currentPair);
					setReadScore(currentPair);
					processMultiMappedReads(currentPair,locations2write,multiLocations2write,brForMultiMappings,chr2writer);
				}
				
				else {
					printLocationToSamFormatFromBufferedReader(currentPair.getFirst().getLocations().get(0).getFilePointer(),-1,-1,-1,-1,-1,-1,-1,br,false,chr2writer,true);
				}
				
				overallParsedReads++;
				prevReadId = readId;
				prevReadIdPrefix = readIdPrefix;
				
				sparseRead = new SparseRead();
				currentPair = new ReadPair<SparseRead,SparseRead>(readIdPrefix);
				currentPair.setFirst(sparseRead);
				
				
				currentSparseLocation = new MultiSparseReadLocation(coordinates,strand);
				currentSparseLocation.setFilePointer(filePointer);
				currentSparseLocation.setScore(locationScore);
				currentSparseReadLocations.clear();
				currentSparseReadLocations.add(currentSparseLocation);
				
				chrA.clear();
				chrB.clear();
				currentChrList = chrA;
				currentChrList.add(chr);
				
				
			}
			
			filePointer = br.getFilePointer();
		}
		
		//process the last bunch
		for(int i = 0; i < currentSparseReadLocations.size(); i++) {
			currentSparseReadLocations.get(i).setAmbiguityFactor(currentSparseReadLocations.size());
			sparseRead.addLocation(currentSparseReadLocations.get(i));
			overallAddedLocations++;
		}
		
		
		
		//both mates found. check for valid pairs
		if(currentPair.getSecond() != null) {
			
			getValidReadPairs(currentPair,chrPair);

			//only found one valid pair. no need for further calculations, print mates.
			if(currentPair.getValidPairs().size() == 1) {
				printLocationToSamFormatFromBufferedReader(currentPair.getFirst().getLocations().get(currentPair.getValidPairs().get(0).getFirst()).getFilePointer(),-1,currentPair.getSecond().getLocations().get(currentPair.getValidPairs().get(0).getSecond()).getFilePointer(),-1,-1,-1,-1,-1,br,false,chr2writer,false);
			}
			
			//at least two valid pairs found. add mates to the multi mapping fraction 
			else if(currentPair.getValidPairs().size() > 1) {
				//this.reads.add(currentPair);
				setReadScore(currentPair);
				processMultiMappedReads(currentPair,locations2write,multiLocations2write,brForMultiMappings,chr2writer);
			}
			
		}

		//found only one mate, but with several mapping positions.
		else if(currentPair.getFirst().getLocations().size() > 1) {
			//reads.add(currentPair);
			setReadScore(currentPair);
			processMultiMappedReads(currentPair,locations2write,multiLocations2write,brForMultiMappings,chr2writer);
		}
		
		else {
			printLocationToSamFormatFromBufferedReader(currentPair.getFirst().getLocations().get(0).getFilePointer(),-1,-1,-1,-1,-1,-1,-1,br,false,chr2writer,true);
		}
		
		
		//write last bunch of reads
		Collections.sort(locations2write, this.locationContainerComparator);
		for(LocationContainer container : locations2write) {
			printLocationToSamFormatFromBufferedReader(container.getPointerToBestLocationFirstMate(),container.getPointerToSecondBestLocationFirstMate(),container.getPointerToBestLocationSecondMate(),container.getPointerToSecondBestLocationSecondMate(), 
					   container.getScoreOfBestLocationFirstMate(),container.getScoreOfSecondBestLocationFirstMate(),container.getScoreOfBestLocationSecondMate(),container.getScoreOfSecondBestLocationSecondMate(),
					   brForMultiMappings,true,chr2writer,container.isNextSegmentUnmapped());
		}
		
		locations2write.clear();
		br.close();
		brForMultiMappings.close();
		
	}
	
	private void processMultiMappedReads(ReadPair<SparseRead,SparseRead> highScoringReadPair,ArrayList<LocationContainer> locations2write,ArrayList<LocationContainer> multiLocations2write,BufferedRandomAccessFile bufferedReader,HashMap<String,BufferedWriter> chr2writer) throws Exception {
		boolean nextSegmentUnmapped;
		MultiSparseReadLocation topLocationFirstMate;
		MultiSparseReadLocation topLocationSecondMate;
		MultiSparseReadLocation secondBestLocationFirstMate;
		MultiSparseReadLocation secondBestLocationSecondMate;
		//get top location of first mate
		if(highScoringReadPair.getTopScoringPair() != null) {
			topLocationFirstMate = highScoringReadPair.getFirst().getLocations().get(highScoringReadPair.getTopScoringPair().getFirst());
		}
		else {
			topLocationFirstMate = highScoringReadPair.getFirst().getLocations().get(highScoringReadPair.getFirst().getLocations().size() - 1);
		}
		
		//get second best location of first mate	
		secondBestLocationFirstMate = null;
		if(highScoringReadPair.getFirst().getLocations().size() > 1) {
			if(highScoringReadPair.getTopScoringPair() != null && highScoringReadPair.getValidPairs().size() > 1) {
				secondBestLocationFirstMate = highScoringReadPair.getFirst().getLocations().get(highScoringReadPair.getValidPairs().get(highScoringReadPair.getValidPairs().size() - 2).getFirst());
			}
			else {
				secondBestLocationFirstMate = highScoringReadPair.getFirst().getLocations().get(highScoringReadPair.getFirst().getLocations().size() - 2);
			}
		}
		
		//get top location of second mate
		nextSegmentUnmapped = true;
		topLocationSecondMate = null;
		secondBestLocationSecondMate = null;
		if(highScoringReadPair.getSecond() != null) {
			
			nextSegmentUnmapped = false;
			
			if(highScoringReadPair.getTopScoringPair() != null) {
				topLocationSecondMate = highScoringReadPair.getSecond().getLocations().get(highScoringReadPair.getTopScoringPair().getSecond());
			}
			else {
				topLocationSecondMate = highScoringReadPair.getSecond().getLocations().get(highScoringReadPair.getSecond().getLocations().size() - 1);
			}
			
			//get second best location of second mate	
			secondBestLocationSecondMate = null;
			if(highScoringReadPair.getSecond().getLocations().size() > 1) {
				if(highScoringReadPair.getTopScoringPair() != null && highScoringReadPair.getValidPairs().size() > 1) {
					secondBestLocationSecondMate = highScoringReadPair.getSecond().getLocations().get(highScoringReadPair.getValidPairs().get(highScoringReadPair.getValidPairs().size() - 2).getSecond());
				}
				else {
					secondBestLocationSecondMate = highScoringReadPair.getSecond().getLocations().get(highScoringReadPair.getSecond().getLocations().size() - 2);
				}
			}
		}
		
		
		processReadPair(highScoringReadPair.getFirst(),highScoringReadPair.getSecond(),highScoringReadPair.getValidPairs(),topLocationFirstMate,secondBestLocationFirstMate,topLocationSecondMate,secondBestLocationSecondMate,bufferedReader,locations2write,multiLocations2write,chr2writer,nextSegmentUnmapped);
	
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
	
	public void getValidReadPairs(ReadPair<SparseRead,SparseRead> readPair, Pair<ArrayList<String>,ArrayList<String>> chrPair) {
		try {
			ArrayList<MultiSparseReadLocation> locationsA;
			ArrayList<String> chromosomesA;
			ArrayList<MultiSparseReadLocation> locationsB;
			ArrayList<String> chromosomesB;
			locationsA = readPair.getFirst().getLocations();
			chromosomesA = chrPair.getFirst();
			int currentStart;
			int currentEnd;
			int currentDistance;
			int secondLocationEnd;
			int coordinatesIndexBeforeClippingA;
			int coordinatesIndexBeforeClippingB;
			if(readPair.getSecond() != null) {
				locationsB = readPair.getSecond().getLocations();
				chromosomesB = chrPair.getSecond();
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
						
						if(!chromosomesA.get(i).equals(chromosomesB.get(j)))
							continue;
						
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
								//distance criterion fulfilled
								currentStart = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getFirst();
								currentEnd = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getSecond();
/*								if(locationsA.get(i).getCoordinates().size() > 1 && locationsA.get(i).getCoordinates().get(1).getFirst() > 0) {
									currentStart = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size()-1).getFirst();
									currentEnd = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size()-1).getSecond();
								}
*/								
								currentDistance = locationsB.get(j).getCoordinates().get(0).getFirst() - currentEnd;
								
								//distance criterion fulfilled
								if(currentDistance >= 0 && currentDistance <= this.maxContextSize) {
									readPair.addValidPair(new Pair<Integer,Integer>(i,j));
								}
								else if(currentDistance < 0) {
									secondLocationEnd = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getSecond();
/*									if(locationsB.get(j).getCoordinates().size() > 1 && locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size() - 1).getSecond() > 0)
										secondLocationEnd = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size() - 1).getSecond();
*/									
									if(currentEnd < secondLocationEnd && secondLocationEnd - currentEnd <= this.maxContextSize) {
										readPair.addValidPair(new Pair<Integer,Integer>(i,j));
									}
								}
							}
							
							//B starts upstream of A
							else if(locationsA.get(i).getCoordinates().get(0).getFirst() > locationsB.get(j).getCoordinates().get(0).getFirst()) {
								//distance criterion fulfilled
								currentStart = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getFirst();
								currentEnd = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getSecond();
/*								if(locationsB.get(j).getCoordinates().size() > 1 && locationsB.get(j).getCoordinates().get(1).getFirst() > 0) {
									currentStart = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size()-1).getFirst();
									currentEnd = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size()-1).getSecond();
								}
*/								
								currentDistance = locationsA.get(i).getCoordinates().get(0).getFirst() - currentEnd;
								
								if(currentDistance >= 0 && currentDistance <= this.maxContextSize) {
									readPair.addValidPair(new Pair<Integer,Integer>(i,j));
								}
								else if(currentDistance < 0) {
									secondLocationEnd = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getSecond();
/*									if(locationsA.get(i).getCoordinates().size() > 1 && locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size() - 1).getSecond() > 0)
										secondLocationEnd = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size() - 1).getSecond();
*/									
									if(currentEnd < secondLocationEnd && secondLocationEnd - currentEnd <= this.maxContextSize) {
										readPair.addValidPair(new Pair<Integer,Integer>(i,j));
									}
								}
							}
						}
					}
				}
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
		
	
}
