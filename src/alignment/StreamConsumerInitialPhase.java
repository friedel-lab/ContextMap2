package alignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import main.Pair;
import main.Triplet;

public class StreamConsumerInitialPhase extends SamStreamConsumer {
	
	private String outputFilePath;
	private String splitCandidatesOutputDir;
	private String unalignedReadsOutputFilePath;
	private String multiSplitCandidatesOutputDir;
	private int maxMismatches;
	private int seedMismatches;
	private int[] splitSeedSizes;
	private int seedLength;
	private int threads;
	private final int linesPerThread = 25000;
	private final double splitCandidateWindowSizeRate = 0.1;
	private final double splitCandidateMismatchRate = 0.5; 
	
	private boolean filterAlignments;
	private boolean skipSplitDetection;
	private boolean skipMultiSplitDetection;
	private boolean verbose;
	private boolean mdFlagPreprocessed;
	
	

	public StreamConsumerInitialPhase(InputStream inputStream, StreamType streamType,int maxMismatches, int seedMismatches, int[] splitSeedSizes, int seedLength, int threads, String outputFilePath, String unalignedReadsOutputFilePath, String splitCandidatesOutputDir, String multiIntronCandidatesOutputDir, boolean filterAlignments, boolean verbose, boolean mdFlagPreprocessed, boolean skipSplitDetection, boolean skipMultiSplitDetection, boolean useAllThreads) {
	super(inputStream,streamType);
	this.outputFilePath = outputFilePath;
	this.unalignedReadsOutputFilePath = unalignedReadsOutputFilePath;
	this.splitCandidatesOutputDir = splitCandidatesOutputDir;
	this.multiSplitCandidatesOutputDir = multiIntronCandidatesOutputDir;
	this.filterAlignments = filterAlignments;
	this.verbose = verbose;
	this.mdFlagPreprocessed = mdFlagPreprocessed;
	this.skipSplitDetection = skipSplitDetection;
	this.skipMultiSplitDetection = skipMultiSplitDetection;
	this.maxMismatches = maxMismatches;
	this.seedMismatches = seedMismatches;
	this.seedLength = seedLength;
	this.splitSeedSizes = splitSeedSizes;
	this.threads = Math.max(1, threads/3);
	if(useAllThreads)
		this.threads = threads;
	}
	
	@Override
	public void run() {
		try {
			InputStreamReader inputStreamReader = new InputStreamReader(this.inputStream);
			BufferedReader br = new BufferedReader(inputStreamReader, 1024 * 1024 * 10);
			while(!br.ready() && this.processRunning.get()) {
				this.sleep(100);
			}
			
			//in case the stream is ready, we output its content now
			LineWriter lineWriter = null;
			BufferedWriter bw = null;
			BufferedWriter unalignedBw = null;
			ExecutorService lineProcessingExecutor = Executors.newFixedThreadPool(this.threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			ArrayList<String> lines = new ArrayList<String>();
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			String readId;
			
			if(br.ready()) {
				if(this.streamType.equals(StreamType.STDOUT)) {
					bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFilePath),true)),1024 * 1024 * 10);
				}
				
				if(this.unalignedReadsOutputFilePath != null)
					unalignedBw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.unalignedReadsOutputFilePath),true)),1024 * 1024 * 10);
				
				String line = null;
				lineWriter = new LineWriter(splitCandidatesOutputDir,multiSplitCandidatesOutputDir, bw,unalignedBw);
				
				int processedLines = 0;
				Date date;
				while ((line = br.readLine()) != null) {
					if(this.streamType.equals(StreamType.ERROR)) {
						System.err.println(line);
					}
					else {
						lines.add(line);
						if(++processedLines >= this.linesPerThread) {
							
							splittedLine = tabPattern.split(line);
							readId = splittedLine[0];
							while((line = br.readLine()) != null) {
								splittedLine = tabPattern.split(line);
								if(readId.equals(splittedLine[0]))
									lines.add(line);
								else
									break;
							}
							
							if(line == null) {
								break;
							}
							
							futures.add(lineProcessingExecutor.submit(new LineProcessor(lines, lineWriter, maxMismatches,seedLength,seedMismatches, splitSeedSizes,filterAlignments,mdFlagPreprocessed,skipSplitDetection,skipMultiSplitDetection)));
							lines = new ArrayList<String>();
							lines.add(line);
							processedLines = 1;
							
							if(futures.size() >= 8 * this.threads) {
								for(Future future : futures)
									future.get();
								futures.clear();
								date = new Date();
								
							}
						}
						
					}
				}
			}
			br.close();
			if(bw != null) {
				futures.add(lineProcessingExecutor.submit(new LineProcessor(lines, lineWriter, maxMismatches,seedLength, seedMismatches, splitSeedSizes,filterAlignments,mdFlagPreprocessed,skipSplitDetection,skipMultiSplitDetection)));
				for(Future future : futures)
					future.get();
				futures.clear();
				lineProcessingExecutor.shutdown();
				bw.close();
				lineWriter.closeSplitWriters();
				lineWriter.closeMultiSplitWriters();
			}
			
			if(unalignedBw != null)
				unalignedBw.close();
			
			synchronized(this) {
				this.notifyAll();
			}
		}
		
		catch (Exception ioe) {
			ioe.printStackTrace();
			synchronized(this) {
				this.notifyAll();
			}
		}
	}
	
	

	
	/**
	 * already adds full AND partial alignments to the all_reads.rmap file
	 * @author bonfert
	 *
	 */
	
	private class LineProcessor extends Thread {
		
		private ArrayList<String> lines;
		private int maxMismatches;
		private int seedLength;
		private int seedMismatches;
		private int[] splitSeedSizes;
		private boolean filterAlignments;
		private boolean mdFlagPreprocessed;
		private boolean skipSplitDetection;
		private boolean skipMultiSplitDetection;
		
		
		private LineWriter lineWriter;
		
		public LineProcessor(ArrayList<String> lines, LineWriter lineWriter, int maxMismatches, int seedLength, int seedMismatches, int[] splitSeedSizes, boolean filterAlignments, boolean mdFlagPreprocessed, boolean skipSplitDetection, boolean skipMultiSplitDetection) {
			this.lines = lines;
			this.lineWriter = lineWriter;
			
			this.maxMismatches = maxMismatches;
			this.seedLength = seedLength;
			this.seedMismatches = seedMismatches;
			this.splitSeedSizes = splitSeedSizes;
			this.filterAlignments = filterAlignments;
			this.mdFlagPreprocessed = mdFlagPreprocessed;
			this.skipSplitDetection = skipSplitDetection;
			this.skipMultiSplitDetection = skipMultiSplitDetection;
		}
		
		public void run() {
			try {
				StringTokenizer st;
				StringBuilder matches = new StringBuilder();
				Pair<String,ArrayList<Integer>> alignmentTypeAndMismatchCounts = new Pair<String,ArrayList<Integer>>();
				Triplet<String,String,Integer> mismatchInfoReadIdStart = new Triplet<String,String,Integer>();
				ArrayList<Integer> mismatchCounts = new ArrayList<Integer>();
				String alignmentType = "";
				int mismatchCount = -1;
				String chr;
				int start;
				String cigar;
				String sequence;
				char strand;
				int readLength;
				boolean isReverseComplement;
				boolean softClipped;
				boolean softClippedAtTheStart;
				String mdField = null;
				String mismatchInfo;
				int mdFieldPosition = -1;
				String readId;
				String tmpReadId;
				String flag;
				String mappingQuality;
				String mateRefName;
				String mateStart;
				String fragmentLength;
				StringBuilder lineToWrite = new StringBuilder();
				int tmpMismatchCount;
				int fwdSeedMismatches;
				String[] splittedSamLine;
				Pattern tabPattern = Pattern.compile("\t");
				
				String prevSequence = null;
				char prevStrand = '+';
				StringBuilder sb  = new StringBuilder();
				
				for(String line : this.lines) {
					
					if(line.charAt(0) == '@')
						continue;
					
					st = new StringTokenizer(line,"\t");
					readId = st.nextToken();
					isReverseComplement = false;
					if(readId.length() >= 3 && readId.substring(readId.length() - 3).equals("/rc"))
						isReverseComplement = true;
					flag = st.nextToken();
					strand = getStrandFromSamFlag(Integer.valueOf(flag));
					chr = st.nextToken();
					
					start = Integer.valueOf(st.nextToken());
					mappingQuality = st.nextToken();
					cigar = st.nextToken();
					readLength = getReadLength(cigar);
					mateRefName = st.nextToken();
					mateStart = st.nextToken();
					fragmentLength = st.nextToken();
					sequence = st.nextToken();
					
					if(chr.equals("*")) {
						this.lineWriter.writeUnalignedRead(readId, sequence);
						continue;
					}
					
					if(sequence.equals("*") || cigar.contains("H")) {
						sb.setLength(0);
						sb.append(prevSequence);
						if(prevStrand != strand) {
							sb.reverse();
							for(int i = 0; i < sb.length(); i++) {
								sb.setCharAt(i, substitute(sb.charAt(i)));
							}
						}
						
						sequence = sb.toString();
					}
					
					else {
						prevSequence = sequence;
						prevStrand = strand;
					}
					
					
					st.nextToken();
					
					mdFieldPosition = 11;
					while(st.hasMoreTokens()) {
						mdField = st.nextToken();
						if(mdField.length() > 1 && mdField.substring(0,2).equals("MD"))
							break;
						mdFieldPosition++;
					}
					
					st = new StringTokenizer(mdField, ":");
					st.nextToken();
					st.nextToken();
					mismatchInfo = st.nextToken();
					
					/**
					 * clipped read
					 */
					softClipped = false;
					softClippedAtTheStart = false;
					if(cigar.contains("S") || cigar.contains("H")) {
						softClipped = true;
						softClippedAtTheStart = modifyMDFlagReadIdStart(mismatchInfo,readId,start, cigar,mismatchInfoReadIdStart,strand,this.mdFlagPreprocessed);
						mismatchInfo = mismatchInfoReadIdStart.getFirst();
						tmpReadId = mismatchInfoReadIdStart.getSecond();
						start = mismatchInfoReadIdStart.getThird();
						
						
						//in case a read is clipped at the beginning as well as at the end with more than minAnchorSize bps (set to 10 per default) 
						//we will consider this read as a multi intron spanning read
						if(tmpReadId == null) {
							
							if(this.skipMultiSplitDetection)
								continue;
							
							lineToWrite.setLength(0);
							if(this.mdFlagPreprocessed)
								lineToWrite.append(line);
							else {
								
								splittedSamLine = tabPattern.split(line);
								splittedSamLine[mdFieldPosition] = "MD:Z:" + mismatchInfo;
								splittedSamLine[9] = sequence;
								lineToWrite.append(splittedSamLine[0]);
								for(int i = 1; i < splittedSamLine.length; i++) {
									lineToWrite.append("\t").append(splittedSamLine[i]);
								}
							}
							this.lineWriter.writeMultiSplitAlignmentLine(chr,lineToWrite.toString());
							continue;
						}
					}
					
					
					setAlignmentTypeAndMismatchCounts(mismatchInfo,cigar,alignmentTypeAndMismatchCounts, mismatchCounts,this.maxMismatches,this.seedLength,this.splitSeedSizes,strand,softClipped,softClippedAtTheStart);
					alignmentType = alignmentTypeAndMismatchCounts.getFirst();
					mismatchCounts = alignmentTypeAndMismatchCounts.getSecond();
					
					if(alignmentType.equals("F") || alignmentType.equals("P")) {
						
						/*
						 * in case we filter the alignments (only for the backward alignment step), we discard full
						 * read alignments with <= # seed mismatches in the whole read or in the fwd seed region (such alignments were already
						 * found in the forward alignment step
						 */
						mismatchCount = mismatchCounts.get(mismatchCounts.size() - 1);
						if(this.filterAlignments && alignmentType.equals("F")) {
							if(mismatchCount <= this.seedMismatches)
								continue;
							
							fwdSeedMismatches = 0;
							tmpMismatchCount = mismatchCounts.get(mismatchCounts.size() - 1);
							for(int i = mismatchCounts.size() - 1; i >= mismatchCounts.size() - this.seedLength; i--) {
								if(mismatchCounts.get(i) < tmpMismatchCount) {
									fwdSeedMismatches++;
									tmpMismatchCount--;
								}
							}
							
							if(fwdSeedMismatches <= this.seedMismatches)
								continue;
							
						}
						
						else if(alignmentType.equals("P")) {
							
							if(this.skipSplitDetection)
								continue;
							
							mismatchCount = this.maxMismatches + 1;
						}
						
						
						
						
						if(isReverseComplement) {
							readId = readId.substring(0, readId.length() - 3);
							if(strand == '+')
								strand = '-';
							else
								strand = '+';
						}
						
						lineToWrite.setLength(0);
						lineToWrite.append(readId).append("\t").append(alignmentType).append("\t").append(chr).append("\t").append(start).append("\t.\t").append(strand).append("\t").append(mismatchCount).append("\t").append(readLength);
						this.lineWriter.writeFullAlignmentLine(lineToWrite.toString());
					}
					
					else {
						if(softClipped) {
							readId = mismatchInfoReadIdStart.getSecond();
							if(readId.length() >= 3 && readId.substring(readId.length() - 3).equals("/rc")) {
								if(strand == '+')
									strand = '-';
								else
									strand = '+';
							}
								
						}
						
						lineToWrite.setLength(0);
						lineToWrite.append(readId).append("\t").append(strand).append("\t").append(start).append("\t").append(sequence).append("\t").append("MD:Z:").append(mismatchInfo).append("\t").append(readLength).append("\t").append(alignmentType.split("_")[1]).append("\t").append(cigar);
						this.lineWriter.writeSplitAlignmentLine(chr,lineToWrite.toString());
					}
				}
				this.lines.clear();
		
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		private char substitute(char n) {
			if(n == 'A' || n == 'a') return 'T';
			if(n == 'T' || n == 't') return 'A';
			if(n == 'C' || n == 'c') return 'G';
			if(n == 'G' || n == 'g') return 'C';
			else return 'N';
		}
	}

	
	private class LineWriter {
		
		private HashMap<String, BufferedWriter> chr2splitWriter;
		private HashMap<String, BufferedWriter> chr2multiSplitWriter;
		private BufferedWriter bw;
		private BufferedWriter unalignedBw;
		private String splitCandidatesOutputDir;
		private String multiSplitCandidatesOutputDir;
		
		public LineWriter(String splitCandidatesOutputDir, String multiSplitCandidatesOutputDir, BufferedWriter bw, BufferedWriter unalignedBw) {
			this.splitCandidatesOutputDir = splitCandidatesOutputDir;
			this.multiSplitCandidatesOutputDir = multiSplitCandidatesOutputDir;
			this.chr2splitWriter = new HashMap<String,BufferedWriter>();
			this.chr2multiSplitWriter = new HashMap<String,BufferedWriter>();
			this.bw = bw;
			this.unalignedBw = unalignedBw;
		}
		
		
		
		public synchronized void writeSplitAlignmentLine(String key, String line) {
			try {
				if(!this.chr2splitWriter.containsKey(key))
					chr2splitWriter.put(key, new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.splitCandidatesOutputDir + "/" + key + ".rmap"),true)), 1024 * 512));
				
				this.chr2splitWriter.get(key).write(line);
				this.chr2splitWriter.get(key).newLine();
				
				//if more than 500 filehandles are opened, we will close all writers
				checkWriterCount(this.chr2splitWriter,500);
			}
			catch(Exception e) {
				e.printStackTrace();
				System.err.println(line);
			}
		}
		
		
		public synchronized void writeMultiSplitAlignmentLine(String key, String line) {
			try {
				if(!this.chr2multiSplitWriter.containsKey(key))
					chr2multiSplitWriter.put(key, new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.multiSplitCandidatesOutputDir + "/" + key + ".sam"),true)), 1024 * 512));
				
				this.chr2multiSplitWriter.get(key).write(line);
				this.chr2multiSplitWriter.get(key).newLine();
				
				checkWriterCount(this.chr2multiSplitWriter,500);
			}
			catch(Exception e) {
				e.printStackTrace();
				System.err.println(line);
			}
		}
		
		private void checkWriterCount(HashMap<String,BufferedWriter> chr2writer, int maxFileHandles) throws Exception {
			if(chr2writer.keySet().size() > maxFileHandles) {
				for(BufferedWriter bw : chr2writer.values()) {
					bw.close();
				}
				chr2writer.clear();
			}
		}
		
		public synchronized void writeFullAlignmentLine(String line) {
			try {
				this.bw.write(line);
				this.bw.newLine();
			}
			catch(Exception e) {
				e.printStackTrace();
				System.err.println(line);
			}
		}
		
		public synchronized void writeUnalignedRead(String readId, String sequence) {
			try {
				if(this.unalignedBw != null) {
					this.unalignedBw.write(">" + readId);
					this.unalignedBw.newLine();
					this.unalignedBw.write(sequence);
					this.unalignedBw.newLine();
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		private char substitute(char n) {
			if(n == 'A' || n == 'a') return 'T';
			if(n == 'T' || n == 't') return 'A';
			if(n == 'C' || n == 'c') return 'G';
			if(n == 'G' || n == 'g') return 'C';
			else return 'N';
		}
		
		public void closeSplitWriters() {
			try {
				for(BufferedWriter bw : this.chr2splitWriter.values())
					bw.close();
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		public void closeMultiSplitWriters() {
			try {
				for(BufferedWriter bw : this.chr2multiSplitWriter.values())
					bw.close();
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		
	}
	
	
	private char getStrandFromSamFlag(int samFlag) {
		if(samFlag < 16)
			return '+';
		
		String binaryFlag = Integer.toBinaryString(samFlag);
		if(binaryFlag.charAt(binaryFlag.length() - 5) == '1')
			return '-';
		
		return '+';
			
	}
	
	private boolean isSecondaryAlignment(int samFlag) {
		if(samFlag < 256)
			return false;
		
		String binaryFlag = Integer.toBinaryString(samFlag);
		if(binaryFlag.charAt(binaryFlag.length() - 9) == '1')
			return true;
		
		return false;
			
	}
	
	private int getMDFieldPosition(String samLine) throws Exception {
		String[] splittedLine;
		splittedLine = samLine.split("\t");
		for(int i = 11; i < splittedLine.length; i++) {
			if(splittedLine[i].substring(0,2).equals("MD"))
				return i;
		}
		return -1;
	}
	
	private String modifyReadSequence(String readSequence, String cigar) {
		boolean hardClippedAtTheStart = false;
    	boolean hardClippedAtTheEnd = false;
    	ArrayList<String> clippings = new ArrayList<String>();
    	Pattern pattern = Pattern.compile("[0-9]+[H]");
		Matcher matcher = pattern.matcher(cigar);
    	
		while(matcher.find()) {
			if(matcher.start() == 0) {
				hardClippedAtTheStart = true;
			}
			else
				hardClippedAtTheEnd = true;
			
			clippings.add(matcher.group());
		}
		
		int clippingLengthStart = 0;
		int clippingLengthEnd = 0;
		if(hardClippedAtTheStart) {
			clippingLengthStart = Integer.valueOf(clippings.get(0).substring(0,clippings.get(0).length() - 1));
			
			for(int i = 0; i < clippingLengthStart; i++) {
				readSequence = "N" + readSequence;
			}
		}
		
		if(hardClippedAtTheEnd) {
			clippingLengthEnd = Integer.valueOf(clippings.get(clippings.size() - 1).substring(0,clippings.get(clippings.size() - 1).length() - 1));
			for(int i = 0; i < clippingLengthEnd; i++) {
				readSequence += "N";
			}
		}
		return readSequence;
	}
	
	private boolean modifyMDFlagReadIdStart(String mismatchInfo,String readId, int start, String cigar, Triplet<String,String,Integer> mismatchInfoReadIdStart, char strand, boolean mdFlagPreprocessed) {
		boolean clippedAtTheStart = false;
    	boolean clippedAtTheEnd = false;
    	ArrayList<String> clippings = new ArrayList<String>();
    	Pattern pattern = Pattern.compile("[0-9]+[S|H]");
		Matcher matcher = pattern.matcher(cigar);
    	
		while(matcher.find()) {
			if(matcher.start() == 0) {
				clippedAtTheStart = true;
			}
			else
				clippedAtTheEnd = true;
			
			clippings.add(matcher.group());
		}
		
		int clippingLengthStart = 0;
		int clippingLengthEnd = 0;
		if(clippedAtTheStart) {
			clippingLengthStart = Integer.valueOf(clippings.get(0).substring(0,clippings.get(0).length() - 1));
			start = start - clippingLengthStart;
			
			if(!mdFlagPreprocessed) {
				for(int i = 0; i < clippingLengthStart; i++) {
					mismatchInfo = "0N" + mismatchInfo;
				}
			}
		}
		
		if(clippedAtTheEnd) {
			clippingLengthEnd = Integer.valueOf(clippings.get(clippings.size() - 1).substring(0,clippings.get(clippings.size() - 1).length() - 1));
			
			if(!mdFlagPreprocessed) {
				for(int i = 0; i < clippingLengthEnd; i++) {
					mismatchInfo += "N0";
				}
			}
		}
		
		
		if((strand == '+' && clippedAtTheStart) || (strand == '-' && clippedAtTheEnd)) {
			readId = readId + "/rc";
		}
		
		
		
		if(clippingLengthStart > 1 && clippingLengthEnd > 1)
		//if(clippingLengthStart > 0 && clippingLengthEnd > 0)
			readId = null;
		
		mismatchInfoReadIdStart.setFirst(mismatchInfo);
		mismatchInfoReadIdStart.setSecond(readId);
		mismatchInfoReadIdStart.setThird(start);
		return (clippingLengthStart > clippingLengthEnd);
	}
	
		
	
	private void setAlignmentTypeAndMismatchCounts(String mismatchInfo,String cigar, Pair<String,ArrayList<Integer>> alignmentTypeAndMismatchCounts, ArrayList<Integer> mismatchCounts,int maxAllowedMismatches,int seedLength, int[] splitSeedSizes,char strand, boolean softClipped, boolean softClippedAtTheStart) {
		try {
			mismatchCounts.clear();
			Pattern pattern = Pattern.compile("[0-9]+|[A-Z,#]+");
			String mismatchInfoMod = mismatchInfo.replaceAll("\\^[A-Z]+", "#");
			Matcher matcher = pattern.matcher(mismatchInfoMod);
			ArrayList<String> chunks = new ArrayList<String>();
	    	while (matcher.find()) {
    	        chunks.add( matcher.group() );
    	    }
	    	
	    	//if(chunks.size() > 1 && ((!softClipped && strand == '-') || (softClipped && softClippedAtTheStart && strand == '+') || (softClipped && !softClippedAtTheStart && strand == '-'))) {
	    	if(chunks.size() > 1 && ((!softClipped && strand == '-') || (softClipped && softClippedAtTheStart))) {
	    		Collections.reverse(chunks);
	    	}
			
	    	
	    	int readLength = 0;
			int mismatchCount = 0;
	    	if(chunks.size() == 1) {
	    		readLength = Integer.valueOf(chunks.get(0));
	    		for(int i = 0; i < readLength; i++)
	    			mismatchCounts.add(mismatchCount);
	    	}
	    	
	    	else {
	    		for(int i = 0; i < chunks.size() - 1; i+=2) {
	    			for(int j = readLength; j < readLength + Integer.valueOf(chunks.get(i)); j++) {
						mismatchCounts.add(mismatchCount);
					}
					
					readLength += Integer.valueOf(chunks.get(i));
					if(!chunks.get(i+1).equals("#")) {
						for(int j = 0; j < chunks.get(i+1).length(); j++) {
						mismatchCount++;	
						mismatchCounts.add(mismatchCount);
						readLength++;
						}
					}
	    		}
				for(int j = readLength; j < readLength + Integer.valueOf(chunks.get(chunks.size() - 1)); j++) {
					mismatchCounts.add(mismatchCount);
				}
				readLength += Integer.valueOf(chunks.get(chunks.size() - 1));
	    	}
	    	
			alignmentTypeAndMismatchCounts.setSecond(mismatchCounts);
			if(mismatchCount <= maxAllowedMismatches) {
				alignmentTypeAndMismatchCounts.setFirst("F");
				return;
			}
			
			int lastAllowedMismatchPosition = 0;
			for(int i = 0; i < mismatchCounts.size(); i++) {
				if(mismatchCounts.get(i) <= maxAllowedMismatches)
					lastAllowedMismatchPosition = i;
				else
					break;
			}
			int windowLength = (int) ((double)(readLength) * this.splitCandidateWindowSizeRate);
			if(windowLength < 4)
				windowLength = 4;
			
			if(windowLength > 20) {
				windowLength = 20;
			}
			
			
			if(readLength - lastAllowedMismatchPosition < splitSeedSizes[0] || mismatchCount < this.splitCandidateMismatchRate * windowLength) {
				alignmentTypeAndMismatchCounts.setFirst("P");
				return;
			}
			
			
			int stepSize = windowLength/3;
			if(stepSize == 0)
				stepSize = 1;
			
			if(stepSize > 5)
				stepSize = 5;
			
			boolean foundGap = false;
			int windowOffset = 0;
			for(int i = seedLength; i < mismatchCounts.size() - windowLength; i+= stepSize) {
				if(mismatchCounts.get(i + windowLength - 1) - mismatchCounts.get(i-1) > this.splitCandidateMismatchRate * windowLength) {
					foundGap = true;
					windowOffset = readLength - i;
					break;
				}
			}
			//check the last window
			if(!foundGap && (mismatchCounts.get(mismatchCounts.size()-1) - mismatchCounts.get(mismatchCounts.size() - windowLength - 1) > this.splitCandidateMismatchRate * windowLength)) {
				foundGap = true;
				windowOffset = windowLength;
			}
			
			if(foundGap) {
				for(int i = splitSeedSizes.length - 1; i >= 0; i--) {
					if(windowOffset >= splitSeedSizes[i]) {
						alignmentTypeAndMismatchCounts.setFirst(String.format("S_%s",i));
						return;
					}
				}
			}
			
			
			alignmentTypeAndMismatchCounts.setFirst("P");
			
		}
		catch(Exception e) {
			e.printStackTrace();
			System.err.println(mismatchInfo);
		}
		
	}
	
	
	private int getReadLength(String cigar) {
    	Pattern pattern = Pattern.compile("[0-9]+[M|I|X|=|S|H]");
		Matcher matcher = pattern.matcher(cigar);
    	int readLength = 0;
    	String currentMatch;
		while(matcher.find()) {
			currentMatch = matcher.group();
			readLength += Integer.valueOf(currentMatch.substring(0,currentMatch.length() - 1));
		}
		return readLength;
	}
}
