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

public class StreamConsumerMissingFields extends SamStreamConsumer {
	
	private String outputDirPath;
	private String unalignedReadsOutputFilePath;
	
	private int threads;
	private final int linesPerThread = 25000;
	
	private boolean skipMultiSplitDetection;
	private boolean filterAlignments;
	private boolean verbose;
	
	

	public StreamConsumerMissingFields(InputStream inputStream, StreamType streamType, int threads, String outputDirPath, String unalignedReadsOutputFilePath, boolean filterAlignments, boolean skipMultiSplitDetection, boolean verbose) {
	super(inputStream,streamType);
	this.outputDirPath = outputDirPath;
	this.unalignedReadsOutputFilePath = unalignedReadsOutputFilePath;
	this.filterAlignments = filterAlignments;
	this.skipMultiSplitDetection = skipMultiSplitDetection;
	this.verbose = verbose;
	this.threads = Math.max(1, threads/3);
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
			BufferedWriter unalignedBw = null;
			ExecutorService lineProcessingExecutor = Executors.newFixedThreadPool(this.threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			ArrayList<String> lines = new ArrayList<String>();
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			String readId;
			String prevReadId;
			if(br.ready()) {
				if(this.unalignedReadsOutputFilePath != null)
					unalignedBw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.unalignedReadsOutputFilePath),true)),1024 * 1024 * 10);
				
				lineWriter = new LineWriter(this.outputDirPath, unalignedBw);
				String line = null;
				
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
							
							futures.add(lineProcessingExecutor.submit(new LineProcessor(lines, lineWriter,this.skipMultiSplitDetection)));
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
			futures.add(lineProcessingExecutor.submit(new LineProcessor(lines, lineWriter,this.skipMultiSplitDetection)));
			for(Future future : futures)
				future.get();
			futures.clear();
			lineProcessingExecutor.shutdown();
			if(lineWriter != null)
				lineWriter.closeAlignmentWriters();
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
	
	
	private class LineProcessor extends Thread {
		private ArrayList<String> lines;
		private LineWriter lineWriter;
		private boolean skipMultiSplitDetection;
		
		public LineProcessor(ArrayList<String> lines,LineWriter lineWriter,boolean skipMultiSplitDetection) {
			this.lines = lines;
			this.lineWriter = lineWriter;
			this.skipMultiSplitDetection = skipMultiSplitDetection;
		}
		
		public void run() {
			try {
				String[] splittedLine;
				Pattern tabPattern = Pattern.compile("\t");
				String readId;
				int flag;
				String chr;
				String cigar;
				String sequence;
				String prevSequence = null;
				char currentStrand;
				char prevStrand = '+';
				boolean secondaryAlignment = false;
				StringBuilder sb  = new StringBuilder();
				for(String line : this.lines) {
					if(line.charAt(0) == '@')
						continue;
					
					splittedLine = tabPattern.split(line);
					readId = splittedLine[0];
					currentStrand = getStrandFromSamFlag(Integer.valueOf(splittedLine[1]));
					flag = Integer.valueOf(splittedLine[1]);
					
					if((flag & (1L << 6)) != 0) {
						readId += "/1";
					}
					else if((flag & (1L << 7)) != 0) {
						readId += "/2";
					}
					splittedLine[0] = readId;
					
					chr = splittedLine[2];
					cigar = splittedLine[5];
					sequence = splittedLine[9];
					
					if(chr.equals("*")) {
						this.lineWriter.writeUnalignedRead(readId, sequence);
						continue;
					}
					
					if(this.skipMultiSplitDetection && isMultiSplitCandidate(cigar)) {
						if(!sequence.equals("*") && !cigar.contains("H")) {
							prevSequence = sequence;
							prevStrand = getStrandFromSamFlag(Integer.valueOf(splittedLine[1]));
						}
						continue;
					}
					
					if(sequence.equals("*") || cigar.contains("H")) {
						sb.setLength(0);
						sb.append(prevSequence);
						if(prevStrand != currentStrand) {
							sb.reverse();
							for(int i = 0; i < sb.length(); i++) {
								sb.setCharAt(i, substitute(sb.charAt(i)));
							}
						}
						
						splittedLine[9] = sb.toString();
						sb.setLength(0);
						sb.append(splittedLine[0]);
						for(int i = 1; i < splittedLine.length; i++) {
							sb.append("\t").append(splittedLine[i]);
						}
						line = sb.toString();
					}
					
					else {
						
						if(cigar.contains("H")) {
							if(!isSecondaryAlignment(flag)) {
								sequence = modifyReadSequence(sequence,cigar);
							}
						}
						
						prevSequence = sequence;
						prevStrand = getStrandFromSamFlag(Integer.valueOf(splittedLine[1]));
					}
					
					this.lineWriter.writeAlignment(chr,line);
					
				}
				this.lines.clear();
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
	}
	
	
	private boolean isMultiSplitCandidate(String cigar) {
		boolean clippedAtTheStart = false;
    	boolean clippedAtTheEnd = false;
    	Pattern pattern = Pattern.compile("[0-9]+[S|H]");
		Matcher matcher = pattern.matcher(cigar);
    	
		while(matcher.find()) {
			if(matcher.start() == 0) {
				clippedAtTheStart = true;
			}
			else
				clippedAtTheEnd = true;
		}
		return (clippedAtTheStart && clippedAtTheEnd);
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
	
	private boolean isSecondaryAlignment(int samFlag) {
		if(samFlag < 256)
			return false;
		
		String binaryFlag = Integer.toBinaryString(samFlag);
		if(binaryFlag.charAt(binaryFlag.length() - 9) == '1')
			return true;
		
		return false;
			
	}
	
	private char substitute(char n) {
		if(n == 'A' || n == 'a') return 'T';
		if(n == 'T' || n == 't') return 'A';
		if(n == 'C' || n == 'c') return 'G';
		if(n == 'G' || n == 'g') return 'C';
		else return 'N';
	}
	
	
	

	private class LineWriter {
		
		private HashMap<String, BufferedWriter> chr2alignmentWriter;
		private BufferedWriter unalignedBw;
		private String alignmentOutputDir;
		
		
		public LineWriter(String alignmentOutputDir, BufferedWriter unalignedBw) {
			this.alignmentOutputDir = alignmentOutputDir;
			this.chr2alignmentWriter = new HashMap<String,BufferedWriter>();
			this.unalignedBw = unalignedBw;
		}
		
		
		public synchronized void writeAlignment(String chr, String line) {
			try {
				if(!this.chr2alignmentWriter.containsKey(chr))
					this.chr2alignmentWriter.put(chr, new BufferedWriter(new FileWriter(new File(this.alignmentOutputDir + "/" + chr + ".sam"),true)));
				
				this.chr2alignmentWriter.get(chr).write(line);
				this.chr2alignmentWriter.get(chr).newLine();
			}
			catch(Exception e) {
				e.printStackTrace();
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
		
		public void closeAlignmentWriters() {
			try {
				for(BufferedWriter bw : this.chr2alignmentWriter.values())
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
	
	private int getMDFieldPosition(String samLine) throws Exception {
		String[] splittedLine;
		splittedLine = samLine.split("\t");
		for(int i = 11; i < splittedLine.length; i++) {
			if(splittedLine[i].substring(0,2).equals("MD"))
				return i;
		}
		return -1;
	}
	
	private boolean modifyMDFlagReadIdStart(String mismatchInfo,String readId, int start, String cigar, Triplet<String,String,Integer> mismatchInfoReadIdStart, char strand) {
		boolean softClippedAtTheStart = false;
    	boolean softClippedAtTheEnd = false;
    	ArrayList<String> clippings = new ArrayList<String>();
    	Pattern pattern = Pattern.compile("[0-9]+[S]");
		Matcher matcher = pattern.matcher(cigar);
    	
		while(matcher.find()) {
			if(matcher.start() == 0) {
				softClippedAtTheStart = true;
			}
			else
				softClippedAtTheEnd = true;
			
			clippings.add(matcher.group());
		}
		
		int clippingLengthStart = 0;
		int clippingLengthEnd = 0;
		if(softClippedAtTheStart) {
			clippingLengthStart = Integer.valueOf(clippings.get(0).substring(0,clippings.get(0).length() - 1));
			start = start - clippingLengthStart;
			for(int i = 0; i < clippingLengthStart; i++) {
				mismatchInfo = "0N" + mismatchInfo;
			}
		}
		
		if(softClippedAtTheEnd) {
			clippingLengthEnd = Integer.valueOf(clippings.get(clippings.size() - 1).substring(0,clippings.get(clippings.size() - 1).length() - 1));
			
			for(int i = 0; i < clippingLengthEnd; i++) {
				mismatchInfo += "N0";
			}
		}
		
		
		if((strand == '+' && softClippedAtTheStart) || (strand == '-' && softClippedAtTheEnd)) {
			readId = readId + "/rc";
		}
		
		mismatchInfoReadIdStart.setFirst(mismatchInfo);
		mismatchInfoReadIdStart.setSecond(readId);
		mismatchInfoReadIdStart.setThird(start);
		return (clippingLengthStart > clippingLengthEnd);
	}
	
		
	
	private void setAlignmentTypeAndMismatchCounts(String mismatchInfo,String cigar, Pair<String,ArrayList<Integer>> alignmentTypeAndMismatchCounts, ArrayList<Integer> mismatchCounts,int maxAllowedMismatches,int[] splitSeedSizes,char strand, boolean softClipped, boolean softClippedAtTheStart) {
		try {
			mismatchCounts.clear();
			Pattern pattern = Pattern.compile("[0-9]+|[A-Z,#]+");
			String mismatchInfoMod = mismatchInfo.replaceAll("\\^[A-Z]+", "#");
			Matcher matcher = pattern.matcher(mismatchInfoMod);
			ArrayList<String> chunks = new ArrayList<String>();
	    	while (matcher.find()) {
    	        chunks.add( matcher.group() );
    	    }
	    	
	    	if(chunks.size() > 1 && ((!softClipped && strand == '-') || (softClipped && softClippedAtTheStart && strand == '+') || (softClipped && !softClippedAtTheStart && strand == '-'))) {
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
			
			if(readLength - lastAllowedMismatchPosition <= splitSeedSizes[0]) {
				alignmentTypeAndMismatchCounts.setFirst("P");
				
			}
			else {
				alignmentTypeAndMismatchCounts.setFirst("S");
				
			}
	    	
			
		}
		catch(Exception e) {
			e.printStackTrace();
			System.err.println(mismatchInfo);
		}
	}
}
