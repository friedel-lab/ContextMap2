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

public class StreamConsumerSlidingWindow extends SamStreamConsumer {
	
	private String outputFilePath;
	private int maxMismatches;


	private int threads;
	private final int linesPerThread = 100000;
	
	private boolean filterAlignments;
	private boolean mdFlagPreprocessed;
	private boolean verbose;
	

	public StreamConsumerSlidingWindow(InputStream inputStream, StreamType error,int maxMismatches, int threads, String outputFilePath, boolean filterAlignments, boolean verbose, boolean mdFlagPreprocessed) {
		super(inputStream,error);
		this.outputFilePath = outputFilePath;
		this.filterAlignments = filterAlignments;
		this.verbose = verbose;
		this.mdFlagPreprocessed = mdFlagPreprocessed;
		this.maxMismatches = maxMismatches;
		this.threads = Math.max(1, threads/3);
	}
	
	@Override
	public void run() {
		try {
			InputStreamReader inputStreamReader = new InputStreamReader(this.inputStream);
			BufferedReader br = new BufferedReader(inputStreamReader, 1024 * 1024 * 100);
			while(!br.ready() && this.processRunning.get()) {
				this.sleep(100);
			}
			
			//in case the stream is ready, we output its content now
			BufferedWriter bw = null;
			
			ArrayList<String> lines = new ArrayList<String>();
			if(br.ready()) {
				if(this.streamType.equals(StreamType.STDOUT)) {
					bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFilePath),true)),1024 * 1024 * 10);
				}
				
				String tmpLine = null;
				int processedLines = 0;
				while ((tmpLine = br.readLine()) != null) {
					if(this.streamType.equals(StreamType.ERROR)) {
						System.err.println(tmpLine);
					}
					else {
						lines.add(tmpLine);
						if(++processedLines >= this.linesPerThread) {
							writeLines(lines,bw);
							lines.clear();
							processedLines = 0;
						}
						
					}
				}
			}
			br.close();
			if(bw != null) {
				// write last lines here
				writeLines(lines,bw);
				lines.clear();
				bw.flush();
				bw.close();
			}
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
	
	
	
	private void writeLines(ArrayList<String> lines, BufferedWriter bw) throws Exception {
		StringTokenizer st;
		Pair<ArrayList<Integer>,Integer> mismatchPositionsAndReadLength = new Pair<ArrayList<Integer>,Integer>();
		ArrayList<Integer> mismatchPositions = new ArrayList<Integer>();
		int readLength = -1;
		String alignmentType = "";
		int mismatchCount = -1;
		String chr;
		int start;
		String cigar;		
		char strand;
		String mdField = null;
		String mismatchInfo = null;
		int mdFieldPosition = -1;
		String readId;
		StringBuilder lineToWrite = new StringBuilder();
		boolean writeLine = true;
		int currentMismatchCount;
		int currentMismatchPosition;
		Pair<String,Integer> mdFlagAndStart = new Pair<String,Integer>();
		String sequence;
		for(String line : lines) {
			
			if(line.charAt(0) == '@')
				continue;
			
			if(mdFieldPosition == -1) {
				mdFieldPosition = getMDFieldPosition(line);
			}
			st = new StringTokenizer(line,"\t");
			readId = st.nextToken();
			strand = getStrandFromSamFlag(Integer.valueOf(st.nextToken()));
			//skip unaligned read lines
			chr = st.nextToken();
			if(chr.equals("*"))
				continue;
			start = Integer.valueOf(st.nextToken());
			
			st.nextToken();
			cigar = st.nextToken();
			st.nextToken();
			st.nextToken();
			st.nextToken();
			sequence = st.nextToken();
			
			//for(int i = 10; i < mdFieldPosition; i++)
			//	st.nextToken();
			//mdField = st.nextToken();
			while(st.hasMoreTokens()) {
				mdField = st.nextToken();
				if(mdField.length() > 1 && mdField.substring(0,2).equals("MD"))
					break;
			}
			
			st = new StringTokenizer(mdField, ":");
			st.nextToken();
			st.nextToken();
			mismatchInfo = st.nextToken();
			if(cigar.contains("S") || cigar.contains("H")) {
				modifyMDFlag(mismatchInfo,start,mdFlagAndStart,cigar,strand,this.mdFlagPreprocessed);
				mismatchInfo = mdFlagAndStart.getFirst();
				start = mdFlagAndStart.getSecond();
			}
			
			if(writeLine) {
				lineToWrite.setLength(0);
				lineToWrite.append(readId).append("\t").append(strand).append("\t").append(start).append("\t").append("MD:Z:").append(mismatchInfo).append("\t").append(readLength);
				bw.write(lineToWrite.toString());
				bw.newLine();
			}
		
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
	
	private void modifyMDFlag(String mismatchInfo,int start, Pair<String,Integer> mdFlagAndStart, String cigar, char strand,boolean mdFlagPreprocessed) {
		boolean softClippedAtTheStart = false;
    	boolean softClippedAtTheEnd = false;
    	ArrayList<String> clippings = new ArrayList<String>();
    	Pattern pattern = Pattern.compile("[0-9]+[S|H]");
		Matcher matcher = pattern.matcher(cigar);
    	
		while(matcher.find()) {
			if(matcher.start() == 0) {
				softClippedAtTheStart = true;
			}
			else
				softClippedAtTheEnd = true;
			
			clippings.add(matcher.group());
		}
		
	/*	if(strand == '-') {
			boolean tmpBool = softClippedAtTheEnd;
			softClippedAtTheEnd = softClippedAtTheStart;
			softClippedAtTheStart = tmpBool;
			Collections.reverse(clippings);
		}
    */
		
		int clippingLength = 0;
		if(softClippedAtTheStart) {
			clippingLength = Integer.valueOf(clippings.get(0).substring(0,clippings.get(0).length() - 1));
			start = start - clippingLength;
			
			if(!mdFlagPreprocessed) {
				for(int i = 0; i < clippingLength; i++) {
					mismatchInfo = "0N" + mismatchInfo;
				}
			}
		}
		
		if(softClippedAtTheEnd) {
			clippingLength = Integer.valueOf(clippings.get(clippings.size() - 1).substring(0,clippings.get(clippings.size() - 1).length() - 1));
			
			if(!mdFlagPreprocessed) {
				for(int i = 0; i < clippingLength; i++) {
					mismatchInfo += "N0";
				}
			}
		}
		
		mdFlagAndStart.setFirst(mismatchInfo);
		mdFlagAndStart.setSecond(start);
	}
	
	
	private void setMismatchPositionsAndReadLength(String mismatchInfo,Pair<ArrayList<Integer>,Integer> mismatchPositionsAndReadLength, ArrayList<Integer> mismatchPositions) {
		Pattern pattern = Pattern.compile("[0-9]+|[A-Z,#]+");
		mismatchInfo = mismatchInfo.replaceAll("\\^[A-Z]+", "#");
		Matcher matcher = pattern.matcher(mismatchInfo);
		ArrayList<String> chunks = new ArrayList<String>();
    	while (matcher.find()) {
	        chunks.add( matcher.group());
	    }
    	
    	mismatchPositions.clear();
		int readLength = 0;
		if(chunks.size() == 1) {
			readLength = Integer.valueOf(chunks.get(0));
		}
		
		else {
	    	for(int i = 0; i < chunks.size() - 1; i+=2) {
	    		readLength += Integer.valueOf(chunks.get(i));
	    		if(!chunks.get(i+1).equals("#")) {
	    			for(int j = 0; j < chunks.get(i+1).length(); j++) {
	    				mismatchPositions.add(readLength);
	    				readLength++;
	    			}
	    		}
	    		
	    	}
	    	readLength += Integer.valueOf(chunks.get(chunks.size()-1));
		}
		
		mismatchPositionsAndReadLength.setFirst(mismatchPositions);
		mismatchPositionsAndReadLength.setSecond(readLength);
	}
}
