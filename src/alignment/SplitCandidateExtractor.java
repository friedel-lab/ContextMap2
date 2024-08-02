package alignment;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import main.Pair;
import main.ReadLocation;

import tools.BufferedRandomAccessFile;
import tools.UnixSort;

public class SplitCandidateExtractor implements ActionListener {

	private ReadAligner readAligner;
	private String alignerBinPath;
	private String alignerIndexerPath;
	private String referenceSequencesDir;
	private String tmpOutputDir;
	
	
	private int maxAllowedMismatches;
	private int[] splitSeedSizes;
	private int[] splitSeedMismatches;
	private int maxContextSize;
	private int contextBufferSize;
	
	private int maxGapSize;
	private int minGapSize;
	private int maxDelSize;
	
	private int maxHits;
	private int maxReadLength;
	private int minReadLength;
	private int threads;
	
	private AtomicInteger windowsNeedingBufferedReference;
	private static final double splitCandidateWindowSizeRate = 0.1;
	private static final double splitCandidateMismatchRate = 0.5; 
	
	
	public SplitCandidateExtractor(ReadAligner readAligner, String alignerBinPath, String alignerIndexerPath, String referenceSequencesDir, String tmpOutputDir, int maxAllowedMismatches, int[] splitSeedSizes, int[] splitSeedMismatches, int maxContextSize, int contextBufferSize, int maxGapSize, int minGapSize, int maxDelSize, int maxHits, int maxReadLength, int minReadLength, int threads) {
		this.readAligner = readAligner;
		this.alignerBinPath = alignerBinPath;
		this.alignerIndexerPath = alignerIndexerPath;
		this.referenceSequencesDir = referenceSequencesDir;
		this.maxAllowedMismatches = maxAllowedMismatches;
		this.splitSeedSizes = splitSeedSizes;
		this.splitSeedMismatches = splitSeedMismatches;
		this.maxContextSize = maxContextSize;
		this.contextBufferSize = contextBufferSize;
		
		this.maxGapSize = maxGapSize;
		this.minGapSize = minGapSize;
		this.maxDelSize = maxDelSize;
		
		this.maxHits = maxHits;
		this.maxReadLength = maxReadLength;
		this.minReadLength = minReadLength;
		this.threads = threads;
		this.tmpOutputDir = tmpOutputDir;
		this.windowsNeedingBufferedReference = new AtomicInteger(0);
	}
	
	
	public static void processMultiSplitCandidates(String multiSplitCandidatesDir, int minExonLength, int splitSeedSize, int maxMismatches, int maxHits, String outputDir, String partialOutputFilePath, String sequenceWriterPath, boolean pairedEnd, boolean spaceInHeaderRenamed, boolean hasRenamedPairedEndHeader) {
		try {
			BufferedWriter candidateWriter;
			BufferedWriter partialWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(partialOutputFilePath),true)), 1024 * 1024 * 10);
			PrintWriter sequenceWriter = new PrintWriter(new FileWriter(new File(sequenceWriterPath)));
			BufferedReader br;
			File[] splitCandidateFiles = new File(multiSplitCandidatesDir).listFiles();
			ArrayList<String> lines = new ArrayList<String>();
			String currentLine;
			String tmpLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			StringTokenizer st;
			String currentReadId;
			String prevReadId = "";
			String tmpReadId;
			String tmpReadIdForSequenceWriter;
			String chr;
			char strand;
			char tmpStrand;
			int currentStart;
			int prevStart;
			int tmpStart; 
			int currentEnd;
			int prevEnd;
			String cigar;
			int[] clippingAndMatchBlocks;
			String sequenceForAlignmentLine;
			String sequenceForReadLine;
			StringBuffer sequenceBuffer = new StringBuffer();
			String tmpSequenceForAlignmentLine;
			String tmpSequenceForReadLine;
			int idCounter;
			int readLength;
			int offset;
			String mdField;
			String tmpMdField;
			StringBuffer lineToWrite = new StringBuffer();
			
			
			int maxAllowedMSCs = (maxHits/2) + 1;
			if(maxAllowedMSCs > 100)
				maxAllowedMSCs = 100;
			
			else if(maxAllowedMSCs < 5) {
				maxAllowedMSCs = 5;
			}
			
			
			for(File splitCandidateFile : splitCandidateFiles) {
				candidateWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outputDir + "/" + splitCandidateFile.getName().substring(0, splitCandidateFile.getName().lastIndexOf('.')) + ".rmap"),true)), 1024 * 1024 * 10);
				br = new BufferedReader(new FileReader(splitCandidateFile));
				while((currentLine = br.readLine()) != null) {
					st = new StringTokenizer(currentLine,"\t");
					currentReadId = st.nextToken();
					
					if(!currentReadId.equals(prevReadId)) {
						//process previous read here
						if(!lines.isEmpty() && lines.size() <= maxAllowedMSCs) {
							writeOutCandidates(lines,sequenceBuffer,lineToWrite,minExonLength,splitSeedSize,maxMismatches,prevReadId,sequenceWriter,candidateWriter,partialWriter, pairedEnd,spaceInHeaderRenamed,hasRenamedPairedEndHeader);
						}
						lines.clear();
						prevReadId = currentReadId;
					}
					lines.add(currentLine);
				}
				
				//write out the last found lines
				if(!lines.isEmpty() && lines.size() <= maxAllowedMSCs) {
					writeOutCandidates(lines,sequenceBuffer,lineToWrite,minExonLength,splitSeedSize,maxMismatches,prevReadId,sequenceWriter,candidateWriter, partialWriter, pairedEnd,spaceInHeaderRenamed,hasRenamedPairedEndHeader);
					lines.clear();
				}
				
				br.close();
				candidateWriter.close();
				
			}
			sequenceWriter.close();
			partialWriter.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void writeOutCandidates(ArrayList<String> lines, StringBuffer sequenceBuffer, StringBuffer lineToWrite, int minExonLength, int splitSeedSize, int maxMismatches, String prevReadId, PrintWriter sequenceWriter, BufferedWriter candidateWriter, BufferedWriter partialWriter, boolean pairedEnd, boolean spaceInHeaderRenamed, boolean hasRenamedPairedEndHeader) throws Exception {
		
		Pattern tabPattern = Pattern.compile("\t");
		String[] splittedLine;
		String tmpReadId;
		String tmpReadIdForSequenceWriter;
		char strand;
		char tmpStrand;
		String chr;
		int currentStart;
		String cigar;
		int[] clippingAndMatchBlocks;
		String sequenceForAlignmentLine;
		String sequenceForReadLine;
		String tmpSequenceForAlignmentLine;
		String tmpSequenceForReadLine;
		String mdField;
		String tmpMdField;
		
		int prevStart = -1;
		int tmpStart;
		int currentEnd;
		int prevEnd = -1;
		int idCounter = 0;
		int readLength;
		int offset;
		
		for(int i = 0; i < lines.size(); i++) {
			
			splittedLine = tabPattern.split(lines.get(i));
			strand = getStrandFromSamFlag(Integer.valueOf(splittedLine[1]));
			chr = splittedLine[2];
			currentStart = Integer.valueOf(splittedLine[3]);
			cigar = splittedLine[5];
			clippingAndMatchBlocks = getClippingAndMatchBlocks(cigar);
			currentEnd = currentStart + clippingAndMatchBlocks[1] - 1;
			sequenceForAlignmentLine = splittedLine[9];
			sequenceForReadLine = sequenceForAlignmentLine;
			
			mdField = null;
			for(int j = 10; j < splittedLine.length; j++) {
				mdField = splittedLine[j];
				if(mdField.length() > 1 && mdField.substring(0,2).equals("MD"))
					break;
			}

			//first block
			if(i == 0) {
				offset = clippingAndMatchBlocks[0];
			}
			//center block
			else {
				offset = Math.min(minExonLength, clippingAndMatchBlocks[0]);
			}
			
			tmpStart = currentStart - offset;
			readLength = offset + clippingAndMatchBlocks[1];
			tmpSequenceForAlignmentLine = sequenceForAlignmentLine.substring(clippingAndMatchBlocks[0] - offset,(clippingAndMatchBlocks[0] - offset) + readLength);
			
			
			if(pairedEnd && prevReadId.length() > 2 && (prevReadId.substring(prevReadId.length() - 2).equals("/1") || prevReadId.substring(prevReadId.length() - 2).equals("/2"))) {  
				tmpReadId = prevReadId + "::MSC::" + chr + "::" + (idCounter) + "::" + prevReadId.substring(prevReadId.length() - 2);
				tmpReadIdForSequenceWriter = tmpReadId;
				tmpReadIdForSequenceWriter = tmpReadIdForSequenceWriter.replace("#*#"," ");
			}
			
			else if(pairedEnd && prevReadId.length() > 2 && (prevReadId.substring(prevReadId.length() - 2).equals("#1") || prevReadId.substring(prevReadId.length() - 2).equals("#2"))) {
				tmpReadId = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.substring(prevReadId.length() - 1) + "::MSC::" + chr + "::" + (idCounter) + "::" + prevReadId.substring(prevReadId.length() - 2);
				tmpReadIdForSequenceWriter = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.substring(prevReadId.length() - 1) + "::MSC::" + chr + "::" + (idCounter) + "::/" + prevReadId.substring(prevReadId.length() - 1);
				tmpReadIdForSequenceWriter = tmpReadIdForSequenceWriter.replace("#*#"," ");
			}
			else {
				tmpReadId = prevReadId + "::MSC::" + chr + "::" + idCounter;
				
				if(!spaceInHeaderRenamed && !hasRenamedPairedEndHeader)
					tmpReadIdForSequenceWriter = tmpReadId;
				
				else {
					tmpReadIdForSequenceWriter = tmpReadId;
					if(hasRenamedPairedEndHeader) {
						tmpReadId = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.substring(prevReadId.length() - 1) + "::MSC::" + chr + "::" + (idCounter);
						tmpReadIdForSequenceWriter = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.charAt(prevReadId.length() - 1) + "::MSC::" + chr + "::" + idCounter;
					}
					if(spaceInHeaderRenamed)
						tmpReadIdForSequenceWriter = tmpReadIdForSequenceWriter.replace("#*#"," ");
				}
			}
			
			sequenceBuffer.setLength(0);
			sequenceBuffer.append(tmpSequenceForAlignmentLine);
			if(strand == '-') {
				sequenceBuffer.reverse();
				for(int j = 0; j < sequenceBuffer.length(); j++) {
					sequenceBuffer.setCharAt(j, substitute(sequenceBuffer.charAt(j)));
				}
			}
			sequenceWriter.println(tmpReadIdForSequenceWriter + "\t" + sequenceBuffer.toString());
			
			
			if(offset < splitSeedSize) {
				lineToWrite.setLength(0);
				lineToWrite.append(tmpReadId).append("\t").append("P").append("\t").append(chr).append("\t").append(tmpStart).append("\t.\t").append(strand).append("\t").append(maxMismatches + 1).append("\t").append(readLength);
				partialWriter.write(lineToWrite.toString());
				partialWriter.newLine();
			}
			
			
			else {
				tmpStrand = strand;
				if(strand == '+') {
					tmpReadId = tmpReadId + "/rc";
					tmpStrand = '-';
				}
				
				cigar = offset + "S" + clippingAndMatchBlocks[1] + "M";
				tmpMdField = modifyMdTag(mdField.split(":")[2],clippingAndMatchBlocks[0] - offset + 1,(clippingAndMatchBlocks[0] - offset) + readLength);
				
				
				lineToWrite.setLength(0);
				lineToWrite.append(tmpReadId).append("\t").append(tmpStrand).append("\t").append(tmpStart).append("\t").append(tmpSequenceForAlignmentLine).append("\t").append("MD:Z:").append(tmpMdField).append("\t").append(readLength).append("\t").append(0).append("\t").append(cigar);
				candidateWriter.write(lineToWrite.toString());
				candidateWriter.newLine();
													
			}
			idCounter++;

						
			if(i == lines.size() - 1) {
				offset = clippingAndMatchBlocks[2];
			}
			else
				offset = Math.min(minExonLength, clippingAndMatchBlocks[2]);
			
				
			tmpStart = currentStart;
			readLength = clippingAndMatchBlocks[1] + offset;
			tmpSequenceForAlignmentLine = sequenceForAlignmentLine.substring(clippingAndMatchBlocks[0], clippingAndMatchBlocks[0] + readLength);
			
			if(pairedEnd && prevReadId.length() > 2 && (prevReadId.substring(prevReadId.length() - 2).equals("/1") || prevReadId.substring(prevReadId.length() - 2).equals("/2"))) {
				//tmpReadId = prevReadId.substring(0,prevReadId.length() - 2) + "::MSC::" + chr + "::" + (idCounter) + "::" + prevReadId.substring(prevReadId.length() - 2);
				tmpReadId = prevReadId + "::MSC::" + chr + "::" + (idCounter) + "::" + prevReadId.substring(prevReadId.length() - 2);
				tmpReadIdForSequenceWriter = tmpReadId;
				tmpReadIdForSequenceWriter = tmpReadIdForSequenceWriter.replace("#*#"," ");
			} 
			
			else if(pairedEnd && prevReadId.length() > 2 && (prevReadId.substring(prevReadId.length() - 2).equals("#1") || prevReadId.substring(prevReadId.length() - 2).equals("#2"))) {
				//tmpReadId = prevReadId.substring(0,prevReadId.length() - 2) + "::MSC::" + chr + "::" + (idCounter) + "::" + prevReadId.substring(prevReadId.length() - 2);
				tmpReadId = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.substring(prevReadId.length() - 1) + "::MSC::" + chr + "::" + (idCounter) + "::" + prevReadId.substring(prevReadId.length() - 2);
				//tmpReadIdForSequenceWriter = prevReadId.substring(0,prevReadId.length() - 2) + "::MSC::" + chr + "::" + (idCounter) + "::/" + prevReadId.substring(prevReadId.length() - 1);
				tmpReadIdForSequenceWriter = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.substring(prevReadId.length() - 1) + "::MSC::" + chr + "::" + (idCounter) + "::/" + prevReadId.substring(prevReadId.length() - 1);
				tmpReadIdForSequenceWriter = tmpReadIdForSequenceWriter.replace("#*#"," ");
				
			}
			else {
				tmpReadId = prevReadId + "::MSC::" + chr + "::" + idCounter;
				
				if(!spaceInHeaderRenamed && !hasRenamedPairedEndHeader)
					tmpReadIdForSequenceWriter = tmpReadId;
				
				else {
					tmpReadIdForSequenceWriter = tmpReadId;
					if(hasRenamedPairedEndHeader) {
						tmpReadId = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.substring(prevReadId.length() - 1) + "::MSC::" + chr + "::" + (idCounter);
						tmpReadIdForSequenceWriter = prevReadId.substring(0,prevReadId.length() - 2) + "/" + prevReadId.charAt(prevReadId.length() - 1) + "::MSC::" + chr + "::" + idCounter;
					}
					if(spaceInHeaderRenamed)
						tmpReadIdForSequenceWriter = tmpReadIdForSequenceWriter.replace("#*#"," ");
				}
				
			}

			sequenceBuffer.setLength(0);
			sequenceBuffer.append(tmpSequenceForAlignmentLine);
			if(strand == '-') {
				sequenceBuffer.reverse();
				for(int j = 0; j < sequenceBuffer.length(); j++) {
					sequenceBuffer.setCharAt(j, substitute(sequenceBuffer.charAt(j)));
				}
			}
			sequenceWriter.println(tmpReadIdForSequenceWriter + "\t" + sequenceBuffer.toString());
			
			
			if(offset < splitSeedSize) {
				lineToWrite.setLength(0);
				lineToWrite.append(tmpReadId).append("\t").append("P").append("\t").append(chr).append("\t").append(tmpStart).append("\t.\t").append(strand).append("\t").append(maxMismatches + 1).append("\t").append(readLength);
				partialWriter.write(lineToWrite.toString());
				partialWriter.newLine();
			}
			
			else {
				
				tmpStrand = strand;
				if(strand == '-') {
					tmpReadId = tmpReadId + "/rc";
					tmpStrand = '+';
				}
				
				cigar = clippingAndMatchBlocks[1] + "M" + offset + "S";
				
				tmpMdField = modifyMdTag(mdField.split(":")[2],clippingAndMatchBlocks[0] + 1,clippingAndMatchBlocks[0] + readLength);
				lineToWrite.setLength(0);
				lineToWrite.append(tmpReadId).append("\t").append(tmpStrand).append("\t").append(tmpStart).append("\t").append(tmpSequenceForAlignmentLine).append("\t").append("MD:Z:").append(tmpMdField).append("\t").append(readLength).append("\t").append(0).append("\t").append(cigar);
				
				
				candidateWriter.write(lineToWrite.toString());
				candidateWriter.newLine();
			}
			idCounter++;
			
		}
	}
	
	private static char substitute(char n) {
		if(n == 'A' || n == 'a') return 'T';
		if(n == 'T' || n == 't') return 'A';
		if(n == 'C' || n == 'c') return 'G';
		if(n == 'G' || n == 'g') return 'C';
		else return 'N';
	}
	
	private static String modifyMdTag(String mdTag, int startLength, int stopLength) {
		Pattern pattern = Pattern.compile("[0-9]+|[A-Z,#]+");
		String mismatchInfoMod = mdTag.replaceAll("\\^[A-Z]+", "#");
		Matcher matcher = pattern.matcher(mismatchInfoMod);
		ArrayList<String> chunks = new ArrayList<String>();
    	while (matcher.find()) {
	        chunks.add( matcher.group() );
	    }
    	
    	if(chunks.size() == 1) {
    		return (String.valueOf(stopLength - startLength + 1));
    	}
    	
    	int readLength = 0;
    	int currentMatchLength;
    	String returnValue = "";
    	boolean inInterval = false;
    	int chunkIndex = 0;
    	while(!inInterval && chunkIndex < chunks.size()) {
    		currentMatchLength = Integer.valueOf(chunks.get(chunkIndex));
			if(readLength + currentMatchLength >= startLength) {
				if((readLength + currentMatchLength) - startLength + 1 <= (stopLength - startLength + 1)) 
					returnValue += String.valueOf(((readLength + currentMatchLength) - startLength) + 1);
				
				else
					return (String.valueOf(stopLength - startLength + 1));
				
				inInterval = true;
			}
			
			readLength += currentMatchLength;
			chunkIndex++;
			
			if(chunkIndex == chunks.size()) {
				inInterval = false;
    			break;
			}
			
			if(!chunks.get(chunkIndex).equals("#")) {
				if(readLength + 1 >= startLength) {
					if(readLength + 1 > stopLength) {
						inInterval = false;
						break;
					}
					
					else if(readLength + 1 == stopLength) {
						returnValue += chunks.get(chunkIndex) + "0";
						inInterval = false;
						break;
					}
					
					
					if(readLength + 1 == startLength)
						returnValue += "0";
					
					returnValue += chunks.get(chunkIndex);
					inInterval = true;
				}
				readLength++;
			}
			chunkIndex++;
		}
    	
    	while(inInterval) {
    		currentMatchLength = Integer.valueOf(chunks.get(chunkIndex));
    		if(readLength + currentMatchLength <= stopLength) {
    			returnValue += String.valueOf(currentMatchLength);
    		}
    		else {
    			returnValue += String.valueOf((stopLength - readLength));
    			break;
    		}
    		
    		readLength += currentMatchLength;
    		chunkIndex++;
    		
    		if(chunkIndex == chunks.size())
    			break;
    		
    		if(!chunks.get(chunkIndex).equals("#")) {
				if(readLength + 1 > stopLength) {
					break;
				}
				
				else if(readLength + 1 == stopLength) {
					returnValue += chunks.get(chunkIndex) + "0";
					break;
				}
					
				readLength++;
				returnValue += chunks.get(chunkIndex);
			}
			chunkIndex++;
    	}
    	
    	
    	return returnValue;
	}
	
	private static void setAlignmentTypeAndMismatchCounts(String mismatchInfo, Pair<String,ArrayList<Integer>> alignmentTypeAndMismatchCounts, ArrayList<Integer> mismatchCounts,int maxAllowedMismatches,int seedLength, int[] splitSeedSizes,char strand, boolean softClipped, boolean softClippedAtTheStart) {
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
			int windowLength = (int) ((double)(readLength) * splitCandidateWindowSizeRate);
			if(windowLength < 4)
				windowLength = 4;
			
			
			if(readLength - lastAllowedMismatchPosition < splitSeedSizes[0] || mismatchCount < splitCandidateMismatchRate * windowLength) {
			//if(mismatchCount < this.splitCandidateMismatchRate * windowLength) {
				alignmentTypeAndMismatchCounts.setFirst("P");
				return;
			}
			
			
			int stepSize = windowLength/3;
			if(stepSize == 0)
				stepSize = 1;
			boolean foundGap = false;
			int windowOffset = 0;
			for(int i = seedLength; i < mismatchCounts.size() - windowLength; i+= stepSize) {
				if(mismatchCounts.get(i + windowLength - 1) - mismatchCounts.get(i-1) > splitCandidateMismatchRate * windowLength) {
					foundGap = true;
					windowOffset = readLength - i;
					break;
				}
			}
			//check the last window
			if(!foundGap && (mismatchCounts.get(mismatchCounts.size()-1) - mismatchCounts.get(mismatchCounts.size() - windowLength - 1) > splitCandidateMismatchRate * windowLength)) {
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
	
	private static char getStrandFromSamFlag(int samFlag) {
		if(samFlag < 16)
			return '+';
		
		String binaryFlag = Integer.toBinaryString(samFlag);
		if(binaryFlag.charAt(binaryFlag.length() - 5) == '1')
			return '-';
		
		return '+';
			
	}
	
	private static int[] getClippingAndMatchBlocks(String cigar) {
		int[] blocks = new int[3];
		ArrayList<String> clippings = new ArrayList<String>();
    	Pattern pattern = Pattern.compile("[0-9]+[S|H|M]");
		Matcher matcher = pattern.matcher(cigar);
    	
		int blockIndex = 0;
		String currentBlock;
		while(matcher.find()) {
			currentBlock = matcher.group();
			blocks[blockIndex] = Integer.valueOf(currentBlock.substring(0,currentBlock.length() - 1));
			blockIndex++;
		}
		
		return blocks;
	}
	
	private HashMap<String,String> mapReferenceName2Path(String referenceSequencesDir) throws Exception {
		File[] referenceFiles = new File(referenceSequencesDir).listFiles();
		HashMap<String,String> refName2Path = new HashMap<String,String>();
		for(File reference : referenceFiles) {
			if(reference.isFile() && reference.getName().contains("."))
				refName2Path.put(reference.getName().substring(0,reference.getName().lastIndexOf('.')),reference.getAbsolutePath());
		}
		return refName2Path;
	}
	
	
	private void bufferNextReferenceSequence(String referenceSequencePath, StringBuilder reference) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(referenceSequencePath)));
		String currentLine;
		//skip the header
		br.readLine();
		reference.setLength(0);
		while((currentLine = br.readLine()) != null) {
			reference.append(currentLine);
		}
		br.close();
	}
	
	public void processAlignments(String inputFolderPath, String outFilePath, boolean verbose) {
		try {
			BufferedRandomAccessFile br;
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outFilePath),true)));
			LineWriter splitCandidateWriter = new LineWriter(bw);
			File[] candidateFiles = new File(inputFolderPath).listFiles();
			String chr;
			String currentLine;
			Pattern tabPattern = Pattern.compile("\t");
			String[] splittedLine;
			String readId;
			int contextStart;
			int contextEnd;
			int nextContextStart;
			int nextContextEnd;
			String strand;
			int currentStart;
			int currentEnd;
			long prevFilePointer;
			long currentFilePointer;
			StringBuilder referenceSequence = new StringBuilder();
			HashMap<String,String> refName2Path = mapReferenceName2Path(this.referenceSequencesDir);
			boolean readReverseComplemented;
			ExecutorService windowProcessingExecutor = Executors.newFixedThreadPool(this.threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			String currentOutputDir;
			for(File candidateFile : candidateFiles) {
				chr = candidateFile.getName().substring(0,candidateFile.getName().lastIndexOf('.'));
				br = new BufferedRandomAccessFile(candidateFile,"r",1024 * 1024);
				currentFilePointer = br.getFilePointer();
				prevFilePointer = currentFilePointer;
				
				if(!refName2Path.containsKey(chr)) {
					//print warning here
					br.close();
					continue;
				}
			
				
				bufferNextReferenceSequence(refName2Path.get(chr),referenceSequence);
				contextStart = Integer.MAX_VALUE;
				contextEnd = Integer.MIN_VALUE;
				while((currentLine = br.getNextLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					readId = splittedLine[0];
					strand = splittedLine[1];
					currentStart = Integer.valueOf(splittedLine[2]);
					currentEnd = currentStart + this.contextBufferSize - 1;
					
					readReverseComplemented = false;
					if(readId.length() >= 3 && readId.substring(readId.length() -3).equals("/rc"))
						readReverseComplemented = true;
					
					if(strand.equals("-")) {
						currentEnd = currentStart;
						currentStart = currentStart - this.contextBufferSize + 1;
					}
				
			        currentStart = Math.max(1, currentStart);
			        currentEnd = Math.min(referenceSequence.length(), currentEnd);
			        currentEnd = Math.max(1, currentEnd);
					
					
					if(contextStart == Integer.MAX_VALUE) {
						contextStart = currentStart;
						contextEnd = currentEnd;
					}
					
					
					nextContextStart = (currentStart < contextStart)? currentStart:contextStart;
					nextContextEnd = (currentEnd > contextEnd)? currentEnd:contextEnd;
					
					//process current window and open a new one
					if(nextContextEnd - nextContextStart + 1 >=  this.maxContextSize) {
						currentOutputDir = String.format("%s/%s_%s_%s",this.tmpOutputDir,chr,contextStart,contextEnd);
						SlidingContextProcessor contextProcessor = new SlidingContextProcessor(this.readAligner,candidateFile,splitCandidateWriter, prevFilePointer, currentFilePointer, contextStart, contextEnd, referenceSequence, this.alignerBinPath,this.alignerIndexerPath,currentOutputDir,this.maxAllowedMismatches,this.splitSeedSizes, this.maxHits, this.maxReadLength, this.minReadLength, this.maxGapSize, this.minGapSize, this.maxDelSize, this.threads,verbose);
						contextProcessor.addListener(this);
						futures.add(windowProcessingExecutor.submit(contextProcessor));
						this.windowsNeedingBufferedReference.incrementAndGet();
						contextStart = currentStart;
						contextEnd = currentEnd;
						prevFilePointer = currentFilePointer;
					}
					
					else {
						contextStart = nextContextStart;
						contextEnd = nextContextEnd;
					}
					
					currentFilePointer = br.getFilePointer();
				}
				br.close();
				//adding the last window
				if(contextStart != Integer.MAX_VALUE) {
					currentOutputDir = String.format("%s/%s_%s_%s",this.tmpOutputDir,chr,contextStart,contextEnd);
					SlidingContextProcessor contextProcessor = new SlidingContextProcessor(this.readAligner,candidateFile,splitCandidateWriter, prevFilePointer, currentFilePointer, contextStart, contextEnd, referenceSequence, this.alignerBinPath,this.alignerIndexerPath,currentOutputDir,this.maxAllowedMismatches,this.splitSeedSizes, this.maxHits, this.maxReadLength, this.minReadLength, this.maxGapSize,this.minGapSize, this.maxDelSize, this.threads,verbose);
					contextProcessor.addListener(this);
					futures.add(windowProcessingExecutor.submit(contextProcessor));
					this.windowsNeedingBufferedReference.incrementAndGet();
				}
				
				while(this.windowsNeedingBufferedReference.get() != 0) {
					Thread.sleep(500);
				}
				
			}
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			windowProcessingExecutor.shutdown();
			bw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	public void actionPerformed(ActionEvent event) {
		if(event.getActionCommand().equals("unlocked"))
			this.windowsNeedingBufferedReference.decrementAndGet();
	}
	
	
	private class SlidingContextProcessor extends Thread {
		
		private ReadAligner readAligner;
		private File alignmentFile;
		private long windowStartPointer;
		private long windowEndPointer;
		private int windowStart;
		private int windowEnd;
		private int maxHits;
		private int maxReadLength;
		private int minReadLength;
		
		private int maxGapSize;
		private int minGapSize;
		private int maxDelSize;
		
		private StringBuilder referenceSequence;
		private String alignerBinPath;
		private String alignerIndexerPath;
		private String outputDir;
		
		private LineWriter splitCandidateWriter;
		
		private int maxMismatches;
		private int[] splitSeedSizes;
		private int threads;
		private final int linesToWriteSize = 50000;
		
		private ArrayList<SplitCandidate> candidates;
		private ArrayList<ActionListener> listeners;
		
		private boolean verbose;
		
		public SlidingContextProcessor(ReadAligner readAligner, File alignmentFile, LineWriter splitCandidateWriter, long windowStartPointer, long windowEndPointer, int windowStart, int windowEnd, StringBuilder referenceSequence, String alignerBinPath, String alignerIndexerPath, String outputDir, int maxMismatches, int[] splitSeedSizes, int maxHits, int maxReadLength, int minReadLength, int maxGapSize,int minGapSize, int maxDelSize, int threads, boolean verbose) {
			this.readAligner = readAligner;
			this.alignmentFile = alignmentFile;
			this.splitCandidateWriter = splitCandidateWriter;
			this.windowStartPointer = windowStartPointer;
			this.windowEndPointer = windowEndPointer;
			this.windowStart = windowStart;
			this.windowEnd = windowEnd;
			this.referenceSequence = referenceSequence;
			this.alignerBinPath = alignerBinPath;
			this.alignerIndexerPath = alignerIndexerPath;
			this.outputDir = outputDir;
			
			this.maxMismatches = maxMismatches;
			this.splitSeedSizes = splitSeedSizes;
			this.maxHits = maxHits;
			this.maxReadLength = maxReadLength;
			this.minReadLength = minReadLength;
			
			this.maxGapSize = maxGapSize;
			this.minGapSize = minGapSize;
			this.maxDelSize = maxDelSize;
			
			this.threads = threads;
			
			//TODO check memory consumption for this array. otherwise write to disk, sort by read id and process each read sequentially.
			this.candidates = new ArrayList<SplitCandidate>();
			this.listeners = new ArrayList<ActionListener>();
			
			this.verbose = verbose;
		}
		
		
		public void run() {
			try {
				
				if(!new File(this.outputDir).exists())
					new File(this.outputDir).mkdirs();
				
				// build fasta sequence of current contig
				Date date = new Date();
				long prevTimePoint = System.currentTimeMillis();
				long currentTimePoint;
				double usedTime;
				double overallUsedTime = 0.0;
				String chrName = this.alignmentFile.getName().substring(0,alignmentFile.getName().lastIndexOf('.'));
				
				String contigFilePath = String.format("%s/%s_%s.fa",this.outputDir,this.windowStart - 1,this.windowEnd - 1);
				createContigSequence(contigFilePath,this.referenceSequence,this.windowStart - 1,this.windowEnd - 1);
				fireAction(new ActionEvent(this,0,"unlocked"));
				currentTimePoint = System.currentTimeMillis();
				usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
				overallUsedTime += usedTime;
				if(this.verbose)
					System.out.println(String.format("%s\t Created contig -> %s,%s,%s (took me: %s s)", this.toString(),chrName,this.windowStart,this.windowEnd,usedTime));
				
				
				//generate read sequences and parse the candidate split alignments
				prevTimePoint = System.currentTimeMillis();
				String[] readSequencesFwdFilePaths = new String[splitSeedSizes.length];
				String[] readSequencesRevFilePaths = new String[splitSeedSizes.length];
				for(int i = 0; i < splitSeedSizes.length; i++) {
					readSequencesFwdFilePaths[i] = String.format("%s/fwd_reads_%s_%s_%s_%s.fa",this.outputDir,this.alignmentFile.getName().substring(0,alignmentFile.getName().lastIndexOf('.')),this.windowStart - 1,this.windowEnd - 1,i);
					readSequencesRevFilePaths[i] = String.format("%s/rev_reads_%s_%s_%s_%s.fa",this.outputDir,this.alignmentFile.getName().substring(0,alignmentFile.getName().lastIndexOf('.')),this.windowStart - 1,this.windowEnd - 1,i);
				}
				createReadSequencesAndParseCandidates(readSequencesFwdFilePaths,readSequencesRevFilePaths, this.alignmentFile, this.windowStartPointer, this.windowEndPointer);
				Collections.sort(this.candidates, new SplitCandidateReadIdComparator());
				currentTimePoint = System.currentTimeMillis();
				usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
				overallUsedTime += usedTime;
				if(this.verbose)
					System.out.println(String.format("%s\t Generated read sequences and parsed %s candidates (took me: %s s)",this.toString(), candidates.size(), usedTime));
				
				
				// build index for current contig and align reads
				prevTimePoint = System.currentTimeMillis();
				String bowtieBasePath = String.format("%s/%s_%s",this.outputDir,this.windowStart - 1,this.windowEnd - 1);
				this.readAligner.buildIndex(this.alignerIndexerPath, contigFilePath, bowtieBasePath);
				currentTimePoint = System.currentTimeMillis();
				usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
				overallUsedTime += usedTime;
				if(this.verbose)
					System.out.println(String.format("%s\t Builded index of current contig (took me %s s)",this.toString(),usedTime));
				
				
				prevTimePoint = System.currentTimeMillis();
				String tmpAlignmentFilePath;
				String alignmentFilePath = String.format("%s/all_alignments_%s_%s_%s.fa",this.outputDir,this.alignmentFile.getName().substring(0,alignmentFile.getName().lastIndexOf('.')),this.windowStart - 1,this.windowEnd - 1);
				for(int i = 0; i < this.splitSeedSizes.length; i++) {
					if(new File(readSequencesFwdFilePaths[i]).length() > 0)
						this.readAligner.alignReadsSlidingWindow(this.alignerBinPath, readSequencesFwdFilePaths[i], bowtieBasePath, alignmentFilePath, splitSeedSizes[i], splitSeedMismatches[i], this.maxMismatches, this.maxHits,this.maxReadLength,this.minReadLength,1, false, true, false, false, false);
					if(new File(readSequencesRevFilePaths[i]).length() > 0)
						this.readAligner.alignReadsSlidingWindow(this.alignerBinPath, readSequencesRevFilePaths[i], bowtieBasePath, alignmentFilePath, splitSeedSizes[i], splitSeedMismatches[i], this.maxMismatches, this.maxHits,this.maxReadLength,this.minReadLength,1, false, false, true, false, false);
					
				}
				currentTimePoint = System.currentTimeMillis();
				usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
				overallUsedTime += usedTime;
				if(this.verbose)
					System.out.println(String.format("%s\t Aligned read sequences (took me: %s s)",this.toString(), usedTime));
				
				
				//sort alignments by read id
				String sortedAlignmentFilePath = alignmentFilePath + ".sorted";
				if(!new File(alignmentFilePath).exists()) {
					date = new Date();
					if(this.verbose)
						System.out.println(String.format("%s\t Did not find any candidate alignments.",this.toString()));
				}
				
				else {
					prevTimePoint = System.currentTimeMillis();
					UnixSort unixSorter = new UnixSort(alignmentFilePath,sortedAlignmentFilePath,outputDir + String.format("/tmp_%s_%s_%s",chrName,windowStart,windowEnd),"\t",1,100,false,false,false);
					unixSorter.start();
					synchronized(unixSorter) {
						unixSorter.wait();
					}
					currentTimePoint = System.currentTimeMillis();
					usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
					overallUsedTime += usedTime;
					if(this.verbose)
						System.out.println(String.format("%s\t Sorted alignments by read id (took me: %s s)",this.toString(),usedTime));
					
					
					//combine alignments to obtain split candidates
					prevTimePoint = System.currentTimeMillis();
					combineSplitCandidateAlignments(candidates,sortedAlignmentFilePath, splitCandidateWriter, this.windowStart - 1, chrName);
					currentTimePoint = System.currentTimeMillis();
					usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
					overallUsedTime += usedTime;
					if(this.verbose)
						System.out.println(String.format("%s\t Obtained split candidates (took me: %s s)",this.toString(), usedTime));
					
				}
				
				//remove tmp files
				this.candidates.clear();
				prevTimePoint = System.currentTimeMillis();
				File[] tmpFiles = new File(this.outputDir).listFiles();
				for(File f : tmpFiles)
					f.delete();
				new File(this.outputDir).delete();
				currentTimePoint = System.currentTimeMillis();
				usedTime = (double)(currentTimePoint - prevTimePoint)/1000.0;
				overallUsedTime += usedTime;
				if(this.verbose)
					System.out.println(String.format("%s\t Everything done. Removing temporary files (overall time used for contig processing: %s s)", this.toString(),overallUsedTime));
				
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
		
		private void modifySamOutput(String inputFilePath, String outputFilePath, int seedLength) throws Exception {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			BufferedWriter bw = new BufferedWriter(new PrintWriter(new FileWriter(new File(outputFilePath))));
			StringTokenizer st;
			Pair<ArrayList<Integer>,Integer> mismatchPositionsAndReadLength = new Pair<ArrayList<Integer>,Integer>();
			ArrayList<Integer> mismatchPositions = new ArrayList<Integer>();
			int readLength;
			String chr;
			String start;
			char strand;
			String mdField;
			String mismatchInfo;
			int mdFieldPosition = -1;
			String readId;
			StringBuilder lineToWrite = new StringBuilder();
			boolean writeLine = true;
			int currentMismatchCount;
			int currentMismatchPosition;
			String line;
			while((line = br.readLine()) != null) {
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
				start = st.nextToken();
				
				st.nextToken();
				st.nextToken();
				st.nextToken();
				st.nextToken();
				st.nextToken();
				st.nextToken();
				
				for(int i = 10; i < mdFieldPosition; i++)
					st.nextToken();
				
				mdField = st.nextToken();
				
				
				writeLine = false;
				if(seedLength == splitSeedSizes[splitSeedSizes.length - 1])
					writeLine = true;
				
				if(!writeLine) {
					st = new StringTokenizer(mdField, ":");
					st.nextToken();
					st.nextToken();
					mismatchInfo = st.nextToken();
					setMismatchPositionsAndReadLengthUsingRegex(mismatchInfo, mismatchPositionsAndReadLength,mismatchPositions);
					mismatchPositions = mismatchPositionsAndReadLength.getFirst();
					readLength = mismatchPositionsAndReadLength.getSecond();
					
					if(strand == '-') {
						for(int i = 0; i < mismatchPositions.size(); i++)
							mismatchPositions.set(i, readLength - mismatchPositions.get(i));
					}
					
					int seedSizeIndex = 0;
					for(int i = 0; i < splitSeedSizes.length; i++) {
						if(splitSeedSizes[i] == seedLength) {
							seedSizeIndex = i + 1;
							break;
						}
					}
					
					currentMismatchCount = 0;
					currentMismatchPosition = 0;
					writeLine = true;
					
						
					while(seedSizeIndex < splitSeedSizes.length && writeLine) {
						for(int j = currentMismatchPosition; j < mismatchPositions.size(); j++) {
							if(mismatchPositions.get(j) < splitSeedSizes[seedSizeIndex]) {
								currentMismatchCount++;
								currentMismatchPosition++;
							}
							
							else {
								if(currentMismatchCount <= seedSizeIndex) {
									writeLine = false;
								}
								break;
							}
						}
						seedSizeIndex++;
					}
				}
				
				if(writeLine) {
					lineToWrite.setLength(0);
					lineToWrite.append(readId).append("\t").append(strand).append("\t").append(start).append("\t").append(mdField);
					bw.write(lineToWrite.toString());
					bw.newLine();
				}
			
			}
			br.close();
			bw.close();
		}
		
		/**
		 * returns 0-based mismatch positions and the read length
		 * 
		 * @param mdField
		 * @param maxAllowedMismatches
		 * @param splitSeedSizes
		 * @return
		 */
		
		private void setMismatchPositionsAndReadLength(String mismatchInfo,Pair<ArrayList<Integer>,Integer> mismatchPositionsAndReadLength, ArrayList<Integer> mismatchPositions) {
			
			mismatchPositions.clear();
			int readLength = 0;
			
			String matches = "0";
			int tmpNumber;
			int mismatchCount = 0;
			boolean considerFollowingLetters = true;
			for(int i = 0; i < mismatchInfo.length(); i++) {
				try {
					tmpNumber = Integer.valueOf(mismatchInfo.substring(i,i+1));
					matches += mismatchInfo.charAt(i);
					considerFollowingLetters = true;
				}
				catch(Exception e) {
					readLength += Integer.valueOf(matches);
					if(mismatchInfo.charAt(i) == '^')
						considerFollowingLetters = false;
					
					else if(considerFollowingLetters) {
						mismatchCount++;
						mismatchPositions.add(readLength);
						readLength++;
						
					}
					matches = "0";
				}
			}
			readLength += Integer.valueOf(matches);
			mismatchPositionsAndReadLength.setFirst(mismatchPositions);
			mismatchPositionsAndReadLength.setSecond(readLength);
		}
		
		
		private void setMismatchPositionsAndReadLengthUsingRegex(String mismatchInfo,Pair<ArrayList<Integer>,Integer> mismatchPositionsAndReadLength, ArrayList<Integer> mismatchPositions) {
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
		
		
		private void combineSplitCandidateAlignments(ArrayList<SplitCandidate> candidates, String sortedAlignmentFilePath, LineWriter splitCandidateWriter, int contigStart, String chrName) {
			String currentLine = null;
			try {
				BufferedRandomAccessFile candidateReader = new BufferedRandomAccessFile(new File(sortedAlignmentFilePath),"r",1024 * 10);
				ArrayList<SplitCandidate> fwdCandidatesFwdStrand = new ArrayList<SplitCandidate>();
				ArrayList<SplitCandidate> fwdCandidatesRevStrand = new ArrayList<SplitCandidate>();
				ArrayList<SplitCandidate> bwdCandidatesFwdStrand = new ArrayList<SplitCandidate>();
				ArrayList<SplitCandidate> bwdCandidatesRevStrand = new ArrayList<SplitCandidate>();
				
				int idComparisonResult;
				SplitCandidate candidateFromRealignment;
				SplitCandidate candidateFromInitialAlignment;
				ArrayList<SplitCandidate> candidatesFromRealignment;
				ArrayList<Integer> mismatchCountsA = new ArrayList<Integer>();
				ArrayList<Integer> mismatchCountsB = new ArrayList<Integer>();
				String readBaseNameFromRealignment;
				String readBaseNameFromInitialAlignment;
				String currentBaseName;
				
				String[] splittedLine;
				Pattern tabPattern = Pattern.compile("\t");
				long prevFilePointer;
				int prevCandidatesIndex;
				boolean isBwdAlignment;
				boolean foundSplitCandidate;
				String tmpId;
				String tmpStrand;
				int minMismatchValue;
				int readLength;
				int splitStart = -1;
				int splitEnd = -1;
				StringBuilder linesToWrite = new StringBuilder();
				int lineCounter = 0;
				HashSet<String> addedCandidates = new HashSet<String>();
				ArrayList<SplitCandidateLine> splitCandidateLines = new ArrayList<SplitCandidateLine>();
				for(int i = 0; i < candidates.size(); i++) {
					candidateFromInitialAlignment = candidates.get(i);
					
					
					if((currentLine = candidateReader.getNextLine()) == null)
						break;
					
					splittedLine = tabPattern.split(currentLine);
					readBaseNameFromInitialAlignment = candidateFromInitialAlignment.getReadId();
					if(readBaseNameFromInitialAlignment.length() >= 3 && readBaseNameFromInitialAlignment.substring(readBaseNameFromInitialAlignment.length() - 3).equals("/rc")) {
						readBaseNameFromInitialAlignment = readBaseNameFromInitialAlignment.substring(0,readBaseNameFromInitialAlignment.length() - 3);
					}
					readBaseNameFromRealignment = splittedLine[0];
					isBwdAlignment = false;
					if(readBaseNameFromRealignment.length() >= 3 && readBaseNameFromRealignment.substring(readBaseNameFromRealignment.length() - 3).equals("/rc")) {
						readBaseNameFromRealignment = readBaseNameFromRealignment.substring(0,readBaseNameFromRealignment.length() - 3);
						isBwdAlignment = true;
					}
					
					idComparisonResult = readBaseNameFromInitialAlignment.compareTo(readBaseNameFromRealignment);
					while(idComparisonResult < 0 && i < candidates.size()) {
						
						tmpId = candidateFromInitialAlignment.getReadId();
						tmpStrand = candidateFromInitialAlignment.getStrand();
						
						if(candidateFromInitialAlignment.getReadId().length() >= 3 && candidateFromInitialAlignment.getReadId().substring(candidateFromInitialAlignment.getReadId().length() -3).equals("/rc")) {
							tmpId = tmpId.substring(0,tmpId.length() - 3);
							if(tmpStrand.equals("+"))
								tmpStrand = "-";
							else
								tmpStrand = "+";
									
						}
						splitCandidateWriter.writeSplitCandidate(String.format("%s\tP\t%s\t%s\t%s\t%s\t%s\t%s\n",tmpId,chrName,candidateFromInitialAlignment.getStart(), ".",tmpStrand,maxMismatches + 1,getReadLength(candidateFromInitialAlignment.getMdField())));
						
						
						candidateFromInitialAlignment = candidates.get(++i);
						readBaseNameFromInitialAlignment = candidateFromInitialAlignment.getReadId();
						if(readBaseNameFromInitialAlignment.length() >= 3 && readBaseNameFromInitialAlignment.substring(readBaseNameFromInitialAlignment.length() - 3).equals("/rc")) {
							readBaseNameFromInitialAlignment = readBaseNameFromInitialAlignment.substring(0,readBaseNameFromInitialAlignment.length() - 3);
						}
						idComparisonResult = readBaseNameFromInitialAlignment.compareTo(readBaseNameFromRealignment);
					}
					
					while(idComparisonResult > 0 && (currentLine = candidateReader.getNextLine()) != null) {
						splittedLine = tabPattern.split(currentLine);
						readBaseNameFromRealignment = splittedLine[0];
						isBwdAlignment = false;
						if(readBaseNameFromRealignment.length() >= 3 && readBaseNameFromRealignment.substring(readBaseNameFromRealignment.length() - 3).equals("/rc")) {
							readBaseNameFromRealignment = readBaseNameFromRealignment.substring(0,readBaseNameFromRealignment.length() - 3);
							isBwdAlignment = true;
						}
						idComparisonResult = readBaseNameFromInitialAlignment.compareTo(readBaseNameFromRealignment);
					}
					
					if(idComparisonResult == 0) {
						fwdCandidatesFwdStrand.clear();
						fwdCandidatesRevStrand.clear();
						bwdCandidatesFwdStrand.clear();
						bwdCandidatesRevStrand.clear();
						prevFilePointer = candidateReader.getFilePointer();
						prevCandidatesIndex = i;
						
						while(currentLine  != null) {
							if(isBwdAlignment) {
								if(splittedLine[1].equals("+"))
									bwdCandidatesFwdStrand.add(new SplitCandidate(splittedLine[0],splittedLine[1],Integer.valueOf(splittedLine[2]),splittedLine[3]));
								else 
									bwdCandidatesRevStrand.add(new SplitCandidate(splittedLine[0],splittedLine[1],Integer.valueOf(splittedLine[2]),splittedLine[3]));
							}
							else {
								if(splittedLine[1].equals("+"))
									fwdCandidatesFwdStrand.add(new SplitCandidate(splittedLine[0],splittedLine[1],Integer.valueOf(splittedLine[2]),splittedLine[3]));
								else 
									fwdCandidatesRevStrand.add(new SplitCandidate(splittedLine[0],splittedLine[1],Integer.valueOf(splittedLine[2]),splittedLine[3]));
							}
							
							currentLine = candidateReader.getNextLine();
							if(currentLine == null)
								break;
							
							splittedLine = tabPattern.split(currentLine);
							readBaseNameFromRealignment = splittedLine[0];
							isBwdAlignment = false;
							if(readBaseNameFromRealignment.length() >= 3 && readBaseNameFromRealignment.substring(readBaseNameFromRealignment.length() - 3).equals("/rc")) {
								readBaseNameFromRealignment = readBaseNameFromRealignment.substring(0,readBaseNameFromRealignment.length() - 3);
								isBwdAlignment = true;
							}
							
							idComparisonResult = readBaseNameFromInitialAlignment.compareTo(readBaseNameFromRealignment);
							if(idComparisonResult != 0) {
								candidateReader.seek(prevFilePointer);
								break;
							}
							prevFilePointer = candidateReader.getFilePointer();
						}
						
						int j = i;
						splitCandidateLines.clear();
						for(; j < candidates.size(); j++) {
							candidateFromInitialAlignment = candidates.get(j);
							
							currentBaseName = candidateFromInitialAlignment.getReadId();
							isBwdAlignment = false;
							if(currentBaseName.length() >= 3 && currentBaseName.substring(currentBaseName.length() - 3).equals("/rc")) {
								currentBaseName = currentBaseName.substring(0,currentBaseName.length() - 3);
								isBwdAlignment = true;
							}
							
							if(!currentBaseName.equals(readBaseNameFromInitialAlignment)) {
								break;
							}
							
							foundSplitCandidate = false;
							if(isBwdAlignment) {
								if(candidateFromInitialAlignment.getStrand().equals("+")) {
									
									
								
									//search for: downstream hit of original read on - strand
									for(SplitCandidate currentCandidate : fwdCandidatesRevStrand) {
										if(currentCandidate.getStart() + contigStart > candidateFromInitialAlignment.getStart() - this.maxMismatches &&
											(currentCandidate.getReadId().length() < 3 || !currentCandidate.getReadId().substring(currentCandidate.getReadId().length() -3).equals("/rc")) &&
											currentCandidate.getStrand().equals("-")) {
											
											setMismatchCountsUsingRegex(candidateFromInitialAlignment.getMdField(),mismatchCountsA,false);
											setMismatchCountsUsingRegex(currentCandidate.getMdField(),mismatchCountsB, true);
											//candidate from initial alignment always has the full read length
											readLength = mismatchCountsA.size();
											splitStart = candidateFromInitialAlignment.getStart();
											splitEnd = currentCandidate.getStart() + contigStart + mismatchCountsB.size() - 1;
											minMismatchValue = getMinMismatchValueForSplitCandidate(mismatchCountsA,mismatchCountsB,splitStart,splitEnd,readLength,maxMismatches);
											
											if((splitEnd - splitStart - readLength <= this.maxGapSize) && ((splitEnd - splitStart - readLength >= this.minGapSize) || (splitEnd - splitStart - readLength <= this.maxDelSize)) && (splitStart + readLength - 1) != splitEnd && minMismatchValue <= this.maxMismatches && !addedCandidates.contains(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd))) {
												
												addedCandidates.add(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd));
												linesToWrite.append((String.format("%s\tS\t%s\t%s\t%s\t%s\t%s\t%s\n",currentCandidate.getReadId(),chrName,splitStart,splitEnd,"-",minMismatchValue,readLength)));
												lineCounter++;
												foundSplitCandidate = true;
											}
										}
									}
								}
								
								else {
									
									//search for: upstream hit of original read on + strand
									for(SplitCandidate currentCandidate : fwdCandidatesFwdStrand) {
										if(currentCandidate.getStart() + contigStart < candidateFromInitialAlignment.getStart() + this.maxMismatches &&
											(currentCandidate.getReadId().length() < 3 || !currentCandidate.getReadId().substring(currentCandidate.getReadId().length() -3).equals("/rc")) &&
											currentCandidate.getStrand().equals("+")) {
											setMismatchCountsUsingRegex(currentCandidate.getMdField(),mismatchCountsA, false);
											setMismatchCountsUsingRegex(candidateFromInitialAlignment.getMdField(),mismatchCountsB,true);
											
											//candidate from initial alignment always has the full read length
											readLength = mismatchCountsB.size();
											splitStart = currentCandidate.getStart() + contigStart;
											splitEnd= candidateFromInitialAlignment.getStart() + readLength - 1;
											minMismatchValue = getMinMismatchValueForSplitCandidate(mismatchCountsA,mismatchCountsB,splitStart,splitEnd,readLength,maxMismatches);
											
											if((splitEnd - splitStart - readLength <= this.maxGapSize) && ((splitEnd - splitStart - readLength >= this.minGapSize) || (splitEnd - splitStart - readLength <= this.maxDelSize)) && (splitStart + readLength - 1) != splitEnd && minMismatchValue <= this.maxMismatches && !addedCandidates.contains(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd))) {
												addedCandidates.add(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd));
												linesToWrite.append(String.format("%s\tS\t%s\t%s\t%s\t%s\t%s\t%s\n",currentCandidate.getReadId(),chrName,splitStart, splitEnd,"+",minMismatchValue,readLength));
												lineCounter++;
												foundSplitCandidate = true;
											}
										}
									}
								}
							}
							
							else {
								
								if(candidateFromInitialAlignment.getStrand().equals("+")) {
									
									//search for: downstream hit of reverse complemented read on - strand
									for(SplitCandidate currentCandidate : bwdCandidatesRevStrand) {
										if(currentCandidate.getStart() + contigStart > candidateFromInitialAlignment.getStart() - this.maxMismatches &&
											(currentCandidate.getReadId().length() >= 3 && currentCandidate.getReadId().substring(currentCandidate.getReadId().length() -3).equals("/rc")) &&
											currentCandidate.getStrand().equals("-")) {
											setMismatchCountsUsingRegex(candidateFromInitialAlignment.getMdField(),mismatchCountsA,false);
											setMismatchCountsUsingRegex(currentCandidate.getMdField(),mismatchCountsB, true);
											
											//candidate from initial alignment always has the full read length
											readLength = mismatchCountsA.size();
											splitStart = candidateFromInitialAlignment.getStart();
											splitEnd = currentCandidate.getStart() + contigStart + mismatchCountsB.size() - 1;
											minMismatchValue = getMinMismatchValueForSplitCandidate(mismatchCountsA,mismatchCountsB,splitStart,splitEnd,readLength,maxMismatches);
											
											if((splitEnd - splitStart - readLength <= this.maxGapSize) && ((splitEnd - splitStart - readLength >= this.minGapSize) || (splitEnd - splitStart - readLength <= this.maxDelSize)) && (splitStart + readLength - 1) != splitEnd && minMismatchValue <= this.maxMismatches && !addedCandidates.contains(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd))) {
												addedCandidates.add(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd));
												linesToWrite.append(String.format("%s\tS\t%s\t%s\t%s\t%s\t%s\t%s\n",candidateFromInitialAlignment.getReadId(),chrName,splitStart,splitEnd,"+",minMismatchValue,readLength));
												lineCounter++;
												foundSplitCandidate = true;
											}
										}
									}
								}
								else {
									
									// search for: upstream hit of reverse complemented read on + strand
									for(SplitCandidate currentCandidate : bwdCandidatesFwdStrand) {
										if(currentCandidate.getStart() + contigStart < candidateFromInitialAlignment.getStart() + this.maxMismatches &&
											(currentCandidate.getReadId().length() >= 3 && currentCandidate.getReadId().substring(currentCandidate.getReadId().length() -3).equals("/rc")) &&
											currentCandidate.getStrand().equals("+")) {
											setMismatchCountsUsingRegex(currentCandidate.getMdField(),mismatchCountsA, false);
											setMismatchCountsUsingRegex(candidateFromInitialAlignment.getMdField(),mismatchCountsB,true);
											
											//candidate from initial alignment always has the full read length
											readLength = mismatchCountsB.size();
											splitStart =  currentCandidate.getStart() + contigStart;
											splitEnd = candidateFromInitialAlignment.getStart() + readLength - 1;
											minMismatchValue = getMinMismatchValueForSplitCandidate(mismatchCountsA,mismatchCountsB,splitStart,splitEnd,readLength,maxMismatches);
											
											if((splitEnd - splitStart - readLength <= this.maxGapSize) && ((splitEnd - splitStart - readLength >= this.minGapSize) || (splitEnd - splitStart - readLength <= this.maxDelSize)) && (splitStart + readLength - 1) != splitEnd && minMismatchValue <= this.maxMismatches && !addedCandidates.contains(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd))) {
												addedCandidates.add(String.format("%s_%s_%s",currentBaseName,splitStart,splitEnd));
												linesToWrite.append(String.format("%s\tS\t%s\t%s\t%s\t%s\t%s\t%s\n",candidateFromInitialAlignment.getReadId(),chrName,splitStart, splitEnd,"-",minMismatchValue,readLength));
												lineCounter++;
												foundSplitCandidate = true;
											}
										}
									}
								}
							}
							
							if(lineCounter >= this.linesToWriteSize) {
								splitCandidateWriter.writeSplitCandidate(linesToWrite.toString());
								linesToWrite.setLength(0);
								lineCounter = 0;
							}
							
							if(!foundSplitCandidate) {
								
								tmpId = candidateFromInitialAlignment.getReadId();
								tmpStrand = candidateFromInitialAlignment.getStrand();
								
								if(candidateFromInitialAlignment.getReadId().length() >= 3 && candidateFromInitialAlignment.getReadId().substring(candidateFromInitialAlignment.getReadId().length() -3).equals("/rc")) {
									tmpId = tmpId.substring(0,tmpId.length() - 3);
									if(tmpStrand.equals("+"))
										tmpStrand = "-";
									else
										tmpStrand = "+";
											
								}
								
								splitCandidateWriter.writeSplitCandidate(String.format("%s\tP\t%s\t%s\t%s\t%s\t%s\t%s\n",tmpId,chrName,candidateFromInitialAlignment.getStart(), ".",tmpStrand,maxMismatches + 1,getReadLength(candidateFromInitialAlignment.getMdField())));
							}
						}
						i = j - 1;
					}
				}
				//write out all candidates at once (all other threads have to wait for this operation!)
				splitCandidateWriter.writeSplitCandidate(linesToWrite.toString());
				candidateReader.close();
			}
			catch(Exception e) {
				e.printStackTrace();
				System.out.println(currentLine);
				System.exit(1);
			}
		}
		
		private class SplitCandidateLine {
			private String line;
			private int minMismatchValue;
			
			public SplitCandidateLine(String line, int minMismatchValue) {
				this.line = line;
				this.minMismatchValue = minMismatchValue;
			}
			
			
			public String getLine() {
				return this.line;
			}
			
			public int getMismatchValue() {
				return this.minMismatchValue;
			}
		}
		
		/**
		 * 
		 * @param mismatchCountsA -> split part A, which begins with the first base of the read
		 * @param mismatchCountsB -> split part B, which ends with the last base of the read
		 * @return
		 */
		
		/*private int getMinMismatchValueForSplitCandidate(ArrayList<Integer> mismatchCountsA, ArrayList<Integer> mismatchCountsB) {
			int minMismatches = Integer.MAX_VALUE;
			int currentMismatches;
			
			int offset = Math.min(mismatchCountsA.size(), mismatchCountsB.size()) - 1;
			int offsetA = 0;
			int offsetB = 2;
			
			if(mismatchCountsA.size() > mismatchCountsB.size()) {
				offsetA = offset;
				offsetB = 1;
			}
			
			else if(mismatchCountsA.size() == mismatchCountsB.size()) {
				offset = mismatchCountsA.size() - 2;
				
			}
			
			
						
			for(int i = 0; i <= offset; i++) {
				currentMismatches = mismatchCountsA.get(i + offsetA) + mismatchCountsB.get(mismatchCountsB.size() - (i + offsetB));
				if(currentMismatches == 0)
					return 0;
				
				
				else if(currentMismatches < minMismatches)
					minMismatches = currentMismatches;
			}
			
			return minMismatches;
		}
		
		*/
		
		/**
		 * 
		 * @param mismatchCountsA -> split part A, which begins with the first base of the read
		 * @param mismatchCountsB -> split part B, which ends with the last base of the read
		 * @return
		 */
		
		
		private int getMinMismatchValueForSplitCandidate(ArrayList<Integer> mismatchCountsA, ArrayList<Integer> mismatchCountsB, int splitStart, int splitEnd, int readLength, int maxMismatches) {
			
			//check if we are processing an insertion candidate
			int insertionSize = readLength - (splitEnd - splitStart + 1);
			if(insertionSize > this.maxDelSize)
				return Integer.MAX_VALUE;
			
			
			int minMismatches = Integer.MAX_VALUE;
			int currentMismatches;
			int offset = Math.min(mismatchCountsA.size(), mismatchCountsB.size()) - 1;
			int offsetA = 0;
			int offsetB = 2;
			
			if(mismatchCountsA.size() > mismatchCountsB.size()) {
				offsetA = offset;
				offsetB = 1;
			}
			
			else if(mismatchCountsA.size() == mismatchCountsB.size()) {
				offset = mismatchCountsA.size() - 2;
				
			}
			
			
						
			for(int i = 0; i <= offset; i++) {
				
				if(splitEnd - splitStart + 1 < readLength) {
					
					if(mismatchCountsB.size() - (i + offsetB) - insertionSize < 0)
						break;
					
					currentMismatches = mismatchCountsA.get(i + offsetA) + mismatchCountsB.get(mismatchCountsB.size() - (i + offsetB) - insertionSize);
				}
				
				else
					currentMismatches = mismatchCountsA.get(i + offsetA) + mismatchCountsB.get(mismatchCountsB.size() - (i + offsetB));
				
				if(currentMismatches == 0)
					return 0;
				
				
				else if(currentMismatches < minMismatches)
					minMismatches = currentMismatches;
			}
			
			return minMismatches;
		}
		
		private void setMismatchCounts(String mdField,ArrayList<Integer> mismatchCounts, boolean reverse) {
			try {
				mismatchCounts.clear();
				int readLength = 0;
				
				String matches = "0";
				int tmpNumber;
				int mismatchCount = 0;
				
				StringTokenizer st = new StringTokenizer(mdField,":");
				st.nextToken();
				st.nextToken();
				String mismatchInfo = st.nextToken();
				
				boolean considerFollowingLetters = true;
				for(int i = 0; i < mismatchInfo.length(); i++) {
					try {
						tmpNumber = Integer.valueOf(mismatchInfo.substring(i,i+1));
						matches += mismatchInfo.charAt(i);
						considerFollowingLetters = true;
					}
					catch(Exception e) {
						
						for(int j = readLength; j < readLength + Integer.valueOf(matches); j++) {
							mismatchCounts.add(mismatchCount);
						}
						
						readLength += Integer.valueOf(matches);
						if(mismatchInfo.charAt(i) == '^')
							considerFollowingLetters = false;
						
						
						else if(considerFollowingLetters) {
							mismatchCount++;
							mismatchCounts.add(mismatchCount);
							readLength++;
							
						}
						matches = "0";
					}
				}
				for(int j = readLength; j < readLength + Integer.valueOf(matches); j++) {
					mismatchCounts.add(mismatchCount);
				}
				
				if(reverse) {
					int maxMismatchCount = mismatchCounts.get(mismatchCounts.size()-1);
					for(int i = 0; i < mismatchCounts.size(); i++) {
						mismatchCounts.set(i, maxMismatchCount - mismatchCounts.get(i));
					}
				}
			}
			catch(Exception e) {
				e.printStackTrace();
				System.err.println(mdField);
			}
		}
		
		private int getReadLength(String mdField) {
			int readLength = 0;
			StringTokenizer st = new StringTokenizer(mdField,":");
			st.nextToken();
			st.nextToken();
			String mismatchInfo = st.nextToken();
			
			Pattern pattern = Pattern.compile("[0-9]+|[A-Z,#]+");
			mismatchInfo = mismatchInfo.replaceAll("\\^[A-Z]+", "#");
			Matcher matcher = pattern.matcher(mismatchInfo);
			ArrayList<String> chunks = new ArrayList<String>();
	    	while (matcher.find()) {
    	        chunks.add( matcher.group() );
    	    }
	    	
	    	if(chunks.size() == 1) {
	    		readLength = Integer.valueOf(chunks.get(0));
	    	}
	    	
	    	else {
	    		for(int i = 0; i < chunks.size() - 1; i+=2) {
					readLength += Integer.valueOf(chunks.get(i));
					if(!chunks.get(i+1).equals("#")) {
						for(int j = 0; j < chunks.get(i+1).length(); j++) {
							readLength++;
						}
					}
	    		}
	    		readLength += Integer.valueOf(chunks.get(chunks.size() - 1));
	    	}
	    	return readLength;
		}
		
		
		private void setMismatchCountsUsingRegex(String mdField,ArrayList<Integer> mismatchCounts, boolean reverse) {
			try {
				mismatchCounts.clear();
				
				
				StringTokenizer st = new StringTokenizer(mdField,":");
				st.nextToken();
				st.nextToken();
				String mismatchInfo = st.nextToken();
				
				Pattern pattern = Pattern.compile("[0-9]+|[A-Z,#]+");
				mismatchInfo = mismatchInfo.replaceAll("\\^[A-Z]+", "#");
				Matcher matcher = pattern.matcher(mismatchInfo);
				ArrayList<String> chunks = new ArrayList<String>();
		    	while (matcher.find()) {
	    	        chunks.add( matcher.group() );
	    	    }
		    	
		    	if(reverse && chunks.size() > 1) {
		    		Collections.reverse(chunks);
		    	}
				
		    	int readLength = 0;
				int mismatchCount = 0;
		    	if(chunks.size() == 1) {
		    		readLength = Integer.valueOf(chunks.get(0));
		    		for(int i = 0; i < readLength; i++)
		    			mismatchCounts.add(0);
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
		    	}
				
			}
			catch(Exception e) {
				e.printStackTrace();
				System.err.println(mdField);
			}
		}
		
		
		private void createContigSequence(String outputFilePath, StringBuilder referenceSequence, int windowStart, int windowEnd) throws Exception {
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			pw.println(String.format(">%s_%s", windowStart,windowEnd));
			pw.println(referenceSequence.substring(windowStart, windowEnd));
			pw.close();
		}
		
		
		private void createReadSequencesAndParseCandidates(String[] outputFilePathsFwdStrand, String[] outputFilePathsRevStrand, File alignmentFile, long windowStartPointer, long windowEndPointer) throws Exception {
			PrintWriter[] fwdPws = new PrintWriter[outputFilePathsFwdStrand.length];
			PrintWriter[] revPws = new PrintWriter[outputFilePathsFwdStrand.length];
			for(int i = 0; i < fwdPws.length; i++) {
				fwdPws[i] = new PrintWriter(new FileWriter(new File(outputFilePathsFwdStrand[i])));
				revPws[i] = new PrintWriter(new FileWriter(new File(outputFilePathsRevStrand[i])));
			}
			
			BufferedRandomAccessFile br = new BufferedRandomAccessFile(alignmentFile,"r",1024 * 1024);
			br.seek(windowStartPointer);
			String currentLine;
			StringTokenizer st;
			String readId;
			String strand;
			int start;
			StringBuilder sequence = new StringBuilder();
			StringBuilder matches = new StringBuilder();
			String mdField;
			int readLength;
			int splitSizeIndex = -1;
			String cigar;
			boolean halfReadMapping;
			HashSet<String> alreadyAddedReads = new HashSet<String>();
			while((currentLine = br.getNextLine()) != null) {
				st = new StringTokenizer(currentLine, "\t");
				readId = st.nextToken();
				strand = st.nextToken();
				start = Integer.valueOf(st.nextToken());
				sequence.setLength(0);
				sequence.append(st.nextToken());
				mdField = st.nextToken();
				readLength = Integer.valueOf(st.nextToken());
				splitSizeIndex = Integer.valueOf(st.nextToken());
				cigar = st.nextToken();
				
				halfReadMapping = isHalfReadMappingValid(cigar,readLength);
				
				
				this.candidates.add(new SplitCandidate(readId,strand,start,mdField));
				
				//backward alignment
				if(readId.length() >= 2 && readId.substring(readId.length() - 2).equals("rc")) {
					if(strand.equals("+")) {
						readId = readId.substring(0,readId.lastIndexOf('/'));
						if(!alreadyAddedReads.contains(readId + "-")) {
							//sequence field: reverse complement of the original read
							//search for: downstream hit of original read on - strand
							alreadyAddedReads.add(readId + "-");
							sequence.reverse();
							for(int i = 0; i < sequence.length(); i++) {
								sequence.setCharAt(i, substitute(sequence.charAt(i)));
							}
							
							revPws[splitSizeIndex].println(">" + readId);
							if(halfReadMapping)
								revPws[splitSizeIndex].println(sequence.toString().substring(0,sequence.length()/2));
							else
								revPws[splitSizeIndex].println(sequence.toString());
						}
						
					}
					else {
						readId = readId.substring(0,readId.lastIndexOf('/'));
						if(!alreadyAddedReads.contains(readId + "+")) {
							//sequence field: original read
							//search for: upstream hit of original read on + strand
							alreadyAddedReads.add(readId + "+");
							fwdPws[splitSizeIndex].println(">" + readId);
							if(halfReadMapping)
								fwdPws[splitSizeIndex].println(sequence.toString().substring(0,sequence.length()/2));
							else
								fwdPws[splitSizeIndex].println(sequence.toString());
						}
						
					}
				}
				
				//forward alignment
				else {
					if(strand.equals("+")) {
						readId += "/rc";
						if(!alreadyAddedReads.contains(readId + "-")) {
							alreadyAddedReads.add(readId + "-");
							// sequence field: original read
							// search for: downstream hit of reverse complemented read on - strand
							sequence.reverse();
							for(int i = 0; i < sequence.length(); i++) {
								sequence.setCharAt(i, substitute(sequence.charAt(i)));
							}
							
							revPws[splitSizeIndex].println(">" + readId);
							if(halfReadMapping)
								revPws[splitSizeIndex].println(sequence.toString().substring(0,sequence.length()/2));
							else
								revPws[splitSizeIndex].println(sequence.toString());
						}
					}
					
					else {
						readId += "/rc";
						if(!alreadyAddedReads.contains(readId + "+")) {
							alreadyAddedReads.add(readId + "+");
							// sequence field: reverse complement of original read
							// search for: upstream hit of reverse complemented read on + strand
							fwdPws[splitSizeIndex].println(">" + readId);
							if(halfReadMapping)
								fwdPws[splitSizeIndex].println(sequence.toString().substring(0,sequence.length()/2));
							else
								fwdPws[splitSizeIndex].println(sequence.toString());
						}
						
					}
				}
				
				
				if(br.getFilePointer() == windowEndPointer)
					break;
			}
			br.close();
			for(int i = 0; i < fwdPws.length; i++) {
				fwdPws[i].close();
				revPws[i].close();
			}
			
		}
		
		
		
		private boolean isHalfReadMappingValid(String cigar, int readLength) {
			boolean clippedAtTheStart = false;
	    	boolean clippedAtTheEnd = false;
	    	ArrayList<String> clippings = new ArrayList<String>();
	    	Pattern pattern = Pattern.compile("[0-9]+[H|S]");
			Matcher matcher = pattern.matcher(cigar);
	    	
			while(matcher.find()) {
				if(matcher.start() == 0) {
					clippedAtTheStart = true;
				}
				else
					clippedAtTheEnd = true;
				
				clippings.add(matcher.group());
			}
			
			if(clippings.isEmpty())
				return false;
			
			//either clipped at the start or end (if there is clipping at both ends, one clipping will be much shorter than the other).
			int clippingLengthStart = 0;
			int clippingLengthEnd = 0;
			if(clippedAtTheStart) {
				clippingLengthStart = Integer.valueOf(clippings.get(0).substring(0,clippings.get(0).length() - 1));
			}
			
			if(clippedAtTheEnd) {
				clippingLengthEnd = Integer.valueOf(clippings.get(clippings.size() - 1).substring(0,clippings.get(clippings.size() - 1).length() - 1));
			}
			
			
			int longestClip = clippingLengthStart;
			if(clippingLengthEnd > clippingLengthStart) longestClip = clippingLengthEnd;
			
			if(longestClip <= readLength/2)
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
		
		
		
		
	}
	
	
	private class SplitCandidate {
		
		private String readId;
		private String strand;
		private int start;
		private String mdField;
		private boolean readReverseComplemented;
		
		public SplitCandidate(String readId, String strand, int start, String mdField) {
			this.readId = readId;
			this.readReverseComplemented = false;
			if(readId.length() >= 2 && readId.substring(readId.length() - 2).equals("rc"))
				readReverseComplemented = true;
			
			this.strand = strand;
			this.start = start;
			this.mdField = mdField;
		}

		public String getReadId() {
			return readId;
		}

		
		public String getStrand() {
			return this.strand;
		}

		public int getStart() {
			return start;
		}

		public String getMdField() {
			return mdField;
		}

		public boolean isReadReverseComplemented() {
			return readReverseComplemented;
		}
		
		
	}
	
	private class SplitCandidateReadIdComparator implements Comparator<SplitCandidate> {

		@Override
		public int compare(SplitCandidate o1, SplitCandidate o2) {
			// TODO Auto-generated method stub
			return o1.getReadId().compareTo(o2.getReadId());
		}
		
	}
	
	private class SplitCandidateStartComparator implements Comparator<SplitCandidate> {

		@Override
		public int compare(SplitCandidate o1, SplitCandidate o2) {
			// TODO Auto-generated method stub
			return Integer.valueOf(o1.getStart()).compareTo(o2.getStart());
		}
		
	}
	
	
	private class LineWriter {
		
		private BufferedWriter bw;
		
		public LineWriter(BufferedWriter bw) {
			this.bw = bw;
		}
		
		
		public synchronized void writeSplitCandidate(String line) {
			try {
				this.bw.write(line);
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	

	public void concatenateFilesWithNIO(ArrayList<String> filePaths, String outputPath) {
		try {
			File checkFile = new File(outputPath);
			if(checkFile.exists())
				checkFile.delete();
			
			//move largest file to output path
			int indexOfLargestFile = getIndexOfLargestFile(filePaths);
			if(indexOfLargestFile == -1) {
				return;
			}
			
			new File(filePaths.get(indexOfLargestFile)).renameTo(new File(outputPath));
			filePaths.remove(indexOfLargestFile);
			
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
				new File(filePath).delete();
			}
			fos.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private int getIndexOfLargestFile(ArrayList<String> rmapFilePaths) {
		int index = -1;
		String path;
		long currentSize;
		long largestSize = Long.MIN_VALUE;
		for(int i = 0; i < rmapFilePaths.size(); i++) {
			path = rmapFilePaths.get(i);
			
			if(!new File(path).exists())
				continue;
			
			currentSize = new File(path).length();
			if(currentSize > largestSize) {
				index = i;
				largestSize = currentSize;
			}
		}
		return index;
	}
	
}
