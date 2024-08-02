package tools;

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
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import main.Pair;
import main.Statistic;


public class SamProcessor {
	
	private final int bunchSize;

	public SamProcessor() {
		this.bunchSize = 10000000;
	}
	
	public void splitByChromosome(String inputFilePath, String outputDirPath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			HashMap<String,PrintWriter> chr2writer = new HashMap<String,PrintWriter>();
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				if(!chr2writer.containsKey(splittedLine[2]))
					chr2writer.put(splittedLine[2], new PrintWriter(new FileWriter(new File(outputDirPath + "/" + splittedLine[2] + ".sam"))));
				
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
	
	
	public void getStartPositionsOnly(String samFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			String readId;
			int alignmentStart;
			String cigar;
			ArrayList<Integer> startPositions;
			int idCounter;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '@') {
					pw.println(currentLine);
					continue;
				}
				
				splittedLine = tabPattern.split(currentLine);
				readId = splittedLine[0];
				alignmentStart = Integer.valueOf(splittedLine[3]);
				cigar = splittedLine[5];
				splittedLine[5] = "1M";
				pw.print(splittedLine[0]);
				for(int i = 1; i < splittedLine.length; i++) {
					pw.print("\t" + splittedLine[i]);
				}
				pw.println();
				
				if(cigar.contains("N")) {
					startPositions = getSplitStartPositions(alignmentStart,cigar);
					idCounter = 1;
					for(int startPos : startPositions) {						
						pw.print(readId + "_" + ++idCounter);
						splittedLine[3] = String.valueOf(startPos);
						for(int i = 1; i < splittedLine.length; i++) {
							pw.print("\t" + splittedLine[i]);
						}
						pw.println();
					}
				}
					
				
			
			}
			br.close();
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private ArrayList<Integer> getSplitStartPositions(int alignmentStart, String cigar) {
		
		ArrayList<Integer> startPositions = new ArrayList<Integer>();
		
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
			
			if(currentOp.equals("M") || currentOp.equals("X") || currentOp.equals("=")) {
				currentPosition += currentLength - 1;
			}
			
			else if(currentOp.equals("N") || currentOp.equals("D")) {
				currentPosition += currentLength + 1;
				if(currentOp.equals("N"))
					startPositions.add(currentPosition);
			}
			
			else if(currentOp.equals("I") && i < matches.size() - 1 && !matches.get(i+1).substring(matches.get(i+1).length() - 1).equals("N") && !matches.get(i+1).substring(matches.get(i+1).length() - 1).equals("D"))
				currentPosition++;
			
		}
		return startPositions;
	}
	
	
	public void removeMateInformationFromId(String inputFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			StringBuilder sb = new StringBuilder();
			 while((currentLine = br.readLine()) != null) {
				 
				 if(currentLine.charAt(0) == '@') {
					 pw.println(currentLine);
					 continue;
				 }
				 
				 splittedLine = tabPattern.split(currentLine);
				 if(splittedLine[0].charAt(splittedLine[0].length() - 2) == '/') {
					 splittedLine[0] = splittedLine[0].substring(0,splittedLine[0].length() - 2);
				 }
				 
				 else if(splittedLine[0].charAt(splittedLine[0].length() - 1) == 'a' || splittedLine[0].charAt(splittedLine[0].length() - 1) == 'b') {
					 splittedLine[0] = splittedLine[0].substring(0,splittedLine[0].length() - 1);
				 }
				 
				 sb.setLength(0);
				 sb.append(splittedLine[0]);
				 for(int i = 1; i < splittedLine.length; i++) {
					 sb.append("\t").append(splittedLine[i]);
				 }
				 pw.println(sb.toString());
			 }
			 br.close();
			 pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	public void addMateInformationToFlag(String inputFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			int flag;
			StringBuilder sb = new StringBuilder();
			
			while((currentLine = br.readLine()) != null) {
				
				if(currentLine.charAt(0) == '@') {
					pw.println(currentLine);
					continue;
				}
				
				splittedLine = tabPattern.split(currentLine);
				flag = Integer.valueOf(splittedLine[1]);
				if(splittedLine[0].charAt(splittedLine[0].length() - 2) == '/' && ((flag & (1L << 6)) == 0) && ((flag & (1L << 7)) == 0)) {
					if(splittedLine[0].charAt(splittedLine[0].length() - 1) == '1') {
						flag += Math.pow(2, 6);
					}
					else if(splittedLine[0].charAt(splittedLine[0].length() - 1) == '2') {
						flag += Math.pow(2, 7);
					}
					
					sb.setLength(0);
					sb.append(splittedLine[0]).append("\t").append(flag);
					for(int i = 2; i < splittedLine.length; i++) {
						sb.append("\t").append(splittedLine[i]);
					}
					pw.println(sb.toString());
				}
				else {
					pw.println(currentLine);
				}
			}
			br.close();
			pw.close();
			
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	public void addMateInformationToId(String inputFilePath, String outputFilePath) {
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			StringBuilder sb = new StringBuilder();
			
			int flag = 0;
			while((currentLine = br.readLine()) != null) {
				
				if(currentLine.charAt(0) == '@')
					continue;
				
				splittedLine = tabPattern.split(currentLine);
				flag = Integer.valueOf(splittedLine[1]);
				if((flag & (1L << 6)) != 0) {
					splittedLine[0] += "/1";
				}
				else if((flag & (1L << 7)) != 0) {
					splittedLine[0] += "/2";
				}
				
				
				sb.setLength(0);
				sb.append(splittedLine[0]).append("\t").append(flag);
				for(int i = 2; i < splittedLine.length; i++) {
					sb.append("\t").append(splittedLine[i]);
				}
				pw.println(sb.toString());
			}
			
			br.close();
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	
	
	/**
	 * this method only works on sam files sorted by read ids
	 * it extracts unmapped and multi mapped reads as well as reads which were truncated
	 * @param samFilePath
	 * @param readFilePath
	 * @param readFormat
	 * @param readLength
	 * @param outputDir
	 */
	public void extractUnmappedAndMultimappedReads(String samFilePath, String readFilePath, String readFormat, int readLength, String outputDir) {
		try {
			BufferedRandomAccessFile readsReader = new BufferedRandomAccessFile(new File(readFilePath), "r",5000000);
			long prevFilePointer = readsReader.getFilePointer();
			HashSet<String> currentReadIds = parseNextBunchOfReads(readsReader,readFormat);
			long currentFilePointer = readsReader.getFilePointer();
			
			String unmappedReadsOutputPath = outputDir + "/unmapped_reads.fa";
			String uniquelyMappedReadsOutputPath = outputDir + "/uniquely_mapped_reads.sam";
			
			PrintWriter unmappedReadsPw = new PrintWriter(new FileWriter(new File(unmappedReadsOutputPath)));
			PrintWriter uniquelyMappedReadsPw = new PrintWriter(new FileWriter(new File(uniquelyMappedReadsOutputPath)));
			while(currentReadIds.size() > 0) {
				//filter for unmapped and multi mapped reads
				removeUniquelyMappedReads(currentReadIds,samFilePath,readLength,uniquelyMappedReadsPw);
				//move the reader back to previous position and write the unmapped reads of the current read set
				readsReader.seek(prevFilePointer);
				writeNextBunchOfReads(currentReadIds,readsReader, unmappedReadsPw, readFormat, currentFilePointer);
				
				//move the reader back and parse the next bunch of reads
				readsReader.seek(currentFilePointer);
				prevFilePointer = currentFilePointer;
				currentReadIds = parseNextBunchOfReads(readsReader,readFormat);
				currentFilePointer = readsReader.getFilePointer();
			}
			readsReader.close();
			unmappedReadsPw.close();
			uniquelyMappedReadsPw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	public void getReverseComplementOfUnalignedReads(String samFilePath, String readFilePath, String readFormat, String outputFilePath) {
		try {
			BufferedRandomAccessFile readsReader = new BufferedRandomAccessFile(new File(readFilePath), "r",5000000);
			long prevFilePointer = readsReader.getFilePointer();
			HashSet<String> currentReadIds = parseNextBunchOfReads(readsReader,readFormat);
			long currentFilePointer = readsReader.getFilePointer();
			
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			while(currentReadIds.size() > 0) {
				//filter for unaligned reads
				removeAlignedReads(currentReadIds,samFilePath);
				//move the reader back to previous position and write the unaligned reads of the current read set
				readsReader.seek(prevFilePointer);
				writeNextBunchOfReverseComplementedReads(currentReadIds,readsReader, pw, readFormat, currentFilePointer);
				
				//move the reader back and parse the next bunch of reads
				readsReader.seek(currentFilePointer);
				prevFilePointer = currentFilePointer;
				currentReadIds = parseNextBunchOfReads(readsReader,readFormat);
				currentFilePointer = readsReader.getFilePointer();
			}
			readsReader.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void convertSamToRmap(String samFilePath, String outputPath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			String readId;
			int bitFlag;
			String binaryString;
			String strand;
			String chr;
			String mappingType;
			int startPosition;
			int endPosition;
			String mappingPosition;
			String cigString;
			int missmatches;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				readId = splittedLine[0];
				bitFlag = Integer.valueOf(splittedLine[1]);
				if(bitFlag == 0)
					strand = "+";
				else { 
					binaryString = Integer.toBinaryString(bitFlag);
					//read is unmapped
					if(binaryString.charAt(binaryString.length() - 3) == '1') {
						continue;
					}
					//read is from the reverse strand
					if(binaryString.charAt(binaryString.length() - 5) == '1')
						strand = "-";
					else
						strand = "+";
				}
				chr = splittedLine[2];
				startPosition = Integer.valueOf(splittedLine[3]);
				cigString = splittedLine[5];
				endPosition = getEndPosition(startPosition,cigString);
				missmatches = getMissmatchCount(splittedLine);
				if(!cigString.contains("D") && !cigString.contains("N")) {
					mappingType = "F";
					mappingPosition = String.valueOf(startPosition) + "\t.";
				}
				else {
					mappingType = "S";
					mappingPosition = String.format("%s\t%s",startPosition,endPosition);
				}
						
				pw.println(String.format("%s\t%s\t%s\t%s\t%s\t%s",readId,mappingType,chr,mappingPosition,strand,missmatches));
			}
			br.close();
			pw.close();
			       
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private int getMissmatchCount(String[] splittedSamLine) {
		int missmatchCount = Integer.MAX_VALUE;
		//sam optional fields start at field 12
		String[] splittedTag;
		Pattern doublePointPattern = Pattern.compile(":");
		for(int i = 11; i < splittedSamLine.length; i++) {
			splittedTag = doublePointPattern.split(splittedSamLine[i]);
			if(splittedTag[0].equals("NM")) {
				missmatchCount = Integer.valueOf(splittedTag[splittedTag.length - 1]);
				break;
			}
		}
		return missmatchCount;
	}
	
	private int getEndPosition(int startPosition, String cigString){
		int endPosition = startPosition;
		int tmpNumber;
		String tmpNumberString = "";
		for(int i = 0; i < cigString.length(); i++) {
			try {
				tmpNumber = Integer.valueOf(cigString.substring(i,i+1));
				tmpNumberString += cigString.charAt(i);
			}
			catch(Exception e) {
				tmpNumber = Integer.valueOf(tmpNumberString);
				if(cigString.charAt(i) == 'N' || cigString.charAt(i) == 'D' || cigString.charAt(i) == 'M' || cigString.charAt(i) == '=' || cigString.charAt(i) == 'X') {
					endPosition += tmpNumber;
				}
				tmpNumberString = "";
			}
		}
		return endPosition - 1;
	}

	private HashSet<String> parseNextBunchOfReads(BufferedRandomAccessFile readsReader, String readFormat) throws Exception {
		HashSet<String> readIds = new HashSet<String>();
		int parsedReads = 0;
		String currentLine;
		while((currentLine = readsReader.getNextLine()) != null) {
			readIds.add(currentLine.substring(1));
	
			if(readFormat.equals("fasta") || readFormat.equals("fa")) {
				readsReader.getNextLine();
			}
			//in case we have fastq reads we have to skip three lines here
			else {
				readsReader.getNextLine();
				readsReader.getNextLine();
				readsReader.getNextLine();
			}
			
			if(++parsedReads == this.bunchSize)
				return readIds;
		}
		
		return readIds;
	}
	
	//it is guaranteed that the reader always points to the header of a read
	private void writeNextBunchOfReads(HashSet<String> unmappedReads,BufferedRandomAccessFile readsReader, PrintWriter pw, String readFormat, long stopFilePointer) throws Exception {
		String currentLine;
		while((currentLine = readsReader.getNextLine()) != null) {
			if(unmappedReads.contains(currentLine.substring(1))) {
				pw.println(">" + currentLine.substring(1));
				pw.println(readsReader.getNextLine());
				
				if(readFormat.equals("fastq")) {
					readsReader.getNextLine();
					readsReader.getNextLine();
				}
			}
			if(readsReader.getFilePointer() == stopFilePointer)
				break;
		}
	}
	
	private void writeNextBunchOfReverseComplementedReads(HashSet<String> unmappedReads,BufferedRandomAccessFile readsReader, PrintWriter pw, String readFormat, long stopFilePointer) throws Exception {
		String currentLine;
		StringBuilder sb = new StringBuilder();
		while((currentLine = readsReader.getNextLine()) != null) {
			if(unmappedReads.contains(currentLine.substring(1))) {
				pw.println(currentLine + "/rc");
				sb.setLength(0);
				sb.append(readsReader.readLine());
				sb.reverse();
				for(int i = 0; i < sb.length(); i++) {
					sb.setCharAt(i, substitute(sb.charAt(i)));
				}
				pw.println(sb.toString());
				
				if(readFormat.equals("fastq")) {
					pw.println(readsReader.readLine());
					sb.setLength(0);
					sb.append(readsReader.readLine());
					pw.println(sb.reverse().toString());
				}
			}
			if(readsReader.getFilePointer() == stopFilePointer)
				break;
		}
	}
	
	private char substitute(char n) {
		if(n == 'A' || n == 'a') return 'T';
		if(n == 'T' || n == 't') return 'A';
		if(n == 'C' || n == 'c') return 'G';
		if(n == 'G' || n == 'g') return 'C';
		else return 'N';
	}
	
	private void removeAlignedReads(HashSet<String> readIds, String samFilePath) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
		String currentLine = br.readLine();
		while(currentLine.charAt(0) == '@') {
			currentLine = br.readLine();
		}
		
		Pattern tabPattern = Pattern.compile("\t");
		String[] splittedLine = tabPattern.split(currentLine);
		if(!splittedLine[2].equals("*"))
			readIds.remove(splittedLine[0]);
		while((currentLine = br.readLine()) != null) {
			splittedLine = tabPattern.split(currentLine);
			if(!splittedLine[2].equals("*"))
				readIds.remove(splittedLine[0]);
		}
		br.close();
	}
	
	
	private void removeUniquelyMappedReads(HashSet<String> readIds, String samFilePath, int readLength, PrintWriter pw) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
		String currentLine = br.readLine();
		
		//skipping the header
		while(currentLine.charAt(0) == '@' && br.ready()) {
			currentLine = br.readLine();
		}
		if(!br.ready())
			return;
		
		String prevLine = currentLine;
		Pattern tabPattern = Pattern.compile("\t");
		String[] splittedLine = tabPattern.split(currentLine);
		String readId = splittedLine[0];
		String chr = splittedLine[2];
		String prevReadId = readId;
		String cigString = splittedLine[5];
		boolean clipped = false;
		if(cigString.contains("S") || cigString.contains("H"))
			clipped = true;
		
		int currentReadLength = getReadLength(cigString); 
		boolean multiMapped = false;
		
		while(br.ready()) {
			currentLine = br.readLine();
			
			if(currentLine.charAt(0) == '@')
				continue;
			
			splittedLine = tabPattern.split(currentLine);
			readId = splittedLine[0];
			chr = splittedLine[2];
			cigString = splittedLine[5];
			if(readId.equals(prevReadId))
				multiMapped = true;
			
			else {
				//if(!multiMapped && !clipped && currentReadLength == readLength && readIds.contains(prevReadId)) {
				//we do not check for the 'correct' read length anymore since we allow for variable read lengths
				if(!multiMapped && !clipped && readIds.contains(prevReadId)) {
					readIds.remove(prevReadId);
					pw.println(prevLine);
				}
				
				prevReadId = readId;
				prevLine = currentLine;
				multiMapped = false;
				if(!chr.equals("*")) {
					currentReadLength = getReadLength(cigString);
				}
				else
					currentReadLength = -1;
				
				clipped = false;
				if(cigString.contains("S") || cigString.contains("H"))
					clipped = true;
			}
		}
	}
	
	private int getReadLength(String cigString) {
		int readLength = 0;
		int tmpNumber;
		String tmpNumberString = "";
		for(int i = 0; i < cigString.length(); i++) {
			try {
				tmpNumber = Integer.valueOf(cigString.substring(i,i+1));
				tmpNumberString += cigString.charAt(i);
			}
			catch(Exception e) {
				if(tmpNumberString.length() > 0) {
					tmpNumber = Integer.valueOf(tmpNumberString);
					if(cigString.charAt(i) == 'M' || cigString.charAt(i) == '=' || cigString.charAt(i) == 'X') {
						readLength += tmpNumber;
					}
					tmpNumberString = "";
				}
			}
		}
		return readLength;
	}
	
	/**
	 * adds the sequence field to a given sam file
	 * ATTENTION: Ordering of the output will not be consistent with the input ordering!
	 * If you wish a consistent ordering please prebuffer all read sequence positions
	 * @param samFilePath
	 * @param outputFilePath
	 * @param readsFilePath
	 * @param readFormat
	 */
	
	
	public void addSequenceField(String samFilePath, String outputFilePath, String readsFilePath, String readFormat, boolean prebuffer, boolean addQuality) {
		try {
			int bufferSize = 40000000;
			if(prebuffer)
				bufferSize = Integer.MAX_VALUE;
			
			HashMap<String,Long> id2sequencePosition = new HashMap<String,Long>();
			BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(readsFilePath),"r",1024 * 8);
			
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			String currentLine;
			long filePointer;
			
			//write the header
			writeHeader(samFilePath,pw);
			
			String[] splittedLine;
			while((currentLine = braf.getNextLine()) != null) {
				filePointer = braf.getFilePointer();
				if(currentLine.charAt(0) == '>' || currentLine.charAt(0) == '@') {
					if(currentLine.contains(" ")) {
						splittedLine = currentLine.split(" ")[1].split(":");
						if(splittedLine.length == 4) {
							currentLine = currentLine.split(" ")[0] + "/" + splittedLine[0];
						}
					}
					id2sequencePosition.put(currentLine.substring(1),filePointer);
					//skip the sequence
					braf.getNextLine();
					if(readFormat.equals("fastq")) {
						braf.getNextLine();
						braf.getNextLine();
					}
				}
				
				if(id2sequencePosition.size() >= bufferSize) {
					printSequences(id2sequencePosition, samFilePath, readsFilePath, pw, prebuffer, addQuality);
					id2sequencePosition.clear();
				}
			}
			
			printSequences(id2sequencePosition, samFilePath, readsFilePath, pw, prebuffer, addQuality);
			braf.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	
	private void writeHeader(String samFilePath, PrintWriter pw) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
		String currentLine;
		while((currentLine = br.readLine()) != null) {
			if(currentLine.charAt(0) == '@') {
				pw.println(currentLine);
			}
			else
				break;
		}
		br.close();
	}
	
	private void printSequences(HashMap<String,Long> id2sequencePosition, String samFilePath, String readsFilePath, PrintWriter pw, boolean prebuffer, boolean addQuality) throws Exception {
		BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(readsFilePath),"r",1024 * 8);
		BufferedReader samReader = new BufferedReader(new FileReader(new File(samFilePath)));
		String currentLine;
		String[] splittedLine;
		char strand;
		StringBuffer sb = new StringBuffer();
		Pattern tabPattern = Pattern.compile("\t");
		String readId;
		int flag;
		while((currentLine = samReader.readLine()) != null) {
			if(currentLine.charAt(0) == '@')
				continue;
			
			splittedLine = tabPattern.split(currentLine);
			readId = splittedLine[0];
			flag = Integer.valueOf(splittedLine[1]);
			if((flag & (1L << 6)) == 64) {
				readId = splittedLine[0] + "/1";
			}
			else if((flag & (1L << 7)) == 128) {
				readId = splittedLine[0] + "/2";
			}
			
			if(id2sequencePosition.containsKey(readId)) {
				braf.seek(id2sequencePosition.get(readId));
				strand = getStrandFromSamFlag(Integer.valueOf(splittedLine[1]));
				if(strand == '+') {
					splittedLine[9] = braf.getNextLine();
					if(addQuality) {
						braf.getNextLine();
						splittedLine[10] = braf.getNextLine();
					}
				}
				else {
					sb.setLength(0);
					sb.append(braf.getNextLine());
					sb.reverse();
					for(int i = 0; i < sb.length(); i++) {
						sb.setCharAt(i, substitute(sb.charAt(i)));
					}
					splittedLine[9] = sb.toString();
					
					if(addQuality) {
						braf.getNextLine();
						sb.setLength(0);
						sb.append(braf.getNextLine());
						sb.reverse();
						splittedLine[10] = sb.toString();
					}
				}
				
				//add template length information
				//splittedLine[8] = String.valueOf(splittedLine[9].length());
				
				pw.print(splittedLine[0]);
				for(int i = 1; i < splittedLine.length; i++) {
					pw.print("\t" + splittedLine[i]);
				}
				pw.println();
			}
			
			//if we have prebuffered all sequences and cannot find the current id, we will print the line as it is.
			else if(prebuffer) {
				pw.println(currentLine);
			}
		}
		braf.close();
		samReader.close();
	}
	
	private char getStrandFromSamFlag(int samFlag) {
		if(samFlag < 16)
			return '+';
		
		String binaryFlag = Integer.toBinaryString(samFlag);
		if(binaryFlag.charAt(binaryFlag.length() - 5) == '1')
			return '-';
		
		return '+';
			
	}
	
	
	public void addMdFlag(String inputDirPath, String outputFilePath, File[] refFiles, boolean considerClippedRegions, int threads) {
		try {
			
			//map reference name 2 file
			HashMap<String, File> refName2File = new HashMap<String, File>();
			for(File refFile : refFiles) {
				if(refFile.isDirectory() || refFile.getName().indexOf('.') == -1)
					continue;
				
				refName2File.put(refFile.getName().substring(0,refFile.getName().lastIndexOf('.')), refFile);
			}
			
			if(new File(inputDirPath).isFile()) {
				Random random = new Random();
				String tmpOutputDirPath = inputDirPath.substring(0,inputDirPath.lastIndexOf('/')) + "/tmp_mdflag_" + random.nextInt();
				new File(tmpOutputDirPath).mkdirs();
				splitByChromosome(inputDirPath, tmpOutputDirPath);
				inputDirPath = tmpOutputDirPath;
			}
			
			File[] samFiles = new File(inputDirPath).listFiles();
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			String refName;
			File refFile;
			StringBuilder sb = new StringBuilder();
			BufferedReader br;
			ArrayList<Future> futures = new ArrayList<Future>();
			
			BufferedRandomAccessFile rf;
			long fileSize;
			long incRate;
			long prevPosition;
			long currentPosition;
			
			int chunkIndex;
			String tmpOutputFilePath;
			ArrayList<String> tmpPaths = new ArrayList<String>();
			FileConcatenator fileConcatenator = new FileConcatenator();
			
			for(File samFile : samFiles) {
				
				refName = samFile.getName().substring(0,samFile.getName().lastIndexOf('.'));
				refFile = null;
				if(refName2File.containsKey(refName)) {
					refFile = refName2File.get(refName);
				}
				else {
					continue;
				}
				
				//buffer current reference
				sb.setLength(0);
				br = new BufferedReader(new FileReader(refFile));
				br.readLine();
				while(br.ready()) {
					sb.append(br.readLine());
				}
				
				rf = new BufferedRandomAccessFile(samFile,"r",1024);
				fileSize = rf.getChannel().size();
				incRate = fileSize/threads;
				prevPosition = 0;
				currentPosition = 0;
				chunkIndex = 0;
				while(currentPosition < fileSize) {
					currentPosition += incRate;
					rf.seek(currentPosition);
					rf.readLine();
					currentPosition = rf.getFilePointer();
					
					tmpOutputFilePath = inputDirPath + "/" + samFile.getName() + ".mod" + chunkIndex;
					futures.add(executor.submit(new MDFlagGenerator(samFile,sb, tmpOutputFilePath,prevPosition,currentPosition,considerClippedRegions)));
					tmpPaths.add(tmpOutputFilePath);
					prevPosition = currentPosition;
					chunkIndex++;
				}
				rf.close();
				
				for(Future future : futures) {
					future.get();
				}
				futures.clear();
			}
			
			executor.shutdown();
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			
			//concatenate tmp files
			fileConcatenator.concatenate(tmpPaths,outputFilePath);
			for(String tmpPath : tmpPaths) {
				new File(tmpPath).delete();
			}
			
		}
	
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private class MDFlagGenerator extends Thread {
		
		private File samFile;
		private StringBuilder reference;
		private String outputFilePath;
		private long startPos;
		private long stopPos;
		
		private boolean considerClippedRegions;
		
		
		public MDFlagGenerator(File samFile, StringBuilder reference, String outputFilePath, long startPos, long stopPos, boolean considerClippedRegions) {
			this.samFile = samFile;
			this.reference = reference;
			this.outputFilePath = outputFilePath;
			this.startPos = startPos;
			this.stopPos = stopPos;
			
			this.considerClippedRegions = considerClippedRegions;
		}
		
		public void run() {
			try {
				BufferedRandomAccessFile br = new BufferedRandomAccessFile(samFile,"r",1024 * 8);
				BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.outputFilePath)));
				String currentLine;
				StringTokenizer st;
				String readId;
				char strand;
				String chr = null;
				int start;
				String cigar;
				String readSequence;
				String completeReadSequence;
				String refSequence;
				StringBuffer sb = new StringBuffer();
				ArrayList<String> clippings = new ArrayList<String>();
		    	Pattern softclippedPattern = Pattern.compile("[0-9]+[S|H]");
		    	Pattern matchingPattern = Pattern.compile("[0-9]+[M]");
		    	Pattern intronPattern = Pattern.compile("[0-9]+[N|D]");
		    	ArrayList<Integer> matches = new ArrayList<Integer>();
		    	ArrayList<Integer> gaps = new ArrayList<Integer>();
				Matcher matcher;
				boolean softClippedAtTheStart;
				boolean softClippedAtTheEnd;
				int clippingLengthStart;
				int clippingLengthEnd;
				
				String nmField;
				String mdField;
				String tmpMdField;
				int currentMatchCount;
				int firstMatchCount;
				int lastMatchCount;
				int firstMatchOffset;
				int lastMatchOffset;
				
				int currentMismatchCount;
				int trueMismatchCount;
				int skippedReads = 0;
				br.seek(this.startPos);
				while((currentLine = br.getNextLine()) != null) {
					
					st = new StringTokenizer(currentLine, "\t");
					readId = st.nextToken();
					strand = getStrandFromSamFlag(Integer.valueOf(st.nextToken()));
					chr = st.nextToken();
					
					if(chr.equals("*"))
						continue;
					
					start = Integer.valueOf(st.nextToken());
					st.nextToken();
					cigar = st.nextToken();
					st.nextToken();
					st.nextToken();
					st.nextToken();
					completeReadSequence = st.nextToken();
					
					nmField = null;
					while(st.hasMoreTokens()) {
						nmField = st.nextToken();
						if(nmField.length() > 1 && nmField.substring(0,2).equals("NM"))
							break;
					}
					st = new StringTokenizer(nmField, ":");
					st.nextToken();
					st.nextToken();
					trueMismatchCount = Integer.valueOf(st.nextToken());
					
					/*
					 * check for soft clipped alignment and trim the read sequence
					 */
					clippings.clear();
					matcher = softclippedPattern.matcher(cigar);
					softClippedAtTheStart = false;
					softClippedAtTheEnd = false;
					while(matcher.find()) {
						if(matcher.start() == 0) {
							softClippedAtTheStart = true;
						}
						else
							softClippedAtTheEnd = true;
						
						clippings.add(matcher.group());
					}
					
					clippingLengthStart = 0;
					clippingLengthEnd = 0;
					if(softClippedAtTheStart) {
						clippingLengthStart = Integer.valueOf(clippings.get(0).substring(0,clippings.get(0).length() - 1));
					}
					
					if(softClippedAtTheEnd) {
						clippingLengthEnd = Integer.valueOf(clippings.get(clippings.size() - 1).substring(0,clippings.get(clippings.size() - 1).length() - 1));
					}
					readSequence = completeReadSequence.substring(clippingLengthStart,completeReadSequence.length() - clippingLengthEnd);
					
					
					/*
					 * check for gaps in the alignment
					 */
					matches.clear();
					gaps.clear();
					matcher = intronPattern.matcher(cigar);
					while(matcher.find()) {
						gaps.add(Integer.valueOf(matcher.group().substring(0,matcher.group().length() - 1)));
					}
					
					if(!gaps.isEmpty()) {
						matcher = matchingPattern.matcher(cigar);
						while(matcher.find()) {
							matches.add(Integer.valueOf(matcher.group().substring(0,matcher.group().length() - 1)));
						}
					}

					/*
					 * now determine the md flag
					 */
					if(gaps.isEmpty())
						refSequence = this.reference.substring(start - 1,start + readSequence.length() - 1).toUpperCase();
					else {
						refSequence = "";
						for(int i = 0; i < gaps.size(); i++) {
							refSequence += this.reference.substring(start - 1, start + matches.get(i) - 1).toUpperCase();
							start += matches.get(i) + gaps.get(i);
						}
						refSequence += this.reference.substring(start - 1, start + matches.get(matches.size() - 1) - 1).toUpperCase();
					}
					
					mdField = "";
					currentMatchCount = 0;
					firstMatchCount = -1;
					firstMatchOffset = -1;
					currentMismatchCount = 0;
					
					for(int i = 0; i < readSequence.length(); i++) {
						if(readSequence.charAt(i) == refSequence.charAt(i)) {
							currentMatchCount++;
						}
						else {
							mdField += currentMatchCount;
							mdField += refSequence.charAt(i);
							
							if(firstMatchCount == -1) {
								firstMatchCount = currentMatchCount;
								firstMatchOffset = mdField.length() - 1;
							}
							currentMatchCount = 0;
							currentMismatchCount++;
						}
					}
					
					lastMatchCount = currentMatchCount;
					lastMatchOffset = String.valueOf(currentMatchCount).length();
					mdField += currentMatchCount;
					
					if(firstMatchCount == -1) {
						firstMatchCount = currentMatchCount;
						firstMatchOffset = mdField.length();
					}
					
					//consistency check
					if(currentMismatchCount == trueMismatchCount) {
						
						//also determine the md field for clipped regions.
						if(this.considerClippedRegions) {
							if(softClippedAtTheStart) {
								if(start - clippingLengthStart - 1 >= 0) {
									readSequence = completeReadSequence.substring(0,clippingLengthStart);
									refSequence = this.reference.substring(start - clippingLengthStart - 1,start - 1).toUpperCase();
									currentMatchCount = 0;
									currentMismatchCount = 0;
									tmpMdField = "";
									for(int i = 0; i < readSequence.length(); i++) {
										if(readSequence.charAt(i) == refSequence.charAt(i)) {
											currentMatchCount++;
										}
										else {
											tmpMdField += currentMatchCount;
											tmpMdField += refSequence.charAt(i);
											currentMatchCount = 0;
											currentMismatchCount++;
										}
									}
									currentMatchCount += firstMatchCount;
									tmpMdField += currentMatchCount;
									
									mdField = tmpMdField + mdField.substring(firstMatchOffset);
								}
								else {
									for(int i = 0; i < clippingLengthStart; i++) {
										mdField = "0N" + mdField;
									}
								}
								
							}
							
							
							if(softClippedAtTheEnd) {
								if(start + completeReadSequence.length() - 1 <= this.reference.length()) {
									
									readSequence = completeReadSequence.substring(completeReadSequence.length() - clippingLengthEnd,completeReadSequence.length());
									refSequence = this.reference.substring(start + completeReadSequence.length() - clippingLengthEnd - 1, start + completeReadSequence.length() - 1).toUpperCase();
									
									currentMatchCount = 0;
									firstMatchCount = -1;
									currentMismatchCount = 0;
									tmpMdField = "";
									for(int i = 0; i < readSequence.length(); i++) {
										if(readSequence.charAt(i) == refSequence.charAt(i)) {
											currentMatchCount++;
										}
										else {
											
											if(firstMatchCount == -1) {
												firstMatchCount = currentMatchCount;
												currentMatchCount += lastMatchCount;
											}
											
											tmpMdField += currentMatchCount;
											tmpMdField += refSequence.charAt(i);
											
											
											currentMatchCount = 0;
											currentMismatchCount++;
										}
									}
									
									if(firstMatchCount == -1) {
										currentMatchCount += lastMatchCount;
									}
									
									tmpMdField += currentMatchCount;
									
									
									mdField = mdField.substring(0,mdField.length() - lastMatchOffset) + tmpMdField;
								}
								
								else {
									for(int i = 0; i < clippingLengthEnd; i++) {
										mdField += "N0";
									}
								}
								
								
								
							}
							
						}
						
						bw.write(currentLine + "\tMD:Z:" + mdField);
						bw.newLine();
					}
					
					else {
						skippedReads++;
					}
					
					if(br.getFilePointer() == this.stopPos) {
						break;
					}
				}
				br.close();
				bw.close();
				if(skippedReads > 0) {
					//System.err.println(String.format("Warning @ md flag generation: Skipped %s reads. Reason: Annotated mismatch count in nmField differed from generated mdField.", skippedReads));
				}
			}
			
			catch(Exception e) {
				e.printStackTrace();
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
		
		private char substitute(char n) {
			if(n == 'A' || n == 'a') return 'T';
			if(n == 'T' || n == 't') return 'A';
			if(n == 'C' || n == 'c') return 'G';
			if(n == 'G' || n == 'g') return 'C';
			else return 'N';
		}
	}
	
	private class FileConcatenator {
		
		
		public FileConcatenator() {
			
		}
		
		public void concatenate(ArrayList<String> filePaths, String outputPath) {
			try {
				FileOutputStream fos = new FileOutputStream(outputPath,true);
				FileChannel writeChannel = fos.getChannel();
				RandomAccessFile rf;
				FileChannel readChannel;
				long currentChannelSize;
				long transferedBytes;
				for(String filePath : filePaths) {
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
		
		private int getIndexOfLargestFile(ArrayList<String> filePaths) {
			int index = -1;
			String path;
			long currentSize;
			long largestSize = Long.MIN_VALUE;
			for(int i = 0; i < filePaths.size(); i++) {
				path = filePaths.get(i);
				currentSize = new File(path).length();
				if(currentSize > largestSize) {
					index = i;
					largestSize = currentSize;
				}
			}
			return index;
		}
	}
	
	/**
	 * generates a fasta file containing reads not aligned in the given sam file
	 * @param inputFilePath
	 * @param outputFilePath
	 */
	
	public void getUnalignedReads(String inputFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			String readId;
			String chr;
			String sequence;
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				readId = splittedLine[0];
				chr = splittedLine[2];
				sequence = splittedLine[9];
				
				if(chr.equals("*")) {
					bw.write(">" + readId);
					bw.newLine();
					bw.write(sequence);
					bw.newLine();
				}
			}
			br.close();
			bw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	//sam file should NOT contain any multi mappings
	public void getMicrobialContentFromUnconcatenatedGenomes(String samFilePath, String indexFolderPath, String refseqCatalogFilePath, String customSpeciesNamesFilePath, boolean mergeContigs) {
		try {
			File[] indexFiles = new File(indexFolderPath).listFiles();
			BufferedReader indexReader;
			String currentIndexName;
			String currentLine;
			StringTokenizer st;
			String indexedId;
			String microbeId;
			String trivialName;
			int genomeLength;
			//hash index files
			HashMap<String,Microbe> idx2microbes = new HashMap<String,Microbe>();
			HashMap<String,Microbe> microbe2mergedContigs = new HashMap<String,Microbe>();
			HashMap<String,String> accession2name = new HashMap<String,String>();
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			for(File f : indexFiles) {
				if(f.getName().substring(f.getName().lastIndexOf('.') + 1).equals("idx")) {
					currentIndexName = f.getName().substring(0,f.getName().lastIndexOf('.'));
					indexReader = new BufferedReader(new FileReader(f));
					while(indexReader.ready()) {
						currentLine = indexReader.readLine();
						st = new StringTokenizer(currentLine,"\t");
						indexedId = st.nextToken();
						microbeId = st.nextToken();
						genomeLength = Integer.valueOf(st.nextToken());
						splittedId = barPattern.split(microbeId);
						trivialName = "NA";
						if(st.hasMoreTokens()) {
							trivialName = st.nextToken();
						}
						if(mergeContigs && microbeId.contains("NZ_")) {
							microbe2mergedContigs.put(splittedId[3].substring(0,9),new Microbe(splittedId[3].substring(0,9) + "000000",null,genomeLength,true));
							accession2name.put(splittedId[3].substring(0,9) + "000000", "NA");
						}
						
						else if(mergeContigs && splittedId.length >= 4 && !splittedId[3].contains("_") && splittedId[3].charAt(12) == '.') {
							microbe2mergedContigs.put(splittedId[3].substring(0,6),new Microbe(splittedId[3].substring(0,6) + "000000",null,genomeLength,true));
							accession2name.put(splittedId[3].substring(0,6) + "000000", "NA");
						}
						
						
						//here we have genbank and refseq ids, we remove the genbank info
						if(splittedId.length >= 4) {
							idx2microbes.put(indexedId, new Microbe(splittedId[3],null,genomeLength,false));
							accession2name.put(splittedId[3],trivialName);
						}
						
						//here we only have chromosome names from the reference (human,mouse or yeast)
						else {
							idx2microbes.put(microbeId, new Microbe(microbeId,null,genomeLength,false));
							accession2name.put(microbeId,"NA");
						}
					}
					indexReader.close();
				}
			}
			
			//map to every accession number the species name
			//we prefer trivial names from index files, so do not overwrite already found names
			Pattern tabPattern = Pattern.compile("\t");
			String[] splittedLine;
			if(refseqCatalogFilePath != null) {
				BufferedReader catalogReader = new BufferedReader(new FileReader(new File(refseqCatalogFilePath)));
				while(catalogReader.ready()) {
					currentLine = catalogReader.readLine();
					splittedLine = tabPattern.split(currentLine);
					if(accession2name.containsKey(splittedLine[2]) && accession2name.get(splittedLine[2]).equals("NA"))
						accession2name.put(splittedLine[2],splittedLine[1]);
					if(mergeContigs && splittedLine[2].contains("NZ_") && accession2name.containsKey(splittedLine[2].substring(0,9) + "000000"))
						accession2name.put(splittedLine[2].substring(0,9) + "000000",splittedLine[1]);
					if(mergeContigs && splittedLine[2].length() > 12 && splittedLine[2].charAt(2) != '_' && splittedLine[2].charAt(12) == '.' && accession2name.containsKey(splittedLine[2].substring(0,6) + "000000")) {
						accession2name.put(splittedLine[2].substring(0,6) + "000000",splittedLine[1]);
					}
				}
				catalogReader.close();
			}
			
			if(customSpeciesNamesFilePath != null) {
				BufferedReader br = new BufferedReader(new FileReader(new File(customSpeciesNamesFilePath)));
				while((currentLine = br.readLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					if(accession2name.containsKey(splittedLine[0]) && accession2name.get(splittedLine[0]).equals("NA")) {
						accession2name.put(splittedLine[0],splittedLine[1]);
					}
				}
				br.close();
			}
			
			//species without trivialname get their accession numbers as trivial names
			for(String speciesName : accession2name.keySet()) {
				if(accession2name.get(speciesName).equals("NA"))
					accession2name.put(speciesName, speciesName);
			}
			
			//process sam file
			int overallContainedReads = 0;
			BufferedReader samReader = new BufferedReader(new FileReader(new File(samFilePath)));
			Microbe currentMicrobe;
			while(samReader.ready()) {
				currentLine = samReader.readLine();
				st = new StringTokenizer(currentLine,"\t");
				//read id bitflag
				st.nextToken();
				st.nextToken();
				currentIndexName = st.nextToken();
				if(idx2microbes.containsKey(currentIndexName)) {
					currentMicrobe = idx2microbes.get(currentIndexName);
					microbeId = currentMicrobe.getId();
					
					if(mergeContigs && microbeId.contains("NZ_")) {
						if(microbe2mergedContigs.containsKey(microbeId.substring(0,9))) {
							microbe2mergedContigs.get(microbeId.substring(0,9)).incrementContainedReads();
						}
					}
					else if(mergeContigs && microbeId.length() >= 6 && microbeId.charAt(2) != '_' && microbe2mergedContigs.containsKey(microbeId.substring(0,6))) {
						microbe2mergedContigs.get(microbeId.substring(0,6)).incrementContainedReads();
					}
					
					else
						currentMicrobe.incrementContainedReads();
					
					overallContainedReads++;
				}
			}
			samReader.close();
			
			
			//output organisms with at least 1 contained read
			HashMap<String,Integer> microbe2genomeLength = new HashMap<String,Integer>();
			TreeMap<Integer,ArrayList<String>> sortedMappings = new TreeMap<Integer,ArrayList<String>>();
			for(Microbe microbe : idx2microbes.values()) {
				if(microbe.containedReads > 0) {
					if(sortedMappings.containsKey(microbe.getContainedReads()))
						sortedMappings.get(microbe.getContainedReads()).add(microbe.getId());
					else {
						ArrayList<String> microbeIds = new ArrayList<String>();
						microbeIds.add(microbe.getId());
						sortedMappings.put(microbe.getContainedReads(),microbeIds);
					}
					microbe2genomeLength.put(microbe.getId(), microbe.getLength());
					
				}
				microbeId = microbe.getId();
				if(mergeContigs && microbe.getId().contains("NZ_")) {
					if(microbe2genomeLength.containsKey(microbeId.substring(0,9))) {
						microbe2genomeLength.put(microbeId.substring(0,9),microbe2genomeLength.get(microbeId.substring(0,9)) + microbe.getLength());
					}
					else
						microbe2genomeLength.put(microbeId.substring(0,9),microbe.getLength());
				}
				else if(mergeContigs && microbeId.length() > 6 && microbe2mergedContigs.containsKey(microbeId.substring(0,6))) {
					
					if(microbe2genomeLength.containsKey(microbeId.substring(0,6))) {
						microbe2genomeLength.put(microbeId.substring(0,6),microbe2genomeLength.get(microbeId.substring(0,6)) + microbe.getLength());
					}
					else
						microbe2genomeLength.put(microbeId.substring(0,6),microbe.getLength());
				}
			}
			
			
			
			//in case we have merged contigs we add those contigs now to the treemap
			for(Microbe microbe : microbe2mergedContigs.values()) {
				if(microbe.getContainedReads() > 0) {
					if(sortedMappings.containsKey(microbe.getContainedReads()))
						sortedMappings.get(microbe.getContainedReads()).add(microbe.getId());
					else {
						ArrayList<String> microbeIds = new ArrayList<String>();
						microbeIds.add(microbe.getId());
						sortedMappings.put(microbe.getContainedReads(),microbeIds);
					}
				}
			}
			
			for(int containedReads : sortedMappings.keySet()) {
				for(String microbe : sortedMappings.get(containedReads)) {
					if(mergeContigs && microbe.contains("NZ_"))
						System.out.println(String.format("%s\t%s\t%s\t%s\t%s",microbe,accession2name.get(microbe),containedReads,((double)containedReads/(double)overallContainedReads),microbe2genomeLength.get(microbe.substring(0,9))));
					else if(mergeContigs && microbe.length() > 6 && microbe2mergedContigs.containsKey(microbe.substring(0,6)))
						System.out.println(String.format("%s\t%s\t%s\t%s\t%s",microbe,accession2name.get(microbe),containedReads,((double)containedReads/(double)overallContainedReads),microbe2genomeLength.get(microbe.substring(0,6))));
					else
						System.out.println(String.format("%s\t%s\t%s\t%s\t%s",microbe,accession2name.get(microbe),containedReads,((double)containedReads/(double)overallContainedReads),microbe2genomeLength.get(microbe)));
				}
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	//sam file should NOT contain any multi mappings
	public void getMicrobialContentFromConcatenatedGenomes(String samFilePath, String indexFolderPath, String refseqCatalogFilePath, String customSpeciesNamesFilePath, boolean mergeContigs, boolean useMDflag, HashSet<String> referenceChromosomes) {
		try {
			File[] indexFiles = new File(indexFolderPath).listFiles();
			BufferedReader indexReader;
			String currentIndexName;
			String currentLine;
			StringTokenizer st;
			String microbeId;
			String trivialName;
			int microbeStart;
			int microbeEnd;
			Date date;
			//hash index files
			HashMap<String,TreeMap<Integer,Microbe>> idx2microbes = new HashMap<String,TreeMap<Integer,Microbe>>();
			HashMap<String,Microbe> microbe2mergedContigs = new HashMap<String,Microbe>();
			HashMap<String,String> accession2name = new HashMap<String,String>();
			TreeMap<Integer,Microbe> microbe2start;
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			date = new Date();
			System.err.println(String.format("[%s]\tDetermining read counts of species contained in the given sam file",date.toLocaleString()));
			
			for(File f : indexFiles) {
				if(f.getName().substring(f.getName().lastIndexOf('.') + 1).equals("idx")) {
					currentIndexName = f.getName().substring(0,f.getName().lastIndexOf('.'));
					
					microbe2start = new TreeMap<Integer,Microbe>();
					idx2microbes.put(currentIndexName, microbe2start);
					indexReader = new BufferedReader(new FileReader(f));
					while(indexReader.ready()) {
						currentLine = indexReader.readLine();
						st = new StringTokenizer(currentLine,"\t");
						microbeId = st.nextToken();
						splittedId = barPattern.split(microbeId);
						microbeStart = Integer.valueOf(st.nextToken());
						microbeEnd = Integer.valueOf(st.nextToken());
						trivialName = "NA";
						if(st.hasMoreTokens()) {
							trivialName = st.nextToken();
						}
						if(mergeContigs && microbeId.contains("NZ_")) {
							microbe2mergedContigs.put(splittedId[3].substring(0,9),new Microbe(splittedId[3].substring(0,9) + "000000",0,0,true));
							accession2name.put(splittedId[3].substring(0,9) + "000000", "NA");
						}
						
						else if(mergeContigs && splittedId.length >= 4 && !splittedId[3].contains("_") && splittedId[3].charAt(12) == '.') {
							microbe2mergedContigs.put(splittedId[3].substring(0,6),new Microbe(splittedId[3].substring(0,6) + "000000",0,0,true));
							accession2name.put(splittedId[3].substring(0,6) + "000000", "NA");
						}
						
						
						//here we have genbank and refseq ids, we remove the genbank info
						if(splittedId.length >= 4) {
							microbe2start.put(microbeStart, new Microbe(splittedId[3],microbeStart,microbeEnd,false));
							accession2name.put(splittedId[3],trivialName);
						}
						
						//here we only have chromosome names from the reference (human,mouse or yeast)
						else {
							microbe2start.put(microbeStart, new Microbe(microbeId,microbeStart,microbeEnd,false));
							accession2name.put(microbeId,"NA");
						}
					}
					indexReader.close();
				}
			}
			
			//map to every accession number the species name
			//we prefer trivial names from index files, so do not overwrite already found names
			Pattern tabPattern = Pattern.compile("\t");
			String[] splittedLine;
			if(refseqCatalogFilePath != null) {
				BufferedReader catalogReader = new BufferedReader(new FileReader(new File(refseqCatalogFilePath)));
				while(catalogReader.ready()) {
					currentLine = catalogReader.readLine();
					splittedLine = tabPattern.split(currentLine);
					if(accession2name.containsKey(splittedLine[2]) && accession2name.get(splittedLine[2]).equals("NA"))
						accession2name.put(splittedLine[2],splittedLine[1]);
					if(mergeContigs && splittedLine[2].contains("NZ_") && accession2name.containsKey(splittedLine[2].substring(0,9) + "000000"))
						accession2name.put(splittedLine[2].substring(0,9) + "000000",splittedLine[1]);
					if(mergeContigs && splittedLine[2].length() > 12 && splittedLine[2].charAt(2) != '_' && splittedLine[2].charAt(12) == '.' && accession2name.containsKey(splittedLine[2].substring(0,6) + "000000")) {
						accession2name.put(splittedLine[2].substring(0,6) + "000000",splittedLine[1]);
					}
				}
				catalogReader.close();
			}
			
			if(customSpeciesNamesFilePath != null) {
				BufferedReader br = new BufferedReader(new FileReader(new File(customSpeciesNamesFilePath)));
				while((currentLine = br.readLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					if(accession2name.containsKey(splittedLine[0]) && accession2name.get(splittedLine[0]).equals("NA")) {
						accession2name.put(splittedLine[0],splittedLine[1]);
					}
				}
				br.close();
			}
			
			//species without trivialname get their accession numbers as trivial names
			for(String speciesName : accession2name.keySet()) {
				if(accession2name.get(speciesName).equals("NA"))
					accession2name.put(speciesName, speciesName);
			}
			
			//process sam file
			int overallContainedReads = 0;
			BufferedReader samReader = new BufferedReader(new FileReader(new File(samFilePath)));
			Microbe currentMicrobe;
			while(samReader.ready()) {
				currentLine = samReader.readLine();
				st = new StringTokenizer(currentLine,"\t");
				//read id bitflag
				st.nextToken();
				st.nextToken();
				currentIndexName = st.nextToken();
				if(idx2microbes.containsKey(currentIndexName)) {
					
					microbe2start = idx2microbes.get(currentIndexName);
					microbeStart = Integer.valueOf(st.nextToken());
					currentMicrobe = microbe2start.get(microbe2start.floorKey(microbeStart));
					microbeId = currentMicrobe.getId();
					
					if(mergeContigs && microbeId.contains("NZ_")) {
						if(microbe2mergedContigs.containsKey(microbeId.substring(0,9))) {
							microbe2mergedContigs.get(microbeId.substring(0,9)).incrementContainedReads();
							microbe2mergedContigs.get(microbeId.substring(0,9)).addCoveredPosition(microbeStart);
						}
					}
					else if(mergeContigs && microbeId.length() >= 6 && microbeId.charAt(2) != '_' && microbe2mergedContigs.containsKey(microbeId.substring(0,6))) {
						microbe2mergedContigs.get(microbeId.substring(0,6)).incrementContainedReads();
						microbe2mergedContigs.get(microbeId.substring(0,6)).addCoveredPosition(microbeStart);
					}
					
					else {
						microbe2start.get(microbe2start.floorKey(microbeStart)).incrementContainedReads();
						microbe2start.get(microbe2start.floorKey(microbeStart)).addCoveredPosition(microbeStart);
					}
					
					overallContainedReads++;
				}
			}
			samReader.close();
			
			
			//output organisms with at least 1 contained read
			HashMap<String,Integer> microbe2genomeLength = new HashMap<String,Integer>();
			HashMap<String,Integer> microbe2coveredPositions = new HashMap<String,Integer>();
			TreeMap<Integer,ArrayList<String>> sortedMappings = new TreeMap<Integer,ArrayList<String>>();
			for(TreeMap<Integer,Microbe> tmpMap : idx2microbes.values()) {
				for(Microbe microbe : tmpMap.values()) {
					if(microbe.containedReads > 0) {
						if(sortedMappings.containsKey(microbe.getContainedReads()))
							sortedMappings.get(microbe.getContainedReads()).add(microbe.getId());
						else {
							ArrayList<String> microbeIds = new ArrayList<String>();
							microbeIds.add(microbe.getId());
							sortedMappings.put(microbe.getContainedReads(),microbeIds);
						}
						microbe2genomeLength.put(microbe.getId(), microbe.getLength());
						microbe2coveredPositions.put(microbe.getId(), microbe.getCoveredPositions().size());
						
					}
					microbeId = microbe.getId();
					if(mergeContigs && microbe.getId().contains("NZ_")) {
						if(microbe2genomeLength.containsKey(microbeId.substring(0,9))) {
							microbe2genomeLength.put(microbeId.substring(0,9),microbe2genomeLength.get(microbeId.substring(0,9)) + microbe.getLength());
						}
						else
							microbe2genomeLength.put(microbeId.substring(0,9),microbe.getLength());
					}
					else if(mergeContigs && microbeId.length() > 6 && microbe2mergedContigs.containsKey(microbeId.substring(0,6))) {
						
						if(microbe2genomeLength.containsKey(microbeId.substring(0,6))) {
							microbe2genomeLength.put(microbeId.substring(0,6),microbe2genomeLength.get(microbeId.substring(0,6)) + microbe.getLength());
						}
						else
							microbe2genomeLength.put(microbeId.substring(0,6),microbe.getLength());
					}
				}
			}
			
			
			//in case we have merged contigs we add those contigs now to the treemap
			for(Microbe microbe : microbe2mergedContigs.values()) {
				if(microbe.getContainedReads() > 0) {
					if(sortedMappings.containsKey(microbe.getContainedReads()))
						sortedMappings.get(microbe.getContainedReads()).add(microbe.getId());
					else {
						ArrayList<String> microbeIds = new ArrayList<String>();
						microbeIds.add(microbe.getId());
						sortedMappings.put(microbe.getContainedReads(),microbeIds);
					}
					
					
					microbe2coveredPositions.put(microbe.getId(), microbe.getCoveredPositions().size());
				}
			}
			
			date = new Date();
			System.err.println(String.format("[%s]\tCalculating confidence values for identified species",date.toLocaleString()));
			
			/**
			 * Now get the confidence values
			 */
			HashMap<String,String> id2trivialName = new HashMap<String,String>();
			for(int containedReads : sortedMappings.keySet()) {
				for(String microbe : sortedMappings.get(containedReads)) {
					id2trivialName.put(microbe,accession2name.get(microbe));
				}
			}
			HashMap<String,Double> species2confidence = new Statistic().generateConfidenceScoresFromConcatenatedGenomes(samFilePath, null, indexFolderPath, false, id2trivialName);
			
			
			date = new Date();
			System.err.println(String.format("[%s]\tCalculating normalized Jenson Shannon divergence values",date.toLocaleString()));
			
			/**
			 * Now get the mismatch rates and calculate the normalized JS values
			 */
			HashMap<String,HashMap<Integer,Double>> species2errorRate =  new Statistic().generateErrorRateTableFromConcatenatedGenomes(samFilePath, null,indexFolderPath, useMDflag, referenceChromosomes,id2trivialName);
			HashMap<String,Double> species2divergence = new HashMap<String,Double>();
			HashMap<String,Double> species2normalizedDivergence = new HashMap<String,Double>();
			
			//obtaining the reference species
			String reference;
			
			if(referenceChromosomes.isEmpty()) {
				reference = getMostAbundantSpecies(species2errorRate);
				date = new Date();
				System.err.println(String.format("[%s]\tThe most abundant species is taken as the reference species for the Jenson-Shannon Divergence calcluation",date.toLocaleString()));
			}
			else if(referenceChromosomes.size() == 1) {
				reference = referenceChromosomes.iterator().next();
				if(!species2errorRate.containsKey(reference)) {
					reference = getMostAbundantSpecies(species2errorRate);
					date = new Date();
					System.err.println(String.format("[%s]\tWARNING: Could not find the given species (%s) in the set of identified microbes",date.toLocaleString(),reference));
					System.err.println(String.format("[%s]\tTherefore, the most abundant species(%s) is taken as the reference species for the Jensen-Shannon Divergence calcluation",date.toLocaleString(),reference));
				}
				date = new Date();
				System.err.println(String.format("[%s]\tThe user defined species(%s) is taken as the reference species for the Jensen-Shannon divergence calcluation",date.toLocaleString(),reference));
			}
			else {
				reference = "Reference";
				date = new Date();
				System.err.println(String.format("[%s]\tThe user defined chromosome names were taken as reference for the Jensen-Shannon Divergence calcluation",date.toLocaleString()));
			}
			
			HashMap<Integer,Double> Q = species2errorRate.get(reference);
			double readCount = 0.0;
			for(double i : Q.values()) {
				readCount += i;
			}
			for(String species : species2errorRate.keySet()) {
				species2divergence.put(species,getKullbackLeiblerDivergence(species2errorRate.get(species),Q,readCount, 1.0));
				species2normalizedDivergence.put(species,getJensenShannonDivergence(species2errorRate.get(species),Q,readCount, 1.0));
			}
			
			
			System.out.println("id\ttrivial_name\tgenome_length\tread_count\tcoverage\tconfidence\tSqrt_of_JS_divergence");
			for(int containedReads : sortedMappings.keySet()) {
				for(String microbe : sortedMappings.get(containedReads)) {
					if(mergeContigs && microbe.contains("NZ_"))
						System.out.println(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s",microbe,accession2name.get(microbe),microbe2genomeLength.get(microbe.substring(0,9)),containedReads,(double)microbe2coveredPositions.get(microbe.substring(0,9))/(double)microbe2genomeLength.get(microbe.substring(0,9)),species2confidence.get(microbe),species2normalizedDivergence.get(microbe)));
					else if(mergeContigs && microbe.length() > 6 && microbe2mergedContigs.containsKey(microbe.substring(0,6)))
						System.out.println(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s",microbe,accession2name.get(microbe),microbe2genomeLength.get(microbe.substring(0,6)),containedReads,(double)microbe2coveredPositions.get(microbe.substring(0,6))/(double)microbe2genomeLength.get(microbe.substring(0,6)),species2confidence.get(microbe),species2normalizedDivergence.get(microbe)));
					else {
						System.out.println(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s",microbe,accession2name.get(microbe),microbe2genomeLength.get(microbe),containedReads,(double)microbe2coveredPositions.get(microbe)/(double)microbe2genomeLength.get(microbe),species2confidence.get(microbe),species2normalizedDivergence.get(microbe)));
					}
				}
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static String getMostAbundantSpecies(HashMap<String,HashMap<Integer,Double>> species2errorRate) {
		String mostAbundantSpecies = "";
		int abundance;
		int maxAbundance = Integer.MIN_VALUE;
		HashMap<Integer,Double> errorRates;
		for(String species : species2errorRate.keySet()) {
			abundance = 0;
			errorRates = species2errorRate.get(species);
			for(int i : errorRates.keySet()) {
				abundance += errorRates.get(i);
			}
			if(abundance > maxAbundance) {
				mostAbundantSpecies = species;
				maxAbundance = abundance;
			}
		}
		
		return mostAbundantSpecies;
	}
	
	/**
	 * 
	 * @param P -> species under consideration
	 * @param Q -> reference species
	 * @return D_kl
	 */
	private static double getKullbackLeiblerDivergence(HashMap<Integer,Double> P,HashMap<Integer,Double> Q, double readCountQ,  double pseudoCount) {
		double readCountP = 0.0;
		for(int i : P.keySet()) {
			readCountP += P.get(i);
		}
		readCountP += (Q.keySet().size() * pseudoCount); 
		readCountQ += (Q.keySet().size() * pseudoCount);
		
		double divergence = 0.0;
		double p_i;
		double q_i;
		for(int i : Q.keySet()) {
			p_i = pseudoCount;
			if(P.containsKey(i))
				p_i += P.get(i);
			
			q_i = Q.get(i) + pseudoCount;
			
			divergence += (Math.log((p_i/readCountP)/(q_i/readCountQ)) * (p_i/readCountP));
		}
			
			
		return divergence;
	}
	
	private static double getJensenShannonDivergence(HashMap<Integer,Double> P,HashMap<Integer,Double> Q,double readCountQ, double pseudoCount) {
		double readCountP = 0.0;
		for(int i : P.keySet()) {
			readCountP += P.get(i);
		}
		readCountP += (Q.keySet().size() * pseudoCount); 
		readCountQ += (Q.keySet().size() * pseudoCount);
		
		
		HashMap<Integer,Double> M = new HashMap<Integer,Double>();
		HashMap<Integer,Double> P_pseudo_counts = new HashMap<Integer,Double>();
		HashMap<Integer,Double> Q_pseudo_counts = new HashMap<Integer,Double>();
		double p_i;
		double q_i;
		double readCountM = 0.0;
		for(int i : Q.keySet()) {
			q_i = Q.get(i) + pseudoCount;
			p_i = pseudoCount;
			if(P.containsKey(i))
				p_i += P.get(i);
			
			M.put(i, 0.5 * ((q_i/readCountQ) + (p_i/readCountP)));
			
			//M.put(i,0.5 * (q_i + p_i));
			//readCountM += 0.5 * (q_i + p_i);
			
			Q_pseudo_counts.put(i, q_i);
			P_pseudo_counts.put(i,p_i);
		}
		
		double divergence = (0.5 * getKullbackLeiblerDivergence(P_pseudo_counts,M,1,0)) + (0.5 * getKullbackLeiblerDivergence(Q_pseudo_counts,M,1,0));
		return Math.sqrt(divergence/Math.log(2));
	}
	
	
	//sam file should NOT contain any multi mappings
	public void separateMappingBySpecies(String samFilePath, String indexFolderPath, String speciesListFilePath, String outputDirPath,boolean mergeContigs,boolean trivialFileNames) {
		try {
			File[] indexFiles = new File(indexFolderPath).listFiles();
			BufferedReader indexReader;
			String currentIndexName;
			String currentLine;
			StringTokenizer st;
			String microbeId;
			String trivialName;
			int microbeStart;
			int microbeEnd;
			//hash index files
			HashMap<String,TreeMap<Integer,Microbe>> idx2microbes = new HashMap<String,TreeMap<Integer,Microbe>>();
			HashMap<String,Microbe> microbe2mergedContigs = new HashMap<String,Microbe>();
			TreeMap<Integer,Microbe> microbe2start;
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			for(File f : indexFiles) {
				if(f.getName().substring(f.getName().lastIndexOf('.') + 1).equals("idx")) {
					currentIndexName = f.getName().substring(0,f.getName().lastIndexOf('.'));
					
					microbe2start = new TreeMap<Integer,Microbe>();
					idx2microbes.put(currentIndexName, microbe2start);
					indexReader = new BufferedReader(new FileReader(f));
					while(indexReader.ready()) {
						currentLine = indexReader.readLine();
						st = new StringTokenizer(currentLine,"\t");
						microbeId = st.nextToken();
						splittedId = barPattern.split(microbeId);
						microbeStart = Integer.valueOf(st.nextToken());
						microbeEnd = Integer.valueOf(st.nextToken());
						trivialName = null;
						if(st.hasMoreTokens())
							trivialName = st.nextToken();
						//here we have genbank and refseq ids, we remove genbank info
						if(splittedId.length >= 4) {
							microbe2start.put(microbeStart, new Microbe(splittedId[3],trivialName,0,false));
							
						}
						//here we only have chromosome names from the reference (human,mouse or yeast)
						else {
							microbe2start.put(microbeStart, new Microbe(microbeId,trivialName,0,false));
							
						}
						
						
						if(mergeContigs && microbeId.contains("NZ_")) {
							microbe2mergedContigs.put(splittedId[3].substring(0,9),new Microbe(splittedId[3].substring(0,9) + "000000",trivialName,0,true));
						}
						
						else if(mergeContigs && splittedId.length >= 4 && !splittedId[3].contains("_") && splittedId[3].charAt(12) == '.') {
							microbe2mergedContigs.put(splittedId[3].substring(0,6),new Microbe(splittedId[3].substring(0,6) + "000000",trivialName,0,true));
						}
						
					}
					indexReader.close();
				}
			}
			
			HashMap<String,String> speciesList = new HashMap<String,String>();
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			if(speciesListFilePath != null) {
				BufferedReader br = new BufferedReader(new FileReader(new File(speciesListFilePath)));
				while(br.ready()) {
					splittedLine = tabPattern.split(br.readLine());
					speciesList.put(splittedLine[0],splittedLine[1]);
				}
				br.close();
			}
			
			//process sam file
			BufferedReader samReader = new BufferedReader(new FileReader(new File(samFilePath)));
			HashMap<String,PrintWriter> species2writer = new HashMap<String,PrintWriter>();
			String currentFileName;
			Microbe microbe;
			while(samReader.ready()) {
				currentLine = samReader.readLine();
				st = new StringTokenizer(currentLine,"\t");
				//read id bitflag
				st.nextToken();
				st.nextToken();
				currentIndexName = st.nextToken();
				if(idx2microbes.containsKey(currentIndexName)) {
					
					microbe2start = idx2microbes.get(currentIndexName);
					microbeStart = Integer.valueOf(st.nextToken());
					microbe = microbe2start.get(microbe2start.floorKey(microbeStart));
					microbeId = microbe.getId();
					trivialName = microbe.getTrivialName();
					if(mergeContigs && microbeId.contains("NZ_")) {
						if(speciesList.isEmpty() || speciesList.containsKey(microbeId.substring(0,9) + "000000")) {
							if(!species2writer.containsKey(microbeId.substring(0,9))) {
								if(speciesList.isEmpty() || !trivialFileNames)
									currentFileName = microbeId.substring(0,9) + "000000.sam";
								else
									currentFileName = speciesList.get(microbeId.substring(0,9) + "000000").replace(" ", "_").replace("/", ".") + ".sam";
									
								species2writer.put(microbeId.substring(0,9),new PrintWriter(new FileWriter(new File(outputDirPath + "/" + currentFileName))));
							}
							
							species2writer.get(microbeId.substring(0,9)).println(currentLine.replace("NM:i:", ""));
						}
					}
					else if(mergeContigs && microbe2mergedContigs.containsKey(microbeId.substring(0,6))) {
						if(speciesList.isEmpty() || speciesList.containsKey(microbeId.substring(0,6) + "000000")) {
							if(!species2writer.containsKey(microbeId.substring(0,6))) {
								if(speciesList.isEmpty() || !trivialFileNames)
									currentFileName = microbeId.substring(0,6) + "000000.sam";
								else
									currentFileName = speciesList.get(microbeId.substring(0,6) + "000000").replace(" ", "_").replace("/", ".") + ".sam";
									
								species2writer.put(microbeId.substring(0,6),new PrintWriter(new FileWriter(new File(outputDirPath + "/" + currentFileName))));
							}
							
							species2writer.get(microbeId.substring(0,6)).println(currentLine.replace("NM:i:", ""));
						}
					}
					
					else {
						if(speciesList.isEmpty() || speciesList.containsKey(microbeId)) {
							if(!species2writer.containsKey(microbeId)) {
								
								if((speciesList.isEmpty() && trivialName == null) || !trivialFileNames)
									currentFileName = microbeId + ".sam";
								else {
									if(!speciesList.isEmpty())
										currentFileName = speciesList.get(microbeId).replace(" ", "_").replace("/", ".") + ".sam";
									else
										currentFileName = trivialName.replace(" ", "_").replace("/", ".") + ".sam";
								}
								
								species2writer.put(microbeId,new PrintWriter(new FileWriter(new File(outputDirPath + "/" + currentFileName))));
							}
							species2writer.get(microbeId).println(currentLine.replace("NM:i:", ""));
						}
					}					
				}
			}
			samReader.close();
			for(PrintWriter pw : species2writer.values())
				pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * checks if:	1) a line begins with a header symbol (@ as char 0). if so, this line will be skipped
	 * 				2) chr column contains a * symbol (unmapped read). if so, this line will be skipped
	 * 				3) the chromsome column contains the prefix "chr". if not, it will be added 
	 * 				4) matches in the cigar string add up to 'readlen'. if not, the last match will be extended
	 */
	
	public void fixSam(String samFilePath, int readlen , String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			while((currentLine = br.readLine()) != null) {
				
				//check (1)
				if(currentLine.charAt(0) == '@')
					continue;
				
				splittedLine = tabPattern.split(currentLine);
				
				//check (2)
				if(splittedLine[2].equals("*"))
					continue;
				
				//check (3)
				if(splittedLine[2].length() < 3 || !splittedLine[2].substring(0,3).equals("chr"))
					splittedLine[2] = "chr" + splittedLine[2];
				
				//check (4)
				splittedLine[5] = checkAndFixCigar(splittedLine[5],readlen);
				
				//print
				pw.print(splittedLine[0]);
				for(int i = 1; i < splittedLine.length; i++) {
					pw.print("\t" + splittedLine[i]);
				}
				pw.println();
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
		
	private String checkAndFixCigar(String cigString, int referenceLength) {
		int readLength = 0;
		int tmpNumber;
		int positionOfLastMatch = -1;
		int lastMatchCount = -1;
		String tmpNumberString = "";
		for(int i = 0; i < cigString.length(); i++) {
			try {
				tmpNumber = Integer.valueOf(cigString.substring(i,i+1));
				tmpNumberString += cigString.charAt(i);
			}
			catch(Exception e) {
				if(tmpNumberString.length() > 0) {
					tmpNumber = Integer.valueOf(tmpNumberString);
					if(cigString.charAt(i) == 'M') {
						readLength += tmpNumber;
						lastMatchCount = tmpNumber;
						positionOfLastMatch = i - String.valueOf(lastMatchCount).length();
					}
					tmpNumberString = "";
				}
			}
		}
		if(readLength == referenceLength)
			return cigString;
		
		return (cigString.substring(0,positionOfLastMatch) + String.valueOf((lastMatchCount + (referenceLength - readLength))) + "M");
	}
	
	private class Microbe {
		
		private String id;
		private String trivialName;
		private int start;
		private int end;
		private int genomeLength;
		private int containedReads;
		private boolean mergedContigs;
		
		private HashSet<Integer> coveredPositions;
		
		public Microbe(String id,String trivialName, int genomeLength, boolean mergedContigs) {
			this.id = id;
			this.trivialName = trivialName;
			this.genomeLength = genomeLength;
			this.containedReads = 0;
			this.mergedContigs = mergedContigs;
			this.coveredPositions = new HashSet<Integer>();
		}
		
		public Microbe(String id,int start, int end, boolean mergedContigs) {
			this.id = id;
			this.start = start;
			this.end = end;
			this.genomeLength = this.end - this.start;
			this.containedReads = 0;
			this.mergedContigs = mergedContigs;
			this.coveredPositions = new HashSet<Integer>();
		}
		
		public void addCoveredPosition(int pos) {
			this.coveredPositions.add(pos);
		}
		
		public HashSet<Integer> getCoveredPositions() {
			return this.coveredPositions;
		}
		
		public String getId() {
			return this.id;
		}
		
		public String getTrivialName() {
			return this.trivialName;
		}
		
		public void setTrivialName(String trivialName) {
			this.trivialName = trivialName;
		}
		
		public int getLength() {
			return this.genomeLength;
		}
		
		public int getContainedReads() {
			return this.containedReads;
		}
		
		public void setContainedReads(int reads) {
			this.containedReads = reads;
		}
		
		public void incrementContainedReads() {
			this.containedReads++;
		}
		
		public boolean hasMergedContigs() {
			return this.mergedContigs;
		}
	}
	
	
	public void getCoverage(String samFilePath, String outputFilePath, String indexFolderPath, int readLength) {
		try {
			TreeMap<Integer,Integer> coverage = new TreeMap<Integer,Integer>();
			TreeMap<Integer,Integer> mismatchCounts = new TreeMap<Integer,Integer>();
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			int currentStart;
			int currentEnd;
			int relativeGenomeStart = getRelativeGenomeStart(samFilePath,indexFolderPath);
			int mismatches;
			String meanMismatches;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				currentStart = Integer.valueOf(splittedLine[3]) - relativeGenomeStart;
				currentEnd = currentStart + readLength - 1;
				mismatches = Integer.valueOf(splittedLine[11]);
				for(int i =  currentStart; i <= currentEnd; i++) {
					
					//update coverages
					if(coverage.containsKey(i))
						coverage.put(i, coverage.get(i) + 1);
					else
						coverage.put(i, 1);
					
					//update mismatch counts
					if(mismatchCounts.containsKey(i))
						mismatchCounts.put(i, mismatchCounts.get(i) + mismatches);
					else
						mismatchCounts.put(i, mismatches);
				}
			}
			br.close();
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			pw.println("position\tcoverage\tmean_mismatches");
			for(int i : coverage.keySet()) {
				meanMismatches = trimDouble((double)mismatchCounts.get(i)/(double)coverage.get(i));
				pw.println(String.format("%s\t%s\t%s",i,coverage.get(i),meanMismatches));
			}
			pw.close();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private int getRelativeGenomeStart(String samFilePath, String indexFolderPath) {
		try {
			
			
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			String firstLine = br.readLine();
			String[] splittedLine = firstLine.split("\t");
			String chrName = splittedLine[2];
			int start = Integer.valueOf(splittedLine[3]);
			br.close();
			
			//get indexFile
			File[] indexFiles = new File(indexFolderPath).listFiles();
			File indexFile = null;
			for(File tmpFile : indexFiles) {
				if(tmpFile.getName().equals(chrName + ".idx")) {
					indexFile = tmpFile;
					break;
				}
			}
			
			if(indexFile == null)
				return Integer.MIN_VALUE;
			
			br = new BufferedReader(new FileReader(indexFile));
			String currentLine;
			Pattern tabPattern = Pattern.compile("\t");
			int currentStart;
			int currentEnd;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				currentStart = Integer.valueOf(splittedLine[1]);
				currentEnd = Integer.valueOf(splittedLine[2]);
				if(start >= currentStart && start <= currentEnd)
					return currentStart;
			}
			br.close();
			return Integer.MIN_VALUE;
		}
		catch(Exception e) {
			e.printStackTrace();
			return Integer.MIN_VALUE;
		}
	}
	
	private String trimDouble(double inValue){
		DecimalFormat twoDec = new DecimalFormat("0.00",new DecimalFormatSymbols(Locale.US));
		twoDec.setGroupingUsed(false);
		return twoDec.format(inValue);
		}
	
	
	public void extractSequences(String fastaFilePath, String samFilePath) {
		try {
			BufferedReader samReader = new BufferedReader(new FileReader(new File(samFilePath)));
			HashSet<String> readIds = new HashSet<String>();
			while(samReader.ready()) {
				readIds.add(samReader.readLine().split("\t")[0]);
			}
			samReader.close();
			
			BufferedReader fastaReader = new BufferedReader(new FileReader(new File(fastaFilePath)));
			String currentLine;
			while(fastaReader.ready()) {
				currentLine = fastaReader.readLine();
				if(currentLine.charAt(0) == '>' && readIds.contains(currentLine.substring(1))) {
					System.out.println(currentLine);
					System.out.println(fastaReader.readLine());
				}
			}
			fastaReader.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	public void modifyPairedEndStrandInformation(String samFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			String currentLine;
			Pattern spacePattern = Pattern.compile(" ");
			Pattern doublePointPattern = Pattern.compile(":");
			Pattern tabPattern = Pattern.compile("\t");
			int mateNumber;
			int flag;
			String[] splittedLine;
			while((currentLine = br.readLine()) != null) {
				mateNumber = Integer.valueOf(doublePointPattern.split(spacePattern.split(currentLine)[1])[0]);
				if(mateNumber == 2) {
					splittedLine = tabPattern.split(currentLine);
					flag = Integer.valueOf(splittedLine[1]);
					if(flag == 256)
						splittedLine[1] = String.valueOf(272);
					else
						splittedLine[1] = String.valueOf(256);
					
					pw.print(splittedLine[0]);
					for(int i = 1; i < splittedLine.length; i++) {
						pw.print("\t" + splittedLine[i]);
					}
					pw.print("\n");
				}
				else
					pw.println(currentLine);
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}
