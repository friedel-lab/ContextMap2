package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;


public class SequenceConverter {

	
	private static final int qualityValue = 40;
	
	
	public static void main(String[] args) {
		String inputPath = null;
		String inputPath2 = null;
		String outputPath = null;
		String inputType = "fastq";
		int trimSize = -1;
		if(args[0].equals("convertSingleEnd")) {
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("-inputformat")) {
					inputType = args[++i];
					continue;
				}
				
			}
			if(inputType.equals("fa") || inputType.equals("fasta")) 
				convertToFastq(inputPath,outputPath);
			
			else if(inputType.equals("fastq"))
				convertToFasta(inputPath,outputPath);
			
			else
				System.err.println("unknown input type: " + inputType);
		}
		else if(args[0].equals("convertPairedEnd")) {
			boolean removeSpaces = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i1")) {
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-i2")) {
					inputPath2 = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("--rmspaces")) {
					removeSpaces = true;
				}
				
			}
			convertToPairedEnd(inputPath,inputPath2,outputPath,removeSpaces);
		}
		
		else if(args[0].equals("trimThreePrime")) {
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("-trimSize")) {
					trimSize = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-format")) {
					inputType = args[++i];
					continue;
				}
			}
			trimThreePrimeEnd(inputPath,outputPath,trimSize,inputType);
		}
		
		else if(args[0].equals("trimFivePrime")) {
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("-trimSize")) {
					trimSize = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-format")) {
					inputType = args[++i];
					continue;
				}
			}
			trimFivePrimeEnd(inputPath,outputPath,trimSize,inputType);
		}

		else if(args[0].equals("extractSequencesFromSample")) {
			String sample = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("-sample")) {
					sample = args[++i];
					continue;
				}
			}
			System.out.println(String.format("input path: %s", inputPath));
			System.out.println(String.format("output path: %s", outputPath));
			System.out.println(String.format("sample: %s", sample));
			extractSequencesFromSample(inputPath,outputPath,sample);
		}
		
		else if(args[0].equals("extractSequencesFromLane")) {
			String lane = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("-lane")) {
					lane = args[++i];
					continue;
				}
			}
			extractSequencesFromLane(inputPath,outputPath,lane);
		}
		
		else if(args[0].equals("splitMultiFasta")) {
			String fastaFilePath = null;
			String outputDirPath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					fastaFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputDirPath = args[++i];
					continue;
				}
			}
			splitMultiFasta(fastaFilePath,outputDirPath);
		}
		else if(args[0].equals("filterMultiFasta")) {
			String fastaFilePath = null;
			String outputFilePath = null;
			String pattern = null;
			boolean positiveFilter = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					fastaFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-p")) {
					pattern = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("--positive")) {
					positiveFilter = true;
					continue;
				}
			}
			filterMultiFasta(fastaFilePath,pattern,positiveFilter,outputFilePath);
		}
		
		else if(args[0].equals("filterFastaByIds")) {
			String fastaFilePath = null;
			String outputFilePath = null;
			String idsFilePath = null;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					fastaFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-ids")) {
					idsFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
					continue;
				}
				
			}
			filterFastaByIds(fastaFilePath,idsFilePath, outputFilePath);
		}
	}
	
	/**
	 * expects a fasta file of paired end reads in the input of the form: xxx/1 and xxx/2
	 * outputs read headers of the form xxx#1 and xxx#2 
	 * 
	 * @param inputFilePath
	 * @param outputFilePath
	 */
	
	
	private static class HeaderRenamer extends Thread {
		private String inputFilePath;
		private String outputFilePath;
		
		private long startPointer;
		private long stopPointer;
		
		private boolean hasPairedEndHeader;
		
		public HeaderRenamer(String inputFilePath, long startPointer, long stopPointer, boolean hasPairedEndHeader, String outputFilePath) {
			this.inputFilePath = inputFilePath;
			this.startPointer = startPointer;
			this.stopPointer = stopPointer;
			this.hasPairedEndHeader = hasPairedEndHeader;
			this.outputFilePath = outputFilePath;
		}
		
		public void run() {
			try {
				BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(inputFilePath),"r",10000);
				braf.seek(this.startPointer);
				PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
				
				String currentLine;
				StringBuffer sb = new StringBuffer();
				while((currentLine = braf.getNextLine()) != null) {
					if(currentLine.charAt(0) == '>') {
						
						//bwa strips spaces in the header line...
						currentLine = currentLine.replace(" ", "#*#");
						
						sb.setLength(0);
						sb.append(currentLine);
						
						if(this.hasPairedEndHeader)
							sb.setCharAt(sb.length()-2, '#');
						
						pw.println(sb.toString());
						pw.println(braf.readLine());
					}
					
					if(braf.getFilePointer() == this.stopPointer)
						break;
				}
				braf.close();
				pw.close();
			}
			
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	
	public static void renameReadHeaderForBwa(String inputFilePath, String outputDirPath, String outputFilePath, boolean hasPairedEndHeader, int threads) {
		try {
			BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(inputFilePath),"r",1024);
			long fileSize = braf.getChannel().size();
			long incRate = fileSize/threads;
			long prevPosition = 0;
			long currentPosition = 0;
			String tmpOutputFilePath;
			
			int chunkIndex = 0;
			
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			String currentLine;
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				braf.seek(currentPosition);
				braf.getNextLine();
				currentLine = braf.getNextLine();
				if(currentLine != null && currentLine.charAt(0) == '>')
					braf.getNextLine();
				
				currentPosition = braf.getFilePointer();
				
				tmpOutputFilePath = outputDirPath + "/renamed_ids_" + chunkIndex + ".fa";
				futures.add(executor.submit(new HeaderRenamer(inputFilePath, prevPosition, currentPosition, hasPairedEndHeader, tmpOutputFilePath)));
				
				prevPosition = currentPosition;
				chunkIndex++;
			}
			braf.close();
			
			executor.shutdown();
			for(Future f : futures) {
				f.get();
			}
			
			concatenateFiles(outputDirPath,outputFilePath);
			File[] tmpFiles = new File(outputDirPath).listFiles();
			for(File f : tmpFiles)
				f.delete();
		
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void concatenateFiles(String inputDirPath, String outputFilePath) throws Exception {
		File outputFile = new File(outputFilePath);
		if(outputFile.isFile())
			outputFile.delete();
		
		FileOutputStream fos = new FileOutputStream(outputFilePath,true);
		FileChannel writeChannel = fos.getChannel();
		RandomAccessFile rf;
		FileChannel readChannel;
		File[] tmpFiles = new File(inputDirPath).listFiles();
		long currentChannelSize;
		long transferedBytes;
		for(File tmpFile : tmpFiles) {
			rf = new RandomAccessFile(tmpFile,"r");
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
	
	public static void renameReadHeaderForBwa(String inputFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			StringBuffer sb = new StringBuffer();
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					sb.setLength(0);
					sb.append(currentLine);
					sb.setCharAt(sb.length()-2, '#');
					pw.println(sb.toString());
					pw.println(br.readLine());
				}
			}
			br.close();
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * input: (fasta || fastq file)
	 * output: length of the first read sequence
	 */
	
	public static int getReadLength(String readFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(readFilePath)));
			//skip header
			br.readLine();
			return((br.readLine()).length());
		}
		catch(Exception e) {
			e.printStackTrace();
			return -1;
		}
	}
	
	public static String getReadFormat(String readFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(readFilePath)));
			String header = br.readLine(); 
			if(header.charAt(0) == '>')
				return "fasta";
			else if(header.charAt(0) == '@')
				return "fastq";
			else
				return null;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private static void splitMultiFasta(String fastaFilePath, String outputDirPath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			PrintWriter pw = null;
			String currentLine;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					if(pw != null)
						pw.close();
					
					pw = new PrintWriter(new FileWriter(new File(outputDirPath + "/" + currentLine.substring(1) + ".fa")));
					pw.println(currentLine);
				}
				else if(pw != null)
					pw.println(currentLine);
			}
			br.close();
			if(pw != null)
				pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void filterMultiFasta(String fastaFilePath,String pattern,  boolean positiveFilter, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			String currentLine;
			boolean writeSequence = false;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					writeSequence = false;
					if(positiveFilter && currentLine.contains(pattern)) {
						pw.println(currentLine);
						pw.println(br.readLine());
						writeSequence = true;
					}
					else if(!positiveFilter && !currentLine.contains(pattern)) {
						pw.println(currentLine);
						pw.println(br.readLine());
						writeSequence = true;
					}
				}
				else if(writeSequence)
					pw.println(currentLine);
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void filterFastaByIds(String fastaFilePath,String idsFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(idsFilePath)));
			HashSet<String> ids = new HashSet<String>();
			String currentLine;
			while((currentLine = br.readLine()) != null) {
				ids.add(currentLine);
			}
			br.close();
			
			br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			boolean writeSequence = false;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					writeSequence = false;
					if(ids.contains(currentLine.substring(1))) {
						pw.println(currentLine);
						pw.println(br.readLine());
						writeSequence = true;
					}
				}
				else if(writeSequence)
					pw.println(currentLine);
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	//assumption: every read block consists of exactly four lines.
	//if the assumption holds, a '@' symbol in the quality line will be ignored...
	
	public static void convertToFasta(String inputPath, String outputPath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			
			String currentLine;
			String header;
			String sequence;
			
			while(br.ready()) {
				currentLine = br.readLine();
				
				if(currentLine.isEmpty())
					continue;
				
				if(currentLine.charAt(0) == '@') {
					header = String.format(">%s",currentLine.substring(1));
					sequence = br.readLine();
					pw.println(header);
					pw.println(sequence);
					//skip the next two fastq lines
					br.readLine();
					br.readLine();
				}
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void convertToFastq(String inputPath,String outputPath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			
			String currentLine;
			String header;
			String sequence;
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.charAt(0) == '>') {
					header = currentLine.replace('>', '@');
					sequence = br.readLine();
					pw.println(header);
					pw.println(sequence);
					pw.println("+");
					for(int i = 0; i < sequence.length(); i++)
						pw.print((char)qualityValue);
					pw.print("\n");
					
				}
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void convertToPairedEnd(String inputPath1, String inputPath2, String outputPath, boolean removeSpaces) {
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			convertHeader(inputPath1,"1",pw,removeSpaces);
			convertHeader(inputPath2,"2",pw,removeSpaces);
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * expects to fasta/fastq files, both in the same ordering of mates and with the same number of mates
	 * outputs a single fasta file, containing one mate after the other
	 * 
	 * @param inputFilePath1
	 * @param inputFilePath2
	 * @param outputFilePath
	 */
	
	public static void convertPairedToSingleFile(String inputFilePath1, String inputFilePath2, String outputFilePath) {
		try {
			//check if we are working with fasta or fastq files
			boolean isFastq = false;
			BufferedReader br1 = new BufferedReader(new FileReader(new File(inputFilePath1)));
			if(br1.ready()) {
				isFastq = (br1.readLine().charAt(0) == '@');
			}
			br1.close();
			
			
			br1 = new BufferedReader(new FileReader(new File(inputFilePath1)));
			BufferedReader br2 = new BufferedReader(new FileReader(new File(inputFilePath2)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			while(br1.ready() && br2.ready()) {
				//write first mate
				pw.println(">" + br1.readLine().substring(1));
				pw.println(br1.readLine());
				if(isFastq) {
					br1.readLine();
					br1.readLine();
				}
				
				//write second mate
				pw.println(">" + br2.readLine().substring(1));
				pw.println(br2.readLine());
				if(isFastq) {
					br2.readLine();
					br2.readLine();
				}
			}
			br1.close();
			br2.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}

	private static void convertHeader(String inputPath, String index, PrintWriter pw, boolean removeSpaces) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
		String currentLine = br.readLine();
		Pattern spacePattern = Pattern.compile(" ");
		String[] splittedLine = spacePattern.split(currentLine);
		StringBuffer sb = new StringBuffer();
		while(br.ready()) {
			currentLine = br.readLine();
			if(currentLine.charAt(0) == '>') {
				pw.print(splittedLine[0] + "/" + index);
				if(!removeSpaces) {
					for(int i = 1; i < splittedLine.length; i++) {
						pw.print(" " + splittedLine[i]);
					}
				}
				pw.print("\n");
				pw.println(sb.toString());
				
				splittedLine = spacePattern.split(currentLine);
				sb.setLength(0);
			}
			else
				sb.append(currentLine);
		}
		br.close();
	}
	
	private static void trimFivePrimeEnd(String inputPath, String outputPath, int trimSize, String readFormat) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			String currentLine;
			String currentSequence;
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.charAt(0) == '>' || currentLine.charAt(0) == '@') {
					currentSequence = br.readLine();
					
					if(currentSequence.length() < trimSize) {
						if(readFormat.equals("fastq")) {
							br.readLine();
							br.readLine();
						}
						continue;
					}
					
					pw.println(currentLine);
					pw.println(currentSequence.substring(0,trimSize));
					
					if(readFormat.equals("fastq")) {
						pw.println(br.readLine());
						pw.println(br.readLine().substring(0,trimSize));
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
	
	
	private static void trimThreePrimeEnd(String inputPath, String outputPath, int trimSize, String readFormat) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			String currentLine;
			String currentSequence;
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.charAt(0) == '>' || currentLine.charAt(0) == '@') {
					currentSequence = br.readLine();
					if(currentSequence.length() < trimSize) {
						if(readFormat.equals("fastq")) {
							br.readLine();
							br.readLine();
						}
						continue;
					}
					
					pw.println(currentLine);
					pw.println(currentSequence.substring(trimSize));
					
					if(readFormat.equals("fastq")) {
						pw.println(br.readLine());
						pw.println(br.readLine().substring(trimSize));
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
	
	public static void extractSequencesFromSample(String inputPath, String outputPath, String sample) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			String currentLine;
			String[] splittedLine;
			Pattern spacePattern = Pattern.compile(" ");
			
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.charAt(0) == '>') {
					splittedLine = spacePattern.split(currentLine);
					if(splittedLine[1].equals(sample)) {
						pw.println(currentLine);
						pw.println(br.readLine());
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
	
	public static void extractSequencesFromLane(String inputPath, String outputPath, String lane) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			String currentLine;
			String[] splittedLine;
			Pattern spacePattern = Pattern.compile(" ");
			
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.charAt(0) == '>') {
					splittedLine = spacePattern.split(currentLine);
					if(splittedLine[4].equals(String.format("lane=%s",lane))) {
						pw.println(currentLine);
						pw.println(br.readLine());
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
}
