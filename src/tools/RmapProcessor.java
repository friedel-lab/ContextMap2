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
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

import main.UnsynchronizedBufferedWriter;

public class RmapProcessor {

	public RmapProcessor() {
		
	}
	
	/**
	 * expects an rmap directory which contains rmap files splitted by chromosome and sorted by start positions.
	 * 
	 * @param inputFilePath
	 * @param outputFilePath
	 * @param annotationFilePath
	 */
	
	
	public void removeAlignmentsToPseudogenes(String inputDirPath, String outputFilePath, String annotationFilePath) {
		try {
			
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * expects a rmap file of paired end reads with headers of the form xxx#1 and xxx#2
	 * outputs read headers of the form xxx/1 and xxx/2 
	 * 
	 * @param inputFilePath
	 * @param outputFilePath
	 */
	
	
	private static class HeaderRenamer extends Thread {
		private String inputFilePath;
		private String outputFilePath;
		
		private long startPointer;
		private long stopPointer;
		
		private boolean pairedEnd;
		private boolean hasPairedEndHeader;
		
		public HeaderRenamer(String inputFilePath, long startPointer, long stopPointer, String outputFilePath, boolean pairedEnd, boolean hasPairedEndHeader) {
			this.inputFilePath = inputFilePath;
			this.startPointer = startPointer;
			this.stopPointer = stopPointer;
			this.outputFilePath = outputFilePath;
			this.pairedEnd = pairedEnd;
			this.hasPairedEndHeader = hasPairedEndHeader;
		}
		
		public void run() {
			try {
				BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(inputFilePath), "r", 10000);
				braf.seek(this.startPointer);
				PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
				
				String currentLine;
				String[] splittedLine;
				Pattern tabPattern = Pattern.compile("\t");
				StringBuffer sb = new StringBuffer();
				String readId;
				int replaceIndex;
				while((currentLine = braf.getNextLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					readId = splittedLine[0];
					readId = readId.replace("#*#"," ");
					
					replaceIndex = readId.length() - 2;
					
					sb.setLength(0);
					sb.append(readId);
					//if(this.hasPairedEndHeader && !readId.contains("::MSC::"))
					if((!readId.contains("::MSC::") && (this.pairedEnd || this.hasPairedEndHeader)) || (readId.contains("::MSC::") && this.pairedEnd))
						sb.setCharAt(replaceIndex, '/');
					
					for(int i = 1; i < splittedLine.length; i++) {
						sb.append("\t" + splittedLine[i]);
					}
					
					pw.println(sb.toString());
					
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
	
	public static void undoHeaderRenamingFromBwa(String inputFilePath, String outputDirPath, String outputFilePath, int threads, boolean pairedEnd, boolean hasPairedEndHeader) {
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
			
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				braf.seek(currentPosition);
				braf.getNextLine();
				currentPosition = braf.getFilePointer();
				
				
				tmpOutputFilePath = outputDirPath + "/renamed_ids_" + chunkIndex + ".rmap";
				futures.add(executor.submit(new HeaderRenamer(inputFilePath, prevPosition, currentPosition, tmpOutputFilePath,pairedEnd,hasPairedEndHeader)));
				
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
	
	
	public static void undoHeaderRenamingFromBwa(String inputFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			StringBuffer sb = new StringBuffer();
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				sb.setLength(0);
				sb.append(splittedLine[0]);
				sb.setCharAt(sb.length() - 2, '/');
				for(int i = 1; i < splittedLine.length; i++) {
					sb.append("\t" + splittedLine[i]);
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
	 * splits for a read of type 'S' (split read) the start and end coordinates and generates therefore
	 * an additional column.
	 * Further the SAM format is 1 based, so we also generate a 1 based rmap file here.
	 * The input file will be deleted
	 * @param inputPath
	 * @param outputPath
	 */
	public void modifyRmapAnnotation(String inputPath, String outputPath, boolean skipPartialHits) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)),1024 * 1024);
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(outputPath)),1024 * 1024);
			String currentLine;
			StringTokenizer st;
			String readId;
			String mappingType;
			String chr;
			String mappingPosition;
			String[] splittedMappingPosition;
			int start;
			int end;
			Pattern commaPattern = Pattern.compile(",");
			StringBuilder currentLineBuilder = new StringBuilder();
			while(br.ready()) {
				currentLine = br.readLine();
				st = new StringTokenizer(currentLine,"\t");
				readId = st.nextToken();
				mappingType = st.nextToken();
				if(skipPartialHits && mappingType.equals("P"))
					continue;
				chr = st.nextToken();
				//read_id	mapping_type	chr
				currentLineBuilder.setLength(0);
				currentLineBuilder.append(readId).append("\t").append(mappingType).append("\t").append(chr).append("\t");
				pw.write(currentLineBuilder.toString());
				//pw.write(String.format("%s\t%s\t%s\t",readId,mappingType,chr));
				mappingPosition = st.nextToken();
				splittedMappingPosition = commaPattern.split(mappingPosition);
				start = Integer.valueOf(splittedMappingPosition[0]) + 1;
				pw.write(String.valueOf(start));
				if(splittedMappingPosition.length == 1)
					pw.write("\t.");
				else
					pw.write("\t" + splittedMappingPosition[1]);
				
				
				currentLineBuilder.setLength(0);
				currentLineBuilder.append("\t").append(st.nextToken()).append("\t").append(st.nextToken());
				pw.write(currentLineBuilder.toString());
				//pw.write(String.format("\t%s\t%s",st.nextToken(),st.nextToken()));
				pw.newLine();
			}
			br.close();
			pw.close();
			new File(inputPath).delete();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void modifyBowtieOutput(String inputPath, String outputPath, int maxMismatches) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)),1024 * 1024);
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(outputPath)),1024 * 1024);
			String currentLine;
			StringTokenizer st;
			String readId;
			String strand;
			String chr;
			int start;
			int end;
			Pattern commaPattern = Pattern.compile(",");
			int mismatches;
			StringBuilder currentLineBuilder = new StringBuilder();
			while(br.ready()) {
				currentLine = br.readLine();
				st = new StringTokenizer(currentLine,"\t");
				readId = st.nextToken();
				strand = st.nextToken();
				chr = st.nextToken();
				start = Integer.valueOf(st.nextToken()) + 1;
				if(st.hasMoreTokens())
					mismatches = commaPattern.split(st.nextToken()).length;
				else
					mismatches = 0;
				
				if(mismatches <= maxMismatches) {
					currentLineBuilder.setLength(0);
					currentLineBuilder.append(readId).append("\tF\t").append(chr).append("\t").append(start).append("\t.\t").append(strand).append("\t").append(mismatches);
					pw.write(currentLineBuilder.toString());
					pw.newLine();
				}
			}
			br.close();
			pw.close();
			new File(inputPath).delete();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public void modifyRmapAnnotation(String inputPath, String outputPath, boolean skipPartialHits, int numberOfThreads) {
		try {
			
			if(numberOfThreads == 1) {
				modifyRmapAnnotation(inputPath, outputPath, skipPartialHits);
				return;
			}
			
			//determine offset
			BufferedRandomAccessFile rf = new BufferedRandomAccessFile(new File(inputPath),"r",1024);
			long fileSize = rf.getChannel().size();
			
			//in case we have an empty input file we just create an empty output file
			if(fileSize == 0) {
				new File(outputPath).createNewFile();
				new File(inputPath).delete();
				return;
			}
			
			int maxThreads = numberOfThreads;
			if(numberOfThreads > 10) {
				maxThreads = 10;
			}
			long incRate = fileSize/maxThreads;
			
			
			long prevPosition = 0;
			long currentPosition = 0;
			ExecutorService executor = Executors.newFixedThreadPool(maxThreads);
			ArrayList<Future> futures = new ArrayList<Future>();
			int chunkIndex = 0;
			String currentOutputPath;
			ArrayList<String> rmapFilePaths = new ArrayList<String>();
			
			
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				rf.seek(currentPosition);
				rf.readLine();
				currentPosition = rf.getFilePointer();
			
				currentOutputPath = outputPath + chunkIndex;
				futures.add(executor.submit(new RmapModifier(inputPath,currentOutputPath,skipPartialHits,prevPosition,currentPosition)));
				
				rmapFilePaths.add(currentOutputPath);
				prevPosition = currentPosition;
				chunkIndex++;
			}
			rf.close();
			
			
			executor.shutdown();
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			
			//delete unmodified file
			new File(inputPath).delete();
			
			if(!rmapFilePaths.isEmpty()) 
				concatenateFilesWithNIO(rmapFilePaths, outputPath);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private class RmapModifier extends Thread {
		
		private String inputPath;
		private String outputPath;
		boolean skipPartialHits;
		long start;
		long stop;
		
		public RmapModifier(String inputPath, String outputPath, boolean skipPartialHits, long start,long stop) {
			this.inputPath = inputPath;
			this.outputPath = outputPath;
			this.skipPartialHits = skipPartialHits;
			this.start = start;
			this.stop = stop;
		}
		
		public void run() {
			try {
				BufferedRandomAccessFile brf = new BufferedRandomAccessFile(new File(this.inputPath),"r",1024 * 1024);
				brf.seek(this.start);
				
				UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(outputPath)),1024 * 1024);
				String currentLine;
				StringTokenizer st;
				String readId;
				String mappingType;
				String chr;
				String mappingPosition;
				String[] splittedMappingPosition;
				int start;
				int end;
				Pattern commaPattern = Pattern.compile(",");
				StringBuilder currentLineBuilder = new StringBuilder();
				while((currentLine = brf.getNextLine()) != null) {
					st = new StringTokenizer(currentLine,"\t");
					readId = st.nextToken();
					mappingType = st.nextToken();
					if(skipPartialHits && mappingType.equals("P"))
						continue;
					chr = st.nextToken();
					//read_id	mapping_type	chr
					currentLineBuilder.setLength(0);
					currentLineBuilder.append(readId).append("\t").append(mappingType).append("\t").append(chr).append("\t");
					pw.write(currentLineBuilder.toString());
					//pw.write(String.format("%s\t%s\t%s\t",readId,mappingType,chr));
					mappingPosition = st.nextToken();
					splittedMappingPosition = commaPattern.split(mappingPosition);
					start = Integer.valueOf(splittedMappingPosition[0]) + 1;
					pw.write(String.valueOf(start));
					if(splittedMappingPosition.length == 1)
						pw.write("\t.");
					else
						pw.write("\t" + splittedMappingPosition[1]);
					
					
					currentLineBuilder.setLength(0);
					currentLineBuilder.append("\t").append(st.nextToken()).append("\t").append(st.nextToken());
					pw.write(currentLineBuilder.toString());
					//pw.write(String.format("\t%s\t%s",st.nextToken(),st.nextToken()));
					pw.newLine();
					
					if(brf.getFilePointer() == this.stop)
						break;
				}
				brf.close();
				pw.flush();
				pw.close();
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	
	
	/**
	 * simply splits a rmap file by chromosome id and finally deletes the input file.
	 * @param inputPath
	 * @param outputDir
	 */
	public void splitByChromosome(String inputPath, String outputDir) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)),1024 * 1024);
			String currentLine;
			StringTokenizer st;
			String chr;
			HashMap<String,UnsynchronizedBufferedWriter> writers = new HashMap<String,UnsynchronizedBufferedWriter>();
			UnsynchronizedBufferedWriter currentWriter;
			while(br.ready()) {
				currentLine = br.readLine();
				st = new StringTokenizer(currentLine,"\t");
				//skip read_id (context_id) and mapping_type (read_id)
				st.nextToken();
				st.nextToken();
				chr = st.nextToken();
				if(!writers.containsKey(chr))
					writers.put(chr,new UnsynchronizedBufferedWriter(new FileWriter(new File(String.format("%s/%s.rmap",outputDir,chr))),1024 * 1024));
				
				currentWriter = writers.get(chr);
				currentWriter.write(currentLine);
				currentWriter.newLine();
			}

			br.close();
			for(UnsynchronizedBufferedWriter pw : writers.values())
				pw.close();
			
			new File(inputPath).delete();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void splitByChromosome(String inputPath, String outputDir, int numberOfThreads) {
		try {
			if(numberOfThreads == 1) {
				splitByChromosome(inputPath,outputDir);
				return;
			}
			
			//determine offset
			BufferedRandomAccessFile rf = new BufferedRandomAccessFile(new File(inputPath),"r",1024);
			long fileSize = rf.getChannel().size();
			
			//here we restrict the maximum pool size, because every thread opens for each chr a file handle.
			int maxThreads = numberOfThreads;
			if(numberOfThreads > 6) {
				maxThreads = 6;
			}
			long incRate = fileSize/maxThreads;
			
			
			long prevPosition = 0;
			long currentPosition = 0;
			ExecutorService executor = Executors.newFixedThreadPool(maxThreads);
			ArrayList<Future> futures = new ArrayList<Future>();
			int chunkIndex = 0;
			String currentOutputDir;
			ArrayList<String> rmapFolderPaths = new ArrayList<String>();
			
			
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				rf.seek(currentPosition);
				rf.readLine();
				currentPosition = rf.getFilePointer();
			
				currentOutputDir = outputDir + chunkIndex;
				futures.add(executor.submit(new ChromosomeSplitter(inputPath, currentOutputDir,prevPosition,currentPosition)));
				
				rmapFolderPaths.add(currentOutputDir);
				prevPosition = currentPosition;
				chunkIndex++;
			}
			rf.close();
			
			
			executor.shutdown();
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			
			//delete unsplitted file
			new File(inputPath).delete();
			
			
			//finally we concatenate each chromosome file to a single file
			ArrayList<TreeMap<String,String>> chrFilePaths = new ArrayList<TreeMap<String,String>>();
			HashSet<String> availableChrNames = new HashSet<String>();
			for(int i = 0; i < rmapFolderPaths.size(); i++) {
				File[] files = new File(rmapFolderPaths.get(i)).listFiles();
				TreeMap<String,String> currentMap = new TreeMap<String,String>();
				for(int j = 0; j < files.length; j++) {
					currentMap.put(files[j].getName(), files[j].getAbsolutePath());
					availableChrNames.add(files[j].getName());
				}
				chrFilePaths.add(currentMap);
			}
			
			maxThreads = numberOfThreads;
			if(numberOfThreads > 10) {
				maxThreads = 10;
			}
			executor = Executors.newFixedThreadPool(maxThreads);
			futures.clear();
			ArrayList<String> currentChromosomeFilePaths ;
			for(String chr : availableChrNames) {
				currentChromosomeFilePaths = new ArrayList<String>();
				for(int i = 0; i < chrFilePaths.size(); i++) {
					if(chrFilePaths.get(i).containsKey(chr))
						currentChromosomeFilePaths.add(chrFilePaths.get(i).get(chr));
				}
				
				futures.add(executor.submit(new FileConcatenator(currentChromosomeFilePaths, String.format("%s/%s",outputDir,chr))));
			}
			
			executor.shutdown();
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			
			File[] tmpFiles;
			for(String rmapFolderPath : rmapFolderPaths) {
				tmpFiles = new File(rmapFolderPath).listFiles();
				for(File tmpFile : tmpFiles) {
					tmpFile.delete();
				}
				new File(rmapFolderPath).delete();
			}
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private class ChromosomeSplitter extends Thread {
		
		private String inputPath;
		private String outputDir;
		private long start;
		private long stop;
		
		public ChromosomeSplitter(String inputPath, String outputDir, long start, long stop) {
			this.inputPath = inputPath;
			this.outputDir = outputDir;
			this.start = start;
			this.stop = stop;
		}
		
		public void run() {
			try {
				File checkFile = new File(this.outputDir);
				if(!checkFile.isDirectory())
					checkFile.mkdirs();
				else {
					File[] oldFiles = new File(this.outputDir).listFiles();
					for(File oldFile : oldFiles) {
						oldFile.delete();
					}
				}
				
				BufferedRandomAccessFile reader = new BufferedRandomAccessFile(new File(this.inputPath),"r",1024 * 1024);
				reader.seek(this.start);
				
				String currentLine;
				StringTokenizer st;
				String chr;
				HashMap<String,UnsynchronizedBufferedWriter> writers = new HashMap<String,UnsynchronizedBufferedWriter>();
				UnsynchronizedBufferedWriter currentWriter;
				while((currentLine = reader.getNextLine()) != null) {
					st = new StringTokenizer(currentLine,"\t");
					//skip read_id (context_id) and mapping_type (read_id)
					st.nextToken();
					st.nextToken();
					chr = st.nextToken();
					if(!writers.containsKey(chr))
						writers.put(chr,new UnsynchronizedBufferedWriter(new FileWriter(new File(String.format("%s/%s.rmap",outputDir,chr))),1024 * 1024));
					
					currentWriter = writers.get(chr);
					currentWriter.write(currentLine);
					currentWriter.newLine();
					
					if(reader.getFilePointer() == this.stop)
						break;
				}

				reader.close();
				for(UnsynchronizedBufferedWriter pw : writers.values()) {
					pw.flush();
					pw.close();
				}
			}
			
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	
	
	
	/**
	 * concatenates inputed rmap files and subsequently deletes them.
	 * 
	 */	
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
	
	
	private class FileConcatenator extends Thread {
		
		private ArrayList<String> filePaths;
		private String outputPath;
		
		public FileConcatenator(ArrayList<String> filePaths, String outputPath) {
			this.filePaths = filePaths;
			this.outputPath = outputPath;
		}
		
		public void run() {
			try {
				File checkFile = new File(this.outputPath);
				if(checkFile.exists())
					checkFile.delete();
				
				//move largest file to output path
				int indexOfLargestFile = getIndexOfLargestFile(this.filePaths);
				new File(filePaths.get(indexOfLargestFile)).renameTo(new File(this.outputPath));
				this.filePaths.remove(indexOfLargestFile);
				
				
				FileOutputStream fos = new FileOutputStream(outputPath,true);
				FileChannel writeChannel = fos.getChannel();
				RandomAccessFile rf;
				FileChannel readChannel;
				long currentChannelSize;
				long transferedBytes;
				for(String filePath : this.filePaths) {
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
	

	
	
	
}
