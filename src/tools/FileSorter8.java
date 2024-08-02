package tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileChannel.MapMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import main.UnsynchronizedBufferedReader;
import main.UnsynchronizedBufferedWriter;



public class FileSorter8 extends Thread {

	
	private String inputPath;
	private String outputPath;
	private String tmpFolderPath;
	
	private HashSet<Integer> columnsOfInterest;
	private int linesPerChunkFile;
	//size in MegaByte
	private int chunkSize;
	private final int maxChunks = 100;
	private String delimiter;
	private boolean numericalSort;
	
	public static void main(String args[]) {
		int[] columnsToSortFor = null;
		int chunkSize = 50;
		String delimiter = "\t";
		String inputPath = null;
		String outputPath = null;
		String tmpFolderPath = null;
		boolean numericalSort = false;
		
		if(args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
			System.out.println("Sorts a file by specific column(s).");
			System.out.println("Available parameters:");
			System.out.println("-c\tcomma separated list of columns to sort for (0-based)");
			System.out.println("-d\tdelimiter which is used in the input file <default: tab>");
			System.out.println("-i\tpath to input file");
			System.out.println("-o\tpath to output file");
			System.out.println("-tmp\tpath to tmp folder for chunk files <default: outputpath/tmp>");
			System.out.println("-chunksize\tsize of chunk files in megabyte <default: 50>");
			System.out.println("--numerical\tnumerical sort <default: dictionary sort>");
			System.exit(0);
		}
		
		
		for(int i = 0; i < args.length; i++) {
			if(args[i].contains("-c")) {
				String[] tmpColumns = args[++i].split(",");
				columnsToSortFor = new int[tmpColumns.length];
				for(int j = 0; j < tmpColumns.length; j++) {
					columnsToSortFor[j] = Integer.valueOf(tmpColumns[j]);
				}
				continue;
			}
			if(args[i].contains("-d")) {
				delimiter = args[++i];
				continue;
			}
			if(args[i].contains("-i")) {
				inputPath = args[++i];
				continue;
			}
			if(args[i].contains("-o")) {
				outputPath = args[++i];
				continue;
			}
			if(args[i].contains("-tmp")) {
				tmpFolderPath = args[++i];
				continue;
			}
			
			if(args[i].contains("-chunksize")) {
				chunkSize = Integer.valueOf(args[++i]);
				continue;
			}
			if(args[i].contains("--numerical")) {
				numericalSort = true;
				continue;
			}
		}
		FileSorter8 fileSorter = new FileSorter8(inputPath,outputPath,tmpFolderPath,columnsToSortFor,chunkSize,delimiter,numericalSort);
		fileSorter.sortFile(inputPath, outputPath, tmpFolderPath);
	}
	
	//constructor for sorting files
	public FileSorter8(String inputPath, String outputPath, String tmpFolderPath, int[] columnsToSortFor,int chunkSize, String delimiter,boolean numericalSort) {
		this.inputPath = inputPath;
		this.outputPath = outputPath;
		this.tmpFolderPath = tmpFolderPath;
		this.columnsOfInterest = new HashSet<Integer>();
		for(int i = 0; i < columnsToSortFor.length; i++ ) {
			this.columnsOfInterest.add(columnsToSortFor[i]);
		}
		
		//first get the size of a line in bytes
		long lineSize = getLineSize(this.inputPath);
		this.chunkSize = chunkSize;
		long chunkSizeInBytes = chunkSize * 1024 * 1024;
		if(lineSize > 0)
			this.linesPerChunkFile = (int)(chunkSizeInBytes/lineSize);
		else
			this.linesPerChunkFile = 0;
		
		this.delimiter = delimiter;
		this.numericalSort = numericalSort;
	}
	
	
	
	//constructor for merging bunch of files
	public FileSorter8(File folderToMerge, int[] columnsToSortFor,int chunkSize, String delimiter,boolean numericalSort) {
		this.columnsOfInterest = new HashSet<Integer>();
		for(int i = 0; i < columnsToSortFor.length; i++ ) {
			this.columnsOfInterest.add(columnsToSortFor[i]);
		}
		File[] files = folderToMerge.listFiles();
		long lineSize = getLineSize(files[0].getAbsolutePath());
		this.chunkSize = chunkSize;
		long chunkSizeInBytes = chunkSize * 1024 * 1024;
		if(lineSize > 0)
			this.linesPerChunkFile = (int)(chunkSizeInBytes/lineSize);
		else
			this.linesPerChunkFile = 0;
		this.delimiter = delimiter;
		this.numericalSort = numericalSort;
	}
	
	
	public void run() {
		try {
			sortFile(this.inputPath,this.outputPath,this.tmpFolderPath);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
		
	//sorts a single file
	public void sortFile(String inputPath, String outputPath, String tmpFolderPath) {
		try {
			BufferedRandomAccessFile br = new BufferedRandomAccessFile(new File(inputPath),"r",10000);
			FileChannel brChannel = br.getChannel();
			
			FileOutputStream fos;
			FileChannel writeChannel;
			
			String chunkFolderPath = tmpFolderPath;
			if(tmpFolderPath == null)
				chunkFolderPath = outputPath.substring(0,outputPath.lastIndexOf(System.getProperty("file.separator"))) + "/tmp";
			File chunkFolder = new File(chunkFolderPath);
			if(!chunkFolder.isDirectory())
				chunkFolder.mkdirs();
			
			ArrayList<SparseLine> currentSparseLines = new ArrayList<SparseLine>();
			ArrayList<long[]> lineInfo;
			Comparator lineComparator;
			
			if(!this.numericalSort) {
				lineComparator = new LineComparatorForStringKeys();
			}
			else
				lineComparator = new LineComparatorForIntegerKeys();
			
	
			int processedLines = 0;
			int chunkCounter = 0;
			String currentLine;
			StringTokenizer tokenizer;
			String sortKey;
			int columnCounter;
			int addedTokens;
			String currentToken;
			StringBuilder tmpKey = new StringBuilder();
			long filePointer = br.getFilePointer();
			
			ArrayList<SparseLine> linePool = new ArrayList<SparseLine>();
			SparseLine tmpLine;
			ExecutorService executor = Executors.newFixedThreadPool(1);
			ArrayList<Future> futures = new ArrayList<Future>();
			while((currentLine = br.getNextLine()) != null) {
				if(this.delimiter.equals("\\t"))
					tokenizer = new StringTokenizer(currentLine, "\t");
				else
					tokenizer = new StringTokenizer(currentLine, this.delimiter);
				
				columnCounter = 0;
				addedTokens = 0;
				tmpKey.setLength(0);
				
				while(tokenizer.hasMoreTokens()) {
					currentToken = tokenizer.nextToken();
					if(this.columnsOfInterest.contains(columnCounter)) {
						tmpKey.append(currentToken);
						addedTokens++;
					}
					columnCounter++;
					if(addedTokens == this.columnsOfInterest.size())
						break;
				}
				sortKey = tmpKey.toString();
				tokenizer = null;
				
				if(currentSparseLines.size() < linePool.size()) {
					tmpLine = linePool.get(currentSparseLines.size());
					tmpLine.update(sortKey, filePointer,currentLine.length() + 1);
				}
				else {
					tmpLine = new SparseLine(sortKey,filePointer,currentLine.length() + 1);
					linePool.add(tmpLine);
				}
				
				currentSparseLines.add(tmpLine);
				
				if(++processedLines == this.linesPerChunkFile) {
					
					if(futures.size() >= 2) {
						for(Future future : futures) {
							future.get();
						}
						futures.clear();
					}
					
					//sort and write to chunk
					Collections.sort(currentSparseLines,lineComparator);
					lineInfo = new ArrayList<long[]>();
					for(SparseLine sparseLine : currentSparseLines) {
						lineInfo.add(new long[]{sparseLine.getFilePointer(),sparseLine.getLineLength()});
						
					}
					fos = new FileOutputStream(String.format(chunkFolderPath + "/chunk_%s.txt",chunkCounter),false);
					
					futures.add(executor.submit(new ConcurrentLineWriter(fos,inputPath,lineInfo)));
					
					
					processedLines = 0;
					currentSparseLines.clear();
					chunkCounter++;
				}
				
				filePointer = br.getFilePointer();
			}
			
			executor.shutdown();
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			
			//writing last chunk file
			Collections.sort(currentSparseLines,lineComparator);
			fos = new FileOutputStream(String.format(chunkFolderPath + "/chunk_%s.txt",chunkCounter),false);
			writeChannel = fos.getChannel();
			for(SparseLine sparseLine : currentSparseLines) {
				brChannel.transferTo(sparseLine.getFilePointer(), sparseLine.getLineLength(), writeChannel);
			}
			brChannel.close();
			br.close();
			br = null;
			
			writeChannel.close();
			fos.close();
			fos = null;
			
			currentSparseLines.clear();
			currentSparseLines = null;
			linePool.clear();
			linePool = null;
			
			//merging chunk files
			if(chunkCounter > 0)
				mergeFiles(chunkFolder,outputPath,this.maxChunks);
			else {
				File outputFile = new File(outputPath);
				if(outputFile.exists()) outputFile.delete();
				File mvChunkFile = new File(chunkFolderPath + "/chunk_0.txt");
				mvChunkFile.renameTo(outputFile);
			}
			
			//deleting tmp folder
			deleteFolderWithContent(chunkFolder);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private long getLineSize(String inputPath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			if(!br.ready()) {
				br.close();
				return 0;
			}
			int lineCounter = 0;
			String currentLine;
			long lineLength = 0;
			while(br.ready() && lineCounter < 10) {
				currentLine = br.readLine();
				lineLength += currentLine.length();
				lineCounter++;
			}
			
			br.close();
			return (lineLength/lineCounter);
		}
		
		catch(Exception e) {
			e.printStackTrace();
			return 0;
		}
	}
	
	private class SparseLine {
		
		private String sortKey;
		private long filePointer;
		private long lineLength;
		
		public SparseLine(String sortKey, long filePointer,long lineLength) {
			this.sortKey = sortKey;
			this.filePointer = filePointer;
			this.lineLength = lineLength;
		}
		
		public void update(String sortKey, long filePointer,long lineLength) {
			this.sortKey = sortKey;
			this.filePointer = filePointer;
			this.lineLength = lineLength;
		}
		
		public String getSortKey() {
			return this.sortKey;
		}
		
		public long getFilePointer() {
			return this.filePointer;
		}
		
		public long getLineLength() {
			return this.lineLength;
		}
	}
	

	private class ConcurrentLineWriter extends Thread {
		
		private FileOutputStream fos;
		private RandomAccessFile rf;
		private ArrayList<long[]> lines;
		
		public ConcurrentLineWriter(FileOutputStream fos, String filePath, ArrayList<long[]> lines) {
			super();
			try {
				this.fos = fos;
				this.rf = new RandomAccessFile(new File(inputPath),"r");
				this.lines = lines;
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		
		public void run() {
			try {
				FileChannel writeChannel = this.fos.getChannel();
				FileChannel readerChannel = this.rf.getChannel();
				
				long currentStart = this.lines.get(0)[0];
				long currentLength = this.lines.get(0)[1];
				for(int i = 1; i < this.lines.size(); i++) {
					if(this.lines.get(i)[0] == (currentStart + this.lines.get(i)[1]))
						currentLength += this.lines.get(i)[1];
					else {
						readerChannel.transferTo(currentStart, currentLength, writeChannel);
						currentStart = this.lines.get(i)[0];
						currentLength = this.lines.get(i)[1];
					}
				}
				readerChannel.transferTo(currentStart,currentLength, writeChannel);
				
				
				/*for(long[] lineInfo : this.lines) {
					readerChannel.transferTo(lineInfo[0],lineInfo[1], writeChannel);
				}*/
				
				
				readerChannel.close();
				writeChannel.close();
				this.fos.close();
				this.rf.close();
				this.lines.clear();
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}

	public void mergeFiles(File folderToMerge, String outputPath, int maxChunks) {
		try {
			File outputFile = new File(outputPath);
			ArrayList<File> filesToMerge = new ArrayList<File>(Arrays.asList(folderToMerge.listFiles()));
			
			//in case there are too many chunk files we process them here again in chunks
			int fileCount = getFileCount(filesToMerge);
			int subsetIndex = 0;
			int deletionIndex = 0;
			String subsetOutputPath;
			while(fileCount > maxChunks) {
				subsetOutputPath = String.format("%s/merged_subset_%s.txt",folderToMerge.getAbsolutePath(),subsetIndex);
				mergeFileSubset(filesToMerge,subsetOutputPath,maxChunks);
				filesToMerge.add(new File(subsetOutputPath));
				fileCount++;
				deletionIndex = 0;
				for(int i = 0; i < maxChunks;i++) {
					if(filesToMerge.get(deletionIndex).isFile()) {
						filesToMerge.get(deletionIndex).delete();
						filesToMerge.remove(deletionIndex);
						fileCount--;
					}
					else
						deletionIndex++;
				}
				subsetIndex++;
			}
			
			ArrayList<BufferedRandomAccessFile> readers = new ArrayList<BufferedRandomAccessFile>();
			ArrayList<FileChannel> channels = new ArrayList<FileChannel>();
			ArrayList<Line> currentRows = new ArrayList<Line>();
			long filePointer;
			String currentLine;
			for(int i = 0; i < filesToMerge.size(); i++) {
				if(filesToMerge.get(i).isFile() && !filesToMerge.get(i).getName().equals(outputFile.getName())) {
					BufferedRandomAccessFile currentReader = new BufferedRandomAccessFile(filesToMerge.get(i),"r",10240);
					filePointer= currentReader.getFilePointer();
					FileChannel currentChannel = currentReader.getChannel();
					if((currentLine = currentReader.getNextLine()) != null) {
						readers.add(currentReader);
						channels.add(currentChannel);
						Line tmpLine = new Line(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
						tmpLine.setReaderIndex(readers.size() - 1);
						currentRows.add(tmpLine);
						
					}
					else {
						currentChannel.close();
						currentReader.close();
					}
				}
			}
			MinPriorityQueue queue = new MinPriorityQueue(currentRows);
			FileOutputStream fos = new FileOutputStream(outputFile);
			FileChannel writerChannel = fos.getChannel();
			int readerIndex;
			int prevReaderIndex;
			long currentStart;
			long currentLength;
			Line maxLine;
			int processedLines = 0;
			
			if((maxLine = (Line)queue.extractMinimum()) != null) {
				readerIndex = maxLine.getReaderIndex();
				prevReaderIndex = readerIndex;
				currentStart = maxLine.getFilePointer();
				currentLength = maxLine.getLineLength();
				filePointer = readers.get(readerIndex).getFilePointer();
				if((currentLine = readers.get(readerIndex).getNextLine()) != null) {
					maxLine.updateLine(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
					maxLine.setReaderIndex(readerIndex);
					queue.insert(maxLine);
				}
				
				
				while((maxLine = (Line)queue.extractMinimum()) != null) {
					readerIndex = maxLine.getReaderIndex();
					if(readerIndex == prevReaderIndex) {
						currentLength += maxLine.getLineLength();
					}
					else {
						channels.get(prevReaderIndex).transferTo(currentStart, currentLength, writerChannel);
						prevReaderIndex = readerIndex;
						currentStart = maxLine.getFilePointer();
						currentLength = maxLine.getLineLength();
					}
					
					filePointer = readers.get(readerIndex).getFilePointer();
					if((currentLine = readers.get(readerIndex).getNextLine()) != null) {
						maxLine.updateLine(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
						maxLine.setReaderIndex(readerIndex);
						queue.insert(maxLine);
					}
					processedLines++;
					
				}
				channels.get(prevReaderIndex).transferTo(currentStart, currentLength, writerChannel);
			}
			queue = null;
			
			for(int i = 0; i < readers.size(); i++) {
				channels.get(i).close();
				readers.get(i).close();
			}
			channels.clear();
			readers.clear();
			
			
			writerChannel.close();
			fos.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	
	private int getFileCount(ArrayList<File> files) {
		int fileCounter = 0;
		for(File f : files) {
			if(f.isFile())
				fileCounter++;
		}
		return fileCounter;
	}
	
	private void mergeFileSubset(ArrayList<File> filesToMerge, String outputPath,int maxChunks) {
		try {
			File outputFile = new File(outputPath);
			ArrayList<BufferedRandomAccessFile> readers = new ArrayList<BufferedRandomAccessFile>();
			ArrayList<FileChannel> channels = new ArrayList<FileChannel>();
			ArrayList<Line> currentRows = new ArrayList<Line>();
			long filePointer;
			String currentLine;
			for(int i = 0; i < filesToMerge.size() && i < maxChunks; i++) {
				if(filesToMerge.get(i).isFile() && !filesToMerge.get(i).getName().equals(outputFile.getName())) {
					BufferedRandomAccessFile currentReader = new BufferedRandomAccessFile(filesToMerge.get(i),"r",10240);
					filePointer = currentReader.getFilePointer();
					FileChannel currentChannel = currentReader.getChannel();
					
					if((currentLine = currentReader.getNextLine()) != null) {
						readers.add(currentReader);
						channels.add(currentChannel);
						Line tmpLine = new Line(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
						tmpLine.setReaderIndex(readers.size() - 1);
						currentRows.add(tmpLine);
						
					}
					else {
						currentChannel.close();
						currentReader.close();
						readers.remove(currentReader);
					}
				}
			}
			
			MinPriorityQueue queue = new MinPriorityQueue(currentRows);
			FileOutputStream fos = new FileOutputStream(outputFile);
			FileChannel writerChannel = fos.getChannel();
			int readerIndex;
			int prevReaderIndex;
			long currentStart;
			long currentLength;
			Line maxLine;
			int processedLines = 0;
			
			if((maxLine = (Line)queue.extractMinimum()) != null) {
				readerIndex = maxLine.getReaderIndex();
				prevReaderIndex = readerIndex;
				currentStart = maxLine.getFilePointer();
				currentLength = maxLine.getLineLength();
				filePointer = readers.get(readerIndex).getFilePointer();
				if((currentLine = readers.get(readerIndex).getNextLine()) != null) {
					maxLine.updateLine(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
					maxLine.setReaderIndex(readerIndex);
					queue.insert(maxLine);
				}
				
				while((maxLine = (Line)queue.extractMinimum()) != null) {
					readerIndex = maxLine.getReaderIndex();
					if(readerIndex == prevReaderIndex) {
						currentLength += maxLine.getLineLength();
					}
					else {
						channels.get(prevReaderIndex).transferTo(currentStart, currentLength, writerChannel);
						prevReaderIndex = readerIndex;
						currentStart = maxLine.getFilePointer();
						currentLength = maxLine.getLineLength();
					}
					
					filePointer = readers.get(readerIndex).getFilePointer();
					if((currentLine = readers.get(readerIndex).getNextLine()) != null) {
						maxLine.updateLine(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
						maxLine.setReaderIndex(readerIndex);
						queue.insert(maxLine);
					}
					processedLines++;
					
				}
				channels.get(prevReaderIndex).transferTo(currentStart, currentLength, writerChannel);
			}
			queue = null;
			
			for(int i = 0; i < readers.size(); i++) {
				channels.get(i).close();
				readers.get(i).close();
			}
			channels.clear();
			readers.clear();
			
			
			writerChannel.close();
			fos.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	

	private class Line implements Comparable  {
		private String sortKey;
		private int readerIndex;
		private long filePointer;
		private long lineLength;
		private boolean numericalSort;
		
		public Line(String completeLine, long filePointer, HashSet<Integer> columnsToSortFor, String splitPattern, boolean numericalSort) {
			this.filePointer = filePointer;
			this.lineLength = completeLine.length() + 1;
			StringTokenizer tokenizer;
			if(splitPattern.equals("\\t"))
				tokenizer = new StringTokenizer(completeLine, "\t");
			else
				tokenizer = new StringTokenizer(completeLine, splitPattern);
			
			this.numericalSort = numericalSort;
			this.sortKey = "";
			int columnCounter = 0;
			int addedTokens = 0;
			String currentToken;
			while(tokenizer.hasMoreTokens()) {
				currentToken = tokenizer.nextToken();
				if(columnsToSortFor.contains(columnCounter)) {
					sortKey += currentToken;
					addedTokens++;
				}
				columnCounter++;
				if(addedTokens == columnsToSortFor.size())
					break;
			}
			
		}
		
		public void updateLine(String completeLine, long filePointer, HashSet<Integer> columnsToSortFor, String splitPattern, boolean numericalSort) {
			this.filePointer = filePointer;
			this.lineLength = completeLine.length() + 1;
			StringTokenizer tokenizer;
			if(splitPattern.equals("\\t"))
				tokenizer = new StringTokenizer(completeLine, "\t");
			else
				tokenizer = new StringTokenizer(completeLine, splitPattern);
			
			this.numericalSort = numericalSort;
			this.sortKey = "";
			int columnCounter = 0;
			int addedTokens = 0;
			String currentToken;
			while(tokenizer.hasMoreTokens()) {
				currentToken = tokenizer.nextToken();
				if(columnsToSortFor.contains(columnCounter)) {
					sortKey += currentToken;
					addedTokens++;
				}
				columnCounter++;
				if(addedTokens == columnsToSortFor.size())
					break;
			}
			
		}
		
		public long getFilePointer() {
			return this.filePointer;
		}
		
		public long getLineLength() {
			return this.lineLength;
		}
		
		public String getSortKey() {
			return this.sortKey;
		}
		
		public void setReaderIndex(int index) {
			this.readerIndex = index;
		}
		
		public int getReaderIndex() {
			return this.readerIndex;
		}

		@Override
		public int compareTo(Object o) {
			if(this.numericalSort)
				return(Double.valueOf(this.sortKey).compareTo(Double.valueOf(((Line)o).getSortKey())));
			else 
				return(this.sortKey.compareTo(((Line)o).getSortKey()));
		}
	}
	
	private class LineComparatorForStringKeys implements Comparator<SparseLine> {
		public int compare(SparseLine l1, SparseLine l2) {
			return(l1.getSortKey().compareTo(l2.getSortKey()));
		}

	}
	
	private class LineComparatorForIntegerKeys implements Comparator<SparseLine> {
		public int compare(SparseLine l1, SparseLine l2) {
			return(Double.valueOf(l1.getSortKey()).compareTo(Double.valueOf(l2.getSortKey())));
		}

	}
	
	
	private void deleteFolderWithContent(File f) {
		File[] files = f.listFiles();
		for(File tmpFile : files) {
			if(tmpFile.isDirectory())
				deleteFolderWithContent(tmpFile);
			else
				tmpFile.delete();
			}
		f.delete();
	}
	
}
