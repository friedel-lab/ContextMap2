package tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
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



public class FileSorter7 extends Thread {

	
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
		FileSorter7 fileSorter = new FileSorter7(inputPath,outputPath,tmpFolderPath,columnsToSortFor,chunkSize,delimiter,numericalSort);
		fileSorter.sortFile(inputPath, outputPath, tmpFolderPath);
	}
	
	//constructor for sorting files
	public FileSorter7(String inputPath, String outputPath, String tmpFolderPath, int[] columnsToSortFor,int chunkSize, String delimiter,boolean numericalSort) {
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
	public FileSorter7(File folderToMerge, int[] columnsToSortFor,int chunkSize, String delimiter,boolean numericalSort) {
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
			String chunkFolderPath = tmpFolderPath;
			if(tmpFolderPath == null)
				chunkFolderPath = outputPath.substring(0,outputPath.lastIndexOf(System.getProperty("file.separator"))) + "/tmp";
			File chunkFolder = new File(chunkFolderPath);
			if(!chunkFolder.isDirectory())
				chunkFolder.mkdirs();
			
			
			ArrayList<SparseLine> currentSparseLines = new ArrayList<SparseLine>();
			Comparator lineComparator;
			
			if(!this.numericalSort) {
				lineComparator = new LineComparatorForStringKeys();
			}
			else
				lineComparator = new LineComparatorForIntegerKeys();
			
	
			int processedLines = 0;
			int chunkCounter = 0;
			String currentLine;
			UnsynchronizedBufferedWriter chunkPw = null;
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
			ArrayList<UnsynchronizedBufferedWriter> chunkPws = new ArrayList<UnsynchronizedBufferedWriter>();
			while((currentLine = br.getNextLine()) != null) {
				if(this.delimiter.equals("\\t"))
					tokenizer = new StringTokenizer(currentLine, "\t");
				else
					tokenizer = new StringTokenizer(currentLine, this.delimiter);
				
				currentLine = null;
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
					tmpLine.update(sortKey, filePointer);
				}
				else {
					tmpLine = new SparseLine(sortKey,filePointer);
					linePool.add(tmpLine);
				}
				currentSparseLines.add(tmpLine);
				filePointer = br.getFilePointer();
				if(++processedLines == this.linesPerChunkFile) {
					//wait until last chunk is written to hdd
					if(futures.size() >= 1) {
						for(int i = 0; i < futures.size(); i++) {
							futures.get(i).get();
							chunkPws.get(i).flush();
							chunkPws.get(i).close();
						}
						futures.clear();
						chunkPws.clear();
					}
					
					//sort and write to chunk
					Collections.sort(currentSparseLines,lineComparator);
					chunkPw = new UnsynchronizedBufferedWriter(new FileWriter(new File(String.format(chunkFolderPath + "/chunk_%s.txt",chunkCounter))),(int)(this.chunkSize * 1.05) * 1024 * 1024);
					chunkPws.add(chunkPw);
					for(SparseLine sparseLine : currentSparseLines) {
						br.seek(sparseLine.getFilePointer());
						chunkPw.write(br.getNextLine());
						chunkPw.newLine();
					}
					br.seek(filePointer);
					futures.add(executor.submit(new ConcurrentFlusher(chunkPw)));
					processedLines = 0;
					currentSparseLines.clear();
					chunkCounter++;
				}
			}
			br.close();
			executor.shutdown();
		
			//writing last chunk file
			br = new BufferedRandomAccessFile(new File(inputPath),"r",500);
			chunkPw = new UnsynchronizedBufferedWriter(new FileWriter(new File(String.format(chunkFolderPath + "/chunk_%s.txt",chunkCounter))));
			Collections.sort(currentSparseLines,lineComparator);
			for(SparseLine sparseLine : currentSparseLines) {
				br.seek(sparseLine.getFilePointer());
				chunkPw.write(br.getNextLine());
				chunkPw.newLine();
			}
			br.close();
			br = null;
			chunkPw.flush();
			chunkPw.close();
			chunkPw = null;
			
			currentSparseLines.clear();
			currentSparseLines = null;
			linePool.clear();
			linePool = null;
			
			//waiting for last write operations
			for(Future future : futures) {
				future.get();
			}
			futures.clear();
			futures = null;
			executor = null;
			
			for(UnsynchronizedBufferedWriter tmpPw : chunkPws) {
				tmpPw.flush();
				tmpPw.close();
			}
			
			chunkPws.clear();
			chunkPws = null;
			
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
		
		public SparseLine(String sortKey, long filePointer) {
			this.sortKey = sortKey;
			this.filePointer = filePointer;
		}
		
		public void update(String sortKey, long filePointer) {
			this.sortKey = sortKey;
			this.filePointer = filePointer;
		}
		
		public String getSortKey() {
			return this.sortKey;
		}
		
		public long getFilePointer() {
			return this.filePointer;
		}
	}
	
	private class ConcurrentFlusher extends Thread {
		
		private UnsynchronizedBufferedWriter pw;
		
		public ConcurrentFlusher(UnsynchronizedBufferedWriter pw) {
			super();
			this.pw = pw;
		}
		
		
		public void run() {
			try {
				this.pw.flush();
				this.pw = null;
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	

private class ConcurrentLineWriter extends Thread {
		
		private UnsynchronizedBufferedWriter pw;
		private ArrayList<String> lines;
		
		public ConcurrentLineWriter(UnsynchronizedBufferedWriter pw,ArrayList<String> lines) {
			super();
			this.pw = pw;
			this.lines = lines;
		}
		
		
		public void run() {
			try {
				for(String line : this.lines) {
					this.pw.write(line);
					this.pw.newLine();
				}
				this.pw.flush();
				this.pw = null;
				this.lines.clear();
				this.lines = null;
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	}


	public void mergeFiles(File folderToMerge, String outputPath, int maxChunks) {
		ArrayList<File> filesToMerge = new ArrayList<File>(Arrays.asList(folderToMerge.listFiles()));
		mergeFiles(filesToMerge,outputPath,maxChunks);
	}

	public void mergeFiles(ArrayList<File> filesToMerge, String outputPath, int maxChunks) {
		try {
			File outputFile = new File(outputPath);
			//in case there are too many chunk files we process them here again in chunks
			int fileCount = getFileCount(filesToMerge);
			int subsetIndex = 0;
			int deletionIndex = 0;
			String subsetOutputPath;
			while(fileCount > maxChunks) {
				subsetOutputPath = String.format("%s/merged_subset_%s.txt",filesToMerge.get(0).getParent(),subsetIndex);
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
			ArrayList<Line> currentRows = new ArrayList<Line>();
			long filePointer;
			String currentLine;
			for(int i = 0; i < filesToMerge.size(); i++) {
				if(filesToMerge.get(i).isFile() && !filesToMerge.get(i).getAbsolutePath().equals(outputFile.getAbsolutePath())) {
					BufferedRandomAccessFile currentReader = new BufferedRandomAccessFile(filesToMerge.get(i),"r",10240);
					readers.add(currentReader);
					filePointer= currentReader.getFilePointer();
					if((currentLine = currentReader.getNextLine()) != null) {
						Line tmpLine = new Line(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
						tmpLine.setReaderIndex(readers.size() - 1);
						currentRows.add(tmpLine);
						
					}
					else {
						currentReader.close();
						readers.remove(currentReader);
					}
				}
			}
			MinPriorityQueue queue = new MinPriorityQueue(currentRows);
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(outputFile),(int)(this.chunkSize * 1.0) * 1024 * 1024);
			int readerIndex;
			Line maxLine;
			int processedLines = 0;
			ExecutorService executor = Executors.newFixedThreadPool(1);
			ArrayList<Future> futures = new ArrayList<Future>();
			ArrayList<String> currentLines = new ArrayList<String>();
			while((maxLine = (Line)queue.extractMinimum()) != null) {
				readerIndex = maxLine.getReaderIndex();
				filePointer = readers.get(readerIndex).getFilePointer();
				readers.get(readerIndex).seek(maxLine.getFilePointer());
				currentLines.add(readers.get(readerIndex).getNextLine());
				readers.get(readerIndex).seek(filePointer);
				if((currentLine = readers.get(readerIndex).getNextLine()) != null) {
					maxLine.updateLine(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
					maxLine.setReaderIndex(readerIndex);
					queue.insert(maxLine);
				}
				processedLines++;
				
				if(currentLines.size() >= (this.linesPerChunkFile/10)) {
					if(!futures.isEmpty() && futures.get(0).isDone())
						futures.remove(0);
					futures.add(executor.submit(new ConcurrentLineWriter(pw,currentLines)));
					currentLines = new ArrayList<String>();
				}
				
				if(futures.size() >=10) {
					for(Future future : futures)
						future.get();
					
					futures.clear();
				}
			}
			futures.add(executor.submit(new ConcurrentLineWriter(pw,currentLines)));
			executor.shutdown();
			for(Future future : futures)
				future.get();
			
			futures.clear();
			futures = null;
			executor = null;
			queue = null;
			
			for(BufferedRandomAccessFile reader : readers) {
				reader.close();
				reader = null;
			}
			readers.clear();
			readers = null;
			pw.flush();
			pw.close();
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
			ArrayList<Line> currentRows = new ArrayList<Line>();
			long filePointer;
			String currentLine;
			for(int i = 0; i < filesToMerge.size() && i < maxChunks; i++) {
				if(filesToMerge.get(i).isFile() && !filesToMerge.get(i).getName().equals(outputFile.getName())) {
					BufferedRandomAccessFile currentReader = new BufferedRandomAccessFile(filesToMerge.get(i),"r",10240);
					readers.add(currentReader);
					filePointer = currentReader.getFilePointer();
					if((currentLine = currentReader.getNextLine()) != null) {
						Line tmpLine = new Line(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
						tmpLine.setReaderIndex(readers.size() - 1);
						currentRows.add(tmpLine);
						
					}
					else {
						currentReader.close();
						readers.remove(currentReader);
					}
				}
			}
			MinPriorityQueue queue = new MinPriorityQueue(currentRows);
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(outputFile),(int)(this.chunkSize * 1.0) * 1024 * 1024);
			int readerIndex;
			Line maxLine;
			int processedLines = 0;
			ExecutorService executor = Executors.newFixedThreadPool(1);
			ArrayList<Future> futures = new ArrayList<Future>();
			ArrayList<String> currentLines = new ArrayList<String>();
			while((maxLine = (Line)queue.extractMinimum()) != null) {
				readerIndex = maxLine.getReaderIndex();
				filePointer = readers.get(readerIndex).getFilePointer();
				readers.get(readerIndex).seek(maxLine.getFilePointer());
				currentLines.add(readers.get(readerIndex).getNextLine());
				readers.get(readerIndex).seek(filePointer);
				if((currentLine = readers.get(readerIndex).getNextLine()) != null) {
					maxLine.updateLine(currentLine,filePointer,this.columnsOfInterest,this.delimiter,this.numericalSort);
					maxLine.setReaderIndex(readerIndex);
					queue.insert(maxLine);
				}
				processedLines++;
				
				if(currentLines.size() >= (this.linesPerChunkFile/10)) {
					if(!futures.isEmpty() && futures.get(0).isDone())
						futures.remove(0);
					
					futures.add(executor.submit(new ConcurrentLineWriter(pw,currentLines)));
					currentLines = new ArrayList<String>();
				}
				
				if(futures.size() >=10) {
					for(Future future : futures)
						future.get();
					
					futures.clear();
				}
			}
			futures.add(executor.submit(new ConcurrentLineWriter(pw,currentLines)));
			executor.shutdown();
			for(Future future : futures)
				future.get();
			
			futures.clear();
			futures = null;
			queue = null;
			for(BufferedRandomAccessFile reader : readers) {
				reader.close();
				reader = null;
			}
			readers.clear();
			readers = null;
			pw.flush();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	

	private class Line implements Comparable  {
		private String sortKey;
		private int readerIndex;
		private long filePointer;
		private boolean numericalSort;
		
		public Line(String completeLine, long filePointer, HashSet<Integer> columnsToSortFor, String splitPattern, boolean numericalSort) {
			this.filePointer = filePointer;
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
