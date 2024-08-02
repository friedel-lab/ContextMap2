package tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.StringTokenizer;


public class FileSorter extends Thread {

	
	private String inputPath;
	private String outputPath;
	private String tmpFolderPath;
	
	private HashSet<Integer> columnsOfInterest;
	private int linesPerChunkFile;
	private final int maxChunks = 50;
	private String delimiter;
	private boolean numericalSort;
	
	public static void main(String args[]) {
		int[] columnsToSortFor = null;
		int linesPerChunkFile = 500000;
		String delimiter = "\t";
		String inputPath = null;
		String outputPath = null;
		boolean numericalSort = false;
		
		if(args.length == 0 || args[0].equals("-h") || args[0].equals("--help")) {
			System.out.println("Sorts a file by specific column(s).");
			System.out.println("Available parameters:");
			System.out.println("-c\tcomma separated list of columns to sort for (0-based)");
			System.out.println("-d\tdelimiter which is used in the input file <default: tab>");
			System.out.println("-i\tpath to input file");
			System.out.println("-o\tpath to output file");
			System.out.println("-linesperchunk\tnumber of lines per chunk file <default: 500000>");
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
			if(args[i].contains("-linesperchunk")) {
				linesPerChunkFile = Integer.valueOf(args[++i]);
				continue;
			}
			if(args[i].contains("--numerical")) {
				numericalSort = true;
				continue;
			}
		}
		FileSorter fileSorter = new FileSorter(columnsToSortFor,linesPerChunkFile,delimiter,numericalSort);
		fileSorter.sortFile(inputPath, outputPath, null);
	}
	
	
	public FileSorter(int[] columnsToSortFor,int linesPerChunkFile, String delimiter,boolean numericalSort) {
		this.columnsOfInterest = new HashSet<Integer>();
		for(int i = 0; i < columnsToSortFor.length; i++ ) {
			this.columnsOfInterest.add(columnsToSortFor[i]);
		}
		this.linesPerChunkFile = linesPerChunkFile;
		this.delimiter = delimiter;
		this.numericalSort = numericalSort;
	}
	
	
	public FileSorter(String inputPath, String outputPath, String tmpFolderPath, int[] columnsToSortFor,int linesPerChunkFile, String delimiter,boolean numericalSort) {
		this.inputPath = inputPath;
		this.outputPath = outputPath;
		this.tmpFolderPath = tmpFolderPath;
		this.columnsOfInterest = new HashSet<Integer>();
		for(int i = 0; i < columnsToSortFor.length; i++ ) {
			this.columnsOfInterest.add(columnsToSortFor[i]);
		}
		this.linesPerChunkFile = linesPerChunkFile;
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
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			String chunkFolderPath = tmpFolderPath;
			if(tmpFolderPath == null)
				chunkFolderPath = outputPath.substring(0,outputPath.lastIndexOf(System.getProperty("file.separator"))) + "/tmp";
			File chunkFolder = new File(chunkFolderPath);
			if(!chunkFolder.isDirectory())
				chunkFolder.mkdirs();
			
			
			Line[] currentChunk = new Line[this.linesPerChunkFile];
			
			Comparator lineComparator;
			if(!this.numericalSort)
				lineComparator = new LineComparatorForStringKeys();
			else
				lineComparator = new LineComparatorForIntegerKeys();
			int processedLines = 0;
			int chunkCounter = 0;
			BufferedWriter chunkPw;
			while(br.ready()) {
				if(currentChunk[processedLines] == null)
					currentChunk[processedLines] = new Line(br.readLine(),this.columnsOfInterest,this.delimiter,this.numericalSort);
				else
					currentChunk[processedLines].updateLine(br.readLine(),this.columnsOfInterest,this.delimiter,this.numericalSort);
				
				if(++processedLines == this.linesPerChunkFile) {
					//sort and write to chunk
					Arrays.sort(currentChunk, lineComparator);
					chunkPw = new BufferedWriter(new FileWriter(new File(String.format(chunkFolderPath + "/chunk_%s.txt",chunkCounter))));
					for(int i = 0; i < currentChunk.length; i++) {
						chunkPw.write(currentChunk[i].getCompleteLine());
						chunkPw.newLine();
					}
					chunkPw.close();
					processedLines = 0;
					chunkCounter++;
				}
			}
			br.close();
			//writing last chunk file
			Line[] lastChunk = Arrays.copyOf(currentChunk, processedLines);
			Arrays.sort(lastChunk, lineComparator);
			chunkPw = new BufferedWriter(new FileWriter(new File(String.format(chunkFolderPath + "/chunk_%s.txt",chunkCounter))));
			for(int i = 0; i < lastChunk.length; i++) {
				chunkPw.write(lastChunk[i].getCompleteLine());
				chunkPw.newLine();
			}
			chunkPw.close();
			currentChunk = null;
			lastChunk = null;
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
			
			
			
			ArrayList<BufferedReader> readers = new ArrayList<BufferedReader>();
			ArrayList<Line> currentRows = new ArrayList<Line>();
			for(int i = 0; i < filesToMerge.size(); i++) {
				if(filesToMerge.get(i).isFile() && !filesToMerge.get(i).getName().equals(outputFile.getName())) {
					BufferedReader currentReader = new BufferedReader(new FileReader(filesToMerge.get(i)));
					readers.add(currentReader);
					if(currentReader.ready()) {
						Line tmpLine = new Line(currentReader.readLine(),this.columnsOfInterest,this.delimiter,this.numericalSort);
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
			BufferedWriter pw = new BufferedWriter(new FileWriter(outputFile));
			int readerIndex;
			Line maxLine;
			while((maxLine = (Line)queue.extractMinimum()) != null) {
				pw.write(maxLine.getCompleteLine());
				pw.newLine();
				readerIndex = maxLine.getReaderIndex();
				if(readers.get(readerIndex).ready()) {
					maxLine.updateLine(readers.get(readerIndex).readLine(),this.columnsOfInterest,this.delimiter,this.numericalSort);
					maxLine.setReaderIndex(readerIndex);
					queue.insert(maxLine);
				}
			}
			for(BufferedReader reader : readers)
				reader.close();
			
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
			
			
			ArrayList<BufferedReader> readers = new ArrayList<BufferedReader>();
			ArrayList<Line> currentRows = new ArrayList<Line>();
			for(int i = 0; i < filesToMerge.size() && i < maxChunks; i++) {
				if(filesToMerge.get(i).isFile() && !filesToMerge.get(i).getName().equals(outputFile.getName())) {
					BufferedReader currentReader = new BufferedReader(new FileReader(filesToMerge.get(i)));
					readers.add(currentReader);
					if(currentReader.ready()) {
						Line tmpLine = new Line(currentReader.readLine(),this.columnsOfInterest,this.delimiter,this.numericalSort);
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
			BufferedWriter pw = new BufferedWriter(new FileWriter(outputFile));
			int readerIndex;
			Line maxLine;
			while((maxLine = (Line)queue.extractMinimum()) != null) {
				pw.write(maxLine.getCompleteLine());
				pw.newLine();
				readerIndex = maxLine.getReaderIndex();
				if(readers.get(readerIndex).ready()) {
					maxLine.updateLine(readers.get(readerIndex).readLine(),this.columnsOfInterest,this.delimiter,this.numericalSort);
					maxLine.setReaderIndex(readerIndex);
					queue.insert(maxLine);
				}
			}
			
			for(BufferedReader reader : readers)
				reader.close();
			
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	

	private class Line implements Comparable  {
		private String completeLine;
		private String sortKey;
		private int readerIndex;
		private boolean numericalSort;
		
		public Line(String completeLine, HashSet<Integer> columnsToSortFor, String splitPattern, boolean numericalSort) {
			this.completeLine = completeLine;
			StringTokenizer tokenizer;
			if(splitPattern.equals("\\t"))
				tokenizer = new StringTokenizer(this.completeLine, "\t");
			else
				tokenizer = new StringTokenizer(this.completeLine, splitPattern);
			
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
		
		public void updateLine(String completeLine, HashSet<Integer> columnsToSortFor, String splitPattern, boolean numericalSort) {
			this.completeLine = completeLine;
			StringTokenizer tokenizer;
			if(splitPattern.equals("\\t"))
				tokenizer = new StringTokenizer(this.completeLine, "\t");
			else
				tokenizer = new StringTokenizer(this.completeLine, splitPattern);
			
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
		
		public String getCompleteLine() {
			return this.completeLine;
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
				return(Integer.valueOf(this.sortKey).compareTo(Integer.valueOf(((Line)o).getSortKey())));
			else 
				return(this.sortKey.compareTo(((Line)o).getSortKey()));
		}
	}
	
	private class LineComparatorForStringKeys implements Comparator<Line> {
		public int compare(Line l1, Line l2) {
			return(l1.getSortKey().compareTo(l2.getSortKey()));
		}

	}
	
	private class LineComparatorForIntegerKeys implements Comparator<Line> {
		public int compare(Line l1, Line l2) {
			return(Integer.valueOf(l1.getSortKey()).compareTo(Integer.valueOf(l2.getSortKey())));
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
