package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import main.StreamConsumer;

public class UnixSort extends Thread {
	
	private String inputFilePath;
	private String outputFilePath;
	private String tmpFolderPath;
	private String delimiter;
	private int column;
	private int maxMem;
	private boolean numericalSort;
	private boolean deleteInputFile;
	private boolean verbose;
	
	private int[] columns;
	private String[] orderingOptions;
	
	public static void main(String args[]) {
		UnixSort sorter = new UnixSort(args[0],args[1],args[2],args[3],Integer.valueOf(args[4]),Integer.valueOf(args[5]),Boolean.valueOf(args[6]),false,Boolean.valueOf(args[7]));
		sorter.start();
	}

	public UnixSort(String inputFilePath, String outputFilePath, String tmpFolderPath, String delimiter, int column, int maxMem, boolean numericalSort, boolean deleteInputFile, boolean verbose) {
		this.inputFilePath = inputFilePath;
		this.outputFilePath = outputFilePath;
		this.tmpFolderPath = tmpFolderPath;
		this.delimiter = delimiter;
		this.column = column;
		this.maxMem = maxMem;
		this.numericalSort = numericalSort;
		this.deleteInputFile = deleteInputFile;
		this.verbose = verbose;
		
		this.columns = null;
		this.orderingOptions = null;
	}
	
	
	//TODO change everything to this constructor and remove the other one....
	public UnixSort(String inputFilePath, String outputFilePath, String tmpFolderPath, String delimiter, int[] columns,String[] orderingOptions, int maxMem,boolean deleteInputFile, boolean verbose) {
		this.inputFilePath = inputFilePath;
		this.outputFilePath = outputFilePath;
		this.tmpFolderPath = tmpFolderPath;
		this.delimiter = delimiter;
		this.columns = columns;
		this.maxMem = maxMem;
		this.orderingOptions = orderingOptions;
		this.deleteInputFile = deleteInputFile;
		this.verbose = verbose;
	}
	
	public UnixSort() {
		
	}
	
	public void run() {
		try {
			if(this.columns == null)
				sort(inputFilePath,outputFilePath,tmpFolderPath,delimiter, column, maxMem, numericalSort,deleteInputFile, verbose);
			else
				sort(inputFilePath,outputFilePath,tmpFolderPath,delimiter, columns, orderingOptions, maxMem,deleteInputFile, verbose);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void sort(String inputFilePath, String outputFilePath, String tmpFolderPath, String delimiter, int[] columns, String[] orderingOptions, int maxMem,boolean deleteInputFile, boolean verbose) {
		try {
			//first check if the tmpFolder exists
			File tmpFolder = new File(tmpFolderPath);
			if(!tmpFolder.isDirectory())
				tmpFolder.mkdirs();
			
			
			Process sort;
			List<String> command = new LinkedList<String>();
			command.add("sort");
			for(int i = 0; i < columns.length; i++) {
				command.add("-k" + columns[i] + orderingOptions[i] + "," + columns[i] + orderingOptions[i]);
			}
			command.add("-o" + outputFilePath);
			command.add("-t" + delimiter);
			command.add("-T" + tmpFolderPath);
			command.add("-S" + maxMem + "M");
			command.add(inputFilePath);
			sort = new ProcessBuilder(command).start();
			
			//handle output
			StreamConsumer errorOut = new StreamConsumer(sort.getErrorStream(),"err",true);
			StreamConsumer stdOut = new StreamConsumer(sort.getInputStream(),"stdout",true);
			errorOut.start();
			stdOut.start();
			int returnValue = sort.waitFor();
			
			//notify stream threads that the process terminated
			errorOut.processTerminated();
			stdOut.processTerminated();
			
			//close everything
			sort.getErrorStream().close();
			sort.getOutputStream().close();
			sort.getInputStream().close();
			sort.destroy();
			
			//remove tmp folder
			new File(tmpFolderPath).delete();
			
			//remove inputFile
			if(deleteInputFile)
				new File(inputFilePath).delete();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void sort(String inputFilePath, String outputFilePath, String tmpFolderPath, String delimiter, int column, int maxMem, boolean numericalSort,boolean deleteInputFile, boolean verbose) {
		try {
			//first check if the tmpFolder exists
			File tmpFolder = new File(tmpFolderPath);
			if(!tmpFolder.isDirectory())
				tmpFolder.mkdirs();
			
			
			Process sort;
			List<String> command = new LinkedList<String>();
			command.add("sort");
			command.add("-k" + column + "," + column);
			command.add("-o" + outputFilePath);
			command.add("-t" + delimiter);
			command.add("-T" + tmpFolderPath);
			command.add("-S" + maxMem + "M");
			if(numericalSort)
				command.add("-n");
			
			command.add(inputFilePath);
			sort = new ProcessBuilder(command).start();
			
			//handle output
			StreamConsumer errorOut = new StreamConsumer(sort.getErrorStream(),"err",true);
			StreamConsumer stdOut = new StreamConsumer(sort.getInputStream(),"stdout",true);
			errorOut.start();
			stdOut.start();
			int returnValue = sort.waitFor();
			
			//notify stream threads that the process terminated
			errorOut.processTerminated();
			stdOut.processTerminated();
			
			//close everything
			sort.getErrorStream().close();
			sort.getOutputStream().close();
			sort.getInputStream().close();
			sort.destroy();
			
			//remove tmp folder
			new File(tmpFolderPath).delete();
			
			//remove input file
			if(deleteInputFile)
				new File(inputFilePath).delete();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public void merge(String inputDirPath, String outputFilePath, String tmpFolderPath, String delimiter, int column, int maxMem, boolean numericalSort,boolean verbose) {
		try {
			//first check if the tmpFolder exists
			File tmpFolder = new File(tmpFolderPath);
			if(!tmpFolder.isDirectory())
				tmpFolder.mkdirs();
			
			
			Process sort;
			List<String> command = new LinkedList<String>();
			command.add("sort");
			command.add("-m");
			command.add("-k" + column + "," + column);
			command.add("-o" + outputFilePath);
			command.add("-t" + delimiter);
			command.add("-T" + tmpFolderPath);
			command.add("-S" + maxMem + "M");
			if(numericalSort)
				command.add("-n");
			
			//File[] files = new File(inputDirPath).listFiles();
			//for(File f : files) {
			//	command.add(f.getAbsolutePath());
			//}
			
			File[] files = new File(inputDirPath).listFiles();
			PrintWriter pw = new PrintWriter(new FileWriter(new File(tmpFolder.getAbsolutePath() + "/files_to_sort.txt")));
			for(File f : files) {
				pw.print(f.getAbsolutePath() + "\0");
			}
			pw.close();
			command.add("--files0-from");
			command.add(tmpFolder.getAbsolutePath() + "/files_to_sort.txt");
			
			
			sort = new ProcessBuilder(command).start();
			
			//handle output
			StreamConsumer errorOut = new StreamConsumer(sort.getErrorStream(),"err",true);
			StreamConsumer stdOut = new StreamConsumer(sort.getInputStream(),"stdout",true);
			errorOut.start();
			stdOut.start();
			int returnValue = sort.waitFor();
			
			//notify stream threads that the process terminated
			errorOut.processTerminated();
			stdOut.processTerminated();
			
			//close everything
			sort.getErrorStream().close();
			sort.getOutputStream().close();
			sort.getInputStream().close();
			sort.destroy();
			
			//remove tmp folder
			new File(tmpFolderPath).delete();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
}
