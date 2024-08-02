package alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.util.ArrayList;

import tools.SamProcessor;
import main.StreamConsumer;

public class ReadAlignerBWA implements ReadAligner {

	private final boolean useAllThreads = true;
	private final boolean mdFlagPreprocessed = false;
	private String tmpOutputDirPath;
	private String referencesDirPath;
	
	public ReadAlignerBWA(String tmpOutputDirPath, String referencesDirPath) {
		this.tmpOutputDirPath = tmpOutputDirPath;
		this.referencesDirPath = referencesDirPath;
		
	}
	
	public void alignReadsInitialPhase(String bwaExecutablePath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, String unalignedReadsOutputFilePath, String splitCandidatesOutputDir,String multiSplitCandidatesOutputDir, int seedLength, int seedMismatches, int[] splitSeedSizes, int maxMismatches,
			int maxHits, int maxReadLength, int minReadLength, int threads, boolean filterFullAlignments, boolean skipSplitDetection, boolean skipMultiSplitDetection, boolean verbose) {
		
		try {
			
			Process bwa;
			int returnValue;
			StreamConsumerInitialPhase errorOut;
			StreamConsumerInitialPhase stdOut;
			String optionalParameters = "";
			
			/**
			 * Here comes the external alignment program call
			 */
			bwa = Runtime.getRuntime().exec(String.format("%s mem -t %s -k %s -w 0 -a -c %s -A 1 -B 4 -O 10000 -E 10000 -L %s -r 1.0 -v 0 -T %s %s %s %s",bwaExecutablePath,threads,seedLength,maxHits,4,seedLength,optionalParameters,genomeIndexBasePath,fastaFilePath));
			
			// COPY AND PASTE FOR THE INTEGRATION OF OTHER ALIGNERS:
			errorOut = new StreamConsumerInitialPhase(bwa.getErrorStream(),StreamType.ERROR,maxMismatches, seedMismatches, splitSeedSizes,seedLength,threads,outputFilePath,unalignedReadsOutputFilePath,splitCandidatesOutputDir,multiSplitCandidatesOutputDir,filterFullAlignments,verbose,mdFlagPreprocessed, skipSplitDetection,skipMultiSplitDetection,useAllThreads);
			stdOut = new StreamConsumerInitialPhase(bwa.getInputStream(),StreamType.STDOUT,maxMismatches, seedMismatches, splitSeedSizes,seedLength,threads,outputFilePath,unalignedReadsOutputFilePath,splitCandidatesOutputDir,multiSplitCandidatesOutputDir,filterFullAlignments,verbose,mdFlagPreprocessed, skipSplitDetection,skipMultiSplitDetection,useAllThreads);
								
			errorOut.start();
			stdOut.start();
			returnValue = bwa.waitFor();
			//notify stream threads that the process terminated
			errorOut.processTerminated();
			stdOut.processTerminated();
			
			synchronized(errorOut) {
				while(errorOut.isAlive()) {
					errorOut.wait();
				}
			}
			
			synchronized(stdOut) {
				while(stdOut.isAlive()) {
					stdOut.wait();
				}
			}
			
			bwa.getErrorStream().close();
			bwa.getOutputStream().close();
			bwa.getInputStream().close();
			bwa.destroy();
			
			if(returnValue != 0) {
			throw new Exception(String.format("Something went wrong, bwa return value: %s.",returnValue));
			}
			
	}
		
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	
	public void alignReadsSlidingWindow(String bwaExecutablePath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, int seedLength, int seedMismatches, int maxMismatches,
			int maxHits, int maxReadLength, int minReadLength, int threads, boolean filterFullAlignments, boolean skipSplitDetection, boolean alignToFwdStrand, boolean alignToRevStrand, boolean verbose) {
		
		try {
			
			if(!new File(this.tmpOutputDirPath).exists()) {
				new File(this.tmpOutputDirPath).mkdir();
			}
			
			Process bwa;
			int returnValue;
			StreamConsumerSlidingWindow errorOut;
			StreamConsumerSlidingWindow stdOut;
			String currentTmpOutputDirPath = this.tmpOutputDirPath;
			String optionalParameters = "";
			
			currentTmpOutputDirPath = this.tmpOutputDirPath + "/" + genomeIndexBasePath.substring(genomeIndexBasePath.lastIndexOf('/'));
			new File(currentTmpOutputDirPath).mkdirs();
			bwa = Runtime.getRuntime().exec(String.format("%s mem -t %s -k %s -w 0 -a -c %s -A 1 -B 4 -O 10000 -E 10000 -r 1.0 -v 0 -T %s -L %s %s %s %s",bwaExecutablePath,threads,seedLength,maxHits,Integer.MIN_VALUE,5,optionalParameters,genomeIndexBasePath,fastaFilePath));
			errorOut = new StreamConsumerSlidingWindow(bwa.getErrorStream(),StreamType.ERROR,maxMismatches,threads,outputFilePath,filterFullAlignments,verbose,mdFlagPreprocessed);
			stdOut = new StreamConsumerSlidingWindow(bwa.getInputStream(),StreamType.STDOUT,maxMismatches,threads,outputFilePath,filterFullAlignments,verbose,mdFlagPreprocessed);
			
			errorOut.start();
			stdOut.start();
			returnValue = bwa.waitFor();
			//notify stream threads that the process terminated
			errorOut.processTerminated();
			stdOut.processTerminated();
			
			synchronized(errorOut) {
				while(errorOut.isAlive()) {
					errorOut.wait();
				}
			}
			
			synchronized(stdOut) {
				while(stdOut.isAlive()) {
					stdOut.wait();
				}
			}
			
			bwa.getErrorStream().close();
			bwa.getOutputStream().close();
			bwa.getInputStream().close();
			bwa.destroy();
			
			if(returnValue != 0) {
			throw new Exception(String.format("Something went wrong, bwa return value: %s.",returnValue));
			}
		}
		
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public void buildIndex(String bwaExecutablePath, String fastaFilePath, String outputFileBase) {
		try {
			//Process bwaBuild = Runtime.getRuntime().exec(String.format("%s index -p %s -a bwtsw %s",bwaExecutablePath,outputFileBase,fastaFilePath));
			Process bwaBuild = Runtime.getRuntime().exec(String.format("%s index -p %s -a is %s",bwaExecutablePath,outputFileBase,fastaFilePath));
			StreamConsumer errorOut = new StreamConsumer(bwaBuild.getErrorStream(),"err",false);
			StreamConsumer stdOut = new StreamConsumer(bwaBuild.getInputStream(),"stdout",false);
			errorOut.start();
			stdOut.start();
			int returnValue = bwaBuild.waitFor();
			
			//notify stream threads that the process terminated
			errorOut.processTerminated();
			stdOut.processTerminated();
			synchronized(errorOut) {
				while(errorOut.isAlive()) {
					errorOut.wait();
				}
			}
			synchronized(stdOut) {
				while(stdOut.isAlive()) {
					stdOut.wait();
				}
			}
			
			//close everything
			bwaBuild.getErrorStream().close();
			bwaBuild.getOutputStream().close();
			bwaBuild.getInputStream().close();
			bwaBuild.destroy();
			
			if(returnValue != 0) {
			throw new Exception(String.format("Something went wrong, bowtie return value: %s.",returnValue));
			} 
		}
		catch(Exception e) {
			e.printStackTrace();
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
