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


/**
 * This class demonstrates how a read alignment program that is not determining a MD field in the SAM output can be integrated into ContextMap
 * @author bonfert
 *
 */

public class ReadAlignerMissingMDField implements ReadAligner {

	private final boolean useAllThreads = true;
	private boolean mdFlagPreprocessed;
	private String tmpOutputDirPath;
	private String referencesDirPath;
	
	public ReadAlignerMissingMDField(String tmpOutputDirPath, String referencesDirPath) {
		this.tmpOutputDirPath = tmpOutputDirPath;
		this.referencesDirPath = referencesDirPath;
		this.mdFlagPreprocessed = true;
		
	}
	
	public void alignReadsInitialPhase(String bwaExecutablePath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, String unalignedReadsOutputFilePath, String splitCandidatesOutputDir,String multiSplitCandidatesOutputDir, int seedLength, int seedMismatches, int[] splitSeedSizes, int maxMismatches,
			int maxHits, int maxReadLength, int minReadLength, int threads, boolean filterFullAlignments, boolean skipSplitDetection, boolean skipMultiSplitDetection, boolean verbose) {
		
		try {
			
			if(!new File(this.tmpOutputDirPath).exists()) {
				new File(this.tmpOutputDirPath).mkdir();
			}
			
			Process bwa;
			int returnValue;
			StreamConsumerMissingFields errorOut;
			StreamConsumerMissingFields stdOut;
			String currentTmpOutputDirPath = this.tmpOutputDirPath;
			String optionalParameters = "";
			
			
			/**
			 * Here comes the external alignment program call
			 */
			bwa = Runtime.getRuntime().exec(String.format("%s mem -t %s -k %s -w 0 -a -c %s -A 1 -B 4 -O 10000 -E 10000 -L %s -r 1.0 -v 0 -T %s %s %s %s",bwaExecutablePath,threads,seedLength,maxHits,4,seedLength,optionalParameters,genomeIndexBasePath,fastaFilePath));

			
			
			/**
			 * If the underlying alignment program does not contain the MD-field (old versions of BWA), 
			 * we use the following StreamConsumer to directly stream the alignments to disk
			 * The missing MD field will be determined afterwards (see a few lines below). 
			 */
	
			// COPY AND PASTE FOR THE INTEGRATION OF OTHER ALIGNERS which do hard clipping or do not determine the MD field:
			errorOut = new StreamConsumerMissingFields(bwa.getErrorStream(),StreamType.ERROR,threads,currentTmpOutputDirPath,unalignedReadsOutputFilePath,filterFullAlignments,skipMultiSplitDetection,verbose);
			stdOut = new StreamConsumerMissingFields(bwa.getInputStream(),StreamType.STDOUT,threads,currentTmpOutputDirPath,unalignedReadsOutputFilePath,filterFullAlignments,skipMultiSplitDetection,verbose);

			
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
			
			
			

			String tmpOutputFilePath = outputFilePath + ".tmp";
			
			
			
			
			/**
			 * Determining MD field. Can be used for every underlying alignment program, that is not outputting it.
			 */
			
			boolean hasMdFlag = false;
			if(new File(currentTmpOutputDirPath).listFiles().length > 0) {
				BufferedReader checkReader = new BufferedReader(new FileReader(new File(currentTmpOutputDirPath).listFiles()[0]));
				String currentLine;
				while((currentLine = checkReader.readLine()) != null) {
					if(currentLine.charAt(0) == '@')
						continue;
					
					if(currentLine.split("\t")[2].equals("*"))
						continue;
					
					if(currentLine.contains("MD:Z:"))
						hasMdFlag = true;
					
					break;
				}
				checkReader.close();
		    }
			
			
			if(hasMdFlag) {
				this.mdFlagPreprocessed = false;
				ArrayList<String> filePaths = new ArrayList<String>();
				for(File f : new File(currentTmpOutputDirPath).listFiles()) {
					filePaths.add(f.getAbsolutePath());
				}
				concatenateFilesWithNIO(filePaths,tmpOutputFilePath);
			}
			
			else {
				SamProcessor samProcessor = new SamProcessor();
				samProcessor.addMdFlag(currentTmpOutputDirPath, tmpOutputFilePath, new File(this.referencesDirPath).listFiles(),true, threads);
				
			}
			
			
			
			
			/**
			 * Now we are ready to start the same StreamConsumer as already used for Bowtie 1 and Bowtie 2 (this time parsing from disk):
			 */
			SamStreamConsumer samStreamConsumer;
			InputStream is = new FileInputStream(new File(tmpOutputFilePath));
			//unaligned reads have already be written to disk
			unalignedReadsOutputFilePath = null;
			samStreamConsumer = new StreamConsumerInitialPhase(is, StreamType.STDOUT,maxMismatches, seedMismatches, splitSeedSizes, seedLength, threads, outputFilePath, unalignedReadsOutputFilePath, splitCandidatesOutputDir,multiSplitCandidatesOutputDir, filterFullAlignments, verbose,mdFlagPreprocessed, skipSplitDetection,skipMultiSplitDetection,useAllThreads);
			
			
			
			samStreamConsumer.start();
			samStreamConsumer.processTerminated();
			
			synchronized(samStreamConsumer) {
				while(samStreamConsumer.isAlive()) {
					samStreamConsumer.wait();
				}
			}
			
			is.close();
			new File(tmpOutputFilePath).delete();
			File[] samFiles = new File(currentTmpOutputDirPath).listFiles();
			for(File samFile : samFiles)
				samFile.delete();
			
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
			StreamConsumerMissingFields errorOut;
			StreamConsumerMissingFields stdOut;
			String currentTmpOutputDirPath = this.tmpOutputDirPath;
			String optionalParameters = "";
			
			currentTmpOutputDirPath = this.tmpOutputDirPath + "/" + genomeIndexBasePath.substring(genomeIndexBasePath.lastIndexOf('/'));
			new File(currentTmpOutputDirPath).mkdirs();
			bwa = Runtime.getRuntime().exec(String.format("%s mem -t %s -k %s -w 0 -a -c %s -A 1 -B 4 -O 10000 -E 10000 -r 1.0 -v 0 -T %s -L %s %s %s %s",bwaExecutablePath,threads,seedLength,maxHits,Integer.MIN_VALUE,5,optionalParameters,genomeIndexBasePath,fastaFilePath));
			errorOut = new StreamConsumerMissingFields(bwa.getErrorStream(),StreamType.ERROR,threads,currentTmpOutputDirPath,null,filterFullAlignments,false,verbose);
			stdOut = new StreamConsumerMissingFields(bwa.getInputStream(),StreamType.STDOUT,threads,currentTmpOutputDirPath,null,filterFullAlignments,false,verbose);
		
			
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
			
			
			String tmpOutputFilePath = outputFilePath + ".tmp";
			
			
			boolean hasMdFlag = false;
			if(new File(currentTmpOutputDirPath).listFiles().length > 0) {
				BufferedReader checkReader = new BufferedReader(new FileReader(new File(currentTmpOutputDirPath).listFiles()[0]));
				String currentLine;
				while((currentLine = checkReader.readLine()) != null) {
					if(currentLine.charAt(0) == '@')
						continue;
					
					if(currentLine.split("\t")[2].equals("*"))
						continue;
					
					if(currentLine.contains("MD:Z:"))
						hasMdFlag = true;
					
					break;
				}
				checkReader.close();
		    }
			
			
			if(hasMdFlag) {
				this.mdFlagPreprocessed = false;
				ArrayList<String> filePaths = new ArrayList<String>();
				for(File f : new File(currentTmpOutputDirPath).listFiles()) {
					filePaths.add(f.getAbsolutePath());
				}
				concatenateFilesWithNIO(filePaths,tmpOutputFilePath);
			}
			
			else {
				SamProcessor samProcessor = new SamProcessor();
				File[] refFiles = new File[1];
				refFiles[0] = new File(genomeIndexBasePath.substring(0,genomeIndexBasePath.lastIndexOf('/') + 1) + genomeIndexBasePath.substring(genomeIndexBasePath.lastIndexOf('/') + 1) + ".fa");
				samProcessor.addMdFlag(currentTmpOutputDirPath, tmpOutputFilePath, refFiles,true, threads);
				
				
			}
			
			
			
			SamStreamConsumer samStreamConsumer;
			InputStream is = new FileInputStream(new File(tmpOutputFilePath));
			samStreamConsumer = new StreamConsumerSlidingWindow(is,StreamType.STDOUT,maxMismatches,threads,outputFilePath,filterFullAlignments,verbose,mdFlagPreprocessed);
			
			samStreamConsumer.start();
			samStreamConsumer.processTerminated();
			
			synchronized(samStreamConsumer) {
				while(samStreamConsumer.isAlive()) {
					samStreamConsumer.wait();
				}
			}
			
			is.close();
			new File(tmpOutputFilePath).delete();
			File[] samFiles = new File(currentTmpOutputDirPath).listFiles();
			for(File samFile : samFiles)
				samFile.delete();
			
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
