package alignment;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;

import main.StreamConsumer;



public class ReadAlignerBowtie1 implements ReadAligner {
	
	private final int baseQuality = 40;
	private final boolean mdFlagPreprocessed = false;
	private final boolean useAllThreads = false;

	public ReadAlignerBowtie1() {
		
	}
	
	/**
	 * 
	 * BOWTIE 1 calls
	 * 
	 */

		
	public void alignReadsInitialPhase(String bowtieExecutablePath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, String unalignedReadsOutputFilePath, String splitCandidatesOutputDir,String multiSplitCandidatesOutputDir, int seedLength, int seedMissmatches, int[] splitSeedSizes, int maxMismatches,
			int maxHits, int maxReadLength, int minReadLength, int threads, boolean filterFullAlignments, boolean skipSplitDetection, boolean skipMultiSplitDetection, boolean verbose) {
		
		try {
			Process bowtie;
			int returnValue;
			SamStreamConsumer errorOut;
			SamStreamConsumer stdOut;
			String optionalParameters = "";
			if(unalignedReadsOutputFilePath != null) {
				optionalParameters = String.format("--un %s",unalignedReadsOutputFilePath);
			}
			
			
			/**
			 * Here comes the external alignment program call:
			 */
			bowtie = Runtime.getRuntime().exec(String.format("%s -f -n %s -l %s -a -m %s %s -S --sam-nohead -e %s --maxbts 1250000 --nomaqround --quiet -p %s %s %s",bowtieExecutablePath,seedMissmatches,seedLength,maxHits,optionalParameters,((maxReadLength/2) + 1 + maxMismatches) * this.baseQuality + 1,threads,genomeIndexBasePath,fastaFilePath));
			
			
			/**
			 * The remaining work is done by the StreamConsumer class for the 'initial' alignment phase,
			 * which works for every underlying alignment program fulfilling the requirements mentioned in the ContextMap 2.0 paper.
			 * Here: Bowtie directly streams the unaligned reads to a separate file (see optional parameter), therefore we set the unalignedReadsOutputFilePath to null for the StreamConsumer
			 */
			
			// COPY AND PASTE FOR THE INTEGRATION OF OTHER ALIGNERS:
			unalignedReadsOutputFilePath = null;
			errorOut = new StreamConsumerInitialPhase(bowtie.getErrorStream(),StreamType.ERROR,maxMismatches, seedMissmatches, splitSeedSizes,seedLength,threads,outputFilePath,unalignedReadsOutputFilePath,splitCandidatesOutputDir,multiSplitCandidatesOutputDir,filterFullAlignments,verbose,mdFlagPreprocessed,skipSplitDetection,skipMultiSplitDetection,useAllThreads);
			stdOut = new StreamConsumerInitialPhase(bowtie.getInputStream(),StreamType.STDOUT,maxMismatches, seedMissmatches, splitSeedSizes,seedLength,threads,outputFilePath,unalignedReadsOutputFilePath,splitCandidatesOutputDir,multiSplitCandidatesOutputDir,filterFullAlignments,verbose,mdFlagPreprocessed,skipSplitDetection,skipMultiSplitDetection,useAllThreads);
		
			errorOut.start();
			stdOut.start();
			returnValue = bowtie.waitFor();
			
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
						
			bowtie.getErrorStream().close();
			bowtie.getOutputStream().close();
			bowtie.getInputStream().close();
			bowtie.destroy();
			
			if(returnValue != 0) {
			throw new Exception(String.format("Something went wrong, bowtie return value: %s.",returnValue));
			}	
		}
		
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public void alignReadsSlidingWindow(String bowtieExecutablePath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, int seedLength, int seedMissmatches, int maxMismatches,
			int maxHits, int maxReadLength, int minReadLength, int threads, boolean filterFullAlignments, boolean skipSplitDetection, boolean alignToFwdStrand, boolean alignToRevStrand, boolean verbose) {
		
		try {
			Process bowtie;
			int returnValue;
			SamStreamConsumer errorOut;
			SamStreamConsumer stdOut;
			
			String optionalParameters = "";
			
			if(alignToFwdStrand && !alignToRevStrand)
				optionalParameters = "--norc";
			
			else if(!alignToFwdStrand && alignToRevStrand)
				optionalParameters = "--nofw";
			
			/**
			 * Here comes the external alignment program call:
			 */
			bowtie = Runtime.getRuntime().exec(String.format("%s -f %s -n %s -l %s -a -m %s -S --sam-nohead -e %s --maxbts 1250000 --nomaqround --quiet -p %s %s %s",bowtieExecutablePath,optionalParameters,seedMissmatches,seedLength,maxHits,10000,threads,genomeIndexBasePath,fastaFilePath));
			
			
			/**
			 * The remaining work is done by the StreamConsumer class for the 'sliding window' alignment phase,
			 * which works for every underlying alignment program fulfilling the requirements mentioned in the ContextMap 2.0 paper:
			 */
			// COPY AND PASTE FOR THE INTEGRATION OF OTHER ALIGNERS:
			errorOut = new StreamConsumerSlidingWindow(bowtie.getErrorStream(),StreamType.ERROR,maxMismatches,threads,outputFilePath,filterFullAlignments,verbose,mdFlagPreprocessed);
			stdOut = new StreamConsumerSlidingWindow(bowtie.getInputStream(),StreamType.STDOUT,maxMismatches,threads,outputFilePath,filterFullAlignments,verbose,mdFlagPreprocessed);
		
			errorOut.start();
			stdOut.start();
			returnValue = bowtie.waitFor();
			
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
						
			bowtie.getErrorStream().close();
			bowtie.getOutputStream().close();
			bowtie.getInputStream().close();
			bowtie.destroy();
			
			if(returnValue != 0) {
			throw new Exception(String.format("Something went wrong, bowtie return value: %s.",returnValue));
			}	
		}
		
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public void buildIndex(String bowtieBuildExecutablePath, String fastaFilePath, String outputFileBase) {
		try {
			Process bowtieBuild = Runtime.getRuntime().exec(String.format("%s --noauto --bmaxdivn 1 --dcv 64 --noref -o 1 -f %s %s",bowtieBuildExecutablePath,fastaFilePath,outputFileBase));
			StreamConsumer errorOut = new StreamConsumer(bowtieBuild.getErrorStream(),"err",false);
			StreamConsumer stdOut = new StreamConsumer(bowtieBuild.getInputStream(),"stdout",false);
			errorOut.start();
			stdOut.start();
			int returnValue = bowtieBuild.waitFor();
			
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
			bowtieBuild.getErrorStream().close();
			bowtieBuild.getOutputStream().close();
			bowtieBuild.getInputStream().close();
			bowtieBuild.destroy();
			
			if(returnValue != 0) {
			throw new Exception(String.format("Something went wrong, bowtie return value: %s.",returnValue));
			} 
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	
}
