package alignment;

import main.StreamConsumer;

public class ReadAlignerBowtie2 implements ReadAligner {
	
	private final boolean mdFlagPreprocessed = false;
	private final boolean useAllThreads = false;
	
	public ReadAlignerBowtie2()  {
		
	}
	
	
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
			if(maxHits > 1)
				bowtie = Runtime.getRuntime().exec(String.format("%s -f -k %s -D 20 -R 2 -N %s -L %s --local --ignore-quals --ma 1 --mp 1,1 --np 1 --score-min L,-%s,0.2 --rdg 10000,10000 --rfg 10000,10000 %s -p %s --quiet --no-hd --no-sq -x %s %s",bowtieExecutablePath,maxHits,seedMissmatches,seedLength,maxMismatches,optionalParameters,threads,genomeIndexBasePath,fastaFilePath));
			else
				bowtie = Runtime.getRuntime().exec(String.format("%s -f -D 20 -R 2 -N %s -L %s --local --ignore-quals --ma 1 --mp 1,1 --np 1 --score-min L,-%s,0.2 --rdg 10000,10000 --rfg 10000,10000 %s -p %s --quiet --no-hd --no-sq -x %s %s",bowtieExecutablePath,seedMissmatches,seedLength,maxMismatches,optionalParameters,threads,genomeIndexBasePath,fastaFilePath));
			
			
			
			/**
			 * The remaining work is done by the StreamConsumer class for the 'initial' alignment phase,
			 * which works for every underlying alignment program fulfilling the requirements mentioned in the ContextMap 2.0 paper.
			 * Here: Bowtie2 directly streams the unaligned reads to a separate file (see optional parameter), therefore we set the unalignedReadsOutputFilePath to null for the StreamConsumer
			 */
			// COPY AND PASTE FOR THE INTEGRATION OF OTHER ALIGNERS:
			unalignedReadsOutputFilePath = null;
			errorOut = new StreamConsumerInitialPhase(bowtie.getErrorStream(),StreamType.ERROR,maxMismatches, seedMissmatches, splitSeedSizes,seedLength,threads,outputFilePath,unalignedReadsOutputFilePath,splitCandidatesOutputDir,multiSplitCandidatesOutputDir,filterFullAlignments,verbose,mdFlagPreprocessed, skipSplitDetection,skipMultiSplitDetection,useAllThreads);
			stdOut = new StreamConsumerInitialPhase(bowtie.getInputStream(),StreamType.STDOUT,maxMismatches, seedMissmatches, splitSeedSizes,seedLength,threads,outputFilePath,unalignedReadsOutputFilePath,splitCandidatesOutputDir,multiSplitCandidatesOutputDir,filterFullAlignments,verbose,mdFlagPreprocessed, skipSplitDetection,skipMultiSplitDetection,useAllThreads);
		
			
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
			if(maxHits > 1)
				bowtie = Runtime.getRuntime().exec(String.format("%s -f -k %s --end-to-end --ignore-quals --mp 1,1 --score-min L,0,-1 -i S,10000,1 --rdg 10000,10000 --rfg 10000,10000 -D 20 -R 1 -N %s -L %s %s -p %s --quiet --no-hd --no-sq -x %s %s",bowtieExecutablePath,maxHits,seedMissmatches,seedLength,optionalParameters,threads,genomeIndexBasePath,fastaFilePath));
			else
				bowtie = Runtime.getRuntime().exec(String.format("%s -f --end-to-end --ignore-quals --mp 1,1 --score-min L,0,-1 -i S,10000,1 --rdg 10000,10000 --rfg 10000,10000 -D 20 -R 1 -N %s -L %s %s -p %s --quiet --no-hd --no-sq -x %s %s",bowtieExecutablePath,seedMissmatches,seedLength,optionalParameters,threads,genomeIndexBasePath,fastaFilePath));
			
			
			/**
			 * The remaining work is done by the StreamConsumer class for the 'sliding window' alignment phase,
			 * which works for every underlying alignment program fulfilling the requirements mentioned in the ContextMap 2.0 paper.
			 * 
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
			Process bowtieBuild = Runtime.getRuntime().exec(String.format("%s --noauto --bmax 10000000 --dcv 64 -o 1 -f %s %s",bowtieBuildExecutablePath,fastaFilePath,outputFileBase));
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
			throw new Exception(String.format("Something went wrong, bowtie2 return value: %s.",returnValue));
			} 
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
}
