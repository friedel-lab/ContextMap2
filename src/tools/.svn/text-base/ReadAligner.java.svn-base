package tools;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;

import main.StreamConsumer;

/**
 * Aligns reads with "rmap", which is a modified bowtie implementation 
 * rmap command line usage: -i <fastafile> -index <bowtie genome index base path> -o <outfile> -tmpdir <tmpdir> -btbuild <path  to bowtie-build executable> -mm <max missmatches per read>
 *
 */
public class ReadAligner {
	
	private final int baseQuality = 40;

	public ReadAligner() {
		
	}

	public void alignReads(String rmapExecutablePath, String fastaFilePath, String genomeIndexBasePath,
							String outputFilePath, String tmpDirPath, String bowtieBuildExecutablePath,
							String seedLength, String seedMissmatches, String maxMissmatches, 
							String splitSeedSizes, String maxIntronLength, String maxIntronCount,
							String maxInsertionSize, int maxHits, String skippedReadsPath,
							String unalignedReadsPath, boolean prebufferReads, String threads,
							String sortMemory, boolean skipFull, boolean skipSplitDetection,boolean verbose) {
		try {
			
			String optionalParameters = String.format("-mm %s -l %s -n %s -splitseedsizes %s -maxintronlength %s -maxintroncount %s -maxdelsize %s -t %s -sMEM %s",maxMissmatches,seedLength,seedMissmatches,splitSeedSizes,maxIntronLength,maxIntronCount,maxInsertionSize,threads,sortMemory);
			if(maxHits != -1) optionalParameters += " -maxHits " + maxHits;
			if(skippedReadsPath != null) optionalParameters += " -skipped " + skippedReadsPath;
			if(unalignedReadsPath != null) optionalParameters += " -nohit " + unalignedReadsPath;
			if(prebufferReads) optionalParameters += " -prebuffer";
			if(skipFull) optionalParameters += " -skipfull";
			if(skipSplitDetection) optionalParameters += " -pass 0 -single";
			
			if(verbose) System.out.println(String.format("./%s -i %s -index %s -o %s -tmpdir %s -btbuild %s %s",rmapExecutablePath,fastaFilePath,genomeIndexBasePath,outputFilePath,tmpDirPath,bowtieBuildExecutablePath, optionalParameters));
			Process rmap = Runtime.getRuntime().exec(String.format("%s -i %s -index %s -o %s -tmpdir %s -btbuild %s %s",rmapExecutablePath,fastaFilePath,genomeIndexBasePath,outputFilePath,tmpDirPath,bowtieBuildExecutablePath, optionalParameters));
			
			//NEW STREAM HANDLING
			//handle output
			StreamConsumer errorOut = new StreamConsumer(rmap.getErrorStream(),"err",verbose);
			StreamConsumer stdOut = new StreamConsumer(rmap.getInputStream(),"stdout",verbose);
			errorOut.start();
			stdOut.start();
			int returnValue = rmap.waitFor();
			
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
			rmap.getErrorStream().close();
			rmap.getOutputStream().close();
			rmap.getInputStream().close();
			rmap.destroy();
			
						
			if(returnValue != 0) {
				throw new Exception(String.format("Something went wrong,rmap return value: %s.",returnValue));
			}	
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	public void alignReadsWithBowtie(String bowtieExecutablePath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, String seedLength, String seedMissmatches, String maxMismatches,
			int maxHits, String threads, boolean verbose) {
		
		try {
			
			String optionalParameters = String.format("-mm %s -l %s -n %s",maxMismatches,seedLength,seedMissmatches);
			if(maxHits != -1) optionalParameters += " -maxHits " + maxHits;
			
			
			
			
			
			if(verbose) System.out.println(String.format("%s -f -n %s -l %s -a -m %s --suppress 5,6,7 -e %s --tryhard --nomaqround -p %s %s %s %s",bowtieExecutablePath,seedMissmatches,seedLength,maxHits,(this.baseQuality * Integer.valueOf(maxMismatches)) + 1,threads,genomeIndexBasePath,fastaFilePath,outputFilePath)) ;
			Process bowtie = Runtime.getRuntime().exec(String.format("%s -f -n %s -l %s -a -m %s --suppress 5,6,7 -e %s --tryhard --nomaqround -p %s %s %s %s",bowtieExecutablePath,seedMissmatches,seedLength,maxHits,(this.baseQuality * Integer.valueOf(maxMismatches)) + 1,threads,genomeIndexBasePath,fastaFilePath,outputFilePath));
			
			
			
			//NEW STREAM HANDLING
			//TODO check changes to old handling
			//handle output
			StreamConsumer errorOut = new StreamConsumer(bowtie.getErrorStream(),"err",verbose);
			StreamConsumer stdOut = new StreamConsumer(bowtie.getInputStream(),"stdout",verbose);
			errorOut.start();
			stdOut.start();
			int returnValue = bowtie.waitFor();
			
			//notify stream threads that the process terminated
			errorOut.processTerminated();
			stdOut.processTerminated();
			
			//close everything
			bowtie.getErrorStream().close();
			bowtie.getOutputStream().close();
			bowtie.getInputStream().close();
			bowtie.destroy();
			
			if(returnValue != 0) {
			throw new Exception(String.format("Something went wrong,rmap return value: %s.",returnValue));
			}	
		}
		
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
}
