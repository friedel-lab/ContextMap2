package alignment;

import java.util.ArrayList;

interface ReadAligner {

	public void alignReadsInitialPhase(String alignerBinPath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, String unalignedReadsOutputFilePath, 
			String splitCandidatesOutputDir, String multiSplitCandidatesOutputDir, int seedLength, int seedMissmatches, int[] splitSeedSizes,
			int maxMismatches, int maxHits, int maxReadLength, int minReadLength, int threads, boolean filterFullAlignments, 
			boolean skipSplitDetection, boolean skipMultiSplitDetection, boolean verbose);
	
	
	public void alignReadsSlidingWindow(String alignerBinPath, String fastaFilePath, String genomeIndexBasePath,
			String outputFilePath, int seedLength, int seedMissmatches, int maxMismatches, int maxHits, int maxReadLength, 
			int minReadLength, int threads, boolean filterFullAlignments, boolean skipSplitDetection, boolean alignToFwdStrand, 
			boolean alignToRevStrand, boolean verbose);
	
	
	public void buildIndex(String indexerBinPath, String fastaFilePath, String outputFileBasePath);
	
	/** to implement this: Simply call: return TestSamHeader.checkIndex(alignerBinPath, genomeIndexBasePath, referenceDir, tmpFolder, this); */
//	public ArrayList<String> checkIndex(String alignerBinPath, String genomeIndexBasePath, String referenceDir, String tmpFolder);
	
	/** command for getting SAM header because of index **/
//	public String getSamHeaderCommand(String alignerBinPath, String genomeIndexBasePath, String path2dummyFastaFile);
	
	
}
