package alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


import tools.RmapProcessor;
import tools.SequenceConverter;
import tools.UnixSort;
import main.Pair;

public class AlignmentCoordinator {

	
	private String alignerName;
	private String alignerBinPath;
	private String alignerIndexerBinPath;
	private String workingDir;
	private String referencesDir;
	
	private int maxAllowedMismatches;
	private int seedSize;
	private int seedMismatches;
	private int maxHits;
	private int[] splitSeedSizes;
	private int[] splitSeedMismatches;
	private int maxIntronLength;
	private int maxIntronCount;
	private int maxContextSize;
	private int contextBufferSize;
	
	private int maxGapSize;
	private int minGapSize;
	private int maxDelSize;
	
	private int threads;
	
	private boolean pairedEnd;
	private boolean skipPass1;
	private boolean skipBackwardAlignment;
	private boolean autosetAlignmentOptions;
	private boolean verbose;
	
	
	private ReadAligner readAligner;
	
	private final int minExonLength = 20;
	
	
	
	public AlignmentCoordinator(String alignerName, String alignerBinPath, String alignerIndexerPath, String workingDir, String referencesDir,
						   int maxAllowedMismatches, int seedSize, int seedMismatches, int maxHits, int[] splitSeedSizes, int[] splitSeedMismatches, int maxContextSize, int contextBufferSize, int maxGapSize, int minGapSize, int maxDelSize, int threads,
						   boolean pairedEnd, boolean skipPass1, boolean autosetAlignmentOptions, boolean verbose) {
		
		
		
		this.alignerName = alignerName;
		this.alignerBinPath = alignerBinPath;
		this.alignerIndexerBinPath = alignerIndexerPath;
		this.workingDir = workingDir;
		this.referencesDir = referencesDir;
		
		this.maxAllowedMismatches = maxAllowedMismatches;
		this.seedSize = seedSize;
		this.seedMismatches = seedMismatches;
		this.maxHits = maxHits;
		this.splitSeedSizes = splitSeedSizes;
		this.splitSeedMismatches = splitSeedMismatches;
		this.maxContextSize = maxContextSize;
		this.contextBufferSize = contextBufferSize;
		this.maxGapSize = maxGapSize;
		this.minGapSize = minGapSize;
		this.maxDelSize = maxDelSize;
		this.threads = threads;
		
		this.pairedEnd = pairedEnd;
		this.skipPass1 = skipPass1;
		this.autosetAlignmentOptions = autosetAlignmentOptions;
		this.verbose = verbose;
		
		
		this.readAligner = initAligner();
		
	}
	
	
	/**
	 * initializes the aligner with specific settings
	 */
	private ReadAligner initAligner() {
		
		/**
		 * 
		 * BOWTIE 1 SETTINGS
		 * 
		 */
		if(this.alignerName.equals("bowtie") || this.alignerName.equals("bowtie1")) {
			if(this.splitSeedSizes == null || this.splitSeedMismatches == null) {
				this.splitSeedSizes = new int[]{10,15};
				this.splitSeedMismatches = new int[]{0,1};
			}
			if(this.seedSize == -1) {
				this.seedSize = 30;
			}
			if(this.seedMismatches == -1)
				this.seedMismatches = 1;

			if(this.maxHits == -1)
				this.maxHits = 20;
			
			this.skipBackwardAlignment = false;
			
			return new ReadAlignerBowtie1();
		}
		
		
		/**
		 * 
		 * BOWTIE 2 SETTINGS
		 * 
		 */
		else if(this.alignerName.equals("bowtie2")) {
			if(this.splitSeedSizes == null || this.splitSeedMismatches == null) {
				this.splitSeedSizes = new int[]{10, 15};
				this.splitSeedMismatches = new int[]{0,1};
			}
			if(this.seedSize == -1) {
				this.seedSize = 20;
			}
			if(this.seedMismatches == -1)
				this.seedMismatches = 0;

			if(this.maxHits == -1)
				this.maxHits = 3;
			
			this.skipBackwardAlignment = true;
			
			return new ReadAlignerBowtie2();
		}
		
		
		/**
		 * 
		 * BWA SETTINGS
		 * 
		 */
		else if(alignerName.equals("bwa") || alignerName.equals("BWA")) {
			if(this.splitSeedSizes == null || this.splitSeedMismatches == null) {
				this.splitSeedSizes = new int[]{15};
				this.splitSeedMismatches = new int[]{0};
			}
			if(this.seedSize == -1) {
				this.seedSize = 20;
			}
			if(this.seedMismatches == -1)
				this.seedMismatches = 0;

			if(this.maxHits == -1)
				this.maxHits = 10;
			
			this.skipBackwardAlignment = true;
			
			return new ReadAlignerBWA(this.workingDir + "/bwa_tmp",this.referencesDir);
		}
		
		
		/**
		 *
		 * ADD CUSTOM ALIGNER SETTINGS HERE
		 * 
		 */
		else {
			
		}
		
		return null;
	}
	
	
	
	
	public void start(String inputFilePath, String outputFilePath, String multiSplitSeqsOutputPath, String bwindexPath, boolean skipSplitDetection, boolean skipMultiSplitDetection) {
		try {
			
			/**
			 * re-naming header for bwa in paired end mode (otherwise bwa removes the extensions /1 and /2 from the read names)
			 */
			String tmpInputFilePath = null;
			
			//here we already know that the input file is not empty and that the first line is a read header
			boolean hasPairedEndHeader = false;
			boolean hasSpacesInHeader = false;
			boolean spaceInHeaderRenamed = false;
			boolean hasRenamedPairedEndHeader = false;
			BufferedReader tmpReader = new BufferedReader(new FileReader(new File(inputFilePath)));
			String tmpHeader = tmpReader.readLine();
			tmpReader.close();
			if(tmpHeader.charAt(tmpHeader.length() - 2) == '/')
				hasPairedEndHeader = true;
			if(tmpHeader.contains(" "))
				hasSpacesInHeader = true;
			
			
			if(((this.alignerName.equals("bwa") || this.alignerName.equals("BWA")) && (this.pairedEnd || hasPairedEndHeader || hasSpacesInHeader)) || ((this.alignerName.equals("bowtie") || this.alignerName.equals("bowtie1")) && hasSpacesInHeader)) {
				String newInputFilePath = this.workingDir + "/all_reads.fa";
				//SequenceConverter.renameReadHeaderForBwa(inputFilePath, newInputFilePath);
				String tmpDirPath = this.workingDir + "/renaming_tmp";
				new File(tmpDirPath).mkdir();
				SequenceConverter.renameReadHeaderForBwa(inputFilePath,tmpDirPath, newInputFilePath, (this.pairedEnd || hasPairedEndHeader), this.threads);
				tmpInputFilePath = inputFilePath;
				inputFilePath = newInputFilePath;
				spaceInHeaderRenamed = hasSpacesInHeader;
				hasRenamedPairedEndHeader = (this.pairedEnd || hasPairedEndHeader);
			}
			
			
			Pair<Integer,Integer> minAndMaxReadLength = getMinAndMaxReadLength(inputFilePath);
			int minReadLength = minAndMaxReadLength.getFirst();
			int maxReadLength = minAndMaxReadLength.getSecond();
			
			
			//automatically sets seed size, splitseedsizes splitseed mismatches dependend on the maximum read length
			if(this.autosetAlignmentOptions) {
				if(this.alignerName.equals("bowtie") || this.alignerName.equals("bowtie1")) {
					if(maxReadLength >= 100) {
						this.splitSeedSizes = new int[]{15,20};
						this.splitSeedMismatches = new int[]{0,1};
						this.seedSize = 30;
						this.seedMismatches = 1;
						this.maxHits = 20;
					}
					
					else if(maxReadLength >= 76) {
						this.splitSeedSizes = new int[]{10,15};
						this.splitSeedMismatches = new int[]{0,1};
						this.seedSize = 30;
						this.seedMismatches = 1;
						this.maxHits = 20;
					}
					
					else {
						this.splitSeedSizes = new int[]{10};
						this.splitSeedMismatches = new int[]{0};
						this.seedSize = 18;
						this.seedMismatches = 1;
						this.maxHits = 50;
					}
				}
				
				else if(this.alignerName.equals("bowtie2")) {
					if(maxReadLength >= 100) {
						this.splitSeedSizes = new int[]{15, 20};
						this.splitSeedMismatches = new int[]{0,1};
						this.seedSize = 25;
						this.seedMismatches = 0;
						this.maxHits = 3;
					}
					
					else if(maxReadLength >= 76) {
						this.splitSeedSizes = new int[]{10, 15};
						this.splitSeedMismatches = new int[]{0,1};
						this.seedSize = 20;
						this.seedMismatches = 0;
						this.maxHits = 3;
					}
					
					else {
						this.splitSeedSizes = new int[]{10};
						this.splitSeedMismatches = new int[]{0};
						this.seedSize = 18;
						this.seedMismatches = 0;
						this.maxHits = 10;
					}
					
					
					
				}
				
				else if(alignerName.equals("bwa") || alignerName.equals("BWA")) {
					if(maxReadLength >= 100) {
						this.splitSeedSizes = new int[]{20};
						this.splitSeedMismatches = new int[]{0};
						this.seedSize = 25;
						this.seedMismatches = 0;
						this.maxHits = 10;
					}
					
					else if(maxReadLength >= 76) {
						this.splitSeedSizes = new int[]{15};
						this.splitSeedMismatches = new int[]{0};
						this.seedSize = 20;
						this.seedMismatches = 0;
						this.maxHits = 10;
					}
					
					else {
						this.splitSeedSizes = new int[]{14};
						this.splitSeedMismatches = new int[]{0};
						this.seedSize = 18;
						this.seedMismatches = 0;
						this.maxHits = 20;
					}
					
					
				}
			}
			
			
			/**
			 * determine forward alignments
			 */
			Date date = new Date();
			ArrayList<String> filePaths = new ArrayList<String>();
			filePaths.add(outputFilePath);
			String splitSeedSizesAsString = Integer.valueOf(this.splitSeedSizes[0]).toString();
			String splitSeedMismatchesAsString = Integer.valueOf(this.splitSeedMismatches[0]).toString();
			for(int i = 1; i < this.splitSeedSizes.length;i++) {
				splitSeedSizesAsString += "," + Integer.valueOf(this.splitSeedSizes[i]).toString();
				splitSeedMismatchesAsString += "," + Integer.valueOf(this.splitSeedMismatches[i]).toString();
			}
			
			System.out.println(String.format("[%s]\tDetermining forward alignments (seed: %s, seed-mm: %s, splitseedsizes: %s, splitseed-mm: %s, skipsplit: %s, skipmultisplit: %s, max_readlen: %s, autoset_options: %s)", date.toLocaleString(),this.seedSize,this.seedMismatches,splitSeedSizesAsString,splitSeedMismatchesAsString,skipSplitDetection,skipMultiSplitDetection,maxReadLength,this.autosetAlignmentOptions));
			String unalignedReadsFirstPassOutputPath = this.workingDir + "/fwd_unaligned_reads_pass_1.fa";
			filePaths.add(unalignedReadsFirstPassOutputPath);
			
			
			String unalignedReadsSecondPassOutputPath = this.workingDir + "/fwd_unaligned_reads_pass_2.fa";;
			filePaths.add(unalignedReadsSecondPassOutputPath);
			
			File tmpFile;
			for(String tmpPath : filePaths) {
				tmpFile = new File(tmpPath);
				if(tmpFile.exists()) {
					tmpFile.delete();				
				}
			}
			filePaths.clear();
			
			String splitCandidatesOutputDir = this.workingDir + "/split_candidates";
			filePaths.add(splitCandidatesOutputDir);
			String multiSplitCandidatesOutputDir = this.workingDir + "/multi_split_candidates";
			filePaths.add(multiSplitCandidatesOutputDir);
			
			for(String tmpPath : filePaths) {
				tmpFile = new File(tmpPath);
				if(!tmpFile.isDirectory())
					tmpFile.mkdir();
				else {
					File[] tmpFiles = tmpFile.listFiles();
					for(File tmp : tmpFiles)
						tmp.delete();
				}
			}
			filePaths.clear();
			
			/**
			 *  starting the initial alignment phase
			 */
			this.readAligner.alignReadsInitialPhase(alignerBinPath, inputFilePath, bwindexPath, outputFilePath,unalignedReadsSecondPassOutputPath, splitCandidatesOutputDir,multiSplitCandidatesOutputDir, seedSize, seedMismatches, splitSeedSizes, maxAllowedMismatches, maxHits, maxReadLength, minReadLength, threads,false,skipSplitDetection,skipMultiSplitDetection,true);
			
			
			if(!skipBackwardAlignment && new File(unalignedReadsSecondPassOutputPath).exists()) {
				/**
				 *  get the reverse complement of unaligned reads reads
				 */
				date = new Date();
				System.out.println(String.format("[%s]\tPreparing read sequences for the backward alignmenpt step", date.toLocaleString()));
				ReadReverser readReverser = new ReadReverser();
				File f = new File(unalignedReadsSecondPassOutputPath);
				String rcReadFilePath = this.workingDir + "/" + f.getName().substring(0,f.getName().lastIndexOf('.')) + "_rc.fa";
				readReverser.getReverseComplement(unalignedReadsSecondPassOutputPath, rcReadFilePath, "fasta", true);
				
				/**
				 * determine backward alignments
				 */
				date = new Date();
				new File(unalignedReadsSecondPassOutputPath).delete();
				System.out.println(String.format("[%s]\tDetermining backward alignments", date.toLocaleString()));
				this.readAligner.alignReadsInitialPhase(this.alignerBinPath, rcReadFilePath, bwindexPath, outputFilePath,null, splitCandidatesOutputDir,multiSplitCandidatesOutputDir, seedSize, seedMismatches, splitSeedSizes, maxAllowedMismatches, maxHits,maxReadLength, minReadLength, threads,true,skipSplitDetection,skipMultiSplitDetection,true);
				
			}
			
			else {
				date = new Date();
				System.out.println(String.format("[%s]\tSkipping backward alignment", date.toLocaleString()));
			}
			
			
			if(!skipSplitDetection) {
				UnixSort unixSorter;
				String sortedSplitCandidatesOutputDir;
				ExecutorService sortExecutor = Executors.newFixedThreadPool(this.threads);
				ArrayList<Future> futures = new ArrayList<Future>();
				int threadIndex;
				File[] candidateFiles;
				
				if( new File(multiSplitCandidatesOutputDir).listFiles().length > 0) {
					/**
					 * Sort multi split candidates by read id and start positions
					 */
					date = new Date();
					System.out.println(String.format("[%s]\tSorting multi split candidate alignments by id and start positions", date.toLocaleString()));
					
					sortedSplitCandidatesOutputDir = this.workingDir + "/sorted_multi_split_candidates";
					tmpFile = new File(sortedSplitCandidatesOutputDir);
					if(!tmpFile.isDirectory())
						tmpFile.mkdir();
					else {
						for(File file : tmpFile.listFiles())
							file.delete();
					}
					
					threadIndex = 0; 
					candidateFiles = new File(multiSplitCandidatesOutputDir).listFiles();
					int[] columns = {1,4};
					//empty option means default (dictionary) order.
					String[] orderingOptions = {"","n"};
					for(File rmapFile : candidateFiles) {
						
						futures.add(sortExecutor.submit(new UnixSort(rmapFile.getAbsolutePath(),sortedSplitCandidatesOutputDir + "/" + rmapFile.getName(),multiSplitCandidatesOutputDir + "/tmp/" + threadIndex,"\t",columns,orderingOptions,500,true,true)));
						threadIndex++;
					}
					sortExecutor.shutdown();
		
					for(Future future : futures) {
						future.get();
					}
					futures.clear();
					for(File rmapFile : candidateFiles) {
						rmapFile.delete();
					}
					
					
					/**
					 * add multi split candidates to the rmap files
					 */
					date = new Date();
					System.out.println(String.format("[%s]\tProcessing multi split candidates", date.toLocaleString()));
					SplitCandidateExtractor.processMultiSplitCandidates(sortedSplitCandidatesOutputDir, this.minExonLength, splitSeedSizes[0],maxAllowedMismatches, this.maxHits, splitCandidatesOutputDir,outputFilePath, multiSplitSeqsOutputPath,pairedEnd,spaceInHeaderRenamed,hasRenamedPairedEndHeader);
				}
				
				
				
				
				
				/**
				 * Sort split candidates by start positions
				 */
				Thread.sleep(500);
				date = new Date();
				System.out.println(String.format("[%s]\tSorting split candidate alignments by start positions", date.toLocaleString()));
				sortedSplitCandidatesOutputDir = this.workingDir + "/sorted_split_candidates";
				tmpFile = new File(sortedSplitCandidatesOutputDir);
				if(!tmpFile.isDirectory())
					tmpFile.mkdir();
				else {
					for(File file : tmpFile.listFiles())
						file.delete();
				}
				
				sortExecutor = Executors.newFixedThreadPool(this.threads);
				futures.clear();
				threadIndex = 0; 
				candidateFiles = new File(splitCandidatesOutputDir).listFiles();
				for(File rmapFile : candidateFiles) {
					
					futures.add(sortExecutor.submit(new UnixSort(rmapFile.getAbsolutePath(),sortedSplitCandidatesOutputDir + "/" + rmapFile.getName(),splitCandidatesOutputDir + "/tmp/" + threadIndex,"\t",3,500,true,true,true)));
					threadIndex++;
				}
				sortExecutor.shutdown();
	
				for(Future future : futures) {
					future.get();
				}
				futures.clear();
				for(File rmapFile : candidateFiles) {
					rmapFile.delete();
				}
				
				
				/**
				 * Determining split candidates
				 */
				date = new Date();
				System.out.print(String.format("[%s]\tDetermining split candidates (mingapsize: %s, maxgapsize: %s, maxdelsize: %s)...", date.toLocaleString(),this.minGapSize,this.maxGapSize,this.maxDelSize));
				String splitDetectionDir = this.workingDir + "/split_extraction";
				tmpFile = new File(splitDetectionDir);
				if(!tmpFile.isDirectory())
					tmpFile.mkdir();
				else {
					File[] tmpFiles = tmpFile.listFiles();
					for(File tmp : tmpFiles)
						tmp.delete();
				}
				SplitCandidateExtractor sce = new SplitCandidateExtractor(readAligner,alignerBinPath,alignerIndexerBinPath, referencesDir, splitDetectionDir, maxAllowedMismatches, splitSeedSizes, splitSeedMismatches, maxContextSize, contextBufferSize,maxGapSize, minGapSize, maxDelSize, maxHits,maxReadLength, minReadLength, threads);
				sce.processAlignments(sortedSplitCandidatesOutputDir, outputFilePath,this.verbose);
				System.out.println();
			}
			
			
			/**
			 * undo header re-naming
			 */
			
			if(((this.alignerName.equals("bwa") || this.alignerName.equals("BWA")) && (this.pairedEnd || hasPairedEndHeader || hasSpacesInHeader)) || ((this.alignerName.equals("bowtie") || this.alignerName.equals("bowtie1")) && hasSpacesInHeader)) {
				inputFilePath = tmpInputFilePath;
				String tmpOutputFilePath = outputFilePath + ".tmp";
				//RmapProcessor.undoHeaderRenamingFromBwa(outputFilePath, tmpOutputFilePath);
				String tmpDirPath = this.workingDir + "/renaming_tmp";
				new File(tmpDirPath).mkdir();
				RmapProcessor.undoHeaderRenamingFromBwa(outputFilePath,tmpDirPath, tmpOutputFilePath,this.threads, this.pairedEnd, hasPairedEndHeader);
				new File(outputFilePath).delete();
				new File(tmpOutputFilePath).renameTo(new File(outputFilePath));
			}
			
			removeTmpFiles(this.workingDir);
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private void removeTmpFiles(String dirPath) {
		File[] files = new File(dirPath).listFiles();
		for(File f : files) {
			if(f.isDirectory())
				removeTmpFiles(f.getAbsolutePath());
			
			f.delete();
		}
		
	}
	
	private Pair<Integer,Integer> getMinAndMaxReadLength(String fastaFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			String currentLine;
			int minReadLength = Integer.MAX_VALUE;
			int maxReadLength = Integer.MIN_VALUE;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					currentLine = br.readLine();
					if(currentLine.length() > maxReadLength)
						maxReadLength = currentLine.length();
					
					if(currentLine.length() < minReadLength)
						minReadLength = currentLine.length();
				}
			}
			br.close();
			
			Pair<Integer,Integer> lengths = new Pair<Integer,Integer>();
			lengths.setFirst(minReadLength);
			lengths.setSecond(maxReadLength);
			return lengths;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}


	public int getSeedSize() {
		return seedSize;
	}


	public void setSeedSize(int seedSize) {
		this.seedSize = seedSize;
	}


	public int getSeedMismatches() {
		return seedMismatches;
	}


	public void setSeedMismatches(int seedMismatches) {
		this.seedMismatches = seedMismatches;
	}


	public int getMaxHits() {
		return maxHits;
	}


	public void setMaxHits(int maxHits) {
		this.maxHits = maxHits;
	}


	public int[] getSplitSeedSizes() {
		return splitSeedSizes;
	}


	public void setSplitSeedSizes(int[] splitSeedSizes) {
		this.splitSeedSizes = splitSeedSizes;
	}


	public int[] getSplitSeedMismatches() {
		return splitSeedMismatches;
	}


	public void setSplitSeedMismatches(int[] splitSeedMismatches) {
		this.splitSeedMismatches = splitSeedMismatches;
	}
	
	
	
	
}
