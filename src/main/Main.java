package main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import assembly.AssemblyProcessor;
import context.ContextExtractor;
import context.MappingProcessor;
import tools.BamConverter;
import tools.FileHandler;
import tools.FileSorter;
import tools.GenomeIndexer;
import tools.ReadAligner;
import tools.ReadSampler;
import tools.RmapProcessor;
import tools.SamProcessor;
import tools.SequenceConverter;

public class Main {

	private static final String contextmapVersion = "v2.7.9";
	
	public static void main(String args[]) {
		
		if(args.length == 0 || args[0].equals("--help") || args[0].equals("-h")) {
			System.out.println();
			System.out.println(String.format("Usage: java -jar ContextMap_%s.jar <tool name> <tool arguments> [tool options]*",contextmapVersion));
			System.out.println();
			System.out.println("Available tools:");
			System.out.println("\t 1. 'mapper' - The ContextMap mapping tool");
			System.out.println("\t 2. 'indexer' - Prepares large sets of genomes for being indexed with Bowtie");
			System.out.println("\t 3. 'inspector' - Determines read counts, confidence values, genome coverages and the square root of the Jensen-Shannon divergence for species contained in a sam file");
			System.out.println("\t 4. 'separator' - Splits up a sam file containing a mixture of mappings to a set of different species");
			//System.out.println("\t 5. 'error_estimator' - Determines the mismatch rates for a set of species contained in a given sam file");
			//System.out.println("\t 6. 'confidence_estimator' - Determines confidence values for a set of species contained in a given sam file");
			//System.out.println("\t 7. 'distance_estimator' - Determines distance values for a set of species contained in a given sam file");
			System.out.println();
			System.out.println("Start a tool without any argument to get a specific help message and a list of available options for this tool");
			System.out.println();
			System.exit(1);
		}
		
		if(args[0].equals("mapper")) {
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar mapper <arguments> [options]*",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-reads\t\t<A comma separated list of file paths to reads in fasta/fastq/zip/gz format. A single path for single-end mapping and two paths (#1 mates and #2 mates) for paired-end mapping>");
				System.out.println("-aligner_name\t<The following alignment tools are supported: \"bwa\", \"bowtie1\" or \"bowtie2\">");
				System.out.println("-aligner_bin\t<The absolute path to the executable of the chosen aligner tool>");
				System.out.println("-indexer_bin\t<The absolute path to the executable of the aligner's indexing tool (not needed for BWA)>");
				System.out.println("-indices\t<A comma separated list of paths to basenames of indices, which can be used by the chosen aligner>");
				System.out.println("-genome\t\t<The path to a directory with genome sequences in fasta format (for each chromosome a separate file)>");
				System.out.println("-o\t\t<The path to the output directory>");
								
				System.out.println();
				System.out.println();
				System.out.println("Options:");
				System.out.println();
				
				
// DEPRECATED				//System.out.println("-sam\t\t\t<The path to an existing mapping in sam format. It is required that the sam file is sorted by read id and that the NM field is set.>");
				System.out.println("--mining\t\t<Set this to mine for infections or contaminations. Changes of ContextMap's mapping behaviour due to this option are described in the manual.>");
				System.out.println("-skipsplit\t\t<A comma separated list of booleans, each element refers to a given aligner index (same ordering). 'true' for no split detection, 'false' otherwise (REQ. in mining mode).>");
				System.out.println("-skipmultisplit\t\t<A comma separated list of booleans, each element refers to a given aligner index (same ordering). 'true' for no multi split detection, 'false' otherwise (REQ. in mining mode).>");
				System.out.println("-speciesindex\t\t<The path to a directory containing index files created with the 'indexer' tool (REQ. in mining mode)>");
				System.out.println("-aligner_tmp\t\t<The path to a directory for temporary alignment files>");
				System.out.println("-seed\t\t\t<The seed length for the alignment> (default: Bwt1: 30, BWA/Bwt2: 20)");
				System.out.println("-seedmismatches\t\t<Allowed mismatches in the seed region> (default: Bwt1: 1, BWA/Bwt2: 0)");
				System.out.println("-mismatches\t\t<Allowed mismatches in the whole read> (default: 4)");
				System.out.println("-maxindelsize\t\t<The maximum allowed size of insertions or deletions> (default: 10)");				
				System.out.println("-mmdiff\t\t\t<The maximum allowed mismatch difference between the best and second best alignment of the same read> (default: 0)");
				System.out.println("-maxhits\t\t<The maximum number of candidate alignments per read. Reads with more hits are skipped (bwa/bwt1) or the already found hits are reported (bwt2) > (default for bwa/bwt1:10, bwt2: 3)");
				System.out.println("-minsize\t\t<The minimum number of reads a genomic region has to contain for being regarded as a local context. (default:10)");
				System.out.println("-annotation\t\t<The path to an annotation file in our own format>");
				System.out.println("-gtf\t\t\t<The path to an annotation file in gtf format. This option is mutually exclusive with the -annotation option>");
				System.out.println("-t\t\t\t<The number of threads used for the run> (default: 1)");
				System.out.println("--autosetalioptions\t<Automatically sets the alignment options dependent on the maximum read length>");
				System.out.println("--noclipping\t\t<Disables the calculation of clipped alignments> (default: not set)");
				System.out.println("--noncanonicaljunctions\t<Enables the prediction of non-canonical splice sites> (default: not set)");	
				System.out.println("--strandspecific\t<Enables strand specific mapping (default: off)>");
				//System.out.println("--pairedend\t\t<Enables mapping of paired-end reads. Nomenclature for mates from the same fragment: base_name/1 and base_name/2, respectively. (default: off)>");
				System.out.println("--polyA\t\t\t<Enables the search for polyA-tails. This option is mutually exclusive with the --noclipping option. (default: off)>");
				System.out.println("--sequenceDB\t\t<Writes a readId -> sequence mapping to disk. Recommended for very large data sets. (default: off)>");
				System.out.println("--verbose\t\t<verbose mode> (default: not set)");
				System.out.println("--keeptmp\t\t<does not delete some temporary files> (default: not set)");
				System.out.println();
				System.exit(1);
			}
			
			Date date = new Date();
			System.out.println(String.format("[%s]\tStarting ContextMap (%s) run.",date.toLocaleString(),contextmapVersion));
			
			/*
			 * global stuff
			 */
			ArrayList<String> filePaths = new ArrayList<String>();
			String readFilePath = null;
			String samFilePath = null;
			String outputDirPath = null;
			String readFormat;
			int readLength = -1;
			boolean verbose = false;
			boolean developer = false;
			boolean keepTmp = false;
			boolean clipping = true;
			boolean polyA = false;
			boolean strandedPolyA = false;
			
			/*
			 * alignment stuff
			 */
			String alignerName = null;
			String alignerBinPath = null;
			String alignerIndexerPath = null;
			String alignerTmpDirPath = null;
			ArrayList<String> genomeIndexBasePaths = new ArrayList<String>();
			
			
			//optional fields (all unset default values are set in the SlidingContext class)
			int seedSize = -1;
			int seedMismatches = -1;
			int maxMismatches = 4;
			int[] splitSeedSizes = null;
			int[] splitSeedMismatches = null;
			int maxIntronLength = 200000;
			int maxIntronCount = 10;
			int maxContextSize = maxIntronCount * maxIntronLength;
			boolean autosetAlignmentOptions = false;
			
			
			int minGapSize = 50;
			int maxGapSize = 300000;
			int maxIndelSize = 10;
			
			int maxHits = -1;
			int threads = 1;
			boolean skipPass1 = true;
			ArrayList<Boolean> skipSplitDetection = new ArrayList<Boolean>();
			ArrayList<Boolean> skipMultiSplitDetection = new ArrayList<Boolean>();
			boolean skipRealignment = false;
			
			/*
			 * context processing stuff
			 */
			String referencesDir = null;
			String annotationFilePath = null;
			String gtfFilePath = null;
			String indexDirPath = null;
			int maxMismatchDifference = 0;
			int minDistanceBetweenContext = 10000;
			int minNumberOfReadsInContext = 10;
			boolean pairedEnd = false;
			boolean strandSpecific = false;
			boolean preferExtensionsWithKnownSpliceSignals = true;
			boolean skipDenovoJunctions = false;
			boolean skipNonCanonicalJunctions = true;
			boolean updateQueue = false;
			boolean printMultiMappings = false;
			boolean printSecondBestChr = false;
			boolean lowCoverage = true;
			int updateInterval = 3;
			
			int minPolyALength = 6;
			int minPolyAReadCount = 3;
			int maxConsideredClippingLength = 30;
			double upperPolyACutoff = 1.0;
			double lowerPolyACutoff = 0.7;
			
			
			boolean writeSequenceDB = false;
			int databaseSortBatchSize = 3000000;
				
			
			
			for(int i = 1; i < args.length; i++) {
				
				//global stuff
				if(args[i].equals("-reads")) {
					readFilePath = args[++i];
					filePaths.addAll(Arrays.asList(readFilePath.split(",")));
					
					if(filePaths.size() == 1)
						pairedEnd = false;
					
					else if(filePaths.size() == 2)
						pairedEnd = true;
					
					else if(filePaths.size() > 2) {
						System.err.println(String.format("[%s]\tOption -reads accepts a maximum of two file paths only. Please check the help by calling ContextMap without any parameters.",date.toLocaleString()));
						System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
						System.exit(1);
					}
					
					continue;
				}
				
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				
				if(args[i].equals("-o")) {
					outputDirPath = args[++i];
					continue;
				}
				if(args[i].equals("-t")) {
					threads = Integer.valueOf(args[++i]);
					continue;
				}
				
				if(args[i].equals("--verbose")) {
					verbose = true;
					continue;
				}
				
				if(args[i].equals("--developer")) {
					developer = true;
					continue;
				}
				
				if(args[i].equals("--keeptmp")) {
					keepTmp = true;
					continue;
				}
				
				//alignment stuff
				if(args[i].equals("-aligner_name")) {
					alignerName = args[++i];
					continue;
				}
				
				if(args[i].equals("-aligner_bin")) {
					alignerBinPath = args[++i];
					if(alignerName.equals("bwa") || alignerName.equals("BWA"))
						alignerIndexerPath = alignerBinPath;
					
					continue;
				}
				
				if(args[i].equals("-indexer_bin")) {
					alignerIndexerPath = args[++i];
					continue;
				}
				
				if(args[i].equals("-aligner_tmp")) {
					alignerTmpDirPath = args[++i];
					continue;
				}
				
				if(args[i].equals("-indices")) {
					String[] tmpGenomeIndexBasePaths = args[++i].split(",");
					genomeIndexBasePaths.addAll(Arrays.asList(tmpGenomeIndexBasePaths));
					continue;
				}				
				
				if(args[i].equals("-seed")) {
					seedSize = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-mismatches")) {
					maxMismatches = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-seedmismatches")) {
					seedMismatches = Integer.valueOf(args[++i]);
					continue;
				}
				
				if(args[i].equals("-splitseedsizes")) {
					String[] tmpSizes = args[++i].split(",");
					splitSeedSizes = new int[tmpSizes.length];
					for(int j = 0; j < tmpSizes.length; j++) {
						splitSeedSizes[j] = Integer.valueOf(tmpSizes[j]);
					}
					continue;
				}
				
				if(args[i].equals("-splitseedmismatches")) {
					String[] tmpSizes = args[++i].split(",");
					splitSeedMismatches = new int[tmpSizes.length];
					for(int j = 0; j < tmpSizes.length; j++) {
						splitSeedMismatches[j] = Integer.valueOf(tmpSizes[j]);
					}
					continue;
				}
				
				if(args[i].equals("--autosetalioptions")) {
					autosetAlignmentOptions = true;
					continue;
				}
				
				
				if(args[i].equals("-mmdiff")) {
					maxMismatchDifference = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-maxhits")) {
					maxHits = Integer.valueOf(args[++i]);
					continue;
				}
				
				if(args[i].equals("--pass1")) {
					skipPass1 = false;
					continue;
				}
				
				if(args[i].equals("-skipsplit")) {
					String[] tmpSkipSplitDetection = args[++i].split(",");
					for(int j = 0; j < tmpSkipSplitDetection.length; j++) {
						skipSplitDetection.add(Boolean.valueOf(tmpSkipSplitDetection[j]));
					}
					continue;
				}
				
				if(args[i].equals("-skipmultisplit")) {
					String[] tmpSkipSplitDetection = args[++i].split(",");
					for(int j = 0; j < tmpSkipSplitDetection.length; j++) {
						skipMultiSplitDetection.add(Boolean.valueOf(tmpSkipSplitDetection[j]));
					}
					continue;
				}
				
				if(args[i].equals("--skipdenovojunctions")) {
					skipDenovoJunctions = true;
					continue;
				}
				
				if(args[i].equals("--noncanonicaljunctions")) {
					skipNonCanonicalJunctions = false;
					continue;
				}
				
				if(args[i].equals("-skiprealignment")) {
					skipRealignment = true;
					continue;
				}
				
				//context processing stuff
				if(args[i].equals("-genome")) {
					referencesDir = args[++i];
					continue;
				}
				if(args[i].equals("-annotation")) {
					annotationFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-gtf")) {
					gtfFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-idx") || args[i].equals("-speciesindex")) {
					indexDirPath = args[++i];
					continue;
				}
				if(args[i].equals("-mindist")) {
					minDistanceBetweenContext = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-minsize")) {
					minNumberOfReadsInContext = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-maxlength")) {
					maxContextSize = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-maxgapsize")) {
					maxGapSize = Integer.valueOf(args[++i]);
					maxIntronLength = maxGapSize;
					continue;
				}
				if(args[i].equals("-mingapsize")) {
					minGapSize = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-maxindelsize")) {
					maxIndelSize = Integer.valueOf(args[++i]);
					continue;
				}
				
				if(args[i].equals("-minPolyALength")) {
					minPolyALength = Integer.valueOf(args[++i]);
				}
				
				if(args[i].equals("-minPolyAReadCount")) {
					minPolyAReadCount = Integer.valueOf(args[++i]);
				}
				
				if(args[i].equals("-polyAUpperWindowCutoff")) {
					upperPolyACutoff = Double.valueOf(args[++i]);
				}
				if(args[i].equals("-polyALowerWindowCutoff")) {
					lowerPolyACutoff = Double.valueOf(args[++i]);
				}
				
				if(args[i].equals("-maxConsideredClippingLength")) {
					maxConsideredClippingLength = Integer.valueOf(args[++i]);
				}
				
				
				
				
				if(args[i].equals("--nospliceprio")) {
					preferExtensionsWithKnownSpliceSignals = false;
					continue;
				}
				
				if(args[i].equals("--noclipping")) {
					clipping = false;
					continue;
				}
				if(args[i].equals("--polyA")) {
					polyA = true;
					continue;
				}
				
				if(args[i].equals("--stranded_polyA")) {
					polyA = true;
					strandedPolyA = true;
					continue;
				}
				if(args[i].equals("--strandspecific")) {
					strandSpecific = true;
					continue;
				}
				if(args[i].equals("--pairedend")) {
					System.err.println(String.format("[%s]\tOption --pairedend deprecated. Paired-end runs are started by passing two file paths to the -reads argument. Please check the help by calling ContextMap without any parameters.",date.toLocaleString()));
					System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
					System.exit(1);
					continue;
				}
				if(args[i].equals("--printmultimappings")) {
					printMultiMappings = true;
					continue;
				}
				if(args[i].equals("--secondbestchr") || args[i].equals("--mining")) {
					printSecondBestChr = true;
					continue;
				}
				
				if(args[i].equals("--highcoverage")) {
					lowCoverage = false;
					continue;
				}
				
				if(args[i].equals("--sequenceDB")) {
					writeSequenceDB = true;
					continue;
				}
				
				if(args[i].equals("-dbSortBatchSize")) {
					databaseSortBatchSize = Integer.valueOf(args[++i]);
					continue;
				}
				
				
			}
			
			/*
			 * The run starts.
			 * First check the parameters
			 */

			if(alignerName == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-aligner_name"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			if(!alignerName.equals("bwa") && !alignerName.equals("bowtie1") && !alignerName.equals("bowtie") && !alignerName.equals("bowtie2")) {
				date = new Date();
				System.err.println(String.format("[%s]\tUnknown aligner name given. Please check the help by calling ContextMap without any parameters",date.toLocaleString()));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			//checking -reads argument
			if(readFilePath == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-reads"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			
			FileHandler fileHandler = new FileHandler();
			fileHandler.checkInput(filePaths);
			
			/*
			 * checking -o argument
			 */
			if(outputDirPath == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-o"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			if(!new File(outputDirPath).isDirectory()) {
				boolean createdDir = new File(outputDirPath).mkdirs();
				if(!createdDir) {
					date = new Date();
					System.err.println(String.format("[%s]\tValue of %s is not a valid directory and cannot be created. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-o"));
					System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
					System.exit(1);
				}
			}
			

			/*
			 * checking -annotation and -gtf options
			 * gtf option is mutually exclusive with annotation option
			 */
			if(gtfFilePath != null && annotationFilePath != null) {
				date = new Date();
				System.err.println(String.format("[%s]\tOptions -annotation and -gtf are mutually exclusive. Either use -annotation or -gtf, but not both.",date.toLocaleString()));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			if(!clipping && polyA) {
				date = new Date();
				System.err.println(String.format("[%s]\tOptions --noclipping and --polyA are mutually exclusive. Either use --noclipping or --polyA, but not both.",date.toLocaleString()));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			if(!pairedEnd && strandedPolyA) {
				date = new Date();
				System.err.println(String.format("[%s]\tThe option --stranded_polyA is not yet available for single-end runs. Please re-start ContextMap with --polyA.",date.toLocaleString()));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			
			/*
			 * alignment stuff
			 */
			if(alignerBinPath == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-aligner"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			if(!new File(alignerBinPath).isFile()) {
				date = new Date();
				System.err.println(String.format("[%s]\tThe path to the aligner executable is not a valid file path. Please check the help by calling ContextMap without any parameters",date.toLocaleString()));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}

			if(genomeIndexBasePaths.size() == 0) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-indices"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			
			if(indexDirPath != null && !printSecondBestChr) {
				date = new Date();
				System.err.println(String.format("[%s]\tWARNING: Found an indexed species list provided via the -speciesindex option, but --mining is not set. If you would like to mine for contaminations and infections please set the --mining option. Otherwise ignore this message.",date.toLocaleString()));
				System.err.println(String.format("[%s]\tContinuing ContextMap run.",date.toLocaleString()));
			}
			
			
			if(skipSplitDetection.size() > 0 && !printSecondBestChr) {
				date = new Date();
				System.err.println(String.format("[%s]\tWARNING: You have used the -skipsplit option, but --mining is not set. If you would like to mine for contaminations and infections please set the --mining option. Otherwise ignore this message.",date.toLocaleString()));
				System.err.println(String.format("[%s]\tContinuing ContextMap run.",date.toLocaleString()));
			}
			
			
			if(printSecondBestChr && indexDirPath == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter, which is required in mining mode. Please check the help by calling ContextMap without any parameters or see the manual",date.toLocaleString(), "-speciesindex"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			if(printSecondBestChr && skipSplitDetection.size() == 0) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter, which is required in mining mode. Please check the help by calling ContextMap without any parameters or see the manual",date.toLocaleString(), "-skipsplit"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			
			if(alignerName.equals("bowtie2") && printSecondBestChr && maxHits <= 3) {
				date = new Date();
				System.err.println(String.format("[%s]\tWARNING: When using bowtie2 as the underlying alignment tool in 'mining' mode, we recommend to increase the maximum number of candidate alignments per read (-maxhits)",date.toLocaleString()));
				System.err.println(String.format("[%s]\tContinuing ContextMap run.",date.toLocaleString()));
			}
			
			//check if split detection is enabled
			if(skipSplitDetection.size() == 0) {
				for(int i = 0; i < genomeIndexBasePaths.size(); i++) {
					skipSplitDetection.add(false);
				}
			}
			
			if(skipMultiSplitDetection.size() == 0) {
				for(int i = 0; i < genomeIndexBasePaths.size(); i++) {
					skipMultiSplitDetection.add(false);
				}
			}
			
			
			if(alignerIndexerPath == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing argument of the %s parameter. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-aligner"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			if(!new File(alignerIndexerPath).isFile()) {
				date = new Date();
				System.err.println(String.format("[%s]\tValue of the aligner's indexing program is not a valid file path. Please check the help by calling ContextMap without any parameters",date.toLocaleString()));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			/*
			 * context processing stuff (default settings are defined in the constructor)
			 */
			if(referencesDir == null) {
				date = new Date();
				System.err.println(String.format("[%s]\tMissing %s parameter. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-genome"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			if(!new File(referencesDir).isDirectory()) {
				date = new Date();
				System.err.println(String.format("[%s]\tValue of %s is not a valid directory. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), "-genome"));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			
			
			/*
			 * Currently, ContextMap only works with fasta files, therefore we convert fastq input here
			 */
			fileHandler.processInput(filePaths, outputDirPath + "/input_tmp", outputDirPath + "/tmp_reads.fa");
			readFilePath = outputDirPath + "/tmp_reads.fa";
			
			if(pairedEnd) {
				
				if(strandSpecific) {
					date = new Date();
					System.err.println(String.format("[%s]\tCurrently, ContextMap is not supporting strand-specific paired-end reads. Please restart ContextMap either in the single-end mode or without the '--strandspecific' option",date.toLocaleString()));
					System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
					System.exit(1);
				}
				
				if(printSecondBestChr) {
					date = new Date();
					System.err.println(String.format("[%s]\tCurrently, ContextMap is not able to print information of the second best chromosome/genome for paired-end reads. Please restart ContextMap either in the single-end mode or without the '--mining' option",date.toLocaleString()));
					System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
					System.exit(1);
				}
			}
			
			
			if(maxGapSize > maxContextSize) {
				maxGapSize = maxContextSize - 1;
			}
			
			
			
			ContextMap contextMap = new ContextMap(contextmapVersion,readFilePath, samFilePath,alignerName, alignerBinPath, alignerTmpDirPath, alignerIndexerPath, genomeIndexBasePaths, referencesDir,annotationFilePath,gtfFilePath, indexDirPath, outputDirPath, "fasta", readLength,
					seedSize,seedMismatches, maxMismatches,maxMismatchDifference, splitSeedSizes,splitSeedMismatches,maxIntronLength, maxIntronCount,maxHits, threads,
					skipPass1,skipSplitDetection,skipMultiSplitDetection,skipRealignment, minDistanceBetweenContext, maxContextSize,maxGapSize, minGapSize, maxIndelSize, minNumberOfReadsInContext,strandSpecific,pairedEnd,preferExtensionsWithKnownSpliceSignals,skipDenovoJunctions,skipNonCanonicalJunctions,updateQueue,updateInterval,printMultiMappings,printSecondBestChr,lowCoverage,clipping,polyA,strandedPolyA,minPolyALength,minPolyAReadCount,upperPolyACutoff,lowerPolyACutoff,maxConsideredClippingLength,verbose,developer,keepTmp,writeSequenceDB,autosetAlignmentOptions,databaseSortBatchSize);
			contextMap.start();
		}
		
		
		
		else if(args[0].equals("indexer")) {
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println("Concatenates the entries of a given multi-fasta file (e.g. all microbial genomes from refseq) and inserts a sequence of N's between two ");
				System.out.println("entries to omit read alignments overlapping subsequent entries.");
				System.out.println("The output is a set of fasta files, each containing at most a pre-defined maximum number of bases (default: 250000000).");
				System.out.println("Furthermore, \"*.idx\" files will be created, which allow to reconstruct the original sequences of the input.");
				System.out.println("This tool was created to be able to index a huge amount of genomes with Bowtie (which is limited to 10.000 genomes per index).");
				System.out.println("Alignments to the concatenated sequences can then be analyzed with other tools contained in the ContextMap package.");
				
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar indexer <arguments> [options]*",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-fasta\t\t<The path to a multi fasta file, which should be indexed>");
				System.out.println("-prefix\t\t<The basename of the index (e.g. 'microbes')>");
				System.out.println("-o\t\t<The path to the output directory>");
				System.out.println();
				System.out.println();
				System.out.println("Options:");
				System.out.println();
				System.out.println("-gapsize\t<Defines the number of N's inserted between two subsequent entries of the input. See the help for an explanation <default: 100>");
				System.out.println("-indexsize\t<The maximum number of bases a index file contains. <default: 250000000>");
				System.out.println();
				System.exit(1);
			}
			
			String genomesFilePath = null;
			String outputDirPath = null;
			String prefix = null;
			int gapSize = 100;
			int maxBasePairsPerIndex = 250000000;
			//in case the index should be restricted to a specific species set, 
			//this set will be extracted from the second column (trivial name) of a tab separated file.
			//this was added for testing, can be deleted later....
			String microbialContentFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-fasta")) { 
					genomesFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-prefix")) { 
					prefix = args[++i];
					continue;
				}
				if(args[i].equals("-microbialcontent")) { 
					microbialContentFilePath = args[++i];
					continue;
				}
				
				if(args[i].equals("-gapsize")) { 
					gapSize = Integer.valueOf(args[++i]);
					continue;
				}
				
				if(args[i].equals("-indexsize")) { 
					maxBasePairsPerIndex = Integer.valueOf(args[++i]);
					continue;
				}
				
				if(args[i].equals("-o")) {
					outputDirPath = args[++i];
					continue;
				}
			}
			
			GenomeIndexer genomeIndexer = new GenomeIndexer(genomesFilePath, outputDirPath, gapSize, maxBasePairsPerIndex);
			genomeIndexer.indexWithConcatenation(prefix,microbialContentFilePath);
		}
		
		
		else if(args[0].equals("inspector")) {
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println();
				System.out.println("Determines all species with mappings in a given sam file and outputs read counts, genome coverages, confidence as well as normalized Jensen-Shannon divergence values for every identified species.");
				
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar inspector <arguments> [options]*",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-sam\t\t<The path to a mapping file in sam format>");
				System.out.println("-idx\t\t<The path to a directory containing index files, which were generated with the 'indexer' tool>");
				System.out.println();
				System.out.println();
				System.out.println("Options:");
				System.out.println();
				System.out.println("-reference\t<A comma seperated list of chr names or a single species id, which will be used as reference for the JS divergence calculation. If not set, the most abundant species in the sample will be used.");
				System.out.println("-refseqcatalog\t<The path to a RefSeq *.catalog file. Trivial species names will be extracted in case they are not available in the index files.>");
				System.out.println("--mergecontigs\t<Mappings to different contigs of the same species will be merged>");
				System.out.println("--mdflag\t<Uses the MD field of the sam file to evaluate mismatch counts. Per default, the NM field is used.>");
				System.out.println();
				System.exit(1);
			}
			
			
			
			String samFilePath = null;
			String indexFolderPath = null;
			String refseqCatalogFilePath = null;
			String customSpeciesNamesFilePath = null;
			boolean mergeContigs = false;
			boolean useMDflag = false;
			HashSet<String> referenceChromosomes = new HashSet<String>();
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-idx")) {
					indexFolderPath = args[++i];
					continue;
				}
				if(args[i].equals("-reference")) {
					String[] references = args[++i].split(",");
					for(String ref : references)
						referenceChromosomes.add(ref);
					continue;
				}
				if(args[i].equals("-customnames")) {
					customSpeciesNamesFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-refseqcatalog")) {
					refseqCatalogFilePath = args[++i];
					continue;
				}
				if(args[i].equals("--mergecontigs")) {
					mergeContigs = true;
					continue;
				}
				if(args[i].equals("--mdflag")) {
					useMDflag = true;
				}
			}
			
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.getMicrobialContentFromConcatenatedGenomes(samFilePath, indexFolderPath,refseqCatalogFilePath,customSpeciesNamesFilePath,mergeContigs,useMDflag,referenceChromosomes);
		}
		
		
		else if(args[0].equals("separator")) {
			
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println("Splits a given sam file into mappings to individual species.");
				
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar separator <arguments> [options]*",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-sam\t\t<The path to a mapping file in sam format>");
				System.out.println("-idx\t\t<The path to a directory containing index files, which were generated with the 'indexer' tool>");
				System.out.println("-o\t\t<The path to the output directory>");
				System.out.println();
				System.out.println();
				System.out.println("Options:");
				System.out.println();
				System.out.println("-specieslist\t<The path to a tab seperated file containing for each species its id and trivial name. The output will be limited to the set of given species.>");
				System.out.println("--mergecontigs\t<Mappings to different contigs of the same species will be merged>");
				System.out.println("--trivialnames\t<The generated files will be named by trivial names of the species>");
				System.out.println();
				System.exit(1);
			}
			
			
			
			String samFilePath = null;
			String indexFolderPath = null;
			String speciesListFilePath = null;
			String outputDirPath = null;
			boolean mergeContigs = false;
			boolean trivialFileNames = false;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-idx")) {
					indexFolderPath = args[++i];
					continue;
				}
				if(args[i].equals("-specieslist")) {
					speciesListFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputDirPath = args[++i];
					continue;
				}
				if(args[i].equals("--mergecontigs")) {
					mergeContigs = true;
					continue;
				}
				if(args[i].equals("--trivialnames")) {
					trivialFileNames = true;
					continue;
				}
			}
			
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.separateMappingBySpecies(samFilePath,indexFolderPath,speciesListFilePath,outputDirPath,mergeContigs,trivialFileNames);
		}
		
		else if(args[0].equals("error_estimator")) {
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println("Determines the mismatch distribution for a set of species contained in a given sam file (requires the NM or MD field).");
				
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar errorestimator <arguments> [options]*",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-sam\t\t<The path to a mapping file in sam format>");
				System.out.println("-idx\t\t<The path to a directory containing index files, which were generated with the 'indexer' tool>");
				System.out.println();
				System.out.println();
				System.out.println("Options:");
				System.out.println();
				System.out.println("-specieslist\t<The path to a tab seperated file containing for each species its id and trivial name. The output will be limited to the set of given species.>");
				System.out.println("--mdflag\t<Uses the MD field of the sam file to evaluate mismatch counts. Per default, the NM field is used.>");
				System.out.println();
				System.exit(1);
			}
			
			
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			boolean useMDflag = false;
			HashSet<String> referenceChromosomes = new HashSet<String>();
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-specieslist")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("-reference")) {
					String[] references = args[++i].split(",");
					for(String ref : references)
						referenceChromosomes.add(ref);
					continue;
				}
				if(args[i].equals("--mdflag")) {
					useMDflag = true;
				}
				
			}
			new Statistic().generateErrorRateTableFromConcatenatedGenomes(samFilePath,microbialContentFilePath,genomeIndexDirPath,useMDflag,referenceChromosomes,null);
		}
		
		
		else if(args[0].equals("confidence_estimator")) {
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println("Determines confidence values for a set of species contained in a given sam file, which was produced with ContextMap");
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar confidence_estimator <arguments>",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-sam\t\t<The path to a mapping file in sam format>");
				System.out.println("-idx\t\t<The path to a directory containing index files, which were generated with the 'indexer' tool>");
				System.out.println("-specieslist\t<The path to a tab seperated file containing for each species its id and trivial name. The output will be limited to the set of given species.>");
				System.out.println();
				System.exit(1);
			}
			
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			boolean mergeYeastChr = true;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-specieslist")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("--noyeastmerge")) {
					mergeYeastChr = false;
				}
			}
			new Statistic().generateConfidenceScoresFromConcatenatedGenomes(samFilePath, microbialContentFilePath,genomeIndexDirPath,mergeYeastChr,null);
		}
		
		
		else if(args[0].equals("distance_estimator")) {
			
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println("Determines distance values for a set of species contained in a given sam file, which was produced with ContextMap");
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar distance_estimator <arguments>",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-sam\t\t<The path to a mapping file in sam format>");
				System.out.println("-idx\t\t<The path to a directory containing index files, which were generated with the 'indexer' tool>");
				System.out.println("-specieslist\t<The path to a tab seperated file containing for each species its id and trivial name. The output will be limited to the set of given species.>");
				System.out.println();
				System.exit(1);
			}
			
			
			
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			boolean mergeYeastChr = true;
			HashSet<String> speciesFilterNames = new HashSet<String>();
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-specieslist")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("-filter")) {
					String[] names = args[++i].split(",");
					speciesFilterNames.addAll(Arrays.asList(names));
				}
				if(args[i].equals("--noyeastmerge")) {
					mergeYeastChr = false;
				}
			}
			new Statistic().generateDistanceTableFromConcatenatedGenomes(samFilePath, microbialContentFilePath,genomeIndexDirPath,speciesFilterNames,mergeYeastChr);
		}
		
		
		else if(args[0].equals("mdfield_generator")) {
			
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println("Determines the md fields for a given sam file. Assumes that the sequence field of the sam file is set.");
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar mdfield_generator <arguments>",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-sam\t\t<The path to a sam file>");
				System.out.println("-ref\t\t<The path to a directory containing reference fasta files>");
				System.out.println("-o\t\t<The output file path>");
				System.out.println("-t\t\t<The number of threads to use for processing>");
				System.out.println();
				System.exit(1);
			}
			
			
			
			String samFilePath = null;
			String referencesDirPath = null;
			String outputFilePath = null;
			int threads  = 1;
			boolean considerClippedRegions = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-ref")) {
					referencesDirPath = args[++i];
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
				}
				if(args[i].equals("-t")) {
					threads = Integer.valueOf(args[++i]);
				}
			}
			new SamProcessor().addMdFlag(samFilePath, outputFilePath, new File(referencesDirPath).listFiles(),considerClippedRegions, threads);
		}
		
		else if(args[0].equals("seqfield_generator")) {
			
			if(args.length == 1 || args[1].equals("--help") || args[1].equals("-h")) {
				System.out.println();
				System.out.println("Determines the sequence fields for a given sam file");
				System.out.println();
				System.out.println(String.format("Usage: java -jar ContextMap_%s.jar seqfield_generator <arguments>",contextmapVersion));
				System.out.println();
				System.out.println("Required arguments:");
				System.out.println();
				System.out.println("-sam\t\t<The path to a sam file>");
				System.out.println("-reads\t\t<The path to the reads file, either in fasta or fastq format>");
				System.out.println("-o\t\t<The output file path>");
				System.out.println();
				System.out.println();
				System.out.println("Options:");
				System.out.println();
				System.out.println("--fastq\t\t<Must be set if the input is given in fastq format>");
				System.out.println("--quality\t\t<The quality field will also be set>");
				System.out.println("--prebuffer\t<Enables prebuffering of all read sequences>");
				
				System.out.println();
				System.exit(1);
			}
			
			
			
			String samFilePath = null;
			String readsFilePath = null;
			String outputFilePath = null;
			String readFormat = "fasta";
			boolean prebuffer = false;
			boolean addQuality = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-reads")) {
					readsFilePath = args[++i];
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
				}
				if(args[i].equals("--fastq")) {
					readFormat = "fastq";
				}
				if(args[i].equals("--quality")) {
					addQuality = true;
				}
				if(args[i].equals("--prebuffer")) {
					prebuffer = true;
				}
			}
			
			
			if(addQuality && !readFormat.equals("fastq")) {
				System.err.println("Quality field can only be set if the underlying read data are given in fastq format (--fastq option).");
				System.err.println("Aborting.");
				System.exit(1);
				
			}
			
			new SamProcessor().addSequenceField(samFilePath, outputFilePath, readsFilePath, readFormat, prebuffer,addQuality);
		}
		
		
		
		else if(args[0].equals("addmateinformation")) {
			boolean toId = false;
			String inputFilePath = null;
			String outputFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					inputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("--id")) { 
					toId = true;
					continue;
				}
			}
			
			SamProcessor sp = new SamProcessor();
			if(!toId)
				sp.addMateInformationToFlag(inputFilePath, outputFilePath);
			
			else
				sp.addMateInformationToId(inputFilePath, outputFilePath);
		}
		
		else if(args[0].equals("removematefromid")) {
			String inputFilePath = null;
			String outputFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					inputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
			}
			
			SamProcessor sp = new SamProcessor();
			sp.removeMateInformationFromId(inputFilePath, outputFilePath);
		}
		
		
		
		
		else if(args[0].equals("unmappedReadsExtractor")) {
			
			String samFilePath = null;
			String readFilePath = null;
			String readFormat = "fasta";
			int readLength = -1;
			String outputDir = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-reads")) {
					readFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-format")) {
					readFormat = args[++i];
					continue;
				}
				if(args[i].equals("-readlen")) {
					readLength = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-o")) {
					outputDir = args[++i];
					continue;
				}
			}
			System.out.println("sam file path:\t" + samFilePath);
			System.out.println("read file path:\t" + readFilePath);
			System.out.println("read format:\t" + readFormat);
			System.out.println("read length:\t" + readLength);
			System.out.println("output dir:\t" + outputDir);
			SamProcessor samFileProcessor = new SamProcessor();
			samFileProcessor.extractUnmappedAndMultimappedReads(samFilePath, readFilePath, readFormat, readLength, outputDir);
		}
		
		else if(args[0].equals("aligner")) {
			
			//required fields
			String rmapExecutablePath = null;
			String fastaFilePath = null;
			String genomeIndexBasePath = null;
			String outputPath = null;
			String tmpDirPath = null;
			String bowtieBuildExecutablePath = null;
			
			//optional fields (defalut settings are defined here)
			String seedLength = "40";
			String seedMissmatches = "2";
			String maxMissmatches = "6";
			String splitSeedSizes = "10,15,20";
			String maxIntronLength = "300000";
			String maxIntronCount = "10";
			String maxInsertionSize = "5";
			//reads with more than maxHits mappings will be skipped and printed into a seperate file
			int maxHits = -1;
			String skippedReadsPath = null;
			String unalignedReadsPath = null;
			boolean prebufferReads = false;
			String threads = "1";
			String sortMemory = "500";
			boolean skipFull = false;
			boolean skipSplitDetection = false;
			boolean verbose = true;
			//rmap command line usage: -i <fastafile> -index <bowtie genome index> -o <outfile> -tmpdir <tmpdir> -btbuild <path  to bowtie-build executable>
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-rmap")) { 
					rmapExecutablePath = args[++i];
					continue;
				}
				if(args[i].equals("-i")) {
					fastaFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-index")) {
					genomeIndexBasePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("-tmpdir")) {
					tmpDirPath = args[++i];
					continue;
				}
				if(args[i].equals("-btbuild")) {
					bowtieBuildExecutablePath = args[++i];
					continue;
				}
				if(args[i].equals("-t")) {
					threads = args[++i];
					continue;
				}
				
				if(args[i].equals("-mm")) {
					maxMissmatches = args[++i];
					continue;
				}
				
				if(args[i].equals("-seed")) {
					seedLength = args[++i];
					continue;
				}
				if(args[i].equals("-skipsplit")) {
					skipSplitDetection = Boolean.valueOf(args[++i]);
					continue;
				}
				
			}
			System.out.println("rmap executable:\t" + rmapExecutablePath);
			System.out.println("read file path:\t" + fastaFilePath);
			System.out.println("genome index:\t" + genomeIndexBasePath);
			System.out.println("output path:\t" + outputPath);
			System.out.println("tmp dir path:\t" + tmpDirPath);
			System.out.println("bowtie build executable:\t" + bowtieBuildExecutablePath);
			System.out.println("max missmatches per read:\t" + maxMissmatches);
			System.out.println("seed length:\t" + seedLength);
			ReadAligner readAligner = new ReadAligner();
			readAligner.alignReads(rmapExecutablePath, fastaFilePath, genomeIndexBasePath, outputPath, tmpDirPath, bowtieBuildExecutablePath,seedLength,seedMissmatches,maxMissmatches,splitSeedSizes,maxIntronLength,maxIntronCount,maxInsertionSize,maxHits,skippedReadsPath,unalignedReadsPath,prebufferReads,threads,sortMemory,skipFull,skipSplitDetection,verbose);
		}
		
		
		else if(args[0].equals("sam2rmap")) {
			
			String samFilePath = null;
			String outputPath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
			}
			System.out.println("sam file path:\t" + samFilePath);
			System.out.println("output path:\t" + outputPath);
			SamProcessor samFileProcessor = new SamProcessor();
			samFileProcessor.convertSamToRmap(samFilePath, outputPath);
		}
		
		else if(args[0].equals("rmap2modifiedRmap")) {
			
			String inputPath = null;
			String outputPath = null;
			boolean skipPartialHits = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("--skippartial")) {
					skipPartialHits = true;
					continue;
				}
			}
			System.out.println("sam file path:\t" + inputPath);
			System.out.println("output path:\t" + outputPath);
			
			RmapProcessor rmapProcessor = new RmapProcessor();
			rmapProcessor.modifyRmapAnnotation(inputPath,outputPath,skipPartialHits,1);
		}
		
		else if(args[0].equals("concatenateRmap")) {
			ArrayList<String> rmaps = new ArrayList<String>();
			String outputPath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					String[] tmpRmaps = args[++i].split(",");
					rmaps.addAll(Arrays.asList(tmpRmaps));
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
			}
			RmapProcessor rmapProcessor = new RmapProcessor();
			rmapProcessor.concatenateFilesWithNIO(rmaps, outputPath);
		}
		
		else if(args[0].equals("splitRmap")) {
			
			String inputPath = null;
			String outputDirPath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputDirPath = args[++i];
					continue;
				}
			}
			System.out.println("sam file path:\t" + inputPath);
			System.out.println("output dir path:\t" + outputDirPath);
			
			RmapProcessor rmapProcessor = new RmapProcessor();
			rmapProcessor.splitByChromosome(inputPath,outputDirPath);
		}
		

		
		
		else if(args[0].equals("extractBestMatchings")) {
			String multiMappingFilePath = null;
			String outputPath = null;
			int maxMissmatchDifference = 0;
			boolean localContext = true;
			boolean preferExtensionsWithKnownSpliceSignals = true;
			boolean skipDenovoJunctions = false;
			boolean skipNonCanonicalJunctions = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					multiMappingFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				if(args[i].equals("-mmdiff")) {
					maxMissmatchDifference = Integer.valueOf(args[++i]);
					continue;
				}
			}
			
			MappingProcessor mappingProcessor = new MappingProcessor(multiMappingFilePath, outputPath, maxMissmatchDifference,preferExtensionsWithKnownSpliceSignals,skipDenovoJunctions,skipNonCanonicalJunctions);
			mappingProcessor.extractBestMatchingAlignmentsInLocalResolution(null,localContext);
			
		}
		
		
		else if(args[0].equals("extractMicrobialContentFromUnconcatenatedGenomes")) {
			String samFilePath = null;
			String indexFolderPath = null;
			String refseqCatalogFilePath = null;
			String customSpeciesNamesFilePath = null;
			boolean mergeContigs = false;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-idx")) {
					indexFolderPath = args[++i];
					continue;
				}
				if(args[i].equals("-customnames")) {
					customSpeciesNamesFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-refseqcatalog")) {
					refseqCatalogFilePath = args[++i];
					continue;
				}
				if(args[i].equals("--mergecontigs")) {
					mergeContigs = true;
					continue;
				}
			}
			
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.getMicrobialContentFromUnconcatenatedGenomes(samFilePath, indexFolderPath,refseqCatalogFilePath,customSpeciesNamesFilePath,mergeContigs);
		}
		
		else if(args[0].equals("extractSequences")) {
			String samFilePath = null;
			String fastaFilePath = null;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-fasta")) {
					fastaFilePath = args[++i];
					continue;
				}
			}
			
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.extractSequences(fastaFilePath, samFilePath);
		}
		
		
		
		
		
		/*
		 * expects a sam file of a single species. uses the index file and the first line of the sam file to
		 * determine relative genome start of the species and outputs the absolute genome coordinates
		 */
		else if(args[0].equals("getCoverage")) {
			String samFilePath = null;
			String outputFilePath = null;
			int readLength = -1;
			String indexFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-rl")) {
					readLength = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-idx")) {
					indexFilePath = args[++i];
					continue;
				}
				
			}
			
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.getCoverage(samFilePath, outputFilePath, indexFilePath, readLength);
		}
		
		else if(args[0].equals("indexGenomesWithoutConcatenation")) {
			String genomesFilePath = null;
			String outputDirPath = null;
			String prefix = null;
			int gapSize = 100;
			int maxBasePairsPerIndex = 250000000;
			//in case the index should be restricted to a specific species set, 
			//this set will be extracted from the second column (trivial name) of a tab separated file.
			//this was added for testing, can be deleted later....
			String microbialContentFilePath = null;
			long maxBasesPerIndex = Long.MAX_VALUE;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-genomes")) { 
					genomesFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-prefix")) { 
					prefix = args[++i];
					continue;
				}
				if(args[i].equals("-basesperindex")) { 
					maxBasesPerIndex = Long.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-microbialcontent")) { 
					microbialContentFilePath = args[++i];
					continue;
				}
				
				if(args[i].equals("-o")) {
					outputDirPath = args[++i];
					continue;
				}
			}
			
			GenomeIndexer genomeIndexer = new GenomeIndexer(genomesFilePath, outputDirPath, gapSize, maxBasePairsPerIndex);
			genomeIndexer.indexWithoutConcatenation(prefix,maxBasesPerIndex,microbialContentFilePath);
		}
		
		
		
		
		else if(args[0].equals("filterAssemblyByLength")) {
			String assemblyFilePath = null;
			int minLength = 0;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-assembly")) { 
					assemblyFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-minlength")) { 
					minLength = Integer.valueOf(args[++i]);
					continue;
				}
			}
			AssemblyProcessor assemblyProcessor = new AssemblyProcessor();
			assemblyProcessor.filterByLength(assemblyFilePath, minLength);
		}
		
		else if(args[0].equals("assemblyStatistic")) {
			String assemblyFilePath = null;
			int minLength = 0;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-assembly")) { 
					assemblyFilePath = args[++i];
					continue;
				}
			}
			AssemblyProcessor assemblyProcessor = new AssemblyProcessor();
			assemblyProcessor.getAssemblyStatistic(assemblyFilePath);
		}
		
		else if(args[0].equals("sampleReads")) {
			String readFilePath = null;
			String outputFilePath = null;
			int sampleSize = 0;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					readFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-samplesize")) { 
					sampleSize = Integer.valueOf(args[++i]);
					continue;
				}
			}
			ReadSampler sampler = new ReadSampler();
			sampler.sampleReads(readFilePath, sampleSize, outputFilePath);
		}
		
		else if(args[0].equals("-samplingstats")) {
			String samplingFolderPath = null;
			double minCov = 0.0;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					samplingFolderPath = args[++i];
					continue;
				}
				if(args[i].equals("-mincov")) { 
					minCov = Double.valueOf(args[++i]);
					continue;
				}
				
			}
			ReadSampler sampler = new ReadSampler();
			sampler.getSampleStatistic(samplingFolderPath,minCov);
		}
		
		
		else if(args[0].equals("addmismatchpositions")) {
			String bamFilePath = null;
			String indexFilePath = null;
			String genomeDirPath = null;
			String fastaFilePath = null;
			String outputFilePath = null ;
			boolean checkModification = false;
			boolean prebuffer = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-bam")) { 
					bamFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-index")) { 
					indexFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-genome")) { 
					genomeDirPath = args[++i];
					continue;
				}
				if(args[i].equals("-reads")) { 
					fastaFilePath = args[++i];
					continue;
				}
				if(args[i].equals("--check")) { 
					checkModification = true;
					continue;
				}
				if(args[i].equals("--prebuffer")) { 
					prebuffer = true;
					continue;
				}
			}
			
			BamConverter converter = new BamConverter();
			converter.addMismatchPositionsToCigarString(bamFilePath,indexFilePath,genomeDirPath,fastaFilePath,outputFilePath,checkModification,prebuffer);
		}
		
		else if(args[0].equals("getseedmismatches")) {
			String bamFilePath = null;
			String indexFilePath = null;
			int seedlen = Integer.MIN_VALUE;
			int maxMismatches = Integer.MAX_VALUE;
			String outputFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-bam")) { 
					bamFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-index")) { 
					indexFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-seedlen")) { 
					seedlen = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-mismatches")) { 
					maxMismatches = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
				
			}
			
			BamConverter converter = new BamConverter();
			converter.getSeedMismatchStats(bamFilePath, indexFilePath, seedlen,maxMismatches,outputFilePath);
		}
		
		else if(args[0].equals("modifyPairedEndStrandInformation")) {
			String samFilePath = null;
			String outputFilePath = null ;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
				
			}
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.modifyPairedEndStrandInformation(samFilePath, outputFilePath);
		}
		
		else if(args[0].equals("fixsam")) {
			String samFilePath = null;
			String outputFilePath = null ;
			int readlen = Integer.MIN_VALUE;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) { 
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-readlen")) { 
					readlen = Integer.valueOf(args[++i]);
					continue;
				}
				
			}
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.fixSam(samFilePath,readlen, outputFilePath);
		}
		
		
		
		else if(args[0].equals("sortFile")) {
			int[] columnsToSortFor = null;
			int linesPerChunkFile = 2000000;
			String delimiter = "\t";
			boolean numericalSort = false;
			String inputPath = null;
			String outputPath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					inputPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputPath = args[++i];
					continue;
				}
				
				if(args[i].equals("-c")) {
					String[] tmpColumns = args[++i].split(",");
					columnsToSortFor = new int[tmpColumns.length];
					for(int j = 0; j < tmpColumns.length; j++) {
						columnsToSortFor[j] = Integer.valueOf(tmpColumns[j]);
					}
					continue;
				}
				if(args[i].equals("-d")) {
					delimiter = args[++i];
					continue;
				}
				
				if(args[i].contains("-l")) {
					linesPerChunkFile = Integer.valueOf(args[++i]);
					continue;
				}
				if(args[i].contains("--numerical")) {
					numericalSort = true;
					continue;
				}
			}
			
			FileSorter fileSorter = new FileSorter(columnsToSortFor,linesPerChunkFile, delimiter,numericalSort);
			fileSorter.sortFile(inputPath, outputPath,null);
		}
		
		else if(args[0].equals("getStartPositionsOnly")) {
			String inputFilePath = null;
			String outputFilePath = null;
		
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) { 
					inputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) { 
					outputFilePath = args[++i];
					continue;
				}
			}
			
			SamProcessor samProcessor = new SamProcessor();
			samProcessor.getStartPositionsOnly(inputFilePath, outputFilePath);
		}
		
		
		else {
			Date date = new Date();
			System.out.println();
			System.err.println(String.format("[%s]\tUnknown argument: '%s'. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), args[0]));
			System.exit(1);
		}
		
	}
	
	/**
	 * should only be called in paired end mode
	 * 
	 **/
	
	private static boolean hasNewIlluminaHeader(String fastaFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			String header = br.readLine().substring(1);
			br.close();
			if(header.contains(" ")) {
				String[] splittedHeader = header.split(" ")[1].split(":");
				if(splittedHeader.length == 4 && (splittedHeader[0].equals("1") || splittedHeader[0].equals("2")) && (splittedHeader[1].equals("Y") || splittedHeader[1].equals("N"))) {
					return true;
				}
			}
			
			return false;
		}
		
		catch(Exception e) {
			e.printStackTrace();
			return false;
		}
	}
	
	/**
	 * here we already know that the input contains the new illumina header
	 * @param inputFilePath
	 * @param outputFilePath
	 */
	
	private static void renameNewIlluminaHeader(String inputFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedHeader;
			Pattern spacePattern = Pattern.compile(" ");
			Pattern doublePointPattern = Pattern.compile(":");
			
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					splittedHeader = spacePattern.split(currentLine);
					splittedHeader[0] += "/" + doublePointPattern.split(splittedHeader[1])[0];
					pw.println(splittedHeader[0]);
				}
				else {
					pw.println(currentLine);
				}
			}
			br.close();
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
}
