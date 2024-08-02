---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


ContextMap is an open source read mapping tool for RNA-Seq data and released under the GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007.

For a description of the ContextMap algorithm and implementation please read the following publications:

ContextMap 2: Fast and accurate context-based RNA-seq mapping. BMC Bioinformatics
Thomas Bonfert, Evelyn Kirner, Gergely Csaba, Ralf Zimmer, Caroline C. Friedel. 2015, accepted, BMC Bioinformatics. 


A context-based approach to identify the most likely mapping for RNA-Seq experiments.
Thomas Bonfert, Gergely Csaba, Ralf Zimmer and Caroline C. Friedel. 2012, BMC Bioinformatics.


For using ContextMap to mine for contaminations and infections please read the following publication:

Mining RNA-seq Data for Infections and Contaminations.
Thomas Bonfert, Gergely Csaba, Ralf Zimmer, Caroline C. Friedel. 2013, PloS one.


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


- PLEASE READ THE MANUAL (available in the ContextMap package, 'manual' directory) to get an overview of ContextMap and to learn more about its usage.

- EXAMPLE DATA AND SHELL SCRIPTS TO DEMONSTRATE THE USAGE OF CONTEXTMAP are included in the ContextMap package ('example_call' and 'mining_example' directories).

- For questions or bug reports please mail to: friedel@bio.ifi.lmu.de


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


- Current release version: 2.7.9

- Usage: java -jar ContextMap_v2.*.jar <tool name> <tool arguments> [tool options]*

Available tools:
1. 'mapper' - The ContextMap mapping tool
2. 'indexer' - Prepares large sets of genomes for being indexed with the desired unspliced Aligner (Bowtie1, Bowtie2 or BWA).
3. 'inspector' - Determines read counts, confidence values, genome coverages and the square root of the Jensen-Shannon divergence for species contained in a sam file

Start a tool without any argument to get a specific help message and a list of available options for this tool


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


- Usage of the mapper tool: 
    java -jar ContextMap_v2.*.jar mapper <arguments> [options]*

Required arguments:

    -reads          <A comma separated list of file paths to reads in fasta/fastq/zip/gz format. A single path for single-end mapping and two paths (#1 mates and #2 mates) for paired-end mapping>
    -aligner_name   <The following alignment tools are supported: "bwa", "bowtie1" or "bowtie2">
    -aligner_bin    <The absolute path to the executable of the chosen aligner tool>
    -indexer_bin    <The absolute path to the executable of the aligner's indexing tool (not needed for BWA)>
    -indices        <A comma separated list of paths to basenames of indices, which can be used by the chosen aligner>
    -genome         <The path to a directory with genome sequences in fasta format (for each chromosome a separate file)>
    -o              <The path to the output directory>

Options:

    --mining                <Set this to mine for infections or contaminations. Changes of ContextMap's mapping behaviour due to this option are described in the manual.>
    -skipsplit              <A comma separated list of booleans, each element refers to a given aligner index (same ordering). 'true' for no split detection, 'false' otherwise (REQ. in mining mode).>
    -speciesindex           <The path to a directory containing index files created with the 'indexer' tool (REQ. in mining mode)>
    -aligner_tmp            <The path to a directory for temporary alignment files>
    -seed                   <The seed length for the alignment> (default: Bwt1: 30, BWA/Bwt2: 20)
    -seedmismatches         <Allowed mismatches in the seed region> (default: Bwt1: 1, BWA/Bwt2: 0)
    -mismatches             <Allowed mismatches in the whole read> (default: 4)
    -mmdiff                 <The maximum allowed mismatch difference between the best and second best alignment of the same read> (default: 0)
    -maxhits                <The maximum number of candidate alignments per read. Reads with more hits are skipped (bwa/bwt1) or the already found hits are reported (bwt2) > (default for bwa/bwt1:10, bwt2: 3)
    -minsize                <The minimum number of reads a genomic region has to contain for being regarded as a local context. (default:10)
    -annotation             <The path to an annotation file in our own format>
    -gtf                    <The path to an annotation file in gtf format. This option is mutually exclusive with the -annotation option>
    -t                      <The number of threads used for the run> (default: 1)
    --noclipping		<Disables the calculation of clipped alignments> (default: not set)";
    -updateinterval         <The number of performed queue updates during global context resolution (default: 3)>
    --updatequeue           <Enables queue updates for the local and global context resolution step> (default: not set)
    --noncanonicaljunctions <Enables the prediction of non-canonical splice sites> (default: not set)
    --strandspecific        <Enables strand specific mapping (default: not set)>
    --polyA			<Enables the search for polyA-tails and writes found cleavage sites to a bed file. (default: off)>
    --sequenceDB		<Writes a readId -> sequence mapping to disk. Recommended for very large data sets. (default: off)>
    --verbose               <verbose mode> (default: not set)
    --keeptmp               <does not delete some temporary files> (default: not set)


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


- Usage of the indexer tool: 
    java -jar ContextMap_v2.*.jar indexer <arguments> [options]*

Required arguments:

    -fasta		<The path to a multi fasta file, which should be indexed>
    -prefix		<The basename of the index (e.g. 'microbes')>
    -o		<The path to the output directory>


Options:

    -gapsize	<Defines the number of N's inserted between two subsequent entries of the input. See the help for an explanation <default: 100>
    -indexsize	<The maximum number of bases a index file contains. <default: 250000000>


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

- Usage of the inspector tool: 
    java -jar ContextMap_v2.*.jar inspector <arguments> [options]*

Required arguments:

    -sam		<The path to a mapping file in sam format>
    -idx		<The path to a directory containing index files, which were generated with the 'indexer' tool>


Options:

    -reference	<A comma seperated list of chr names or a single species id, which will be used as reference for JS divergence calculations. If not set, the most abundant species will be used.
    -refseqcatalog	<The path to a RefSeq *.catalog file. Trivial species names will be extracted in case they are not available in the index files.>
    --mergecontigs	<Mappings to different contigs of the same species will be merged>
    --mdflag	<Uses the MD field of the sam file to evaluate mismatch counts. Per default, the NM field is used.>


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Building ContextMap from source

See the manual for ContextMap.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Example call:

The example call for a usual ContextMap run can be found in the directory 'example_call'. Simply execute the script example_call.sh to try ContextMap. To use ContextMap with your own data, change the paths within the script.

The package also contains the directory 'mining_example', which itself contains data for demonstrating the mining process. By executing the 'run_example.sh' script, you can try the mining process. For using it with your own data, change the paths within the script.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

External libraries:

ContextMap makes use of a unmodified version of the MapDB library (http://www.mapdb.org), which is released under the Apache 2 license. You can find the source code of MapDB in the 'libs' directory of the ContextMap package.

