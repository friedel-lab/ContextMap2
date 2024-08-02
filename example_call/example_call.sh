#definition of input paths
#
read -p "Which unspliced aligner do you wish to use? (bwa,bowtie1,bowtie2) " aligner_name
    case $aligner_name in
    	bwa )  ;;
	bowtie1 ) ;;
	bowtie2 ) ;;
	* ) echo "Please answer bwa, bowtie1 or bowtie2";
    esac

if [ "$aligner_name" = bwa ]
	then read -p "Please enter the path for the bwa binary " aligner_bin
	indexer_bin=$aligner_bin
elif [ "$aligner_name" = bowtie1 ]
	then read -p "Please enter the path for the bowtie1 binary " aligner_bin
	read -p "Please enter the path for the bowtie1 indexer " indexer_bin
else read -p "Please enter the path for the bowtie2 binary " aligner_bin
	read -p "Please enter the path for the bowtie2 indexer " indexer_bin
fi

#output path is the current working directory. the jar file ashould be located in the parent directory.
#
output_dir=`pwd`
parent_dir="$(dirname "$output_dir")"

#
# Simply modify the following paths to start a ContextMap run with your own data
#
contextmap_jar=$parent_dir/ContextMap_v2.7.9.jar
reads=$output_dir/human_demo_set.fa
genomes_dir=$output_dir/reference_sequences
aligner_index=$output_dir/aligner_index/demo_index

#
# building indices with the original indexer
#
mkdir $output_dir/aligner_index/ >> log.txt 2>&1

if [ "$aligner_name" = bowtie1 -o "$aligner_name" = bowtie2 ]
then ${indexer_bin} -f $genomes_dir/chr1_ENSG00000198746.fa ${aligner_index} >> log.txt 2>&1
else ${indexer_bin} index -p ${aligner_index} $genomes_dir/chr1_ENSG00000198746.fa >> log.txt 2>&1
fi

#
# starting contextmap
#

java -Xms1000m -Xmx1000m -jar $contextmap_jar mapper -aligner_name $aligner_name -aligner_bin $aligner_bin -indexer_bin $indexer_bin -reads $reads -o $output_dir -indices $aligner_index -genome $genomes_dir >> log.txt 2>&1
echo ""


#
# for starting a multi threaded run, you have to add the '-t' option. Together with some performance tuning for the virtual machine the call will look like this (using 8 threads):
#
# java -Xms1000m -Xmx4000m -XX:+UseConcMarkSweepGC -XX:NewSize=300M -XX:MaxNewSize=300M -jar $contextmap_jar mapper -reads $reads -readlen $readlen -o $output_dir -rmap $rmap -bwtbuild $bwtbuild -bwtindex $bwtindex -genome $genome_dir -t 8
