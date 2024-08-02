output_dir=`pwd`
parent_dir="$(dirname "$output_dir")"

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
elif [ "$aligner_name" = bowtie2 ] 
	then read -p "Please enter the path for the bowtie2 binary " aligner_bin
	read -p "Please enter the path for the bowtie2 indexer " indexer_bin
fi

#
# Simply modify the following paths to mine for contaminations and infections with your own data
#
contextmap_jar=$parent_dir/ContextMap_v2.7.9.jar
reads=$output_dir/reads.fa
readlen=100
aligner_index=$output_dir/aligner_index/demo_index
reference_sequences_basename=$output_dir/reference_sequences/genomes
indexer_output_dir=$output_dir/indexed_references
genomes_dir=$output_dir/indexed_references/fasta
indices_dir=$output_dir/indexed_references/indices


touch log.txt
true > log.txt
#
# 1. Indexing the given genome sequences with the 'indexer' tool
#
mkdir $output_dir/indexed_references >> log.txt 2>&1
unzip -o $output_dir/reference_sequences/genomes.zip -d $output_dir/reference_sequences >> log.txt 2>&1
echo "#"
echo "# Indexing given reference sequences with the 'indexer' tool"
echo "#"

java -jar $contextmap_jar indexer -fasta ${reference_sequences_basename}_0.fa -prefix microbesA -o $indexer_output_dir >> log.txt 2>&1
java -jar $contextmap_jar indexer -fasta ${reference_sequences_basename}_1.fa -prefix microbesB -o $indexer_output_dir >> log.txt 2>&1


echo ""

#
# 2. Building an index for the desired aligner on the newly generated fasta sequences with the original indexer
#

echo "#"
echo "# Building indices on the newly generated fasta sequences with the original indexer"
echo "#"
mkdir $output_dir/aligner_index/ >> log.txt 2>&1

if [ "$aligner_name" = bowtie1 -o "$aligner_name" = bowtie2 ]
then ${indexer_bin} -f $genomes_dir/microbesA_0.fa ${aligner_index}_A >> log.txt 2>&1
     ${indexer_bin} -f $genomes_dir/microbesB_0.fa ${aligner_index}_B >> log.txt 2>&1
else ${indexer_bin} index -p ${aligner_index}_A $genomes_dir/microbesA_0.fa >> log.txt 2>&1
     ${indexer_bin} index -p ${aligner_index}_B $genomes_dir/microbesB_0.fa >> log.txt 2>&1
fi


echo ""

#
# 3. Running ContextMap
#
echo "#"
echo "# Running ContextMap"
echo "#"
unzip -o $output_dir/reads.zip -d $output_dir >> log.txt 2>&1
java -Xms1000m -Xmx1000m -jar $contextmap_jar mapper --mining -aligner_name $aligner_name -aligner_bin $aligner_bin -indexer_bin $indexer_bin -reads $reads -o $output_dir/mapping -indices ${aligner_index}_A,${aligner_index}_B -genome $genomes_dir -skipsplit true,true -speciesindex $indices_dir -maxhits 50 >> log.txt 2>&1
echo ""

#
# 4. Extracting contained species in the generated SAM file with the 'inspector' tool
#
echo "#"
echo "# Extracting contained species in the generated SAM file with the 'inspector' tool"
echo "#"
java -Xms1000m -Xmx1000m -jar $contextmap_jar inspector -sam $output_dir/mapping/mapping.sam -idx $indices_dir > $output_dir/contained_species.txt 2>> log.txt
echo ""
echo "#"
echo "# Everything done. Results were printed to the 'contained_species.txt' file"
echo "#"
