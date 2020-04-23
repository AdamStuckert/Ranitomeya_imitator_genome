#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=24
#SBATCH --mem 110Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH --output star-htseq.log


#### Purpose: align reads to a genome and count them
### Inputs: path the genome fasta, path to genome gff file, 
## Note: if files are not gzipped the "--readFilesCommand zcat" flag in STAR needs to be changed.

DIR=$(pwd)
#GENOME=$1
#GFF=$2
#READ_DIR=$3
#SUFFIX=$4

# parse input

while getopts g:f:d:s: flag
do
    case "${flag}" in
        -g) GENOME=${OPTARG};;
        -f) GFF=${OPTARG};;
        -d) READ_DIR=${OPTARG};;
        -s) SUFFIX=${OPTARG};; 
    esac
done 
echo "Genome: $GENOME"
echo "gff: $GFF"
echo "Reads directory: $READ_DIR"
echo "Reads suffix: $SUFFIX"

# make a directory for star genome index
mkdir $DIR/STAR_index

# make a directory for star read mapping
mkdir $DIR/STAR_mapped

# Index genome with STAR
STAR --runMode genomeGenerate \
--genomeDir $DIR/STAR_index \
--genomeFastaFiles $GENOME \
--runThreadN 24 \
--limitGenomeGenerateRAM 70744733397 \
--sjdbOverhang 149 \
--genomeChrBinNbits 18 \
--sjdbGTFfile $GFF \
--sjdbGTFtagExonParentTranscript Parent

#####################
##### map reads #####
#####################

# list all samples to map
cd $READ_DIR
samples=$(basename -a $SUFFIX | sed 's/R"$SUFFIX"//g") # this lists each sample identifier without the filename suffix.

echo Mapping these samples: $samples

# map reads with STAR
for sample in $samples
do
STAR --runMode alignReads \
--genomeDir $DIR/STAR_index \
--readFilesIn ${READ_DIR}/${sample}R1${SUFFIX} ${READ_DIR}/${sample}R2${SUFFIX} \
--readFilesCommand zcat \
--runThreadN 16 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--quantMode GeneCounts \
--outFileNamePrefix ${DIR}/STAR_mapped/${sample} \
--outSAMtype BAM Unsorted SortedByCoordinate
done


#######################
##### count reads #####
#######################


# count reads with htseq-count
echo Counting genes with htseq-count for these samples: $samples

for sample in $samples
do
htseq-count -s no -f bam -t exon -i Parent $DIR/STAR_mapped/${sample}Aligned.sortedByCoord.out.bam $GFF > $DIR/STAR_mapped/${sample}.gene.counts
done

