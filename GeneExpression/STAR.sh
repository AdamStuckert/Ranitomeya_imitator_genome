#!/bin/bash STAR.sh
# USAGE: sh STAR.sh -a GENOME -g GFF -r READ_DIR -s SUFFIX
# USAGE: sh /project/stuckert/users/Stuckert/scripts/STAR50.sh -a /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/R_imi_1.0.fa -g /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/Maker_round3cat/R_imi.UNK.AED_0.7.functional.gff3 -r /project/stuckert/users/Stuckert/MultiSpeciesDevSeries/imitator_reads -s .fastq.gz

while getopts a:g:r:s: option
do
case "${option}"
in
a) ASSEMBLY=${OPTARG};;
g) GFF=${OPTARG};;
r) READ_DIR=${OPTARG};;
s) SUFFIX=${OPTARG};;
esac
done


DIR=$(pwd)
GENOME=$(echo $ASSEMBLY)
GENOME_NAME=$(basename $GENOME | sed "s/.fasta//" | sed "s/.fa//")
mkdir ${GENOME_NAME}_index
mkdir ${DIR}/STAR_mapped50

printf "Assembly is $ASSEMBLY\n"
printf "GFF is $GFF\n"
printf "Reads are in $READ_DIR\n"
printf "Read suffixes are $SUFFIX\n"

if [ -f ${DIR}/${GENOME_NAME}_index/SAindex ]
then
  printf "Genome index built\n\n"
else
  printf "Building genome index\n\n"
  
  cat << EOF > STAR_idx.job
#!/bin/bash
#SBATCH -p general
#sbatch -J STAR_idx
#SBATCH --cpus-per-task=48
#SBATCH -t 5-0
#SBATCH --mem=180G
#SBATCH -o STAR_idx_%j.out

    STAR --runMode genomeGenerate \
    --genomeDir $DIR/${GENOME_NAME}_index \
    --genomeFastaFiles $GENOME \
    --runThreadN 48 \
    --limitGenomeGenerateRAM 70744733397 \
    --sjdbOverhang 150 \
    --sjdbGTFfile $GFF
    
EOF

sbatch STAR_idx.job  | tee STAR_idx.job.txt

  # get dependencies
  IDXJOBID=$(tail -n1 STAR_idx.job.txt | cut -f 4 -d " ") 
  #rm STAR_idx.job.txt
  MAPJOBIDDEP=$(printf "#SBATCH --dependency=afterok:$IDXJOBID")
  echo $MAPJOBIDDEP
fi


#####################
##### map reads #####
#####################
# list all samples to map
samples=$(basename -a ${READ_DIR}/*1$SUFFIX | sed "s/1${SUFFIX}//g")

# for loop to submit jobs
for sample in $samples
do
    if [ -f ${DIR}/STAR_mapped50/${sample}_gene.counts ]
    then
        printf "$sample has been aligned and counted, moving to next sample\n"
    else
        printf "Aligning $sample\n"

    cat << EOF > ${sample}.mapping.job
#!/bin/bash
#SBATCH -p general
#sbatch -J $sample.map
#SBATCH --cpus-per-task=16
#SBATCH -t 3-0
#SBATCH --mem=100G
#SBATCH -o STAR_${sample}_%j.out
$MAPJOBIDDEP


        printf "Mapping $sample with STAR\n"

        # map with star
        STAR --runMode alignReads \
        --genomeDir ${DIR}/${GENOME_NAME}_index \
        --readFilesIn ${READ_DIR}/${sample}1${SUFFIX} ${READ_DIR}/${sample}2${SUFFIX} \
        --readFilesCommand zcat \
        --runThreadN 16 \
        --outFilterMultimapNmax 5 \
        --outFilterMismatchNmax 10 \
        --outFilterScoreMinOverLread 0.5 \
        --outFilterMatchNminOverLread 0.5 \
        --quantMode GeneCounts \
        --outFileNamePrefix ${DIR}/STAR_mapped50/${sample} \
        --outSAMtype BAM Unsorted

        #######################
        ##### count reads #####
        #######################
        # count reads with htseq-count
        # load conda source script!
        . /project/stuckert/software/anaconda3/etc/profile.d/conda.sh
        module purge
        conda activate gene_expression

        echo Counting $sample
        # run htseq-count
        htseq-count -s no -f bam -t exon -i Parent $DIR/STAR_mapped50/${sample}Aligned.out.bam $GFF > $DIR/STAR_mapped50/${sample}_gene.counts
        # rm bam files to free up disk space
        # rm $DIR/STAR_mapped/${sample}Aligned.out.bam
EOF
sbatch ${sample}.mapping.job

fi
done
