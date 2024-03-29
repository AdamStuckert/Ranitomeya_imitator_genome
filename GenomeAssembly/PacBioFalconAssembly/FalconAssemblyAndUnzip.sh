#!/bin/bash
#SBATCH --job-name=pacbio
#SBATCH --output=pacbio.log
#SBATCH --cpus-per-task=24
#SBATCH --partition=macmanes,shared
# echo commands to stdout
set -x

DIR=$(pwd)

# load environments
module purge
module load linuxbrew/colsa

# create fasta files
cd ${DIR}/raw_PacBio_data/
PBdat=$(ls */*subreads.bam | sed "s/bam//g")
for data in PBdat
do
samtools fasta ${data}bam > ${data}fasta
done

cd ${DIR}

# change environments
module purge
module load anaconda/colsa
source activate pacbio-20190801


### create FOFN files
mkdir falcon_assembly
# fasta files first
touch falcon_assembly/PacBioFastaFiles.fofn
for fasta in $(ls raw_PacBio_data/*/*subreads.fasta)
do
(printf '%s/%s \n' "$DIR" "$fasta") >> falcon_assembly/PacBioFastaFiles.fofn
done

# bam files second
touch falcon_assembly/PacBioBamFiles.fofn
for fasta in $(ls raw_PacBio_data/*/*subreads.bam)
do
(printf '%s/%s \n' "$DIR" "$fasta") >> falcon_assembly/PacBioBamFiles.fofn
done

# run Falcon assembly
fc_run fc_run.cfg  &> FalconRun.log &
# follow along
tail -f FalconRun.log

# run Falcon unzip
mv all.log all0.log
fc_unzip.py fc_unzip.cfg &> FalconUnzip.std &


# Follow up steps:
# 1. Phase genomics
# 2. Gap fill with LRGap
# 3. Pilon polish via 10x data
