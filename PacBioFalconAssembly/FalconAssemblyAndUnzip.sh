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
cd ${DIR}/raw_PacBio_data/3_C01/
bam2fastx -o m64019_190803_173458.subreads m64019_190803_173458.subreads.bam

cd ${DIR}

# change environments
module purge
module load anaconda/colsa
source activate pacbio-20190801


### create FOFN files
# fasta files first
touch PacBioFastaFiles.fofn
for fasta in $(ls raw_PacBio_data/*fasta)
do
(printf '%s/%s \n' "$DIR" "$fasta") >> PacBioFastaFiles.fofn
done

# bam files second
touch PacBioBamFiles.fofn
for fasta in $(ls raw_PacBio_data/*subreads.bam)
do
(printf '%s/%s \n' "$DIR" "$fasta") >> PacBioBamFiles.fofn
done

# run Falcon assembly
fc_run fc_run.cfg  &> FalconRun.log &

# run Falcon unzip
mv all.log all0.log
fc_unzip.py fc_unzip.cfg &> FalconUnzip.std &
