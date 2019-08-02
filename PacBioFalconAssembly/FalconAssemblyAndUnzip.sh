#/bin/bash

module purge
module load anaconda/colsa
source activate pacbio-20190801

DIR=$(pwd)

### create FOFN files
# fasta files first
touch PacBioFastaFiles.fofn
for fasta in $(ls raw_PacBio_data/*fasta)
do
(printf '%s/%s \n' "$DIR" "$fasta") >> PacBioFastaFiles.fofn
done

# bam files second
touch PacBioBamFiles.fofn
for fasta in $(ls raw_PacBio_data/*bam)
do
(printf '%s/%s \n' "$DIR" "$fasta") >> PacBioBamFiles.fofn
done

# run Falcon assembly
fc_run fc_run.cfg  &> FalconRun.log &

# run Falcon unzip
mv all.log all0.log
fc_unzip.py fc_unzip.cfg &> FalconUnzip.std &
