## New assembly details

First, assemble with wtdbg and polish with racon

```sh
#!/usr/bin/sh
#SBATCH --mem=700GB
#SBATCH --job-name="imi_assembly"
#SBATCH --output="imi_assembly.log"
#SBATCH --partition=macmanes
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118


DIR=$(pwd)
READ1="${HOME}/ratshit/data/m64019_190803_173458.subreads.fa"
READ2="${HOME}/ratshit/data/m64019_190913_173453.subreads.fasta"
READ3="${HOME}/ratshit/data/m64019_190915_200441.subreads.fasta"
TENX1="${HOME}/ratshit/data/trimmed10x.R1.fastq.gz"
TENX2="${HOME}/ratshit/data/trimmed10x.R2.fastq.gz"
REF="imitator.wtdbg.ctg.fa"
REF_PATH="${DIR}/wtdbg/${REF}"


module purge
# load conda source script!
. ${HOME}/anaconda3/etc/profile.d/conda.sh


mkdir wtdbg
cd wtdbg

wtdbg2 \
-o imitator.wtdbg \
-x sq \
-g 6.7g \
-L 2500 \
-t 40 \
-i ${READ1} \
-i ${READ2} \
-i ${READ3}

wtpoa-cns -t 40 -i imitator.wtdbg.ctg.lay.gz -fo $REF

# rename fasta headers!
awk '{print $1}' $REF > new.fasta
mv new.fasta $REF

#### Calculate genome metrics
printf "Submitting genome metrics job\n\n\n"
sbatch $HOME/scripts/genomemetrics.job $REF $(echo $REF | sed "s/.fa//")

#### polish with Pilon
cat << EOF > bwa.job
#!/bin/bash
#SBATCH --partition=macmanes
#SBATCH -J bwa
#SBATCH --output bwa.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude node117,node118


#REF=$REF

species="Ranitomeya_imitator"
prefix="Rimi"


module purge
module load linuxbrew/colsa

echo Running BWA on $REF


mkdir pilon
cd pilon

bwa index ${REF_PATH}

rm -f "$prefix".sorted.bam
rm -f "$prefix".sorted.bam.bai

bwa mem -t 40 ${REF_PATH} ${Tenx1} ${Tenx2} \
| samtools view -@20 -Sb - \
| samtools sort -T "$prefix" -O bam -@20 -l9 -m2G -o "$prefix".sorted.bam -
samtools index "$prefix".sorted.bam


### Split assembly into chunks to speed up mapping

# make folder
mkdir chunks

# remove any chunks and list from previous iteration
rm chr.list
rm chunks/*

# genome headers
grep ">" ${REF_PATH} | sed 's_>__' | shuf | tee -a chr.list


# slit into 80 chunks
split -d -n l/80 chr.list chunks/genomechunk.

cd chunks/
rename genomechunk.0 genomechunk. genomechunk.0*
printf "Submitting pilon array job\n\n\n"
sbatch ../pilon.job
EOF

printf "Submitting bwa alignment and chunk splitting job\n\n\n"
sbatch bwa.job


#### polish with racon
cat << EOF > racon.job
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J racon
#SBATCH --output racon.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=700Gb


module load linuxbrew/colsa

# preparation
mkdir racon
cd racon

ln -s $REF_PATH

### First align reads with minimap2
echo aligning with minimap2
minimap2 -I10G -t 40 -xmap-pb $REF $READ1 $READ2 $READ3 | gzip -c - > Imi.PB.paf.gz

## merge reads into a single file for racon...
cat $READ1 $READ2 $READ3 > AllReads.fa

### Run racon
echo Polishing with racon
racon -t 40 AllReads.fa  Imi.PB.paf.gz $REF > ${REF}.raconpolished.fa

## genome metrics
sbatch $HOME/scripts/genomemetrics.job ${REF}.raconpolished.fa ${REF}.raconpolished

EOF

printf "Submitting pilon polishing job\n\n\n"
sbatch racon.job

```

Next, scaffold with prna_scaffolder


```sh
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J PRNAscaf
#SBATCH --cpus-per-task=40
#SBATCH --output PRNAscaf.log
#SBATCH --mem 300Gb
#SBATCH --exclude=node117,node118


# input files
ASSEMBLY="${HOME}/ratshit/imitator.pacbio.wtdbg.ctg.fa.raconpolished.fa"
READ1="$HOME/ratshit/data/RNA.R1.fq"
READ2="$HOME/ratshit/data/RNA.R2.fq"

mkdir prna
cd prna

# This uses only the RNA seq reads that went in to the transcriptome assembly

##### Maping with BWA
printf "#####################################################################\n"
printf "#####################################################################\n"
printf "################## Indexing genome with BWA #########################\n"
printf "#####################################################################\n"
printf "#####################################################################\n"

bwa index $ASSEMBLY

printf "#####################################################################\n"
printf "#####################################################################\n"
printf "################# Alighing RNA reads with BWA    ####################\n"
printf "#####################################################################\n"
printf "#####################################################################\n"

bwa mem -t 40 $ASSEMBLY $READ1 $READ2 > input.sam


printf "#####################################################################\n"
printf "#####################################################################\n"
printf "#################### Scaffolding with P_RNA #########################\n"
printf "#####################################################################\n"
printf "#####################################################################\n"

P_RNA_scaffolder.sh -d $HOME/software/P_RNA_scaffolder -i input.sam -j $ASSEMBLY -F $READ1 -R $READ2  -t 40

echo Donezos
```

Scaffold with ARCS

```
#!/bin/bash
#SBATCH --partition=macmanes
#SBATCH --ntasks=40
#### #SBATCH --mem 700Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH -J arcs
#SBATCH --output arcs.log


#### USAGE:
## scaffold an assembly ($1) with 10x data, output to a new scaffolded assembly ($2)


# longranger basic --id=imitator10x --fastqs=${HOME}/ratshit/data/10x --sample=Rimi --localcores=40 --localmem=290


source /mnt/lustre/macmaneslab/ams1236/.bash_profile


conda activate base
module load linuxbrew/colsa


# Local variables
INPUTFASTA="/mnt/lustre/macmaneslab/ams1236/ratshit/wtdbg/prna/imitator.wtdbg.ctg.raconpolished.prnascaf.fa"
OUTPUTFASTA="imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.fa"
READS="${HOME}/ratshit/imitator10x/outs/barcoded.fastq.gz"
READNAME=$(echo $READS | sed "s/.fastq.gz//g")
ASSNAME=$(basename $INPUTFASTA)
DIR=$(pwd)
NOSUFFIXASS=$(echo $ASSNAME | sed "s/.fasta//g" | sed "s/.fa//g" )

# copy assembly into arcs folder
mkdir arcs_run
# cp $INPUTFASTA $DIR/arcs_run/$ASSNAME


$HOME/software/arcs/Examples/arcs-make arcs draft=$DIR/arcs_run/$NOSUFFIXASS reads=$READNAME

# rename scaffolded assembly
NEWASS=$(ls -lhtr arcs_run/*fa | awk '{print $9}' | tail -n1)
cp $NEWASS $DIR/arcs_run/$OUTPUTFASTA

# Rename headers so they aren't too long for BUSCO, etc
awk -F',' '{print $1}' $OUTPUTFASTA > tmp.fa
mv tmp.fa $OUTPUTFASTA

# cleanup
rm arcs_run/$INPUTFASTA
rm arcs_run/${NOSUFFIXASS}.renamed*
rm arcs_run/${NOSUFFIXASS}.sorted.bam
```

PRNA
```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J PRNAscaf
#SBATCH --cpus-per-task=40
#SBATCH --output PRNAscaf.log
#SBATCH --mem 300Gb
#SBATCH --exclude=node117,node118


# input files
ASSEMBLY="${HOME}/ratshit/imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.fa"
READ1="$HOME/ratshit/data/RNA.R1.fq"
READ2="$HOME/ratshit/data/RNA.R2.fq"

mkdir prna
cd prna

# This uses only the RNA seq reads that went in to the transcriptome assembly

##### Maping with BWA
printf "#####################################################################\n"
printf "#####################################################################\n"
printf "################## Indexing genome with BWA #########################\n"
printf "#####################################################################\n"
printf "#####################################################################\n"

bwa index $ASSEMBLY

printf "#####################################################################\n"
printf "#####################################################################\n"
printf "################# Alighing RNA reads with BWA    ####################\n"
printf "#####################################################################\n"
printf "#####################################################################\n"

bwa mem -t 40 $ASSEMBLY $READ1 $READ2 > input.sam


printf "#####################################################################\n"
printf "#####################################################################\n"
printf "#################### Scaffolding with P_RNA #########################\n"
printf "#####################################################################\n"
printf "#####################################################################\n"

P_RNA_scaffolder.sh -d $HOME/software/P_RNA_scaffolder -i input.sam -j $ASSEMBLY -F $READ1 -R $READ2  -t 40

echo Donezos
```

Racon 

```
REF=${HOME}/ratshit/imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.fa
DIR=$(pwd)
READ1="${HOME}/ratshit/data/m64019_190803_173458.subreads.fa"
READ2="${HOME}/ratshit/data/m64019_190913_173453.subreads.fasta"
READ3="${HOME}/ratshit/data/m64019_190915_200441.subreads.fasta"
cat << EOF > racon.job
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J racon
#SBATCH --output racon.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=700Gb


module load linuxbrew/colsa

# preparation
mkdir racon
cd racon

ln -s $REF_PATH

### First align reads with minimap2
echo aligning with minimap2
minimap2 -I10G -t 40 -xmap-pb $REF $READ1 $READ2 $READ3 | gzip -c - > Imi.PB.paf.gz

## merge reads into a single file for racon...
cat $READ1 $READ2 $READ3 > AllReads.fa

### Run racon
echo Polishing with racon
racon -t 40 AllReads.fa  Imi.PB.paf.gz $REF > ${REF}.raconpolished.fa

## genome metrics
sbatch $HOME/scripts/genomemetrics.job ${REF}.raconpolished.fa ${REF}.raconpolished

EOF

printf "Submitting pilon polishing job\n\n\n"
sbatch racon.job
```

Scaffold with RAILS

```
#!/bin/bash

#SBATCH --partition=macmanes
#SBATCH -J rails
#SBATCH --cpus-per-task=40
#SBATCH --output rails.log
#SBATCH --mem 700Gb
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

# set variables
ASSEMBLY=$HOME/ratshit/genome_versions/imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.fa
ASSNAME=$(basename $ASSEMBLY)
READ_PARAMETERS="nil"
ASSEMBLY_OUT=imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.fa
READS="$HOME/ratshit/data/combined_PacBio_Nanopore_reads.fa"
READNAME=$(basename $READS)
DIR=$(pwd)


##### prep raw data:
# convert ONT to fasta
seqtk seq -a $HOME/ratshit/data/allcombinedNanoporedataJune2019.fastq > $HOME/ratshit/data/allcombinedNanoporedataJune2019.fasta
cat $HOME/ratshit/data/allcombinedNanoporedataJune2019.fasta $HOME/ratshit/wtdbg/racon/AllReads.fa > $HOME/ratshit/data/combined_PacBio_Nanopore_reads.fa

mkdir rails
cd rails

# make sure the symlink doesn't exist...
rm $ASSNAME
ln -s $ASSEMBLY $ASSNAME
ln -s $READS ${DIR}/rails/$READNAME

runRAILSminimap.sh \
$ASSNAME \
$READNAME \
90 0.90 90 1 $READ_PARAMETERS \
~/linuxbrew/.linuxbrew/bin/samtools \
40 # threads

# copy assembly into new name
cp *rails.scaffolds.fa $ASSEMBLY_OUT

# clean up big files
rm *bam
rm *formatted.fa
```

ARCS
```
#!/bin/bash
#SBATCH --partition=macmanes
#SBATCH --ntasks=40
#### #SBATCH --mem 700Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH -J arcs
#SBATCH --output arcs.log


#### USAGE:
## scaffold an assembly ($1) with 10x data, output to a new scaffolded assembly ($2)


# longranger basic --id=imitator10x --fastqs=${HOME}/ratshit/data/10x --sample=Rimi --localcores=40 --localmem=290


source /mnt/lustre/macmaneslab/ams1236/.bash_profile


conda activate base
module load linuxbrew/colsa


# Local variables
INPUTFASTA="/mnt/lustre/macmaneslab/ams1236/ratshit/rails/imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.fa"
OUTPUTFASTA="imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.fa"
READS="${HOME}/ratshit/imitator10x/outs/barcoded.fastq.gz"
READNAME=$(echo $READS | sed "s/.fastq.gz//g")
ASSNAME=$(basename $INPUTFASTA)
DIR=$(pwd)
NOSUFFIXASS=$(echo $ASSNAME | sed "s/.fasta//g" | sed "s/.fa//g" )

# copy assembly into arcs folder
mkdir arcs_run
cp $INPUTFASTA arcs_run/


$HOME/software/arcs/Examples/arcs-make arcs draft=$DIR/arcs_run/$NOSUFFIXASS reads=$READNAME

# rename scaffolded assembly
NEWASS=$(ls -lhtr arcs_run/*fa | awk '{print $9}' | tail -n1)
cp $NEWASS $DIR/arcs_run/$OUTPUTFASTA

# Rename headers so they aren't too long for BUSCO, etc
awk -F',' '{print $1}' $OUTPUTFASTA > tmp.fa
mv tmp.fa $OUTPUTFASTA

# cleanup
rm arcs_run/$INPUTFASTA
rm arcs_run/${NOSUFFIXASS}.renamed*
rm arcs_run/${NOSUFFIXASS}.sorted.bam
```

RACON
```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J racon
#SBATCH --output racon.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=700Gb


module load linuxbrew/colsa

DIR=$(pwd)
# preparation
mkdir racon
cd racon

ln -s ${DIR}/imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.fa

### First align reads with minimap2
echo aligning with minimap2
minimap2 -I10G -t 40 -xmap-pb imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.fa /mnt/lustre/macmaneslab/ams1236/ratshit/data/m64019_190803_173458.subreads.fa /mnt/lustre/macmaneslab/ams1236/ratshit/data/m64019_190913_173453.subreads.fasta /mnt/lustre/macmaneslab/ams1236/ratshit/data/m64019_190915_200441.subreads.fasta | gzip -c - > Imi.PB.paf.gz

## merge reads into a single file for racon...
#cat /mnt/lustre/macmaneslab/ams1236/ratshit/data/m64019_190803_173458.subreads.fa /mnt/lustre/macmaneslab/ams1236/ratshit/data/m64019_190913_173453.subreads.fasta /mnt/lustre/macmaneslab/ams1236/ratshit/data/m64019_190915_200441.subreads.fasta > AllReads.fa

### Run racon
echo Polishing with racon
racon -t 40 ${HOME}/ratshit/wtdbg/racon/AllReads.fa  Imi.PB.paf.gz imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.fa > imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.raconpolished.fa

## genome metrics
sbatch /mnt/lustre/macmaneslab/ams1236/scripts/genomemetrics.job imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.raconpolished.fa imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.raconpolished
```

symlink assembly

```bash
cd ${HOME}/ratshit/
ln -s /mnt/lustre/macmaneslab/ams1236/ratshit/rails/arcs_run/racon/imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.raconpolished.fa .
/mnt/lustre/macmaneslab/ams1236/ratshit/rails/arcs_run/racon/imitator.wtdbg.ctg.raconpolished.prnascaf.arcsscaf.prnascaf.raconpolished.rails.arcsscaf.raconpolished.fa
```
