Aligning and mapping RNA seq reads from what were going to be a few different projects.

_Ranitomeya imitator_ first:

```bash
sbatch AlignmentReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/devseries/reads_from_enrique .fastq.gz
```

Since I have already indexed the genome, I'll submit a slightly modified script without the indexing step for the other species.

_R fantastica run:
```bash
# note submitted from: /mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/test
# I want to override the output from the sbatch header in the script as well
sbatch --output RNAseqReadCountFantastica.log ReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MultispeciesDevSeries/readfiles/fantastica_reads .fq.gz
```

_R variabilis_ run:
```bash
# note submitted from: /mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/test
# I want to override the output from the sbatch header in the script as well
sbatch --output RNAseqReadCountVariabilis.log ReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MultispeciesDevSeries/readfiles/variabilis_reads .fq.gz
```

Rename all the populations of imitator reads.

```bash
for i in $(ls *counts | grep -v "R_"); do mv -- "$i" "R_imitator_$i"; done

#### 
# Huallaga
for i in $(ls *H*.counts); do mv -- "$i" "${i/H/striped_}"; done

## Sauce
for i in $(ls *S*.counts); do mv -- "$i" "${i/S/banded_}"; done

## Tarapoto
for i in $(ls *T*.counts); do mv -- "$i" "${i/T/spotted_}"; done

## Varadero
for i in $(ls *V*.counts); do mv -- "$i" "${i/V/redheaded_}"; done

# change - to _
for i in $(ls R_imitator*); do mv -- "$i" "${i/-/_}"; done 
```

Get mapping data:

```bash
printf "Sample\tNumber_reads\tUniquely_mapping_reads\tMapped_to_multiple_loci\tMapped_to_too_many_loci\tUnmapped_reads_too_short\n" > mappingdata.tab
files=$(ls *Log.final.out)

for file in $files
do
sample=$(echo $file | sed "s/_Log.final.out//g" | sed "s/.Log.final.out//g")
reads=$(grep "Number of input reads" $file | cut -f 2)
unique=$(grep "Uniquely mapped reads %" $file | cut -f 2) 
multi=$(grep "% of reads mapped to multiple loci" $file | cut -f 2)
toomany=$(grep "% of reads mapped to too many loci" $file | cut -f 2)
unmapped=$(grep "reads unmapped: too short" $file | cut -f 2)
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$reads" "$unique" "$multi" "$toomany" "$unmapped" >> mappingdata.tab
done
```


Quick test of trimming reads to see how much better STAR does after trimming.

Imitator trimming script:
```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=24
#SBATCH --mem 110Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH --output trimmomatic.log

DIR=$(pwd)

mkdir $DIR/trimmed_reads
cd /mnt/lustre/macmaneslab/ams1236/devseries/reads_from_enrique

samples=$(ls *R1.fastq.gz | sed "s/.R1.fastq.gz//g")
for sample in $samples
do
trimmomatic PE -threads 24 \
-baseout $DIR/trimmed_reads/$sample.fastq.gz $sample.R1.fastq.gz $sample.R2.fastq.gz \
LEADING:3 TRAILING:3 ILLUMINACLIP:barcodes.fa:2:30:10 MINLEN:25
done
```

Fantastica reads:

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=24
#SBATCH --mem 110Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH --output trimmomatic.fant.log

DIR=$(pwd)

mkdir $DIR/trimmed_reads
cd /mnt/lustre/macmaneslab/ams1236/MultispeciesDevSeries/readfiles/fantastica_reads

cp $DIR/barcodes.fa .

samples=$(ls *R1.fq.gz | sed "s/_R1.fq.gz//g")
for sample in $samples
do
printf "Read 1: %s_R1.fq.gz" "$sample"
trimmomatic PE -threads 24 \
-baseout $DIR/trimmed_reads/$sample.fq.gz ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
LEADING:3 TRAILING:3 ILLUMINACLIP:barcodes.fa:2:30:10 MINLEN:25
done
```

Variabilis reads:
```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=24
#SBATCH --mem 110Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH --output trimmomatic.var.log

DIR=$(pwd)

mkdir $DIR/trimmed_reads
cd /mnt/lustre/macmaneslab/ams1236/MultispeciesDevSeries/readfiles/variabilis_reads

cp $DIR/barcodes.fa .

samples=$(ls *R1.fq.gz | sed "s/_R1.fq.gz//g")
for sample in $samples
do
printf "Read 1: %s_R1.fq.gz" "$sample"
trimmomatic PE -threads 24 \
-baseout $DIR/trimmed_reads/$sample.fq.gz ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
LEADING:3 TRAILING:3 ILLUMINACLIP:barcodes.fa:2:30:10 MINLEN:25
done
```

For this test, I will just use the reads for which there is a forward and reverse read. I'll rename those trimmed reads so my script works.

```bash
cd trimmed_reads
for i in $(ls *1P.fastq.gz); do mv -- "$i" "${i/1P/R1}"; done 
for i in $(ls *2P.fastq.gz); do mv -- "$i" "${i/2P/R2}"; done 
for i in $(ls *1P.fq.gz); do mv -- "$i" "${i/1P/R1}"; done  
for i in $(ls *2P.fq.gz); do mv -- "$i" "${i/2P/R2}"; done  

Imitator reads all end in `.fastq.gz`:

```bash
sbatch --output RNAseqReadCountTrimmedImitator.log ReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/trimmed_reads .fastq.gz
```

Variabilis/fantastica reads all end in `.fq.gz`:

```bash
sbatch --output RNAseqReadCountTrimmedModels.log ReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/trimmed_reads .fq.gz
```

Trimming helped decrease the proportion of reads tossed because they were considered too short by STAR. Now I want to see if a small adjustment in the length requirements will help recover some more reads (without hurting accuracy). 

Imitator reads all end in `.fastq.gz`:

```bash
sbatch --output RNAseqReadCountTrimmed_50percentlength_Imitator.log ReadCount_50percent_readlength.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/trimmed_reads .fastq.gz
```

Variabilis/fantastica reads all end in `.fq.gz`:

```bash
sbatch --output RNAseqReadCountTrimmed_50percentlength_Models.log ReadCount_50percent_readlength.job  \
$HOME/imitator_genome/imitator.1.3.6.fa \
$HOME/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
$HOME/MimicryGeneExpression/trimmed_reads .fq.gz
```

Tests concluded. Rerun with finalized version of the genome.

Imitator reads all end in `.fastq.gz`:

```bash
sbatch --output FINAL.RNAseqReadCountTrimmed_50percentlength_Imitator.log AlignmentReadCount_50percent_readlength.job  \
$HOME/imitator_genome/imitator.1.3.6.fa \
$HOME/imitator_genome/maker_1.3.6.masked_24June/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
$HOME/MimicryGeneExpression/trimmed_reads .fastq.gz
```

Variabilis/fantastica reads all end in `.fq.gz`:

```bash
sbatch --output FINAL.RNAseqReadCountTrimmed_50percentlength_Models.log ReadCount_50percent_readlength.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_24June/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/trimmed_reads .fq.gz
```
