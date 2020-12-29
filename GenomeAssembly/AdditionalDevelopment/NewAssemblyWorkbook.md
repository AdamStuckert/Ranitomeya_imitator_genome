## Starting fresh to eliminate duplicates!

Trying a revised version of the assembly using `wtdbg2`.

```bash
#!/bin/bash
#SBATCH --partition=shared,macmanes
#SBATCH -J imigenome
#SBATCH --output imi.axolotlparameters.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=700000

module load linuxbrew/colsa

# convert subreads.bam to fasta file
# samtools fasta raw_data/1_A01/m64019_200422_015111.subreads.bam > raw_data/S_parvus_smrtcell_1.fasta

DIR=$(pwd)

module purge

mkdir wtdbg_imitator_axolotlparameters
cd $HOME/imitator_genome/wtdbg_imitator_axolotlparameters


# run assembler (genome size is a total guess here, but placed on the high end)
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtdbg2 \
-x sq \
-o imitator_axolotlparameters \
-g 6.8g \
-L 5000 \
-p 21 \
-S 2 \
--aln-noskip \
--rescue-low-cov-edges \
--tidy-reads 2500 \
-i ${DIR}/raw_PacBio_data/m64019_190918_221316.subreads.fa \
-i ${DIR}/raw_PacBio_data/m64019_190803_173458.subreads.fa \
-i ${DIR}/raw_PacBio_data/m64019_190912_004210.subreads.fa

# run consensus
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns \
-t 40 \
-i imitator_axolotlparameters.ctg.lay.gz \
-fo imitator_axolotlparameters.ctg.fa
```

### Results, compared to previous attempts

Assembly | Genome Size (GB) | Contig N50 | Number of contigs | BUSCO 
--- | --- | --- | --- | ---
Initial wtdbg2 assembly, polished (imi_wtdbg.ctg.polished.fa) | 6.77 | 198779 | *add* | C:92.3%[S:75.4%,D:16.9%],F:4.6%,M:3.1%,n:3950
imitator.1.3.6 (final version in bioRxiv submission) | 6.79 | 301327 | *add* | C:92.7%[S:73.6%,D:19.1%],F:4.3%,M:3.0%,n:3950
imitator_axolotlparameters.ctg.fa (not polished) | 7.1 | 175,973 | 92,773 | C:90.6%[S:88.1%,D:2.5%],F:4.8%,M:4.6%,n:3950

This led to a pretty dramatic decrease in the overall number of duplicated BUSCO orthologs, even though it led to slightly higher fragmented and "missing" gene content. Some polishing might improve this. Guess I need to pursue this more. Polishing (with pilon and racon) and Hi-C scaffolding are the next steps to test things out.

### Polishing:

I tried a round of Pilon polishing.

```bash
# first map Illumina (10x) reads to the genome
sbatch bwa.job imitator_axolotlparameters.ctg.fa
# second, split up the genome into 80 chunks so pilon polishing doesn't take forever
./chunks.sh imitator_axolotlparameters.ctg.fa
# submit pilon polishing array
sbatch --dependency=afterok:3830 pilon.job imitator_axolotlparameters.ctg.fa imitator_axolotlparameters.ctg.pilonpolished.fa
```

I also tried a round of Racon polishing.

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J racon
#SBATCH --output racon.log
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node117,node118
#SBATCH --mem=700000
set -x

DIR=$(pwd)
ASSEMBLY="$HOME/imitator_genome/imitator_axolotlparameters.ctg.fa"
genome=$(basename $ASSEMBLY)
READS=$"$HOME/imitator_genome/raw_PacBio_data/R_imitator_PacBio.fa"


module load linuxbrew/colsa

# preparation
mkdir racon
cd racon

cp $ASSEMBLY .

awk '{print $1}' $genome > new.fasta
mv new.fasta $genome

### First align reads with minimap2
echo aligning with minimap2
minimap2 -I10G -t 40 -xmap-pb $genome $READS | gzip -c - > Rimi.PB.paf.gz

### Run racon
echo Polishing with racon
racon -t 40 $READS Rimi.PB.paf.gz $genome > imitator_axolotlparameters.ctg.raconpolished.fa
```


Results:

Assembly | Genome Size (GB) | Contig N50 | Number of contigs | BUSCO 
--- | --- | --- | --- | ---
imitator_axolotlparameters.ctg.fa (not polished) | 7.1 | 175,973 | 92,773 | C:90.6%[S:88.1%,D:2.5%],F:4.8%,M:4.6%,n:3950
imitator_axolotlparameters.ctg.pilonpolished.fa | 7.1 | 176,010 | 92,773 | C:92.4%[S:79.4%,D:13.0%],F:4.2%,M:3.4%,n:3950
imitator_axolotlparameters.ctg.raconpolished.fa | 7.1 | 179,075 | 85,427 | C:91.8%[S:83.8%,D:8.0%],F:4.5%,M:3.7%,n:3950

Alright, so this Racon polished assembly looks pretty good! Next I'll try purging duplicates...

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J purge_dups
#SBATCH --output purge_dups_slurm.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=300G
set -x

DIR=$(pwd)
ASSEMBLY="$HOME/imitator_genome/racon/imitator_axolotlparameters.ctg.raconpolished.fa"
genome=$(basename $ASSEMBLY)
READS=$"$HOME/imitator_genome/reads/PacBio_reads.fa"


module load linuxbrew/colsa

ln -s $ASSEMBLY

### run pipeline step by step
echo aligning with minimap2
minimap2 -I50G -t 40 -xmap-pb $genome $READS | gzip -c - > RimiPB.paf.gz

echo initial calculations
pbcstat RimiPB.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

echo Splitting assembly
split_fa $genome > $assembly.split
minimap2 -I50G -t 40 -xasm5 -DP $assembly.split $assembly.split | gzip -c - > $assembly.split.self.paf.gz

echo Purging haplotigs with overlaps
purge_dups -2 -T cutoffs -c PB.base.cov $assembly.split.self.paf.gz > dups.bed 2> purge_dups.log

echo Getting purged sequences...
get_seqs dups.bed $genome
```

## Results

Assembly | Genome Size (GB) | Contig N50 | Number of contigs | BUSCO 
--- | --- | --- | --- | ---
imitator_axolotlparameters.ctg.fa (not polished) | 7.1 | 175,973 | 92,773 | C:90.6%[S:88.1%,D:2.5%],F:4.8%,M:4.6%,n:3950
imitator_axolotlparameters.ctg.pilonpolished.fa | 7.1 | 176,010 | 92,773 | C:92.4%[S:79.4%,D:13.0%],F:4.2%,M:3.4%,n:3950
imitator_axolotlparameters.ctg.raconpolished.fa | 7.1 | 179,075 | 85,427 | C:91.8%[S:83.8%,D:8.0%],F:4.5%,M:3.7%,n:3950
imitator_axolotlparameters.ctg.raconpolished.purged.fa | 6.96 | 185,798 | 76,258 | C:91.8%[S:83.9%,D:7.9%],F:4.5%,M:3.7%,n:3950


I also did a second round of racon polishing. Just out of curiousity. And ran pilon by just polishing bases and ignoring gaps. These were run separately.

Assembly | Genome Size (GB) | Contig N50 | Number of contigs | BUSCO 
--- | --- | --- | --- | ---
imitator_axolotlparameters.ctg.fa (not polished) | 7.1 | 175,973 | 92,773 | C:90.6%[S:88.1%,D:2.5%],F:4.8%,M:4.6%,n:3950
imitator_axolotlparameters.ctg.pilonpolished.fa | 7.1 | 176,010 | 92,773 | C:92.4%[S:79.4%,D:13.0%],F:4.2%,M:3.4%,n:3950
imitator_axolotlparameters.ctg.raconpolished.fa | 7.1 | 179,075 | 85,427 | C:91.8%[S:83.8%,D:8.0%],F:4.5%,M:3.7%,n:3950
imitator_axolotlparameters.ctg.raconpolished2x.fa | 7.1 | 180,375 | 83,067 | C:91.9%[S:81.8%,D:10.1%],F:4.6%,M:3.5%,n:3950

Given no real improvement in the 2x racon polished assembly, I ran arcs to scaffold the 1x racon polished assembly:

Assembly | Genome Size (GB) | Contig N50 | Number of contigs | Scaffold N50 | Number of Scaffolds | BUSCO 
--- | --- | --- | --- | ---
imitator_axolotlparameters.ctg.fa (not polished) | 7.1 | 175,973 | 92,773 | NA | NA | C:90.6%[S:88.1%,D:2.5%],F:4.8%,M:4.6%,n:3950
imitator_axolotlparameters.ctg.raconpolished.fa | 7.1 | 179,075 | 85,427 | NA | NA | C:91.8%[S:83.8%,D:8.0%],F:4.5%,M:3.7%,n:3950
imitator_axolotlparameters.ctg.raconpolished.arcsscaff.fa | 7.1 | 179,075 | 85,427 | 303,634 | 74312 | C:91.9%[S:82.1%,D:9.8%],F:2.4%,M:5.7%,n:5310
