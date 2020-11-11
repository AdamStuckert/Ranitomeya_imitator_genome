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

