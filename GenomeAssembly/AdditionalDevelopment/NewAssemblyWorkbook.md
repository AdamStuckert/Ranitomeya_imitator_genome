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

