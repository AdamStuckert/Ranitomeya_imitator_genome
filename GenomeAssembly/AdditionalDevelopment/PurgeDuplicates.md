# Additional genome development

This documents additional ongoing development of the genome. I'm especially interested in 1) attempting to remove regions of uncollapsed heterozygosity and 2) scaffold this further if possible.

## Regions of heterozygosity

I am going to attempt to use both the programs `purge_dups` and `purge_haplotigs` to do this. The latter is a bit more annoying, as it keeps getting hung and requires X11 support and both of those are, frankly, annoying as shit.

### Purge dups run

The first program is [purge dups](https://github.com/dfguan/purge_dups).

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J purge_dups
#SBATCH --output purge_dups_slurm.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
set -x

DIR=$(pwd)
ASSEMBLY="$HOME/imitator_genome/imitator.1.3.6.fa"
genome=$(basename $ASSEMBLY)
READS=$"$HOME/imitator_genome/reads/PacBio_reads.fa"


module load linuxbrew/colsa

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

### UNKNOWN FROM GIT: Step 4. Merge hap.fa and $hap_asm and redo the above steps to get a decent haplotig set.
## I have to figure out precisely what this means.
```

After running this I am left with an assembly that has been purged of duplicates, which is called `purged.fa`. I then check the quality with my genome metrics assessment script.

> sbatch ~/scripts/genomemetrics.job purged.fa purged.genome

