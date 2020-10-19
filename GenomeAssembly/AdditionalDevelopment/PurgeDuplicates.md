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
```

After running this I am left with an assembly that has been purged of duplicates, which is called `purged.fa`. I then check the quality with my genome metrics assessment script.

> sbatch ~/scripts/genomemetrics.job purged.fa purged.genome

Results:

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.3.6 | 6.79 | 301,327 | 397,629 | 0.01 | C:92.7%[S:73.6%,D:19.1%],F:4.3%,M:3.0%,n:3950
purged.fa | 6.64 | 318,047 | 416,770 | 0.01 | C:92.6%[S:73.5%,D:19.1%],F:4.3%,M:3.1%,n:3950

The purged genome assembly has 65,005 contigs in 60,563 scaffolds. As you can see however, basically no improvement in genic content, particularly duplicated orthologs.

I then ran SALSA on this purged genome. This purged, SALSA'd genome had 65,072 contigs in 59,843 scaffolds.

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.3.6 | 6.79 | 301,327 | 397,629 | 0.01 | C:92.7%[S:73.6%,D:19.1%],F:4.3%,M:3.0%,n:3950
purged.fa | 6.64 | 318,047 | 416,770 | 0.01 | C:92.6%[S:73.5%,D:19.1%],F:4.3%,M:3.1%,n:3950
purged.SALSAd | 6.64 | 315,462 | 455,728 | 0.01 | C:92.7%[S:73.6%,D:19.1%],F:4.3%,M:3.0%,n:3950

really not much improvement any which way.

So, I went back to basically the drawing board and ran `purge_dups` on the initial assembly:


Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | --- 
imitator.1.0.fa | 6.77 | 198,779 | NA | 0.00% | C:92.3%[S:75.4%,D:16.9%],F:4.6%,M:3.1%,n:3950
imitator.1.0.purged.fa | 6.63 | 207,656 | NA | 0.00% | pending

The histogram looks like this:
![alt text](SupplementalFiles/PB.stat.png?raw=true "Why are you hovering over a dang histogram?")

The automated pipeline used these as cutoffs:

> 5       -8      36      37      57      141

These refer to....????
