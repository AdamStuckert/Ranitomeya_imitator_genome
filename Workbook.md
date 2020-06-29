# *Ranitomeya imitator* assembly working notebook

This document is meant to be a working notebook of the insanity that is attempting this project. Much of this is experimental computational work, so I will document what we did in order to keep track of it for my own benefit.

#### A note about genome assembly assessment
Genomes were assessed with a script (`genomemetrics.job`) which runs `assemblathon_stats.pl` (script from [the assemblathon2 GitHub rep](https://github.com/ucdavis-bioinformatics/assemblathon2-analysis) and BUSCO 3.0.2 against the tetrapod core gene set. 

#### A note about scripts used
Scripts used in this notebook/project are found in the scripts directory. All `*.job` scripts are slurm submission scripts (our HPC here at UNH uses slurm) and all `*.sh` scripts are just generic bash scripts that run quickly on the head node.

## Supernova assembly
This took a while to run, but eventually completed. However, the assembly was small relative to what we expected (25% of tetrapod core genes) and crappy. Scaffolding with ONT reads helped, but it was still crappy. Abondoned in favor of assemblies with PacBio data.

## Falcon assembly
This seemed very promising. However, it produced > 16 TB of intermediate files and we hit quota. I aim to resume this eventually, if possible.

## Masurca assembly
This ran for about 35 days, was making almost no progress, and so I eventually gave up after reaching quota issues due to the Falcon assembly.

## wtdbg2 assembly
This is the assembler that was used to produce the axolotl genome. This 1) works and 2) is fast. Much of this documentation is for this approach, largely because we were able to actually assemble a genome from it. The "good" `subreads.bam` from the sequencer had already been converted to fasta files using `samtools fasta` before this for our attempt using Falcon.

We took 4 approaches to the initial assembly:

1. PacBio data with PacBio specs
2. Nanopore data with ONT specs
3. Assembly with PacBio and Nanopore data with PacBio specs
4. Assembly with PacBio and Nanopore data with ONT specs

ONT only data was not good--it was shorter than expected and had poor gene coverage. Abandoned without polishing.

The assemblies with ONT data were larger (~0.5 GB) than the PacBio only assembly. However, after we Pilon polished these assemblies, we found that there was both a lower overall contiguity (Contig N50) and a marked increase in gene duplicates (from a BUSCO analysis to the tetrapod gene set).

Genome polishing (both these assemblies and subsequent iterations of this genome) is done in this fashion.

```bash
# first map Illumina (10x) reads to the genome
sbatch bwa.job $ASSEMBLY $PILON_DIRECTORY
# second, split up the genome into 80 chunks so pilon polishing doesn't take forever
chunks.sh $ASSEMBLY
# submit pilon polishing array
sbatch pilon.job $ASSEMBLY $ASSEMBLY_OUT # directory within --frags flag has to be $PILON_DIRECTORY
```

Proof from polished assemblies:

Assembly | Genome Size (GB) | Contig N50 | BUSCO 
--- | --- | --- | ---
Pacbio only (imi_wtdbg.ctg.polished.fa) | 6.77 | 198779 | C:92.3%[S:75.4%,D:16.9%],F:4.6%,M:3.1%,n:3950
Nanopore only | | | C:0.1%[S:0.1%,D:0.0%],F:1.7%,M:98.2%,n:3950 
All data, PacBio specs (imitator.alldata.pacbiospec.wtdbg.ctg.polished.fa) | 7.82 | 149205 | C:91.9%[S:72.6%,D:19.3%],F:4.3%,M:3.8%,n:3950
All data, ONT specs (imitator.alldata.nanoporespec.wtdbg.ctg.polished.fa) | 7.83 | 149048 | C:88.1%[S:69.4%,D:18.7%],F:7.3%,M:4.6%,n:3950

I then quickmerged the PacBio only assembly with the all data PacBio specs assembly. 

Code to merge:
```
merge_wrapper.py -pre imitator.1.0 \
$HOME/imitator_genome/imitator.alldata.pacbiospec.wtdbg.ctg.polished.fa \
$HOME/imitator_genome/imi_wtdbg.ctg.polished.fa
```

This assembly was probably not great.

Assembly | Genome Size (GB) | Contig N50 | BUSCO 
--- | --- | --- | ---
Merged assembly | 8.28 | 186687 | C:92.5%[S:72.4%,D:20.1%],F:4.1%,M:3.4%,n:3950

In particular, I thought that the duplicated genes were inflated along with overall size. This may in fact be real, but it seems more so to be an artifact of assembly + merge. 

## Improving the assembly via scaffolding, gap-filling, polishing, etc

I then elected to move forward with just the initial PacBio assembly that I Pilon polished. That seemed to be the better approach. I then scaffolded this assembly with 10x data using `arcs`. I used the provided `arcs.mk` file and ran the `arcs` processes (ie `arcs.mk arcs`). Submission script for arcs:

```bash
sbatch arcs.job imitator.1.0.fa imitator.1.1.fa
```
(I probably have to rerun longranger to recreate this barcoded file...)
Some notes here:
1. My attempt to automate changing the produced assembly name did not initially work (and needs to be fixed still). Error: `cp: target ‘/mnt/lustre/macmaneslab/ams1236/imitator_genome/arcs_run/imitator.1.1.fa’ is not a directory`
2. When I attempted to run my genome metrics script (`genomemetrics.job`, which runs an Assemblathon script + BUSCO), it turned out that the header names from the `arcs` program were too long. I fixed with:

```bash
awk -F',' '{print $1}' imitator.1.1.fa > tmp.fa
mv tmp.fa imitator.1.1.fa 
```

**This has been fixed, presumably. Leaving this here in case it does not work, so I have a starting point. I also added some cleanup code to remove files that aren't needed**

Genome metrics from assembly scaffolded with 10x data (which is imitator.1.1.fa):

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.1.fa | 6.78 | 198779 | 330635 | 0.02 | C:92.3%[S:74.9%,D:17.4%],F:4.4%,M:3.3%,n:3950

The next step is to run RAILS to fill gaps. We have PacBio data as well as ONT data that we can run RAILS with. RAILS requires sequencing platform details (as it changes the algorithm of minimap2), and is specified as `pb`, `ont`, or `nil`. It is unclear the best way to run this.

So we are taking an experimental approach to this.

1. PacBio data only (submission: `sbatch --output rails.pb.log rails.pb.job imitator.1.1.fa pb`)
2. ONT data only (submission: `sbatch --output rails.ont.log rails.ont.job imitator.1.1.fa ont`)
3. PacBio + ONT data, with `pb` specs (submission: `sbatch --output rails.alldat.pbspecs.log rails.job imitator.1.1.fa pb`)
4. PacBio + ONT data, with `ont` specs (submission: `sbatch --output rails.alldat.ontspecs.log rails.job imitator.1.1.fa ont`) # this one is the only one not currently running, other 4 currently are
5. PacBio + ONT data, with `nil` specs (submission: `sbatch --output rails.alldat.nilspecs.log rails.job imitator.1.1.fa nil`)

Who knows what is best??? Verify this with the assemblathon script and BUSCO.

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
PB data only | 6.78 | 201277 | 336960 | 0.02 | C:92.5%[S:75.1%,D:17.4%],F:4.4%,M:3.1%,n:3950 |
ONT data only | 6.79 | 201935 | 332101 | | C:92.3%[S:74.9%,D:17.4%],F:4.4%,M:3.3%,n:3950 | 
PB + ONT data, PB specs | 6.79 | 204382 | 339129 | 0.01 | C:92.5%[S:75.1%,D:17.4%],F:4.4%,M:3.1%,n:3950 |
PB + ONT data, ONT specs | 6.78 | 201935 | 332101 | 0.02 | C:92.3%[S:74.9%,D:17.4%],F:4.4%,M:3.3%,n:3950 |
PB + ONT data, NIL specs | 6.79 | 204382 | 339129 | 0.01 | C:92.5%[S:75.1%,D:17.4%],F:4.4%,M:3.1%,n:3950 |

Based off of this, I'm going to use the assembly that had all of the data with the "nil" specs.

Using this assembly, renaming it `imitator.1.2.fa`.

`reads/combined_PacBio_Nanopore_reads.fa_vs_imitator.1.1.fa_90_0.90_rails.scaffolds.fa`

**note: After this round of runs I modified the script to 1) accept an output assembly name and 2) remove large extraneous files (bam files, formatted data files)**


Following this: pilon polishing and gapfilling with LRgap.


**been repeatedly bitten by the damn bwa/pilon naming of contigs/scaffold issue. add a workaround to one script or another!**

Workaround, probably should go in `bwa.job` at the beginning of script:
```bash
awk '{print $1}' $ASSEMBLY > new.fasta
mv new.fasta $ASSEMBLY 
```

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.2.1.fa | 6.79 | 211576 | 339195 | 0.01 | C:92.7%[S:74.5%,D:18.2%],F:4.3%,M:3.0%,n:3950

Ran a round of gapfilling. First I filled with ONT data and then I fileed that assembly with PacBio data.

Code to submit these two attempts:
```bash
sbatch lrgap.ont.job imitator.1.2.1.fa
sbatch lrgap.pb.job imitator.1.2.2.fa
```

Gap-filled assembly metrics:

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.2.2 (ONT filled) | 6.79 | 247642 | 339195 | 0.01 | C:92.6%[S:74.5%,D:18.1%],F:4.3%,M:3.1%,n:3950
imitator.1.2.3 (PacBio filled) | 6.79 | 272070 | 339195 | 0.00 | C:92.7%[S:74.5%,D:18.2%],F:4.3%,M:3.0%,n:3950

Next I will polish this gap-filled assembly with Pilon using 10x Illumina reads.

**quick note to self: fix pilon.job and bwa.job so that the pilon directory is hard coded**

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.3 (1.2.3, pilon polished) | 6.79 | 275328 | 339195 | 0.00 | C:92.7%[S:74.0%,D:18.7%],F:4.3%,M:3.0%,n:3950

Scaffolded again with 10X data using arcs.

```bash
sbatch arcs.job imitator.1.3.fa imitator.1.3.1.fa
```

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.3.1 (10X scaffolded) | 6.79 | 275328 | 397353 | 0.01 | C:92.6%[S:73.6%,D:19.0%],F:4.3%,M:3.1%,n:3950-

Our cluster went down for routine maintenance/upgrade at this point. Everythings location went from `/mnt/lustre/` to `/mnt/home/`. I changed all the job submission scripts to represent this:

```bash
jobs=$(ls *job)
for job in $jobs
do
sed -i "s/mnt\/lustre/mnt\/home/g" $job
done
```

Anything that was a random script elsewhere also got updated in a similar fashion. Hopefully things will work with relative ease.

Scaffold + gapfill with cobbler/rails.

```bash
sbatch --output rails.alldat.nilspecs.log rails.job imitator.1.3.1.fa nil imitator.1.3.2.fa
```

New assembly metrics:

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.3.2 | 6.79 | 275704 | 397634 | 0.01 | C:92.6%[S:73.6%,D:19.0%],F:4.3%,M:3.1%,n:3950

Next we did a round of pilon polishing.

```bash
# first map Illumina (10x) reads to the genome
sbatch bwa.job imitator.1.3.2.fa
# second, split up the genome into 80 chunks so pilon polishing doesn't take forever
chunks.sh imitator.1.3.2.fa
# submit pilon polishing array
sbatch pilon.job imitator.1.3.2.fa imitator.1.3.3.fa
```

Now on to LR gapfilling with the Nanopore data first, followed by the PacBio data. 

```bash
sbatch lrgap.ont.job imitator.1.3.3.fa imitator.1.3.4.fa
sbatch lrgap.pb.job imitator.1.3.4.fa imitator.1.3.5.fa 
```

Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.3.4 | 6.79 | 292624 | 397633 | 0.01 | C:92.7%[S:73.6%,D:19.1%],F:4.3%,M:3.0%,n:3950
imitator.1.3.5 | 6.79 | 300673 | 397633 | 0.01 | C:92.7%[S:73.6%,D:19.1%],F:4.3%,M:3.0%,n:3950

Next we did a round of pilon polishing.

```bash
# first map Illumina (10x) reads to the genome
sbatch bwa.job imitator.1.3.5.fa
# second, split up the genome into 80 chunks so pilon polishing doesn't take forever
chunks.sh imitator.1.3.5.fa
# submit pilon polishing array
sbatch pilon.job imitator.1.3.5.fa imitator.1.3.6.fa
```


Assembly | Genome Size (GB) | Contig N50 | Scaffold N50 | %Ns | BUSCO 
--- | --- | --- | --- | --- | ---
imitator.1.3.6 | 6.79 | 301327 | 397629 | 0.01 | C:92.7%[S:73.6%,D:19.1%],F:4.3%,M:3.0%,n:3950


The arrival of our Hi-C data is imminent. Rather than continuing to incrementally eek out improvements, I will wait for this data. Then, I will use the Hi-C data to scaffold the current assembly (1.3.5), as well as scaffold the original assembly. The rationale behind this is that while we have definitely improved the quality of the assembly, we have also had an increase of ~3% of duplicated genes. I would like to make sure that this is a "real" thing and not an assembly artifact.

### Update from the pandemic 

I'm beginning to think that our Hi-C data will never materialize. Thanks coronavirus! We are moving forward without the Hi-C data. I will run Maker to annotate the genome, as well as run BUSCO against the new version of the database/software. I've done some preliminary tests that indicate that running Maker alone to annotate will have issues due to repeat regions. So I'm running Maker thrice: once without transcript evidence, once with transcript evidence, and finally once with transcript evidence after running Repeat Modeler and Repeat Masker.

Stand alone Maker code is currently in the Maker directory. 

Repeat Modeler code:

```#!/bin/bash
#SBATCH --job-name=repeatmod
#SBATCH --output=repeatmodeler.log
#SBATCH --partition=macmanes
#SBATCH --ntasks=40
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118


module load linuxbrew/colsa

# CPU=$(lscpu | grep CPU\(s\) | head -n1 | awk '{print $2}')

BuildDatabase -name imitator.1.3.6.repeatmodeler_db -engine ncbi imitator.1.3.6.fa

RepeatModeler -database imitator.1.3.6.repeatmodeler_db -pa 40
```

Repeat Masker code:

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J repeatmask
#SBATCH --output repeatmask_imi2.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

mkdir repeatmask_postModeler
cd repeatmask_postModeler
PATH=/mnt/lustre/macmaneslab/macmanes/ncbi-blast-2.7.1+/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/shared/augustus_config/config


# input requires two fasta files. 1) the masked fasta file from RepeatMasker; 2) the genome assembly

RepeatMasker -pa 40 -gff -lib /mnt/lustre/macmaneslab/ams1236/imitator_genome/consensi.fa.classified -q /mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa
```

I checked the annotations with a few quality metrics. Total number of contigs, total unique contigs, number of proteins of unknown function, BUSCO transcriptome score.

Annotation | Total transcripts | Unique transcripts | # unknown proteins  | BUSCO 
--- | --- | --- | --- | ---
No transcript evidence | 30803 | 16360 | 1610 |  C:77.7%[S:60.9%,D:16.8%],F:9.4%,M:12.9%,n:3950
Transcript evidence | 144683 | 108705 | 18051 | C:84.1%[S:65.0%,D:19.1%],F:9.3%,M:6.6%,n:3950
Transcript evidence + masked | 52336 | 26336 | 16929 | C:82.4%[S:64.1%,D:18.3%],F:10.8%,M:6.8%,n:3950

At this point I was able to get my hands on the newest version of the RepBase database (link to website here). This is version `RepBase25.05.2020`. So I downloaded the newest version of RepeatMasker, linked in the newest version of dfam (3.1) and then used the newest RepBase database to rerun RepeatMasker. I concatenated all the RepBase fastas (`cat *ref > RepBase25.05.fasta`) and specified that as the library in Repeat Masker.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J repeatmask
#SBATCH --output repeatmasker.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

DIR=$(pwd)
PATH=/mnt/lustre/macmaneslab/macmanes/ncbi-blast-2.7.1+/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/shared/augustus_config/config

# cp $masked .

mkdir repeatmask_adamsinstall_newRepBaseDB
cd repeatmask_adamsinstall_newRepBaseDB

which perl


$HOME/software/RepeatMasker/RepeatMasker -pa 40 -gff -lib /mnt/lustre/macmaneslab/ams1236/software/RepeatMasker/Libraries/RepBase25.05.fasta  -q ${DIR}/imitator.1.3.6.fa
```

This actually produced less informative masking than with just the dfam 3.1 database run on the Repeat Modeled imitator genome. So we will use the dfam 3.1/repeat modeled imitator genome. I then ran Maker on this assembly, and calculated some statistics.

Annotation | Total transcripts | Unique transcripts | # unknown proteins  | BUSCO 
--- | --- | --- | --- | ---
imitator.1.3.6 repeat modeled, masked, and with transcript evidence | 52325 | 24862 | 16926 |  C:82.5%[S:63.9%,D:18.6%],F:10.6%,M:6.9%,n:3950

For clarification, to caluclate "unique transcripts" I used anything that should be considered a unique gene. This means all unknown proteins count as 1, and two transcripts with a top hit to the same gene (eg A1BG Alpha-1B-glycoprotein (Homo sapiens OX=9606) and A1bg Alpha-1B-glycoprotein (Rattus norvegicus OX=10116) would be counted once in total). Code for that:

> cat annotations.txt | sed "s/(.*)//g" | tr '[:upper:]' '[:lower:]' | sort | uniq | wc -l
