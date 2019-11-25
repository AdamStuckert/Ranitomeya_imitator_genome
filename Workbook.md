# *Ranitomeya imitator* assembly working notebook

This document is meant to be a working notebook of the insanity that is attempting this project. Much of this is experimental computational work, so I will document what we did in order to keep track of it for my own benefit.

### Supernova assembly
This took a while to run, but eventually completed. However, the assembly was small relative to what we expected (25% of tetrapod core genes) and crappy. Scaffolding with ONT reads helped, but it was still crappy. Abondoned in favor of assemblies with PacBio data.

### Falcon assembly
This seemed very promising. However, it produced > 16 TB of intermediate files and we hit quota. I aim to resume this eventually, if possible.

### Masurca assembly
This ran for about 35 days, was making almost no progress, and so I eventually gave up after reaching quota issues due to the Falcon assembly.

### wtdbg2 assembly
This is the assembler that was used to produce the axolotl genome. This 1. works and 2. is fast. Much of this documentation is for this approach, largely because we were able to actually assemble a genome from it.

We took 4 approaches to the initial assembly:

1. PacBio data with PacBio specs
2. Nanopore data with ONT specs
3. Assembly with PacBio and Nanopore data with PacBio specs
4. Assembly with PacBio and Nanopore data with ONT specs

ONT only data was not good--it was shorter than expected and had poor gene coverage. Abandoned without polishing.

The assemblies with ONT data were larger (~0.5 GB) than the PacBio only assembly. However, after we Pilon polished these assemblies, we found that there was both a lower overall contiguity (Contig N50) and a marked increase in gene duplicates (from a BUSCO analysis to the tetrapod gene set).

Proof from polished assemblies:

Assembly | Genome Size (GB) | Contig N50 | BUSCO 
--- | --- | --- | ---
Pacbio only | 6.77 | 198779 | C:92.3%[S:75.4%,D:16.9%],F:4.6%,M:3.1%,n:3950
Nanopore only | | | C:0.1%[S:0.1%,D:0.0%],F:1.7%,M:98.2%,n:3950 
All data, PacBio specs | 7.82 | 149205 | C:91.9%[S:72.6%,D:19.3%],F:4.3%,M:3.8%,n:3950
All data, ONT specs | 7.83 | 149048 | C:88.1%[S:69.4%,D:18.7%],F:7.3%,M:4.6%,n:3950

I then quickmerged the PacBio only assembly with the all data PacBio specs assembly. 

Code to merge:
```
merge_wrapper.py -pre imitator.1.0 \
$HOME/imitator_genome/imitator.alldata.pacbiospec.wtdbg.ctg.polished.fa \
$HOME/imitator_genome/imi_wtdbg.ctg.polished.fa
```

This assembly was probably not great (Contig N50 = 186687; total size = 8.28; BUSCO = C:92.5%[S:72.4%,D:20.1%],F:4.1%,M:3.4%,n:3950). In particular, I thought that the duplicated genes were inflated along with overall size. This may in fact be real, but it seems more so to be an artifact of assembly + merge. 




