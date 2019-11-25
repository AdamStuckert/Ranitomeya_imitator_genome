# *Ranitomeya imitator* assembly working notebook

This document is meant to be a working notebook of the insanity that is attempting this project. Much of this is experimental computational work, so I will document what we did in order to keep track of it for my own benefit.

### Supernova assembly
This took a while to run, but eventually completed. However, the assembly was small relative to what we expected and crappy. Scaffolding with ONT reads helped, but it was still crappy. Abondoned in favor of assemblies with PacBio data.

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

ONT only data was not good--it was shorter than expected and had poor gene coverage. Abondoned.

The assemblies with ONT data were larger (~0.5 GB) than the PacBio only assembly. However, after we Pilon polished these assemblies, we found that there was both a lower overall contiguity (Contig N50) and a marked increase in gene duplicates (from a BUSCO analysis to the tetrapod gene set).
