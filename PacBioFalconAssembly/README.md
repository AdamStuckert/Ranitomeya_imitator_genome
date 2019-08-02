# Kind of a placeholder for my notes on how to do a Falcon assembly of  _Ranitomeya imitator_

#############################################################################
####################### FALCON DE NOVO ASSEMBLY #############################
#############################################################################

### INFORMATION FROM: https://github.com/PacificBiosciences/pb-assembly

# Falcon uses only a single input, a config file. Up to date info on config is here: https://github.com/PacificBiosciences/FALCON-integrate/wiki/Configuring-Unzip

# Assemble
```bash
fc_run fc_run.cfg
```

# Unzip and polish
```bash
fc_unzip.py fc_unzip.cfg
```

# Extended phasing with HiC # not currently applicable.
```bash
fc_phase.py fc_phase.cfg
```


##### General information on configs: 
```bash
[General]
input_fofn=input.fofn # carriage return separated list of input fasta files with specified paths
input_type=raw # `raw` or `preads` raw will invoke `0-rawreads` pre-assembly phase, `preads` will skip
pa_DBdust_option=true # default is dusting is on and run after generating the raw read database. Can be modified with flag `pa_DBdust_option`
pa_fasta_filter_option=streamed-internal-median # default, uses median-length subread for each sequencing reaction well. Most users will not change...

[Data Partitioning]
# large genomes
pa_DBsplit_option=-x500 -s200 # the `-x` flag filters reads smaller than what's specified while the -s flag controls the size of DB blocks.
ovlp_DBsplit_option=-x500 -s200 # not sure what we should specify here....

# small genomes (<10Mb)
pa_DBsplit_option = -x500 -s50
ovlp_DBsplit_option = -x500 -s50

[Repeat Masking]
pa_HPCTANmask_option= # I don't know anything about this, but I guess I should consult this website:  https://dazzlerblog.wordpress.com/2016/04/01/detecting-and-soft-masking-repeats/
pa_REPmask_code=0,300;0,300;0,300

[Pre-assmbly]
genome_size=1000000000 # size of genome in base pairs
seed_coverage=30 # PacBio generally recommends 20-40x seed coverage.
length_cutoff=-1    # using -1 forces seed coverage auto-calculation, otherwise set `length_cutoff` to manually set limit
		# `length_cutoff` makes assembly faster, but anything smaller than that can't be used to phase in unzip. 
		# TRADEOFF: faster compute time in preassembly, but a smaller dataset for phasing.
pa_HPCdaligner_option=-v -B128 -M24 # confusing, for information see https://dazzlerblog.wordpress.com/2014/07/10/dalign-fast-and-sensitive-detection-of-all-pairwise-local-alignments/ and https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide/
pa_daligner_option=-e0.8 -l2000 -k18 -h480  -w8 -s100 # -e=average sequence identity (0.70, low quality data; 0.8, high quality data), higher values help prevent haplotype collapse; -l=minimum length of overlap (1000 [short library] - 5000 [longer library]); -k=kmer size (14 [low quality data] - 18 [high quality data])
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 3 --max-n-read 400 # set minimum alignment identity, minimum coverage, and max number of reads for consensus to make the preads
falcon_sense_greedy=False

[Pread overlapping]
ovlp_daligner_option=-e.96 -s1000 -h60 # -e=average sequence identity (0.93 inbread - 0.96 outbred)
ovlp_HPCdaligner_option=-v -M24 -l500  # -l=minimum length of overlap (800 [poor preassembly, short/low quality library] - 6000 [long, high quality library])

[Final Assembly]
overlap_filtering_setting=--max-diff 100 --max-cov 100 --min-cov 2 # `overlap_filter_setting` allows setting criteria for filtering pread overlaps; `--max-diff` filters overlaps that have coverage differences between the 5' and 3' ends; `--max-cov` filters highly represented overalps caused by contaminants or repeats; `--min-cov` allows specification of a minimum overlap coverage--setting this too low allows more overlaps to be detected at the expencse of additional chimeric/mis-assemblies
fc_ovlp_to_graph_option=
length_cutoff_pr=1000 # minimum length of pre-assembled preads used for final assmebly. Typically set to 15-30 fold coverage of corrected reads

[Miscellaneous configuration options]
target=assembly
skip_checks=False
LA4Falcon_preload=false

[job.defaults]
job_type=slurm
pwatcher_type=blocking # `fs_based` is the default and relies on the pipeline polling the file system periodically to determine whether a sentinel file has appeared that would signal the pipeline to continue; `blocking` process watcher which can help with systems that have issues with filesystem latency. In this case, the end of the job is determined by the finishing of the system call, rather than by file system polling
JOB_QUEUE = default
MB = 32768
NPROC = 6
njobs = 32
submit = qsub -S /bin/bash -sync y -V  \
  -q ${JOB_QUEUE}     \
  -N ${JOB_NAME}      \
  -o "${JOB_STDOUT}"  \
  -e "${JOB_STDERR}"  \
  -pe smp ${NPROC}    \
  -l h_vmem=${MB}M    \
  "${JOB_SCRIPT}"

[job.step.da]
NPROC=4
MB=49152
njobs=240
```


##########################################
############### To Submit ################
##########################################
```bash
srun --wait=0 -p myqueue -J ${JOB_NAME} -o ${JOB_STDOUT} -e ${JOB_STDERR} --mem-per-cpu=${MB}M --cpus-per-task=${NPROC} ${JOB_SCRIPT}
```

### FALCON-Unzip Configuration

```[General]
max_n_open_files = 1000 # can produce way too many `sam` files at a time, so they say maybe 300 files at a time...
[Unzip]
input_fofn=input.fofn # redundant with `fc_run.cfg` config file
input_bam_fofn=input_bam.fofn # for polishing genome, specified list of bam files.```
