#!/usr/bin/make -f
# Pipeline for the ARCS program
# Written by Jeffrey Tse
#Default Parameters

# Input Names
draft=draft
reads=reads

# tigmint Parameters
minsize=2000
as=0.65
nm=5
dist=50000
mapq=0
trim=0
span=20
window=1000

# ARCS ARKS Common Parameters
barcodeMultiplicityFile=barcodeMultiplicityArcs
c=5
m=50-10000
z=500
s=98
r=0.05
e=30000
D=false
dist_upper=false
d=0
gap=100
B=20

# ARCS Specific Paramters
s=98

# ARKS Specific Parameters
j=0.55
k=30
t=40

# LINKS Parameters
l=5
a=0.3
bin=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))

SHELL=bash -e -o pipefail
ifeq ($(shell zsh -e -o pipefail -c 'true' 2>/dev/null; echo $$?), 0)
#Set pipefail to ensure that all commands of a pipe succeed.
SHELL=zsh -e -o pipefail
# Report run time and memory usage with zsh.
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
endif

# Record run time and memory usage in a file using GNU time.
ifdef time
ifneq ($(shell command -v gtime),)
gtime=command gtime -v -o $@.time
else
gtime=command time -v -o $@.time
endif
endif

.PHONY: all version help clean tigmint arcs arcs-tigmint arcs-with-tigmint arks arks-tigmint arks-with-tigmint
.DELETE_ON_ERROR:
.SECONDARY:


all: help
# Help
help:
        @echo "Usage: ./arcs-make [COMMAND] [OPTION=VALUE]..."
        @echo "    Commands:"
        @echo ""
        @echo " arcs            run arcs in default mode only, skipping tigmint"
        @echo " arcs-tigmint    run tigmint, and run arcs in default mode with the output of tigmint"
        @echo " arks            run arcs in kmer mode only, skipping tigmint"
        @echo " arks-tigmint    run tigmint, and run arcs in kmer mode with the output of tigmint"
        @echo " help            display this help page"
        @echo " version         display the software version"
        @echo " clean           remove intermediate files"
        @echo ""
        @echo "    General Options:"
        @echo ""
        @echo " draft           draft name [draft]. File must have .fasta or .fa extension"
        @echo " reads           read name [reads]. File must have .fastq.gz or .fq.gz extension"
        @echo " time            logs time and memory usage to file for main steps (Set to 1 to enable logging)"
        @echo ""
        @echo "    bwa Options:"
        @echo ""
        @echo " t               number of threads used [8]"
        @echo ""
        @echo "    Tigmint Options:"
        @echo ""
        @echo " minsize         minimum molecule size [2000]"
        @echo " as              minimum AS/read length ratio [0.65]"
        @echo " nm              maximum number of mismatches [5]"
        @echo " dist            max dist between reads to be considered the same molecule [50000]"
        @echo " mapq            mapping quality threshold [0]"
        @echo " trim            bp of contigs to trim after cutting at error [0]"
        @echo " span            min number of spanning molecules to be considered assembled [20]"
        @echo " window          window size for checking spanning molecules [1000]"
        @echo ""
        @echo "    Common Options:"
        @echo ""
        @echo " c               minimum aligned read pairs per barcode mapping [5]"
        @echo " m               barcode multiplicity range [50-10000]"
        @echo " z               minimum contig length [500]"
        @echo " r               p-value for head/tail assigment and link orientation [0.05]"
        @echo " e               contig head/tail length for masking aligments [30000]"
        @echo " D               enable distance estimation [false]"
        @echo " dist_upper      use upper bound distance over median distance [false]"
        @echo " B               estimate distance using N closest Jaccard scores [20]"
        @echo " d               max node degree in scaffold graph [0]"
        @echo " gap             fixed gap size for dist.gv file [100]"
        @echo " barcode-counts  name of output barcode multiplicity TSV file [barcodeMultiplicityArcs]"
        @echo ""
        @echo " ARCS Specific Options:"
        @echo " s               minimum sequence identity [98]"
        @echo ""
        @echo " ARKS Specific Options:"
        @echo " j               minimum fraction of read kmers matching a contigId [0.55]"
        @echo " k               size of a k-mer [30]"
        @echo " t               number of threads [8]"
        @echo ""
        @echo "    LINKS Options:"
        @echo ""
        @echo " l               minimum number of links to compute scaffold [5]"
        @echo " a               maximum link ratio between two best contig pairs [0.3]"
        @echo ""
        @echo "Example: To run tigmint and arcs with myDraft.fa, myReads.fq.gz, and a custom multiplicty range, run:"
        @echo " ./arcs-make arcs-tigmint draft=myDraft reads=myReads m=[User defined multiplicty range]"
        @echo "To ensure that the pipeline runs correctly, make sure that the following tools are in your PATH: bwa, tigmint, samtools, arcs (>= v1.1.0), LINKS (>= v1.8.6)"

clean:
        rm -f *.amb *.ann *.bwt *.pac *.sa *.dist.gv *.fai *.bed *.molecule.tsv *.sortbx.bam
        @echo "Clean Done"

version:
        @echo "arcs-make v1.1.0"

#Preprocessing

# Create a .fa file that is soft linked to .fasta
%.fa: %.fasta
        ln -s $^ $@

# Create a .fq.gz file that is soft linked to .fastq.gz
%.fq.gz: %.fastq.gz
        ln -s $^ $@


#Run Tigmint
arcs-tigmint: tigmint arcs-with-tigmint
arks-tigmint: tigmint arks-with-tigmint

# Main
tigmint: $(draft).tigmint.fa
# Run tigmint
$(draft).tigmint.fa: $(draft).fa $(reads).fq.gz
        $(gtime) tigmint tigmint draft=$(draft) reads=$(reads) minsize=$(minsize) as=$(as) nm=$(nm) dist=$(dist) mapq=$(mapq) trim=$(trim) span=$(span) window=$(window) t=$t

#Run ARCS
arcs: $(draft)_c$c_m$m_s$s_r$r_e$e_z$z_l$l_a$a.scaffolds.fa
arcs-with-tigmint: $(draft).tigmint_c$c_m$m_s$s_r$r_e$e_z$z_l$l_a$a.scaffolds.fa
arks: $(draft)_c$c_m$m_k$k_r$r_e$e_z$z_l$l_a$a.scaffolds.fa
arks-with-tigmint: $(draft).tigmint_c$c_m$m_k$k_r$r_e$e_z$z_l$l_a$a.scaffolds.fa

# Convert Scaffold Names into Numerical Numbers
%.renamed.fa: %.fa
        perl -ne 'chomp; if(/>/){$$ct+=1; print ">$$ct\n";}else{print "$$_\n";} ' < $^ > $@

# Make bwa index from Draft Assembly
%.renamed.fa.bwt: %.renamed.fa
        $(gtime) bwa index $^

# Use bwa mem to Align Reads to Draft Assembly and Sort it
%.sorted.bam: %.renamed.fa $(reads).fq.gz %.renamed.fa.bwt
        $(gtime) sh -c 'bwa mem -t$t -C -p $< $(reads).fq.gz | samtools view -@$t -Sb - | samtools sort -@$t -l9 -m2G -n - -o $@'

# Create an fof File Containing the bam File
%_bamfiles.fof: %.sorted.bam
        echo $^ > $@

# Run ARCS Program
%_c$c_m$m_s$s_r$r_e$e_z$z_original.gv: %.renamed.fa %_bamfiles.fof
ifneq ($D, true)
        $(gtime) arcs -f $< -a $(word 2,$^) -v -c $c -m $m -s $s -r $r -e $e -z $z -d $d --gap $(gap) -b $(patsubst %_original.gv,%,$@) --barcode-counts $(barcodeMultiplicityFile).tsv
else ifneq ($(dist_upper), true)
        $(gtime) arcs -D -B $B -v -f $< -a $(word 2,$^) -c $c -m $m -s $s -r $r -e $e -z $z -d $d --gap $(gap) -b $(patsubst %_original.gv,%,$@) --barcode-counts $(barcodeMultiplicityFile).tsv
else
        $(gtime) arcs -D -B $B --dist_upper -v -f $< -a $(word 2,$^) -c $c -m $m -s $s -r $r -e $e -z $z -d $d --gap $(gap) -b $(patsubst %_original.gv,%,$@) --barcode-counts $(barcodeMultiplicityFile).tsv
endif

# kmer
%_c$c_m$m_k$k_r$r_e$e_z$z_original.gv: %.renamed.fa $(reads).fq.gz
ifneq ($D, true)
        $(gtime) arcs --arks -v -f $< -c $c -m $m -r $r -e $e -z $z -j $j -k $k -t $t -d $d --gap $(gap) -b $(patsubst %_original.gv,%,$@) $(word 2,$^) --barcode-counts $(barcodeMultiplicityFile).tsv

else ifneq ($(dist_upper), true)
        $(gtime) arcs --arks -v -D -B $B -f $< -c $c -m $m -r $r -e $e -z $z -j $j -k $k -t $t -d $d --gap $(gap) -b $(patsubst %_original.gv,%,$@) $(word 2,$^) --barcode-counts $(barcodeMultiplicityFile).tsv
else
        $(gtime) arcs --arks -v -D -B $B --dist_upper -f $< -c $c -m $m -r $r -e $e -z $z -j $j -k $k -t $t -d $d --gap $(gap) -b $(patsubst %_original.gv,%,$@) $(word 2,$^) --barcode-counts $(barcodeMultiplicityFile).tsv
endif

# Generate TSV from ARCS
# default
%_c$c_m$m_s$s_r$r_e$e_z$z.tigpair_checkpoint.tsv: %_c$c_m$m_s$s_r$r_e$e_z$z_original.gv %.renamed.fa
        python /mnt/lustre/macmaneslab/ams1236/scripts/makeTSVfile.py $< $@ $(word 2,$^)
# kmer
%_c$c_m$m_k$k_r$r_e$e_z$z.tigpair_checkpoint.tsv: %_c$c_m$m_k$k_r$r_e$e_z$z_original.gv %.renamed.fa
        python /mnt/lustre/macmaneslab/ams1236/scripts/makeTSVfile.py $< $@ $(word 2,$^)

# Adds a and l paramters to the filename
%_z$z_l$l_a$a.tigpair_checkpoint.tsv: %_z$z.tigpair_checkpoint.tsv
        ln -s $^ $@

# Make an Empty fof File
empty.fof:
        touch $@

# Run LINKS
# default
%_c$c_m$m_s$s_r$r_e$e_z$z_l$l_a$a.scaffolds.fa: %.renamed.fa empty.fof %_c$c_m$m_s$s_r$r_e$e_z$z_l$l_a$a.tigpair_checkpoint.tsv
        $(gtime) LINKS -f $< -s empty.fof -b $(patsubst %.scaffolds.fa,%,$@) -l $l -a $a -z $z
# kmer
%_c$c_m$m_k$k_r$r_e$e_z$z_l$l_a$a.scaffolds.fa: %.renamed.fa empty.fof %_c$c_m$m_k$k_r$r_e$e_z$z_l$l_a$a.tigpair_checkpoint.tsv
        $(gtime) LINKS -f $< -s empty.fof -b $(patsubst %.scaffolds.fa,%,$@) -l $l -a $a -z $z
