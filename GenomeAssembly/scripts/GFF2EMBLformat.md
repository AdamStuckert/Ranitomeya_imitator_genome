## Convert genome annotations into the EMBL format for submission

This is the code I used for getting my genome annotations from Maker into the correct format for deposition in the European Nucleotide Archive.

First step, install [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3). I did this via conda using a `yml` file I wrote to eliminate issues with software versions that were automatically installed and don't work with EMBLmyGFF3.

```bash
conda env create -f myEMBLmyGFF3.yml
```



I then copied everything into a new directory, juuuuuust in case.

```bash
mkdir gff2embl
cp imitator.1.3.6.fa gff2embl
cp Ranitomeya_imitator.imitator.1.3.6.functional.gff3 gff2embl
```

I also cleaned up the scaffold names a bit.

```bash
cd gff2embl
sed "s/_pilon//g" imitator.1.3.6.fa > Ranitomeya_imitator_genomeassembly_1.0.fa
sed "s/_pilon//g" Ranitomeya_imitator.imitator.1.3.6.functional.gff3 > Ranitomeya_imitator_genomeassembly_1.0.gff3
```

Create the flat file. For some reason I had issues with this just submitting it as regular code. Using their example maker shell script seemed to be a functional work around.

```bash
module purge
conda activate embl
./maker_imitator.sh
```

Contents of `maker_imitator.sh`:

```
#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#PATH to the FASTA file used to produce the annotation
GENOME=`dirname "$0"`"/Ranitomeya_imitator_genomeassembly_1.0.fa"

#PATH to the ANNOTATION in gff3 FORMAT
ANNOTATION=`dirname "$0"`"/Ranitomeya_imitator_genomeassembly_1.0.gff3"

#PROJECT name registered on EMBL
PROJECT="PRJEB28312 "

#Locus tag registered on EMBL
LOCUS_TAG="RIMI"

# species name
SPECIES="Ranitomeya imitator"

# Taxonomy
#TAXONOMY="VERT"

#The working groups/consortia that produced the record. No default value
#REFERENCE_GROUP="XXX"

#Translation table
TABLE="1"

#Molecule type of the sample.
MOLECULE="genomic DNA"

myCommand="EMBLmyGFF3 -i $LOCUS_TAG -p $PROJECT -m \"$MOLECULE\" -r $TABLE -t linear -s \"$SPECIES\" -o Ranitomeya_imitator_genomeassembly_1.0.embl $ANNOTATION $GENOME $@"
echo -e "Running the following command:\n$myCommand"

#execute the command
eval $myCommand

```
        
        
And finally, the manifest file to submit!

```
STUDY   PRJEB28312
SAMPLE   ERS5098337
ASSEMBLYNAME   Ranitomeya_imitator_genomeassembly_1.0
ASSEMBLY_TYPE	isolate
COVERAGE   35X
PROGRAM   wtdbg2
PLATFORM   PacBio Sequel RSII, Oxford Nanopore Technologies, Illumina HiSeq 2500
MINGAPLENGTH   50
MOLECULETYPE   genomic DNA
FASTA   Ranitomeya_imitator_genomeassembly_1.0.fa.gz
```
