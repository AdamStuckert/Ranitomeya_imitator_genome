### Analyses of duplicated orthologs

This document details our approach to examining why so many of the expected orthologs are present in duplicated copies.

`Minimap2` was used to align PacBio reads and `samtools depth` was used to calculate depth: 

```
#!/bin/bash
#SBATCH --partition=macmanes
#SBATCH -J bwaMini
#SBATCH --output PBMiniMap.log
#SBATCH --cpus-per-task=40
#SBATCH --mem=300000
#SBATCH --exclude=node117,node118

DIR="/mnt/lustre/macmaneslab/tml1019/seniorThesis/genomeFiles"
ASSEMBLY="/mnt/lustre/macmaneslab/tml1019/seniorThesis/genomeFiles/imitator.1.3.6.fa"
READS="/mnt/lustre/macmaneslab/ams1236/imitator_genome/reads/PacBio_reads.fa"

module purge
module load linuxbrew/colsa

cd /mnt/lustre/macmaneslab/tml1019/seniorThesis/BWAandDepth/PacBio_aligned_reads

# Run MiniMap becuase this is long read

minimap2 -x map-pb -I50g -N 10 -a -t 40 $ASSEMBLY $READS | samtools view -Sb - | samtools sort -T PacBio -O bam -@40 -l9 -m2G -o PacBioMiniMap.sorted.bam -

samtools index PacBioMiniMap.sorted.bam

# samtools depth
samtools depth -aa PacBioMiniMap.sorted.bam > PacBioMiniMapDepth.tsv
```

We calculated summary statistics (mean, median, mode, sd) on read depth:

```
#! /usr/bin/env python3

import statistics

# empty list
stats_list = []

try:
    with open("PacBioMiniMapDepth.tsv", "r") as input:
        for line in input:
            fields = line.split("\t")
            stats_list.append(int(fields[2]))
except IOError:
    print("problem reading file")

print("Mean:  ", statistics.mean(stats_list))
print("Median:  ", statistics.median(stats_list))
print("Mode:  ", statistics.mode(stats_list))
print("Stdv:  ", statistics.stdev(stats_list))
```

Pull out depth at specific duplicated regions:

```
#! /usr/bin/env python3

# *****************************************************
#   Script returns the depth at the bases inside the regions of intrest indicated in the BUSCO Output (duplicate/single copy)
# *****************************************************

#Initialize dictionary that will store duplicates
dupDictionary = dict()

#Read in file with list of genes from BUSCO output
with open("BUSCO4_duplicatedGenes.tsv", "r") as duplicatedGenes:
    for duplicatedLine in duplicatedGenes:
        duplicatedLine_stripped = duplicatedLine.strip()
        duplicate = duplicatedLine_stripped.split("\t")

        #Set the key to be the sacffold value
        key = duplicate[2]

        #If key doesnt exist in dictionary, make it
        dupDictionary.setdefault(key, [])

        #Get the BUSCO ID, the start and the end of the gene from the BUSCO output
        valuesToAppend = [duplicate[0], duplicate[3], duplicate[4]]

        #Add those values to the dictionary
        dupDictionary[key].append(valuesToAppend)

#Read in genome depth file

with open("PacBioMiniMapDepth.tsv", "r") as genomeDepth:
    for genomeLine in genomeDepth:
        genomeLine_stripped = genomeLine.strip()
        genome = genomeLine_stripped.split("\t")

        scaffold_name = (genome[0])
        coordinates = int(genome[1])

        #If the current line in the genome depth file contains a scaffold from the BUSCO file, report the depth at each base in that region of the genome, output to a TSV
        if scaffold_name in dupDictionary:

            for value in dupDictionary[scaffold_name]:
                endofDup = value[2]
                endofDup = int(endofDup)

                startofDup = value[1]
                startofDup = int(startofDup)

                # Get +- 5kb on either side of the gene of intrest
                lower = int(startofDup) - 50000
                upper = int(endofDup) + 50000

                if coordinates >= (lower) and coordinates <= (upper):
                    print(value[0], "\t", scaffold_name, "\t", coordinates, "\t", genome[2])
```

Calculate number of duplicated orthologs per scaffold:

```
#Load Packages

library(plyr)
library(dplyr)

#Read in File
duplicated_genes <- read.table("BUSCO4_duplicated.tsv", sep = "\t", header = T, stringsAsFactors = F)

head(duplicated_genes)
colnames(duplicated_genes) <- c("Busco_id", "Status", "Scaffold", "Start", "End", "Score", "Length", "OrthoDB", "Gene")

# *****************************************************
#   Count the number of duplicates on each Scaffold
# *****************************************************

scaffold_count <- count(duplicated_genes, Scaffold)

write.csv(scaffold_count, "scaffold_count.csv")
```

Calculate number of duplicated orthologs per scaffold normalized by scaffold length:

```
#! /usr/bin/env python3

# *****************************************************
#   Number of duplicates on each Scaffold Normalized by Scaffold Length
# *****************************************************

comboDictionary = dict()

#Read in file with the scaffold name and its length; input file is the index of the genome from samtools
with open("imitator.1.3.6.fa.fai", "r") as duplicatedGenes:
    for duplicatedLine in duplicatedGenes:
        duplicatedLine_stripped = duplicatedLine.strip()
        duplicate = duplicatedLine_stripped.split(",")

        #Key = scaffold name
        key = duplicate[0]

        #If key doesnt exist it makes one
        comboDictionary.setdefault(key, [])

        #Adds the scaffold length to the dictionary
        valuesToAppend = [duplicate[1]]
        comboDictionary[key].append(valuesToAppend)

#Reads in a file with the number of duplicates per scaffold
with open("scaffold_count.csv", "r") as scaffoldcount:
    for line in scaffoldcount:
        line_stripped = line.strip()
        scaffold = line_stripped.split(",")

        key = scaffold[0]

        comboDictionary.setdefault(key, [])

        #Adds the number of duplicates by scaffold
        valuesToAppend = [scaffold[1]]
        comboDictionary[key].append(valuesToAppend)

percentagePerScafDict = dict()

for key in comboDictionary:

    #takes the number of duplicates and normalizes it by scaffold length
    entry1and2 = []
    for entry in comboDictionary[key]:
        for line in entry:
            entry1and2.append(line)

    try:
        print(key, "\t", int(entry1and2[1])/int(entry1and2[0]))
    except:
        continue
```

Do both copies of orthologs ever appear on the same scaffold?

```
#! /usr/bin/env python3

# *****************************************************
#   If the duplicate and its complement are on the same scaffold it will print the BUSCO ID and the scaffold that they appear on
# *****************************************************

dupDictionary = dict()

#Read in BUSCO Output
with open("BUSCO4_duplicated.tsv", "r") as duplicatedGenes:
    for duplicatedLine in duplicatedGenes:
        duplicatedLine_stripped = duplicatedLine.strip()
        duplicate = duplicatedLine_stripped.split("\t")

        #key = BUSCO ID
        key = duplicate[0]

        #If key doesnt already exist, make it
        dupDictionary.setdefault(key, [])

        #Gets the scaffold of the duplicate
        valuesToAppend = [duplicate[2]]

        #Adds it to dictionary
        dupDictionary[key].append(valuesToAppend)

    for key in dupDictionary:

        #Initialize a list to capture the two values
        entry1and2 = []
        for entry in dupDictionary[key]:
            for line in entry:
                entry1and2.append(line)


        #If the duplicate and its complement are on the same scaffold it will print the BUSCO ID and the scaffold it is on
        #If not, it will break the loop and move on

        item = entry1and2[0]
        check = True
        for i in entry1and2:
            if item != i:
                check = False
                break
        if (check == True):
            print(key, dupDictionary[key])
```

We then calculated sliding window statistics for the data using an [R script that we ran on the cluster](ADD LINK HERE). The results of this script were then used to compare orthologs present in single copies and duplicated copies in a [subsequent R script which we ran locally](ADD LINK HERE).
