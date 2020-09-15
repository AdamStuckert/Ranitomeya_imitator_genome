#!/bin/bash

### pull out star mapping data:
printf "Sample\tNumber_reads\tUniquely_mapping_reads\tMapped_to_multiple_loci\tMapped_to_too_many_loci\tUnmapped_reads_too_short\n" > mappingdata.tab
files=$(ls *Log.final.out)

for file in $files
do
sample=$(echo $file | sed "s/_Log.final.out//g" | sed "s/.Log.final.out//g")
reads=$(grep "Number of input reads" $file | cut -f 2)
unique=$(grep "Uniquely mapped reads %" $file | cut -f 2) 
multi=$(grep "% of reads mapped to multiple loci" $file | cut -f 2)
toomany=$(grep "% of reads mapped to too many loci" $file | cut -f 2)
unmapped=$(grep "reads unmapped: too short" $file | cut -f 2)
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$reads" "$unique" "$multi" "$toomany" "$unmapped" >> mappingdata.tab
done


### Calculate read mapping statistics with a custom R script I wrote
STATS=$(which ReadMappingStats.R)
Rscript $STATS
