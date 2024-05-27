Aligning and mapping RNA seq reads from what were going to be a few different projects.

_Ranitomeya imitator_ first:

```bash
sh /project/stuckert/users/Stuckert/scripts/STAR.sh -a /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/R_imi_1.0.fa -g /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/Maker_round3cat/R_imi.UNK.merged.0.7.functional.gff3 -r /project/stuckert/users/Stuckert/MultiSpeciesDevSeries/imitator_reads -s .fastq.gz
```

Since I have already indexed the genome, I'll submit a slightly modified script without the indexing step for the other species.

_R fantastica run:
```bash
sh /project/stuckert/users/Stuckert/scripts/STAR.sh -a /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/R_imi_1.0.fa -g /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/Maker_round3cat/R_imi.UNK.merged.0.7.functional.gff3 -r /project/stuckert/users/Stuckert/MultiSpeciesDevSeries/fantastica_reads -s .fq.gz
```

_R variabilis_ run:
```bash
sh /project/stuckert/users/Stuckert/scripts/STAR.sh -a /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/R_imi_1.0.fa -g /project/stuckert/users/Stuckert/mimicry_genome/R_imi.UNK.annotation/Maker_round3cat/R_imi.UNK.merged.0.7.functional.gff3 -r /project/stuckert/users/Stuckert/MultiSpeciesDevSeries/variabilis_reads -s .fq.gz
```

Rename all the populations of imitator reads.

```bash
for i in $(ls *counts | grep -v "R_"); do mv -- "$i" "R_imitator_$i"; done

#### 
# Huallaga
for i in $(ls *H*.counts); do mv -- "$i" "${i/H/striped_}"; done

## Sauce
for i in $(ls *S*.counts); do mv -- "$i" "${i/S/banded_}"; done

## Tarapoto
for i in $(ls *T*.counts); do mv -- "$i" "${i/T/spotted_}"; done

## Varadero
for i in $(ls *V*.counts); do mv -- "$i" "${i/V/redheaded_}"; done

# change - to _
for i in $(ls R_imitator*); do mv -- "$i" "${i/-/_}"; done 
```

Get mapping data:

```bash
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
```
