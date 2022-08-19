Aligning and mapping RNA seq reads from what were going to be a few different projects.

_Ranitomeya imitator_ first:

```bash
sbatch AlignmentReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/ratshit/R_imi_1.0.fa \
/mnt/lustre/macmaneslab/ams1236/ratshit/annotation2/Maker_round4/Ranitomeya_imitator.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/imitator_data .fastq.gz
```

Since I have already indexed the genome, I'll submit a slightly modified script without the indexing step for the other species.

_R fantastica run:
```bash
sbatch --output RNAseqReadCountVariabilis.log ReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/ratshit/R_imi_1.0.fa \
/mnt/lustre/macmaneslab/ams1236/ratshit/annotation2/Maker_round4/Ranitomeya_imitator.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/variabilis_reads .fq.gz
```

_R variabilis_ run:
```bash
# note submitted from: /mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/test
# I want to override the output from the sbatch header in the script as well
sbatch --output RNAseqReadCountVariabilis.log ReadCount.job  \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa \
/mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6.masked_28April/Ranitomeya_imitator.imitator.1.3.6.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MultispeciesDevSeries/readfiles/variabilis_reads .fq.gz
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

If we decrease read length threshold in STAR:

```bash
sbatch --output RNAseqReadCountImitator_trimmed.log ReadCount_50percent_readlength.job \
/mnt/lustre/macmaneslab/ams1236/ratshit/R_imi_1.0.fa \
/mnt/lustre/macmaneslab/ams1236/ratshit/annotation2/Maker_round4/Ranitomeya_imitator.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/imitator_data .fastq.gz

sbatch --output RNAseqReadCountFantastica_trimmed.log ReadCount_50percent_readlength.job \
/mnt/lustre/macmaneslab/ams1236/ratshit/R_imi_1.0.fa \
/mnt/lustre/macmaneslab/ams1236/ratshit/annotation2/Maker_round4/Ranitomeya_imitator.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/fantastica_data .fq.gz


sbatch --output RNAseqReadCountVariabilis_trimmed.log ReadCount_50percent_readlength.job \
/mnt/lustre/macmaneslab/ams1236/ratshit/R_imi_1.0.fa \
/mnt/lustre/macmaneslab/ams1236/ratshit/annotation2/Maker_round4/Ranitomeya_imitator.functional.gff3 \
/mnt/lustre/macmaneslab/ams1236/MimicryGeneExpression/variabilis_reads .fq.gz
```


Following this, read counts were downloaded to a local machine and differential expression analyses and weighted gene coexpression network analyses (WCGNA) were conducted. Those are R scripts here.
