Aligning and mapping RNA seq reads from what were going to be a few different projects.

_Ranitomeya imitator_ first:

```bash
sbatch AlignmentReadCount.job  /mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa  /mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6_noRepeatModeler_test/Ranitomeya_imitator.imitator.1.3.6.fa.gff3  /mnt/lustre/macmaneslab/ams1236/devseries/readfiles .fastq.gz
```

Since I have already indexed the genome, I'll submit a slightly modified script without the indexing step for the other species.

_R fantastica run:
```bash
sbatch AlignmentReadCount.job  /mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa  /mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6_noRepeatModeler_test/Ranitomeya_imitator.imitator.1.3.6.fa.gff3  /mnt/lustre/macmaneslab/ams1236/MultispeciesDevSeries/readfiles/fantastica_reads .fq.gz
```

_R variabilis_ run:
```bash
sbatch AlignmentReadCount.job  /mnt/lustre/macmaneslab/ams1236/imitator_genome/imitator.1.3.6.fa  /mnt/lustre/macmaneslab/ams1236/imitator_genome/maker_1.3.6_noRepeatModeler_test/Ranitomeya_imitator.imitator.1.3.6.fa.gff3  /mnt/lustre/macmaneslab/ams1236/MultispeciesDevSeries/readfiles/variabilis_reads .fq.gz
```
