#!/usr/bin/sh
#SBATCH --mem=700GB
#SBATCH --job-name="imi_assembly"
#SBATCH --output="imi_assembly.log"
#SBATCH --partition=macmanes
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118


DIR=$(pwd)
READ1="${HOME}/ratshit/data/m64019_190803_173458.subreads.fa"
READ2="${HOME}/ratshit/data/m64019_190913_173453.subreads.fasta"
READ3="${HOME}/ratshit/data/m64019_190915_200441.subreads.fasta"
TENX1="${HOME}/ratshit/data/trimmed10x.R1.fastq.gz"
TENX2="${HOME}/ratshit/data/trimmed10x.R2.fastq.gz"
REF="imitator.wtdbg.ctg.fa"
REF_PATH="${DIR}/wtdbg/${REF}"


module purge
# load conda source script!
. ${HOME}/anaconda3/etc/profile.d/conda.sh


mkdir wtdbg
cd wtdbg

wtdbg2 \
-o imitator.wtdbg \
-x sq \
-g 6.7g \
-L 2500 \
-t 40 \
-i ${READ1} \
-i ${READ2} \
-i ${READ3}

wtpoa-cns -t 40 -i imitator.wtdbg.ctg.lay.gz -fo $REF

# rename fasta headers!
awk '{print $1}' $REF > new.fasta
mv new.fasta $REF

#### Calculate genome metrics
printf "Submitting genome metrics job\n\n\n"
sbatch $HOME/scripts/genomemetrics.job $REF $(echo $REF | sed "s/.fa//")

#### polish with Pilon
cat << EOF > bwa.job
#!/bin/bash
#SBATCH --partition=macmanes
#SBATCH -J bwa
#SBATCH --output bwa.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude node117,node118


#REF=$REF

species="Ranitomeya_imitator"
prefix="Rimi"


module purge
module load linuxbrew/colsa

echo Running BWA on $REF


mkdir pilon
cd pilon

bwa index ${REF_PATH}

rm -f "$prefix".sorted.bam
rm -f "$prefix".sorted.bam.bai

bwa mem -t 40 ${REF_PATH} ${Tenx1} ${Tenx2} \
| samtools view -@20 -Sb - \
| samtools sort -T "$prefix" -O bam -@20 -l9 -m2G -o "$prefix".sorted.bam -
samtools index "$prefix".sorted.bam


### Split assembly into chunks to speed up mapping

# make folder
mkdir chunks

# remove any chunks and list from previous iteration
rm chr.list
rm chunks/*

# genome headers
grep ">" ${REF_PATH} | sed 's_>__' | shuf | tee -a chr.list


# slit into 80 chunks
split -d -n l/80 chr.list chunks/genomechunk.

cd chunks/
rename genomechunk.0 genomechunk. genomechunk.0*
printf "Submitting pilon array job\n\n\n"
sbatch ../pilon.job
EOF

printf "Submitting bwa alignment and chunk splitting job\n\n\n"
sbatch bwa.job


################ pilon array job
printf "Writing array job for pilon polishing\n\n\n"
cat << "EOF" > pilon.job
 #!/bin/bash
 #SBATCH --partition=macmanes,shared
 #SBATCH -J pilon
 #SBATCH --output logs4pilon/pilon.%A_%a.log
 #SBATCH --array=0-79%5
 #SBATCH --exclude=node117,node118
 #SBATCH --dependency=afterany:281601

EOF

cat << EOF >> pilon.job
 REF=$REF
 REF_OUT=$(echo $REF | sed "s/.fa//")
 REF_OUT+=".pilonpolished.fa"
 species="Ranitomeya_imitator"
 prefix="Rimi"

EOF

cat << "EOF" >> pilon.job
 echo "SLURM_JOBID: " $SLURM_JOBID
 echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
 echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

 DIR=$(pwd)
 java -jar -Xmx105G /mnt/lustre/macmaneslab/macmanes/pilon-1.23.jar --genome ${REF_PATH} \
 --frags ${DIR}/pilon/"$prefix".sorted.bam \
 --output ${DIR}/chunks/pilonchunk.$SLURM_ARRAY_TASK_ID \
 --fix bases,gaps \
 --diploid \
 --threads 24 \
 --flank 5 \
 --verbose \
 --mingap 1 \
 --nostrays \
 --targets ${DIR}/wtdbg/chunks/genomechunk.$SLURM_ARRAY_TASK_ID

 cat ${DIR}/chunks/pilon*fasta > ${DIR}/${REF_OUT}

 ## genome metrics
 sbatch $HOME/scripts/genomemetrics.job ${DIR}/${REF_OUT} $(echo ${DIR}/${REF_OUT} | sed "s/.fa//")

EOF





#### polish with racon
cat << EOF > racon.job
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J racon
#SBATCH --output racon.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=700Gb


module load linuxbrew/colsa

# preparation
mkdir racon
cd racon

ln -s $REF_PATH

### First align reads with minimap2
echo aligning with minimap2
minimap2 -I10G -t 40 -xmap-pb $REF $READ1 $READ2 $READ3 | gzip -c - > Imi.PB.paf.gz

## merge reads into a single file for racon...
cat $READ1 $READ2 $READ3 > AllReads.fa

### Run racon
echo Polishing with racon
racon -t 40 AllReads.fa  Imi.PB.paf.gz $REF > ${REF}.raconpolished.fa

## genome metrics
sbatch $HOME/scripts/genomemetrics.job ${REF}.raconpolished.fa ${REF}.raconpolished

EOF

printf "Submitting pilon polishing job\n\n\n"
sbatch racon.job
