## Convert genome annotations into the EMBL format for submission

This is the code I used for getting my genome annotations from Maker into the correct format for deposition in the European Nucleotide Archive.

First step, install [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3). I did this via conda.

```bash
conda create --embl
conda activate embl
conda install -c bioconda emblmygff3
```

I then copied everything into a new directory, juuuuuust in case.

```bash
mkdir gff2embl
cp imitator.1.3.6.fa gff2embl
cp $GFF gff2embl
```

I also cleaned up the scaffold names a bit.

```bash
cd gff2embl
sed "s/_pilon//g" imitator.1.3.6.fa > Ranitomeya_imitator_genomeassembly_1.0.fa
sed "s/_pilon//g" Ranitomeya_imitator.imitator.1.3.6.functional.gff3 > Ranitomeya_imitator_genomeassembly_1.0.gff3
```

Create the flat file

```bash
conda activate embl
EMBLmyGFF3 Ranitomeya_imitator_genomeassembly_1.0.fa Ranitomeya_imitator_genomeassembly_1.0.gff3 \
        --data_class WGS \
        --topology linear \
        --molecule_type "genomic DNA" \
        --transl_table 1  \
        --translate \
        --species 'Ranitomeya imitator' \
        --locus_tag RIMITATOR \
        --project_id PRJEB28312 \
        -o Ranitomeya_imitator_genomeassembly_1.0.embl
```
        
        
