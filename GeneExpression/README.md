### Gene expression analyses

This folder documents work done to analyze our gene expression after assembling the genome. We parameterized aligning and counting RNA seq reads using `STAR` and `htseq-count.` For details on this, please see the [GeneExpressionWorkbook.md](https://github.com/AdamStuckert/Ranitomeya_imitator_genome/blob/master/GeneExpression/GeneExpressionWorkbook.md) file. We used a number of scripts we wrote for this purpose in the parameterization process, notably scripts to [index the genome, align reads, and count expression](https://github.com/AdamStuckert/Ranitomeya_imitator_genome/blob/master/GeneExpression/AlignmentReadCount.job) as well as one that just [alignsreads, and countsexpression](https://github.com/AdamStuckert/Ranitomeya_imitator_genome/blob/master/GeneExpression/ReadCount.job).

Following this we conducted differential expression analyses using the R package DESeq2. We have a [detailed R markdown file of our differential expression analyses](https://github.com/AdamStuckert/Ranitomeya_imitator_genome/blob/master/GeneExpression/MimeticGeneExpressionGeneLevel.Rmd). In addition to doing analyses, this script produces a number of the figures in the publication. Relevant data for this script (e.g., gene counts for each sample, genome annotation documents used) are all in the `./data/` directory. 

Additionally, we used the WGCNA package in R to [examine how genes are coexpressed and their correlation to treatment groups](https://github.com/AdamStuckert/Ranitomeya_imitator_genome/blob/master/GeneExpression/MimeticWGCNA.Rmd).

Finally, this folder has [the script we used to produce various figures describing the genome that are in the publication](https://github.com/AdamStuckert/Ranitomeya_imitator_genome/blob/master/GeneExpression/GenomeLengthFigure.R). Does this really belong here? No probably not. But it is a pandemic and I'm in charge of this repository.