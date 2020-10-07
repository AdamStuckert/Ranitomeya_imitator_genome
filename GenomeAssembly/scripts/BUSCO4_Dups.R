#/usr/bin/env/R

# GenomeDuplicates.R
# Author: Adam Stuckert


## Note: this script uses read depth data from Illumina, PacBio, and Nanopore sequencing. For the submitted manuscript we are only including PacBio data because this was the primary sequencing technology used.

## load libraries
library(ggplot2)
library(tidyverse)
library(data.table)
library(RcppRoll)

## Data import
# import PacBio data
PBdat <- fread("data/duplicates/PacBiodepthAtDuplicatedRegions_BUSCO4_duplicates.tsv", header = FALSE, sep = "\t")
colnames(PBdat) <- c("BUSCO_ID", "scaffold","base","PacBio_depth")

### import ONT data
ONTdat <- fread("data/duplicates/NanoporedepthAtDuplicatedRegions_BUSCO4_duplicates.tsv", header = FALSE, sep = "\t")
colnames(ONTdat) <- c("BUSCO_ID", "scaffold","base","ONT_depth")

# import Illumina data
Illdat <- fread("data/duplicates/IlluminadepthAtDuplicatedRegions_BUSCO4_duplicates.tsv", header = FALSE, sep = "\t")
colnames(Illdat) <- c("BUSCO_ID", "scaffold","base","Illumina_depth")

# merge them all
dat <- cbind(PBdat, ONTdat$ONT_depth, Illdat$Illumina_depth)
colnames(dat) <- c("BUSCO_ID", "scaffold","base","PacBio_depth", "ONT_depth", "Illumina_depth")

# clean up extra data from workspace to save speed/make my life tolerable
rm(PBdat, ONTdat, Illdat)

# Import BUSCO duplicate information
busco <- fread("data/duplicates/BUSCO4_duplicated.tsv", header = FALSE, sep = "\t")
colnames(busco) <- c("BUSCO_ID", "status", "scaffold", "dup_start", "dup_end", "score", "length", "query", "gene_name")

# remove pilon
busco$scaffold <- gsub(pattern="scaffold", replacement = "", busco$scaffold)
busco$scaffold <- gsub(pattern="_pilon_pilon", replacement = "", busco$scaffold)

# merge column names
busco <- busco %>% tidyr::unite("dup_scaf", BUSCO_ID, scaffold, sep = "_", remove = FALSE)

# remove unneeded columns
buscodat <- busco[,c("dup_scaf", "dup_start", "dup_end", "gene_name")]

## Data manipulation
# remove "pilon": probably substantially faster to do in bash, but ok here.
dat$scaffold <- gsub(pattern="scaffold", replacement = "", dat$scaffold)
dat$scaffold <- gsub(pattern="_pilon_pilon", replacement = "", dat$scaffold)

# add a column for merged duplicate by scaffold
dat <- dat %>% tidyr::unite("dup_scaf", BUSCO_ID, scaffold, sep = "_", remove = FALSE)

# add start + stop of duplicate data
dat <- left_join(dat, buscodat, by = "dup_scaf")

## Caculate sliding window statistics
# verify multiple dup by scaf occurrences:
numunique = length(unique(dat$dup_scaf))
cat(paste0("There are ", numunique, " unique combinations of duplicates by scaffolds.\n"))

# adding this for cluster to temporarily supress warnings
oldw <- getOption("warn")
options(warn = -1)

# means and shit
newdf <- data.frame()
regiondata <- data.frame()
for(i in 1:length(unique(dat$dup_scaf))){
  # which one are we on?
  cat(paste0("Calculating sliding window statistics for unique duplicate by scaffold combination ", i, " out of ", numunique, "\n"))
  dupscaf <- unique(dat$dup_scaf)[i]
  tmpdf <-  filter(dat, dup_scaf == dupscaf)

  # PacBio stats
  tmpdf$PB_means <- roll_mean(tmpdf$PacBio_depth, n = 1000, by = 1, align = "right", na.rm = FALSE, fill = "NA")
  tmpdf$PB_sd <- roll_sd(tmpdf$PacBio_depth, n = 1000, by = 1, align = "right", na.rm = FALSE, fill = "NA")

  # Nanopore stats
  tmpdf$ONT_means <- roll_mean(tmpdf$ONT_depth, n = 1000, by = 1, align = "right", na.rm = FALSE, fill = "NA")
  tmpdf$ONT_sd <- roll_sd(tmpdf$ONT_depth, n = 1000, by = 1, align = "right", na.rm = FALSE, fill = "NA")

  # Illumina stats
  tmpdf$Ill_means <- roll_mean(tmpdf$Illumina_depth, n = 1000, by = 1, align = "right", na.rm = FALSE, fill = "NA")
  tmpdf$Ill_sd <- roll_sd(tmpdf$Illumina_depth, n = 1000, by = 1, align = "right", na.rm = FALSE, fill = "NA")
  newdf <- rbind(newdf, tmpdf)

  ### calculate means + sd of just the BUSCO regions of interest
  dup_start <- tmpdf$dup_start[1]
  dup_end <- tmpdf$dup_end[1]

  # pull out just region of interest
  region <- tmpdf %>% filter(base >= dup_start & base <= dup_end) %>%
    summarise(dup_scaf = unique(dup_scaf),
              BUSCO_ID = unique(BUSCO_ID),
              scaffold = unique(scaffold),
              mean_PB_depth = mean(PacBio_depth),
              sd_PB_depth = sd(PacBio_depth),
              mean_ONT_depth = mean(ONT_depth),
              sd_ONT_depth = sd(ONT_depth),
              mean_Ill_depth = mean(Illumina_depth),
              sd_Ill_depth = sd(Illumina_depth))

  # merge together
  regiondata <- rbind(regiondata, region)
}

# end suppression of warnings
options(warn = oldw)

## save output
fwrite(newdf, "slidingwindowstatistics_duplicated_BUSCO4.tsv", sep = "\t", row.names = FALSE)
fwrite(regiondata, "BUSCOregionstatistics_duplicated_BUSCO4.tsv", sep = "\t", row.names = FALSE)
