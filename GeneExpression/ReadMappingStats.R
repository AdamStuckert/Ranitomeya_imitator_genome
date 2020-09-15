#/usr/bin/env/R

# ReadMappingStats
# Author: Adam Stuckert


## Assumptions: You have produced a file called "mappingdata.tab" which is in the same directory that you called this script from.
## This file (mappingdata.tab) is produced by another script I wrote (STARdata.sh), this Rscript is incorporated into that.

## import data
dat <- read.delim("mappingdata.tab", sep = "\t", header = TRUE)
# remove all "%"
dat$Uniquely_mapping_reads <- as.numeric(gsub(pattern = "%", replacement = "", dat$Uniquely_mapping_reads))

## Calculate read summary statistics
ave_reads = round(mean(dat$Number_reads), 2)
sd_reads = round(sd(dat$Number_reads), 2)
med_reads = round(median(dat$Number_reads), 2)

## Print read summary statistics
cat(paste0("Average number of reads: ", ave_reads,"\n"))
cat(paste0("Standard deviation of the number of reads: ", sd_reads,"\n"))
cat(paste0("Median number of reads: ", ave_reads,"\n"))

## Calculate mapping summary statistics
ave_map = round(mean(dat$Uniquely_mapping_reads), 2)
sd_map = round(sd(dat$Uniquely_mapping_reads), 2)
med_map = round(median(dat$Uniquely_mapping_reads), 2)

## Print mapping summary statistics
cat(paste0("Average mapping rate: ", ave_map,"%\n"))
cat(paste0("Standard deviation of the mapping rate: ", sd_map,"%\n"))
cat(paste0("Median mapping rate: ", med_map,"%\n"))

