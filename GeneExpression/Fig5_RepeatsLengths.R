#/usr/bin/env/R

# Fig5_RepeatLengths.Rmd
# Author: Adam Stuckert

library(ggplot2)
library(tidyverse)
library(data.table)
dir.create("results")
dir.create("results4supplemental")
#library(RcppRoll)

### Get in kids, we are looking at some repeats!



# note about these data. They are in a file produced by repeatmasker called `imitator.1.3.6.fa.out`

#!/bin/bash
# awk '{print $5, $6, $7, $11}' R_imi_1.0.fa.out > R_imi_1.0.certaincolumns.out
# I then removed the weird header lines that are funky. 
#wc -l R_imi_1.0.certaincolumns.out | cut -f1 -d " " ## gets total lines
#sed -n '4,13634563p' R_imi_1.0.certaincolumns.out > tmp
#mv tmp R_imi_1.0.certaincolumns.out



# import data
reps <- fread("data/R_imi_1.0.certaincolumns.out")

# Rename columns
colnames(reps) <- c("scaffold", "repeat_begin", "repeat_end", "repeat_type")

# calculate length of repeat
reps <- dplyr::mutate(reps, length = repeat_end - repeat_begin)

# summariiiiiieeeees
reps_sum <- reps %>% group_by(repeat_type) %>% summarise(mean_length = mean(length), sd = sd(length), number_of_repeats = n(), minimum_length = min(length), maximum_length = max(length), total_bp = sum(length), Gbp = sum(length)/1000000000)

# total genomic content in repeats here
print(paste0("Total repeat content in the genome from these major families (in Gbp): ", sum(reps_sum$Gbp)))

# save these summary data
write.csv(reps_sum, "results/RepeatSummaryStats.csv", row.names = FALSE)
write.csv(reps_sum, "results4supplemental/RepeatSummaryStats.csv", row.names = FALSE)

# drop repeats with very little data (here, < 50Mbp )

reps_sum2 <- reps_sum %>% filter(total_bp > 50000000)


# Plot just the summary data!
  
  

sum_plot <- ggplot(reps_sum2, aes(x = repeat_type, y = mean_length)) +
  geom_point() + 
  geom_errorbar(aes(ymin = mean_length - sd, ymax = mean_length + sd)) +
  ylab("Length (bp)") +
  xlab("Repeat type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
sum_plot


# Now pull out only repeats that are one of the 26 types with >50Mbp of repeats

reps_subset <- reps %>% filter(repeat_type %in% reps_sum2$repeat_type)

# violin plot of this...
reps_plot <- ggplot(reps_subset, aes(x = repeat_type, y = length)) +
  geom_point(alpha = 0.05, position = "jitter", size = 0.1, pch = ".") +  geom_violin(size = 1, color = "red", fill = NA) +
  ylab("Length (bp)") +
  xlab("Repeat type") + 
  scale_y_log10(labels = scales::comma, breaks = c(10, 100, 1000, 10000, 100000)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
reps_plot


ggsave("figures4publication/Fig_5RepeatLengths.png", width = 7.4, height = 4, dpi = 600)


