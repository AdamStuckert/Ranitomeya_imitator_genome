#/usr/bin/env/R

# Fig4_Repeats.R
# Author: Adam Stuckert


## Assumptions: You have a single column of data with read lengths (in base pairs).
## Purpose: This script will make nice figures of genome scaffolds + lengths, etc.
## Requirements: ggplot2, cowplot
# load libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)
#library(scattermore)

##### Repeats figure
repdat <- read.delim("data/Repeats.csv", sep = ",", header = TRUE)

# Compute the position of labels for pie chart
repdat <- repdat %>% 
  arrange(desc(Repeat_type)) %>%
  mutate(ypos = cumsum(Proportion)- 0.5*Proportion )

# order factors
DNA_order <- c('LINEs', 'LTRs', 'DNA transposons', 'Simple repeats', 'Other', 'Unclassified', 'Not classified as repeats') #this vector might be useful for other plots/analyses

# make pie chart
pie <- ggplot(repdat, aes(x="", y=Proportion, fill=factor(Repeat_type, level = DNA_order)))+
  geom_bar(width = 1, stat = "identity", color="white") + 
  coord_polar("y", start=0) +
  ylab("Proportion of genome") +
  xlab("") +
  theme_void() + 
  theme(legend.title = element_blank(), legend.text=element_text(size=13)) 
pie  
#geom_text(aes(y = ypos, label = Repeat_type), color = "white", size=4) 
#pie

ggsave("figures4publication/Fig4_RepeatPiechart.png", width = 7, height = 7, dpi = 600)
