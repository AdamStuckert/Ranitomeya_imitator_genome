#/usr/bin/env/R

# ReadMappingStats.R
# Author: Adam Stuckert


## Assumptions: You have a single column of data with read lengths (in base pairs).
## Purpose: This script will make nice figures of genome scaffolds + lengths, etc.
## Requirements: ggplot2, cowplot
# load libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
# import data
dat <- read.delim("data/imitator.1.3.6.fa.sequencelengths.txt", sep = "\t", header = FALSE)

# order by length
dat$V1 <- sort(dat$V1, decreasing = TRUE)

# add column for scaffold number
scaf_num <- seq(1, length(dat$V1), 1)
dat <- cbind(scaf_num, dat)
colnames(dat) <- c("Scaffold_number", "Scaffold_length")

## Cumulative genome size
# First calculate running total
totalbases = 0
numbases = sum(dat$Scaffold_length)
dat$Total_bases = NA
dat$Proportion_of_genome = NA

for (i in 1:length(dat$Scaffold_length)){
  totalbases = totalbases + dat$Scaffold_length[i]
  dat$Total_bases[i] = totalbases
  dat$Proportion_of_genome[i] = totalbases/numbases
}

## calculate N50
tmp <- dat %>% filter(Total_bases > (numbases / 2)) 
N50 <- tmp[1,]

## calculate L50

#### Plotting

## Scaffold length plot
# to make the plot easier to 
scaflen <- ggplot(dat, aes(Scaffold_length)) +
  geom_freqpoly(bins = 150) +
  #geom_area(aes(y = ..count..), stat = "bin") +
  xlab("Scaffold length") +
  ylab("Number of scaffolds") +
  theme_classic()  +
  scale_x_log10(labels = scales::comma, breaks = c(1000,10000,100000,1000000)) + 
  scale_y_continuous(labels = scales::comma) + 
  coord_cartesian(expand = FALSE, xlim = c(1000,10000000), ylim = c(0,2000)) +
  geom_vline(xintercept = N50$Scaffold_length, color = "red") # adds N50
scaflen  



# add a proportion column
#numbases = dat$Total_bases
#dat$Proportion_of_genome = (dat$Total_bases/numbases)

# Plot cumulative total
propgenome <- ggplot(dat, aes(x = Scaffold_number, y = Proportion_of_genome)) +
  geom_line() +
  xlab("Total scaffolds") +
  ylab("Proportion of genome") +
  theme_classic() + 
  scale_x_log10(labels = scales::comma, breaks = c(100,1000,10000,75000)) + 
  scale_y_continuous(labels = scales::comma) + 
  coord_cartesian(expand = FALSE, ylim = c(0,1),xlim = c(1,100000)) +
  geom_vline(xintercept = N50$Scaffold_number, color = "red") # adds L50
propgenome


## Combine plots and save:
#grid.arrange(scaflen, propgenome, ncol=1)
#gridplot <- marrangeGrob(scaflen, propgenome, ncol=1, nrow = 1)
#ggsave("figures4publication/GenomeFig.png", gridplot, width = 8.69, height = 5.9, dpi = 600)



plot_grid(scaflen, propgenome, labels=c("A", "B"), ncol = 1, nrow = 2)


## Violin plot of genome size relative to predicted genes
# modify scaffold df to list type
dat$seq_type <- "scaffold"
scafdat <- dat[,c(2,5)]
colnames(scafdat)[1] <- "length"
# import lengths from predicted transcripts
txdat <- read.delim("data/Ranitomeya_imitator.imitator.1.3.6.functional.transcripts.fasta.sequencelengths.txt", sep = "\t", header = FALSE)
# modify transcript df to list type
colnames(txdat)[1] <- "length"

#### Calculate L50 and N50 for tx
# order by length
txdat$length <- sort(txdat$length, decreasing = TRUE)

# add column for scaffold number
tx_num <- seq(1, length(txdat$length), 1)
txdat <- data.frame(cbind(tx_num, txdat$length))
colnames(txdat)[2] <- "length"
txdat$seq_type = "transcript"
## Cumulative genome size
# First calculate running total
totalbases = 0
numbases = sum(txdat$length)
txdat$Total_bases = NA
txdat$Proportion_of_genome = NA

for (i in 1:length(txdat$length)){
  totalbases = totalbases + txdat$length[i]
  txdat$Total_bases[i] = totalbases
  txdat$Proportion_of_genome[i] = totalbases/numbases
}

## calculate N50
tmp <- txdat %>% filter(Total_bases > (numbases / 2)) 
txN50 <- tmp[1,]


## scaffold violin plot
scafviol <- ggplot(dat, aes(x = seq_type, y = Scaffold_length)) +
  geom_violin(size = 1.5) + geom_point(alpha = 0.1, position = "jitter", size = 0.25) + 
  scale_y_log10(labels = scales::comma, breaks = c(1000,10000,100000,1000000,10000000,10000000,100000000,100000000)) +
  ylab("Scaffold length") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  geom_hline(yintercept = N50$Scaffold_length, color = "red") # adds N50
scafviol


## transcript violin plot
txviol <- ggplot(txdat, aes(x = seq_type, y = length)) +
  geom_violin(size = 1) + geom_point(alpha = 0.1, position = "jitter", size = 0.25) + 
  scale_y_log10(labels = scales::comma, breaks = c(1,10,100,1000,10000)) +
  ylab("Transcript length") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  geom_hline(yintercept = txN50$length, color = "red") # adds N50
txviol


bottom_row <- plot_grid(scafviol, txviol, labels = c('C', 'D'), ncol = 2, nrow = 1)


plot_grid(scaflen, propgenome, bottom_row, labels=c("A", "B", "", ""), ncol = 1, nrow = 3)
ggsave("figures4publication/GenomeLengthFigure.png", width = 8.69, height = 7, dpi = 600)



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

ggsave("figures4publication/Piechart.png", width = 7, height = 7, dpi = 600)

