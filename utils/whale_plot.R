#! /usr/bin/Rscript
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))

args = commandArgs(trailingOnly=TRUE)

A = args[1]
B = args[2]
C = args[3]
D = args[4]
E = args[5]
F = args[6]
out = args[7]
run = args[8]

ii <- read.delim(file = C, sep=',')
ii = ii %>%
  gather('Classification', 'Time', 2:ncol(ii))
ii$type = 'C) Split Read Start'
jj <- read.delim(file = D, sep=',')
jj = jj %>%
  gather('Classification', 'Time', 2:ncol(jj))
jj$type = 'D) Internal Read Start'
kk <- read.delim(file = A, sep=',')
kk = kk %>%
  gather('Classification', 'Time', 2:ncol(kk))
kk$type = 'A) Unique Read Start'

i <- read.delim(file = E, sep=',')
i = i %>%
  gather('Classification', 'Time', 2:ncol(i))
i$type = 'E) Internal Read Ends'
j <- read.delim(file = F, sep=',')
j = j %>%
  gather('Classification', 'Time', 2:ncol(j))
j$type = 'F) Split Read End'
k <- read.delim(file = B, sep=',')
k = k %>%
  gather('Classification', 'Time', 2:ncol(k))
k$type = 'B) Unique Read End'

df = rbind(kk, jj, ii, k, j, i)
df = filter(df, df$Time >= -10 & df$Time <= 10)

if (is.na(run) || is.null(run)){
  title = "Bulk fast5 file classification"
} else {
  title = paste("Bulk fast5 file classification -", run, collapse="")
}


g <- ggplot(subset(df, Classification %in% c("above", "adapter", "pore", "transition",  "unblocking", "unclassified")), aes(Time,fill=Classification)) +
  geom_density(aes(y=..count..),alpha=0.5) +
  facet_grid(Classification ~ type, scales="free_y") +
  scale_x_continuous(breaks = c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  theme_minimal() +
  theme(legend.position='none',plot.title = element_text(family = "Helvetica", face = "bold", size = (15),hjust = 0.5)) +
  ggtitle(title) +
  geom_vline(xintercept=c(0), linetype="dotted") +
  xlab("Time (Seconds)") +
  ylab("Classification Count")

ggsave(filename = out, plot = g, width = 345, height = 272, units = "mm")