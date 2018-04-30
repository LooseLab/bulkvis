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

col_c <- read.delim(file = C, sep=',')
col_c = col_c %>%
  gather('Classification', 'Time', 2:ncol(col_c))
col_c$type = 'C) Split Read Start'
col_d <- read.delim(file = D, sep=',')
col_d = col_d %>%
  gather('Classification', 'Time', 2:ncol(col_d))
col_d$type = 'D) Internal Read Start'
col_a <- read.delim(file = A, sep=',')
col_a = col_a %>%
  gather('Classification', 'Time', 2:ncol(col_a))
col_a$type = 'A) Unique Read Start'

col_e <- read.delim(file = E, sep=',')
col_e = col_e %>%
  gather('Classification', 'Time', 2:ncol(col_e))
col_e$type = 'E) Internal Read Ends'
col_f <- read.delim(file = F, sep=',')
col_f = col_f %>%
  gather('Classification', 'Time', 2:ncol(col_f))
col_f$type = 'F) Split Read End'
col_b <- read.delim(file = B, sep=',')
col_b = col_b %>%
  gather('Classification', 'Time', 2:ncol(col_b))
col_b$type = 'B) Unique Read End'

df = rbind(col_a, col_b, col_c, col_d, col_e, col_f)
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