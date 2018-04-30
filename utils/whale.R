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

c <- read.delim(file = C, sep=',')
c = c %>%
  gather('Classification', 'Time', 2:ncol(c))
c$type = 'C) Split Read Start'
d <- read.delim(file = D, sep=',')
d = d %>%
  gather('Classification', 'Time', 2:ncol(d))
d$type = 'D) Internal Read Start'
a <- read.delim(file = A, sep=',')
a = a %>%
  gather('Classification', 'Time', 2:ncol(a))
a$type = 'A) Unique Read Start'

e <- read.delim(file = E, sep=',')
e = e %>%
  gather('Classification', 'Time', 2:ncol(e))
e$type = 'E) Internal Read Ends'
f <- read.delim(file = F, sep=',')
f = f %>%
  gather('Classification', 'Time', 2:ncol(f))
f$type = 'F) Split Read End'
b <- read.delim(file = B, sep=',')
b = b %>%
  gather('Classification', 'Time', 2:ncol(b))
b$type = 'B) Unique Read End'

df = rbind(a, b, c, d, e, f)
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