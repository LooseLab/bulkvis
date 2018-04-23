library(ggplot2)
library(dplyr)
library(tidyr)

ii <- read.delim('just_fused_dist.csv', sep=',')
ii = ii %>%
  gather('Classification', 'Time', 2:ncol(ii))
ii$type = 'C) Split Read Start'
jj <- read.delim('to_be_fused_dist.csv', sep=',')
jj = jj %>%
  gather('Classification', 'Time', 2:ncol(jj))
jj$type = 'D) Internal Read Start'
kk <- read.delim('un_fused_dist.csv', sep=',')
kk = kk %>%
  gather('Classification', 'Time', 2:ncol(kk))
kk$type = 'A) Unique Read Start'

i <- read.delim('end_just_fused_dist.csv', sep=',')
i = i %>%
  gather('Classification', 'Time', 2:ncol(i))
i$type = 'E) Internal Read Ends'
j <- read.delim('end_to_be_fused_dist.csv', sep=',')
j = j %>%
  gather('Classification', 'Time', 2:ncol(j))
j$type = 'F) Split Read End'
k <- read.delim('end_un_fused_dist.csv', sep=',')
k = k %>%
  gather('Classification', 'Time', 2:ncol(k))
k$type = 'B) Unique Read End'

df = rbind(kk,jj,ii,k, j, i)
df = filter(df, df$Time >= -10 & df$Time <= 10)

g <- ggplot(subset(df, Classification %in% c("above", "adapter", "pore", "transition",  "unblocking", "unclassified")), aes(Time,fill=Classification)) +
  geom_density(aes(y=..count..),alpha=0.5) +
  facet_grid(Classification ~ type, scales="free_y") +
  scale_x_continuous(breaks = c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  theme_minimal() +
  theme(legend.position='none',plot.title = element_text(family = "Helvetica", face = "bold", size = (15),hjust = 0.5)) +
  ggtitle("Bulk File Classifications - RAD004 Sample (97030)") +     geom_vline(xintercept=c(0), linetype="dotted") + xlab("Time (Seconds)") + ylab("Classification Count")

g