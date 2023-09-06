##ggplot2 plotting nuclear data (bold_github environment)
library(ggplot2)
library(tidyverse)
library(ggpubr)

cov<-read.table("mapping/STATS_mem/out.cov", sep="\t")
covmarkdup<-read.table("mapping/STATS_mem/markdup.cov",sep="\t")

p1 <- ggplot(cov, aes(x=cov$V1, y=cov$V2)) + 
xlim(0,1001) +
ylim(0,14000000) +
labs(x = "Coverage",y = "No. of mapped bp") +
ggtitle("Cleaned") +
geom_point(color="blue") +
theme(axis.title.x = element_text(size=14)) +
theme(axis.title.y = element_text(size=14)) +
theme_bw() +
theme(plot.title = element_text(face="italic", hjust = 0.5, size=14),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), 
axis.line = element_line(colour = "black")) #+
#scale_y_log10()

p3 <- ggplot(covmarkdup, aes(x=covmarkdup$V1, y=covmarkdup$V2)) + 
xlim(0,1001) +
ylim(0,14000000) +
labs(x = "Coverage",y = "No. of mapped bp") +
ggtitle("markdup") +
geom_point(color="lightgreen") +
theme(axis.title.x = element_text(size=14)) +
theme(axis.title.y = element_text(size=14)) +
theme_bw() +
theme(plot.title = element_text(face="italic", hjust = 0.5, size=14),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), 
axis.line = element_line(colour = "black")) #+
#scale_y_log10()

#Plotting read lengths
RL<-read.table("mapping/STATS_mem/out.RL", sep="\t")
markdup_RL<-read.table("mapping/STATS_mem/markdup.RL",sep="\t")

p2 <- ggplot(RL, aes(x=RL$V1, y=RL$V2)) + 
xlim(0,200) +
ylim(0,8000000) +
labs(x = "Read length",y = "Frequency") +
ggtitle("Cleaned") +
geom_point(color="blue") +
theme(axis.title.x = element_text(size=14)) +
theme(axis.title.y = element_text(size=14)) +
theme_bw() +
theme(plot.title = element_text(face="italic", hjust = 0.5, size=14),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), 
axis.line = element_line(colour = "black")) #+
#scale_y_log10()

p4 <- ggplot(markdup_RL, aes(x=markdup_RL$V1, y=markdup_RL$V2)) + 
xlim(0,200) +
ylim(0,8000000) +
labs(x = "Read length",y = "Frequency") +
ggtitle("markdup") +
geom_point(color="lightgreen") +
theme(axis.title.x = element_text(size=14)) +
theme(axis.title.y = element_text(size=14)) +
theme_bw() +
theme(plot.title = element_text(face="italic", hjust = 0.5, size=14),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), 
axis.line = element_line(colour = "black")) #+
#scale_y_log10()

ggarrange(p1, p2, p3, p4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
#ggsave("Raw_Coverage.pdf", width=15, height=8.5)
ggsave("Cov_RL_cleaned_and_markdup.pdf", width=15, height=8.5)