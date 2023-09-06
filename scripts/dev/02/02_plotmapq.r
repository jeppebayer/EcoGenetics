library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

mapqfile <- read.table(args[1], header = FALSE, sep = "\t")

p <- ggplot() +
        geom_histogram(binwidth = 1, aes(x = mapqfile$V1)) +
        geom_vline(aes(xintercept = mean(mapqfile$V1)),
                linetype = "dashed", size = 1) +
        scale_x_continuous(limits = c(0, max(mapqfile$V1) + 10),
                breaks = seq(0, max(mapqfile$V1) + 10, 10)) +
        xlab("MapQ") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))

plotdir <- dirname(args[1])
plotname1 <- "/"
plotname2 <- basename(plotdir)
plotname3 <- "_mapq_histogram.png"
plotfile <- paste(plotdir, plotname1, plotname2, plotname3, sep = "")

ggsave(plotfile,
        plot = p, device = "png", height = 8.5, width = 15)
