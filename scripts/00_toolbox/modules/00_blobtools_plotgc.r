library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

file <- read.table(args[1], skip = 11, header = FALSE, sep = "\t")

if (length(args) > 1) {
    p <- ggplot() +
        geom_histogram(binwidth = 0.005, aes(x = file$V3)) +
        geom_vline(aes(xintercept = mean(file$V3)),
                linetype = "dashed", size = 1) +
        geom_vline(aes(xintercept = as.numeric(args[2])),
                linetype = "dashed", size = 1, color = "red") +
        geom_vline(aes(xintercept = as.numeric(args[3])),
                linetype = "dashed", size = 1, color = "red") +
        scale_x_continuous(limits = c(0.1, max(file$V3) + 0.01),
                breaks = seq(0.1, max(file$V3) + 0.1, 0.01)) +
        ylab("Count") +
        xlab("GC Content") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
} else {
    p <- ggplot() +
        geom_histogram(binwidth = 0.005, aes(x = file$V3)) +
        geom_vline(aes(xintercept = mean(file$V3)),
                linetype = "dashed", size = 1) +
        scale_x_continuous(limits = c(0.1, max(file$V3) + 0.01),
                breaks = seq(0.1, max(file$V3) + 0.1, 0.01)) +
        ylab("Count") +
        xlab("GC Content") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
}

plotdir <- dirname(args[1])
plotname1 <- "/"
plotname2 <- "gc_histogram.png"
plotfile <- paste(plotdir, plotname1, plotname2, sep = "")

ggsave(plotfile,
        plot = p, device = "png", height = 8.5, width = 15)