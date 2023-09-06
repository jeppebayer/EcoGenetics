library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

covfile <- read.table(args[1], header = FALSE, sep = "\t")

p <- ggplot() +
        geom_histogram(binwidth = 1, aes(x = covfile$V2)) +
        geom_vline(aes(xintercept = mean(covfile$V2)),
                linetype = "dashed", size = 1) +
        scale_x_continuous(limits = c(0, max(covfile$V2) + 10),
                breaks = seq(0, max(covfile$V2) + 10, 10)) +
        xlab("Coverage") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))

plotdir <- dirname(args[1])
plotname <- "/contig.coverage.png"
plotfile <- paste(plotdir, plotname, sep = "")

ggsave(plotfile,
        plot = p, device = "png", height = 8.5, width = 15)

# cov1 <- read.table("people/Jeppe_Bayer/steps/04_assembly/Lepidocyrtus_sp/blobtools/Lep_sp_blobtools_db.blobDB.table_arthropoda.txt",
#         header = FALSE, sep = "\t")

# p1 <- ggplot() +
#         geom_histogram(binwidth = 1, aes(x = cov1$V5)) +
#         geom_vline(aes(xintercept = mean(cov2$V5)),
#                 linetype = "dashed", size = 1) +
#         scale_x_continuous(limits = c(2, 95), breaks = seq(0, 100, 2)) +
#         scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
#         xlab("Coverage") +
#         theme_bw()

# ggsave("people/Jeppe_Bayer/steps/04_assembly/Lepidocyrtus_sp/blobtools/arthropoda_covplot.png",
#         plot = p1, device = "png", height = 8.5, width = 15)

# cov2 <- read.table("people/Jeppe_Bayer/steps/04_assembly/Isotoma_sp/blobtools/Iso_sp_blobtools_db.blobDB.table_arthropoda.txt",
#         header = FALSE, sep = "\t")

# p2 <- ggplot() +
#         geom_histogram(binwidth = 1, aes(x = cov2$V5)) +
#         geom_vline(aes(xintercept = mean(cov2$V5)),
#                 linetype = "dashed", size = 1) +
#         scale_x_continuous(limits = c(2, 95), breaks = seq(0, 100, 2)) +
#         scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
#         xlab("Coverage") +
#         theme_bw()

# ggsave("people/Jeppe_Bayer/steps/04_assembly/Isotoma_sp/blobtools/arthropoda_covplot.png",
#         plot = p2, device = "png", height = 8.5, width = 15)
