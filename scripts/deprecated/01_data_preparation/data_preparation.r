library(ggplot2)
library(dplyr)

depth <- read.table(file = "people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/intermediate/QC_sort/E100050923_L01_ANIcnqkR228563-642.markdup.depth", header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE, col.names = c("Scaffold", "Locus", "Depth"))

# Compressing the dataframe in windows
depth %>%
mutate(rounded_position = round(Locus, -2)) %>%
    group_by(rounded_position) %>%
        summarize(mean_cov = mean(Depth)) -> compressed

# Plotting the data
p <- ggplot(data =  compressed, aes(x = rounded_position, y = mean_cov)) + geom_area() + theme_classic() + ylim(0, 400)

# Saving your coverage plot
ggsave("people/Jeppe_Bayer/steps/data_preparation/Orchesella_cincta/NYS-F/intermediate/QC_sort/CoveragePlot.png", p)
