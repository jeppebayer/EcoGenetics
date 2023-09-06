library(ggplot2)


sfs <- t(read.csv("people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_villosa?/FÅJ-C5/E100049659_L01_ANIcnqkR228559-607.spectrum", sep = " ", header =FALSE))

ggplot() +
geom_point(aes(x = seq(1, 51), y = sfs[1:51])) +
ylim(0, 0.0001)

ggplot() +
geom_bar(aes(y = seq(1, 51), x = sfs[1:51]),
stat = "identity", color = "black", fill = "black") +
xlim(0.00002, 0.00008)
