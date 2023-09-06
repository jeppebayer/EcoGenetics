args <- commandArgs(trailingOnly = TRUE)

df <- read.table(args[1])

png(paste(args[1], ".png", sep = ""), width = 800, height = 400, units = "px")
plot(df[2:100, ], type = "l")
points(df[2:100, ])
dev.off()

number_of_fields <- (nrow(df) - 1)

total_kmer_count <- (sum(as.numeric(df[2:number_of_fields, 1] * df[2:number_of_fields, 2])))

genome_size <- (sum(as.numeric(df[2:10000, 1] * df[2:10000, 2])) / as.numeric(df$V1[df$V2 == max(df[2:number_of_fields, 2])])) / 1000000

numbers <- data.frame(Number_of_fields = c(number_of_fields),
                        Total_number_of_Kmers = c(total_kmer_count),
                        Genome_size = c(paste(genome_size, "Mb")))

write.table(numbers, paste(args[1], ".stats", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")