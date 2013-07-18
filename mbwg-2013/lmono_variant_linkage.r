setwd("~/Documents/CCQM/lmono_full_length_analysis/")
library(reshape2)
library(ggplot2)


sam_csv <- list.files()[grep("+.csv", c(list.files()), value = F)]

#generating a dataframe combining the pileup parse files
biovars_df <- data.frame()

for(i in sam_csv){
  sam_biovar <- read.csv(i)
  sam_biovar <- sam_biovar[sam_biovar$ref_position < 500,]
  biovars <- dcast(sam_biovar, read_name~ref_position, value.var = "base_call")
  biovars <- biovars[!(is.na(biovars$"418") | is.na(biovars$"174")),]
  biovars$variants <- paste(biovars$"174",biovars$"187",biovars$"418", sep = "")
  biovars$dataset <- i
  biovars_df <- rbind(biovars_df, biovars)
}
remove(i, biovars, sam_biovar,sam_csv)
pdf("Lmono_variant_string_analysis.pdf", width = 14, height = 8)
var_strings <- c("GTA","GTG","GCA","GCG","TCA","TCG")
ggplot(biovars_df[biovars_df$variants %in% var_strings,]) + geom_bar(aes(x = variants, fill = dataset), position = "dodge")
ggplot(biovars_df[biovars_df$variants %in% var_strings,]) + geom_bar(aes(x = variants)) + facet_wrap(~dataset, scale = "free_y")
dev.off()

#modification of dataset for Steve
#need to remove filtered datasets
unfiltered <- unique(biovars_df$dataset)[grep("F.csv", unique(biovars_df$dataset), value = F)]
biovars_df_unfilt <- biovars_df[biovars_df$dataset %in% unfiltered,]
steve <- dcast(biovars_df_unfilt, variants~dataset)

strings <- c("GTA", "GTG", "GCA","GCG","TTA", "TTG", "TCA","TCG")
steve2 <- steve[steve$variants %in% strings,]
write.csv(steve2, "~/Desktop/copy_proportions.csv")
