setwd("~/Documents/CCQM/lmono_full_length_analysis/")
library(reshape2)
library(ggplot2)

# input files are parsed sam files using SAM_parse.py 
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

#modification of dataset for maximum likelihood analysis
unfiltered <- unique(biovars_df$dataset)[grep("F.csv", unique(biovars_df$dataset), value = F)]
biovars_df_unfilt <- biovars_df[biovars_df$dataset %in% unfiltered,]
biovars_cast <- dcast(biovars_df_unfilt, variants~dataset)

strings <- c("GTA", "GTG", "GCA","GCG","TTA", "TTG", "TCA","TCG")
write.csv(biovars_cast[biovars_cast$variants %in% strings,], "~/Desktop/copy_proportions.csv")
