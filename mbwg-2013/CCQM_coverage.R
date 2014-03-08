##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%
##
## Coverage plot for CCQM identity work
##
##
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%


library(ggplot2)
library(reshape2)
library(plyr)

#set working directory for CCQM folder with pileup.csv, .tsv, and vcf files
setwd("/Volumes/UNTITLED/CCQM/ccqm_proportions_analysis/")

#######
## Import of pileup files parsed using PileUpParse.pl
## input file names need to be the following format or change the input loop below
## ORG_LAB_PLATFORM_.pileup
#######

pileup_csv <- list.files()[grep("-1.pileup.csv", c(list.files()), value = F)]

#generating a dataframe combining the pileup parse files
pileup = data.frame()
for(i in 1:length(pileup_csv)){
  fileTable <- read.csv(pileup_csv[i], header = F, stringsAsFactors=F)
  colnames(fileTable) <- c("Control","Location","ref","string","cigar","blank","Coverage","A","C","G","T")
  id <- strsplit(pileup_csv[i], "-")
  fileTable$org <- id[[1]][1]
  fileTable$lab <- id[[1]][3]
  fileTable$platform <- id[[1]][2]
  fileTable$rep <- id[[1]][4]
  fileTable$data_set <- pileup_csv[i]
  pileup <- rbind(pileup, fileTable)
}
remove(fileTable, i, id, pileup_csv)
pileup <- subset(pileup, select = -c(Control))

pileup$Location <- as.numeric(pileup$Location)
ggplot(pileup) + 
  geom_path(aes(x = Location, y = Coverage, color = lab, 
                linetype = org, group = data_set)) + 
  facet_wrap(~platform) + scale_y_log10() + 
  labs(y = "Coverage (log10)", x = "Reference Location (bp)")

ggplot(pileup[pileup$platform == "454",]) + 
  geom_path(aes(x =Location, y= Coverage, color = lab, 
                linetype = org, group = data_set)) +
  scale_y_log10() + 
  labs(y = "Coverage (log10)", x = "Reference Location (bp)")
