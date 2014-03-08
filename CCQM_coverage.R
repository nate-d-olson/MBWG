##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%
##
## Coverage plot for CCQM identity work
##
##
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%


library(ggplot2)
library(reshape2)

#set working directory for CCQM folder with pileup.csv, .tsv, and vcf files
setwd("~/Documents/CCQM/CCQM_Analysis_Results/pileup_parse/")

#######
## Import of pileup files parsed using PileUpParse.pl
## input file names need to be the following format or change the input loop below
## ORG_LAB_PLATFORM_.pileup
#######

pileup_csv <- list.files()[grep("+IRMM-F.pileup.csv", c(list.files()), value = F)]

#generating a dataframe combining the pileup parse files
pileup = data.frame()
for(i in 1:length(pileup_csv)){
  fileTable <- read.csv(pileup_csv[i], header = T)
  id <- strsplit(pileup_csv[i], "-")
  fileTable$org <- id[[1]][1]
  fileTable$lab <- id[[1]][3]
  fileTable$platform <- id[[1]][2]
  fileTable$rep <- id[[1]][4]
  fileTable$strain <- id[[1]][5]
  #fileTable$trim <- id[[1]][6]
  fileTable$data_set <- pileup_csv[i]
  pileup <- rbind(pileup, fileTable)
}
remove(fileTable, i, id, pileup_csv)
pileup <- subset(pileup, select = -c(Control))

ggplot(pileup) + 
  geom_line(aes(x = Location, y = Coverage, color = platform, linetype = lab, group = data_set)) + 
  facet_wrap(~org)
ggplot(pileup[pileup$trim != "T.pileup.csv" & pileup$data_set != "Ecoli-ION-NMIC-1-UNK-F.pileup.csv", ]) + 
  geom_path(aes(x = Location, y = Coverage, color = lab, linetype = lab, group = data_set)) + 
  facet_grid(platform~org) + scale_y_log10() + 
  labs(y = "Coverage (log10)", x = "Reference Location (bp)")

ggplot(pileup[pileup$trim != "T.pileup.csv" & pileup$data_set != "Ecoli-ION-NMIC-1-UNK-F.pileup.csv", ]) + 
  geom_line(aes(x = Location, y = Coverage, color = lab, linetype = lab, group = data_set)) + 
  facet_grid(platform~org) + scale_y_log10() + 
  labs(y = "Coverage (log10)", x = "Reference Location (bp)")



