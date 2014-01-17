##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%
##
## Coverage plot for CCQM identity work
## written by: Nate Olson
## May 2013
## Modified 7/20/2013 - cleanup for release on github
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%

library(ggplot2)
library(reshape2)

# set working directory for CCQM folder with pileup.csv files (generated using ccqm_pipe.sh)
# will need to change to the appropriate location for local computer
setwd("/media/second/mirror/CCQM/Identity-Study-I/CCQM_Analysis_Results/pileup_parse/")

#######
## Import of pileup files parsed using PileUpParse.pl
## input file names need to be the following format or change the input loop below
## ORG-PLATFORM-LAB-REPLICATE.pileup.csv
#######

pileup_csv <- list.files()[grep("F.pileup.csv", c(list.files()), value = F)]

#generating a dataframe combining the pileup parse files
pileup = data.frame()
for(i in 1:length(pileup_csv)){
  fileTable <- read.csv(pileup_csv[i], header = T)
  dataset <- gsub(".pileup.csv", "", i)
  id <- strsplit(pileup_csv[i], "-")
  fileTable$org <- id[[1]][1]
  fileTable$lab <- id[[1]][3]
  fileTable$platform <- id[[1]][2]
  fileTable$data_set <- pileup_csv[i]
  pileup <- rbind(pileup, fileTable)
}
remove(fileTable, i, id, pileup_csv)
pileup <- subset(pileup, select = -c(Control))

ggplot(pileup) + 
  geom_path(aes(x = Location, 
		y = Coverage, 
		color = lab, 
		linetype = org, 
		group = data_set)) + 
  facet_wrap(~platform) + scale_y_log10() + 
  labs(y = "Coverage", x = "Reference Location (bp)", color = "Participants", linetype = "Organisms") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")




