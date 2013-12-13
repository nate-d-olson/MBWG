##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%
##
## Proportion analysis and plots for biologically variant positions
## written by: Nate Olson
## May 2013
## Modified 7/20/2013 - cleanup for release on github
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%

# required packages
library(ggplot2)
library(reshape2)
library(plyr)

# set working directory for variant_analysis folder
# will need to change to the appropriate location for local computer
setwd("~/Documents/CCQM/CCQM_Analysis_Results/pileup_parse/")
pileup_csv <- list.files()[grep("pileup.csv$", c(list.files()), value = F)]

# generating a dataframe combining the pileup parse files
pileup_parse = data.frame()
for(i in pileup_csv){
  fileTable <- read.table(i, header = T, sep = ",")
  dataset <- gsub(".pileup.csv", "", i)
  id <- strsplit(dataset, "-")
  fileTable$dataset <- dataset
  fileTable$org <- id[[1]][1]
  fileTable$lab <- id[[1]][3]
  fileTable$platform <- id[[1]][2]
  fileTable$rep <- id[[1]][4]
  pileup_parse <- rbind(pileup_parse, 
			fileTable[!(fileTable$Ref %in% c("A","T","C","G")) & 
				    fileTable$Location < 1490,])
}

pileup_parse$A <- pileup_parse$A + pileup_parse$a
pileup_parse$C <- pileup_parse$C + pileup_parse$c
pileup_parse$G <- pileup_parse$G + pileup_parse$g
pileup_parse$T <- pileup_parse$T + pileup_parse$t
pileup_parse <- subset(pileup_parse, 
		       select = -c(Control,a,c,g,t,Insertion,Deletion))

# position associated ambiguities
ambigs <- list(
  "118"=c("A","G"),
  "975"=c("G","A"),
  "979"=c("G","C"),
  "983"=c("T","C"),
  "992"=c("A","G"),
  "993"=c("G","A"),
  "994"=c("A","T"),
  "995"=c("A","T"),
  "996"=c("T","G"),
  "1011"=c("C","T"),
  "1395"=c("A","G"),
  "175"=c("T","G"),
  "188"=c("C","T"),
  "419"=c("A","G"))

# bayesian analysis
max_copy = NULL
max_prob = NULL
prop = NULL
for(i in 1:nrow(pileup_parse)){
  y = pileup_parse[i,ambigs[as.character(pileup_parse$Location[i])][[1]][2]] 
  n = pileup_parse[i,ambigs[as.character(pileup_parse$Location[i])][[1]][2]] +
    pileup_parse[i,ambigs[as.character(pileup_parse$Location[i])][[1]][1]] 
  print(c(y,n))
  if(pileup_parse$org[i] == "Ecoli"){
    copy_number = 7
    p = (0:7)/7
  } else {
    copy_number = 6
    p = (0:6)/6
  }

  z = (y-n*p)/sqrt(n*p*(1-p))
  z <- replace(z,is.nan(z),-Inf) 
  post = dnorm(z)/sum(dnorm(z))
  if(y == 0){
    max_prob = c(max_prob, 1)
    max_copy = c(max_copy, 0)
    prop <- c(prop, 0)
  }else if(y == n){
    max_prob = c(max_prob, 1)
    max_copy = c(max_copy, copy_number)
    prop <- c(prop, 1)
  }else{
    max_copy <- c(max_copy, order(-post)[1]-1)
    max_prob <- c(max_prob, max(post))
    prop <- c(prop, y/n)
  }
}

pileup_parse <- data.frame(cbind(pileup_parse, prop, max_copy, max_prob))
pileup_parse$ratio[pileup_parse$Location %in% 
				c(975,979, 983,992:996,1011,175,419)] <- "high"
pileup_parse$ratio[pileup_parse$Location %in% 
				c(118,1395,188)]                      <- "even"

pileup_parse$max_copy <- as.factor(pileup_parse$max_copy)
pileup_parse$Location <- as.factor(pileup_parse$Location)


# code for figure
ggplot(pileup_parse) + 
  geom_jitter(aes(x = platform, y = prop, color = max_copy, shape  = lab),
              position = position_jitter(width = .25)) + 
  facet_grid(ratio~org, scale = "free_x") + 
  labs(x = "Platform", y = "Proportions", shape = "Participants", color = "Abundant\nBase\nCount")
