##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%
##
## subsampled proportion analysis and plots for subsampled biologically variant positions
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
#setwd("~/Documents/mirror/CCQM/Identity-Study-I/ccqm_biovar_analysis/")


#%^)(*&^%$%^&*()(*&^% 
#$%^&*()(*&^&*()*&^%^&*()(*&^   temporary - will change before saving 1/6/2014



setwd("~/Desktop/to_file/ccqm_biovar_analysis/")
pileup_csv <- list.files()[grep(".pileup.csv$", c(list.files()), value = F)]

#%#%#% generating a dataframe combining the pileup parse files

pileup_parse = data.frame()
for(i in pileup_csv){
  print(i)
  fileTable <- read.table(i, header = F, sep = ",")
    colnames(fileTable) <- c("Reference","Location", "Ref","String","Qual", "blank","depth","A","C","G","T") 
    dataset <- gsub(".pileup.csv", "", i)
    id <- strsplit(dataset, "-")
    fileTable$dataset <- dataset
    fileTable$org <- id[[1]][1]
    fileTable$lab <- id[[1]][3]
    fileTable$platform <- id[[1]][2]
    fileTable$rep <- id[[1]][4]
    fileTable$frac <- id[[1]][7]
    fileTable <- subset(fileTable[!(fileTable$Ref %in% c("A","T","C","G")),], select = -c(String,Qual,blank))
    pileup_parse <- rbind(pileup_parse, fileTable)
}


pileup_parse <- pileup_parse[pileup_parse$Ref != "Traversal" & pileup_parse$Location %in% 100:1400,]
pileup_parse$Location <- as.numeric(as.character(pileup_parse$Location))


# list of biologically variant positions and potential alleles
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
  n = pileup_parse$depth[i]
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

pileup_parse$max_copy <- as.factor(pileup_parse$max_copy)

pileup_parse$exp[pileup_parse$Location %in% c(975,979,983,992,993,994,995,996,1011)] <- "6:1" 
pileup_parse$exp[pileup_parse$Location %in% c(118,1395)] <- "4:3"
pileup_parse$exp[pileup_parse$Location %in% c(118,1395)] <- "4:3"
pileup_parse$exp[pileup_parse$Location %in% c(175,419)] <- "5:1"
pileup_parse$exp[pileup_parse$Location == 188] <- "3:3"

high <- c("6:1","5:1")
even <- c("4:3","3:3")

ggplot(pileup_parse[pileup_parse$exp %in% high,]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  scale_x_log10() + geom_vline(aes(xintercept = 100), linetype = 2, alpha = 0.5) + facet_wrap(~exp) +
  labs (x = "log(Depth)", y = "Observed Proportions", shape = "Sequencing\nPlatform", color = "Copy Ratio")

ggplot(pileup_parse[pileup_parse$exp %in% even,]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  scale_x_log10() + geom_vline(aes(xintercept = 200), linetype = 2, alpha = 0.5) + facet_wrap(~exp) +
  labs (x = "log(Depth)", y = "Observed Proportions", shape = "Sequencing\nPlatform", color = "Copy Ratio")

library("Hmisc")
pileup_parse$depth_groups <- as.numeric(cut2(pileup_parse$depth))
library("plyr")
ggplot(pileup_parse) + geom_histogram(aes(x = depth_groups)) + facet_wrap(~exp)
depth_decay <- ddply(pileup_parse,.(exp, depth_groups), summarise, decay = sd(prop), depth = mean(depth), size = length(platform))
ggplot(depth_decay) + geom_point(aes(x = depth, y= decay, color = exp)) + 
  geom_smooth(aes(x = depth, y= decay, color = exp)) + scale_x_log10() + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(y = "Standard Deviation of Observed Proportion", x = "Coverage", color = "Expected Ratio")

ddecay <- ddply(pileup_parse, .(exp), depth_groups = as.numeric(cut(depth, 30))) 
decay_sd <- ddply(pileup_parse, .(exp, depth_groups), summarize, depth = median(depth), decay = sd(prop))
ggplot(decay_sd) + geom_point(aes(x = depth, y= decay)) + 
  geom_smooth(aes(x = depth, y= decay)) + scale_x_log10() + 
  theme_bw() + 
  labs(y = "Standard Deviation of Observed Proportion") +
  facet_wrap(~exp)