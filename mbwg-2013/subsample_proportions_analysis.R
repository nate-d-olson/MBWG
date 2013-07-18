setwd("~/Documents/CCQM/ccqm_proportions_analysis/")
library(ggplot2)
library(reshape2)
library(plyr)

#set working directory for variant_analysis folder

pileup_csv <- list.files()[grep("csv$", c(list.files()), value = F)]
#%#%#% generating a dataframe combining the pileup parse files

pileup_parse = data.frame()
for(i in pileup_csv){
  fileTable <- read.table(i, header = F, sep = ",")
    colnames(fileTable) <- c("Reference","Location", "Ref","String","Qual", "blank","depth","A","C","G","T") 
    dataset <- gsub(".pileup.csv", "", i)
    id <- strsplit(dataset, "-")
    fileTable$dataset <- dataset
    fileTable$org <- id[[1]][1]
    fileTable$lab <- id[[1]][3]
    fileTable$platform <- id[[1]][2]
    fileTable$rep <- id[[1]][4]
    fileTable$strain <- id[[1]][5]
    fileTable$trim <- id[[1]][6]
    fileTable$frac <- id[[1]][7]
    fileTable <- subset(fileTable[!(fileTable$Ref %in% c("A","T","C","G")),], select = -c(String,Qual,blank))
    pileup_parse <- rbind(pileup_parse, fileTable)
  }
}
pileup_parse <- pileup_parse[pileup_parse$Ref != "Traversal" & pileup_parse$Location %in% 100:1400,]
write.csv(pileup_parse, "pileup_counts_variant_positions_only.csv")

# reading proportion count file
pileup_parse <- read.csv("pileup_counts_variant_positions_only.csv")

pileup_parse$Location <- as.numeric(as.character(pileup_parse$Location))


# statistical analysis -bayes
#ambiguities
R <- c("A","G")
S <- c("C","G")
Y <- c("C","T")
W <- c("A","T")
K <- c("G","T")
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


library(ggplot2)
#ggplot(pileup_parse) + geom_jitter(aes(x = Location, y = prop, color = log(depth))) + facet_wrap(~org, scale = "free") 
#ggplot(pileup_parse) + geom_point(aes(depth, prop)) + xlim(0,200) + facet_wrap(~Location)
#ggplot(pileup_parse[pileup_parse$Location %in% 0:1400,]) + geom_point(aes(depth, prop)) + xlim(0,100) + facet_wrap(~Location)
#ggplot(pileup_parse) + geom_smooth(aes(depth, prop)) + scale_x_log10() + facet_wrap(~Location)

pileup_parse$max_copy <- as.factor(pileup_parse$max_copy)
ggplot(pileup_parse) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  scale_x_log10() + facet_wrap(~org*Location)
ggplot(pileup_parse[pileup_parse$depth < 200,]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  #scale_x_log10() + 
  facet_wrap(~org*Location)

ggplot(pileup_parse[pileup_parse$depth < 200 & pileup_parse$org == "Ecoli",]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  facet_wrap(~org*Location)

ggplot(pileup_parse[pileup_parse$depth < 200 & pileup_parse$org == "Lmono",]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  facet_wrap(~org*Location)

pileup_parse$exp[pileup_parse$Location %in% c(975,979,983,992,993,994,995,996,1011)] <- "6:1" 
pileup_parse$exp[pileup_parse$Location %in% c(118,1395)] <- "4:3"
pileup_parse$exp[pileup_parse$Location %in% c(118,1395)] <- "4:3"
pileup_parse$exp[pileup_parse$Location %in% c(175,419)] <- "5:1"
pileup_parse$exp[pileup_parse$Location == 188] <- "3:3"

ggplot(pileup_parse) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  scale_x_log10() +
  facet_wrap(~org*exp)

ggplot(pileup_parse[pileup_parse$depth < 200,]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  facet_wrap(~org*exp)

high <- c("6:1","5:1")
even <- c("4:3","3:3")

ggplot(pileup_parse[pileup_parse$exp %in% high,]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  scale_x_log10() + geom_vline(aes(xintercept = 100), linetype = 2, alpha = 0.5) + facet_wrap(~exp) +
  labs (x = "log(Depth)", y = "Observed Proportions", shape = "Sequencing\nPlatform", color = "Copy Ratio")

ggplot(pileup_parse[pileup_parse$exp %in% even,]) + geom_point(aes(depth, prop, color = max_copy, shape = platform)) + 
  scale_x_log10() + geom_vline(aes(xintercept = 200), linetype = 2, alpha = 0.5) + facet_wrap(~exp) +
  labs (x = "log(Depth)", y = "Observed Proportions", shape = "Sequencing\nPlatform", color = "Copy Ratio")
