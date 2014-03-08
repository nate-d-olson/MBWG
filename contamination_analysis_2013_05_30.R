# analysis of read contamination - plots for generating contamination figures
library(ggplot2)
library(reshape2)
library(plyr)

contam <- read.csv("~/Documents/CCQM//analysis_scripts/contam.csv")


contam$lab_id <- paste(contam$lab, paste("Rep ",contam$rep, sep = ""), sep = "\n")
# meta$percent <- meta$prop * 100
ggplot(contam[contam$template == "IRMM",]) + 
   geom_bar(aes(x = lab_id, fill = id_type), position = "fill") + 
   facet_grid(.~org*platform, scale = "free_x", space = "free") +
  # theme(axis.text.x = element_text(angle = 90)) +
   labs(x = "Datasets", y = "Identification Proportions", fill = "Identificaiton\nType")

ggplot(contam[contam$template == "IRMM" & contam$level != 1,]) + 
  geom_bar(aes(x = platform, fill = id_type), position = "fill") + 
  facet_wrap(~org) +
  labs(x = "Platform", y = "Identification Proportions", fill = "Identification\nType")

contam$level_id <- factor(contam$level, levels = c(1:6), labels = c("Domain","Phylum","Class","Order","Family", "Genus"))
ggplot(contam[contam$template == "IRMM" & contam$level != 1,]) + geom_boxplot(aes(x = platform, y = confidence, color = id_type)) + 
  facet_grid(org~level_id) + 
  labs(x = "Platform", color = "Identification\nType")

ggplot(contam[contam$template == "IRMM" & contam$level != 1,]) + geom_boxplot(aes(x = platform, y = confidence, color = id_type)) + 
  facet_wrap(~org) +
  labs(x = "Platform", color = "Identification\nType", y = "Confidence Score")