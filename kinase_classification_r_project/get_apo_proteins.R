###splits full.data into apo and ligand bound sets

source("train_random_forest.R")
library(tidyverse)
library(plyr)

meta = read.table("../../stdy_kinase.xtal.complete.170601.csv", header = T, sep =",")
meta.ligand =meta[which(meta$ligand!="NaN"),]
meta.apo = meta[which(meta$ligand=="NaN"),]
full.data$row.names = row.names(full.data)
full.data[1,18] = "1ATP_A"
full.data$pdb = ldply(strsplit(full.data$row.names, "_"))$V1

full.data.ligand = full.data[which(full.data$pdb %in% meta.ligand$pdb_id),]
full.data.apo = full.data[which(full.data$pdb %in% meta.apo$pdb_id),]
full.data.ligand$type = rep("ligand", nrow(full.data.ligand))
full.data.apo$type = rep("apo", nrow(full.data.apo))
full.data.2 = as.data.frame(rbind(full.data.apo, full.data.ligand))
head(full.data)

write.table(full.data.ligand, "../../2_predicted_classes/8.29.17.predicted+verified.ligand.csv" , quote =F, eol = "\n" ,col.names = T  )
write.table(full.data.apo, "../../2_predicted_classes/8.29.17.predicted+verified.apo.csv" , quote =F, eol = "\n" ,col.names = T  )
