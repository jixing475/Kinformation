tani = read.table("../kinase_metrics.v5/../../pdb.ligands.smi.tanimoto.fp.txt.uniq.cases", sep = ",", header = F)
lig = read.table("../../pdb.ligand.smi.csv.complete.cases", sep = ",", header = F)
library(qgraph)
library(tidyr)
lig$V4= NULL
lig$V3= NULL
colnames(lig) = c("Row.names", "Ligand")
full.data$Row.names = rownames(full.data)
full.data.lig = merge(full.data, lig, by = "Row.names")
lig.conf = full.data.lig[,c(20,19)]
head(lig.conf)
lig.conf.counted = as.data.frame(table(lig.conf))
lig.conf.counted.wide = spread(lig.conf.counted, V18 ,value =  Freq , fill = 0 )
ligand.dominant.class = data.frame()
for (i in 1:nrow(lig.conf.counted.wide))
{
  rw = lig.conf.counted.wide[i,]
  ligand = as.character(rw$Ligand)
  rw$Ligand = NULL
  tot = sum(rw)
  if (tot == 0)
  {
    next ;
  }
  cols = colnames(rw)
  for (ii in 1:ncol(rw))
  {
   freq = rw[,ii] / tot 
   n = cols[ii]
    if ( freq > 0.5 )
   {
      
     row = as.data.frame(t(as.data.frame(c(ligand, n, freq))))
     ligand.dominant.class  = as.data.frame(rbind(ligand.dominant.class , row))
    }
   else 
   {
     #row = as.data.frame(t(as.data.frame(c(ligand, "nonspecific" , freq))))
     #ligand.dominant.class  = as.data.frame(rbind(ligand.dominant.class , row))
   }
  }
}

head(full.data.lig)
chembl = read.table("../../xtal.chembl.combined.fp.txt.uniq.cases", sep = ",", header = F)
chemid = unique(as.character(chembl$V1))
tani = rbind(tani, chembl)

tani$weight = (tani$V3 *1.2) + (tani$V4*1.2) + (tani$V5 *.6) 
tani$z = (tani$weight - mean(tani$weight))/ sd(tani$weight)
tani$p = pnorm(-abs(tani$z))
tani$padj = p.adjust(tani$p, method = "fdr", n = length(tani$p))
tani.sig = tani[which(tani$padj < .1 ), ]
tani.sig = tani.sig[which(tani.sig$weight > mean(tani$weight)) ,]
library(qgraph)
library(tidyr)
tani.sig.full = tani.sig
tani.sig[,c(3:6,8:ncol(tani.sig))] = NULL
tani.sig$V1  = as.character(tani.sig$V1)
tani.sig$V2 = as.character(tani.sig$V2)
##spread to fill in missing pairs 
#tani.sig.wide = spread(data = tani.sig, key = V2 , value = V3 , fill = 0  )
#row.names(tani.sig.wide) = tani.sig.wide$V1
#tani.sig.wide$V1 = NULL
#cols = colnames(tani.sig.wide)
#tani.sig.wide$ligand1 = row.names(tani.sig.wide)
#tani.sig.complete = gather(data = tani.sig.wide , ligand2 , distance , cols)
colnames(ligand.dominant.class) = c("ligand1", "struct", "freq")
missing.pdb= unique(tani.sig[!(tani.sig$V1 %in% tani.sig$V2),]$V1)
missing.df = as.data.frame(cbind(missing.pdb, missing.pdb, rep(3,length(missing.pdb))))
colnames(missing.df)= colnames(tani.sig)
tani.sig.combined = rbind(tani.sig,missing.df)
head(tani.sig.combined)
for ( i in 1:nrow(tani.sig.combined))
{
  rw = tani.sig.combined[i,]
  lig = rw$V2 
  class = as.character(ligand.dominant.class[which(ligand.dominant.class$ligand1 == lig),2])
  
  if ( length(class) == 0)
  {
    if ( lig %in% chemid) 
    {
      #print(paste("chembl"))
      tani.sig.combined[i,4] = "chembl" 
    }
    else 
    {
      #print(paste("nonspecific"))
      tani.sig.combined[i,4] = "nonspecific" 
    }
    
  }
  else 
  {
    tani.sig.combined[i,4] = class 
    #print(paste(lig, class, sep = " "))
  }
}

write.table(tani.sig.combined, file = "../tani.sig.distance.12.04.2017.no.missing.classes.chembl.combined.csv" ,sep = ",", quote = F, eol =  "\n", row.names = F )
