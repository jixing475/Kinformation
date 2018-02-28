##generates box/violin plot of features 
##supplemental figure 1 
library(tidyverse)
library(readxl)
library(clusterSim)
library(randomForest)
library(ggplot2)
library(ggthemes)
library(plotly)
##set seed for randomforest 
set.seed(10)
##read data 
dat = as.data.frame(read_excel("../../1_manual_classes/stdy_kinase.param.170829.xlsx",sheet = 1,col_names = T ,))
#get training data and relevant data columns 
dat.2 = data.matrix(dat[2:265,c(6,7,8,9,10,11,12,13,15,16,23,24,33,34,35)])
##impute missing data 
dat.2.impute = as.data.frame(rfImpute(
  x = dat.2,
  y = as.factor(dat$Group[2:265]) ,
  ntree = 1000))
dat.2.impute$y = NULL
##normalize data to fit -1 to 1  
dat.2.scaled = data.Normalization(dat.2.impute,type = "n5", normalization = "column" )
##add classes to normalized data 
dat.2.scaled$clust = dat$Group[2:265]
##add pdb id to normalized data 
row.names(dat.2.scaled) = dat$pdb_id[2:265]
##format data for ggplot 
dat.2.scaled.split = split(dat.2.scaled, f = dat.2.scaled$clust)
complete.df = data.frame()
for (i in dat.2.scaled.split) 
{
    c = unique(as.character(i$clust))
    i$clust = NULL
    i.st = stack(i)
    i.st$clust = rep(c, nrow(i.st))
    complete.df.l = as.data.frame(rbind(complete.df, i.st))
    complete.df = complete.df.l
}
##fill colors for classes 
colors = c("#CAF270","#73D487","#30B096","#40607A", "#453B52")
##plot box/violin plot 

colors = c("#443333","#cc5522","#00ccff","#AA96DA", "grey50")
ggplot(complete.df, aes(clust , values, fill = clust))  + geom_boxplot() + geom_violin(color = "lightgrey",scale = "width", trim =T, alpha = .8)  +  scale_fill_manual(values =  colors ) +  facet_wrap(~ind)  + theme_light(base_size = 20) + ylim(-1.1,1.1)  + theme(legend.position = "none", strip.text = element_text(colour = 'black', size = 15) , axis.title.x = element_blank() , strip.background = element_rect(fill = "grey",colour = "lightgrey"), axis.text = element_text(size = 15)) + ylab(label = "Normalized Score")
plot(p)
 ##interactive figure 
ggplotly(p)
       