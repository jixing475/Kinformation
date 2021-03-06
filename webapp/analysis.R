x= read.table("../tani.sig.distance.12.04.2017.no.missing.classes.chembl.combined.csv",sep = ",", header = T)
dist = x[,c(2,1,3)]
library(tidyverse)
d =dist[!duplicated(dist),]
types = x[,c(2,4)]
colnames(types) = c("n", "type")
x.wide = spread(data = d, key = V2 ,value = weight, fill = 0 )
row.names(x.wide) = x.wide$V1
x.wide$V1 = NULL
library(Rtsne)
x.tsne = Rtsne(x.wide , dims = 3, perplexity = 50 , max_iter = 1000, theta = 0, pca = T,  is_distance = F,   check_duplicates = F, verbose = T)
x.tsne.df = as.data.frame(x.tsne$Y)
x.tsne.df$n = row.names(x.wide)
x.tsne.df.m = merge(x.tsne.df, types, by = "n")
library(plotly)
plot_ly(data = x.tsne.df.m, x = ~V1, y = ~V2, z = ~V3,color = ~type, text=~paste('Name: ', n) )
