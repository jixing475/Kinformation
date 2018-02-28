###generates training model for DFG and C-Helix prediction 
###Generates variables: 
####training.dfg.rf to train dfg model  
####training.chelix.rf to train chelix model based on data + dfg classification
####test.chelix.predictions final prediction 
#### predictions in : /2_predicted_classes/8.29.17.predicted.chelix.dfg.conformation.csv"
#### and 2_predicted_classes/8.29.17.predicted+verified.chelix.dfg.conformation.csv
library(tidyverse)
library(randomForest)
library(plotly)
library(clusterSim)
##set seed for randomforest 
set.seed(42)
dat = read.table("../../1_manual_classes/stdy_kinase.param.170829.rrahman.csv",header =T, sep ="," )
##get relevant columns 
dat.rel.col = dat[,c(1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,30,33)]
##partition and format test/training data 
training = dat.rel.col[c(1,3:261),] ## row 1 is duplcate of row 2 
row.names(training) = training$pdb_id
training$dfg_st_2 = factor(as.vector(training$dfg_st_2))
training.class = as.factor(as.character(training$Group))
t.c.df = as.data.frame(training.class)
t.c.sp <- separate(t.c.df, training.class, c(paste0("V",LETTERS[1:2])),sep = c(2,4))
training.dfg = training$dfg_st
training[,c(1,2,3,12)] = NULL
test = dat.rel.col[262:nrow(dat.rel.col),]
row.names(test) = test$pdb_id
test[,c(1,2,3,12)] = NULL
##remove row with NA lose about 139 cases 
test.complete = test[complete.cases(test),]
##impute missing training data using c-helix class 

training.impute = rfImpute(x = training, y = training.class ,   ntree = 1000)
training.impute$training.class = NULL
##normalize data to 0
combined.data = as.data.frame(rbind(training.impute, test.complete))
combined.n = data.Normalization(combined.data,type = "n2", normalization = "column" )
training.n = combined.n[1:260,]
test.n = combined.n[261:nrow(combined.n),]
##fix dfg stat values
training.dfg.2 = c()
for (i in 1:length(training.class))
{
  if (training.class[i] == "other")
  {
    training.dfg.2[i] = 3
  }
  else
  {
    g = substr(training.class[i], 3, 5 )
    if (g == "di") {
      training.dfg.2[i] = 1
    }
    else {
      training.dfg.2[i] = 0              
    }
  }
}
training.class = as.factor(t.c.sp[,1])

training.dfg.2 = as.factor(training.dfg.2)
###train DFG classifier 
set.seed(10)
training.dfg.rf = randomForest( training.dfg ~., data=training.n, 
                               ntree = 1000)
##train Chelix classifer 
training.n.dfg = training.n
training.n.dfg$dfg = as.factor(training.dfg) ##add dfg data to data 
set.seed(10)
training.chelix.rf = randomForest( training.class ~., data=training.n.dfg, 
                                ntree = 2000)
training.chelix.rf
training.dfg.rf
##predict test classes 
test.pred.dfg = as.data.frame(predict(object = training.dfg.rf, newdata = test.n, type = "prob"))
test.pred.dfg.2 = test.n
for (i in 1:nrow(test.pred.dfg)) {
  col = (which.max(test.pred.dfg[i,]))
  dfgstat = as.character(row.names(as.data.frame(col)))
  test.pred.dfg.2[i, 14] = dfgstat
} 
test.dfg =  test.pred.dfg.2$V14
test.pred.dfg.2$dfg  = test.dfg
test.pred.dfg.2$V14 = NULL
test.pred.dfg.2$dfg = as.factor(test.pred.dfg.2$dfg)
test.chelix.pred = as.data.frame(predict( object = training.chelix.rf, newdata =  test.pred.dfg.2, type = "prob"))

test.chelix.predictions = test.pred.dfg.2

##assign class to rows based on class with max probablity 
for ( i in 1:nrow(test.chelix.pred))
{
  col = which.max(test.chelix.pred[i,])
  test.chelix.predictions[i,15] =  as.character(row.names(as.data.frame(col)))
  test.chelix.predictions[i,16] = as.numeric(test.chelix.pred[i,col])
}
#test.chelix.predictions[which(test.chelix.predictions$dfg == "random"),15]  = "random"
write.table(x = test.chelix.predictions, file = "../../2_predicted_classes/9.21.17.predicted.chelix.dfg.conformation.csv",  sep =",", quote = F, eol = "\n", row.names = T, col.names = T )


head(test.chelix.predictions)
head(training.n)
training.n$dfg = training.dfg
training.n$V15 = training.class
#training.n[which(training.n$dfg == "random"),17] = "random"  
training.n$V16 = rep(1, nrow(training))
full.data = as.data.frame(rbind(training.n,test.chelix.predictions))
source = c(rep("training",nrow(training)),rep("test",nrow(test.chelix.predictions)))
full.data$source = source
full.data[,18]  = paste(full.data$V15, full.data$dfg, sep = "-") 
full.data[full.data$dfg =="random", 18] = "omega"

write.table(x = full.data, file = "../../2_predicted_classes/9.21.17.predicted+verified.chelix.dfg.conformation.twostage.classification.csv" ,      sep =",", quote = F, eol = "\n", row.names = T, col.names = T )


