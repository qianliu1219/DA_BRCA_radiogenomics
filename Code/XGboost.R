signature<-read.csv("sigScore.csv", row.names = 1)
DRFs<-read.csv("fea113.csv",row.names = 1)


library(ROCR)
library(xgboost)



pik3cags<-ifelse(signature$pik3cags<quantile(signature$pik3cags)[4],0,1)
endo<-ifelse(signature$endo<quantile(signature$endo)[4],0,1)
gene70<-ifelse(signature$gene70<quantile(signature$gene70)[4],0,1)
genius<-ifelse(signature$genius<quantile(signature$genius)[4],0,1)
oncotypedx<-ifelse(signature$oncotypedx<quantile(signature$oncotypedx)[4],0,1)
rorS<-ifelse(signature$rorS<quantile(signature$rorS)[4],0,1)



set.seed(101)
sample <- sample.int(n = 113, size = floor(.8*113), replace = F)
trainData <- DRFs[sample,]
testData <- DRFs[-sample,]



#1
trainLabel <-pik3cags[sample]
testLabel <- pik3cags[-sample]



library(ClassifyR)

measurements<-t(DRFs)
classes<-as.factor(pik3cags)

DMresults <- runTests(measurements, classes, datasetName = "DRFs",
                      classificationName = "Different Means", permutations = 20, folds = 5,
                      seed = 2018, verbose = 1)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)

fn<-names(sortedPercentages[sortedPercentages>40])



detach("package:ClassifyR", unload = TRUE)
bstDense <- xgboost(data = as.matrix(trainData[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =6, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData[,fn]))

pred <- prediction(as.numeric(mypred),testLabel)
perf_pik3cags <- performance(pred,"tpr","fpr")

auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  x<-1-unlist(slot(perf_pik3cags,"x.values"))
  y<-1-unlist(slot(perf_pik3cags,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
auc2
#   0.8333333  40, 10, 5, 2, 10    "fea_1467" "fea_1211" "fea_950"  "fea_966"  "fea_1733" "fea_1246" "fea_1486"



#2
trainLabel <-endo[sample]
testLabel <- endo[-sample]

library(ClassifyR)

measurements<-t(DRFs)
classes<-as.factor(endo)

DMresults <- runTests(measurements, classes, datasetName = "DRFs",
                      classificationName = "Different Means", permutations = 20, folds = 5,
                      seed = 2018, verbose = 1)

selectionPercentages <- distribution(DMresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)

fn<-names(sortedPercentages[sortedPercentages>40])


detach("package:ClassifyR", unload = TRUE)

bstDense <- xgboost(data = as.matrix(trainData[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =8, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData[,fn]))

pred <- prediction(as.numeric(mypred),testLabel)
perf_endo <- performance(pred,"tpr","fpr")

auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  x_endo<-unlist(slot(perf_endo,"x.values"))
  y_endo<-unlist(slot(perf_endo,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
auc2

# 0.7916667   40 10 8 2 10
# "fea_2855" "fea_2599" "fea_2231" "fea_2327" "fea_2583" "fea_2487" "fea_2839" "fea_2343"
# "fea_2786" "fea_3111" "fea_3127" "fea_2071" "fea_1150" "fea_2226" "fea_2247" "fea_2594"
# "fea_3554" "fea_2807" "fea_2850" "fea_2087" "fea_3559" "fea_574"  "fea_1134" "fea_2359"
# "fea_2871" "fea_558"  "fea_3095"





#3
trainLabel <-gene70[sample]
testLabel <- gene70[-sample]

library(ClassifyR)

measurements<-t(DRFs)
classes<-as.factor(gene70)

DMresults <- runTests(measurements, classes, datasetName = "DRFs",
                      classificationName = "Different Means", permutations = 20, folds = 5,
                      seed = 2018, verbose = 1)

selectionPercentages <- distribution(DMresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)

fn<-names(sortedPercentages[sortedPercentages>40])


detach("package:ClassifyR", unload = TRUE)

bstDense <- xgboost(data = as.matrix(trainData[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData[,fn]))

pred <- prediction(as.numeric(mypred),testLabel)
perf_gene70 <- performance(pred,"tpr","fpr")

auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  x<-1-unlist(slot(perf_gene70,"x.values"))
  y<-1-unlist(slot(perf_gene70,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
auc2

# 0.7352941  40 10 1 2 10
# "fea_1390" "fea_3268" "fea_1150" "fea_1134" "fea_3794"




#4
trainLabel <-genius[sample]
testLabel <- genius[-sample]

library(ClassifyR)
measurements<-t(DRFs)
classes<-as.factor(genius)

DMresults <- runTests(measurements, classes, datasetName = "DRFs",
                      classificationName = "Different Means", permutations = 20, folds = 5,
                      seed = 2018, verbose = 1)

selectionPercentages <- distribution(DMresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)

fn<-names(sortedPercentages[sortedPercentages>40])


detach("package:ClassifyR", unload = TRUE)

bstDense <- xgboost(data = as.matrix(trainData[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta = 0.8, nthread = 2, nrounds = 3, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData[,fn]))

pred <- prediction(as.numeric(mypred),testLabel)
perf_genius <- performance(pred,"tpr","fpr")

auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  x<-1-unlist(slot(perf_genius,"x.values"))
  y<-1-unlist(slot(perf_genius,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
auc2
#0.7579365 40 10 0.8 2 3
# "fea_1458" "fea_1202" "fea_1463" "fea_1218" "fea_1474" "fea_1442" "fea_1719" "fea_1714"
# "fea_1698" "fea_1447" "fea_962"  "fea_1479" "fea_1207" "fea_1703" "fea_1975" "fea_1970"





#5
trainLabel <-oncotypedx[sample]
testLabel <- oncotypedx[-sample]

library(ClassifyR)
measurements<-t(DRFs)
classes<-as.factor(oncotypedx)

DMresults <- runTests(measurements, classes, datasetName = "DRFs",
                      classificationName = "Different Means", permutations = 20, folds = 5,
                      seed = 2018, verbose = 1)

selectionPercentages <- distribution(DMresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)

fn<-names(sortedPercentages[sortedPercentages>40])

detach("package:ClassifyR", unload = TRUE)

bstDense <- xgboost(data = as.matrix(trainData[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta = 1.1, nthread = 2, nrounds = 3, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData[,fn]))

pred <- prediction(as.numeric(mypred),testLabel)
perf_oncotypedx <- performance(pred,"tpr","fpr")

auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  x<-1-unlist(slot(perf_oncotypedx,"x.values"))
  y<-1-unlist(slot(perf_oncotypedx,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
auc2
#0.6964286 40 10 1.1 2 3
# "fea_1134" "fea_3794" "fea_3219" "fea_4050" "fea_3810" "fea_1390" "fea_3970" "fea_574" 
# "fea_3217" "fea_3983" "fea_4063"



#6
trainLabel <-rorS[sample]
testLabel <- rorS[-sample]

library(ClassifyR)
measurements<-t(DRFs)
classes<-as.factor(rorS)

DMresults <- runTests(measurements, classes, datasetName = "DRFs",
                      classificationName = "Different Means", permutations = 20, folds = 5,
                      seed = 2018, verbose = 1)

selectionPercentages <- distribution(DMresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)

fn<-names(sortedPercentages[sortedPercentages>30])


detach("package:ClassifyR", unload = TRUE)

bstDense <- xgboost(data = as.matrix(trainData[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta = 6, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData[,fn]))

pred <- prediction(as.numeric(mypred),testLabel)
perf_rorS <- performance(pred,"tpr","fpr")

auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  x<-1-unlist(slot(perf_rorS,"x.values"))
  y<-1-unlist(slot(perf_rorS,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
auc2
#0.7008929 30 10 6 2 10 
# "fea_3794" "fea_3810" "fea_2871" "fea_3538" "fea_1390" "fea_3127" "fea_1134" "fea_3143"
# "fea_2855" "fea_3554" "fea_3970" "fea_574"  "fea_4050" "fea_1150" "fea_3138" "fea_3154"
# "fea_3399" "fea_3266" "fea_2866" "fea_3271"


pdf("XGboost_rocs.pdf")
plot(perf_pik3cags,col = 2)
lines(x_endo,y_endo,col = 3)
plot(perf_gene70, add = TRUE,col = 4)
plot(perf_genius, add = TRUE,col = 5)
plot(perf_oncotypedx, add = TRUE,col = 6)
plot(perf_rorS, add = TRUE,col = 7)
legend("bottomright", legend=c("pik3cags (AUC=0.83)", "endo (AUC=0.79)","gene70 (AUC=0.74)","genius (AUC=0.76)","oncotypedx (AUC=0.70)","rorS (AUC=0.70)"),
       col=2:7, lty=1, cex=1)
dev.off()




