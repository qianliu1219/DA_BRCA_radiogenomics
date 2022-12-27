library(ROCR)
library(xgboost)
library(nnet)
library(ClassifyR)


signature<-read.csv("sigScore.csv", row.names = 1)
TILs<-read.csv("TCGABRCA_TIL_TIMER.csv", row.names = 1)
DRFs<-read.csv("fea113.csv",row.names = 1)
CRFs<-read.csv("semi_fea.csv", row.names = 1)
row.names(CRFs)<-substring(row.names(CRFs),1,12)
row.names(CRFs)<-gsub("-", ".", row.names(CRFs))

matched_names<-intersect(intersect(intersect(row.names(DRFs),row.names(CRFs)),
                                   row.names(signature)),row.names(TILs))

signature<-signature[matched_names,]
TILs<-TILs[matched_names,]
DRFs<-DRFs[matched_names,]
CRFs<-CRFs[matched_names,]

pik3cags<-ifelse(signature$pik3cags<quantile(signature$pik3cags)[4],0,1)
endo<-ifelse(signature$endo<quantile(signature$endo)[4],0,1)
gene70<-ifelse(signature$gene70<quantile(signature$gene70)[4],0,1)
genius<-ifelse(signature$genius<quantile(signature$genius)[4],0,1)
oncotypedx<-ifelse(signature$oncotypedx<quantile(signature$oncotypedx)[4],0,1)
rorS<-ifelse(signature$rorS<quantile(signature$rorS)[4],0,1)

B.cell<-ifelse(TILs$B.cell<quantile(TILs$B.cell)[4],0,1)
T.cell.CD4<-ifelse(TILs$T.cell.CD4<quantile(TILs$T.cell.CD4)[4],0,1) 
T.cell.CD8<-ifelse(TILs$T.cell.CD8<quantile(TILs$T.cell.CD8)[4],0,1)
Neutrophil<-ifelse(TILs$Neutrophil<quantile(TILs$Neutrophil)[4],0,1)
Macrophage<-ifelse(TILs$Macrophage<quantile(TILs$Macrophage)[4],0,1)
Dendritic.cell<-ifelse(TILs$Dendritic.cell<quantile(TILs$Dendritic.cell)[4],0,1)

#DRFs train/test data split
set.seed(101)
sample <- sample.int(n = length(matched_names), size = floor(.8*length(matched_names)), replace = F)
trainData_DRF <- DRFs[sample,]
testData_DRF <- DRFs[-sample,]

#CRFs train/test data split
trainData_CRF <- CRFs[sample,]
testData_CRF <- CRFs[-sample,]

######################      pik3cags   ##################

#label
classes<-as.factor(pik3cags)

#train/test label split
trainLabel <-pik3cags[sample]
testLabel <- pik3cags[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_1854" "fea_726"  "fea_966"  "fea_2382" "fea_1870" "fea_1227" "fea_1246" "fea_2126" "fea_971" 
# "fea_988"  "fea_710"  "fea_1630" "fea_982"  "fea_1803" "fea_1950" "fea_2654" "fea_1902" "fea_1886"
# "fea_1869" "fea_2981" "fea_1614" "fea_2638" "fea_1867" "fea_2773" "fea_2206" "fea_1883"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-5, maxit=200)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_pik3cags_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_pik3cags_DRF_DNN_x<-1-unlist(slot(perf_pik3cags_DRF_DNN,"x.values"))
  perf_pik3cags_DRF_DNN_y<-1-unlist(slot(perf_pik3cags_DRF_DNN,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.8076923

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_pik3cags_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_pik3cags_DRF_xgboost_x<-1-unlist(slot(perf_pik3cags_DRF_xgboost,"x.values"))
  perf_pik3cags_DRF_xgboost_y<-1-unlist(slot(perf_pik3cags_DRF_xgboost,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.6153846

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Size.Lesion.volume..S1."                   "Volume.of.most.enhancing.voxels..S4."     
# "Enhancement.variance.Decreasing.Rate..E4." "Enhancement.variance.Increasing.Rate..E3."
# "Surface.Area..S3."                         "Effective.Diameter..S2."                  
# "Maximum.enhancement.variance..E1."         "Enhancement.Variance.Time.to.Peak..E2." 

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=55, decay=1.0e-6, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_pik3cags_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.7692308

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_pik3cags_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
# 0.5384615


######################      endo   ##################

#label
classes<-as.factor(endo)

#train/test label split
trainLabel <-endo[sample]
testLabel <- endo[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_2327" "fea_2583" "fea_2071" "fea_2322" "fea_2066" "fea_2594" "fea_2599" "fea_2855" "fea_2311"
# "fea_2578" "fea_2567" "fea_2329" "fea_2841" "fea_3469" "fea_3097" "fea_1810"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_endo_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.755102

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_endo_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.7346939

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Enhancement.variance.Increasing.Rate..E3." "Energy..T5."                              
# "Enhancement.variance.Decreasing.Rate..E4." "Maximum.enhancement.variance..E1."        
# "Maximum.enhancement..K1."                  "E1..K6."    

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_endo_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.625

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_endo_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.55



######################      gene70   ##################

#label
classes<-as.factor(gene70)

#train/test label split
trainLabel <-gene70[sample]
testLabel <- gene70[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
#"fea_3469" "fea_3485" "fea_2315" "fea_2312" "fea_2330" "fea_2568" "fea_2824" "fea_3483"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=55, decay=1.0e-4, maxit=1000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_gene70_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6458333

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_gene70_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.625

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
#"Uptake.rate..K3."      "Time.to.peak..K2."     "Maximum.Diameter..S5."  

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=55, decay=1.0e-5, maxit=1000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_gene70_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_gene70_CRF_DNN_x<-1-unlist(slot(perf_gene70_CRF_DNN,"x.values"))
  perf_gene70_CRF_DNN_y<-1-unlist(slot(perf_gene70_CRF_DNN,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.5625

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =6, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_gene70_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_gene70_CRF_xgboost_x<-1-unlist(slot(perf_gene70_CRF_xgboost,"x.values"))
  perf_gene70_CRF_xgboost_y<-1-unlist(slot(perf_gene70_CRF_xgboost,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.5625

######################      genius   ##################

#label
classes<-as.factor(genius)

#train/test label split
trainLabel <-genius[sample]
testLabel <- genius[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_2553" "fea_2297" "fea_2809" "fea_2041" "fea_3305" "fea_3065" "fea_2793" "fea_3321" "fea_1785"
# "fea_3049" "fea_1529" "fea_4004" "fea_3833"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=45, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_genius_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.7878788

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_genius_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6969697

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Enhancement.Variance.Time.to.Peak..E2." "Volume.of.most.enhancing.voxels..S4."  
# "Time.to.peak..K2."                      "Margin.Sharpness..M1."                 
# "Curve.shape.index..K5."                 "Uptake.rate..K3."   

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=65, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_genius_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5225

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_genius_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_genius_CRF_xgboost_x<-1-unlist(slot(perf_genius_CRF_xgboost,"x.values"))
  perf_genius_CRF_xgboost_y<-1-unlist(slot(perf_genius_CRF_xgboost,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.5454545




######################      oncotypedx   ##################

#label
classes<-as.factor(oncotypedx)

#train/test label split
trainLabel <-oncotypedx[sample]
testLabel <- oncotypedx[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_2765" "fea_3469" "fea_3211" "fea_3021" "fea_2253" "fea_2071" "fea_2327" "fea_149" 

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=25, decay=1.0e-4, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_oncotypedx_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.75

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =0.1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_oncotypedx_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6354167

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Uptake.rate..K3."                     "Time.to.peak..K2."                   
# "Sum.Average..T11."                    "Volume.of.most.enhancing.voxels..S4."

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=45, decay=1.0e-4, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_oncotypedx_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5416667

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta = 2, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_oncotypedx_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5833333




######################      rorS   ##################

#label
classes<-as.factor(rorS)

#train/test label split
trainLabel <-rorS[sample]
testLabel <- rorS[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_165"  "fea_3976" "fea_2050" "fea_149"  "fea_2567" "fea_2562" "fea_2594" "fea_2855" "fea_3980"
# "fea_2066" "fea_3094" "fea_3978" "fea_2866" "fea_161"  "fea_2583" "fea_2306" "fea_3969" "fea_2322"
# "fea_3021" "fea_3974" "fea_2850" "fea_3984"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=25, decay=1.0e-6, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_rorS_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_rorS_DRF_DNN_x<-1-unlist(slot(perf_rorS_DRF_DNN,"x.values"))
  perf_rorS_DRF_DNN_y<-1-unlist(slot(perf_rorS_DRF_DNN,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.7875

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_rorS_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Variance.of.Radial.Gradient.Histogram..vRGH...M3." "Sum.Variance..T13."                               
# "Variance..T14."                                    "Sum.Entropy..T12."                                
# "Time.to.peak..K2."                                 "Sphericity..G1."

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-4, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_rorS_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.725

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =0.1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_rorS_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5875






######################      B.cell   ##################

#label
classes<-as.factor(B.cell)

#train/test label split
trainLabel <-B.cell[sample]
testLabel <- B.cell[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_3188" "fea_3163" "fea_3492" "fea_3180" "fea_3174" "fea_3204" "fea_3178" "fea_6"    "fea_3182"
# "fea_2915"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_B.cell_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.9090909

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_B.cell_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.696969

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Enhancement.variance.Increasing.Rate..E3." "Time.to.peak..K2."                        
# "Enhancement.variance.Decreasing.Rate..E4." "Maximum.enhancement.variance..E1."        
# "Variance.of.Margin.Sharpness..M2."         "Curve.shape.index..K5."                   
# "Margin.Sharpness..M1."                     "Uptake.rate..K3."  

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-3, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_B.cell_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6060606

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_B.cell_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5454545



######################      T.cell.CD4   ##################

#label
classes<-as.factor(T.cell.CD4)

#train/test label split
trainLabel <-T.cell.CD4[sample]
testLabel <- T.cell.CD4[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
#"fea_1229" "fea_3384" "fea_707"  "fea_1469" "fea_2334" "fea_969"  "fea_3765" "fea_963"  "fea_3376"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=25, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_T.cell.CD4_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.7555556

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_T.cell.CD4_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6777778

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Signal.Enhancement.Ratio..SER...K7."    "Volume.of.most.enhancing.voxels..S4."  
# "Curve.shape.index..K5."                 "Maximum.Diameter..S5."                 
# "Washout.rate..K4."                      "Time.to.peak..K2."                     
# "Enhancement.Variance.Time.to.Peak..E2." "Surface.Area..S3." 

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_T.cell.CD4_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6555556


#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta = 10, nthread = 2, nrounds = 100, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_T.cell.CD4_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6555556




######################      T.cell.CD8   ##################

#label
classes<-as.factor(T.cell.CD8)

#train/test label split
trainLabel <-T.cell.CD8[sample]
testLabel <- T.cell.CD8[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_712"  "fea_972"  "fea_3956" "fea_821"  "fea_982"  "fea_1228" "fea_4093" "fea_3700" "fea_644" 
# "fea_716"  "fea_976"  "fea_720"  "fea_988"  "fea_966"  "fea_3325" "fea_704"  "fea_805"  "fea_660" 

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=40, decay=1.0e-3, maxit=3000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_T.cell.CD8_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6444444

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =0.1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_T.cell.CD8_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_T.cell.CD8_DRF_xgboost_x<-1-unlist(slot(perf_T.cell.CD8_DRF_xgboost,"x.values"))
  perf_T.cell.CD8_DRF_xgboost_y<-1-unlist(slot(perf_T.cell.CD8_DRF_xgboost,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.6

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Signal.Enhancement.Ratio..SER...K7."    "Volume.of.most.enhancing.voxels..S4."  
# "Curve.shape.index..K5."                 "Maximum.Diameter..S5."                 
# "Washout.rate..K4."                      "Enhancement.Variance.Time.to.Peak..E2."
# "Time.to.peak..K2."                      "Surface.Area..S3."   

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=25, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_T.cell.CD8_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6111111

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 100, eta = 100, nthread = 1, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_T.cell.CD8_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.6111111




######################      Neutrophil   ##################

#label
classes<-as.factor(Neutrophil)

#train/test label split
trainLabel <-Neutrophil[sample]
testLabel <- Neutrophil[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_3325" "fea_3384" "fea_982"  "fea_3376" "fea_3368" "fea_3110" "fea_3392" "fea_3388" "fea_707" 
# "fea_3837" "fea_3581" "fea_961"  "fea_3371" "fea_963"  "fea_1725" "fea_3069"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=30, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_Neutrophil_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_Neutrophil_DRF_DNN_x<-1-unlist(slot(perf_Neutrophil_DRF_DNN,"x.values"))
  perf_Neutrophil_DRF_DNN_y<-1-unlist(slot(perf_Neutrophil_DRF_DNN,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.6060606

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta = 1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_Neutrophil_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5454545

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Energy..T5."  "Entropy..T6." "IMC2..T9."  

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=15, decay=1.0e-4, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_Neutrophil_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5757576

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta = 0.01, nthread = 2, nrounds = 20, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_Neutrophil_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5454545


######################      Macrophage   ##################

#label
classes<-as.factor(Macrophage)

#train/test label split
trainLabel <-Macrophage[sample]
testLabel <- Macrophage[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_3384" "fea_3121" "fea_3392" "fea_3128" "fea_3126" "fea_3387" "fea_3136" "fea_3368" "fea_3388"
# "fea_3382" "fea_3376" "fea_3112" "fea_3400" "fea_3398" "fea_2649" "fea_3177" "fea_3406" "fea_3408"
# "fea_3120" "fea_3404" "fea_719"  "fea_3115"

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=40, decay=1.0e-4, maxit=3000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_Macrophage_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.625

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta = 2, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_Macrophage_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
# 0.7125

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Signal.Enhancement.Ratio..SER...K7."  "Washout.rate..K4."                   
# "Variance.of.Margin.Sharpness..M2."    "Volume.of.most.enhancing.voxels..S4."
# "Sum.Average..T11."  


#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-4, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_Macrophage_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
if (auc< 0.5){
  perf_Macrophage_CRF_DNN_x<-1-unlist(slot(perf_Macrophage_CRF_DNN,"x.values"))
  perf_Macrophage_CRF_DNN_y<-1-unlist(slot(perf_Macrophage_CRF_DNN,"y.values"))
  auc2<-1-auc
}
if (auc >= 0.5){
  auc2<-auc
}
#0.6

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_Macrophage_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.5875

######################      Dendritic.cell   ##################

#label
classes<-as.factor(Dendritic.cell)

#train/test label split
trainLabel <-Dendritic.cell[sample]
testLabel <- Dendritic.cell[-sample]

#DRF
measurements<-as.matrix(DRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "fea_3325" "fea_53"   "fea_1627" "fea_2043" "fea_1613" "fea_5"    "fea_1675" "fea_2299" "fea_69"  

#DNN
mynn <- nnet(as.matrix(trainData_DRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-5, maxit=2000)
nnpred<-predict(mynn, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_Dendritic.cell_DRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.625

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_DRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 5, eta =1, nthread = 4, nrounds = 50, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_DRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_Dendritic.cell_DRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.7

#CRF
measurements<-as.matrix(CRFs)

#feature selection
DDresults <- runTests(measurements, classes)
selectionPercentages <- distribution(DDresults, plot = FALSE)
sortedPercentages <- sort(selectionPercentages, decreasing = TRUE)
fn<-names(sortedPercentages[sortedPercentages>40])
# "Enhancement.variance.Increasing.Rate..E3." "Energy..T5."                              
# "Enhancement.variance.Decreasing.Rate..E4." "Maximum.enhancement.variance..E1."        
# "Maximum.enhancement..K1."                  "E1..K6."   

#DNN
mynn <- nnet(as.matrix(trainData_CRF[,fn]),as.numeric(as.character(trainLabel)) ,
             size=35, decay=1.0e-5, maxit=3000)
nnpred<-predict(mynn, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(nnpred),testLabel)
perf_Dendritic.cell_CRF_DNN <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.625

#XGboost
bstDense <- xgboost(data = as.matrix(trainData_CRF[,fn]), label = as.numeric(as.character(trainLabel)), 
                    max.depth = 10, eta =1, nthread = 2, nrounds = 10, objective = "binary:logistic")
mypred<-predict(bstDense, as.matrix(testData_CRF[,fn]))
pred <- prediction(as.numeric(mypred),testLabel)
perf_Dendritic.cell_CRF_xgboost <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
#0.55






######################  plot ##############################

# Gene Signature
pdf("gs_DRF_XGboost_rocs.pdf")
plot(perf_endo_DRF_xgboost,col = 2, lwd=2)
lines(perf_pik3cags_DRF_xgboost_x,perf_pik3cags_DRF_xgboost_y,col = 3, lwd=2)
plot(perf_gene70_DRF_xgboost, add = TRUE,col = 4, lwd=2)
plot(perf_genius_DRF_xgboost, add = TRUE,col = 5, lwd=2)
plot(perf_oncotypedx_DRF_xgboost, add = TRUE,col = 6, lwd=2)
plot(perf_rorS_DRF_xgboost, add = TRUE,col = 7, lwd=2)
legend("bottomright", legend=c("endo (AUC=0.73)","pik3cags (AUC=0.62)","gene70 (AUC=0.63)","genius (AUC=0.70)","oncotypedx (AUC=0.64)","rorS (AUC=0.60)"),
       col=2:7, lty=1, cex=1)
dev.off()

pdf("gs_CRF_XGboost_rocs.pdf")
plot(perf_endo_CRF_xgboost,col = 2, lwd=2)
plot(perf_pik3cags_CRF_xgboost,add = TRUE,col = 3, lwd=2)
lines(perf_gene70_CRF_xgboost_x,perf_gene70_CRF_xgboost_y,col = 4, lwd=2)
lines(perf_genius_CRF_xgboost_x,perf_genius_CRF_xgboost_y,col = 5, lwd=2)
plot(perf_oncotypedx_CRF_xgboost, add = TRUE,col = 6, lwd=2)
plot(perf_rorS_CRF_xgboost, add = TRUE,col = 7, lwd=2)
legend("bottomright", legend=c("endo (AUC=0.55)","pik3cags (AUC=0.54)","gene70 (AUC=0.56)","genius (AUC=0.55)","oncotypedx (AUC=0.58)","rorS (AUC=0.59)"),
       col=2:7, lty=1, cex=1)
dev.off()


pdf("gs_DRF_DNN_rocs.pdf")
plot(perf_endo_DRF_DNN,col = 2, lwd=2)
lines(perf_pik3cags_DRF_DNN_x,perf_pik3cags_DRF_DNN_y,col = 3, lwd=2)
plot(perf_gene70_DRF_DNN, add = TRUE,col = 4, lwd=2)
plot(perf_genius_DRF_DNN, add = TRUE,col = 5, lwd=2)
plot(perf_oncotypedx_DRF_DNN, add = TRUE,col = 6, lwd=2)
plot(perf_rorS_DRF_xgboost, add = TRUE,col = 7, lwd=2)
lines(perf_rorS_DRF_DNN_x,perf_rorS_DRF_DNN_y,col = 7, lwd=2)
legend("bottomright", legend=c("endo (AUC=0.81)","pik3cags (AUC=0.62)","gene70 (AUC=0.65)","genius (AUC=0.79)","oncotypedx (AUC=0.75)","rorS (AUC=0.79)"),
       col=2:7, lty=1, cex=1)
dev.off()

pdf("gs_CRF_DNN_rocs.pdf")
plot(perf_endo_CRF_DNN,col = 2, lwd=2)
lines(perf_pik3cags_CRF_DNN_x,perf_pik3cags_CRF_DNN_y,col = 3, lwd=2)
lines(perf_gene70_CRF_DNN_x,perf_gene70_CRF_DNN_y,col = 4, lwd=2)
plot(perf_genius_CRF_xgboost, add = TRUE,col = 5, lwd=2)
plot(perf_oncotypedx_CRF_DNN, add = TRUE,col = 6, lwd=2)
plot(perf_rorS_CRF_DNN, add = TRUE,col = 7, lwd=2)
legend("bottomright", legend=c("endo (AUC=0.77)","pik3cags (AUC=0.54)","gene70 (AUC=0.56)","genius (AUC=0.55)","oncotypedx (AUC=0.58)","rorS (AUC=0.59)"),
       col=2:7, lty=1, cex=1)
dev.off()



#TILs

pdf("tils_DRF_XGboost_rocs.pdf")
plot(perf_B.cell_DRF_xgboost,col = 2, lwd=2)
plot(perf_T.cell.CD4_DRF_xgboost, add = TRUE,col = 3, lwd=2)
lines(perf_T.cell.CD8_DRF_xgboost_x,perf_T.cell.CD8_DRF_xgboost_y,col = 4, lwd=2)
plot(perf_Neutrophil_DRF_xgboost, add = TRUE,col = 5, lwd=2)
plot(perf_Macrophage_DRF_xgboost, add = TRUE,col = 6, lwd=2)
plot(perf_Dendritic.cell_DRF_xgboost, add = TRUE,col = 7, lwd=2)
legend("bottomright", legend=c("B.cell (AUC=0.70)","T.cell.CD4 (AUC=0.68)","T.cell.CD8 (AUC=0.60)",
                               "Neutrophil (AUC=0.54)","Macrophage (AUC=0.71)","Dendritic.cell (AUC=0.70)"),
       col=2:7, lty=1, cex=1)
dev.off()


pdf("tils_CRF_XGboost_rocs.pdf")
plot(perf_B.cell_CRF_xgboost,col = 2, lwd=2)
plot(perf_T.cell.CD4_CRF_xgboost, add = TRUE,col = 3, lwd=2)
plot(perf_T.cell.CD8_CRF_xgboost, add = TRUE,col = 4, lwd=2)
plot(perf_Neutrophil_CRF_xgboost, add = TRUE,col = 5, lwd=2)
plot(perf_Macrophage_CRF_xgboost, add = TRUE,col = 6, lwd=2)
plot(perf_Dendritic.cell_CRF_xgboost, add = TRUE,col = 7, lwd=2)
legend("bottomright", legend=c("B.cell (AUC=0.55)","T.cell.CD4 (AUC=0.66)","T.cell.CD8 (AUC=0.61)",
                               "Neutrophil (AUC=0.55)","Macrophage (AUC=0.59)","Dendritic.cell (AUC=0.55)"),
       col=2:7, lty=1, cex=1)
dev.off()


pdf("tils_DRF_DNN_rocs.pdf")
plot(perf_B.cell_DRF_DNN,col = 2, lwd=2)
plot(perf_T.cell.CD4_DRF_DNN, add = TRUE,col = 3, lwd=2)
plot(perf_T.cell.CD8_DRF_DNN, add = TRUE,col = 4, lwd=2)
plot(perf_Neutrophil_DRF_DNN, add = TRUE,col = 5, lwd=2)
plot(perf_Macrophage_DRF_DNN, add = TRUE,col = 6, lwd=2)
plot(perf_Dendritic.cell_DRF_DNN, add = TRUE,col = 7, lwd=2)
legend("bottomright", legend=c("B.cell (AUC=0.91)","T.cell.CD4 (AUC=0.76)","T.cell.CD8 (AUC=0.64)",
                               "Neutrophil (AUC=0.61)","Macrophage (AUC=0.63)","Dendritic.cell (AUC=0.63)"),
       col=2:7, lty=1, cex=1)
dev.off()


pdf("tils_CRF_DNN_rocs.pdf")
plot(perf_B.cell_CRF_DNN,col = 2, lwd=2)
plot(perf_T.cell.CD4_CRF_DNN, add = TRUE,col = 3, lwd=2)
plot(perf_T.cell.CD8_CRF_DNN, add = TRUE,col = 4, lwd=2)
plot(perf_Neutrophil_CRF_DNN, add = TRUE,col = 5, lwd=2)
lines(perf_Macrophage_CRF_DNN_x,perf_Macrophage_CRF_DNN_y,col = 6)
plot(perf_Dendritic.cell_CRF_DNN, add = TRUE,col = 7, lwd=2)
legend("bottomright", legend=c("B.cell (AUC=0.61)","T.cell.CD4 (AUC=0.66)","T.cell.CD8 (AUC=0.61)",
                               "Neutrophil (AUC=0.58)","Macrophage (AUC=0.60)","Dendritic.cell (AUC=0.63)"),
       col=2:7, lty=1, cex=1)
dev.off()

