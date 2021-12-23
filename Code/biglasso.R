rm(list=ls())

library(biglasso)
library(glmnet)
library(ROCR)
library(MLmetrics)


#clinic
clinic<-read.csv("clinic.csv",header = T,row.names=1)
ER_Status<-ifelse(as.vector(clinic$ER.Status)=="Positive",1,0) 
PR_Status<-ifelse(as.vector(clinic$PR.Status)=="Positive",1,0)
HER2_Status<-ifelse(as.vector(clinic$HER2.Final.Status)=="Positive",1,0)
T_Status<-ifelse(as.vector(clinic$Tumor..T1.Coded)=="Positive",1,0)
N_Status<-ifelse(as.vector(clinic$Node.Coded)=="Positive",1,0)

Status<-ER_Status
#Status<-PR_Status
#Status<-HER2_Status
#Status<-T_Status
#Status<-N_Status

set.seed(101) 
sample <- sample.int(n = length(Status[!is.na(Status)]), size = floor(.7*length(Status[!is.na(Status)])), replace = F)

#RFs

RFs<-read.csv("semi_features.csv",header = T,row.names = 1)

trainData <- RFs[!is.na(Status),][sample,-1]
testData <- RFs[!is.na(Status),][-sample,-1]
trainLabel <-Status[!is.na(Status)][sample]
testLabel <- Status[!is.na(Status)][-sample]


##DRFs
#DRFs<-read.csv("features.csv",header = T,row.names = 1)

#trainData <- DRFs[!is.na(Status),][sample,-1]
#testData <- DRFs[!is.na(Status),][-sample,-1]
#trainLabel <-Status[!is.na(Status)][sample]
#testLabel <- Status[!is.na(Status)][-sample]



cvlasso<-cv.glmnet(as.matrix(trainData),trainLabel,type.measure="auc",family="binomial",nfolds=5)


x<-as.big.matrix(as.matrix(trainData))
is.big.matrix(x)
fit <- biglasso(x,trainLabel,family="binomial",nlambda = 20)



x_test<-as.big.matrix(as.matrix(testData))

lassopred<-predict(fit, x_test, type="response", lambda=fit$lambda)
lassopred<-as.matrix(lassopred)


AUC<-rep(NA,nrow(testData))
for (i in 1:nrow(testData)){
  pred <- prediction(lassopred[,i],testLabel)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  AUC[i] <- unlist(slot(auc, "y.values"))



mybeta<-fit$beta 
beta<-as.matrix(mybeta)[-1,]
beta<-as.data.frame(beta)

nzero<-rep(NA, nrow(testData))
for (i in (1:nrow(testData))){
  nzero[i]<-sum(beta[,i]!=0)
}
non0<-cbind(colnames(beta),nzero)

non0<-cbind(non0,AUC)



lassopred<-predict(fit, x_test, type="class", lambda=fit$lambda)
lassopred<-as.matrix(lassopred)
F1<-rep(NA,nrow(testData))
for (i in 1:nrow(testData)){
  if (length(unique(lassopred[,i]))==1){
    F1[i]<-0
  } else {
    F1[i]<-F1_Score(y_pred = lassopred[,i], y_true = testLabel, positive = "1")
  } 
}


ERtable<-cbind(non0,F1)
ERtable<-as.data.frame(ERtable)

#PRtable<-non0
#PRtable<-as.data.frame(PRtable)

#HERtable<-non0
#HERtable<-as.data.frame(HERtable)

#Ttable<-non0
#Ttable<-as.data.frame(Ttable)

#Ntable<-non0
#Ntable<-as.data.frame(Ntable)


dev.off()

plot(Ttable$nzero,Ttable$AUC,type="l",col=1,ylim=0:1)
lines(Ntable$nzero,Ntable$AUC,col= 2)
lines(ERtable$nzero,ERtable$AUC,col= 3)
lines(PRtable$nzero,PRtable$AUC,col= 4)
lines(HERtable$nzero,HERtable$AUC,col= 5)
legend("bottomright",c("Pathological_T","Pathological_N","ER status","PR status","HER2 status"),lty=1,col = 1:5)

plot(Ttable$nzero,Ttable$F1,type="l",col=1,ylim=0:1)
lines(Ntable$nzero,Ntable$F1,col= 2)
lines(ERtable$nzero,ERtable$F1,col= 3)
lines(PRtable$nzero,PRtable$F1,col= 4)
lines(HERtable$nzero,HERtable$F1,col= 5)
legend("bottomright",c("Pathological_T","Pathological_N","ER status","PR status","HER2 status"),lty=1,col = 1:5)










