rm(list=ls())
library(nlme)

#lme
signature<-read.csv("sigScore.csv",header = T,row.names=1)
DRFs<-read.csv("features.csv",header=T,row.names=1)

 
res<-matrix(NA,4096,12)
row.names(res)<-colnames(Fea)
colnames(res)<-paste(rep(colnames(signature),each=2),c("Pvalue","adj.Pvalue"),sep="_")

start_time <- Sys.time()
for(l in 1:6){
  phe<-cbind(DRFs[,1],signature[,l])
  Fea<-DRFs[,-1]
  for (i in 1:4096){
    tryCatch({
      MyData<-cbind(phe,Fea[,i])
      MyData<-as.data.frame(MyData)
      colnames(MyData)<-c("PatientID","Signature","Feature")
      MyData_group<-groupedData(Feature ~ Signature|PatientID,data = MyData)
      fm<-lme(Feature ~ Signature,data=MyData_group,random=~1|PatientID
              ,control = lmeControl(maxIter = 1,niterEM = 1)
      )
      res[i,(2*l-1)]<-coef(summary(fm))[2,5]
      print(i)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    res[,2*l]<-p.adjust(res[,(2*l-1)],method = "BH")
  }
}
end_time <- Sys.time()
end_time - start_time


### heatmap #####
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


png("lme_heatmap.png", width=2400, height=1800)
result<-res[,c(2,4,6,8,10,12)] #adjusted p-value
result<-ifelse(result<0.05,0,1)
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="ward.D2")}

heatmap.3(result, na.rm = TRUE, scale="none",margins=c(1,12), hclustfun=myclust, distfun=mydist,
          Rowv=F, Colv=T, dendrogram ="column",symbreaks=F, key=F,symkey=F,
          density.info="none", trace="none", 
          labRow=colnames(Signature),
          col=c("red","blue"))
dev.off() 





