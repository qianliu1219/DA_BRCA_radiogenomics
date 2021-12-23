fea113<-read.csv("fea113.csv",header = T,row.names=1)
clinic<-read.csv("clinic.csv",header = T,row.names=1)
clab<-clinic[,2:6]
clab<-ifelse(clab=="Positive","black","gray")

main_title="Image features"
par(cex.main=1)
col_breaks = c(seq(0,0.35,length=10),  
               seq(0.35,0.6,length=10),           
               seq(0.6,1,length=10))    
heatmap.3(t(fea113), na.rm = TRUE, scale="none", dendrogram="both",margins=c(1,1),
          breaks=col_breaks,Rowv=T, Colv=T,ColSideColors=clab, symbreaks=FALSE, 
          key=TRUE,symkey=F,density.info="none", trace="none", main=main_title, 
          labCol=FALSE, col = colorRampPalette(brewer.pal(9,"Spectral"))(29),
          cexRow=0.5,ColSideColorsSize=7, RowSideColorsSize=2)
legend(x=0.85,y=1.0,cex=0.5,
       legend=c("Positive","Negative","NA"),
       fill=c("black","gray","white"), 
       border=FALSE, bty="n", y.intersp = 0.7)
dev.off()


library(Rtsne)
rtsne_out <- Rtsne(as.matrix(fea113),perplexity = 5,max_iter = 3000)
cluster<-as.data.frame(rtsne_out$Y)
row.names(cluster)<-row.names(clinic)
colnames(cluster)<-c("x","y")

cluster1<-cluster[cluster$x < -40,]
cluster2<-cluster[cluster$x > -40 & cluster$x < 0 & cluster$y < 0,]
cluster3<-cluster[cluster$x > -40 & cluster$x < 0 & cluster$y > 0,]
cluster4<-cluster[cluster$x > 0 & cluster$y < 0,]
cluster5<-cluster[cluster$x > 0 & cluster$y > 0,]

image_count<-read.csv("image_count.csv",header = T,row.names=1)
sum(image_count[row.names(cluster1),2])
sum(image_count[row.names(cluster2),2])
sum(image_count[row.names(cluster3),2])
sum(image_count[row.names(cluster4),2])
sum(image_count[row.names(cluster5),2])

tsne60229<-read.csv("tsne60229.csv",header = T,row.names = 1)



plot(cluster1,xlab = "t-SNE_1",ylab = "t-SNE_2",pch=16,cex=0.5,col=2,xlim=c(-65,50),ylim=c(-50,50))
points(cluster2,xlab = "t-SNE_1",ylab = "t-SNE_2",pch=16,cex=0.5,col=3)
points(cluster3,xlab = "t-SNE_1",ylab = "t-SNE_2",pch=16,cex=0.5,col=4)
points(cluster4,xlab = "t-SNE_1",ylab = "t-SNE_2",pch=16,cex=0.5,col=5)
points(cluster5,xlab = "t-SNE_1",ylab = "t-SNE_2",pch=16,cex=0.5,col=6)

plot(tsne60229[substr(row.names(tsne60229),1,12)==row.names(cluster1),],xlab = "t-SNE_1",ylab = "t-SNE_2",
     pch=16,cex=0.5,col=2,xlim=c(-40,40),ylim=c(-50,50))
points(tsne60229[substr(row.names(tsne60229),1,12)==row.names(cluster2),],xlab = "t-SNE_1",ylab = "t-SNE_2",
       pch=16,cex=0.5,col=3)
points(tsne60229[substr(row.names(tsne60229),1,12)==row.names(cluster3),],xlab = "t-SNE_1",ylab = "t-SNE_2",
       pch=16,cex=0.5,col=4)
points(tsne60229[substr(row.names(tsne60229),1,12)==row.names(cluster4),],xlab = "t-SNE_1",ylab = "t-SNE_2",
       pch=16,cex=0.5,col=5)
points(tsne60229[substr(row.names(tsne60229),1,12)==row.names(cluster5),],xlab = "t-SNE_1",ylab = "t-SNE_2",
       pch=16,cex=0.5,col=6)


dev.off()
