####Clustering#####

library(gridGraphics)
library(gridExtra)
library(tidyr)
library(reshape)
library(gplots)

grab_grob <- function(){
  grid.echo()
  grid.grab()
}

Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/Data_Visualization_Analysis/PCA_And_Clustering/";
source(paste(Dir,"myheatmap.txt",sep=""))

###Panel1###############################################
a=read.table(paste(Dir,"9Envs_envParas_DAP200",sep=""),header = T)
a$days=rep(1:200,9)
a=a[,c(1,6,12)];

data_wide <- spread(a, env_code, GDD)
a=t(data_wide[,-1])
gdd.std=rescaler(a)

gdd.dist=dist(gdd.std)
gdd.hc=hclust(gdd.dist,method="ward.D")


mycol=c("#4daf4a","#377eb8","#e41a1c")
names(mycol)=c(1,2,3)
mytree=cutree(gdd.hc,3)
x=mycol[match(mytree,names(mycol))]

dev.new(width=6, height=6)

heatmap.2(as.matrix(gdd.dist),Rowv=as.dendrogram(gdd.hc),Colv=as.dendrogram(gdd.hc),trace="none",main="Temperature",
          revC=T,symm=T,ColSideColors=x,cexRow=1,cexCol=1,srtCol=360,adjCol=c(0.5,0),
          key.title=NA,key.xlab=NA,key.ylab=NA,key=T,keysize=1.0,
          key.par=list(mar=c(3,3,3,0),cex=0.6,bty="n",lwd=0.01))

f1=grab_grob()

###Panel2###############################################
b=read.table(paste(Dir,"LbE_table",sep=""),header = T)
b=t(b[,-1])

gdd.std=rescaler(b)

gdd.dist=dist(gdd.std)
gdd.hc=hclust(gdd.dist,method="ward.D")


mycol=c("#377eb8","#4daf4a","#e41a1c")
names(mycol)=c(1,2,3)
mytree=cutree(gdd.hc,3)
x=mycol[match(mytree,names(mycol))]
dev.new(width=6, height=6)

heatmap.2(as.matrix(gdd.dist),Rowv=as.dendrogram(gdd.hc),Colv=as.dendrogram(gdd.hc),trace="none",main="Flowering time",
          revC=T,symm=T,ColSideColors=x,cexRow=1,cexCol=1,srtCol=360,adjCol=c(0.5,0),
          key.title=NA,key.xlab=NA,key.ylab=NA,key=T,keysize=1.0,
          key.par=list(mar=c(3,3,3,0),cex=0.6,bty="n",lwd=0.01))

f2=grab_grob()

png(paste(Dir,"Clustering.png",sep=""),width=12,height=6,units="in",res=600)

lay <- grid.layout(nrow = 1, ncol=2)
pushViewport(viewport(layout = lay))
grid.draw(editGrob(f1, vp=viewport(layout.pos.row = 1, 
                                   layout.pos.col = 1, clip=TRUE)))
grid.draw(editGrob(f2, vp=viewport(layout.pos.row = 1, 
                                   layout.pos.col = 2, clip=TRUE)))
upViewport(1)


dev.off()








