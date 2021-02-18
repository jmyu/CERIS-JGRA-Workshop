########PCA##########################################
Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/Data_Visualization_Analysis/PCA_And_Clustering/";

png(paste(Dir,"PCA.png",sep=""),width=3.5,height=3.5,pointsize=9,
    units="in",res=600)

mars <- c(3, 2.5, 1.5, 1.5);
par(mar = mars , mgp = c(1, 0.1, 0), tck = -0.0075);

library("plotrix")
ge=read.table(paste(Dir,"GbE_table",sep=""),header=T,sep="\t");
ge.1=ge[,-1]

NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
ge.1[] <- lapply(ge.1, NA2mean)

pr<-prcomp(t(ge.1), center = TRUE, scale = TRUE)

tt=numeric();
for(i in 1:9)
{t=(pr$sdev[i])^2/sum((pr$sdev)^2);
tt=c(tt,t);
}
cumsum(tt)

color=rgb(1,0,0,alpha=0.3) 

plot(pr$x[,1],pr$x[,2],xlim=c(-15,17),ylim=c(-11,8),xlab=paste("PC1(",round(tt[1]*100,2),"%)",sep=""),ylab=paste("PC2(",round(tt[2]*100,2),"%)",sep=""),pch=19,col="brown",bty='l')
text(pr$x[,1]+1,pr$x[,2]+1,gsub("[^0-9A-Za-z///' ]","" , names(ge.1) ,ignore.case = TRUE),cex=0.9)
draw.ellipse(c(-12,-1,10), c(3.5,-7.5,3), c(4.5,5,7),c(5,4,5),angle=c(0,15,25),col=c(rgb(1,0,0,alpha=0.1),rgb(0,1,0,alpha=0.1),rgb(0,0,1,alpha=0.1)))

#mtext('A', side = 3, at = -16, line = 0.5, cex = 1.2);
dev.off()