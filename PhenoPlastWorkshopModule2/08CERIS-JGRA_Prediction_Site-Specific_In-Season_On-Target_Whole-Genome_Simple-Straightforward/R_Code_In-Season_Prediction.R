####Empirical validation####################################
####Untested genotypes in untested environments
library("rrBLUP")
library("RColorBrewer")
Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/08CERIS-JGRA_Prediction_Site-Specific_In-Season_On-Target_Whole-Genome_Simple-Straightforward/";

envir=read.table(paste(Dir,"Envir_Index.txt",sep=""),header=T,sep="")
phenotype=read.table(paste(Dir,"Phenotype.txt",sep=""),header=T,sep="");

phe.1=phenotype[,-1];

inter=numeric();
slo=numeric();

for(j in 1:239)
{
  y1= as.vector(t(phe.1[j,]));
  x1= as.vector(envir[,3]);
  coe=lm(y~x,data=data.frame(x=x1,y=y1));
  intercept=summary(coe)$coefficients[1,1];
  slope=summary(coe)$coefficients[2,1];
  inter=c(inter,intercept);
  slo=c(slo,slope);
}

####2015##########
yha=numeric();
for(i in 1:237)
{
  y_hat=inter[i]+slo[i]*354;
  yha=c(yha,y_hat);
}

Observed=read.table(paste(Dir,"Observed flowering time_2015.txt",sep=""),header=T,sep="")
cor(yha,Observed[,1],use="pairwise.complete.obs")
dat_2015=data.frame(pre=yha,obs=Observed[,1]);

####2016##########
yha=numeric();
for(i in 1:237)
{
  y_hat=inter[i]+slo[i]*416;
  yha=c(yha,y_hat);
}

Observed=read.table(paste(Dir,"Observed flowering time_2016.txt",sep=""),header=T,sep="")
cor(yha,Observed[,1],use="pairwise.complete.obs")
dat=data.frame(pre=yha,obs=Observed[,1]);
dat_2016=dat[!rowSums((is.na(dat))),]

#######################################Generate_Figures##############################################
#####################################################################################################
png(paste(Dir,"In-Season_Prediction.png",sep=""),width=3.5,height=3.5,pointsize=6,units="in",res=600)
mat=matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T);
layout(mat,widths=c(3.5/2,3.5/2),heights=c(3.5/2,3.5/2))
par(mar=c(4.0,4.0,2.0,2.0),mgp=c(2.0,0.5,0),tcl=-0.2)

##Fig1###
mycolor=rep(rev(brewer.pal(8,"Set1")[-6]),each=1)
y=as.numeric(apply(phenotype[,-1],2,mean));
x=envir[,3];
lm(formula = y ~ x)
plot(envir[,3],apply(phenotype[,-1],2,mean),ylim=c(1200,4000),col="white",pch=16,xlim=c(300,600),
     xlab="Photothermal time",ylab="Predicted flowering time",cex.lab=1,cex.axis=1)
#abline(294.979,3.816,lty=2,col="grey")

points(rep(354,237),dat_2015$pre,col="darkgreen",cex=0.3)
points(rep(416,233),dat_2016$pre,col="black",cex=0.3)

fix=3700;
size=1.0;
aj=c(0,12,0,0,4,0,0);
st=3930;ed=4000;
for(i in 1:7)
{
  points(rep(envir[i,3],237),inter[1:237]+slo[1:237]*envir[i,3],col="grey",cex=0.3);
  text(labels=envir[i,1], y = fix, x=envir[i,3]+aj[i],cex=size,srt=90,col="black");
  lines( y = c(st,ed), x=rep(x1[i],2));
}

lines(envir[,3],inter[238]+slo[238]*envir[,3],col="blue")
lines(envir[,3],inter[239]+slo[239]*envir[,3],col="red")
legend(300,3000,col=c("red","blue"),bty="n",legend=c("Parent1","Parent2"),cex=1,lty=1)

text(labels="IA15", y = fix, x=354,cex=size,srt=90)
text(labels="IA16", y = fix, x=416+4,cex=size,srt=90)
lines( y = c(st,ed), x=rep(354,2))
lines( y = c(st,ed), x=rep(416,2))

##Fig2###
plot(dat_2015$pre,dat_2015$obs,col="darkgreen",pch=16,cex=0.8,ylab="Observed flowering time",xlab="Predicted flowering time",
     cex.lab=1.0,cex.axis=1.0,xlim=c(1300,2500),ylim=c(1300,2500))
text(2300,1700,substitute(paste(italic('r'), " = 0.80")),cex=1.2)
text(1500,2400,"IA15",cex=1.2)
lines(x=c(1000,3000),y=c(1000,3000),col="grey")

##Fig3###
plot(dat_2016$pre,dat_2016$obs,col="black",pch=16,cex=0.8,ylab="Observed flowering time",xlab="Predicted flowering time",
     cex.lab=1.0,cex.axis=1.0,xlim=c(1300,2500),ylim=c(1300,2500))
text(2300,1700,substitute(paste(italic('r'), " = 0.88")),cex=1.2)
text(1500,2400,"IA16",cex=1.2)
lines(x=c(1000,3000),y=c(1000,3000),col="grey")
dev.off()





