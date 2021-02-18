#################JGRA_through_Norm_Parameter_Estimation#####################
library("rrBLUP")
library("RColorBrewer")
library("BGLR")
library("yarrr")
Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/JGRA_Joint_Genomic_Regression_Analysis/JGRA_through_Wide_Marker_Effect_Continuum/";

source(paste(Dir,"Function_JGRA_Marker.R",sep=""))

png(paste(Dir,"JGRA_Marker.png",sep=""),width= 6.5,height= 6.5,pointsize=10,units="in",res=600)

mat=matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T);
layout(mat,widths=c(2.5,2.5),heights=c(2.5,2.5))

par(mar=c(3,3,2,0.5),mgp=c(2,0.5,0),tcl=-0.2)
#Panel0#######################################################################################################################
########################################################
##marker effect###
LOC=read.table(paste(Dir,"9Env_meta_table",sep=""),header=T,sep="\t");
pheno=read.table(paste(Dir,"LbE_table",sep=""),header=T,sep="\t");
envir=read.table(paste(Dir,"FTdap_envMeanPara_9_50",sep=""),header=T,sep="");
geno=read.table(paste(Dir,"Genotype.txt",sep=""),header=T,sep="");

name_env=as.character(colnames(pheno)[-1]);
name_line=as.character(pheno$line_code);

envir=envir[match(name_env,envir$env_code),];
LOC=LOC[match(name_env,LOC$env_code),];
geno=geno[match(name_line,geno$RIL),];

no_e=nrow(LOC);
no_l=nrow(geno);
no_m=ncol(geno)-1;

effect=matrix(999,no_m,no_e);
intercept=matrix(999,1,no_e)

for(i in 1:9)
{
  fit=mixed.solve(pheno[,1+i],Z=geno[,-1])
  effect[,i]=(fit$u);
}

size=1.2;
x1=envir[,4];
plot(NULL,xlim=c(min(x1),max(x1)),xlab="Temperature",ylim=c(-3,4.5),
     ylab="Genomewide marker index",cex.lab=1.2,cex.axis=1.2,cex=1.2)

vv=numeric();
for(i in 1:162)
{
  Y.axis=effect[i,];
  reg=data.frame(x=as.numeric(x1),y=Y.axis);
  coe=lm(y~x,data=reg)
  intercept=summary(coe)$coefficients[1,1];
  slope=summary(coe)$coefficients[2,1];
  vv=c(vv,slope);
}

for(i in c(1:162))
{
  Y.axis=effect[i,];
  reg=data.frame(x=as.numeric(x1),y=Y.axis);
  coe=lm(y~x,data=reg)
  intercept=summary(coe)$coefficients[1,1];
  slope=summary(coe)$coefficients[2,1];
  lines(c(10,36),c(intercept+(slope)*10,intercept+(slope)*36),type="l",lwd=1.5,
        col=rev(colorRampPalette(c("Gold","Green"))(162)[cut(slope,breaks=unique(vv))]))
  #points(reg$x,reg$y,col=rev(colorRampPalette(c("Gold","Green"))(162)[cut(slope,breaks=unique(vv))]),cex=0.2,pch=20)
}

text(labels=gsub("[^0-9A-Za-z///' ]","" , envir$env_code[1:9] ,ignore.case = TRUE), y = rep(4.5,9), x=c(x1[1]-1.3,x1[2]-0.7,x1[3:9]),srt = 90, cex = 1, adj = 1, family = "mono", font = 2)
lines(c(0,600),c(0,0),lty=2,col="gray")
#lines(c(25.6,25.6),c(-20,20),lty=2,col="gray")


#Panel1#######################################################################################################################
#1:2
LOC=read.table(paste(Dir,"9Env_meta_table",sep=""),header=T,sep="\t");
pheno=read.table(paste(Dir,"LbE_table",sep=""),header=F,sep="\t");
envir=read.table(paste(Dir,"FTdap_envMeanPara_9_50",sep=""),header=T,sep="");
geno=read.table(paste(Dir,"Genotype.txt",sep=""),header=T,sep="");

tt.line=300;##line quality
tt.e=6;##Envir quality
enp=4;
out=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.E")
out.1=out[[1]];
out.2=round(out[[2]],2);
size=1;size.p=0.2;
plot(out.1$pre,out.1$obs,col=as.character(out.1$col),pch=16,cex=0.8,ylab="Observed flowering time",xlab="Predicted flowering time",
     cex.lab=size,cex.axis=size,xlim=c(50,190),ylim=c(50,190))

coll=as.character(unique(out.1$col));
nameE=gsub("[^0-9A-Za-z///' ]","" , names(out.2) ,ignore.case = TRUE);
op <- par(family = "mono")
legend(137,130,pch=16,col=coll,bty="n",
       legend=c(paste(nameE," (",out.2,")", sep = '')),
       cex=size)
par(op)
r1 <- sprintf("%.2f", cor(out.1$obs, out.1$pre, use = "complete.obs"));
text(160,180,substitute(paste(italic('r'), " = ", R1), list(R1 = r1)),cex=size+0.2)
lines(x=c(40,200),y=c(40,200),col="Gray")
bb=0.15;

#text(2400,3500,"Tested genotypes in untested environment (119)",cex=size)
xl=60;yb=130;xr=100;yu=170; 
rect(xl,yb,xr,yu,border="black",lwd=1)
lines(c(xl,xr),c((yb+yu)/2,(yb+yu)/2),col="black",lwd=1)
lines(c((xl+xr)/2,(xl+xr)/2),c(yb,yu),col="black",lwd=1)

text(xl+(xr-xl)/4,yu-(yu-yb)/4,"1",cex=size+0.5,col="black")
text(xr-(xr-xl)/4,yu-(yu-yb)/4,"2",cex=size+0.5,col="black")
text(xl+(xr-xl)/4,yb+(yu-yb)/4,"3",cex=size+0.5,col="Gray")
text(xr-(xr-xl)/4,yb+(yu-yb)/4,"4",cex=size+0.5,col="Gray")

text(xl+(xr-xl)/4,yu+(yu-yb)/4-5,"Tested",cex=size-bb)
text(xr-(xr-xl)/4,yu+(yu-yb)/4-5,"Untested",cex=size-bb)
text((xr+xl)/2,yu+(yu-yb)/2-8,"Environment",cex=size)

text(xl-(xr-xl)/4+5,yu-(yu-yb)/4,"Tested",cex=size-bb,srt=90)
text(xl-(xr-xl)/4+5,yb+(yu-yb)/4,"Untested",cex=size-bb,srt=90)
text(xl-(xr-xl)/2+10,(yu+yb)/2,"Genotype",cex=size,srt=90)
##adarrow##
x0=xl+(xr-xl)/4
y0=yu-(yu-yb)/4
x1=xr-(xr-xl)/4
y1=yu-(yu-yb)/4
arrows(x0+5,y0,x1-5,y1,angle=15,length = 0.07)
mtext('A', side = 3, at = 47.5, line = 0, cex = 1.2);

#Panel2#######################################################################################################################
###############################################################1:3
out=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.G")

out.1=out[[1]];
out.2=round(out[[2]],2);
size=1;size.p=0.2;
plot(out.1$pre,out.1$obs,col=as.character(out.1$col),pch=16,cex=0.8,ylab="Observed flowering time",xlab="Predicted flowering time",
     cex.lab=size,cex.axis=size,xlim=c(50,190),ylim=c(50,190))

coll=as.character(unique(out.1$col));
nameE=gsub("[^0-9A-Za-z///' ]","" , names(out.2) ,ignore.case = TRUE);
op <- par(family = "mono")
legend(137,130,pch=16,col=coll,bty="n",
       legend=c(paste(nameE," (",out.2,")", sep = '')),
       cex=size)
par(op)
r1 <- sprintf("%.2f", cor(out.1$obs, out.1$pre, use = "complete.obs"));
text(160,180,substitute(paste(italic('r'), " = ", R1), list(R1 = r1)),cex=size+0.2)
lines(x=c(40,200),y=c(40,200),col="Gray")


#text(2400,3500,"Tested genotypes in untested environment (119)",cex=size)
xl=60;yb=130;xr=100;yu=170; 
rect(xl,yb,xr,yu,border="black",lwd=1)
lines(c(xl,xr),c((yb+yu)/2,(yb+yu)/2),col="black",lwd=1)
lines(c((xl+xr)/2,(xl+xr)/2),c(yb,yu),col="black",lwd=1)

text(xl+(xr-xl)/4,yu-(yu-yb)/4,"1",cex=size+0.5,col="black")
text(xr-(xr-xl)/4,yu-(yu-yb)/4,"2",cex=size+0.5,col="Gray")
text(xl+(xr-xl)/4,yb+(yu-yb)/4,"3",cex=size+0.5,col="black")
text(xr-(xr-xl)/4,yb+(yu-yb)/4,"4",cex=size+0.5,col="Gray")

text(xl+(xr-xl)/4,yu+(yu-yb)/4-5,"Tested",cex=size-bb)
text(xr-(xr-xl)/4,yu+(yu-yb)/4-5,"Untested",cex=size-bb)
text((xr+xl)/2,yu+(yu-yb)/2-8,"Environment",cex=size)

text(xl-(xr-xl)/4+5,yu-(yu-yb)/4,"Tested",cex=size-bb,srt=90)
text(xl-(xr-xl)/4+5,yb+(yu-yb)/4,"Untested",cex=size-bb,srt=90)
text(xl-(xr-xl)/2+10,(yu+yb)/2,"Genotype",cex=size,srt=90)
##adarrow##
x0=xl+(xr-xl)/4
y0=yu-(yu-yb)/4
x1=xl+(xr-xl)/4
y1=yb+(yu-yb)/4
arrows(x0,y0-6,x1,y1+6,angle=15,length = 0.07)
mtext('B', side = 3, at = 48, line = 0, cex = 1.2);

#Panel3#######################################################################################################################
###############################################################1:4
out=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE")

out.1=out[[1]];
out.2=round(out[[2]],2);
size=1;size.p=0.2;
plot(out.1$pre,out.1$obs,col=as.character(out.1$col),pch=16,cex=0.8,ylab="Observed flowering time",xlab="Predicted flowering time",
     cex.lab=size,cex.axis=size,xlim=c(50,190),ylim=c(50,190))

coll=as.character(unique(out.1$col));
nameE=gsub("[^0-9A-Za-z///' ]","" , names(out.2) ,ignore.case = TRUE);
op <- par(family = "mono")
legend(137,130,pch=16,col=coll,bty="n",
       legend=c(paste(nameE," (",out.2,")", sep = '')),
       cex=size)
par(op)
r1 <- sprintf("%.2f", cor(out.1$obs, out.1$pre, use = "complete.obs"));
text(160,180,substitute(paste(italic('r'), " = ", R1), list(R1 = r1)),cex=size+0.2)
lines(x=c(40,200),y=c(40,200),col="Gray")


#text(2400,3500,"Tested genotypes in untested environment (119)",cex=size)
xl=60;yb=130;xr=100;yu=170; 
rect(xl,yb,xr,yu,border="black",lwd=1)
lines(c(xl,xr),c((yb+yu)/2,(yb+yu)/2),col="black",lwd=1)
lines(c((xl+xr)/2,(xl+xr)/2),c(yb,yu),col="black",lwd=1)

text(xl+(xr-xl)/4,yu-(yu-yb)/4,"1",cex=size+0.5,col="black")
text(xr-(xr-xl)/4,yu-(yu-yb)/4,"2",cex=size+0.5,col="Gray")
text(xl+(xr-xl)/4,yb+(yu-yb)/4,"3",cex=size+0.5,col="Gray")
text(xr-(xr-xl)/4,yb+(yu-yb)/4,"4",cex=size+0.5,col="black")

text(xl+(xr-xl)/4,yu+(yu-yb)/4-5,"Tested",cex=size-bb)
text(xr-(xr-xl)/4,yu+(yu-yb)/4-5,"Untested",cex=size-bb)
text((xr+xl)/2,yu+(yu-yb)/2-8,"Environment",cex=size)

text(xl-(xr-xl)/4+5,yu-(yu-yb)/4,"Tested",cex=size-bb,srt=90)
text(xl-(xr-xl)/4+5,yb+(yu-yb)/4,"Untested",cex=size-bb,srt=90)
text(xl-(xr-xl)/2+10,(yu+yb)/2,"Genotype",cex=size,srt=90)
##adarrow##
x0=xl+(xr-xl)/4
y0=yu-(yu-yb)/4
x1=xr-(xr-xl)/4
y1=yb+(yu-yb)/4
arrows(x0+3,y0-5,x1-3,y1+5,angle=15,length = 0.070)
mtext('C', side = 3, at = 48, line = 0, cex = 1.2);


dev.off()


