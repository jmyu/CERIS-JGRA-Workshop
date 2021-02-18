#################JGRA_through_Norm_Parameter_Estimation#####################
library("rrBLUP")
library("RColorBrewer")
library("BGLR")
library("yarrr")
Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/0607JGRA_Joint_Genomic_Regression_Analysis/JGRA_through_Norm_Parameter_Estimation/";
source(paste(Dir,"Function_JGRA.R",sep=""))

png(paste(Dir,"JGRA_Norm_CV.png",sep=""),width= 6,height= 6,pointsize=10,units="in",res=600)
mat=matrix(c(1:4),nrow=2,ncol=2,byrow=T);
layout(mat,widths=c(2.5,2.5),heights=c(2.5,2.5))
par(mar=c(3,3,2,0.5),mgp=c(2,0.5,0),tcl=-0.2)

#Panel0#######################################################################################################################
LOC=read.table(paste(Dir,"9Env_meta_table",sep=""),header=T,sep="\t");
pheno=read.table(paste(Dir,"LbE_table",sep=""),header=T,sep="\t");
envir=read.table(paste(Dir,"FTdap_envMeanPara_9_50",sep=""),header=T,sep="");

name_env=as.character(colnames(pheno)[-1]);
name_line=as.character(pheno$line_code);

envir=envir[match(name_env,envir$env_code),];
LOC=LOC[match(name_env,LOC$env_code),];

no_e=nrow(LOC);
no_l=nrow(pheno);

slope=numeric();
for(i in 1:nrow(pheno))
{
  slo=lm(as.numeric(pheno[i,-1])~envir$GDD,na.action=na.exclude)$coefficients[2];
  slope=c(slope,slo);
}

pheno$slope=slope;
pheno$col_order=rank(slope);
line_col_index <- pheno$col_order;
line_colors=colorRampPalette(c("Cyan","White","Magenta"))(180)[line_col_index]
Env_order4=envir$GDD;names(Env_order4)=envir$env_code;

plot(0, 0, col = "white", xlim = range(Env_order4), ylim = c(50,170),  ylab = 'Flowering time (DAP)',  xlab = 'Temperature', cex.axis = 1, cex.lab = 1.1); # , fg = "gray50"
for (i in 1:nrow(pheno)) {
  fit=lm(as.numeric(pheno[i, names(Env_order4)])~Env_order4,na.action=na.exclude);
  inter=coef(fit)[1];
  sloo=coef(fit)[2];
  points(Env_order4, inter+sloo*Env_order4, col = line_colors[i], type = "l", pch = 19, lwd = 1)
}

mtext('A', side = 3, at = 10, line = 0, cex = 1.2);

text(Env_order4-c(1,0.5,rep(0,7)), rep(170, 9), gsub("[^0-9A-Za-z///' ]","" , names(Env_order4),ignore.case = TRUE),  srt = 90, cex = .8, adj = 1, family = "mono", font = 2)

#Panel1#######################################################################################################################

LOC=read.table(paste(Dir,"9Env_meta_table",sep=""),header=T,sep="\t");
pheno=read.table(paste(Dir,"LbE_table",sep=""),header=F,sep="\t");
envir=read.table(paste(Dir,"FTdap_envMeanPara_9_50",sep=""),header=T,sep="");
geno=read.table(paste(Dir,"Genotype.txt",sep=""),header=T,sep="");

tt.line=300;##line quality
tt.e=6;##Envir quality
enp=4; ##column number of environmental index
fold=8;
reshuffle=20;
#1:2
out=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.E",fold,reshuffle)
out.1=out[[1]];
out.2=round(out[[2]][,1],2);
out.3=round(out[[3]],2);
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
r1 <- sprintf("%.2f", out.3);
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
mtext('B', side = 3, at = 47.5, line = 0, cex = 1.2);

#Panel2#######################################################################################################################
###############################################################1:3
out=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.G",fold,reshuffle)

out.1=out[[1]];
out.2=round(apply(out[[2]],2,mean),2);
out.3=round(mean(out[[3]]),2);
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
r1 <- sprintf("%.2f", out.3);
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
mtext('C', side = 3, at = 48, line = 0, cex = 1.2);

#Panel3#######################################################################################################################
###############################################################1:4
out=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE",fold,reshuffle)

out.1=out[[1]];
out.2=round(apply(out[[2]],2,mean),2);
out.3=round(mean(out[[3]]),2);
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
r1 <- sprintf("%.2f", out.3);
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
mtext('D', side = 3, at = 48, line = 0, cex = 1.2);


dev.off()







