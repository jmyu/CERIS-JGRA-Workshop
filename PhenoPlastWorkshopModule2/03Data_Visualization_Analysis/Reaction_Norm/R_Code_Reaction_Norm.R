rm(list=ls())
library(RColorBrewer)
library("plotrix")

Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/03Data_Visualization_Analysis/Reaction_Norm/";

png(paste(Dir,"Reaction_Norm.png",sep=""),width=6,height=4,pointsize=9,
    units="in",res=600)

mat=matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T);
layout(mat,widths=c(3,3),heights=c(2,2))

########Panel 1##########################################
par(mar=c(3,3,1.3,1),mgp=c(1.8,0.25,0),tcl=-0.2)

pheno=read.table(paste(Dir,"LbE_table",sep=""),header=T,sep="\t");
env_name=colnames(pheno)[-1]

meanY=apply(pheno[,-1],2,mean,na.rm=T)
names(meanY)=colnames(pheno)[-1]

slope=numeric();
for(i in 1:nrow(pheno))
{
  slo=lm(as.numeric(pheno[i,-1])~meanY,na.action=na.exclude)$coefficients[2];
  slope=c(slope,slo);
}

pheno$slope=slope;
pheno$col_order=rank(slope);
line_col_index <- pheno$col_order;
line_colors <- rep("gray",180)[line_col_index];

Env_order1=as.character(env_name);

plot(0, 0, col = "white", xlim = c(0.5, 9.5), ylim = c(50,170),  ylab = 'Flowering time (DAP)',  xlab = '', xaxt = "n", cex.lab = 1.1); ##, fg = "gray50"
axis(1, 1:9, labels = gsub("[^0-9A-Za-z///' ]","" , Env_order1 ,ignore.case = TRUE), cex.axis = .9,las=2);

for (i in 1:nrow(pheno)) 
{
  points(c(1:9), as.vector(pheno[i, Env_order1]), col = line_colors[i], type = "l", pch = 19, lwd = .8)
  points(c(1:9), as.vector(pheno[i, Env_order1]), col = line_colors[i],  pch = 19, cex = .3)
}

mtext('A', side = 3, at = 0.5, line = -0.1, cex = 1.2);

########Panel 2##########################################
Env_order2=names(meanY)[order(meanY)];

plot(0, 0, col = "white", xlim = c(0.5, 9.5), ylim = c(50,170),  ylab = 'Flowering time (DAP)',  xlab = '', xaxt = "n", cex.lab = 1.1); ##, fg = "gray50"
axis(1, 1:9, labels = gsub("[^0-9A-Za-z///' ]","" , Env_order2 ,ignore.case = TRUE), cex.axis = .9,las=2);

for (i in 1:nrow(pheno)) {
  points(c(1:9), as.vector(pheno[i, Env_order2]), col = line_colors[i], type = "l", pch = 19, lwd = .8)
  points(c(1:9), as.vector(pheno[i, Env_order2]), col = line_colors[i],  pch = 19, cex = .3)
}
mtext('B', side = 3, at = 0.5, line = -0.1, cex = 1.2);

#######Panel 3############################################
line_colors=colorRampPalette(c("Cyan","White","Magenta"))(180)[line_col_index]
Env_order3=names(meanY)[order(meanY)];

plot(0, 0, col = "white", xlim = range(meanY), ylim = c(50,170),  ylab = 'Flowering time (DAP)',  xlab = 'Environmental mean', cex.axis = 1, cex.lab = 1.1); 
for (i in 1:nrow(pheno)) {
  points(meanY, as.vector(pheno[i, Env_order3]), col = line_colors[i], type = "l", pch = 19, lwd = .8)
  points(meanY, as.vector(pheno[i, Env_order3]), col = line_colors[i],  pch = 19, cex = .3)
}

mtext('C', side = 3, at = meanY[1], line = -0.1, cex = 1.2);

text(meanY, rep(170, 9), gsub("[^0-9A-Za-z///' ]","" , Env_order3 ,ignore.case = TRUE),  srt = 90, cex = .8, adj = 1, family = "mono", font = 2)

########Panel 4##########################################
line_colors=colorRampPalette(c("Cyan","White","Magenta"))(180)[line_col_index]
Env_order4=names(meanY)[order(meanY)];

plot(0, 0, col = "white", xlim = range(meanY), ylim = c(50,170),  ylab = 'Flowering time (DAP)',  xlab = 'Environmental mean', cex.axis = 1, cex.lab = 1.1); # , fg = "gray50"
for (i in 1:nrow(pheno)) {
  fit=lm(as.numeric(pheno[i, Env_order4])~meanY,na.action=na.exclude);
  inter=coef(fit)[1];
  sloo=coef(fit)[2];
  points(meanY, inter+sloo*meanY, col = line_colors[i], type = "l", pch = 19, lwd = 1)
}

mtext('D', side = 3, at = meanY[1], line = -0.1, cex = 1.2);

text(meanY, rep(170, 9), gsub("[^0-9A-Za-z///' ]","" , Env_order4 ,ignore.case = TRUE),  srt = 90, cex = .8, adj = 1, family = "mono", font = 2)


dev.off()

