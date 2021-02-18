#install.packages("agricolae")
library(agricolae)
Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/02Partition_GxE_PCs/";

Pheno=read.table(paste(Dir,"Input_Phenotypes.txt",sep=""),header=T,sep="\t");

result=AMMI(Pheno$Loc, Pheno$Entry, Pheno$Rep, as.numeric(as.character(Pheno$GDD)))
PCA=result$analysis

write.table(PCA,file=paste(Dir,"AMMI_PCA.txt",sep=""),quote=F,row.names = T)
