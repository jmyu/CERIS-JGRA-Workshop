#install.packages("VCA")
library("VCA")
Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/02Partition_GxE_Heterogeneity_of_Genotypic_Variance_Lack_of_Genetic_Correlation/";

Input_VCA=read.table(paste(Dir,"Input_Phenotypes.txt",sep=""),header=T)
Input_VCA[1:4,]

Input_VCA$Loc <- as.factor(Input_VCA$Loc)

A=as.character(Input_VCA$Loc);
B=as.character(Input_VCA$Entry);
C=as.character(Input_VCA$Rep);
y=as.numeric(as.character(Input_VCA$GDD));

mydataframe=data.frame(y,A,B,C);

###For variance component analysis
fit_r <-remlMM(y ~ (A)+(B)+(A/C)+(A*B), mydataframe)
fit_r

GxE_Var=fit_r$aov.tab[5,2];
G_Var=fit_r$aov.tab[3,2];
Error_Var=fit_r$aov.tab[6,2];

###For ANOVA (Analysis of variance)
fit_a <-anovaMM(y ~ A+B+A/C+A*B, mydataframe)
fit_a

#A linear model for describing the phenotypic performance of genotypes in individual environments.
###
EN=unique(A);
GenoVarInd=numeric();
for(i in 1:7)
{
  data1=mydataframe[mydataframe$A==EN[i],];
  data1=data1[which(!(is.na(data1$y))),];
  fit1 <-remlMM(y ~ (B)+(C), data1)
  print(fit1)
  #Variance component for genotypes###
  GenoVar=fit1$aov.tab[2,2];
  GenoVarInd=c(GenoVarInd,GenoVar)
}

#V
V=sum((sqrt(GenoVarInd)-mean(sqrt(GenoVarInd)))^2)/(length(EN)-1)
V.percent=V/GxE_Var;
V
V.percent
#L
L=GxE_Var-V;
L.percent=L/GxE_Var;
L
L.percent

#Line/Entry mean heritability
h2=G_Var/(G_Var+GxE_Var/length(EN)+Error_Var/(length(EN)*2));
h2
#Pooled genetic correlation
rg=G_Var/(G_Var+L)
rg

