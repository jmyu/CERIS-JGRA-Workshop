if (!require(rrBLUP)) { install.packages("rrBLUP", repos = "https://cloud.r-project.org");}
if (!require(BGLR)) { install.packages("BGLR", repos = "https://cloud.r-project.org");}
if (!require(yarrr)) { install.packages("yarrr", repos = "https://cloud.r-project.org");}
if (!require(openxlsx)) { install.packages("openxlsx", repos = "https://cloud.r-project.org");}


Top_dir <- 'C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/04CERIS-Critical_Environmental_Regressor_through_Informed_Search/'; #### Modify this for your own directory
subfunction_file <- paste(Top_dir, 'Sub_functions_Workshop.r', sep = '');
source(subfunction_file);
######################################
experiment <- '2Rice'; ## 1Sorghum; 2Rice; 3Maize
######################################
exp_dir <- paste(Top_dir, experiment, '/', sep = '')
trait <- 'FTdap'; ### modify this correspondently
exp_trait_dir <- paste(exp_dir, trait,  '/',  sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir, recursive= T)};

maxR_dap1 <- 9;
maxR_dap2 <- 50;
kPara_Name <- 'GDD';

Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/0607JGRA_Joint_Genomic_Regression_Analysis/"
#########################################################
############ From Tingting Guo (tguo@iastate.edu) #################
#### 3 prediction scenarios: 1->2; 1->3; 1->4
LOC=read.table(paste(exp_dir,"Env_meta_table.txt",sep=""),header=T,sep="\t");
geno=read.table(paste(exp_dir,"Genotype.txt",sep=""),header=T,sep="\t");
pheno=read.table(paste(exp_trait_dir,"LbE_table.txt",sep=""),header=F,sep="\t");
envir=read.table(paste(exp_trait_dir,trait,'_envMeanPara_', maxR_dap1, '_', maxR_dap2, '.txt',sep=""),header=T,sep="\t");

tt.line=nrow(pheno)*0.5;##remove environment if the number of missing lines > tt.line
tt.e=c(ncol(pheno)-1)-3;##remove line if the number of missing environment > tt.e
enp=which(colnames(envir) == kPara_Name); 
fold=5;
reshuffle=10;

###Throught reaction norm parameter
out1.4=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE",fold,reshuffle)  ## 1->4 prediction
write.xlsx(out1.4, file = paste(Dir,"Result_Norm_1.4.xlsx",sep=""))
### Throught marker effect
out1.4=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE",fold,reshuffle)
write.xlsx(out1.4, file = paste(Dir,"Result_Marker_1.4.xlsx",sep=""))

