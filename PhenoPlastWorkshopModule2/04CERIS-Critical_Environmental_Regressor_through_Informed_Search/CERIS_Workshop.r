###### Two sectioons: the first section is to identify the environmetnal index (by Xianran Li); and the second section is for performance prediction (by Tingting Guo)
## install required packages
if (!require(rnoaa)) { install.packages("rnoaa", repos = "https://cloud.r-project.org");}
if (!require(dplyr)) { install.packages("dplyr", repos = "https://cloud.r-project.org");}
if (!require(R.utils)) { install.packages("R.utils", repos = "https://cloud.r-project.org");}
if (!require(corrgram)) { install.packages("corrgram", repos = "https://cloud.r-project.org");}
if (!require(lubridate)) { install.packages("lubridate", repos = "https://cloud.r-project.org");}
if (!require(geosphere)) { install.packages("geosphere", repos = "https://cloud.r-project.org");}
if (!require(data.table)) { install.packages("data.table", repos = "https://cloud.r-project.org");}
if (!require(colorspace)) { install.packages("colorspace", repos = "https://cloud.r-project.org");}
if (!require(RColorBrewer)) { install.packages("RColorBrewer", repos = "https://cloud.r-project.org");}

if (!require(rrBLUP)) { install.packages("rrBLUP", repos = "https://cloud.r-project.org");}
if (!require(BGLR)) { install.packages("BGLR", repos = "https://cloud.r-project.org");}
if (!require(yarrr)) { install.packages("yarrr", repos = "https://cloud.r-project.org");}
if (!require(openxlsx)) { install.packages("openxlsx", repos = "https://cloud.r-project.org");}

col_wdw <- 25;
col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1)
t_base <- 50; t_max1 <- 100; t_max2 <- 1000; Haun_threshold <- 0; p <- 1 

###
Top_dir <- 'C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/04CERIS-Critical_Environmental_Regressor_through_Informed_Search/'; #### Modify this for your own directory

###### If you modifity some funcition in this file, please make sure to run Line 27 each time to reload the updated functions 
subfunction_file <- paste(Top_dir, 'Sub_functions_Workshop.r', sep = '');
source(subfunction_file);

######################################
experiment <- '3Maize'; ## 1Sorghum; 2Rice; 3Maize
######################################

exp_dir <- paste(Top_dir, experiment, '/', sep = '')
env_meta_file <- paste(exp_dir, 'Env_meta_table.txt', sep = ''); ## make sure the PlantingData formated as 'YYYY-MM-DD'
env_meta_info_0 <- read.table(env_meta_file, header = T, sep = "\t", stringsAsFactors = F);

exp_s_year <- min(env_meta_info_0$TrialYear); 
exp_e_year <- max(env_meta_info_0$TrialYear) + 1; 

searching_daps <- 80; if (experiment == '1Sorghum') { searching_daps <- 122};
###################################################
### 
### Download the year gz file form ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/by_year/ and put the downloaded file (such as 2007.csv.gz) in this directory.
ghcn_year_dir <- paste(Top_dir, '0ghcn/', sep = ''); if (!dir.exists(ghcn_year_dir))  { dir.create(ghcn_year_dir)}; ###
###################################################
ptt_ptr_file <- paste(exp_dir, nrow(env_meta_info_0), 'Envs_envParas_DAP',  searching_daps, '.txt.', sep = ''); ##  
if (!file.exists(ptt_ptr_file)) { Compile_PTT_PTR_local_GHCN(exp_dir, env_meta_info_0, exp_s_year, exp_e_year, searching_daps, ptt_ptr_file, t_base, t_max1, t_max2, Top_dir) };
PTT_PTR <- read.table(ptt_ptr_file, header = T , sep = "\t");
Paras <- c('DL', 'GDD', 'DTR','PTT', 'PTR', 'PTD', 'PTD2', 'PTS');

exp_traits_file <- paste(exp_dir, 'Trait_records.txt', sep = '');
exp_traits <- read.table(exp_traits_file, sep = "\t", header = T, stringsAsFactors = F, na.string = 'NA');

if(!('FTdap' %in% colnames(exp_traits))) {exp_traits$FTdap <- exp_traits$DTA};

all_env_codes <- unique(exp_traits$env_code);
env_cols <- rainbow_hcl(length(all_env_codes), c = 80, l = 60, start = 0, end = 300, fixup = TRUE, alpha = 0.75);

##### For Sorghum: the trait is called as FTgdd; for Rice, the trait is FTdap; for Maize, the traits are DTA, PH, and YLD; 
trait <- 'PH'; ### modify this correspondently
lInd <- which(colnames(exp_traits) == 'line_code'); eInd <- which(colnames(exp_traits) == 'env_code'); tInd <- which(colnames(exp_traits) == trait);
exp_trait_dir <- paste(exp_dir, trait,  '/',  sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir, recursive= T)};
exp_trait <- exp_traits[,c(lInd, eInd, tInd)]; 

colnames(exp_trait)[3] <- 'Yobs';
exp_trait <- aggregate(Yobs ~  line_code + env_code, data = exp_trait, mean) ## To make sure only one phenotype record per each line in each environment
exp_trait <- exp_trait[!is.na(exp_trait$Yobs),];

line_codes <- unique(exp_trait$line_code); 
env_mean_trait_0 <- na.omit(aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code), mean, na.rm = T));
colnames(env_mean_trait_0)[2] <- 'meanY';
env_mean_trait <- env_mean_trait_0[order(env_mean_trait_0$meanY),];

### pairwise correlations among enviroments; trait distribution across enviroments;
### two figures and the correspondent output files will be saved in the trait directory;
try(Pairwise_trait_env_distribution_plot(exp_trait, exp_trait_dir, trait, all_env_codes, env_meta_info_0));

##### searching the critical window with the highest correlation with environmental mean
##### the window can be adjusted based on biological meaning
##### 'FTgdd_7Envs_PTTPTR_0LOO_cor.txt' stores all correlations from all the tested windows and environmental paramters;
##### 'MaxR_FTgdd_7Envs_0LOO.png' is the visulization 
pop_cor_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', 0, 'LOO_cor.txt', sep = '');
Exhaustive_search(env_mean_trait, PTT_PTR, searching_daps, exp_trait_dir, exp_traits$FTdap, trait, 1, searching_daps, searching_daps, 0, Paras, pop_cor_file)#; searching_daps, searching_daps);

#######################################################
##### From Laura Cortes (ltibbs@iastate.edu) ##########
##### To help navigate the strongest correlations identified for each environment paramter. 
##### View the results with 'FTgdd_7Envs_PTTPTR_0LOO_cor.txt' and 'MaxR_FTgdd_7Envs_0LOO.png'
##### R_PTT means the original correlation between PTT and population mean; while nR_PTT means the correlation of 0 - R_PTT. 
search.results <- read.table(pop_cor_file, header = T)
search.results <- search.results %>%
  tidyr::gather(key="Parameter", value="Corr", -Day_x, -Day_y, -window, -pop_code) %>%
  arrange(-Corr)
search.results <- search.results %>%
  group_by(Parameter) %>%
  top_n(5, Corr)    
View(search.results);
#############################################

### change the following three parameters for the window and environmental parameter with the strongest correlation
maxR_dap1 <- 22;
maxR_dap2 <- 35;
kPara_Name <- 'PTS';
#####  
#####
PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); 
#### Visulization of the relationships between environmental mean and environmental parameters from the selected window.  
Plot_Trait_mean_envParas(env_mean_trait, PTT_PTR, maxR_dap1, maxR_dap2, trait, exp_trait_dir, env_cols, Paras);  

#### Output intercept and slope estimation for each line based on environmental mean and environmental parameter
Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, exp_trait, line_codes, exp_trait_dir);
#### LOOCV function for 1 -> 2
obs_prd_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, 'D', maxR_dap1, '_', maxR_dap2, '.txt', sep = '');
LOO_pdf_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, 'D', maxR_dap1, '_', maxR_dap2, '.png', sep = '');
if (!file.exists(obs_prd_file)) { 
  prdM <- LOOCV(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, PTT_PTR_ind, exp_trait, obs_prd_file, p)
}
Plot_prediction_result(obs_prd_file, all_env_code, prdM, kPara_Name, LOO_pdf_file,env_cols);

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
fold=2;
reshuffle=1;

###Throught reaction norm parameter
out1.2=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.E",fold,reshuffle)  ## 1->2 prediction
out1.3=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.G",fold,reshuffle)  ## 1->3 prediction
out1.4=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE",fold,reshuffle) ## 1->4 prediction
write.xlsx(out1.2, file = paste(exp_trait_dir,"Result_Norm_1.2.xlsx",sep=""))
write.xlsx(out1.3, file = paste(exp_trait_dir,"Result_Norm_1.3.xlsx",sep=""))
write.xlsx(out1.4, file = paste(exp_trait_dir,"Result_Norm_1.4.xlsx",sep=""))
### Throught marker effect
out1.2=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.E",fold,reshuffle)
out1.3=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.G",fold,reshuffle)
out1.4=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE",fold,reshuffle)
write.xlsx(out1.2, file = paste(exp_trait_dir,"Result_Marker_1.2.xlsx",sep=""))
write.xlsx(out1.3, file = paste(exp_trait_dir,"Result_Marker_1.3.xlsx",sep=""))
write.xlsx(out1.4, file = paste(exp_trait_dir,"Result_Marker_1.4.xlsx",sep=""))

##For ploting
O1=out1.2[[1]];
n.envir=length(unique(O1$col));
n.line=length(O1$col[which(O1$col==O1$col[1])]);

coloo=piratepal(palette="basel",trans=0)[1:n.envir];
col=rep(coloo,times=rep(n.line,n.envir))

plot(O1$pre,O1$obs,col=col,pch=16,cex=0.8,ylab="Observed flowering time",xlab="Predicted flowering time",xlim=range(O1[,1:2],na.rm=T),
     ylim=range(O1[,1:2],na.rm=T))
abline(coef = c(0,1))
text(mean(par("usr")[1:2]),par("usr")[2]-diff(par("usr")[1:2])*0.1,paste("r_between = ",round(cor(O1$pre,O1$obs,use = "complete.obs"),2),sep=""))





