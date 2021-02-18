library("yarrr")

Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/03Data_Visualization_Analysis/Environmental_Factors/"

png(paste(Dir,"Environmental_Factors.png",sep=""),width=6,height=6/2,pointsize=8,
    units="in",res=600)

par(mfrow = c(1,2),
    mar=c(3,4,1.5,1.5),mgp=c(1.8,0.25,0),tcl=-0.2)

Envir.1=read.table(file=paste(Dir,"9Envs_envParas_DAP200",sep=""),header=T,sep="\t");
pheno=read.table(file=paste(Dir,"Trait_Rice_addFTgdd",sep=""),header=T,sep="\t")

Name=c("TS_07","TS_08E","TS_08L","TS_09","ISA_08","FU_08","ISI_08","TH_08","HA_08");

colo=piratepal(palette="basel",trans=0)

###Daylength###
plot(NULL,xlim=c(0,200),ylim=c(12,16),ylab="Day length (hours)",xlab="Days after planting (DAP)",
     main="Photoperiod");

for(i in 1:9)
{
  E1=Envir.1[Envir.1$env_code%in%Name[i],]
  P1=pheno[pheno$env_code%in%Name[i],];
  lines(1:200,E1$DL,lty=1,col="gray")
  lines(9:50,E1$DL[9:50],lwd=3,col=colo[i])
  points(round(mean(P1$FTdap,na.rm=T),0),E1$DL[round(mean(P1$FTdap,na.rm=T),0)],lwd=1.5,col=colo[i])
}

#legend(130,16,legend=Name,lty=1,col=colo,bty = "n",cex=0.5,y.intersp=0.5)
mtext('A', side = 3, at = -40, line = -0.1, cex = 1.2);
###GDD###

plot(NULL,xlim=c(0,200),ylim=c(0,4000),ylab="Cumulative growing degree days (GDD)",xlab="Days after planting (DAP)",
     main="Temperature");

for(i in 1:9)
{
  E1=Envir.1[Envir.1$env_code%in%Name[i],];
  P1=pheno[pheno$env_code%in%Name[i],];
  day.range=range(P1$FTdap,na.rm = T);
  gdd.range=range(P1$FTgdd,na.rm = T);
  
  miss <- is.na(E1$GDD)
  E1$GDD[miss] <- 0
  cs <- cumsum(E1$GDD)
  cs[miss] <- NA
  
  lines(1:200,cs,lty=1,col="gray")
  lines(9:50,cs[9:50],lwd=1.5,col=colo[i])
  points(round(mean(P1$FTdap,na.rm=T),0),cs[round(mean(P1$FTdap,na.rm=T),0)],lwd=1.5,col=colo[i])
}

legend(130,1900,legend=gsub("[^0-9A-Za-z///' ]","" , Name ,ignore.case = TRUE),lty=1,col=colo,bty = "n",cex=0.9,y.intersp=1,pt.cex=1.5)
mtext('B', side = 3, at = -40, line = -0.1, cex = 1.2);
dev.off()



