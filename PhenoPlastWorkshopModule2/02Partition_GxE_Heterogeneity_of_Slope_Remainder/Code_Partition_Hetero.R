####Heterogeity of slope ###############
Dir="C:/Users/tguo/Desktop/PhenoPlastWorkshopModule2/02Partition_GxE_Heterogeneity_of_Slope_Remainder/";

LBE=read.table(file=paste(Dir,"LBE_Table.txt",sep=""),header=T,sep="\t")
LBE.matrix=as.matrix(LBE[,-1]);

E.filter=apply(LBE.matrix,2,function(x) length(which(!is.na(x))));
G.filter=apply(LBE.matrix,1,function(x) length(which(!is.na(x))));

LBE.matrix.1=LBE.matrix[G.filter>=5,E.filter>=200];

pop.mean=apply(LBE.matrix.1,2,mean,na.rm=T)

vv=numeric();
for(i in 1:nrow(LBE.matrix.1))
{
  v=predict(lm(LBE.matrix.1[i,]~pop.mean))
  vv=rbind(vv,v)
}

mu=mean(LBE.matrix.1,na.rm=T)
ipsen=apply(LBE.matrix.1,2,mean,na.rm=T)-mu
mu.ipsen=apply(LBE.matrix.1,2,mean,na.rm=T)
ipsen.sqr=ipsen^2

d=apply(LBE.matrix.1,1,mean,na.rm=T)-mu
beta.1=apply(LBE.matrix.1,1,function(x) lm(x ~ ipsen)$coefficients[2])
d.sqr=d^2

LBE.matrix.2=sweep(LBE.matrix.1,2,mu.ipsen)
beta.2=apply(LBE.matrix.2,1,function(x) lm(x ~ ipsen)$coefficients[2])
beta.2.sqr=beta.2^2

output=as.data.frame(matrix(NA,5,5));
row.names(output)=c("E","G","GE","error","total");
colnames(output)=c("degrees of freedom","sums of squares","means of squares","F","probability")

t=nrow(LBE.matrix.1)
s=ncol(LBE.matrix.1)


LBE.matrix.3=(LBE.matrix.1-mu)^2
data=LBE.matrix.3;
for(i in 1:ncol(data)){
  data[is.na(data[,i]), i] <- mean(data[,i], na.rm = TRUE)
}




output[1,1]=s-1;
output[2,1]=t-1;
output[3,1]=t-1;
output[5,1]=s*t-1;
output[4,1]=(t-1)*(s-2)

output[1,2]=sum(ipsen.sqr)*t;
output[2,2]=sum(d.sqr)*s;
output[3,2]=sum(beta.2.sqr)*sum(ipsen.sqr);
output[5,2]=sum(data,na.rm=T);
output[4,2]=output[5,2]-sum(output[1:3,2])

output[1,3]=output[1,2]/output[1,1];
output[2,3]=output[2,2]/output[2,1];
output[3,3]=output[3,2]/output[3,1];
output[4,3]=output[4,2]/output[4,1];

output[1,4]=output[1,3]/output[4,3];
output[2,4]=output[2,3]/output[4,3];
output[3,4]=output[3,3]/output[4,3];

output[1,5]=1-pf(output[1,4], df1=output[1,1], df2=output[4,1])
output[2,5]=1-pf(output[2,4], df1=output[2,1], df2=output[4,1])
output[3,5]=1-pf(output[3,4], df1=output[3,1], df2=output[4,1])

write.table(output,file=paste(Dir,"Output_Table.txt",sep=""),quote=F,row.names = T)









