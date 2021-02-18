
######## From Tingting Guo ######
JGRA=function(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets,fold=10,reshuffle=50)
{
  
  colnames(pheno)=as.character(unlist(pheno[1,]));
  envir.name=colnames(pheno)[-1];
  
  pheno=pheno[-1,];
  m=pheno[,-1];
  m=as.numeric(as.character(unlist(m)));m <- matrix(data=m, ncol=dim(pheno)[2]-1, nrow=dim(pheno)[1]);
  colnames(m)=colnames(pheno)[-1];
  pheno_=data.frame(line_code=pheno$line_code,m);
  colnames(pheno_)=c("line_code",envir.name);
  pheno=pheno_;
  
  envir=envir[envir$env_code%in%envir.name,];
  LOC=LOC[LOC$env_code%in%envir.name,];
  
  order=match(envir$env_code,LOC$env_code);
  LOC=LOC[order,];
  #pheno=pheno[,as.character(envir$env_code)];
  
  rm.env=colnames(pheno)[colSums(is.na(pheno))>=tt.line];
  LOC=LOC[which(!(LOC$env_code%in%rm.env)),];
  envir=envir[which(!(envir$env_code%in%rm.env)),];
  pheno=pheno[,which(!(colnames(pheno)%in%rm.env))];
  
  pheno.1=pheno
  keep=dim(pheno.1)[2]-rowSums(is.na(pheno.1))
  pheno=pheno.1[which(keep>=tt.e),];
  n.line=dim(pheno)[1];
  n.envir=dim(pheno)[2]-1;
  
  library("yarrr")
  coloo=piratepal(palette="basel",trans=0)[1:n.envir];
  #coloo=heat.colors(n.envir);
  
  if(mets=="RM.E") 
  {
    pheno.hat=matrix(999,n.line,dim(envir)[1]);
    cor.whole=numeric();
    for(k in 1:n.envir)
    {
      for(j in 1:n.line)
      {
        x1=envir[,enp][-k];
        y1=as.vector(t(pheno[j,-c(1,1+k)]));
        
        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[,enp][k];
        pheno.hat[j,k]=y.hat;
      }
      cor.each=cor(pheno.hat[,k],pheno[,k+1], use = "complete.obs");
      cor.whole=c(cor.whole,cor.each);
    }
    observe=as.vector(as.matrix(pheno[,-1]));
    predict=as.vector(pheno.hat);
    
    r_within=cor.whole;names(r_within)=colnames(pheno)[-1];
    r_within=data.frame(cor_within=r_within,envir=colnames(pheno)[-1]);
    r_across=cor(observe,predict,use = "complete.obs");
    outforfigure=data.frame(obs=observe,pre=predict,col=rep(coloo,times=rep(n.line,n.envir)),envir=rep(colnames(pheno)[-1],times=rep(n.line,n.envir)))
  }
  
  if(mets=="RM.G") 
  {
    intercept=numeric();
    slope=numeric();
    for(j in 1:n.line)
    {
      x1=envir[,enp];
      y1=as.vector(t(pheno[j,-c(1)]));
      
      coe=lm(y~x,data=data.frame(x=x1,y=y1));
      inter=summary(coe)$coefficients[1,1]
      slop=summary(coe)$coefficients[2,1];
      intercept=c(intercept,inter);
      slope=c(slope,slop);
    }
    
    genotype.match=match(pheno[,1],geno[,1])
    genotype=geno[genotype.match,];
    Marker=genotype[,-1];
    
    intercept.hat=numeric();slope.hat=numeric();cor.within=matrix(999,reshuffle,n.envir);cor.all=numeric();
    for(i in 1:reshuffle)
    { 
      
      cross=sample(rep(1:fold,each=ceiling(n.line/fold)),n.line);
      yhat.whole.cross=numeric();yobs.whole.cross=numeric();
      for(f in 1:fold)
      {
        id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
        ##Intercept###
        y0=intercept; 
        ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
        e=as.matrix(ans$u)
        G.pred=Marker[id.V,]
        y_pred=as.matrix(G.pred) %*% e
        GEBV.inter=c(y_pred[,1])+c(ans$beta);
        ##Slope###
        y0=slope; 
        ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
        e=as.matrix(ans$u)
        G.pred=Marker[id.V,]
        y_pred=as.matrix(G.pred) %*% e
        GEBV.slope=c(y_pred[,1])+c(ans$beta);
        ###All the predicted slope and intercept
        yhat.envir=matrix(999,length(id.V),n.envir);yobs.envir=matrix(999,length(id.V),n.envir)
        for(j in 1:n.envir)
        {
          yhat=GEBV.inter+GEBV.slope*envir[j,enp];
          yobs=pheno[id.V,j+1];
          yhat.envir[,j]=yhat;yobs.envir[,j]=yobs;
        }
        yhat.whole.cross=rbind(yhat.whole.cross,yhat.envir);
        yobs.whole.cross=rbind(yobs.whole.cross,yobs.envir);
      }
      for(j in 1:n.envir)
      {cor.within[i,j]=cor(yhat.whole.cross[,j],yobs.whole.cross[,j],use = "complete.obs");}
      
      cor.all=c(cor.all,cor(as.vector(yhat.whole.cross),as.vector(yobs.whole.cross),use = "complete.obs"));
    }
    #Correlation within environment 50 times
    r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
    r_within=data.frame(cor_within=r_within,envir=colnames(pheno)[-1]);
    #Correlation across environment 50 times
    r_across=mean(cor.all);
    #Observation and prediction last time
    outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                            pre=as.vector(yhat.whole.cross),
                            col=rep(coloo,times=rep(n.line,n.envir)),envir=rep(colnames(pheno)[-1],times=rep(n.line,n.envir)))
    colnames(cor.within)=colnames(pheno)[-1];
    r_within=cor.within;
    r_across=cor.all;
  }   
  
  if(mets=="RM.GE") 
  { 
    genotype.match=match(pheno[,1],geno[,1])
    genotype=geno[genotype.match,];
    Marker=genotype[,-1];
    
    cor.within=matrix(999,reshuffle,n.envir);cor.all=numeric();
    for(i in 1:reshuffle)
    {
      obs_matrix=matrix(999,n.line,n.envir);pre_matrix=matrix(999,n.line,n.envir);
      for(k in 1:n.envir)
      {
        intercept=numeric();
        slope=numeric();
        for(j in 1:n.line)
        {
          x1=envir[-k,enp];
          y1=as.vector(t(pheno[j,-c(1,1+k)]));
          
          coe=lm(y~x,data=data.frame(x=x1,y=y1));
          inter=summary(coe)$coefficients[1,1]
          slop=summary(coe)$coefficients[2,1];
          intercept=c(intercept,inter);
          slope=c(slope,slop);
        }
        
        cross=sample(rep(1:fold,each=ceiling(n.line/fold)),n.line);
        yhat.whole=numeric();yobs.whole=numeric();
        
        for(f in 1:fold)
        {
          id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
          ##Intercept###
          y0=intercept; 
          ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
          e=as.matrix(ans$u)
          G.pred=Marker[id.V,]
          y_pred=as.matrix(G.pred) %*% e
          GEBV.inter=c(y_pred[,1])+c(ans$beta);
          ##Slope###
          y0=slope; 
          ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
          e=as.matrix(ans$u)
          G.pred=Marker[id.V,]
          y_pred=as.matrix(G.pred) %*% e
          GEBV.slope=c(y_pred[,1])+c(ans$beta);
          ###All the predicted slope and intercept
          yhat=GEBV.inter+GEBV.slope*envir[k,enp];
          yobs=pheno[id.V,k+1];
          
          yhat.whole=c(yhat.whole,yhat);
          yobs.whole=c(yobs.whole,yobs);
        }
        cor.within[i,k]=cor(yhat.whole,yobs.whole,use = "complete.obs");
        obs_matrix[,k]=yobs.whole;
        pre_matrix[,k]=yhat.whole;
      }
      cor.shuffle=cor(as.vector(obs_matrix),as.vector(pre_matrix),use = "complete.obs")
      cor.all=c(cor.all,cor.shuffle);
    }
    
    yhat.whole.cross=pre_matrix;
    yobs.whole.cross=obs_matrix;
    #Correlation within environment 50 times
    r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
    #Correlation across environment 50 times
    r_across=mean(cor.all);
    #Observation and prediction last time
    outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                            pre=as.vector(yhat.whole.cross),
                            col=rep(coloo,times=rep(n.line,n.envir)),envir=rep(colnames(pheno)[-1],times=rep(n.line,n.envir)))
    colnames(cor.within)=colnames(pheno)[-1];
    r_within=cor.within;
    r_across=cor.all;
    
  }
  
  return(list(outforfigure,r_within,r_across));
}
JGRA.marker=function(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets,fold=10,reshuffle=50)
{
  colnames(pheno)=as.character(unlist(pheno[1,]));
  envir.name=colnames(pheno)[-1];
  
  pheno=pheno[-1,];
  m=pheno[,-1];
  m=as.numeric(as.character(unlist(m)));m <- matrix(data=m, ncol=dim(pheno)[2]-1, nrow=dim(pheno)[1]);
  colnames(m)=colnames(pheno)[-1];
  pheno_=data.frame(line_code=pheno$line_code,m);
  colnames(pheno_)=c("line_code",envir.name);
  pheno=pheno_;
  
  envir=envir[envir$env_code%in%envir.name,];
  LOC=LOC[LOC$env_code%in%envir.name,];
  
  order=match(envir$env_code,LOC$env_code);
  LOC=LOC[order,];
  #pheno=pheno[,as.character(envir$env_code)];
  
  rm.env=colnames(pheno)[colSums(is.na(pheno))>=tt.line];
  LOC=LOC[which(!(LOC$env_code%in%rm.env)),];
  envir=envir[which(!(envir$env_code%in%rm.env)),];
  pheno=pheno[,which(!(colnames(pheno)%in%rm.env))];
  
  pheno.1=pheno
  keep=dim(pheno.1)[2]-rowSums(is.na(pheno.1))
  pheno=pheno.1[which(keep>=tt.e),];
  n.line=dim(pheno)[1];
  n.envir=dim(pheno)[2]-1;
  
  library("yarrr")
  coloo=piratepal(palette="basel",trans=0)[1:n.envir];
  #coloo=heat.colors(n.envir);
  
  genotype.match=match(pheno[,1],geno[,1])
  genotype=geno[genotype.match,];
  Marker=genotype[,-1];
  n.marker=dim(Marker)[2];
  
  if(mets=="RM.E") 
  {
    
    effect=matrix(999,n.marker,n.envir);
    intercept=matrix(999,1,n.envir)
    
    for(i in 1:n.envir)
    {
      fit=mixed.solve(pheno[,1+i],Z=Marker)
      effect[,i]=fit$u
      intercept[,i]=fit$beta
    }
    
    pheno.hat=matrix(999,n.line,dim(envir)[1]);
    cor.whole=numeric();
    for(k in 1:n.envir)
    {
      
      effect.hat=numeric();
      for(j in 1:n.marker)
      {
        x1=envir[,enp][-k];
        y1=effect[j,-k];
        
        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[,enp][k];
        effect.hat=c(effect.hat,y.hat);
      }
      
      ##Environment intercept####
      reg.intercept=data.frame(x=as.numeric(envir[,enp][-k]),y=as.vector(intercept[-k]));
      coe.intercept=lm(y~x,data=reg.intercept)
      y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(envir[,enp][k]);
      
      pheno.hat[,k]=y.intercept+as.matrix(Marker)%*%as.matrix(effect.hat);
      cor.each=cor(pheno.hat[,k],pheno[,k+1], use = "complete.obs");
      cor.whole=c(cor.whole,cor.each);
    }
    observe=as.vector(as.matrix(pheno[,-1]));
    predict=as.vector(pheno.hat);
    r_within=cor.whole;names(r_within)=colnames(pheno)[-1];
    r_within=data.frame(cor_within=r_within,envir=colnames(pheno)[-1]);
    r_across=cor(observe,predict,use = "complete.obs");
    outforfigure=data.frame(obs=observe,pre=predict,col=rep(coloo,times=rep(n.line,n.envir)),envir=rep(colnames(pheno)[-1],times=rep(n.line,n.envir)))
  }
  
  if(mets=="RM.G") 
  {
    cor.within=matrix(999,reshuffle,n.envir);cor.all=numeric();
    for(i in 1:reshuffle)
    {
      cross=sample(rep(1:fold,each=ceiling(n.line/fold)),n.line);
      yhat.whole.cross=numeric();yobs.whole.cross=numeric();
      for(f in 1:fold)
      {
        id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
        ##marker effect###
        effect=matrix(999,n.marker,n.envir);
        intercept=matrix(999,1,n.envir)
        
        for(k in 1:n.envir)
        {
          fit=mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
          effect[,k]=fit$u
          intercept[,k]=fit$beta
        }
        
        ##Slope###
        effect.hat=numeric();
        for(j in 1:n.marker)
        {
          x1=envir[,enp];
          y1=effect[j,];
          coe=lm(y~x,data=data.frame(x=x1,y=y1));
          y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[,enp];
          effect.hat=rbind(effect.hat,y.hat);
        }
        
        ##Environment intercept####
        reg.intercept=data.frame(x=as.numeric(envir[,enp]),y=as.vector(intercept));
        coe.intercept=lm(y~x,data=reg.intercept)
        y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(envir[,enp]);
        
        ###All the predicted slope and intercept
        yhat.envir=matrix(999,length(id.V),n.envir);yobs.envir=matrix(999,length(id.V),n.envir)
        for(j in 1:n.envir)
        {
          yhat=y.intercept[j]+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat))[,j];
          yobs=pheno[id.V,j+1];
          yhat.envir[,j]=yhat;yobs.envir[,j]=yobs;
        }
        yhat.whole.cross=rbind(yhat.whole.cross,yhat.envir);
        yobs.whole.cross=rbind(yobs.whole.cross,yobs.envir);
      }
      
      for(j in 1:n.envir)
      {cor.within[i,j]=cor(yhat.whole.cross[,j],yobs.whole.cross[,j],use = "complete.obs");}
      cor.all=c(cor.all,cor(as.vector(yhat.whole.cross),as.vector(yobs.whole.cross),use = "complete.obs"));
    }
    #Correlation within environment 50 times
    r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
    #Correlation across environment 50 times
    r_across=mean(cor.all);
    #Observation and prediction last time
    outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                            pre=as.vector(yhat.whole.cross),
                            col=rep(coloo,times=rep(n.line,n.envir)),envir=rep(colnames(pheno)[-1],times=rep(n.line,n.envir)))
    colnames(cor.within)=colnames(pheno)[-1];
    r_within=cor.within;
    r_across=cor.all;
  }   
  
  if(mets=="RM.GE") 
  { 
    cor.within=matrix(999,reshuffle,n.envir);cor.all=numeric();
    for(i in 1:reshuffle)
    {
      cross=sample(rep(1:fold,each=ceiling(n.line/fold)),n.line);
      obs_matrix=numeric();pre_matrix=numeric();
      for(f in 1:fold)
      {
        id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
        ##marker effect###
        effect=matrix(999,n.marker,n.envir);
        intercept=matrix(999,1,n.envir)
        
        for(k in 1:n.envir)
        {
          fit=mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
          effect[,k]=fit$u
          intercept[,k]=fit$beta
        }
        
        ##predict marker effect###
        obs.envir=numeric();pre.envir=numeric();
        for(kk in 1:n.envir)
        {
          effect.hat=numeric();
          for(j in 1:n.marker)
          {
            x1=envir[-kk,enp];
            y1=effect[j,-kk];
            coe=lm(y~x,data=data.frame(x=x1,y=y1));
            y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[kk,enp];
            effect.hat=rbind(effect.hat,y.hat);
          }
          ##Environment intercept####
          reg.intercept=data.frame(x=as.numeric(envir[-kk,enp]),y=as.vector(intercept[-kk]));
          coe.intercept=lm(y~x,data=reg.intercept)
          y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(envir[kk,enp]);
          
          ###All the predicted slope and intercept
          yhat=y.intercept+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat));
          yobs=pheno[id.V,kk+1];
          obs.envir=cbind(obs.envir,yobs);
          pre.envir=cbind(pre.envir,yhat);
        }
        obs_matrix=rbind(obs_matrix,obs.envir);pre_matrix=rbind(pre_matrix,pre.envir);
        
      }
      cor.shuffle=cor(as.vector(obs_matrix),as.vector(pre_matrix),use = "complete.obs")
      cor.all=c(cor.all,cor.shuffle);
      for(j in 1:n.envir)
      {cor.within[i,j]=cor(obs_matrix[,j],pre_matrix[,j],use = "complete.obs");}
    } 
    
    yhat.whole.cross=pre_matrix;
    yobs.whole.cross=obs_matrix;
    #Correlation within environment 50 times
    r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
    #Correlation across environment 50 times
    r_across=mean(cor.all);
    #Observation and prediction last time
    outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                            pre=as.vector(yhat.whole.cross),
                            col=rep(coloo,times=rep(n.line,n.envir)),envir=rep(colnames(pheno)[-1],times=rep(n.line,n.envir)))
    colnames(cor.within)=colnames(pheno)[-1];
    r_within=cor.within;
    r_across=cor.all;
    
  }   
  
  return(list(outforfigure,r_within,r_across));
}
