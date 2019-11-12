# rm(list=ls())
library(glmmTMB)
library(tidyr)
library(mvabund)
library(ggplot2)
library(prclust)
library(fields)
data("spider")




a_lam<-function(mbeta,lambda,w=NULL,tau=Inf){
  n_all=length(mbeta)
  cols=ncol(mbeta)
  rows=nrow(mbeta)
  betalong=matrix(as.vector(mbeta),ncol=1)
  blanki=matrix(0,nrow=rows,ncol=1)
  out=matrix(0,n_all,n_all)
  c=1e-8
  if(is.null(w)){
    w=matrix(1,length(betalong),length(betalong))
  }
  #eq 20 in petry 2011
  for(i in 1:(rows-1)){
    for(k in (i+1):rows){
      aj=blanki
      aj[i,1]=w[i,k]
      aj[k,1]=-w[i,k]
      ajmat=kronecker(diag(cols),aj)
      ajbeta=abs(t(ajmat)%*%betalong) #i maj be able to just do aj%*%beta
      if(sqrt(sum(ajbeta^2))<tau){
        out=out+1/sqrt(sum(ajbeta^2)+c)*(ajmat%*%t(ajmat))
      }
    }
  }
  lambda*out
}







equalise<-function(beta,same,rows){
  mbeta=matrix(beta,nrow=rows)
  for(i in 1:nrow(same)){
    mbeta[same[i,1],]=mbeta[same[i,2],]
  }
  c(mbeta)
}


clust_lasso<-function(dat_mat,lambda,tau=Inf,beta_init=NULL,delta=1, weights=TRUE){
  # print(paste("lambda=",lambda))
  cols=ncol(dat_mat)
  rows=nrow(dat_mat)
  y=as.vector(as.matrix(dat_mat))
  x=factor(1:length(y))
  TMBStruc=glmmTMB(y~-1+x,family=poisson,doFit = FALSE)
  obj <- with(TMBStruc, TMB::MakeADFun(data.tmb, parameters, 
                                       map = mapArg, random = randomArg, profile = NULL, 
                                       silent = !verbose, DLL = "glmmTMB"))
  
  
  if(is.null(beta_init)){
    beta_init=log(y+0.001)  
  }
  if(weights){
    wi=dist(log(dat_mat+0.001))
    w=1/wi
    w[wi==0]=0
    w=w/sum(w)*length(w)
    w=as.matrix(w)
  }else(
    w=matrix(1,length(y),length(y))
  )

  
  
  N_beta=length(beta)
  eps=1e-3
  diff=eps+1
  beta=beta_init
  same_track=NULL
  pars=obj$par
  while(diff>eps){
    diffold=diff
    betaold=beta
    

    # }
    
    
    
    
    pars=c(beta)
    score=obj$gr(pars)
    hess=obj$he(pars)
    alam=a_lam(matrix(beta,nrow=rows),lambda,w=w,tau=tau)
    mbeta=matrix(beta,ncol=1)
    hess_inv=try(solve(hess+alam),silent = TRUE)
    if(class(hess_inv)!="try-error"){
      beta=beta-delta*hess_inv%*%(t(score)+alam%*%mbeta)
    } else {
      beta=beta-delta*1/mean(diag(hess+alam))(t(score)+alam%*%mbeta)
    }
    
    
    blank=matrix(NA,rows,rows)
    blank[lower.tri(blank)]=apply(apply(matrix(beta,nrow=rows),2,function(v) c(dist(v))),1,function(x) sqrt(mean(x^2)))  
    # print(sum(blank==0,na.rm=T))
    # blank[same_track]=NA
    # image.plot(blank)
    same=which(blank<0.01,arr.ind = TRUE)
    # same_track=unique(rbind(same_track,same))
    beta=equalise(beta,same,rows)
    # blank[lower.tri(blank)]=apply(apply(matrix(beta,nrow=rows),2,function(v) c(dist(v))),1,function(x) sqrt(mean(x^2))) 
    # print(sum(blank==0,na.rm=T))
    
    diff=max(abs(betaold-beta)/abs(beta))
    # print(diff)
    if(diff>diffold){
      delta=delta*0.9
      beta=betaold
    }
  }
  beta

}





clust_lasso_gauss<-function(dat_mat,lambda,tau=Inf,beta_init=NULL,delta=1, weights=TRUE){
  # print(paste("lambda=",lambda))
  cols=ncol(dat_mat)
  rows=nrow(dat_mat)
  y=log(as.vector(as.matrix(dat_mat))+0.001)
  x=factor(1:length(y))
  TMBStruc=glmmTMB(y~-1+x,doFit = FALSE)
  obj <- with(TMBStruc, TMB::MakeADFun(data.tmb, parameters, 
                                       map = mapArg, random = randomArg, profile = NULL, 
                                       silent = !verbose, DLL = "glmmTMB"))
  
  if(is.null(beta_init)){
    beta_init=y
  }
  if(weights){
    wi=dist(log(dat_mat+0.001))^2
    w=1/wi
    w[wi==0]=0
    w=w/sum(w)*length(w)
    w=as.matrix(w)
  }else(
    w=matrix(1,length(y),length(y))
  )
  
  
  
  N_beta=length(beta)
  eps=1e-3
  diff=eps+1
  beta=beta_init
  same_track=NULL
  
  while(diff>eps){
    diffold=diff
    betaold=beta

    blank=matrix(NA,rows,rows)
    blank[lower.tri(blank)]=apply(apply(matrix(beta,nrow=rows),2,function(v) c(dist(v))),1,function(x) mean(abs(x)))
    blank[same_track]=NA
    same=which(blank<0.01,arr.ind = TRUE)
    same_track=unique(rbind(same_track,same))
    if(nrow(same)>0){
      # print(same)
      beta=equalise(beta,same,rows)
    }

    
    score=obj$gr(c(beta,1))
    last=length(score)
    score=score[-last]
    he=obj$he(c(beta,1))
    he=he[-last,-last]
    alam=a_lam(matrix(beta,nrow=rows),lambda,w=w,tau=tau)
    mbeta=matrix(beta,ncol=1)
    hess=try(solve(he+alam),silent = TRUE)
    if(class(hess)!="try-error"){
      beta=beta-delta*hess%*%(t(t(score))+alam%*%mbeta)
    } else {
      beta=beta-delta*(t(t(score))+alam%*%mbeta)
    }

    
    
    diff=max(abs(betaold-beta)/abs(beta))
    # print(diff)
    if(diff>diffold){
      delta=delta/2
    }
  }
  beta
}


clust_path<-function(dat_mat,lambdas,tau=Inf, weights=TRUE,poisson=TRUE){
  nlam=length(lambdas)
  betas=matrix(NA,prod(dim(dat_mat)),nlam)
  fn=clust_lasso
  if(!poisson){fn=clust_lasso_gauss}
  betas[,1]=fn(dat_mat,lambdas[1],tau)
  for(i in 2:nlam){
    betas[,i]=fn(dat_mat,lambdas[i],beta_init=betas[,(i-1)],tau=tau,weights=weights)
  }
  betas
}



plot_path_2d<-function(dat_mat,beta_out){
  nsites=nrow(dat_mat)
  nspecies=ncol(dat_mat)
  beta_plot=as.data.frame(beta_out)
  beta_plot$species=rep(1:nspecies,each=nsites)
  beta_plot$sites=factor(rep(1:nsites,nspecies))
  
  beta_plot=gather(beta_plot,lambda,beta,1:length(lambdas),factor_key=TRUE)
  beta_plot=spread(beta_plot,species,beta,drop=FALSE)
  beta_plot$lambda=rep(lambdas,nsites)
  colnames(beta_plot)[-(1:2)]=c("sp1","sp2")
  ggplot(beta_plot,aes(x=exp(sp1),y=exp(sp2),color=sites,group=sites))+geom_path()+geom_point()+
    theme(legend.position = "none")
}





plot_path_2d_alt<-function(dat_mat,beta_out,title=""){
  nsites=nrow(dat_mat)
  nspecies=ncol(dat_mat)
  y=c(dat_mat)
  y_beta_out=cbind(log(y+0.01),beta_out)
  beta_plot=as.data.frame(y_beta_out)
  beta_plot$species=rep(1:nspecies,each=nsites)
  beta_plot$sites=factor(rep(1:nsites,nspecies))
  beta_plot=gather(beta_plot,lambda,beta,1:ncol(y_beta_out),factor_key=TRUE)
  beta_plot$lambda=factor(beta_plot$lambda)
  beta_plot=spread(beta_plot,species,beta,drop=FALSE)
  # beta_plot$lambda=rep(lambdas,nsites)
  colnames(beta_plot)[-(1:2)]=c("sp1","sp2")
  ggplot(beta_plot,aes(x=exp(sp1),y=exp(sp2)))+
    geom_point(aes(color=lambda,group=lambda))+geom_path(aes(group=sites))+
    theme(legend.position = "none")+ggtitle(title)
}



plot_path_1d<-function(dat_mat,beta_out){
  nsites=nrow(dat_mat)
  nspecies=ncol(dat_mat)
  beta_plot=as.data.frame(beta_out)
  beta_plot$Y=dat_mat
  beta_plot$species=rep(1:nspecies,each=nsites)
  beta_plot$sites=factor(rep(1:nsites,nspecies))
  beta_plot=gather(beta_plot,lambda,beta,1:length(lambdas),factor_key=TRUE)
  beta_plot$lambda=factor(beta_plot$lambda)
  ggplot(beta_plot,aes(x=Y,y=exp(beta)))+
    geom_point(aes(color=lambda,group=lambda))+geom_path(aes(group=sites))+
    theme(legend.position = "none")
}



same_clust<-function(betavec,cluster){
  out=0
  rows=length(cluster)
  betamat=matrix(betavec,nrow=rows)
  for(i in 1:(rows-1)){
    for(j in i:rows){
      betadiff=dist(betamat[i,]-betamat[j,])
      if((cluster[i]==cluster[j])&(betadiff>0)){
        out=out+1
      }
      if((cluster[i]!=cluster[j])&(betadiff==0)){
        out=out+1
      }
        
    }
  }
  out
  
}

dmpois<-function(y,betamat){
  p=length(y)
  out=1
  for(i in 1:p){
    out=out*dpois(y[i],lambda = exp(betamat[,i]))
  }
  out
}


mix_llik<-function(betavec,dat_mat){
  ## work out pi from the dpois()
  
  ## for each y they sum to 1, then sum across cluster
  out=0
  rows=nrow(dat_mat)
  betamat=matrix(betavec,nrow=rows)
  clustermat=unique(betamat)
  
  
  
  
  for(i in 1:rows){
    dim1=dpois(dat_mat[i,1],lambda = exp(betamat[,1]))
    dim2=dpois(dat_mat[i,2],lambda = exp(betamat[,2]))
    
    out=out+log(sum(dim1*dim2/rows))
  }
  out
}
