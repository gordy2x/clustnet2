rm(list=ls())
library(glmmTMB)
library(tidyr)
library(mvabund)
library(ggplot2)
data("spider")

a_lam<-function(beta,lambda,tau=1){
  n_all=length(beta)
  cols=ncol(beta)
  rows=nrow(beta)
  betalong=matrix(as.vector(beta),ncol=1)
  blanki=matrix(0,nrow=rows,ncol=1)
  out=matrix(0,n_all,n_all)
  c=1e-8
  #eq 20 in petry 2011
  for(i in 1:(rows-1)){
    for(k in (i+1):rows){
      aj=blanki
      aj[i,1]=1
      aj[k,1]=-1
      ajmat=kronecker(diag(cols),aj)
      ajbeta=abs(t(ajmat)%*%betalong) #i maj be able to just do aj%*%beta
      if(sqrt(sum(ajbeta^2))<tau){
        out=out+1/sqrt(sum(ajbeta^2)+c)*(ajmat%*%t(ajmat))
      }
    }
  }
  lambda*out
}



a_lam_w<-function(mbeta,lambda,w){
  n_all=length(mbeta)
  cols=ncol(mbeta)
  rows=nrow(mbeta)
  betalong=matrix(as.vector(mbeta),ncol=1)
  blanki=matrix(0,nrow=rows,ncol=1)
  out=matrix(0,n_all,n_all)
  c=1e-8
  #eq 20 in petry 2011
  for(i in 1:(rows-1)){
    for(k in (i+1):rows){
      aj=blanki
      aj[i,1]=w[i,k]
      aj[k,1]=-w[i,k]
      ajmat=kronecker(diag(cols),aj)
      ajbeta=abs(t(ajmat)%*%betalong) #i maj be able to just do aj%*%beta
      out=out+1/sqrt(sum(ajbeta^2)+c)*(ajmat%*%t(ajmat))
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


clust_lasso<-function(dat_mat,lambda,tau,beta_init=NULL,delta=1){
  print(paste("lambda=",lambda))
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
  wi=dist(log(dat_mat+0.001))^2
  w=1/wi
  w[wi==0]=0
  w=w/sum(w)*length(w)
  w=as.matrix(w)
  
  
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
      print(same)
      beta=equalise(beta,same,rows)
    }

    
    score=obj$gr(beta)
    # print(score)
    hess=obj$he(beta)
  
    alam=a_lam_w(matrix(beta,nrow=rows),lambda,w=w)
    alam2=a_lam(matrix(beta,nrow=rows),lambda,tau=0.5)
    # image.plot(alam)
    # image.plot(alam2)
    mbeta=matrix(beta,ncol=1)
    hess=try(solve(hess+alam),silent = TRUE)
    if(class(hess)!="try-error"){
      beta=beta-delta*hess%*%(t(score)+alam%*%mbeta)
    } else {
      beta=beta-delta*(t(score)+alam%*%mbeta)
    }
    diff=max(abs(betaold-beta)/abs(beta))
    # print(diff)
    if(diff>diffold){
      delta=delta/2
    }
  }
  beta
}


clust_path<-function(dat_mat,lambdas,tau){
  nlam=length(lambdas)
  betas=matrix(NA,prod(dim(dat_mat)),nlam)
  betas[,1]=clust_lasso(dat_mat,lambdas[1],tau)
  for(i in 2:nlam){
    betas[,i]=clust_lasso(dat_mat,lambdas[i],beta_init=betas[,(i-1)],tau=tau)
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





plot_path_2d_alt<-function(dat_mat,beta_out){
  nsites=nrow(dat_mat)
  nspecies=ncol(dat_mat)
  beta_plot=as.data.frame(beta_out)
  beta_plot$species=rep(1:nspecies,each=nsites)
  beta_plot$sites=factor(rep(1:nsites,nspecies))
  beta_plot=gather(beta_plot,lambda,beta,1:length(lambdas),factor_key=TRUE)
  beta_plot$lambda=factor(beta_plot$lambda)
  beta_plot=spread(beta_plot,species,beta,drop=FALSE)
  # beta_plot$lambda=rep(lambdas,nsites)
  colnames(beta_plot)[-(1:2)]=c("sp1","sp2")
  ggplot(beta_plot,aes(x=exp(sp1),y=exp(sp2)))+
    geom_point(aes(color=lambda,group=lambda))+geom_path(aes(group=sites))+
    theme(legend.position = "none")
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


load(file="poisson_bivariate_4_grps")
lambdas=seq(0,1,length.out = 10)

tau=NA  #0.8 is interesting. 5 seems wrong
dat_mat=dat$Y
beta_out<-clust_path(dat_mat,lambdas)
plot_path_2d_alt(dat_mat,beta_out)


set.seed(1)
neach=10
means=expand.grid(rep=1:neach,mu1=c(2,15),mu2=c(2,15))[,-1]
dat=rpois(prod(dim(means)),lambda = as.matrix(means))
dat_mat=matrix(dat,ncol=2)
lambdas=lambdas=seq(0,1,length.out = 10)
tau=NA

beta_out<-clust_path(dat_mat,lambdas,tau)

plot_path_2d_alt(dat_mat,beta_out)







## to do

## more sparcity
## speed
## higher dimension





#1d sim

set.seed(1)
neach=10
means=expand.grid(rep=1:neach,mu1=c(2,15),mu2=c(2,15))[,-1]
dat=rpois(prod(dim(means)),lambda = as.matrix(means))
dat_mat=matrix(dat,ncol=2)
dat_mat=as.matrix(dat_mat[,1])
lambdas=seq(0,1,0.1)
tau=5

beta_out<-clust_path(dat_mat,lambdas,tau)
plot_path_1d(dat_mat,beta_out)



# 
# 
# pen_gr<-function(beta,obj,lambda){
#   alam=a_lam(matrix(beta,nrow=rows),lambda,tau=tau)
#   obj$gr(beta)+matrix(beta,nrow=1)%*%alam
# }
# 
# 
# 
# pen_he<-function(beta,obj,lambda){
#   alam=a_lam(matrix(beta,nrow=rows),lambda,tau=tau)
#   obj$he(beta)+alam
# }
