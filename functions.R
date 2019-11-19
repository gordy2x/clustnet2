# rm(list=ls())
library(glmmTMB)
library(tidyr)
library(mvabund)
library(ggplot2)
library(prclust)
library(fields)
library(psych)
library(Hmisc)
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
  if(length(same)>0){
    for(i in 1:nrow(same)){
      mbeta[same[i,1],]=mbeta[same[i,2],]
    }
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
    beta_init=log(y+1)  
  }
  if(weights){
    wi=dist(log(dat_mat+1))
    w=1/wi
    w[wi==0]=0
    w=w/sum(w)*length(w)
    w=as.matrix(w)
  }else(
    w=matrix(1,length(y),length(y))
  )

  
  
  N_beta=length(beta)
  eps=1e-2
  diff=eps+1
  beta=beta_init
  same_track=NULL
  pars=obj$par
  while(diff>eps){
    diffold=diff
    betaold=beta
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
    same=which(blank<0.1,arr.ind = TRUE)
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




# 
# clust_lasso_gauss<-function(dat_mat,lambda,tau=Inf,beta_init=NULL,delta=1, weights=TRUE){
#   # print(paste("lambda=",lambda))
#   cols=ncol(dat_mat)
#   rows=nrow(dat_mat)
#   y=log(as.vector(as.matrix(dat_mat))+0.001)
#   x=factor(1:length(y))
#   TMBStruc=glmmTMB(y~-1+x,doFit = FALSE)
#   obj <- with(TMBStruc, TMB::MakeADFun(data.tmb, parameters, 
#                                        map = mapArg, random = randomArg, profile = NULL, 
#                                        silent = !verbose, DLL = "glmmTMB"))
#   
#   if(is.null(beta_init)){
#     beta_init=y
#   }
#   if(weights){
#     wi=dist(log(dat_mat+0.001))^2
#     w=1/wi
#     w[wi==0]=0
#     w=w/sum(w)*length(w)
#     w=as.matrix(w)
#   }else(
#     w=matrix(1,length(y),length(y))
#   )
#   
#   
#   
#   N_beta=length(beta)
#   eps=1e-3
#   diff=eps+1
#   beta=beta_init
#   same_track=NULL
#   
#   while(diff>eps){
#     diffold=diff
#     betaold=beta
# 
#     blank=matrix(NA,rows,rows)
#     blank[lower.tri(blank)]=apply(apply(matrix(beta,nrow=rows),2,function(v) c(dist(v))),1,function(x) mean(abs(x)))
#     blank[same_track]=NA
#     same=which(blank<0.05,arr.ind = TRUE)
#     same_track=unique(rbind(same_track,same))
#     if(nrow(same)>0){
#       # print(same)
#       beta=equalise(beta,same,rows)
#     }
# 
#     
#     score=obj$gr(c(beta,1))
#     last=length(score)
#     score=score[-last]
#     he=obj$he(c(beta,1))
#     he=he[-last,-last]
#     alam=a_lam(matrix(beta,nrow=rows),lambda,w=w,tau=tau)
#     mbeta=matrix(beta,ncol=1)
#     hess=try(solve(he+alam),silent = TRUE)
#     if(class(hess)!="try-error"){
#       beta=beta-delta*hess%*%(t(t(score))+alam%*%mbeta)
#     } else {
#       beta=beta-delta*(t(t(score))+alam%*%mbeta)
#     }
# 
#     
#     
#     diff=max(abs(betaold-beta)/abs(beta))
#     # print(diff)
#     if(diff>diffold){
#       delta=delta/2
#     }
#   }
#   beta
# }


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




plot_path_lat<-function(dat_mat,beta_out,title=""){
  nsites=nrow(dat_mat)
  nspecies=ncol(dat_mat)
  nlam=ncol(beta_out)+1
  #add original Y to plot
  y_beta_out=as.data.frame(cbind(log(c(dat_mat)+0.01),beta_out))
  
  
  
  if(nspecies>2){
    lY=scale(log(dat_mat+0.01),scale=FALSE)
    dat_fa=fa(lY,nfactors=2,rotate="none",scores="Thurstone",fm="ml")
    load=as.matrix(dat_fa$loadings)
    sY=t(scale(lY))
    # dat_fa$scores
    score_loc=t(solve(diag(2)+t(load)%*%load)%*%t(load)%*%sY)
    # plot(dat_fa$scores,t(score_loc))
    
    
    
    load_mat=matrix(NA,nsites*2,nlam)
    load_mat[,1]=c(score_loc)
    for(i in 1:nlam){
      beta_loc=matrix(y_beta_out[,i],ncol=nspecies)
      sY=scale(t(beta_loc),scale = FALSE)
      load_loc=solve(diag(2)+t(load)%*%load)%*%t(load)%*%sY

      load_mat[,i]=c(t(load_loc))
    }
    
    beta_plot=as.data.frame(load_mat)
    beta_plot$species=rep(1:2,each=nsites)
    beta_plot$sites=factor(rep(1:nsites,2))
    beta_plot=gather(beta_plot,lambda,beta,1:ncol(load_mat),factor_key=TRUE)
    beta_plot$lambda=factor(beta_plot$lambda)
    beta_plot=spread(beta_plot,species,beta,drop=FALSE)
    colnames(beta_plot)[-(1:2)]=c("factor1","factor2")
    ggplot(beta_plot,aes(x=exp(factor1),y=exp(factor2)))+
      geom_point(aes(color=lambda,group=lambda))+geom_path(aes(group=sites))+
      theme(legend.position = "none")+ggtitle(title)
    

  }else{
    beta_plot=y_beta_out
    beta_plot$species=rep(1:2,each=nsites)
    beta_plot$sites=factor(rep(1:nsites,nspecies))
    beta_plot=gather(beta_plot,lambda,beta,1:ncol(y_beta_out),factor_key=TRUE)
    beta_plot$lambda=factor(beta_plot$lambda)
    beta_plot=spread(beta_plot,species,beta,drop=FALSE)
    colnames(beta_plot)[-(1:2)]=c("sp1","sp2")
    ggplot(beta_plot,aes(x=exp(sp1),y=exp(sp2)))+
      geom_point(aes(color=lambda,group=lambda))+geom_path(aes(group=sites))+
      theme(legend.position = "none")+ggtitle(title)
  }
  
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
  out=0
  for(i in 1:p){
    out=out+dpois(y[i],lambda = exp(betamat[,i]),log = TRUE)
  }
  exp(out)
}


z_opt<-function(betavec,dat_mat){

  out=0
  rows=nrow(dat_mat)
  betamat=matrix(betavec,nrow=rows)
  clustermat=unique(betamat)
  
  #should there be a loop here? 
  #z depends on tau (2.9/2.22) which depends on z(2.10)
  w=apply(dat_mat,1,dmpois,betamat=clustermat)
  z=scale(w, center = FALSE, scale = colSums(w))#assumes tau is 1
  tau=apply(z,1,sum)
  
  #is the even right?
  eps=1e-5
  delta=eps+1
  while(delta>eps){
    tauold=tau
    znum=t(scale(t(w),center = FALSE,scale=1/tau))
    z=scale(znum, center = FALSE, scale = colSums(znum))
    tau=apply(z,1,sum)
    delta=sum((tau-tauold)^2)
  }
  z
}

EM<-function(betavec,dat_mat){
  z=z_opt(betavec,dat_mat)
  zlz=z*log(z)
  -sum(zlz)
}


  
mix_llik<-function(betavec,dat_mat){
  rows=nrow(dat_mat)
  betamat=matrix(betavec,nrow=rows)
  clustermat=unique(betamat)
  w=apply(dat_mat,1,dmpois,betamat=clustermat)
  z=z_opt(betavec,dat_mat)
  tau=apply(z,1,sum)
  znum=t(scale(t(w),center = FALSE,scale=1/tau))
  
  sum(log(apply(znum,2,sum)))
  
}


sim_pois<-function(beta,neach){
  nclust=length(neach)
  index=rep(1:nclust,times=neach)
  lambda=exp(beta[index,])
  dat=rpois(prod(dim(lambda)),lambda = as.matrix(lambda))
  dat_mat=matrix(dat,ncol=ncol(beta))
  out=list(Y=dat_mat,beta=log(lambda))
  out
}


all_dim_equal<-function(x1,x2,mat=dat_mat){
  # x1==x2
  all(mat[x1,]==mat[x2,])
}

find_clust<-function(beta,dat_mat){
  
  cols=ncol(dat_mat)
  betamat=matrix(beta,ncol=cols)
  find.matches(betamat, betamat,tol=rep(0,cols), maxmatch=1)$matches
}



create_merge<-function(beta_out,dat_mat,lambdas){
  
  lambdas=c(0,lambdas)
  out=list()
  clusts=cbind(1:nrow(dat_mat),apply(beta_out,2,find_clust,dat_mat=dat_mat))
  clusts=as.data.frame(clusts)
  cols=ncol(clusts)
  rows=nrow(clusts)
  clusts_ord=arrange_all(clusts[,rev(1:cols)])
  out$order=as.integer(clusts_ord[,cols])
  nclust=apply(clusts,2,function(x) length(unique(x)))
  eps=1e-7
  merge=NULL
  height=NULL
  ct=1
  clust_memb=matrix(NA,rows,nrow(dat_mat))
  clust_memb[,1]=-clusts[,1]
  lamatmetge=rep(NA,rows-1)
  for(i in 2:cols){
    diff=clusts[,i]!=clusts[,i-1]
    if(any(diff)){
      w=which(diff)
      for(j in 1:length(w)){
          in_clust=which(clusts[,i]==clusts[w[j],i]) #find all other members of the cluster
          tomerge=unique(clust_memb[in_clust,ct])
          if(length(tomerge)>1){
            merge=rbind(merge,as.integer(tomerge[1:2]))
            clust_memb[,(ct+1)]=clust_memb[,ct]
            clust_memb[clust_memb[,ct]%in%tomerge[1:2],(ct+1)]=ct
            lamatmetge[ct]=lambdas[i]
            ct=ct+1
            if(j==1){
              height=c(height,lambdas[i]-lambdas[i-1])
            }else{
              height=c(height,eps)
            }
          }
      }
    }

    # print(merge)
    # print(c(i,w))
  }
  out$merge=merge
  out$height=lamatmetge
  # out$lam=lamatmetge
  class(out)="hclust"
  out
}


# 
# dat=gather(data.frame(t(clusts_ord[,rev(1:cols)])))
# 
# dat$lambda=c(-5,log(lambdas))
# dat%>%
#   ggplot(aes(lambda,value,color=key))+
#   geom_line()+theme_classic()
# 
# for(i in 1:40){
#   print(sum(clust_memb[,i]!=clust_memb[,i-1]))
# }

