source("functions.R")

# pdf("poisson_dat_compare.pdf")
set.seed(1)  #1 good
base=c(2,15)
beta_in=log(expand.grid(base,2*base[1],3*base,4*base[1],5*base[2]))
neach=rep(10,nrow(beta_in))
dat=sim_pois(beta_in,neach)
dat_mat=dat$Y
rownames(dat_mat)=
rows=nrow(dat_mat)
# beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=FALSE,poisson = TRUE)
# plot_path_2d_alt(dat_mat,beta_out,title="Unweighted poisson")
# 
# 
# beta_out<-clust_path(dat_mat,lambdas/20,tau=Inf,weights=FALSE,poisson = FALSE)
# plot_path_2d_alt(dat_mat,beta_out,title="Unweighted gaussian log")
lambdas=seq(0.01,3,length.out = 150)


beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=TRUE,poisson = TRUE)
plot_path_lat(dat_mat,beta_out,title="5 species 4 clusters")

form_clust=create_merge(beta_out,dat_mat,lambdas)
dend=as.dendrogram(form_clust)
plot(dend)



# plot_path_2d_alt(dat_mat,beta_out,title="Weighted poisson")
# mse=apply(beta_out,2,function(x) sum((matrix(x,ncol=5)-log(means))^2))
# right_cluster=apply(beta_out,2,same_clust,cluster=rep(1:4,each=10))
# lik=apply(beta_out,2,mix_llik,dat_mat=dat_mat)
# em=apply(beta_out,2,EM,dat_mat=dat_mat)
clusters=apply(beta_out,2,function(x) unique(matrix(x,nrow=rows)))
nclust=unlist(lapply(clusters, function(x) nrow(x)))
# bic_pen=log(rows)*(nclust*cols)
# bic=-2*lik+bic_pen
# icl=bic-em


# plot_path_2d_alt(dat_mat,beta_out[,1:which(bic==min(bic))[1]],
                 # title="Weighted poisson bic")
# plot_path_2d_alt(dat_mat,beta_out[,1:which(icl==min(icl))[1]],
#                  title="Weighted poisson ICL")
# plot_path_2d_alt(dat_mat,beta_out[,1:which(right_cluster==min(right_cluster))[1]],
#                  title="Weighted poisson clust")

# print(nclust[bic==min(bic)])
# print(nclust[icl==min(icl)])
# 
# betamat=matrix(beta_out[,which(icl==min(icl))[1]],nrow=rows)
# image.plot(betamat)


plot(nclust,bic,type="l",ylim=c(0,1000))
lines(nclust,-2*lik,col="blue")
lines(nclust,bic_pen,col="red")
lines(nclust,mse*10,col="orange")
lines(nclust,right_cluster*5,col="purple")



# plot(lambdas,mse,type="l")
# plot(lambdas,right_cluster,type="l")
# plot(lambdas,-2*lik,type="l")
# plot(lambdas,bic,type="l")
# plot(lambdas,icl,type="l")
# # 
# plot_path_2d_alt(dat_mat,beta_out[,1:which(mse==min(mse))],title="Weighted poisson mse")
# plot_path_2d_alt(dat_mat,beta_out[,1:which(right_cluster==min(right_cluster))[1]],title="Weighted poisson min good clust")




## write function to plug in data with cluster centres
## need to do sparsity for gaussian

# 
# beta_out<-clust_path(dat_mat,lambdas/2,tau=Inf,weights=TRUE,poisson = FALSE)
# plot_path_2d_alt(dat_mat,beta_out,title="Weighted gaussian log")
# 
# 
# 
# J=matrix(runif(prod(dim(dat_mat))),nrow = nrow(dat_mat))/5
# Y=log((dat_mat+J))
# lambdas=c(0,seq(0.001,0.035,length.out = 100))
# 
# betasHat=matrix(NA,prod(dim(Y)),length(lambdas))
# for(i in 1:length(lambdas))
# {
#   betasHat[,i] <- as.vector(t(PRclust(t(Y),lambda1=0.4 ,lambda2=lambdas[i],tau=1e10,loss.method="lasso",
#                                       grouping.penalty = "gtlp",algorithm = "ADMM")$mu))
# }
# plot_path_2d_alt(Y,betasHat,title="PRclust gaussian log")
# 
# 
# dev.off()
# 
# 
# load(file="poisson_bivariate_4_grps")
# lambdas=c(0,seq(0.001,2,length.out = 100))
# 
# dat_mat=dat$Y
# beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=FALSE)
# plot_path_2d_alt(dat_mat,beta_out)
# 
# beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=TRUE)
# plot_path_2d_alt(dat_mat,beta_out)
# 
# 
# J=matrix(runif(prod(dim(dat$Y))),nrow = nrow(dat$Y))
# Y=log((dat$Y++J))
# lambdas=c(0,seq(0.001,0.07,length.out = 100))
# 
# betasHat=matrix(NA,prod(dim(Y)),length(lambdas))
# for(i in 1:length(lambdas))
# {
#   betasHat[,i] <- as.vector(t(PRclust(t(Y),lambda1=0.4 ,lambda2=lambdas[i],tau=1e10,loss.method="lasso",algorithm = "ADMM")$mu))
# }
# 
# plot_path_2d_alt(Y,betasHat)
# 
# 
# 
# ## to do
# 
# ## more sparcity
# ## speed
# ## higher dimension
# ## compare to PRclust
# 
# 
# 
# 
# #1d sim
# 
# set.seed(1)
# neach=10
# means=expand.grid(rep=1:neach,mu1=c(2,15),mu2=c(2,15))[,-1]
# dat=rpois(prod(dim(means)),lambda = as.matrix(means))
# dat_mat=matrix(dat,ncol=2)
# dat_mat=as.matrix(dat_mat[,1])
# lambdas=seq(0,1,0.1)
# tau=5
# 
# beta_out<-clust_path(dat_mat,lambdas,tau)
# plot_path_1d(dat_mat,beta_out)
# 
# 
# 
# # 
# # 
# # pen_gr<-function(beta,obj,lambda){
# #   alam=a_lam(matrix(beta,nrow=rows),lambda,tau=tau)
# #   obj$gr(beta)+matrix(beta,nrow=1)%*%alam
# # }
# # 
# # 
# # 
# # pen_he<-function(beta,obj,lambda){
# #   alam=a_lam(matrix(beta,nrow=rows),lambda,tau=tau)
# #   obj$he(beta)+alam
# # }
# 
# 




set.seed(1)  #1 good
base=c(2,15)
beta_in=log(expand.grid(base,2*base[1],3*base,4*base[1],5*base[2]))
neach=rep(2,nrow(beta_in))
dat=sim_pois(beta_in,neach)
dat_mat=dat$Y
rows=nrow(dat_mat)
# beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=FALSE,poisson = TRUE)
# plot_path_2d_alt(dat_mat,beta_out,title="Unweighted poisson")
# 
# 
# beta_out<-clust_path(dat_mat,lambdas/20,tau=Inf,weights=FALSE,poisson = FALSE)
# plot_path_2d_alt(dat_mat,beta_out,title="Unweighted gaussian log")
lambdas=seq(0.01,9,length.out = 100)


beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=TRUE,poisson = TRUE)
plot_path_lat(dat_mat,beta_out,title="5 species 4 clusters")
