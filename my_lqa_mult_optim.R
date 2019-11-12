source("functions.R")


set.seed(2)  #1 good
neach=10
means=expand.grid(rep=1:neach,mu1=c(2,15)*5,mu2=c(2,15))[,-1]
dat=rpois(prod(dim(means)),lambda = as.matrix(means))
dat_mat=matrix(dat,ncol=2)
lambdas=seq(0,3,length.out = 30)

beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=FALSE)
plot_path_2d_alt(dat_mat,beta_out,title="Unweighted poisson")

lambdas=lambdas*5
beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=TRUE)
plot_path_2d_alt(dat_mat,beta_out,title="Weighted poisson")


J=matrix(runif(prod(dim(dat_mat))),nrow = nrow(dat_mat))/5
Y=log((dat_mat+J))
lambdas=c(0,seq(0.001,0.07,length.out = 100))

betasHat=matrix(NA,prod(dim(Y)),length(lambdas))
for(i in 1:length(lambdas))
{
  betasHat[,i] <- as.vector(t(PRclust(t(Y),lambda1=0.4 ,lambda2=lambdas[i],tau=1e10,loss.method="lasso",algorithm = "ADMM")$mu))
}
plot_path_2d_alt(Y,betasHat,title="PRclust default log of data")





load(file="poisson_bivariate_4_grps")
lambdas=c(0,seq(0.001,2,length.out = 100))

dat_mat=dat$Y
beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=FALSE)
plot_path_2d_alt(dat_mat,beta_out)

beta_out<-clust_path(dat_mat,lambdas,tau=Inf,weights=TRUE)
plot_path_2d_alt(dat_mat,beta_out)


J=matrix(runif(prod(dim(dat$Y))),nrow = nrow(dat$Y))
Y=log((dat$Y++J))
lambdas=c(0,seq(0.001,0.07,length.out = 100))

betasHat=matrix(NA,prod(dim(Y)),length(lambdas))
for(i in 1:length(lambdas))
{
  betasHat[,i] <- as.vector(t(PRclust(t(Y),lambda1=0.4 ,lambda2=lambdas[i],tau=1e10,loss.method="lasso",algorithm = "ADMM")$mu))
}

plot_path_2d_alt(Y,betasHat)



## to do

## more sparcity
## speed
## higher dimension
## compare to PRclust




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


