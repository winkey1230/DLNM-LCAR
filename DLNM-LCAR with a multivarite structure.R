################################################################################
# DLNM-LCAR with a multivariate structure                                      #
# The following code is just an example focusing on several interesting DATs  # 
# INLA package may not directly implement the parameter estimation, so we use  #
# the 'jags' software, which may cost much computing resources                 #                                                     
################################################################################

################################################################################
## load function and exposure-response data from the first stage ###############
################################################################################
path <- "your path"
setwd(path)
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%`
multimerge <- function(dat=list(),...){
  if(length(dat)<2)return(as.data.frame(dat))
  mergedat<-dat[[1]]
  dat[[1]]<-NULL
  for(i in dat){
    mergedat<-merge(mergedat,i,...)
  }
  return(mergedat)
}
library(tmap)
library(RColorBrewer)
library(R2jags)
#### load data from the first stage, DLNM without time-varying
# cumfit: cummulated log-RR; cumcov: the covariance for  cummulated log-RR
# ismissing: True indicates the DAT is out of the range
# R: associated with the adjacent matrix, the details can be seen in LCAR model
load("Data for DLNM-LCAR with a multivariate structure.Rdata")

################################################################################
##Multivariate LCAR based on the second-hand datasets with estimation error ####
################################################################################
xvar <- c(-8,-3,2,7,16,20,24,28) 
## likelihood model -----
modelLL <- function(){
  for (i in 1:n) {
    y[i,] ~ dmnorm(theta[i,1:p],invsigma.obs[i,,])
    theta[i,1:p] <- eta + xi[i,]
  }
  xi <- A %*% x %*% t(B)
  # L is set as constant in this example, i.e., V is a structure correlation matrix
  # When a non-structure correlation matrix is needed, L is a lower triangular
  # matrix to be estimated, and for each non-zero element, a weak prior is set
  A <- L%*% S 
  for (i in 1:n) {
    S[i,1:n] <- 1/(tau - tau * rho * lambda)^0.5 * I[i,]
  }
  tau ~ dgamma(0.1,0.1)
  rho ~ dunif(lowrho,uprho)
  for (i in 1:p) {
    eta[i] ~ dnorm(0,0.0001)
    x[1:n,i] ~ dmnorm(rep(0,n),I)
  }
}

## run model -----
nstate <- ncol(ismissing)
invcumcov <- cumcov
for (i in 1:nstate) {
  invcumcov[i,,] <- chol2inv(chol(cumcov[i,,]))
}
y <- cumfit; invsigma.obs <- invcumcov
p <- ncol(y); n <- nrow(y);I <- diag(n); Cmatrix <- diag(n) - R
eig <- eigen(Cmatrix); L <- eig$vectors; lambda <- eig$values
uprho <- min(1/lambda[lambda>0]); lowrho <- max(1/lambda[lambda<0])
vstr <- apply(apply(cumcov, 1, cov2cor), 1, mean) %>% matrix(nrow = length(xvar))
eig <- eigen(vstr); B <- eig$vectors %*% diag(eig$values^0.5)
jags_data <- list("y","n","p","invsigma.obs","I","uprho","lowrho","lambda","L","B")
jags_params <- c("eta","theta","rho")
# system.time({
#   jagsfit <- jags.parallel(data=jags_data,parameters.to.save = jags_params,#inits = jags.inits, 
#                            n.iter=10000,model.file=modelLL,n.burnin = 10000,
#                            n.thin = 10,n.chains = 3)
# })
jagsfit <- jags(data=jags_data,parameters.to.save = jags_params,#inits = jags.inits,
                n.iter=200,model.file=modelLL,n.burnin = 100,
                n.thin = 1,n.chains = 2)

eta0.5 <- apply(jagsfit$BUGSoutput$sims.list$eta, 2, median)
eta0.975 <- apply(jagsfit$BUGSoutput$sims.list$eta, 2, quantile,probs = 0.975)
eta0.025 <- apply(jagsfit$BUGSoutput$sims.list$eta, 2, quantile,probs = 0.025)


theta0.5 <- sapply(1:length(xvar), function(x) apply(jagsfit$BUGSoutput$sims.list$theta[,,x], 2, median)) 
theta0.975 <- sapply(1:length(xvar), function(x) apply(jagsfit$BUGSoutput$sims.list$theta[,,x], 2, quantile,probs = 0.975)) 
theta0.025 <- sapply(1:length(xvar), function(x) apply(jagsfit$BUGSoutput$sims.list$theta[,,x], 2, quantile,probs = 0.025)) 
colnames(theta0.5) <- colnames(theta0.975) <- colnames(theta0.025) <- xvar

rho0.5 <- apply(jagsfit$BUGSoutput$sims.list$rho, 2, quantile,probs = c(0.025,0.5,0.975))
traceplot(jagsfit)
mean((theta_fit0.5-theta)^2)

