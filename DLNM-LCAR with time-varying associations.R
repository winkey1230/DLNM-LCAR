################################################################################
## Use the temperamental LCAR to estimate the modification effects             #
## The time-varying DLNM-LCAR is similar to DLNM-LCAR,so the code is not given #
## in this file                                                                #
## Author: Wang, Wei   2022-07-15                                              #
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
library(INLA)
##############################################################################-#
# The following date is of interest due to the representativeness             -#
# c("2021-01-01","2021-02-01","2021-03-01","2021-04-01","2021-05-01",         -#
#   "2021-06-01","2021-07-01","2021-08-01","2021-09-01","2021-10-01",         -#
#   "2021-11-01","2021-12-01","2022-01-01","2022-02-01","2022-03-01")         -#
# first, adopt time-varying DLNMs to obtain the associations at these times   -#
# for each state.                                                             -#
# cumfit_all: cummulated log-RR ; cumse_all: the standard deviation           -#
# R: associated with the adjacent matrix, the details can be seen in LCAR     -#
# The values in vaccine_PP1V and vaccine_PPFV are ordered by time, state      -#
##############################################################################-#
load("Data for DLNM-LCAR with time-varying associations.Rdata")

################################################################################
####### temperamental LCAR for estimating the modification effect ##############
################################################################################
nstate <- 49; cen <- 12
Cmatrix <- diag(nstate) - R
sdunif="expression:
      logdens=-log_precision/2;
      return(logdens)"
lunif = "expression:
      a = 1;
      b = 1;
      beta = exp(theta)/(1+exp(theta));
      logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
      log_jacobian = log(beta*(1-beta));
      return(logdens+log_jacobian)"
## PP1V ----
poolfitvalues1 <- matrix(NA,nrow = nrow(cumfit_all),ncol = 3,
                         dimnames = list(rownames(cumfit_all),c("beta","betahigh","betalow")))
for (i in 1:nrow(cumse_all)) {
  if(as.numeric(rownames(cumse_all)[i]) == 12){
    poolfitvalues1[i,] <- NA
  } else{
    d <- data.frame(y = cumfit_all[i,],regionID = regionID,timeID = timeID, se = cumse_all[i,])
    formula <- y ~ vaccine_PP1V +
      f(d$regionID, model="generic1", Cmatrix = Cmatrix, constr=TRUE,
        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
      f(d$timeID,model = "rw1")
    
    model <- inla(formula, family="gaussian", data=d,
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"),
                  control.family = list(initial = 0, fixed=TRUE),
                  #control.family = list(hyper = list(prec = prec.prior)),
                  scale = 1/se^2,keep = F)
    poolfitvalues1[i,] <- as.numeric(model$summary.fixed[2,c("0.5quant","0.975quant","0.025quant")]) 
  }
  cat(i," ")
}
## PPFV ----
poolfitvalues2 <- matrix(NA,nrow = nrow(cumfit_all),ncol = 3,
                         dimnames = list(rownames(cumfit_all),c("beta","betahigh","betalow")))
for (i in 1:nrow(cumse_all)) {
  if(as.numeric(rownames(cumse_all)[i]) == 12){
    poolfitvalues1[i,] <- NA
  } else{
    d <- data.frame(y = cumfit_all[i,],regionID = regionID,timeID = timeID, se = cumse_all[i,])
    formula <- y ~ vaccine_PPFV +
      f(d$regionID, model="generic1", Cmatrix = Cmatrix, constr=TRUE,
        hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0))) +
      f(d$timeID,model = "rw1")
    
    model <- inla(formula, family="gaussian", data=d,
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"),
                  control.family = list(initial = 0, fixed=TRUE),
                  #control.family = list(hyper = list(prec = prec.prior)),
                  scale = 1/se^2,keep = F)
    poolfitvalues2[i,] <- as.numeric(model$summary.fixed[2,c("0.5quant","0.975quant","0.025quant")]) 
  }
  cat(i," ")
}

### plot curve-----
tiff("average curve-timevary-vaccine.tiff",width = 16,height = 11,units = "cm",res = 300)
aa <- setdiff(seq(-10,30,by=1),cen) %>% as.character()
poolRR1 <- exp(poolfitvalues1) %>% as.data.frame()
poolRR1 <- poolRR1[aa,]
xvar <- rownames(poolRR1) %>% as.numeric()
poolRR2 <- exp(poolfitvalues2) %>% as.data.frame()
poolRR2 <- poolRR2[aa,]

oldpar <- par(cex.axis=1,cex.lab=1,cex.main=1.2, mar=c(4,4,3.6,0.8))
plot(xvar,y = rep(1,length(xvar)),type="l", ylim=range(cbind(poolRR1,poolRR2),na.rm = T),
     xlab="",ylab="",lty = 2,lwd = 0.5)
mtext("Exp(beta)",side = 2,line = 2,cex = 0.9)
mtext("DAT (degrees Celsius)",side = 1,line = 2.1,cex = 0.9)
points(xvar,poolRR1$beta, type="l", lwd=1, col="blue")
redtrans1 <- rgb(0, 0, 255, 40, maxColorValue=255)
polygon(c(rev(xvar),xvar),c(rev(poolRR1$betahigh),poolRR1$betalow),col=redtrans1, border = NA)

points(xvar,poolRR2$beta, type="l", lwd=1, col="red")
redtrans2 <- rgb(255, 0, 0, 40, maxColorValue=255)  
polygon(c(rev(xvar),xvar),c(rev(poolRR2$betahigh),poolRR2$betalow),col=redtrans2, border = NA)

legend(x="top",inset =0, legend=c("Vaccinated at least one dose (PP1V)", "Fully vaccinated (PPFV)"),
       lwd=1.5, lty=1, col=c("blue", "red"), bty="n",ncol=1, cex=1)
par(oldpar)
dev.off()


