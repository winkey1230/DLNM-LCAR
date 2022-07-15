################################################################################
## This code is to implement the LCAR in the second stage for the DLNM-LCAR    #
## two-stage strategy.                                                         #
## The codes for DLNM have been well detailed in Package "DLNM", so they are   #
## not included in this file. Only the results from the DLNM are provided      #
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
library(tmap)
library(RColorBrewer)
#library(ggplot2)
#### load data from the first stage, DLNM without time-varying
# cumfit: cummulated log-RR; cumse: the standard deviation for  cummulated log-RR
# ismissing: True indicates the DAT is out of the range
# R: associated with the adjacent matrix, the details can be seen in LCAR model
load("Data for DLNM-LCAR without time-varying associations.Rdata")

################################################################################
####### LCAR for the second stage###############################################
################################################################################
nstate <- ncol(cumfit)
cen <- 12
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
fitvalues <- fitvalueslow <- fitvalueshigh <- matrix(NA,nrow = nrow(cumfit),ncol = nstate,
                                                     dimnames = list(rownames(cumfit),colnames(cumfit)))
poolfitvalues <- matrix(NA,nrow = nrow(cumfit),ncol = 3,
                        dimnames = list(rownames(cumfit),c("average","averagehigh","averagelow")))
performance_CAR <- matrix(NA,nrow = nrow(cumse),ncol = 3,
                          dimnames = list(rownames(cumse),c("pD","DIC","LS")))
for (i in 1:nrow(cumse)) {
  if(as.numeric(rownames(cumse)[i]) == cen){
    poolfitvalues[i,] <- 0
    fitvalues[i,] <- 0
    fitvalueshigh[i,] <- 0
    fitvalueslow[i,] <- 0
  } else{
    d <- data.frame(y = cumfit[i,],regions = 1:nstate, se = cumse[i,])
    formula <- y ~ f(d$regions, model="generic1", Cmatrix = Cmatrix, constr=TRUE,
                     hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
    model <- inla(formula, family="gaussian", data=d,
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"),
                  control.family = list(initial = 0, fixed=TRUE),
                  #control.family = list(hyper = list(prec = prec.prior)),
                  scale = 1/se^2,keep = F)
    poolfitvalues[i,] <- as.numeric(model$summary.fixed[,c("0.5quant","0.975quant","0.025quant")]) 
    fitvalues[i,] <- model$summary.fitted.values$`0.5quant`
    fitvalueshigh[i,] <- model$summary.fitted.values$`0.975quant`
    fitvalueslow[i,] <- model$summary.fitted.values$`0.025quant`
    performance_CAR[i,]<- c(model$dic$p.eff,model$dic$dic,
                            -sum(log(model$cpo$cpo),na.rm = T))
  }
  cat(i," ")
}

################################################################################
################# plot the average exposure-response relationship ##############
###############################################################################
tiff("average curve.tiff",width = 14,height = 12,units = "cm",res = 300)
poolRR <- exp(poolfitvalues) %>% as.data.frame()
x1 <- seq(-20,30,by=1) %>% as.character()
poolRR1 <- poolRR[x1,]
xvar <- rownames(poolRR1) %>% as.numeric()

isNA <- t(ismissing)
spatialRR <- exp(t(fitvalues)) %>% as.data.frame()
spatialRR[isNA] <- NA

oldpar <- par(cex.axis=1,cex.lab=1,cex.main=1.2, mar=c(4,4,3.6,0.8))
plot(xvar,y = rep(1,length(xvar)),type="l", ylim=c(0.6,1.35),
     xlab="DAT (℃)",ylab="Relative Risk (RR)", main="")
points(xvar,poolRR1$average, type="l", lwd=1, col="red")
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
polygon(c(rev(xvar),xvar),c(rev(poolRR1$averagehigh),poolRR1$averagelow),col=redtrans, border = NA)
for (i in 1:nstate) {
  lines(xvar,spatialRR[i,],lty = 3,col = "grey50",lwd = 0.5)
}
par(oldpar)
dev.off()
################################################################################
################# plot spatial distribution of association  ####################
################################################################################
spatialRRdata <- cbind(rownames(spatialRR),spatialRR)
colnames(spatialRRdata) <- c("STUSPS","tem"%+%colnames(spatialRR))
spatialRRhigh <- exp(t(fitvalueshigh))
spatialRRlow <- exp(t(fitvalueslow))

sigdata <- (spatialRRhigh < 1 | spatialRRlow >1) %>% as.data.frame()
xx <- as.matrix(sigdata) 
sigdata[xx] <- "*"
sigdata[!xx] <- " "
sigdata <- cbind(rownames(sigdata),sigdata) 
colnames(sigdata) <- c("STUSPS","sig_tem"%+%colnames(spatialRR))

sigcoldata <- matrix(NA,nrow = nrow(spatialRR),ncol = ncol(spatialRR)) %>% as.data.frame()
sigcoldata[spatialRRhigh < 1] <- "black"
sigcoldata[spatialRRlow > 1] <- "black"
sigcoldata <- cbind(rownames(spatialRR),sigcoldata) 
colnames(sigcoldata) <- c("STUSPS","sigcol_tem"%+%colnames(spatialRR))

spatialRRdata <- multimerge(list(spatialRRdata,sigdata,sigcoldata),by = "STUSPS")

plotdata <- merge(mapdata,spatialRRdata)

tmap_mode("plot")
DAT <- c(-8,-3,2,7,16,20,24,28)
objtem <- "tem"%+%DAT

x <- as.data.frame(plotdata)[,objtem] %>% unlist()
x <- round(range(x,na.rm = T) + c(-0.005,0.005),2)
v <- c(x[1],0.7,0.8,0.9,0.95,1.0,1.05,1.1,1.2,1.3,x[2]+0.01)
colorlables <- NULL
for (i in 2:length(v)) {
  colorlables <- c(colorlables,"["%+%v[i-1]%+%"-"%+%v[i]%+%")")
}
paleta <- brewer.pal(length(colorlables),"RdYlGn")[length(colorlables):1]
panelnames <- "DAT = " %+% DAT
pdf("RRmap-CAR.pdf")
tm_shape(plotdata) + 
  tm_polygons(col = objtem,palette=paleta, title="RR referring to 12", 
              legend.show=T, border.alpha=0.1,legend.reverse=T, style="fixed", 
              breaks=v, interval.closure="left",labels=colorlables) + 
  tm_layout(panel.labels=panelnames,legend.outside=T,legend.title.size = 0.8) +
  tm_text(text = "sig_"%+%objtem,col = "sigcol_"%+%objtem) + 
  tm_add_legend(type = "text",labels = "p < 0.05",
                text = "*",col = "black",size = 15)
dev.off()
################################################################################
#################comparison among LCAR, meta, and stratified analysis###########
################################################################################
## LCAR ----
nstate <- ncol(cumfit)
cen <- 12
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
fitvalues <- fitvalueslow <- fitvalueshigh <- matrix(NA,nrow = nrow(cumfit),ncol = nstate,
                                                     dimnames = list(rownames(cumfit),colnames(cumfit)))
poolfitvalues <- matrix(NA,nrow = nrow(cumfit),ncol = 3,
                        dimnames = list(rownames(cumfit),c("average","averagehigh","averagelow")))
# pD is the number of efficient parameters
performance_CAR <- matrix(NA,nrow = nrow(cumse),ncol = 3,
                          dimnames = list(rownames(cumse),c("pD","DIC","LS")))
for (i in 1:nrow(cumse)) {
  if(as.numeric(rownames(cumse)[i]) == cen){
    poolfitvalues[i,] <- 0
    fitvalues[i,] <- 0
    fitvalueshigh[i,] <- 0
    fitvalueslow[i,] <- 0
  } else{
    d <- data.frame(y = cumfit[i,],regions = 1:nstate, se = cumse[i,])
    formula <- y ~ f(d$regions, model="generic1", Cmatrix = Cmatrix, constr=TRUE,
                     hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
    model <- inla(formula, family="gaussian", data=d,
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"),
                  control.family = list(initial = 0, fixed=TRUE),
                  #control.family = list(hyper = list(prec = prec.prior)),
                  scale = 1/se^2,keep = F)
    poolfitvalues[i,] <- as.numeric(model$summary.fixed[,c("0.5quant","0.975quant","0.025quant")]) 
    fitvalues[i,] <- model$summary.fitted.values$`0.5quant`
    fitvalueshigh[i,] <- model$summary.fitted.values$`0.975quant`
    fitvalueslow[i,] <- model$summary.fitted.values$`0.025quant`
    performance_CAR[i,]<- c(model$dic$p.eff,model$dic$dic,
                            -sum(log(model$cpo$cpo),na.rm = T))
  }
  cat(i," ")
}

## meta----
fitvalues <- fitvalueslow <- fitvalueshigh <- matrix(NA,nrow = nrow(cumfit),ncol = nstate,
                                                     dimnames = list(rownames(cumfit),colnames(cumfit)))
poolfitvalues <- matrix(NA,nrow = nrow(cumfit),ncol = 3,
                        dimnames = list(rownames(cumfit),c("average","averagehigh","averagelow")))
performance_meta <- matrix(NA,nrow = nrow(cumse),ncol = 3,
                           dimnames = list(rownames(cumse),c("pD","DIC","LS")))
Cmatrix_meta <- diag(nrow(Cmatrix))
for (i in 1:nrow(cumse)) {
  if(as.numeric(rownames(cumse)[i]) == cen){
    poolfitvalues[i,] <- 0
    fitvalues[i,] <- 0
    fitvalueshigh[i,] <- 0
    fitvalueslow[i,] <- 0
  } else{
    d <- data.frame(y = cumfit[i,],regions = 1:nstate, se = cumse[i,])
    formula <- y ~ f(d$regions, model="generic0", Cmatrix = Cmatrix_meta, constr=TRUE,
                     hyper=list(prec=list(prior=sdunif)))
    model <- inla(formula, family="gaussian", data=d,
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"),
                  control.family = list(initial = 0, fixed=TRUE),
                  #control.family = list(hyper = list(prec = prec.prior)),
                  scale = 1/se^2,keep = F)
    poolfitvalues[i,] <- as.numeric(model$summary.fixed[,c("0.5quant","0.975quant","0.025quant")]) 
    fitvalues[i,] <- model$summary.fitted.values$`0.5quant`
    fitvalueshigh[i,] <- model$summary.fitted.values$`0.975quant`
    fitvalueslow[i,] <- model$summary.fitted.values$`0.025quant`
    performance_meta[i,]<- c(model$dic$p.eff,model$dic$dic,
                             -sum(log(model$cpo$cpo),na.rm = T))
  }
  cat(i," ")
}


## stratified analysis----
fitvalues <- fitvalueslow <- fitvalueshigh <- matrix(NA,nrow = nrow(cumfit),ncol = nstate,
                                                     dimnames = list(rownames(cumfit),colnames(cumfit)))
poolfitvalues <- matrix(NA,nrow = nrow(cumfit),ncol = 3,
                        dimnames = list(rownames(cumfit),c("average","averagehigh","averagelow")))
performance_independent <- matrix(NA,nrow = nrow(cumse),ncol = 3,
                                  dimnames = list(rownames(cumse),c("pD","DIC","LS")))
for (i in 1:nrow(cumse)) {
  if(as.numeric(rownames(cumse)[i]) == cen){
    poolfitvalues[i,] <- 0
    fitvalues[i,] <- 0
    fitvalueshigh[i,] <- 0
    fitvalueslow[i,] <- 0
  } else{
    d <- data.frame(y = cumfit[i,],regions = 1:nstate, se = cumse[i,])
    formula <- y ~ f(d$regions, model="generic1", Cmatrix = Cmatrix, constr=TRUE,
                     hyper=list(theta1=list(prior="loggamma",param=c(1,100000))))
    model <- inla(formula, family="gaussian", data=d,
                  control.compute=list(dic=TRUE, mlik=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.inla=list(strategy="simplified.laplace"),
                  control.family = list(initial = 0, fixed=TRUE),
                  #control.family = list(hyper = list(prec = prec.prior)),
                  scale = 1/se^2,keep = F)
    poolfitvalues[i,] <- as.numeric(model$summary.fixed[,c("0.5quant","0.975quant","0.025quant")]) 
    fitvalues[i,] <- model$summary.fitted.values$`0.5quant`
    fitvalueshigh[i,] <- model$summary.fitted.values$`0.975quant`
    fitvalueslow[i,] <- model$summary.fitted.values$`0.025quant`
    performance_independent[i,]<- c(model$dic$p.eff,model$dic$dic,
                                    -sum(log(model$cpo$cpo),na.rm = T))
  }
  cat(i," ")
}
# get compare results----
obtem <- seq(-20,30,by=1) %>% as.character()
pD <- cbind(performance_CAR[obtem,"pD"],performance_meta[obtem,"pD"],performance_independent[obtem,"pD"])
DIC <- cbind(performance_CAR[obtem,"DIC"],performance_meta[obtem,"DIC"],performance_independent[obtem,"DIC"])
LS <- cbind(performance_CAR[obtem,"LS"],performance_meta[obtem,"LS"],performance_independent[obtem,"LS"])
performance <- cbind(pD,DIC,LS) %>% round(2) %>% as.data.frame()
performance <- cbind(DAT = as.integer(obtem),performance)
names(performance) <- c("DAT","pD-CAR","pD-meta","pD-sub","DIC-CAR",
                        "DIC-meta","DIC-sub","LS-CAR","LS-meta","LS-sub")
apply(performance, 2, sum,na.rm = T)[-1]




