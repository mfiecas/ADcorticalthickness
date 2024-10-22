library(parallel)
library(RcppEigen)
library(BH)
library(withr)
library(crayon)
library(ggplot2)
library(gridExtra)
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
load("datal.RData")
load("roiv.RData")
ids=which(dat$fn %in% c("lh.corpuscallosum","rh.corpuscallosum","lh.unknown","rh.unknown"))
dat=dat[-ids,]
fnlist=as.character(dat$fn)
datal=datal[fnlist]
initl=initl[fnlist]
fc=function(fn,datal,initl,dat)
{
  md=stan_model(file="spa.stan")
  mc=sampling(md,data=datal[[fn]],init=list(initl[[fn]]),iter=1e4+500,warmup=500,chains=1)
  save(mc,file=paste0(fn,".RData"))
}
cl=makeCluster(length(fnlist))
clusterEvalQ(cl, library(rstan))
clusterEvalQ(cl,rstan_options(auto_write = TRUE))
clusterExport(cl, c("datal","initl","fnlist","fc","dat"))
clusterSetRNGStream(cl,iseed=100)
clusterEvalQ(cl, RNGkind())
parLapply(cl,fnlist,fc,datal=datal,initl=initl,dat=dat)
stopCluster(cl)