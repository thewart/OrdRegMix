#source("~/code/brianzhong/Plattlab/master.R")
library(ordinal)
library(rstan)
library(standardize)
library(lme4)
source("~/code/OrdRegMix/messyprep.R")

behaviors = c("SDB","Travel","Feed","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)", "Approach:initiate(partner)","Agg_give","Agg_rec")
d <- length(behaviors)
leftside <- "factor(Year) + Group + SEX + poly(Age,2) + ORD_RANK + (1|FocalID) + (1|Observer)"
#used_obs <- all_obs[Year>2011]
used_obs <- all_obs

#get design matrix
std <- standardize(as.formula(paste0(behaviors[1]," ~ ",leftside)),used_obs)
Xf <- model.matrix(lmer(formula=std$formula,data=std$data))[,-1]
Xr1 <- model.matrix(~ 0 + FocalID,data=used_obs)
Xr2 <- model.matrix(~ 0 + Observer,data=used_obs)
Y <- used_obs[,behaviors,with=F] %>% as.matrix()

fit <- list()
fixed <- matrix(nrow=ncol(Xf),ncol=d)
ranid <- matrix(nrow=ncol(Xr1),ncol=d)
ranobs <- matrix(nrow=ncol(Xr2),ncol=d)
eta <- matrix(nrow=nrow(used_obs),ncol=d)
baseline <- vector(length = d)
cutdiff <- vector(length=d)
for (i in 1:d) {
  std <- standardize(as.formula(paste0("ordered(`",behaviors[i],"`) ~ ",leftside)),used_obs,family=poisson)
  fit[[i]] <- clmm(std$formula,std$data,link="probit")
  fixed[,i] <- fit[[i]]$coefficients[-(1:2)]
  ranid[,i] <- ranef(fit[[i]])$Focal[[1]]
  ranobs[,i] <- ranef(fit[[i]])$Observer[[1]]
  baseline[i] <- fit[[i]]$coefficients[1]
  cutdiff[i] <- fit[[i]]$coeff[2] - fit[[i]]$coeff[1]
  eta[,i] <- Xf %*% fixed[,i] + Xr1 %*% ranid[,i] + Xr2 %*% ranobs[,i]
  #eta[,i] <- eta[,i]-baseline[i]
}

nfold <- 20
fold <- sample(rep(1:nfold,nrow(Y)/nfold))

ordmodel <- stan_model("~/code/OrdRegMix/mordreg_topic.stan")
ordmodel_cv <- stan_model("~/code/OrdRegMix/mordreg_topic_cv.stan")
K <- 5
standat <- list(N=nrow(eta),D=d,K=K,Y=Y,
                eta=eta,L=rep(3,d))
system.time(juh <- foreach(1:8) %dopar% {library(rstan); library(gtools)
  standat$K <- K
  init <- list(gamma=-rep(baseline,K) %>% array(dim = c(d,K)) + rnorm(K*d,sd = 0.5),
               pi=MCMCpack::rdirichlet(1,rep(1,K)) %>% as.vector(),
               cutdiff=cutdiff)
  stanfit <- optimizing(ordmodel,standat,init=init,as_vector=F,verbose=T,tol_obj=0.01)
  return(stanfit)
})
init <- juh[[which.min(sapply(juh,function(x) x$val))]]$par

juh <- foreach(i=1:nfold) %dopar% {library(rstan); library(gtools)
  standat <- list(N=sum(fold!=i),D=d,K=K,Y=Y[fold!=i,],eta=eta[fold!=i,],L=rep(3,d),
                  N_out=sum(fold==i),Y_out=Y[fold==i,],eta_out=eta[fold==i,])
  stanfit <- optimizing(ordmodel_cv,standat,init=init,as_vector=F,verbose=T,tol_obj=0.01)
}

K <- 3
ordmodel <- stan_model("~/code/OrdRegMix/mordreg_topic.stan")
standat <- list(N=nrow(eta),D=d,K=K,Y=Y,
                eta=eta,L=rep(3,d))
init <- list(gamma=-rep(baseline,K) %>% array(dim = c(d,K)) + rnorm(K*d,sd = 0.5),
             pi=MCMCpack::rdirichlet(1,rep(1,K)) %>% as.vector(),
             cutdiff=cutdiff)
system.time(stanfit <- optimizing(ordmodel,standat,init=init,as_vector=F,verbose=T))


standat$K <- K+1

init <- list()
init$gamma <- cbind(stanfit$par$gamma,baseline + rnorm(d,sd=0))
init$pi <- c(stanfit$par$pi,1e-6)
init$pi <- init$pi/sum(init$pi)
init$cutdiff <- stanfit$par$cutdiff

stanfit_kp1 <- optimizing(ordmodel,standat,init=init,as_vector=F,verbose=T)


K <- 5
ordmodel <- stan_model("~/code/OrdRegMix/mordreg_topic_dir.stan")
standat <- list(N=nrow(eta),D=d,K=K,Y=Y,
                eta=eta,L=rep(3,d),alpha=1)
init <- list(gamma=-rep(baseline,K) %>% array(dim = c(d,K)) + rnorm(K*d,sd = 0.1),
             pi=MCMCpack::rdirichlet(1,rep(1,K)) %>% as.vector(),
             cutdiff=cutdiff)
stanfit <- optimizing(ordmodel,standat,init=init,as_vector=F,verbose=T)

init$gamma <- init$gamma[,1:9]
init$pi <- init$pi[1:9]/sum(init$pi[1:9])
standat$K <- 9

sig <- 1/stanfit$par$cutdiff %>% as.vector()
ranobsf <- sweep(ranobs,2,sig,"*")
ranidf <- sweep(ranid,2,sig,"*")
fixedf <- sweep(fixed,2,sig,"*")
alphaf <- t(stanfit$par$gamma)*sig

path <- "~/analysis/OrdRegMix/"
write.matrix(sig,paste0(path,"sigma_0.csv"))
write.matrix(ranid,paste0(path,"ranef1_0.csv"))
write.matrix(ranobs,paste0(path,"ranef2_0.csv"))
write.matrix(fixed,paste0(path,"fixef_0.csv"))
write.matrix(stanfit$par$r,paste0(path,"r_0.csv"))
write.matrix(stanfit$par$gamma,paste0(path,"alpha_0.csv"))

write.matrix(Xf,paste0(path,"Xf.csv"))
write.matrix(Xr1,paste0(path,"Xr1.csv"))
write.matrix(Xr2,paste0(path,"Xr2.csv"))
write.matrix(Y,paste0(path,"Y.csv"))
write.matrix(used_obs[,length(`Observation`),by=FocalID]$V1,paste0(path,"docrng.csv"))

