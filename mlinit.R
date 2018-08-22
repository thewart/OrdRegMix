#source("~/code/brianzhong/Plattlab/master.R")
library(rstan)
library(standardize)
library(lme4)
source("~/code/OrdRegMix/messyprep.R")

behaviors = c("SDB","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)", "Approach:initiate(partner)",
              "NonConAgg_give","NonConAgg_rec","contactAgg:direct'n(give)","contactAgg:direct'n(receive)")
d <- length(behaviors)
#leftside <- "factor(Year) + Group + SEX + poly(Age,2) + ORD_RANK + (1|FocalID) + (1|Observer)"
leftside <- "factor(Year) + poly(Age,2) + ORD_RANK + (1|FocalID)"
#used_obs <- all_obs[Year>2011]
used_obs <- all_obs[Group=="F" & SEX=="f"]

Xdf <- unique(used_obs[,.(Group=names(which.max(table(Group))),SEX,Y=mean(SDB)),by=FocalID])
std <- standardize(Y~Group + SEX,Xdf)
Xin <- model.matrix(lm(std$formula,std$data))[,-1]
#get design matrix
std <- standardize(as.formula(paste0(behaviors[1]," ~ ",leftside)),family = binomial,used_obs)
Xf <- model.matrix(lmer(formula=std$formula,data=std$data))[,-1]
Xr1 <- model.matrix(~ 0 + FocalID,data=used_obs)
# Xr2 <- model.matrix(~ 0 + Observer,data=used_obs)
Y <- used_obs[,behaviors,with=F] %>% as.matrix()

fit <- list()
fixed <- matrix(nrow=ncol(Xf),ncol=d)
ranid <- matrix(nrow=ncol(Xr1),ncol=d)
# ranobs <- matrix(nrow=ncol(Xr2),ncol=d)
eta <- matrix(nrow=nrow(used_obs),ncol=d)
baseline <- vector(length = d)
# cutdiff <- vector(length=d)
for (i in 1:d) {
  # std <- standardize(as.formula(paste0("ordered(`",behaviors[i],"`) ~ ",leftside)),used_obs,family=poisson)
  # fit[[i]] <- clmm(std$formula,std$data,link="probit")
  std <- standardize(as.formula(paste0("`",behaviors[i],"` ~ ",leftside)),used_obs,family=binomial(link="probit"))
  fit[[i]] <- glmer(std$formula,std$data,family=binomial(link="probit"))
  # fixed[,i] <- fit[[i]]$coefficients[-(1:2)]
  fixed[,i] <- fixef(fit[[i]])[-1]
  ranid[,i] <- ranef(fit[[i]])$Focal[[1]]
  # ranobs[,i] <- ranef(fit[[i]])$Observer[[1]]
  # baseline[i] <- fit[[i]]$coefficients[1]
  baseline[i] <- fixef(fit[[i]])[1]
  # cutdiff[i] <- fit[[i]]$coeff[2] - fit[[i]]$coeff[1]
  eta[,i] <- Xf %*% fixed[,i] + Xr1 %*% ranid[,i]
  #eta[,i] <- eta[,i]-baseline[i]
}

#save progress
#save(Xf,Xr1,Xr2,Y,fit,fixed,ranid,ranobs,eta,baseline,cutdiff,file="~/analysis/OrdRegMix/062118fF/clmmfits.Rdat")

# load("~/analysis/OrdRegMix/062518fF/clmmfits.Rdat")

nfold <- 8
fold <- sample(rep(1:nfold,ceiling(nrow(Y)/nfold)),size=nrow(Y))
ordmodel_cv <- stan_model("~/code/OrdRegMix/mordreg_topic_cv.stan")
ordmodel <- stan_model("~/code/OrdRegMix/mordreg_topic.stan")
probmodel <- stan_model("~/code/OrdRegMix/probreg_topic.stan")

K <- 2:9
stanfitlist <- list()
ll <- matrix(nrow=100,ncol=length(K))
for (i in 1:length(K)) { 
  cl <- readyparallel(4)
  cat(paste0(i,"\n"))
  standat <- list(N=nrow(eta),D=d,K=K[i],Y=Y,eta=eta)
  system.time(juh <- foreach(1:100) %dopar% {library(rstan); library(gtools)
    init <- list(gamma=rep(baseline,K[i]) %>% array(dim = c(d,K[i])) +
    rnorm(K[i]*d,sd = 1), pi=MCMCpack::rdirichlet(1,rep(1,K[i])) %>% array())
    #              cutdiff=cutdiff)
    stanfit <- optimizing(probmodel,standat,init=init,as_vector=F,verbose=T,tol_obj=0.01)
    return(stanfit)
  })
  ll[,i] <- sapply(juh,function(x) x$val)
  stanfitlist[[i]] <- juh[[which.max(ll[,i])]]
  stopCluster(cl)
}

#cross-validation
cv <- matrix(nrow=nfold,ncol=length(K))
for (j in 1:length(K)) {
  cl <- readyparallel(4)
  juh <- foreach(i=1:nfold) %dopar% {library(rstan); library(gtools)
    standat <- list(N=sum(fold!=i),D=d,K=K[j],Y=Y[fold!=i,],eta=eta[fold!=i,],L=rep(3,d),
                    N_out=sum(fold==i),Y_out=Y[fold==i,],eta_out=eta[fold==i,])
    stanfit <- optimizing(ordmodel_cv,standat,init=stanfitlist[[j]]$par,as_vector=F,verbose=T,tol_obj=0.01)
  }
  cv[,j] <- sapply(juh,function(x) x$par$loglik_out)
  stopCluster(cl)
}

#save for output
path <- "~/analysis/OrdRegMix/071818fF/"
write.matrix(Xf,paste0(path,"Xf.csv"))
write.matrix(Xr1,paste0(path,"Xr1.csv"))
write.matrix(Y,paste0(path,"Y.csv"))
write.matrix(used_obs[,length(`Observation`),by=FocalID]$V1,paste0(path,"docrng.csv"))
#write.matrix(Xin,paste0(path,"Xin.csv"))

suff <- "_k1.csv"
write.matrix(baseline,paste0(path,"alpha_0",suff))
write.matrix(ranid,paste0(path,"ranef1_0",suff))
write.matrix(fixed,paste0(path,"fixef_0",suff))

for (i in 1:length(K)) {
  stanfit <- stanfitlist[[i]]
  # sig <- 1/stanfit$par$cutdiff %>% as.vector()
  sig <- rep(1,ncol(Y))
  ranidf <- sweep(ranid,2,sig,"*")
  fixedf <- sweep(fixed,2,sig,"*")
  alphaf <- t(t(stanfit$par$gamma)*sig)
  
  suff <- paste0("_k",K[i],".csv")
  # write.matrix(sig,paste0(path,"sigma_0",suff))
  write.matrix(ranidf,paste0(path,"ranef1_0",suff))
  write.matrix(fixedf,paste0(path,"fixef_0",suff))
  write.matrix(stanfit$par$r,paste0(path,"r_0",suff))
  write.matrix(alphaf,paste0(path,"alpha_0",suff))
}
