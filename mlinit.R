library(rstan)
library(standardize)
library(lme4)
source("~/code/OrdRegMix/collectcleanfocaldata.R")

# behaviors = c("SDB","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)", "Approach:initiate(partner)",
#               "NonConAgg_give","NonConAgg_rec","contactAgg:direct'n(give)","contactAgg:direct'n(receive)")
# d <- length(behaviors)
# #leftside <- "factor(Year) + Group + SEX + poly(Age,2) + ORD_RANK + (1|FocalID) + (1|Observer)"
# leftside <- "factor(Year) + poly(Age,2) + ORD_RANK + (1|FocalID)"
# #used_obs <- all_obs[Year>2011]
# used_obs <- all_obs[Group=="F" & SEX=="f"]
# 
# Xdf <- unique(used_obs[,.(Group=names(which.max(table(Group))),SEX,Y=mean(SDB)),by=FocalID])
# std <- standardize(Y~Group + SEX,Xdf)
# Xin <- model.matrix(lm(std$formula,std$data))[,-1]
#get design matrix

combodat <- merge(used_obs,cov_dat,by=c("FocalID","Year"))
std <- standardize(as.formula(paste0(behaviors[1]," ~ ",leftside,"+",leftside_ranef)),family = binomial,combodat)
Xf <- model.matrix(lmer(formula=std$formula,data=std$data))[,-1]
Xr <- model.matrix(~ 0 + FocalID,data=combodat)

fit <- list()
fixed <- matrix(nrow=ncol(Xf),ncol=d)
ranid <- matrix(nrow=ncol(Xr),ncol=d)
eta <- matrix(nrow=nrow(used_obs),ncol=d)
baseline <- vector(length = d)
for (i in 1:d) {
  std <- standardize(as.formula(paste0("`",behaviors[i],"` ~ ",leftside,"+",leftside_ranef)),
                     combodat,family=binomial(link="probit"))
  fit[[i]] <- glmer(std$formula,std$data,family=binomial(link="probit"))
  fixed[,i] <- fixef(fit[[i]])[-1]
  ranid[,i] <- ranef(fit[[i]])$Focal[[1]]
  baseline[i] <- fixef(fit[[i]])[1]
  eta[,i] <- Xf %*% fixed[,i] + Xr %*% ranid[,i]
}

#save progress
#save(Xf,Xr1,Xr2,Y,fit,fixed,ranid,ranobs,eta,baseline,cutdiff,file="~/analysis/OrdRegMix/062118fF/clmmfits.Rdat")

# load("~/analysis/OrdRegMix/062518fF/clmmfits.Rdat")

probmodel <- stan_model("~/code/OrdRegMix/probreg_topic.stan")

K <- 2:9
iter <- 50
stanfitlist <- list()
ll <- matrix(nrow=iter,ncol=length(K))
for (i in 1:length(K)) { 
  cl <- readyparallel(4)
  cat(paste0(i,"\n"))
  standat <- list(N=nrow(eta),D=d,K=K[i],Y=Y,eta=eta)
  system.time(juh <- foreach(1:iter) %dopar% {library(rstan); library(gtools)
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

#save for output
path <- "~/analysis/OrdRegMix/071818fF/"
write.matrix(Xf,paste0(path,"Xf.csv"))
write.matrix(Xr,paste0(path,"Xr.csv"))
write.matrix(Y,paste0(path,"Y.csv"))
write.matrix(n,paste0(path,"docrng.csv"))
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
