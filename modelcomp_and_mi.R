library(loo)
if (!exists("path")) path <- "022119fF/"
source("/home/seth/code/OrdRegMix/loadfromjl.R")
iter <- "101:1000"
l <- julia_eval("length(lp)")

##### waic
waic_rs <- waic_r <- waic_0 <- list()
Ks <- vector()
for (i in 1:l) {
  is <- as.character(i)
  
  waic_rs[[i]] <- waic(julia_eval(paste0("lp[",is,"][:,",iter,"]';")))
  waic_r[[i]] <- waic(julia_eval(paste0("lpu[",is,"][:,",iter,"]';")))
  waic_0[[i]] <- waic(julia_eval(paste0("lp0[",is,"][:,",iter,"]';")))
  
  Ks[i] <- julia_eval(paste0("size(getfield(foof[",is,"][1],fieldname(HYBRIDsample,1)),2);"))
}

##### repeatability
mi_rs <- mi_r <- vector()
sdev_rs <- sdev_r <- vector()
Ks  <- vector()
for (i in 1:l) {
  cat(i,"\r")
  is <- as.character(i)
  
  prs <- julia_eval(paste0("calc_py(foof[",is,"][",iter,"]);"))
  rep_rs <- apply(prs,3,calc_mi,normalize=F)
  rep_rs_marg <- apply(prs,3,calc_mi_marg,normalize=F)
  
  midiff <- (rep_rs-rep_rs_marg)/rep_rs_marg
  mi_rs[i] <- mean(midiff)
  sdev_rs[i] <- sd(midiff)
  
  # pr <- julia_eval(paste0("calc_py(foofu[",is,"][",iter,"]);"))
  # rep_r <- apply(pr,3,calc_mi,normalize=F)
  # rep_r_marg <- apply(pr,3,calc_mi_marg,normalize=F)
  # 
  # midiff <- (rep_r-rep_r_marg)/rep_r_marg
  # mi_r[i] <- mean(midiff)
  # sdev_r[i] <- sd(midiff)
  
  Ks[i] <- julia_eval(paste0("size(getfield(foof[",is,"][1],fieldname(HYBRIDsample,1)),2);"))
}

midat_bsr <- data.table(mi=mi_rs,sdev=sdev_rs,model=paste0("BS",Ks))


