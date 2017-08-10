source("/home/seth/code/LatentSocialPheno/parsefocaldata.R")
source("/home/seth/code/LatentSocialPheno/getAge.R")
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,"F2013/Txtexports_all_processed.csv")
ptetho <- defaultpoint2()[type!="misc" & !(behavior%in%c("Vigilnce","PsCnTerm","GrmTerm","GrmPrsnt"))]
stetho <- defaultstate2()[type!="misc" & state!="Corral"]
Y <- collectfocal(fpath,ptetho,stetho,group = c("F"))
Y <- Y[FocalID %in% Y[,length(Observation),by=FocalID][V1>10,FocalID]]

ncovcols <- 5

#remove behaviors that happen too infrequenlty and baseline state behaviors
filt <- Y[,lapply(.SD,function(x) mean(x>0)),.SD=eventslices(names(Y),ptetho)] > 0.01
filt <- colnames(filt)[!filt] %>% c(stetho[baseline==T,behavior]) %>% unique()
Y <- Y[,-filt,with=F]
Yraw <- copy(Y)

#discritize!!
Y <- cbind(Y[,1:ncovcols],Y[,lapply(.SD,function(x) 
  cut(x,c(0,quantile(x,c(0.025,0.335,0.665,0.975,1.0))+0.5) %>% unique,include.lowest = T)),
  .SD=-(1:ncovcols)])
Y <- cbind(Y[,1:ncovcols],Y[,lapply(.SD,as.numeric),.SD=-(1:ncovcols)])

#collect covariates
library(xlsx)
drank <- read.xlsx("~/Dropbox/Subjects_attributes, dominance, etc/Dominance Hierarchies/DOMINANCE_ALLSUBJECTS_LONGLIST.xlsx",1) %>% as.data.table()
drank[,ID:=toupper(ID) %>% str_replace_all(.,"_","")]
drank <- drank[YEAR %in% c("2012","2013","2014","2015")]
drank[,ORD_RANK:=ordered(ORD_RANK,levels=c("L","M","H"))]
drank <- drank[,.(rank=as.numeric(ORD_RANK) %>% mean()),by=ID]


Xdf <- data.table(FocalID=Y$FocalID,sex=getsex(Y$FocalID),age=getage(Y$FocalID,Y$Year),group=Y$Group)
Xdf[,age:=mean(unique(age)),by=FocalID]
Xdf <- merge(Xdf,drank,by.x = "FocalID",by.y="ID")
#generate covariate design matrix

#Xdf[group %in% c("BB","KK","S","V"),group:="O"]
X <- model.matrix( ~ sex + poly(age,1) + sex*poly(rank,1),Xdf)[,-1]
#+ group*poly(rank,2) + group*poly(age,2)
X[,2:3] <- apply(X[,2:3],2,function(x) (x-mean(x))/sd(x))
X[,4:5] <- X[,2:3] * X[,1]

Y[!(Observer %in% c("JG","JN")),Observer:="O"]
Z <- cbind(Xdf$FocalID %>% as.factor() %>% as.numeric(), Y$Observer %>% as.factor() %>% as.numeric())

Y <- Y[,-c("ScanProx","avoid:winner?(foc)","avoid:winner?(partnr)","Submit:direct'n(give)","Submit:direct'n(receive)")]


library(rstan)
ordregtop <- stan_model("~/code/OrdRegMix/mordreg_topic.stan")
standat <- list(N=nrow(Y),D=ncol(Y)-ncovcols,L=sapply(Y[,-(1:5)],function(x) length(unique(x))),K=1,Y=Y[,-(1:ncovcols)])
initfit <- optimizing(ordregtop,data=standat,verbose=T,as_vector=F)

Kseq <- rep(c(2,4,8),each=10)
fit_noreg <- foreach(K=Kseq) %dopar% { library(gtools); library(rstan)
  standat$K <- K
  moo <- optimizing(ordregtop,standat,verbose=T,as_vector=F,iter=800)
  return(moo)
}

standat <- list(N=nrow(Y),D=ncol(Y)-ncovcols,P=1,L=sapply(Y[,-(1:5)],function(x) length(unique(x))),K=1,Y=Y[,-(1:ncovcols)],X=matrix(0,nrow(Y),1))
nullfit <- optimizing(ordregtop,data=standat,verbose=T,as_vector=F)

Xfull <- cbind(X,model.matrix( ~ Observer,data=Y))
standat <- list(N=nrow(Y),D=ncol(Y)-ncovcols,P=ncol(Xfull),L=sapply(Y[,-(1:5)],function(x) length(unique(x))),K=1,Y=Y[,-(1:ncovcols)],X=Xfull)
initfit <- optimizing(ordregtop,data=standat,verbose=T,as_vector=F)

Kseq <- rep(c(2,4,8),each=10)
fit_both <- foreach(K=Kseq) %dopar% { library(gtools); library(rstan)
  standat$K <- K
  moo <- optimizing(ordregtop,standat,verbose=T,as_vector=F,iter=1200)
  return(moo)
}
bothdat <- data.table(AIC=sapply(fit_both,function(x) -2*x$value + with(x$par,(length(pi)-1+length(gamma)+length(beta))*2)),K=Kseq)
bothdat <- bothdat[,.(AIC=min(AIC),model="rate+state"),by=K]

statedat <- data.table(AIC=sapply(fit_nocov,function(x) -2*x$value + with(x$par,(length(pi)-1+length(gamma)+length(beta))*2)),K=rep(c(2,4,8),each=5))
statedat <- statedat[,.(AIC=min(AIC),model="state"),by=K]

datdat <- rbind(bothdat,statedat)
datdat <- rbind(datdat,data.table(K=1,AIC=c(-2*nullfit$value+2*with(nullfit$par,length(gamma)),
                                            -2*initfit$value+2*with(initfit$par,length(gamma)+length(beta))),
                                  model=c("state","rate+state")))
