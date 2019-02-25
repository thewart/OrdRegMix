if (!exists("path")) path <- "022119fF/"
source("/home/seth/code/OrdRegMix/loadfromjl.R")


alpha <- julia_eval("gf(foof[4],:α)[:,:,101:1000];")
D <- c("SDB","Feed","Travel","passcont","Agg_give","Avoid","Submit",
       "Groom_higher","Groom_lower","Approach_high","Approach_low")
K <- dim(alpha)[2]
niter <- dim(alpha)[3]

pi <- julia_eval("cat(3,[mapslices(softmax,f.tlmm.η,1) for f in foof[4][101:1000]]...);")
n <- dim(pi)[2]
pidat <- data.table(pi=as.vector(pi),K=1:K,FocalID=rep(cov_dat[,unique(FocalID)],each=K),iter=rep(1:niter,each=K*n))
topicord <- pidat[,mean(pi),by=K][,rank(-V1)] %>% ordered()
levels(topicord) <- paste0("S",levels(topicord))
pidat[,K:=NULL]j
pidat[,K:=topicord]

alphdat <- data.table(alpha=pnorm(as.vector(alpha)),
                      behav=D,K=rep(topicord,each=length(D)),iter=rep(1:niter,each=length(D)*K))
alphdat <- alphdat[,.(alpha=mean(alpha),lb=quantile(alpha,0.025),ub=quantile(alpha,0.975)),by=.(behav,K)]
ggplot(alphdat,aes(x=K,color=K,y=alpha,ymin=lb,ymax=ub)) +
  geom_pointrange(position = position_dodge(width=0.5)) +
  facet_wrap(.~ordered(behav,levels=D),scale="free") + scale_y_continuous(limits = c(0,NA))

pihat <- pidat[,mean(pi),by=K]
foo <- pidat[,mean(pi),by=.(K,FocalID)] %>% dcast(FocalID~K) 
ggplot(foo,aes(y=S5,x=S6)) + geom_point()
