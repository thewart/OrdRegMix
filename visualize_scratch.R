if (!exists("path")) path <- "022119fF/"
source("/home/seth/code/OrdRegMix/loadfromjl.R")


alpha <- julia_eval("gf(foof[4],:α)[:,:,101:1000];")
D <- c("SDB","Feed","Travel","PassiveContact","Agg_give","Avoid","Submit",
       "Groom_high","Groom_low","Approach_high","Approach_low")
D <- ordered(D,levels=D[c(1:3,5:7,4,9,11,8,10)])
dtype <- rep(c("Self-directed","Affiliative","Agonistic","Affiliative"),c(3,1,3,4))
K <- dim(alpha)[2]
niter <- dim(alpha)[3]

pi <- julia_eval("cat(3,[mapslices(softmax,f.tlmm.η,1) for f in foof[4][101:1000]]...);")
n <- dim(pi)[2]
pidat <- data.table(pi=as.vector(pi),Strategy=1:K,FocalID=rep(cov_dat[,unique(FocalID)],each=K),iter=rep(1:niter,each=K*n))
topicord <- pidat[,mean(pi),by=Strategy][,rank(-V1)] %>% ordered()
levels(topicord) <- paste0("S",levels(topicord))
pidat[,Strategy:=NULL]
pidat[,Strategy:=topicord]

alphdat <- data.table(alpha=pnorm(as.vector(alpha)),
                      behav=D,btype=dtype,Strategy=rep(topicord,each=length(D)),
                      iter=rep(1:niter,each=length(D)*K))
alphdat1 <- alphdat[,.(prob=mean(alpha),sdev=sd(alpha),lb=quantile(alpha,0.025),
                       ub=quantile(alpha,0.975)),by=.(behav,Strategy)]
ggplot(alphdat1,aes(x=Strategy,color=Strategy,y=prob,ymin=prob-sdev,ymax=prob+sdev)) + geom_pointrange() +
  facet_wrap(.~behav,scale="free") + ylab("Probability") + theme_gray()


alphdat2 <- alphdat[,.(alpha=alpha/max(alpha),Strategy),
                    by=.(iter,behav,btype)][,.(prob=mean(alpha),lb=quantile(alpha,0.025),
                                               ub=quantile(alpha,0.975),sdev=sd(alpha)),
                                            by=.(behav,Strategy,btype)]
ggplot(alphdat2,aes(y=behav,color=btype,x=prob,xmin=prob-sdev,xmax=prob+sdev)) + 
  geom_point() + geom_errorbarh(height=0) + facet_wrap(.~Strategy) + scale_color_discrete(NULL) + 
  ylab("Behavior") + xlab("Relative likelihood") + theme_gray()

pihat <- pidat[,mean(pi),by=Strategy]
foo <- pidat[,mean(pi),by=.(Strategy,FocalID)] %>% dcast(FocalID~Strategy) 


juh <- used_obs[,mean(passcont),by=.(FocalID)]
juh <- merge(foo,juh)
p3 <- ggplot(juh,aes(y=S2,x=S5,color=Hmisc::cut2(V1,g=3))) + geom_point() + 
  scale_colour_ordinal("Passive contact",labels=c("low","middle","high"),position="bottom") + 
  theme_light()
p1 <- ggplot(juh,aes(y=V1,x=S2)) + geom_point() + geom_smooth(method="lm",color="black") +
  theme_light() + ylab("Passive contact")
p2 <- ggplot(juh,aes(y=V1,x=S5)) + geom_point() + geom_smooth(method="lm",color="black") +
  theme_light() + ylab("Passive contact")

library(cowplot)
p12 <- plot_grid(p1,p2,nrow=2,labels=c("a","b"),vjust=1)
plot_grid(p12,p3+theme(legend.position = "bottom",plot.margin=unit(c(30,5.5,30,5.5),"pt")),
          ncol=2,rel_widths = c(1,1.2),labels=c("","c"),hjust = c(0,-3))




juh <- combodat[,mean(Avoid),by=.(FocalID,Year,ORD_RANK)]
juh <- juh[,V1:=lm(V1~Year+ORD_RANK) %>% resid][,mean(V1),by=FocalID]
juh <- merge(foo,juh)
