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
# needs mlinit.R to load cov_dat
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

b <- "passcont"
juh <- used_obs[,mean(get(b)),by=.(FocalID)]
juh <- merge(foo,juh)
juh[,tritile:=Hmisc::cut2(V1,g=3)]
pp <- ggplot(juh,aes(y=S2,x=S4,color=tritile)) + geom_point() + 
  scale_colour_ordinal("Passive \n contact",labels=c("low","middle","high"),position="bottom") + 
  theme_light()
b <- "Feed"
juh <- used_obs[,mean(get(b)),by=.(FocalID)]
juh <- merge(foo,juh)
juh[,tritile:=Hmisc::cut2(V1,g=3)]
pf <- ggplot(juh,aes(y=S1,x=S5,color=tritile)) + geom_point() + 
  scale_colour_ordinal("Feeding",labels=c("low","middle","high"),position="bottom") + 
  theme_light()

plot_grid(pp,pf,nrow=2,labels=c("a","b"))

library(cowplot)
p12 <- plot_grid(p1,p2,nrow=2,labels=c("a","b"),vjust=1)
plot_grid(p12,p3+theme(legend.position = "bottom",plot.margin=unit(c(30,5.5,30,5.5),"pt")),
          ncol=2,rel_widths = c(1,1.2),labels=c("","c"),hjust = c(0,-3))

juh <- combodat[,mean(Avoid),by=.(FocalID,Year,ORD_RANK)]
juh <- juh[,V1:=lm(V1~Year+ORD_RANK) %>% resid][,mean(V1),by=FocalID]
juh <- merge(foo,juh)

mustrat <- pidat[,mean(pi),by=.(iter,Strategy)][,.(mean(V1),quantile(V1,0.025),
                                                   quantile(V1,0.975)),by=Strategy]
ggplot(mustrat,aes(y=V1,x=Strategy,ymin=V2,ymax=V3)) + geom_pointrange()

sdstrat <- pidat[,sd(pi/mean(pi)),by=.(iter,Strategy)][,.(mean(V1),quantile(V1,0.025),
                                                   quantile(V1,0.975)),by=Strategy]
ggplot(sdstrat,aes(y=V1,x=Strategy,ymin=V2,ymax=V3)) + geom_pointrange()

foof <- pidat[,.(mean(pi),sd(pi)),by=.(Strategy,FocalID)]
ggplot(foof,aes(x=reorder_within(FocalID,V1,Strategy),y=V1,ymin=V1-V2,ymax=V1+V2)) +
  geom_pointrange() + facet_wrap(~Strategy,scale="free") +
  scale_x_reordered("FocalID",breaks=NULL) + ylab("Probability")

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}



