library(loo)
if (!exists("path")) path <- "022119fF/"
source("/home/seth/code/OrdRegMix/loadfromjl.R")
iter <- "101:1000"
l <- julia_eval("length(lp)")

baseline <- waic(julia_eval(paste0("lp0[1][:,",iter,"]'")))
waicdat <- data.table(elpd_diff=vector(),se=vector(),K=vector())
rsvrdat <- data.table(elpd_diff=vector(),se=vector(),K=vector())
for (i in 1:l) {
  is <- as.character(i)
  
  wcomp <- compare(baseline,waic(julia_eval(paste0("lp[",is,"][:,",iter,"]';"))))
  wucomp <- compare(baseline,waic(julia_eval(paste0("lpu[",is,"][:,",iter,"]';"))))
  w0comp <- compare(baseline,waic(julia_eval(paste0("lp0[",is,"][:,",iter,"]';"))))
  
  tmp <- rbind(wcomp,wucomp,w0comp) %>% as.data.table()
  tmp$K <- julia_eval(paste0("size(getfield(foof[",is,"][1],fieldname(HYBRIDsample,1)),2);"))
  waicdat <- rbind(waicdat,tmp)
}
waicdat$model <- c("Rates & States","Rates only","Nothing")
ggplot(waicdat[K!=1 | model!="Rates & States"],aes(x=K,y=elpd_diff,color=model)) + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se),position = position_dodge(width=0.2),fatten=2.0) +
  geom_line(position=position_dodge(width=0.2)) + scale_color_discrete("Variability in...") + xlab("Number of states") + ylab(expression(paste(Delta,"WAIC")))

repdat <- data.table(rep=vector(),sd=vector(),K=vector())
for (i in 1:l) {
  cat(i,"\r")
  is <- as.character(i)
  repf <- julia_eval(paste0("calc_py(foof[",is,"][",iter,"]);")) %>% apply(3,calc_repeat)
  repu <- julia_eval(paste0("calc_py(foofu[",is,"][",iter,"]);")) %>% apply(3,calc_repeat)
  tmp <- data.table(rep=c(mean(repf),mean(repu)),sd = c(sd(repf),sd(repu)))
  tmp$K <- julia_eval(paste0("size(getfield(foof[",is,"][1],fieldname(HYBRIDsample,1)),2);"))
  repdat <- rbind(repdat,tmp)
}
repdat$model <- c("Rates & States","Rates only")
# ggplot(repdat[K!=1 | model!="Rates & states"],aes(y=rep,x=K,ymin=rep-sd,ymax=rep+sd,color=model)) +
#   geom_line(position=position_dodge(0.25)) + geom_pointrange(position = position_dodge(0.25)) + 
#   ylim(c(0,0.05)) + xlab("Number of states") + ylab("Repeatability") + scale_color_discrete("Variability in...")
# dcast(repdat,K ~ model,value.var="rep")[,(`Rates & States`- `Rates only`)/ `Rates & States`,by=K]  


