sourcepath <- "/home/seth/code/OrdRegMix/modelcomp.R"
path <- "121518fF_kin/"
source(sourcepath)
repdat_kin <- copy(repdat)
waicdat_kin <- copy(waicdat)
repdat_kin$data <- "Kin"
waicdat_kin$data <- "Kin"

path <- "121518fF_rank/"
source(sourcepath)
repdat_rank <- copy(repdat)
waicdat_rank <- copy(waicdat)
repdat_rank$data <- "Rank"
waicdat_rank$data <- "Rank"

path <- "121518fF/"
source(sourcepath)
repdat$data <- "None"
waicdat$data <- "None"

repdat_tlmm <- rbind(repdat,repdat_rank,repdat_kin)
waicdat_tlmm <- rbind(waicdat,waicdat_rank,waicdat_kin)

wptlmm <- ggplot(waicdat_tlmm[K!=1 | model!="Rates & States"],aes(x=K,y=elpd_diff,color=model)) + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se),position = position_dodge(width=0.2),fatten=2.0) +
  geom_line(position=position_dodge(width=0.2)) + scale_color_discrete("Variability in...") + xlab("Number of states") + ylab(expression(paste(Delta,"WAIC"))) +
  facet_wrap(. ~ data)

rptlmm <- ggplot(repdat_tlmm[K!=1 | model!="Rates & States"],aes(y=rep,x=K,ymin=rep-sd,ymax=rep+sd,color=model)) +
  geom_line(position=position_dodge(0.25)) + geom_pointrange(position = position_dodge(0.25)) +
  ylim(c(0,0.05)) + xlab("Number of states") + ylab("Repeatability") + scale_color_discrete("Variability in...") + 
  facet_wrap(.~data)
dcast(repdat_tlmm,K + data ~ model,value.var="rep")[,(`Rates & States`- `Rates only`)/ `Rates & States`,by=.(K,data)]


comppath <- "/home/seth/code/MultVarBinom/modelcomp.R"
reppath <- "/home/seth/code/MultVarBinom/repeatability.R"
ffile <- "mvbfits_kin.Rdat"
source(comppath)
source(reppath)
repdat_kin <- copy(repdat)
waicdat_kin <- copy(waicdat)
repdat_kin$data <- "Kin"
waicdat_kin$data <- "Kin"

ffile <- "mvbfits_rank.Rdat"
source(comppath)
source(reppath)
repdat_rank <- copy(repdat)
waicdat_rank <- copy(waicdat)
repdat_rank$data <- "Rank"
waicdat_rank$data <- "Rank"

ffile <- "mvbfits.Rdat"
source(comppath)
source(reppath)
repdat$data <- "None"
waicdat$data <- "None"

repdat_mvb <- rbind(repdat,repdat_rank,repdat_kin)
waicdat_mvb <- rbind(waicdat,waicdat_rank,waicdat_kin)

wpmvb <- ggplot(waicdat_mvb,aes(y=elpd_diff,x=data,color=model)) + 
  geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se),position=position_dodge(0.25)) + 
  scale_color_discrete("Variability in...") + xlab("") + ylab(expression(paste(Delta,"WAIC")))

rpmvb <- ggplot(repdat_mvb,aes(y=rep,x=data,ymin=rep-sd,ymax=rep+sd,color=model)) +
  geom_pointrange(position = position_dodge(0.25)) +
  ylim(c(0,0.05)) + xlab("") + ylab("Repeatability") + scale_color_discrete("Variability in...")

dcast(repdat_mvb,data ~ model,value.var="rep")[,(`Rates & Cors`- `Rates only`)/ `Rates & Cors`,by=data]
