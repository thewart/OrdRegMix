#specify behaviors and their modifiers ----
behavior <- c("Scratch","SelfGrm","GroomGIVE","GroomGET","passcont","threat","noncontactAgg","contactAgg","Approach")
modifier <- c(NA,NA,"","","",rep("direct'n\\(\\w+\\)",3),"initiate\\(\\w+\\)")
sb <- !is.na(modifier)
if (modonrank) modifier[sb] <- paste(modifier[sb],"partnerrank\\(\\w+\\)",sep=";")
if (modonsex) modifier[sb] <- paste(modifier[sb],"samesex\\(\\w+\\)",sep=";")
if (modonkin) modifier[sb] <- paste(modifier[sb],"iskin\\(\\w+\\)",sep=";")
modifier <- str_remove(modifier,"^;")
ptetho <- data.table(behavior,modifier)

addmodifier(bdat,iskin,kinmat=A)
addmodifier(bdat,samesex,demo=pedigree)
addmodifier(bdat,partnerrank,rank=dominance)


# extract and count specified behaviors ----
bdat[,Behavior:=eventsplit(Behavior,BehaviorModifier,ptetho)]
all_obs <- countprep(ptetho$behavior,bdat) %>% 
  merge(x=unique(bdat[,c("Observation","FocalID","Observer","Year","Group")]),by="Observation")
setkey(all_obs,"FocalID","Year","Observer")

#filter unwanted columns
all_obs <- all_obs[,-str_subset(names(all_obs),"NA"),with=F]
if (modonsex) all_obs <- all_obs[,-str_subset(names(all_obs),"samesex\\(FALSE\\)"),with=F]

# all_obs[,NonConAgg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)`]
# all_obs[,NonConAgg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)`]
# all_obs[,Agg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)` + `contactAgg:direct'n(give)`]
# all_obs[,Agg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)` + `contactAgg:direct'n(receive)`]
# all_obs[,Avoid_give:=`displace:winner?(foc)` + `avoid:winner?(foc)`]
# all_obs[,Avoid_rec:=`displace:winner?(partnr)` + `avoid:winner?(partnr)`]
all_obs[,SDB:= SelfGrm + Scratch]
all_obs[,Agg_kin:=`threat:direct'n(give):samesex(TRUE):iskin(TRUE)` +
          `noncontactAgg:direct'n(give):samesex(TRUE):iskin(TRUE)` +
          `contactAgg:direct'n(give):samesex(TRUE):iskin(TRUE)`]
all_obs[,Agg_nonkin:=`threat:direct'n(give):samesex(TRUE):iskin(FALSE)` +
          `noncontactAgg:direct'n(give):samesex(TRUE):iskin(FALSE)` +
          `contactAgg:direct'n(give):samesex(TRUE):iskin(FALSE)`]

# loop quantile function for all behaviors ----
blevels <- list()
for (i in 6:ncol(all_obs)) {
  all_obs[[i]] <- cut(all_obs[[i]],c(0,quantile(all_obs[[i]],c(0.5,1.0))+0.5),
                      include.lowest = T,ordered_result = T)
  blevels[[i-5]] <- levels(all_obs[[i]])
  all_obs[[i]] <- as.numeric(all_obs[[i]]) - 1
}

