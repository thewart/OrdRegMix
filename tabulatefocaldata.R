#specify behaviors and their modifiers ----
basebehav <- c("Scratch","SelfGrm","Travel","Feed","GroomGIVE","GroomGET","passcont",
              "Approach","threat","noncontactAgg","contactAgg","displace","avoid")
modifier <- c(rep(NA,4),rep("",3),"initiate\\(\\w+\\)",
              rep("direct'n\\(\\w+\\)",3),rep("winner\\?\\(\\w+\\)",2))
sb <- !is.na(modifier)
if (modonrank) modifier[sb] <- paste(modifier[sb],"partnerrank\\(\\w+\\)",sep=";")
if (modonsex) modifier[sb] <- paste(modifier[sb],"samesex\\(\\w+\\)",sep=";")
if (modonkin) modifier[sb] <- paste(modifier[sb],"iskin\\(\\w+\\)",sep=";")
modifier <- str_remove(modifier,"^;")
modifier[modifier==""] <- NA
ptetho <- data.table(behavior=basebehav,modifier)

addmodifier(bdat,iskin,kinmat=A)
addmodifier(bdat,samesex,demo=pedigree)
addmodifier(bdat,partnerrank,rank=dominance)

# extract and count specified basebehavs ----
bdat[,Behavior:=eventsplit(Behavior,BehaviorModifier,ptetho)]
behav_withmods <- bdat[,eventslices(unique(Behavior),ptetho$behavior)]

all_obs <- addupbehaviors(behav_withmods,bdat) %>% 
  merge(x=unique(bdat[,c("Observation","FocalID","Observer","Year","Group")]),by="Observation")
setkey(all_obs,"FocalID","Year","Observer")

#filter unwanted columns
nacols <- str_subset(names(all_obs),"NA")
if (length(nacols>0)) all_obs <- all_obs[,-nacols,with=F]
if (modonsex) all_obs <- all_obs[,-str_subset(names(all_obs),"samesex\\(FALSE\\)"),with=F]

all_obs[,SDB:= SelfGrm + Scratch]
if (!modonkin & !modonrank) {
  all_obs[,NonConAgg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)`]
  all_obs[,NonConAgg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)`]
  all_obs[,Agg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)` + `contactAgg:direct'n(give)`]
  all_obs[,Agg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)` + `contactAgg:direct'n(receive)`]
  all_obs[,Avoid_give:=`displace:winner?(foc)` + `avoid:winner?(foc)`]
  all_obs[,Avoid_rec:=`displace:winner?(partnr)` + `avoid:winner?(partnr)`]
}

if (modonkin) {
  all_obs[,SDB:= SelfGrm + Scratch]
  all_obs[,Agg_kin:=`threat:direct'n(give):samesex(TRUE):iskin(TRUE)` +
            `noncontactAgg:direct'n(give):samesex(TRUE):iskin(TRUE)` +
            `contactAgg:direct'n(give):samesex(TRUE):iskin(TRUE)`]
  all_obs[,Agg_nonkin:=`threat:direct'n(give):samesex(TRUE):iskin(FALSE)` +
            `noncontactAgg:direct'n(give):samesex(TRUE):iskin(FALSE)` +
            `contactAgg:direct'n(give):samesex(TRUE):iskin(FALSE)`]
}
