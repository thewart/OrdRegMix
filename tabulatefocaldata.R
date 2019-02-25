#specify behaviors and their modifiers ----
basebehav <- c("Scratch","SelfGrm","Travel","Feed","GroomGIVE","GroomGET","passcont",
              "Approach","threat","noncontactAgg","contactAgg","Submit","FearGrm","displace","avoid")
modifier <- c(rep(NA,4),rep("",3),"initiate\\(\\w+\\)",
              rep("direct'n\\(\\w+\\)",5),rep("winner\\?\\(\\w+\\)",2))
sb <- !is.na(modifier)
if (modonsex) modifier[sb] <- paste(modifier[sb],"samesex\\(\\w+\\)",sep=";")
if (modonrank) modifier[sb] <- paste(modifier[sb],"partnerrank\\(\\w+\\)",sep=";")
if (modonkin) modifier[sb] <- paste(modifier[sb],"iskin\\(\\w+\\)",sep=";")
modifier <- str_remove(modifier,"^;")
modifier[modifier==""] <- NA
ptetho <- data.table(behavior=basebehav,modifier)

addmodifier(bdat,iskin,kinmat=A)
addmodifier(bdat,partnerrank,rank=dominance)
addmodifier(bdat,samesex,demo=pedigree)


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
