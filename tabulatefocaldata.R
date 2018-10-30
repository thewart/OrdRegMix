#load info for partner kin ----
pedigree <- read_excel("Dropbox/Pedigree and Life-History Data/PEDIGREE_updated2016.xlsx", sheet = "Demographic Data")
pedigree <- as.data.table(pedigree)
pedigree[is.na(DAM),DAM:=`BEHAVIORAL MOM`]
pedigree[,SEX:=tolower(SEX)]

Afull <- pedigree[,2*kinship2::kinship(Focal_ID,dadid=SIRE,momid=DAM)]
matchup <- match(
  unique(c(bdat$FocalID,str_extract(bdat$PartnerID,"^\\w{3}$") %>% na.omit)),
  rownames(Afull))
A <- Afull[matchup,matchup]

#load info for partner rank ----
dominance <- read_excel("~/Dropbox/Subjects_attributes, dominance, etc/Dominance Hierarchies/DOMINANCE_ALLSUBJECTS_LONGLIST.xlsx", sheet = "DOMINANCE_ALLSUBJECTS") %>% as.data.table()
dominance$YEAR <- as.character(dominance$YEAR)
dominance[,ORD_RANK:=ordered(ORD_RANK,levels=c("L","M","H"))]

# extract and count relevant behaviors ----
ptetho <- defaultpoint2()[type!="misc" & !(behavior%in%c("Vigilnce","PsCnTerm","GrmTerm","GrmPrsnt"))]
ptetho <- rbind(ptetho,list(c("GroomGIVE","GroomGET"),rep(NA,2),rep("Affiliative",2)))
bdat[,Behavior:=eventsplit(Behavior,BehaviorModifier,ptetho)]
all_obs <- countprep(ptetho$behavior,bdat) %>% 
  merge(x=unique(bdat[,c("Observation","FocalID","Observer","Year","Group")]),by="Observation")
setkey(all_obs,"FocalID","Year","Observer")

all_obs[,NonConAgg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)`]
all_obs[,NonConAgg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)`]
all_obs[,Agg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)` + `contactAgg:direct'n(give)`]
all_obs[,Agg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)` + `contactAgg:direct'n(receive)`]
all_obs[,Avoid_give:=`displace:winner?(foc)` + `avoid:winner?(foc)`]
all_obs[,Avoid_rec:=`displace:winner?(partnr)` + `avoid:winner?(partnr)`]
all_obs[,SDB:= SelfGrm + Scratch]

# loop quantile function for all behaviors
blevels <- list()
for (i in 6:ncol(all_obs)) {
  # all_obs[[i]] <- cut(all_obs[[i]],c(0,0.5,quantile(all_obs[[i]][all_obs[[i]]>0],c(0.5,1.0))+0.5) %>% unique(),include.lowest = T)
  all_obs[[i]] <- cut(all_obs[[i]],c(0,quantile(all_obs[[i]],c(0.5,1.0))+0.5),include.lowest = T)  
  blevels[[i-5]] <- levels(all_obs[[i]])
  # all_obs[[i]] <- as.numeric(all_obs[[i]])
  all_obs[[i]] <- as.numeric(all_obs[[i]]) - 1
}