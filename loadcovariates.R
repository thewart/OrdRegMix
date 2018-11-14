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
dominance[,YEAR:= as.character(YEAR)]
dominance[,SEX:=tolower(SEX)]
dominance[,ORD_RANK:=ordered(ORD_RANK,levels=c("L","M","H"))]
#uuh... are you shitting me
dominance <- dominance[paste0(ID,SEX) %in% pedigree[,paste0(ID,SEX)]]

#remove animals with no rank info
bdat <- bdat[paste0(FocalID,Year,Group) %in% dominance[,paste0(ID,YEAR,GROUP)]]

#get kin counts ----
censifiles <- paste0("~/Dropbox/Cayo Census/",c("2013/2013.march_animalels in cs-orig matriline march 2013.xls",
                                                "2014/2014.AUG_animalels in cs-orig matriline August 2014.xls",
                                                "2015/2015.march_animalels in cs-orig matriline mar 2015.xls",
                                                "2016/2016.March_animalels in cs-orig matriline Mar 2016.xls"))
censiyear <- c("2013","2014","2015","2016")

censi <- lapply(censifiles,read_excel)
for (i in 1:length(censi)) censi[[i]]$YEAR <- censiyear[i]
censi <- do.call(rbind,censi) %>% as.data.table()
censi <- censi[paste0(ngroup,YEAR) %in% bdat[,paste0(Group,Year)]]
censi[,SEX:=tolower(SEX)]

# remove animals in focal data who aren't in censi (probably died early in the year)
bdat <- bdat[paste0(FocalID,Year,Group) %in% censi[,paste0(animal_id2,YEAR,ngroup)]]

kincounts <- censi[str_length(animal_id2)==3,countkin(animal_id2,Afull,SEX),by=.(YEAR,ngroup)]
setnames(kincounts,c("YEAR","ngroup"),c("Year","Group"))

kincounts$kin_samesex <- cut(kincounts$kin_samesex,c(0,0.5,quantile(kincounts$kin_samesex,c(0.5,1.0))+0.5),include.lowest = T,ordered_result = T)  
kincounts$kin_diffsex <- cut(kincounts$kin_diffsex,c(0,0.5,quantile(kincounts$kin_diffsex,c(0.5,1.0))+0.5),include.lowest = T,ordered_result = T)  

#extract IDs
id_dat <- unique(bdat[,.(FocalID,Year,Group)])
# merge and calculate age data
id_dat <-  merge(id_dat,pedigree[,.(Focal_ID,SEX,CS.BIRTH.SEASON)], by.x = "FocalID",by.y = "Focal_ID")
id_dat[,Age:=as.numeric(Year) - CS.BIRTH.SEASON]
id_dat[,CS.BIRTH.SEASON:=NULL]

#put it all together
id_dat <- merge(id_dat,dominance[,c("ID","YEAR","GROUP","SEX","ORD_RANK")] %>% unique(),
                by.x=c("FocalID","Year","Group","SEX"),by.y=c("ID","YEAR","GROUP","SEX"))
id_dat <- merge(id_dat,kincounts[,.(ID,Year,kin_samesex)],by.x=c("FocalID","Year"),by.y=c("ID","Year"))
