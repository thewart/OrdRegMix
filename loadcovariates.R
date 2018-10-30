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
# remove animals in focal data who aren't in censi (probably died early in the year)
bdat <- bdat[paste0(FocalID,Year,Group) %in% censi[,paste0(animal_id2,YEAR,ngroup)]]

kincounts <- censi[str_length(animal_id2)==3,countkin(animal_id2,Afull,0.5),by=.(YEAR,ngroup)]
setnames(kincounts,c("YEAR","ngroup"),c("Year","Group"))


#extract IDs
id_dat <- unique(all_obs[,.(FocalID,Year)])

# merge and calculate age data
id_dat <-  merge(id_dat,pedigree[,c(1, 2, 11)], by.x = "FocalID",by.y = "Focal_ID")
id_dat[,Age:=as.numeric(Year) - CS.BIRTH.SEASON]
id_dat[,CS.BIRTH.SEASON:=NULL]
id_dat <- merge(id_dat,dominance[,c("ID","YEAR","ORD_RANK")],by.x=c("FocalID","Year"),by.y=c("ID","YEAR"))

nc1 <- ncol(all_obs)
# remove animals w/ no dominance rank
all_obs <- merge(all_obs,id_dat,by=c("FocalID","Year"))
nc2 <- ncol(all_obs)
all_obs <- all_obs[,c(1:5,(nc2-2):nc2,6:nc1),with=F]


