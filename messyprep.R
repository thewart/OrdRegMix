library(readxl)

mainpath <- "~/code/brianzhong/Plattlab/"
source(paste0(mainpath,"ReadXL.R"), echo=TRUE)
source("~/code/LatentSocialPheno/parsefocaldata.R")
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv",
                           "HH2014/Txtexports_all_processed.csv",
                           "R2015/Txtexports_all_processed.csv"))
foo <- readfocfiles(fpath,group = c("F","HH","R"))

F_2016 <- as.data.table(F_2016)
F_2016[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
F_2016[,Year:=as.character(Year)]
F_2016[,Year:="2016"]
F_2016 <- F_2016[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
F_2016$Group <- "F"

# F_2014 <- as.data.table(F_2014)
# F_2014[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,x,x__1,x__2,sep="; ")]
# F_2014[,Year:=as.character(Year)]
# F_2014[,Year:="2014"]
# F_2014 <- F_2014[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
#                     `Start Time`,`Stop Time`, Duration, Behavior,`Behavior Modifier`,PartnerID)]
# F_2014$Group <- "F"

V_2015 <- as.data.table(V_2015)
V_2015[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
V_2015[,Year:=as.character(Year)]
V_2015[,Year:="2015"]
V_2015 <- V_2015[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
V_2015$Group <- "V"

V_2016 <- as.data.table(V_2016)
V_2016[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
V_2016[,Year:=as.character(Year)]
V_2016[,Year:="2016"]
V_2016 <- V_2016[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
V_2016$Group <- "V"

KK_2015 <- as.data.table(KK_2015)
KK_2015[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
KK_2015[,Year:=as.character(Year)]
KK_2015[,Year:="2015"]
KK_2015 <- KK_2015[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
KK_2015$Group <- "KK"

all_files <- list(KK_2015, V_2015, V_2016, F_2016)
juh <- do.call(rbind,all_files)
juh[,`Behavior Modifier`:=str_replace_all(`Behavior Modifier`,"NA;*","")]
setnames(juh,c("EventName","Observation","FocalID","Observer","Year","StartTime",
               "StopTime","Duration","Behavior","BehaviorModifier","PartnerID","Group"))
juh <- juh[!is.na(Observation) & !is.na(StartTime) & !is.na(FocalID)]
juh[,StartTime:=format(StartTime,format="%H:%M:%S") %>% strptime(format="%H:%M:%S")]
juh[,StopTime:=format(StopTime,format="%H:%M:%S") %>% strptime(format="%H:%M:%S")]
juh[,RelativeEventTime:=as.numeric(StartTime)-min(as.numeric(StartTime)),by=Observation]
juh[,Duration:=as.numeric(StopTime)-as.numeric(StartTime),by=Observation]
juh <- juh[Duration>=0]

juh <- juh[!(Observation %in% juh[RelativeEventTime>1e4,unique(Observation)])]
juh <- juh[RelativeEventTime<630]

foo[File=="JosueWeek3_2013 - Episode Selection 010.txt",Observer:="JN"]
foo[str_length(Observer)==1 | Observer=="04",Observer:="JG"]
foo <- foo[,which(colnames(foo) %in% colnames(juh)),with=F]
setcolorder(juh,colnames(foo))
foo[,StartTime:=strptime(StartTime,format="%H:%M:%S")]
foo[,StopTime:=strptime(StopTime,format="%H:%M:%S")]
bdat <- rbind(foo,juh)
setkey(bdat,"FocalID","Observation","RelativeEventTime")

bdat[,Observer:=toupper(Observer)]
bdat[Observer %in% c("GP2","GV"),Observer:="GP"]
bdat[Observer=="JN5",Observer:="JN"]
bdat <- bdat[!(str_detect(Behavior,"revised") | str_detect(Behavior,"overtime"))]
bdat[,FocalID:=toupper(FocalID)]
bdat[,PartnerID:=toupper(PartnerID)]

ptetho <- defaultpoint2()[type!="misc" & !(behavior%in%c("Vigilnce","PsCnTerm","GrmTerm","GrmPrsnt"))]
stetho <- defaultstate2()[type!="misc" & state!="Corral"]
all_obs <- collectfocal(bdat,ptetho,stetho)

obscount <- all_obs[,length(Observation),by=c("FocalID","Year","Group")]
obskeep <- obscount[,.(FocalID,V1,V1>(mean(V1)-2*sd(V1))),by=c("Year","Group")][V3==T]
all_obs <- all_obs[all_obs[,paste0(FocalID,Year,Group)] %in% obskeep[,paste0(FocalID,Year,Group)]]

all_obs[,NonConAgg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)`]
all_obs[,NonConAgg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)`]
all_obs[,Avoid_give:=`displace:winner?(foc)` + `avoid:winner?(foc)`]
all_obs[,Avoid_rec:=`displace:winner?(partnr)` + `avoid:winner?(partnr)`]
all_obs[,SDB:= SelfGrm + Scratch]


# loop quantile function for all behaviors
blevels <- list()
for (i in 6:ncol(all_obs)) {
  all_obs[[i]] <- cut(all_obs[[i]],c(0,0.5,quantile(all_obs[[i]][all_obs[[i]]>0],c(0.5,1.0))+0.5) %>% unique(),include.lowest = T)
  blevels[[i-5]] <- levels(all_obs[[i]])
  all_obs[[i]] <- as.numeric(all_obs[[i]])
}

#extract IDs
id_dat <- unique(all_obs[,.(FocalID,Year)])

# merge and calculate age data
dominance$YEAR <- as.character(dominance$YEAR)
id_dat <-  merge(id_dat,pedigree[,c(1, 2, 11)], by.x = "FocalID",by.y = "Focal_ID")
id_dat[,Age:=as.numeric(Year) - CS.BIRTH.SEASON]
id_dat[,CS.BIRTH.SEASON:=NULL]
id_dat <- merge(id_dat,dominance[,c("ID","YEAR","ORD_RANK")],by.x=c("FocalID","Year"),by.y=c("ID","YEAR"))
id_dat[,ORD_RANK:=ordered(ORD_RANK,levels=c("L","M","H"))]

nc1 <- ncol(all_obs)
# remove animals w/ no dominance rank
all_obs <- merge(all_obs,id_dat,by=c("FocalID","Year"))
nc2 <- ncol(all_obs)
all_obs <- all_obs[,c(1:5,(nc2-2):nc2,6:nc1),with=F]
