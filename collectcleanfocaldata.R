library(readxl)
source("~/code/LatentSocialPheno/parsefocaldata.R")
basepath <- "~/Dropbox/focaldata_processed/"
behpath <- "~/Dropbox/Behavioural Data for Seth/"

# read in and format data -----
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv",
                           "HH2014/Txtexports_all_processed.csv",
                           "R2015/Txtexports_all_processed.csv"))
foo <- readfocfiles(fpath,group = c("F","HH","R"),minobs = 0)
#manually fix ambiguous parnter IDs
foo[PartnerID=="4B2;40V",PartnerID:="4B2"]
foo[PartnerID=="80T;8B0",PartnerID:="8B0"]
foo[PartnerID=="30L;3A2",PartnerID:="30L"]
foo[PartnerID=="01K;0K1",PartnerID:="01K"]

F_2016 <- read_excel(paste0(behpath,"2016.GrpF_Masterfile.xlsx"), sheet = "2016.GrpF_MASTER")
F_2016 <- as.data.table(F_2016)
F_2016[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
F_2016[,Year:=as.character(Year)]
F_2016[,Year:="2016"]
F_2016 <- F_2016[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
F_2016$Group <- "F"

F2014SDB <- fread(paste0(behpath,"Masterdatafile2014_ALL_GroupF_SDB.csv"))
F2014SDB[,`Start Time`:=as.POSIXct(`Start Time`,format="%H:%M:%S")]
F2014SDB[,`Stop Time`:=as.POSIXct(`Stop Time`,format="%H:%M:%S")]
F2014SDB[,Year:=as.character(Year)]
F2014SDB[,Year:="2014"]
F2014SDB <- F2014SDB[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`, Duration, Behavior,`Behavior Modifier`,PartnerID)]

F_2014 <- fread(paste0(behpath,"Masterdatafile2014_ALL_GroupF.csv"))
F_2014[,`Start Time`:=as.POSIXct(`Start Time`,format="%H:%M:%S")]
F_2014[,`Stop Time`:=as.POSIXct(`Stop Time`,format="%H:%M:%S")]
F_2014[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,x,x_1,x_2,sep="; ")]
F_2014[,Year:=as.character(Year)]
F_2014[,Year:="2014"]
F_2014[,`Start Time`:=as.POSIXct(`Start Time`,format="%H:%M:%S")]
F_2014[,`Stop Time`:=as.POSIXct(`Stop Time`,format="%H:%M:%S")]
F_2014 <- F_2014[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`, Duration, Behavior,`Behavior Modifier`,PartnerID)]
F_2014 <- rbind(F_2014,F2014SDB)
F_2014$Group <- "F"
F_2014[PartnerID %in% c("0A10","0A11","0A12","0A13"),PartnerID:="0A9"]

V_2015 <- read_excel(paste0(behpath,"2015_GroupV_FOCALdata.xlsx"), sheet = "groupV 2015")
V_2015 <- as.data.table(V_2015)
V_2015[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
V_2015[,Year:=as.character(Year)]
V_2015[,Year:="2015"]
V_2015 <- V_2015[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
V_2015$Group <- "V"

V_2016 <- read_excel(paste0(behpath,"2016_GroupV_MASTERFILE.xlsx"), sheet = "Sheet1")
V_2016 <- as.data.table(V_2016)
V_2016[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
V_2016[,Year:=as.character(Year)]
V_2016[,Year:="2016"]
V_2016 <- V_2016[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
V_2016$Group <- "V"

KK_2015 <- read_excel(paste0(behpath,"2015_groupKK_FOCALdata.xlsx"), sheet = "grpKK_2015")
KK_2015 <- as.data.table(KK_2015)
KK_2015[,`Behavior Modifier`:= paste(`Behavior Modifier`,`Behavior Modifier 1`,X,X__1,X__2,sep="; ")]
KK_2015[,Year:=as.character(Year)]
KK_2015[,Year:="2015"]
KK_2015 <- KK_2015[,.(`Event Name`,`Observation Name`,`Focal ID`,Observer,Year,
                    `Start Time`,`Stop Time`,Duration,Behavior,`Behavior Modifier`,PartnerID)]
KK_2015$Group <- "KK"

F_2015 <- fread(paste0(behpath,"/GroupF2015_FocalData.txt"))
F_2015 <- F_2015[constrained.duration!="overtime"]
F_2015[,constrained.duration:=as.numeric(constrained.duration)]
F_2015[,behaviour.starttime:=strptime(behaviour.starttime,format="%H:%M:%S")]
setnames(F_2015,old = c("focal.id","behaviour","constrained.duration","partner.id","observation.name","year"),
         new=c("Focal ID","Event Name","Duration_Revised","PartnerID","Observation Name","Year"))
F_2015[,StopTime:=behaviour.starttime+Duration_Revised]
F_2015[,Year:=rep("2015",nrow(F_2015))]
F_2015[,Group:="F"]
F_2015 <- F_2015[,c("Event Name","Observation Name","Focal ID","Observer","Year",
                    "behaviour.starttime","StopTime","Duration_Revised","Event Name","behaviour.modifers","PartnerID","Group")]
setnames(F_2015,colnames(KK_2015))

all_files <- list(F_2014, F_2015, KK_2015, V_2015, V_2016, F_2016)
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

#get rid of animals w/ suspiciously few observations
obscount <- bdat[,length(unique(Observation)),by=c("FocalID","Year","Group")]
obskeep <- obscount[,.(FocalID,V1,V1>(mean(V1)-2*sd(V1))),by=c("Year","Group")][V3==T]
bdat <- bdat[bdat[,paste0(FocalID,Year,Group)] %in% obskeep[,paste0(FocalID,Year,Group)]]
