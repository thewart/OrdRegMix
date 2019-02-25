source("~/code/OrdRegMix/recodebehaviors.R")
if (modonkin) {
  if (modonsex) {
    behaviors <- c("SDB","Feed","Travel","passcont",
                   "Agg_give_kin","Agg_give_nonkin",
                   "Groom_kin","Groom_nonkin",
                   "Approach:initiate(focal):samesex(TRUE):iskin(TRUE)",
                   "Approach:initiate(focal):samesex(TRUE):iskin(FALSE)")
  } else 
    behaviors <- c("SDB","Feed","Travel","passcont",
                   "Agg_give_kin","Agg_give_nonkin",
                   "Groom_kin","Groom_nonkin",
                   "Approach:initiate(focal):iskin(TRUE)",
                   "Approach:initiate(focal):iskin(FALSE)")
  
  used_obs <- all_obs[Group=="F" & (paste0(FocalID,Year) %in% id_dat[SEX=="f" & kin_samesex>min(kin_samesex),paste0(FocalID,Year)])] %>% copy()
  cov_dat <- id_dat[SEX=="f" & Group=="F" & kin_samesex>min(kin_samesex)]
  leftside <- "Year + kin_samesex"
}

if (modonrank) {
  behaviors <- c("SDB","Feed","Travel","passcont","Agg_give","Avoid","Submit","Groom_higher","Groom_lower",
                 "Approach:initiate(focal):samesex(TRUE):partnerrank(higher)",
                 "Approach:initiate(focal):samesex(TRUE):partnerrank(lower)")
  
  used_obs <- all_obs[Group=="F" & (paste0(FocalID,Year) %in% id_dat[SEX=="f" & ORD_RANK<"H",paste0(FocalID,Year)])] %>% copy()
  cov_dat <- id_dat[SEX=="f" & Group=="F" & ORD_RANK<"H"]
  leftside <- "Year + ORD_RANK"
}

if (!modonkin & !modonrank) {
  behaviors = c("SDB","Feed","Travel","passcont","Groom","Agg_give","Approach:initiate(focal)")
  used_obs <- all_obs[Group=="F" & (FocalID %in% id_dat[SEX=="f",FocalID])] %>% copy()
  cov_dat <- id_dat[SEX=="f" & Group=="F"]
  leftside <- "Year"
}

d <- length(behaviors)
leftside_ranef <- "(1|FocalID)"
used_obs <- used_obs[,c(names(used_obs)[1:5],behaviors),with=F]

# loop quantile function for all behaviors ----
Y <- used_obs[,behaviors,with=F]
Y <- lapply(Y, function(x) cut(x,c(0,quantile(x,c(0.5,1.0))+0.5),
                               include.lowest = T,ordered_result = T))
behav_bins <- sapply(Y,levels)
Y <- lapply(Y,function(x) as.numeric(x) - 1)
used_obs[,(behaviors):=Y]

