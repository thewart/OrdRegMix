if (modonkin) {
  behaviors <- c("SDB","Agg_kin","Agg_nonkin",
               "GroomGIVE:samesex(TRUE):iskin(TRUE)",
               "GroomGIVE:samesex(TRUE):iskin(FALSE)",
               "passcont:samesex(TRUE):iskin(TRUE)",
               "passcont:samesex(TRUE):iskin(FALSE)",
               "Approach:initiate(focal):samesex(TRUE):iskin(TRUE)",
               "Approach:initiate(focal):samesex(TRUE):iskin(FALSE)")
  
  used_obs <- all_obs[Group=="F" & (paste0(FocalID,Year) %in% id_dat[SEX=="f" & kin_samesex>min(kin_samesex),paste0(FocalID,Year)])] %>% copy()
  cov_dat <- id_dat[SEX=="f" & Group=="F" & kin_samesex>min(kin_samesex)]
  leftside <- "Year + kin_samesex"
}

if (!modonkin & !modonrank) {
  behaviors = c("SDB","Feed","Travel","GroomGIVE","GroomGET","passcont",
                "Approach:initiate(focal)","NonConAgg_give","contactAgg:direct'n(give)")
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

