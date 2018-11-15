behaviors <- c("SDB","Agg_kin","Agg_nonkin",
               "GroomGIVE:samesex(TRUE):iskin(TRUE)",
               "GroomGIVE:samesex(TRUE):iskin(FALSE)",
               "passcont:samesex(TRUE):iskin(TRUE)",
               "passcont:samesex(TRUE):iskin(FALSE)",
               "Approach:initiate(focal):samesex(TRUE):iskin(TRUE)",
               "Approach:initiate(focal):samesex(TRUE):iskin(FALSE)")
# behaviors = c("SDB","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)",
#               "Approach:initiate(partner)","NonConAgg_give","NonConAgg_rec",
#               "contactAgg:direct'n(give)","contactAgg:direct'n(receive)")
d <- length(behaviors)
Y <- used_obs[,behaviors,with=F] %>% as.matrix()

used_obs <- all_obs[Group=="F" & paste0(FocalID,Year) %in% id_dat[SEX=="f" & kin_samesex>min(kin_samesex),paste0(FocalID,Year)]]
cov_dat <- id_dat[SEX=="f" & Group=="F" & kin_samesex>min(kin_samesex)]
leftside <- "Year + kin_samesex"
leftside_ranef <- "(1|FocalID)"
