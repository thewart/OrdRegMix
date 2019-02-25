all_obs[,SDB:= SelfGrm + Scratch]
if (!modonkin & !modonrank) {
  all_obs[,NonConAgg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)`]
  all_obs[,NonConAgg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)`]
  all_obs[,Agg_give:=`threat:direct'n(give)` + `noncontactAgg:direct'n(give)` + `contactAgg:direct'n(give)`]
  all_obs[,Agg_rec:=`threat:direct'n(receive)` + `noncontactAgg:direct'n(receive)` + `contactAgg:direct'n(receive)`]
  all_obs[,Avoid_give:=`displace:winner?(foc)` + `avoid:winner?(foc)`]
  all_obs[,Avoid_rec:=`displace:winner?(partnr)` + `avoid:winner?(partnr)`]
  all_obs[,Groom:=GroomGIVE+GroomGET]
}

if (modonkin) {
  if (modonsex) {
    all_obs[,Agg_give_kin:=`threat:direct'n(give):samesex(TRUE):iskin(TRUE)` +
              `noncontactAgg:direct'n(give):samesex(TRUE):iskin(TRUE)` +
              `contactAgg:direct'n(give):samesex(TRUE):iskin(TRUE)`]
    all_obs[,Agg_give_nonkin:=`threat:direct'n(give):samesex(TRUE):iskin(FALSE)` +
              `noncontactAgg:direct'n(give):samesex(TRUE):iskin(FALSE)` +
              `contactAgg:direct'n(give):samesex(TRUE):iskin(FALSE)`]
    
    all_obs[,Groom_kin:=`GroomGIVE:samesex(TRUE):iskin(TRUE)` + 
              `GroomGET:samesex(TRUE):iskin(TRUE)`]
    all_obs[,Groom_nonkin:=`GroomGIVE:samesex(TRUE):iskin(FALSE)` + 
              `GroomGET:samesex(TRUE):iskin(FALSE)`]
    
    all_obs[,passcont:=`passcont:samesex(TRUE):iskin(TRUE)` + 
              `passcont:samesex(TRUE):iskin(FALSE)`]
  } else {
    all_obs[,Agg_give_kin:=`threat:direct'n(give):iskin(TRUE)` +
              `noncontactAgg:direct'n(give):iskin(TRUE)` +
              `contactAgg:direct'n(give):iskin(TRUE)`]
    all_obs[,Agg_give_nonkin:=`threat:direct'n(give):iskin(FALSE)` +
              `noncontactAgg:direct'n(give):iskin(FALSE)` +
              `contactAgg:direct'n(give):iskin(FALSE)`]
    
    all_obs[,Groom_kin:=`GroomGIVE:iskin(TRUE)` + `GroomGET:iskin(TRUE)`]
    all_obs[,Groom_nonkin:=`GroomGIVE:iskin(FALSE)` + `GroomGET:iskin(FALSE)`]
    
    all_obs[,passcont:=`passcont:iskin(TRUE)` + 
              `passcont:iskin(FALSE)`]
  }
}

if (modonrank) {
  all_obs[,Agg_give_higher:=`threat:direct'n(give):samesex(TRUE):partnerrank(higher)` +
            `noncontactAgg:direct'n(give):samesex(TRUE):partnerrank(higher)` +
            `contactAgg:direct'n(give):samesex(TRUE):partnerrank(higher)`]
  all_obs[,Agg_give_lower:=`threat:direct'n(give):samesex(TRUE):partnerrank(lower)` +
            `noncontactAgg:direct'n(give):samesex(TRUE):partnerrank(lower)` +
            `contactAgg:direct'n(give):samesex(TRUE):partnerrank(lower)`]
  all_obs[,Agg_give:=Agg_give_higher+Agg_give_lower]
  
  all_obs[,Avoid:=`avoid:winner?(partnr):samesex(TRUE):partnerrank(higher)` +
            `avoid:winner?(partnr):samesex(TRUE):partnerrank(lower)` + 
            `displace:winner?(partnr):samesex(TRUE):partnerrank(higher)` +
            `displace:winner?(partnr):samesex(TRUE):partnerrank(lower)`]
  all_obs[,Submit:=`Submit:direct'n(give):samesex(TRUE):partnerrank(higher)` +
            `Submit:direct'n(give):samesex(TRUE):partnerrank(lower)` + 
            `FearGrm:direct'n(give):samesex(TRUE):partnerrank(higher)` + 
            `FearGrm:direct'n(give):samesex(TRUE):partnerrank(lower)`]
  
  all_obs[,Groom_higher:=`GroomGIVE:samesex(TRUE):partnerrank(higher)` + 
            `GroomGET:samesex(TRUE):partnerrank(higher)`]
  all_obs[,Groom_lower:=`GroomGIVE:samesex(TRUE):partnerrank(lower)` + 
            `GroomGET:samesex(TRUE):partnerrank(lower)`]
  
  all_obs[,passcont:=`passcont:samesex(TRUE):partnerrank(higher)` + 
            `passcont:samesex(TRUE):partnerrank(lower)`]
}
