library(lubridate)


seas_Loss_func_lead <- function(month_in_seas,vec_obs, date_obs, array_hind, date_init, lead){
  #vec_obs = time
  #array_hind = lead*init
  
  Hind_sel = array_hind[lead,which(month(date_init + lead-1) %in% month_in_seas)]
  
  ref_dates_obs = which((as.logical(date_obs %in% (date_init + lead-1))) & (month(date_obs) %in% month_in_seas))
  Sel_dates_obs = date_obs[ref_dates_obs]
  Obs_sel = vec_obs[ref_dates_obs]
  
  # return((mean((Obs_sel == 1 & Hind_sel == 0) | 
  #                (Obs_sel == 0 & Hind_sel == 1)))/mean(Obs_sel == 1 |
  #                                                        Hind_sel == 1))
  
  return(mean((Obs_sel + Hind_sel) == 1) / mean((Obs_sel + Hind_sel) >= 1))
}# end seas loss func




Week_acc_Seas_Loss_func_lead <- function(lead_start, nb_days_acc, nb_ev, month_in_seas, vec_obs, date_obs, array_hind, date_init){
  #vec_obs = time
  #array_hind = lead*init
  Nb_events_obs_wlead = numeric()
  Nb_events_hind_wlead = numeric()
  
  count <- 0
  for (INIT in 1:length(date_init)) {
    if(month(date_init[INIT]+lead_start-1) %in% month_in_seas){
      count <- count+1
      Nb_events_obs_wlead[count] <- sum(vec_obs[which(date_obs>=(date_init[INIT]+lead_start-1) & date_obs<(date_init[INIT]+lead_start+nb_days_acc-1))])
      Nb_events_hind_wlead[count] <- sum(array_hind[lead_start:(lead_start+nb_days_acc-1),INIT])
    }#end if good seas
  }#end for INIT
  
  Obs_sel = Nb_events_obs_wlead>=nb_ev
  Hind_sel = Nb_events_hind_wlead>=nb_ev
  
  return(mean((Obs_sel + Hind_sel) == 1) / mean((Obs_sel + Hind_sel) >= 1))
}# end seas loss func



Week_acc_Seas_Loss_func_lead <- function(lead_start, nb_days_acc, nb_ev, month_in_seas, vec_obs, date_obs, array_hind, date_init){
  #vec_obs = time
  #array_hind = member*lead*init
  Nb_events_obs_wlead = numeric()
  Nb_events_hind_wlead = numeric()
  
  sub_init = which(month(date_init+lead_start-1) %in% month_in_seas)
  
  for (INIT in 1:length(sub_init)) {
    Nb_events_obs_wlead[INIT] <- sum(vec_obs[which(date_obs>=(date_init[sub_init[INIT]]+lead_start-1) & date_obs<(date_init[sub_init[INIT]]+lead_start+nb_days_acc-1))])
    sum_mem <- apply(array_hind[,lead_start:(lead_start+nb_days_acc-1),sub_init[INIT]], MARGIN = c(1), FUN = sum)
    Nb_events_hind_wlead[INIT] <- median(sum_mem)
  }#end for INIT
  
  Obs_sel = Nb_events_obs_wlead>=nb_ev
  Hind_sel = Nb_events_hind_wlead>=nb_ev
  
  return(mean((Obs_sel + Hind_sel) == 1) / mean((Obs_sel + Hind_sel) >= 1))
}# end seas loss func



Spati_acc_seas_Loss_func_lead <- function(lead_time, nb_ev, month_in_seas, square_obs_neighb, date_obs, square_hind_neighb, date_init){
  
  #square_obs_neighb = lon*lat*time
  #square_hind_neighb = lon*lat*member*lead*init
  
  Nb_events_obs_square = numeric()
  Nb_events_hind_square = numeric()
  
  sub_init = which(month(date_init+lead_time-1) %in% month_in_seas)
  
  for (INIT in 1:length(sub_init)) {
    Nb_events_obs_square[INIT] <- sum(square_obs_neighb[,,which(date_obs==(date_init[sub_init[INIT]]+lead_time-1))])
    sum_mem <- apply(square_hind_neighb[,,,lead_time,sub_init[INIT]], MARGIN = c(3), FUN = sum)
    Nb_events_hind_square[INIT] <- median(sum_mem)
  }#end for INIT
  
  Obs_sel = Nb_events_obs_square>=nb_ev
  Hind_sel = Nb_events_hind_square>=nb_ev
  
  return(mean((Obs_sel + Hind_sel) == 1) / mean((Obs_sel + Hind_sel) >= 1))
}# end spatial accumulation seas loss func




Spati_acc_WR_Loss_func_lead <- function(lead_time, nb_ev, dates_WR, square_obs_neighb, date_obs, square_hind_neighb, date_init){
  
  #square_obs_neighb = lon*lat*time
  #square_hind_neighb = lon*lat*member*lead*init
  
  Nb_events_obs_square = numeric()
  Nb_events_hind_square = numeric()
  
  sub_init = which((date_init+lead_time-1) %in% dates_WR)
  
  for (INIT in 1:length(sub_init)) {
    Nb_events_obs_square[INIT] <- sum(square_obs_neighb[,,which(date_obs==(date_init[sub_init[INIT]]+lead_time-1))])
    sum_mem <- apply(square_hind_neighb[,,,lead_time,sub_init[INIT]], MARGIN = c(3), FUN = sum)
    Nb_events_hind_square[INIT] <- median(sum_mem)
  }#end for INIT
  
  Obs_sel = Nb_events_obs_square>=nb_ev
  Hind_sel = Nb_events_hind_square>=nb_ev
  
  BLS = mean((Obs_sel + Hind_sel) == 1) / mean((Obs_sel + Hind_sel) >= 1)
  
  return(list(BLS=BLS, nb_ev_obs = sum(Obs_sel), nb_ev_hind = sum(Hind_sel)))
}# end spatial accumulation seas loss func

