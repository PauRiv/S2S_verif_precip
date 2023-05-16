library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)


# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""

# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()))
perc=95
alpha = 0.05
boot = 1000

ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO" 
# ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA" 

nb_days_acc = 7

# Useful function --------------------------------------------------------

get_CI_paralelized_counts <- function(nb_days_acc, nb_ev, month_in_seas,vec_obs, bootsize, date_obs, date_init){
  
  Nb_events_obs_wlead = numeric() # count number of events in nb_days_acc days, for observation
  
  sub_init = which(month(date_init) %in% month_in_seas)
  
  for (INIT in 1:length(sub_init)) {
    Nb_events_obs_wlead[INIT] <- sum(vec_obs[which(date_obs>=(date_init[sub_init[INIT]]) & date_obs<(date_init[sub_init[INIT]]+nb_days_acc))])
  }#end for INIT
  
  Obs_sel = Nb_events_obs_wlead>=nb_ev # change to binary data
  
  RES_list = foreach (boot_member = 1:bootsize) %dopar% { #bootstrap, random drawn of number of events during _acc days, for the good season
    
    Nb_events_hind_wlead = numeric()
    #randomly choose length(sub_init) fake init days:
    sample_date_obs = sample(which((month(date_obs) %in% month_in_seas) & ((date_obs+nb_days_acc-1) < (max(date_init)+nb_days_acc))), #carreful with end of time series
                             size = length(sub_init), replace = F)
    
    for (INIT in 1:length(sub_init)) {#count the events during nb_days_acc days after the fake init
      Nb_events_hind_wlead[INIT] <- sum(vec_obs[sample_date_obs[INIT]:(sample_date_obs[INIT]+nb_days_acc-1)])
    }#end for INIT
    
    
    Hind_sel = Nb_events_hind_wlead>=nb_ev # change to binary data
    
    return(mean((Obs_sel + Hind_sel) == 1) / mean((Obs_sel + Hind_sel) >= 1)) # get boot BLF, for one member
    
  }#end foreach boot_member
  
  RES <- numeric()
  
  for (boot_member in 1:bootsize) {
    if(length(RES_list[[boot_member]])==length(RES[boot_member])){
      RES[boot_member] <- RES_list[[boot_member]]
    } else {
      RES[boot_member] <- NA
    }#end if no pb paral
  }#end for boot_m
  
  return(RES) #vector of size bootsize
} # end function try_get_CI




# Compute CI and first time member in CI ---------------------------------------------
CI_climato_5th_perc_1ev <- array(dim=c(130,85))
CI_climato_5th_perc_2ev <- array(dim=c(130,85))
CI_climato_5th_perc_3ev <- array(dim=c(130,85))

for (lon in 1:130) {
  T1<-Sys.time()
  nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_",perc,"thperc_season_dependant.nc"))
  Exceed_obs = ncvar_get(nc_obs, varid = paste0("exceed",perc), start = c(lon,1,1), count = c(1,-1,-1))
  time_obs =  ncvar_get(nc_obs, varid = "time")
  nc_close(nc_obs)
  
  time_obs_days = as.Date(time_obs, origin = "1970-01-01")
  # hindcast file names, one file per init
  File_names_exceed <- list.files(hind_exceed_folder)
  init_dates <- substr(File_names_exceed, 35, 44)
  init_days = as.Date(init_dates)
  
  set.seed(2022)
  
  for (lat in 1:85) {
    CI_BLF_1=get_CI_paralelized_counts(nb_days_acc = nb_days_acc, nb_ev=1,month_in_seas=ext_seas, vec_obs = Exceed_obs[lat,],
                                       bootsize = boot, date_obs = time_obs_days, date_init = init_days)
    CI_climato_5th_perc_1ev[lon,lat] <- quantile(CI_BLF_1, p=alpha, na.rm=T)
    
    CI_BLF_2=get_CI_paralelized_counts(nb_days_acc = nb_days_acc, nb_ev=2,month_in_seas=ext_seas, vec_obs = Exceed_obs[lat,],
                                       bootsize = boot, date_obs = time_obs_days, date_init = init_days)
    CI_climato_5th_perc_2ev[lon,lat] <- quantile(CI_BLF_2, p=alpha, na.rm=T)
    
    CI_BLF_3=get_CI_paralelized_counts(nb_days_acc = nb_days_acc, nb_ev=3,month_in_seas=ext_seas, vec_obs = Exceed_obs[lat,],
                                       bootsize = boot, date_obs = time_obs_days, date_init = init_days)
    CI_climato_5th_perc_3ev[lon,lat] <- quantile(CI_BLF_3, p=alpha, na.rm=T)
  }
  T2 <- Sys.time()
  
  print(paste("done with lon",lon,"out of 130"))
  print(difftime(T2,T1))
}#end for lon



# Save in netcdf file ---------------------------------------------------

nc_hind = nc_open(hind_exceed_folder)
lon_dim =  ncvar_get(nc_hind, varid = "lon")
unit_lon = nc_hind$dim$lon$units
lat_dim =  ncvar_get(nc_hind, varid = "lat")
unit_lat = nc_hind$dim$lat$units
nc_close(nc_hind)

gc()

fillvalue <- -9999
latdim <- ncdim_def("lat",unit_lat,lat_dim)
londim <- ncdim_def("lon",unit_lon,lon_dim)

BLS_clim_lowerbound_1ev <- ncvar_def(name = "BLS_clim_lowerbound_1ev", units = "no", dim = list(londim,latdim), missval = fillvalue,
                                     longname = "BLS_clim_lowerbound_5thperc_accumu_1ev")
BLS_clim_lowerbound_2ev <- ncvar_def(name = "BLS_clim_lowerbound_2ev", units = "no", dim = list(londim,latdim), missval = fillvalue,
                                     longname = "BLS_clim_lowerbound_5thperc_accumu_2ev")
BLS_clim_lowerbound_3ev <- ncvar_def(name = "BLS_clim_lowerbound_3ev", units = "no", dim = list(londim,latdim), missval = fillvalue,
                                     longname = "BLS_clim_lowerbound_5thperc_accumu_3ev")

ncout <- nc_create(filename = paste0(BLI_folder,"CI_",100*(alpha),"thlowerbound_tempACC_",nb_days_acc,
                                     "days_seasonal_BLF_perlead_andmember_",perc,"th_perc_",name_seas,".nc"),
                   vars = list(BLS_clim_lowerbound_1ev,BLS_clim_lowerbound_2ev,BLS_clim_lowerbound_3ev), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLS_clim_lowerbound_1ev, vals = CI_climato_5th_perc_1ev)
ncvar_put(nc = ncout, varid = BLS_clim_lowerbound_2ev, vals = CI_climato_5th_perc_2ev)
ncvar_put(nc = ncout, varid = BLS_clim_lowerbound_3ev, vals = CI_climato_5th_perc_3ev)

nc_close(ncout)
