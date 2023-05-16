library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)




# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2.5))
perc=95
alpha = 0.05
boot = 500

# ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO" # 31h
ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA" #30h


# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""

# Useful functions --------------------------------------------------------

get_CI_paralelized_local <- function(nb_days_acc, nb_ev, month_in_seas,vec_obs, bootsize, leadsize, date_obs, date_init){
  
  RES_list = foreach (boot_member = 1:bootsize) %dopar% {
    
    
    extr_justdata_seas_obs = vec_obs[which(month(date_obs) %in% month_in_seas)]
    
    blf_clim_lead <- numeric()
    
    for (lead in 1:leadsize) {
      ref_dates_obs = which((as.logical(date_obs %in% (date_init + lead-1))) & (month(date_obs) %in% month_in_seas))
      Sel_dates_obs = date_obs[ref_dates_obs]
      Obs_sel = vec_obs[ref_dates_obs]
      
      
      extr_justdata_seas_obs = vec_obs[which(month(date_obs + lead-1) %in% month_in_seas)]
      
      random_obs=sample(extr_justdata_seas_obs, size = length(Obs_sel), replace = F)
      
      blf_clim_lead[lead] = mean((Obs_sel + random_obs) == 1) / mean((Obs_sel + random_obs) >= 1)
    }#end for lead
    return(blf_clim_lead)
  }#end foreach boot_m
  
  RES <- matrix(nrow = bootsize, ncol = leadsize)
  
  for (boot_member in 1:bootsize) {
    if(length(RES_list[[boot_member]]) == length(RES[boot_member,])) {
      RES[boot_member,] <- RES_list[[boot_member]]
    }#end if pb parall?
    
  }#end for boot_m
  
  return(RES)
} # end function try_get_CI




# Compute CI and first time member in CI ---------------------------------------------
CI_BLF <- array(dim=c(130,85,boot,46))
CI_climato_5th_perc <- array(dim=c(130,85,46))

for (lon in 1:130) {
  
  nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_", perc,"thperc_season_dependant.nc"))
  Exceed_obs = ncvar_get(nc_obs, varid = paste0("exceed",perc), start = c(lon,1,1), count = c(1,-1,-1))
  time_obs =  ncvar_get(nc_obs, varid = "time")
  nc_close(nc_obs)
  
  time_obs_days = as.Date(time_obs, origin = "1970-01-01")
  # hindcast file names, one file per init
  File_names_exceed <- list.files(hind_exceed_folder)
  init_dates <- substr(File_names_exceed, 35, 44)
  init_days = as.Date(init_dates)
  Exceed_hind <- array(dim = c(85, 11, 46, 2080))
  
  Exceed_hind_list <- foreach (INI = 1:2080) %dopar% {
    Init = init_dates[INI]
    nc_hind = nc_open(paste0(hind_exceed_folder,"Exceedences_Hincast_",perc ,"thperc_init_",Init,".nc"))
    A = ncvar_get(nc_hind, varid = paste0("exceed_perc",perc),
                  start = c(lon,1,1,1), count = c(1,-1,-1,-1))
    nc_close(nc_hind)
    gc()
    return(A)
  }#end foreach INI
  
  for (INI in 1:2080) {
    Exceed_hind[,,,INI] <- Exceed_hind_list[[INI]]
  }#end for INI
  
  set.seed(2022)
  
  for (lat in 1:85) {
    CI_BLF[lon,lat,,]=get_CI_paralelized_local(month_in_seas=ext_seas, vec_obs = Exceed_obs[lat,], bootsize = boot,
                                               leadsize = 46, date_obs = time_obs_days, date_init = init_days)
    CI_climato_5th_perc[lon,lat,] <- apply(CI_BLF[lon,lat,,], MARGIN = c(2), FUN = quantile, p=alpha, na.rm=T)
  }#end for lat
  
}#end for lon




# Save in netcdf file ---------------------------------------------------

nc_hind = nc_open(paste0(hind_exceed_folder,"Exceedences_Hincast_",perc,"thperc_init_2001-01-04.nc"))
lon_dim =  ncvar_get(nc_hind, varid = "lon")
unit_lon = nc_hind$dim$lon$units
lat_dim =  ncvar_get(nc_hind, varid = "lat")
unit_lat = nc_hind$dim$lat$units
nc_close(nc_hind)

gc()

fillvalue <- -9999
latdim <- ncdim_def("lat",unit_lat,lat_dim)
londim <- ncdim_def("lon",unit_lon,lon_dim)
leaddim <- ncdim_def("leadtime","days",1:46)
bootdim <- ncdim_def("boot","boot",1:boot)

BLS_clim_boot <- ncvar_def(name = "BLS_clim_boot", units = "no", dim = list(londim,latdim,bootdim,leaddim), missval = fillvalue,
                           longname = "Binary_Loss_Function_climato_boot")
BLS_clim_lowerbound <- ncvar_def(name = "BLS_clim_lowerbound", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                                 longname = "BLS_clim_lowerbound_5thperc")

ncout <- nc_create(filename = paste0(BLI_folder,"BLF_seas_",name_seas ,"_confidenceinterval_climato_with_R_", perc,"th_perc.nc"),
                   vars = list(BLS_clim_boot, BLS_clim_lowerbound), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLS_clim_boot, vals = CI_BLF)
ncvar_put(nc = ncout, varid = BLS_clim_lowerbound, vals = CI_climato_5th_perc)

nc_close(ncout)
