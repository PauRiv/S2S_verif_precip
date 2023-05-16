library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)

# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2)) 
perc=95
ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO" #
# ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA" 
nb_days_acc = 7

# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""

source("Binary_loss_index_functions.R")

# load observations ------------------------------------------------------

nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_", perc,"thperc_season_dependant.nc"))
Exceed_obs = ncvar_get(nc_obs, varid = paste0("exceed",perc))
time_obs =  ncvar_get(nc_obs, varid = "time")
nc_close(nc_obs)

time_obs_days = as.Date(time_obs, origin = "1970-01-01")

File_names_exceed <- list.files(hind_exceed_folder) # hindcast file names, one file per init
init_dates <- substr(File_names_exceed, 35, 44)
init_days = as.Date(init_dates)

gc()




# Compute BLF: loop on longitudes -----------------------------------------



BLF_list <- foreach(LON=1:130) %dopar% {
  Exceed_hind_LON <- array(dim = c(85, 11, 46, 2080))
  
  for (INI in 1:length(init_dates)) {
    Init = init_dates[INI]
    nc_hind = nc_open(paste0("/scratch3/pauline/S2S_precip_clustering/Data_05_Europe/Exceedances/perc",
                             perc,"/Exceedences_Hincast_",perc ,"thperc_init_",Init,".nc"))
    Exceed_hind_LON[,,,INI] = ncvar_get(nc_hind, varid = paste0("exceed_perc",perc),
                                        start = c(LON,1,1,1), count = c(1,-1,-1,-1)) #get only one longitude
    nc_close(nc_hind)
    gc()
  }#end for INI
  
  
  
  BLF_1 <- array(dim = c(85,6))
  BLF_2 <- array(dim = c(85,6))
  BLF_3 <- array(dim = c(85,6))
  
  for (LAT in 1:85) {
    BLF_1[LAT,] <- apply(X=as.matrix(seq(1,42,7)), MARGIN=c(1), FUN=Week_acc_Seas_Loss_func_lead,
                         nb_days_acc = nb_days_acc, nb_ev=1, month_in_seas=ext_seas, vec_obs = Exceed_obs[LON,LAT,],
                         array_hind = Exceed_hind_LON[LAT,,,], date_obs = time_obs_days, date_init = init_days)
    
    BLF_2[LAT,] <- apply(X=as.matrix(seq(1,42,7)), MARGIN=c(1), FUN=Week_acc_Seas_Loss_func_lead,
                         nb_days_acc = nb_days_acc, nb_ev=2, month_in_seas=ext_seas, vec_obs = Exceed_obs[LON,LAT,],
                         array_hind = Exceed_hind_LON[LAT,,,], date_obs = time_obs_days, date_init = init_days)
    
    BLF_3[LAT,] <- apply(X=as.matrix(seq(1,42,7)), MARGIN=c(1), FUN=Week_acc_Seas_Loss_func_lead,
                         nb_days_acc = nb_days_acc, nb_ev=3, month_in_seas=ext_seas, vec_obs = Exceed_obs[LON,LAT,],
                         array_hind = Exceed_hind_LON[LAT,,,], date_obs = time_obs_days, date_init = init_days)
    
  }#end for LAT
  
  return(list(BLF_1=BLF_1,BLF_2=BLF_2,BLF_3=BLF_3))
}#end foreach LON #7 min per LON

BLF_array_1 <- array(dim=c(130,85,6))
BLF_array_2 <- array(dim=c(130,85,6))
BLF_array_3 <- array(dim=c(130,85,6))

for (LON in 1:130) {
  BLF_array_1[LON,,] <- BLF_list[[LON]]$BLF_1
  BLF_array_2[LON,,] <- BLF_list[[LON]]$BLF_2
  BLF_array_3[LON,,] <- BLF_list[[LON]]$BLF_3
}#end for LON

# 1h10min



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
leaddim <- ncdim_def("leadweek","week",1:6)

BLF_var_1 <- ncvar_def(name = "BLF_1ev", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_1eventormore")
BLF_var_2 <- ncvar_def(name = "BLF_2ev", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_2eventsormore")
BLF_var_3 <- ncvar_def(name = "BLF_3ev", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_3eventsormore")

ncout <- nc_create(filename = paste0(BLI_folder,"ACC_",nb_days_acc,
                                     "days_seasonal_BLF_perlead_",perc,"th_perc_",name_seas,".nc"),
                   vars = list(BLF_var_1,BLF_var_2,BLF_var_3), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLF_var_1, vals = BLF_array_1)
ncvar_put(nc = ncout, varid = BLF_var_2, vals = BLF_array_2)
ncvar_put(nc = ncout, varid = BLF_var_3, vals = BLF_array_3)

nc_close(ncout)
