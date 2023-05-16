library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)



# to adapt accordingly:
exceed_folder = ""
hind_exceed_folder = ""
Brier_folder = ""


# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2))
perc=95
ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO" 
# ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA" 


# Useful function --------------------------------------------------------

Brier_seas_lead <- function(lead_time, month_in_seas, obs_vec, date_obs, hind_array, date_init){
  
  #square_obs_neighb = time
  #hind_array = member*lead*init
  
  sub_init = which(month(date_init+lead_time-1) %in% month_in_seas)
  
  Obs_sel = numeric()
  Hind_sel = numeric()
  
  for (INIT in 1:length(sub_init)) {
    Obs_sel[INIT] = obs_vec[which(date_obs==(date_init[sub_init[INIT]]+lead_time-1))]
    Hind_sel[INIT] <- mean(hind_array[,lead_time,sub_init[INIT]])
  }#end for INIT
  
  return(mean((Obs_sel-Hind_sel)^2))
}# end spatial accumulation Brier



# Load data ---------------------------------------------------------------

nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_",perc,"thperc_season_dependant.nc"))
Exceed_obs = ncvar_get(nc_obs, varid = paste0("exceed",perc))
time_obs =  ncvar_get(nc_obs, varid = "time")
nc_close(nc_obs)

time_obs_days = as.Date(time_obs, origin = "1970-01-01")

File_names_exceed <- list.files(hind_exceed_folder) # hindcast file names, one file per init
init_dates <- substr(File_names_exceed, 35, 44)
init_days = as.Date(init_dates)

gc()



# Compute BLF: loop on longitudes -----------------------------------------

BLF_array <- array(dim=c(130,85,46))

for (LON in 1:130) {
  Exceed_hind_LON <- array(dim = c(85, 11, 46, 2080))
  
  Exceed_list <- foreach (INI = 1:length(init_days)) %dopar% {
    Init = init_dates[INI]
    nc_hind = nc_open(paste0(hind_exceed_folder,"Exceedences_Hincast_",perc ,"thperc_init_",Init,".nc"))
    A = ncvar_get(nc_hind, varid = paste0("exceed_perc",perc), start = c(LON,1,1,1), count = c(1,-1,-1,-1)) #get only one longitude
    nc_close(nc_hind)
    gc()
    return(A)
  }#end for INI
  
  
  for (INI in 1:length(init_days)) {
    Exceed_hind_LON[,,,INI] <- Exceed_list[[INI]]
  }#end for INI
  
  RES <- foreach (LAT = 1:85) %dopar% {
    BLF_1 = apply(X=as.matrix(1:46), MARGIN=c(1), FUN=Brier_seas_lead, month_in_seas=ext_seas,
                  obs_vec = Exceed_obs[LON,LAT,], hind_array = Exceed_hind_LON[LAT,,,],
                  date_obs = time_obs_days, date_init = init_days)
    
    return(BLF_1)
    
  }#end for LAT
  
  for (LAT in 1:85) {
    BLF_array[LON,LAT,] <- RES[[LAT]]
  }#end for LAT
  
  
}#end foreach LON #4 min per LON


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

BLF_var <- ncvar_def(name = "Brier", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                     longname = "Brier_score_local")

ncout <- nc_create(filename = paste0(Brier_folder,"Brier_score_perlead_",perc,"th_perc_",name_seas,".nc"),
                   vars = list(BLF_var), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLF_var, vals = BLF_array)

nc_close(ncout)
