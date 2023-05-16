library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)

#to adapt accordingly:
exceed_folder = ""
hind_exceed_folder = ""
Brier_folder = ""

# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2.3)) #27 cores
perc=95
ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO" # h 27 cores
# ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA" #h 27 cores





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

BLF_array <- array(dim=c(130,85))

for (LON in 1:130) {
  RES <- foreach (LAT = 1:85) %dopar% {
    
    
    mean_value = mean(Exceed_obs[LON,LAT,which(month(time_obs_days) %in% ext_seas)])
    
    BLF_1 = mean((Exceed_obs[LON,LAT,which(month(time_obs_days) %in% ext_seas)] - mean_value)^2)
    
    return(BLF_1)
    
  }#end for LAT
  
  for (LAT in 1:85) {
    BLF_array[LON,LAT] <- RES[[LAT]]
  }#end for LAT
  # print(LON)
  
}#end foreach LON #35 sec per LON


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

BLF_var <- ncvar_def(name = "Brier", units = "no", dim = list(londim,latdim), missval = fillvalue,
                     longname = "Brier_score_local")

ncout <- nc_create(filename = paste0(Brier_folder,"Brier_score_climato_perlead_",perc,"th_perc_",name_seas,"_CONSTCLIMATO.nc"),
                   vars = list(BLF_var), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLF_var, vals = BLF_array)

nc_close(ncout)
