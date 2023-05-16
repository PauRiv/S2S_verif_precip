library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)

# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2.8))
perc=95
ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO"
#ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA"

Member="med"
#Member="MIN"
#Member="MAX"


# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""


# Useful functions --------------------------------------------------------

source("Binary_loss_index_functions.R")

# Get data -------------------------------------------------------- 

nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_", perc,"thperc_season_dependant.nc"))
Exceed_obs = ncvar_get(nc_obs, varid = paste0("exceed",perc))
time_obs =  ncvar_get(nc_obs, varid = "time")
nc_close(nc_obs)

time_obs_days = as.Date(time_obs, origin = "1970-01-01")

# hindcast file names, one file per init
File_names_exceed <- list.files(hind_exceed_folder)

init_dates <- substr(File_names_exceed, 42, 51)
init_days = as.Date(init_dates)

gc()


BLF_list <- foreach(LON=1:130) %dopar% {
  Exceed_hind_LON <- array(dim = c(85, 46, 2080))
  
  for (INI in 1:length(init_dates)) {
    Init = init_dates[INI]
    if(Member=="med"){nc_hind = nc_open(paste0(hind_median_exceed_folder,"Memb_Median_Exceed_Hincast_",perc ,"thperc_init_",Init,".nc"))}
    if(Member=="MIN"){nc_hind = nc_open(paste0(hind_min_exceed_folder,"Memb_MIN_Exceed_Hincast_95thperc_init_",Init,".nc"))}
    if(Member=="MAX"){nc_hind = nc_open(paste0(hind_min_exceed_folder,"mb_MAX_Exceed_Hincast_95thperc_init_", Init,".nc"))}
    
    Exceed_hind_LON[,,INI] = ncvar_get(nc_hind, varid = paste0("sel_memb_exceed",perc),
                                       start = c(LON,1,1), count = c(1,-1,-1)) #get only one longitude
    nc_close(nc_hind)
    gc()
  }#end for INI
  
  
  BLF <- array(dim = c(85,46))
  for (LAT in 1:85) {
    BLF[LAT,] <- apply(X=as.matrix(1:46), MARGIN = c(1), FUN = seas_Loss_func_lead,
                       vec_obs = Exceed_obs[LON,LAT,], array_hind = Exceed_hind_LON[LAT,,],
                       date_obs = time_obs_days, date_init = init_days, month_in_seas=ext_seas)
  }#end for LAT
  
  return(BLF)
}#end foreach LON #8 min per LON

BLF_array <- array(dim=c(130,85,46))

for (LON in 1:130) {
  BLF_array[LON,,] <- BLF_list[[LON]]
}#end for LON

# 1h10min



# Save in netcdf file ---------------------------------------------------

nc_hind = nc_open(paste0(hind_data_path,"Exceedences_Hincast_",perc,"thperc_init_2001-01-04.nc"))
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

BLF_var <- ncvar_def(name = "BLF", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                     longname = "Binary_Loss_Function")


if(Member=="med"){ncout <- nc_create(filename = paste0(BLI_folder,"seasonal_BLF_perlead_medmemb_", perc,"th_perc_",name_seas,".nc"),
                                     vars = list(BLF_var), force_v4=TRUE)}
if(Member=="MIN"){ncout <- nc_create(filename = paste0(BLI_folder,"seasonal_BLF_perlead_MINmemb_", perc,"th_perc_",name_seas,".nc"),
                                     vars = list(BLF_var), force_v4=TRUE)}
if(Member=="MAX"){ncout <- nc_create(filename = paste0(BLI_folder,"seasonal_BLF_perlead_MAXmemb_", perc,"th_perc_",name_seas,".nc"),
                                     vars = list(BLF_var), force_v4=TRUE)}

ncvar_put(nc = ncout, varid = BLF_var, vals = BLF_array)

nc_close(ncout)
