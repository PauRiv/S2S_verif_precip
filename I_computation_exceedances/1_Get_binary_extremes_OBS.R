library(lubridate);library(ncdf4)
library(foreach);library(iterators);library(parallel);library(doParallel)

if(detectCores()>7){ncores=4}else(ncores=1)

registerDoParallel(cores=ncores) #one core per season

# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""


nc_obs=nc_open(obs_data_path)
lon=ncvar_get(nc_obs, varid = "lon")
lon=lon[1:130]
unit_lon = nc_obs$dim$lon$units
lat=ncvar_get(nc_obs, varid = "lat")
unit_lat = nc_obs$dim$lat$units
precip_obs = ncvar_get(nc_obs, varid = "pr")
precip_obs <- precip_obs[1:130,,]  #adapt to hindcast
# precip_obs <- precip_obs[,length(lat_obs):1,] #reverse latitude, to match with hindcast
time_obs =  ncvar_get(nc_obs, varid = "time")
unit_time = nc_obs$dim$time$units
nc_close(nc_obs)

time_obs_days = as.Date(time_obs, origin="1979-01-01")


init_dates <- substr(list.files(hind_data_path), 88, 97)
ref_sel_prec = which(time_obs_days>= as.Date(min(init_dates)) & time_obs_days<=(as.Date(max(init_dates))+45))
precip_obs = precip_obs[,,ref_sel_prec]
time_obs_days = time_obs_days[ref_sel_prec]

gc()

# Get percentiles ----

perc = 95


percentile = array(data=NA, dim=c(130,85,4)) # lon * lat * season


PERC_L <- foreach(season = 1:4) %dopar% { #1="MAM", 2="JJA",...
  
  months_in_season <- (((c(1,2,3)+(3*season-1))-1)%%12)+1
  
  date_ref <- which(month(time_obs_days) %in% months_in_season)
  
  return(apply(precip_obs[,,date_ref], MARGIN = c(1,2), FUN = quantile, probs = perc/100))
}#end foreach season

for (season in 1:4) {
  percentile[,,season] <- PERC_L[[season]]
}

gc()

# Store results in ncdf ---------------------------------------------------

ncfname <- paste0(exceed_folder, "/Observ_",perc,"thperc_season_dependant.nc")


londim <- ncdim_def(name = "lon", units = unit_lon, vals = lon)
latdim <- ncdim_def(name = "lat", units = unit_lat, vals = lat)
seasdim <- ncdim_def(name = "season", units = "season", vals = 1:4)

# define variables
fillvalue <- 1e32

perc_var <- ncvar_def(name = paste0("perc", perc), units = "mm/day",
                      dim = list(londim, latdim, seasdim), missval = fillvalue, prec = "float")

# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(perc_var),force_v4=TRUE)

# put variables
ncvar_put(ncout,perc_var,percentile)
nc_close(ncout)

# Get exceedences ----

exceedances = array(dim=dim(precip_obs))
for (season in 1:4) {
  months_in_season <- (((c(1,2,3)+(3*season-1))-1)%%12)+1
  
  date_ref <- which(month(time_obs_days) %in% months_in_season)
  
  for (d in date_ref) {
    exceedances[,,d] <- 1*(precip_obs[,,d] > percentile[,,season])
  }#end for d
  
}# end for season

gc()



# Store results in ncdf ---------------------------------------------------

ncfname <- paste0(exceed_folder, "Exceed_Observ_", perc,"thperc_season_dependant.nc")


londim <- ncdim_def(name = "lon", units = unit_lon, vals = lon)
latdim <- ncdim_def(name = "lat", units = unit_lat, vals = lat)
timedim <- ncdim_def(name = "time", units = "days since 1970-01-01", vals = as.numeric(time_obs_days))

# define variables
fillvalue <- 1e32

exceed_var <- ncvar_def(name = paste0("exceed", perc), units = "yon", dim = list(londim, latdim, timedim), prec = "integer")

# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(exceed_var),force_v4=TRUE)

# put variables
ncvar_put(ncout,exceed_var,exceedances)
nc_close(ncout)
