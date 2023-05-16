######### Get the median member for each lead time and inititalisation date #########

library(lubridate);library(ncdf4)
library(foreach);library(iterators);library(parallel);library(doParallel)

registerDoParallel(cores=floor(detectCores()/2.7))

# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""

# Define parameters --------------------------------------------------------

perc = 95

File_names <- list.files(paste0(hind_exceed_folder))


nc_exc <- nc_open(paste0(hind_exceed_folder,File_names[1]))
lon <- ncvar_get(nc_exc, varid = "lon")
unit_lon <- nc_exc$dim$lon$units
lat <- ncvar_get(nc_exc, varid = "lat")
unit_lat <- nc_exc$dim$lat$units
leaddim <- ncvar_get(nc_exc, varid = "time")
unit_lead <- nc_exc$dim$time$units
nc_close(nc_exc)

init_dates = substr(File_names, 35,44)

# Loop on files and lead ---------------------------------

foreach (file_nb = 1:length(File_names)) %dopar% {
  nc_exc <- nc_open(paste0(hind_exceed_folder,File_names[file_nb]))
  EXCEED <- ncvar_get(nc_exc, varid = "exceed_perc95")
  nc_close(nc_exc)
  gc()
  
  Median_exceed = apply(EXCEED, MARGIN = c(1,2,4), FUN = median)
  
  ncfname <- paste0(hind_median_exceed_folder,"Memb_Median_Exceed_Hincast_",perc,"thperc_init_",init_dates[file_nb],".nc")
  
  londim <- ncdim_def(name = "lon", units = unit_lon, vals = lon)
  latdim <- ncdim_def(name = "lat", units = unit_lat, vals = lat)
  timedim <- ncdim_def(name = "time", units = unit_lead, vals = leaddim)
  
  # define variables
  fillvalue <- 1e32
  
  exceed_var <- ncvar_def(name = paste0("med_memb_exceed", perc), units = "yon", dim = list(londim, latdim, timedim), prec = "integer")
  
  # create netCDF file and put arrays
  ncout <- nc_create(ncfname,list(exceed_var),force_v4=TRUE)
  
  # put variables
  ncvar_put(ncout,exceed_var, Median_exceed)
  nc_close(ncout)
}#end foreach file_nb. 27 sec per file: 40 min everyting?


gc()
