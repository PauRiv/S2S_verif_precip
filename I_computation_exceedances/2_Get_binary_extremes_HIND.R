library(lubridate);library(ncdf4)
library(foreach);library(iterators);library(parallel);library(doParallel)

registerDoParallel(cores=floor(detectCores()/2))

# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""

# Define parameters --------------------------------------------------------

perc = 95

LON_chunk_size = 10 ; LAT_chunk_size = 5

File_names <- list.files(hind_data_path) # hindcast file names, one file per init
init_dates <- substr(File_names, 88, 97) # all the init dates, in the order of the files in the folder
SEASONS <- c("MAM", "JJA", "SON", "DJF")


# useful function ---------------------------------------------------------

extract_lead_seas_percentile_chunk <- function(perc, lon_start, lon_chunk_size, lat_start, lat_chunk_size, lead, season, file_names, init_dates){
  SEASONS <- c("MAM", "JJA", "SON", "DJF")
  
  ref_seas <- which(SEASONS==season)
  
  months_in_season <- (((c(1,2,3)+(3*ref_seas-1))-1)%%12)+1
  
  files_ref <- which(month(as.Date(init_dates) + lead - 1) %in% months_in_season)
  
  data_precip <- array(dim=c(lon_chunk_size, lat_chunk_size, 11, length(files_ref)))
  count <- 0
  for (file in file_names[files_ref]) {
    count <- count + 1
    nc_f <- nc_open(paste0(hind_data_path,file))
    data_precip[,,,count] <- ncvar_get(nc_f, varid = "tp", start = c(lon_start, lat_start, 1, lead),
                                       count = c(lon_chunk_size, lat_chunk_size, -1, 1))
    nc_close(nc_f)
    gc()
  }# end for file
  
  return(apply(data_precip, MARGIN = c(1,2), FUN = quantile, probs = perc/100))
  
}# end function extract_data


# Loop on Lon, Lat, season, lead to extract percentile ---------------------------------

percentile = array(data=NA, dim=c(130,85,46,4)) # lon * lat * leadtime * season

for (LON_start in seq(1,130,LON_chunk_size)) {
  
  for (LAT_start in seq(1,85,LAT_chunk_size)) {
    
    for (seas in c("MAM", "JJA", "SON", "DJF")) {
      
      PERC_CHUNK <- foreach(lead=1:46) %dopar% {
        return(extract_lead_seas_percentile_chunk(perc = perc, lon_start = LON_start, lon_chunk_size = LON_chunk_size,
                                                  lat_start =  LAT_start, lat_chunk_size = LAT_chunk_size, lead = lead, season = seas,
                                                  file_names = File_names, init_dates = init_dates))
      }#foreach lead : with 16 cores , 4min, with 23 cores 3.4 min
      gc()
      for (lead in 1:46) {
        percentile[(LON_start:(LON_start+LON_chunk_size-1)), (LAT_start:(LAT_start+LAT_chunk_size-1)), lead, which(SEASONS==seas)] <- PERC_CHUNK[[lead]]
      }#end for lead
      
    }#end for seas
    
  } # end for LON_start
  
} # end for LON_start




# Store results in ncdf ---------------------------------------------------

ncfname <- paste0(hind_exceed_folder,"/Hincast_",perc,"thperc_daylead_season_dependant.nc")

nc_f <- nc_open(hind_data_path,"ecmwf_reforecast_daily_precipitation_Europe_0.5deg_model_date_2021-01-04_hindcast_date_2001-01-04.nc")
lon <- ncvar_get(nc_f, varid = "lon")
unit_lon <- nc_f$dim$lon$units
lat <- ncvar_get(nc_f, varid = "lat")
unit_lat <- nc_f$dim$lat$units
nc_close(nc_f)

londim <- ncdim_def(name = "lon", units = unit_lon, vals = lon[1:130])
latdim <- ncdim_def(name = "lat", units = unit_lat, vals = lat)
leaddim <- ncdim_def(name = "lead", units = "days", vals = 1:46)
seasdim <- ncdim_def(name = "season", units = "season", vals = 1:4)

# define variables
fillvalue <- 1e32

perc_var <- ncvar_def(name = paste0("perc", perc), units = "mm/day", dim = list(londim, latdim, leaddim, seasdim), missval = fillvalue, prec = "float")

# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(perc_var),force_v4=TRUE)

# put variables
ncvar_put(ncout,perc_var,percentile)
nc_close(ncout)



# Get exceedences: Loop on files and lead ---------------------------------

foreach (file_nb = 1:length(File_names)) %dopar% {
  nc_f <- nc_open(paste0(hind_data_path,File_names[file_nb]))
  data_precip <- ncvar_get(nc_f, varid = "tp")
  lon <- ncvar_get(nc_f, varid = "lon")
  unit_lon <- nc_f$dim$lon$units
  lat <- ncvar_get(nc_f, varid = "lat")
  unit_lat <- nc_f$dim$lat$units
  timed <- ncvar_get(nc_f, varid = "time")
  unit_time <- nc_f$dim$time$units
  member <- ncvar_get(nc_f, varid = "member")
  unit_member <- nc_f$dim$member$units
  nc_close(nc_f)
  
  data_exceedances <- array(dim = c(130,dim(data_precip)[2:4]))
  
  for (lead in 1:46) {
    date=as.Date(init_dates[file_nb]) + lead - 1
    data_exceedances_inter <- array(dim = dim(data_exceedances)[1:3])
    seas_ref = floor(((month(date)-3)%%12)/3)+1 #find ref of the corresponding season (c("MAM", "JJA", "SON", "DJF") = 1:4)
    
    for (member in 1:11) {
      data_exceedances_inter[,,member] <- data_precip[1:130,,member,lead] > percentile[,, lead, seas_ref]
    }#end for member
    
    data_exceedances[,,,lead] <- data_exceedances_inter
  }#end for lead
  gc()
  
  ncfname <- paste0(exceed_folder,"Exceedences_Hincast_",perc,"thperc_init_",init_dates[file_nb],".nc")
  
  londim <- ncdim_def(name = "lon", units = unit_lon, vals = lon[1:130])
  latdim <- ncdim_def(name = "lat", units = unit_lat, vals = lat)
  timedim <- ncdim_def(name = "time", units = unit_time, vals = timed)
  ensdim <- ncdim_def(name = "ensemble", units = unit_member, vals = 1:11)
  
  # define variables
  fillvalue <- 1e32
  
  exceed_var <- ncvar_def(name = paste0("exceed_perc", perc), units = "yon", dim = list(londim, latdim, ensdim, timedim), prec = "integer")
  
  # create netCDF file and put arrays
  ncout <- nc_create(ncfname,list(exceed_var),force_v4=TRUE)
  
  # put variables
  ncvar_put(ncout,exceed_var,1*data_exceedances)
  nc_close(ncout)
  
}#end foreach file_nb