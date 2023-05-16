library(ncdf4);library(geosphere)

# to adapt accordingly:
hind_exceed_folder = ""
BLI_folder = ""
last_skillful_day_folder = ""

# For simple local analysis -----------------------------------------------

perc=95
name_seas="MJJASO"
# name_seas="NDJFMA"

nc_CI5 <- nc_open(paste0(BLI_folder,"BLF_seas_",name_seas,"_confidenceinterval_climato_with_R_", perc,"th_perc.nc"))
lower_bound = ncvar_get(nc_CI5, varid = "BLS_clim_lowerbound")
nc_close(nc_CI5)
rm(nc_CI5)

nc_BLF <- nc_open(paste0(BLI_folder,"seasonal_BLF_perlead_medmemb_",perc,"th_perc_",name_seas,".nc"))
BLF_lead_memb = ncvar_get(nc_BLF, varid = "BLF")
nc_close(nc_BLF)

Get_last_day_skill <- function(vec_median_BLF, vec_lower_bound){
  # length of the vectors = leadtime
  return(min(which(vec_median_BLF-vec_lower_bound>0)) - 1)
}


Last_leadday_skill <- matrix(nrow = 130, ncol = 85)
for (LON in 1:130) {
  for (LAT in 1:85) {
    Last_leadday_skill[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb[LON,LAT,], vec_lower_bound = lower_bound[LON,LAT,])
  }#end for LAT
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

Last_leadday_skill_var <- ncvar_def(name = "Last_lead_outCI", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                    longname = "Last_lead_day_median_memb_lower_5thperc_CI")
ncout <- nc_create(filename = paste0(last_skillful_day_folder,"BLF_seas_Last_lead_day_median_lower_5thperc_CI_",perc,"th_perc_",name_seas,".nc"),
                   vars = list(Last_leadday_skill_var), force_v4=TRUE)

ncvar_put(nc = ncout, varid = Last_leadday_skill_var, vals = Last_leadday_skill)

nc_close(ncout)




# For temporal accumulated precip analysis -----------------------------------------------

perc=95
ndays_acc = 7
name_seas="MJJASO"
# name_seas="NDJFMA"

nc_CI5 <- nc_open(paste0(BLI_folder,"CI_5thlowerbound_tempACC_",ndays_acc,"days_seasonal_BLF_perlead_andmember_",perc,"th_perc_",name_seas,".nc"))
lower_bound_1ev = ncvar_get(nc_CI5, varid = "BLS_clim_lowerbound_1ev") #no lead time dimension, it is useless, accumlation over week lead1 covers all the climatology
lower_bound_2ev = ncvar_get(nc_CI5, varid = "BLS_clim_lowerbound_2ev")
lower_bound_3ev = ncvar_get(nc_CI5, varid = "BLS_clim_lowerbound_3ev")
nc_close(nc_CI5)


nc_BLF <- nc_open(paste0(BLI_folder,"ACC_",ndays_acc,"days_seasonal_BLF_perlead_",perc,"th_perc_",name_seas,".nc"))
BLF_lead_memb_1ev = ncvar_get(nc_BLF, varid = "BLF_1ev")
BLF_lead_memb_2ev = ncvar_get(nc_BLF, varid = "BLF_2ev")
BLF_lead_memb_3ev = ncvar_get(nc_BLF, varid = "BLF_3ev")
nc_close(nc_BLF)


Get_last_week_skill <- function(vec_median_BLF, vec_lower_bound){
  # length of the vectors = leadtime
  WHICH_first_in = which(vec_median_BLF-vec_lower_bound>0)
  if(length(WHICH_first_in)>0){
    return(min(WHICH_first_in) - 1)
  } else {
    return(0)
  }#end ifelse no skill
}

Last_leadweek_skill_1ev <- matrix(nrow = 130, ncol = 85)
Last_leadweek_skill_2ev <- matrix(nrow = 130, ncol = 85)
Last_leadweek_skill_3ev <- matrix(nrow = 130, ncol = 85)
for (LON in 1:130) {
  for (LAT in 1:85) {
    Last_leadweek_skill_1ev[LON,LAT] <- Get_last_week_skill(vec_median_BLF = BLF_lead_memb_1ev[LON,LAT,], vec_lower_bound = lower_bound_1ev[LON,LAT])
    Last_leadweek_skill_2ev[LON,LAT] <- Get_last_week_skill(vec_median_BLF = BLF_lead_memb_2ev[LON,LAT,], vec_lower_bound = lower_bound_2ev[LON,LAT])
    Last_leadweek_skill_3ev[LON,LAT] <- Get_last_week_skill(vec_median_BLF = BLF_lead_memb_3ev[LON,LAT,], vec_lower_bound = lower_bound_3ev[LON,LAT])
  }#end for LAT
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

Last_leadweek_skill_1ev_var <- ncvar_def(name = "Last_week_lowerCI_1ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                         longname = "Last_lead_week_lower_5thperc_CI_1ev")
Last_leadweek_skill_2ev_var <- ncvar_def(name = "Last_week_lowerCI_2ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                         longname = "Last_lead_week_lower_5thperc_CI_2ev")
Last_leadweek_skill_3ev_var <- ncvar_def(name = "Last_week_lowerCI_3ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                         longname = "Last_lead_week_lower_5thperc_CI_3ev")

ncout <- nc_create(filename = paste0(last_skillful_day_folder,"BLF_tempACC_",ndays_acc,"days_seasonal_Last_lead_day_median_lower_5thperc_CI_", perc,"th_perc_",name_seas,".nc"),
                   vars = list(Last_leadweek_skill_1ev_var,Last_leadweek_skill_2ev_var,Last_leadweek_skill_3ev_var), force_v4=TRUE)

ncvar_put(nc = ncout, varid = Last_leadweek_skill_1ev_var, vals = Last_leadweek_skill_1ev)
ncvar_put(nc = ncout, varid = Last_leadweek_skill_2ev_var, vals = Last_leadweek_skill_2ev)
ncvar_put(nc = ncout, varid = Last_leadweek_skill_3ev_var, vals = Last_leadweek_skill_3ev)

nc_close(ncout)




# For spatially accumulated precip analysis -----------------------------------------------

perc=95

# name_seas="MJJASO"
name_seas="NDJFMA"

nc_CI5 <- nc_open(paste0(BLI_folder,"CI_5thlowerbound_spatiACC_seasonal_BLF_perlead_andmember_",perc,"th_perc_",name_seas,".nc"))
lower_bound_1ev = ncvar_get(nc_CI5, varid = "BLF_climato_1ev_spa")
lower_bound_2ev = ncvar_get(nc_CI5, varid = "BLF_climato_2ev_spa")
lower_bound_3ev = ncvar_get(nc_CI5, varid = "BLF_climato_3ev_spa")
lower_bound_5ev = ncvar_get(nc_CI5, varid = "BLF_climato_5ev_spa")
lower_bound_10ev = ncvar_get(nc_CI5, varid = "BLF_climato_10ev_spa")
nc_close(nc_CI5)


nc_BLF <- nc_open(paste0(BLI_folder,"spatiACC_seas_BLF_perlead_",perc,"th_perc_",name_seas,".nc"))
BLF_lead_memb_1ev = ncvar_get(nc_BLF, varid = "BLF_1ev_spa")
BLF_lead_memb_2ev = ncvar_get(nc_BLF, varid = "BLF_2ev_spa")
BLF_lead_memb_3ev = ncvar_get(nc_BLF, varid = "BLF_3ev_spa")
BLF_lead_memb_5ev = ncvar_get(nc_BLF, varid = "BLF_5ev_spa")
BLF_lead_memb_10ev = ncvar_get(nc_BLF, varid = "BLF_10ev_spa")
nc_close(nc_BLF)


# plot(median_BLF_memb[1,1,])
# lines(lower_bound[1,1,])


Get_last_day_skill <- function(vec_median_BLF, vec_lower_bound){
  # length of the vectors = leadtime
  if(is.na(vec_median_BLF[1])){
    return(NA)
  } else {
    WHICH_first_in = which(vec_median_BLF-vec_lower_bound>0)
    if(length(WHICH_first_in)>0){
      return(min(WHICH_first_in) - 1)
    } else {
      return(0)
    }#end ifelse no skill
  }#end if no data zone
} #end Get_last_day_skill function

Last_leadday_skill_1ev <- matrix(nrow = 130, ncol = 85)
Last_leadday_skill_2ev = Last_leadday_skill_3ev = Last_leadday_skill_1ev
Last_leadday_skill_5ev = Last_leadday_skill_10ev= Last_leadday_skill_1ev

for (LON in 1:130) {
  for (LAT in 1:85) {
    Last_leadday_skill_1ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_1ev[LON,LAT,], vec_lower_bound = lower_bound_1ev[LON,LAT,])
    Last_leadday_skill_2ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_2ev[LON,LAT,], vec_lower_bound = lower_bound_2ev[LON,LAT,])
    Last_leadday_skill_3ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_3ev[LON,LAT,], vec_lower_bound = lower_bound_3ev[LON,LAT,])
    Last_leadday_skill_5ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_5ev[LON,LAT,], vec_lower_bound = lower_bound_5ev[LON,LAT,])
    Last_leadday_skill_10ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_10ev[LON,LAT,], vec_lower_bound = lower_bound_10ev[LON,LAT,])
  }#end for LAT
}#end for lon



# fill the blanks with value of the center
size_chunks_lat = rep(3,28)
size_1_chunk_lat = 3
lat_start_ref = c(1,cumsum(size_chunks_lat[-length(size_chunks_lat)])+1)

dX <- 3*55.4327 # window size, in km, to get 3 GP in delta latitudes
dLat <- destPoint(c(0,30),0,d=dX*1E3)[2]-30
# dLat %/% 0.5

latitudess=seq(30,72,0.5)[lat_start_ref+1]
dLon <- 0*latitudess
for (j in 1:length(latitudess)){
  dLon[j] <- destPoint(c(0,latitudess[j]),90,d=dX*1E3)[1]
}
size_chunks_lon = dLon %/% 0.5

lon_start_ref_list = list()
for (lat_inter in 1:length(size_chunks_lon)) {
  sum_cumu = cumsum(rep(size_chunks_lon[lat_inter], (130 %/% size_chunks_lon[lat_inter])))
  lon_start_ref_list[[lat_inter]]=c(1,(sum_cumu[-length(sum_cumu)]+1))
}#end for lat_inter


for (lat in lat_start_ref) {
  ref_lat_chunck = which(lat_start_ref == lat) #have an index
  
  for (lon in lon_start_ref_list[[ref_lat_chunck]]) {
    size_1_chunk_lon = size_chunks_lon[ref_lat_chunck] #size of the chunk here
    
    Last_leadday_skill_1ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_1ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    Last_leadday_skill_2ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_2ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    Last_leadday_skill_3ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_3ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    Last_leadday_skill_5ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_5ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    Last_leadday_skill_10ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_10ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                                lat+(floor(size_1_chunk_lat/2))],
                                                                                                 nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
  }#end for lon
}#end for lat










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

Last_leadday_skill_1ev_var <- ncvar_def(name = "Last_lead_lowerCI_1ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_1ev")
Last_leadday_skill_2ev_var <- ncvar_def(name = "Last_lead_lowerCI_2ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_2ev")
Last_leadday_skill_3ev_var <- ncvar_def(name = "Last_lead_lowerCI_3ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_3ev")
Last_leadday_skill_5ev_var <- ncvar_def(name = "Last_lead_lowerCI_5ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_5ev")
Last_leadday_skill_10ev_var <- ncvar_def(name = "Last_lead_lowerCI_10ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                         longname = "Last_lead_day_median_lower_5thperc_CI_10ev")
ncout <- nc_create(filename = paste0(last_skillful_day_folder,"BLF_spatiACC_seasonal_Last_lead_day_median_lower_5thperc_CI_",perc,"th_perc_",name_seas,".nc"),
                   vars = list(Last_leadday_skill_1ev_var,Last_leadday_skill_2ev_var,Last_leadday_skill_3ev_var,Last_leadday_skill_5ev_var,
                               Last_leadday_skill_10ev_var), force_v4=TRUE)

ncvar_put(nc = ncout, varid = Last_leadday_skill_1ev_var, vals = Last_leadday_skill_1ev)
ncvar_put(nc = ncout, varid = Last_leadday_skill_2ev_var, vals = Last_leadday_skill_2ev)
ncvar_put(nc = ncout, varid = Last_leadday_skill_3ev_var, vals = Last_leadday_skill_3ev)
ncvar_put(nc = ncout, varid = Last_leadday_skill_5ev_var, vals = Last_leadday_skill_5ev)
ncvar_put(nc = ncout, varid = Last_leadday_skill_10ev_var, vals = Last_leadday_skill_10ev)

nc_close(ncout)




# For weather regime spatially accumulated precip analysis -----------------------------------------------

perc=95

ref_WR = 3  # NAO-=3, NAO+=6

nc_CI5 <- nc_open(paste0(BLI_folder,"WR",ref_WR,"_CI5thlowerbound_spatiACC_BLF_perlead_andmember_",perc,"th_perc.nc"))
lower_bound_1ev = ncvar_get(nc_CI5, varid = "BLF_climato_1ev_spa")
lower_bound_2ev = ncvar_get(nc_CI5, varid = "BLF_climato_2ev_spa")
lower_bound_3ev = ncvar_get(nc_CI5, varid = "BLF_climato_3ev_spa")
lower_bound_5ev = ncvar_get(nc_CI5, varid = "BLF_climato_5ev_spa")
nc_close(nc_CI5)


nc_BLF <- nc_open(paste0(BLI_folder,"WR",ref_WR,"_spatiACC_seas_BLF_perlead_",perc,"th_perc.nc"))
BLF_lead_memb_1ev = ncvar_get(nc_BLF, varid = "BLF_1ev_spa")
BLF_lead_memb_2ev = ncvar_get(nc_BLF, varid = "BLF_2ev_spa")
BLF_lead_memb_3ev = ncvar_get(nc_BLF, varid = "BLF_3ev_spa")
BLF_lead_memb_5ev = ncvar_get(nc_BLF, varid = "BLF_5ev_spa")
nb_1ev_obs = ncvar_get(nc_BLF, varid = "nb_obs_geq1ev_WR")
nb_2ev_obs = ncvar_get(nc_BLF, varid = "nb_obs_geq2ev_WR")
nb_1ev_hind = ncvar_get(nc_BLF, varid = "nb_hind_geq1ev_WR")
nb_2ev_hind = ncvar_get(nc_BLF, varid = "nb_hind_geq2ev_WR")
nc_close(nc_BLF)

Get_last_day_skill <- function(vec_median_BLF, vec_lower_bound){
  # length of the vectors = leadtime
  if(is.na(vec_median_BLF[1])){
    return(NA)
  } else {
    WHICH_first_in = which(vec_median_BLF-vec_lower_bound>0)
    if(length(WHICH_first_in)>0){
      return(min(WHICH_first_in) - 1)
    } else {
      return(0)
    }#end ifelse no skill
  }#end if no data zone
} #end Get_last_day_skill function

Last_leadday_skill_1ev <- matrix(nrow = 130, ncol = 85)
Last_leadday_skill_2ev = Last_leadday_skill_3ev = Last_leadday_skill_1ev
Last_leadday_skill_5ev = Last_leadday_skill_1ev

for (LON in 1:130) {
  for (LAT in 1:85) {
    Last_leadday_skill_1ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_1ev[LON,LAT,], vec_lower_bound = lower_bound_1ev[LON,LAT,])
    Last_leadday_skill_2ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_2ev[LON,LAT,], vec_lower_bound = lower_bound_2ev[LON,LAT,])
    Last_leadday_skill_3ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_3ev[LON,LAT,], vec_lower_bound = lower_bound_3ev[LON,LAT,])
    Last_leadday_skill_5ev[LON,LAT] <- Get_last_day_skill(vec_median_BLF = BLF_lead_memb_5ev[LON,LAT,], vec_lower_bound = lower_bound_5ev[LON,LAT,])
  }#end for LAT
}#end for lon



# fill the blanks with value of the center
size_chunks_lat = rep(3,28)
size_1_chunk_lat = 3
lat_start_ref = c(1,cumsum(size_chunks_lat[-length(size_chunks_lat)])+1)

dX <- 3*55.4327 # window size, in km, to get 3 GP in delta latitudes
dLat <- destPoint(c(0,30),0,d=dX*1E3)[2]-30
# dLat %/% 0.5

latitudess=seq(30,72,0.5)[lat_start_ref+1]
dLon <- 0*latitudess
for (j in 1:length(latitudess)){
  dLon[j] <- destPoint(c(0,latitudess[j]),90,d=dX*1E3)[1]
}
size_chunks_lon = dLon %/% 0.5

lon_start_ref_list = list()
for (lat_inter in 1:length(size_chunks_lon)) {
  sum_cumu = cumsum(rep(size_chunks_lon[lat_inter], (130 %/% size_chunks_lon[lat_inter])))
  lon_start_ref_list[[lat_inter]]=c(1,(sum_cumu[-length(sum_cumu)]+1))
}#end for lat_inter


nb_1ev_hind_mat <- array(dim=c(130,85,46))
nb_1ev_obs_mat = nb_2ev_obs_mat = nb_2ev_hind_mat = nb_1ev_hind_mat

for (lat in lat_start_ref) {
  ref_lat_chunck = which(lat_start_ref == lat) #have an index
  
  for (lon in lon_start_ref_list[[ref_lat_chunck]]) {
    size_1_chunk_lon = size_chunks_lon[ref_lat_chunck] #size of the chunk here
    
    Last_leadday_skill_1ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_1ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    Last_leadday_skill_2ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_2ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    Last_leadday_skill_3ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_3ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    Last_leadday_skill_5ev[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1)] <- matrix(data = Last_leadday_skill_5ev[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                              lat+(floor(size_1_chunk_lat/2))],
                                                                                                nrow = size_1_chunk_lon,ncol = size_1_chunk_lat)
    
    for (lead in 1:46) {
      
      nb_1ev_obs_mat[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),lead] <- array(data = nb_1ev_obs[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                lat+(floor(size_1_chunk_lat/2)),lead],
                                                                                              dim = c(size_1_chunk_lon,size_1_chunk_lat))
      
      nb_2ev_obs_mat[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),lead] <- array(data = nb_2ev_obs[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                lat+(floor(size_1_chunk_lat/2)),lead],
                                                                                              dim = c(size_1_chunk_lon,size_1_chunk_lat))
      
      nb_1ev_hind_mat[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),lead] <- array(data = nb_1ev_hind[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                  lat+(floor(size_1_chunk_lat/2)),lead],
                                                                                               dim = c(size_1_chunk_lon,size_1_chunk_lat))
      
      nb_2ev_hind_mat[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),lead] <- array(data = nb_2ev_hind[lon+(floor(size_1_chunk_lon/2)),
                                                                                                                  lat+(floor(size_1_chunk_lat/2)),lead],
                                                                                               dim = c(size_1_chunk_lon,size_1_chunk_lat))
    }#end for lead
    
  }#end for lon
}#end for lat










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
leaddim <-  ncdim_def("leaddays","days",1:46)

Last_leadday_skill_1ev_var <- ncvar_def(name = "Last_lead_lowerCI_1ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_1ev")
Last_leadday_skill_2ev_var <- ncvar_def(name = "Last_lead_lowerCI_2ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_2ev")
Last_leadday_skill_3ev_var <- ncvar_def(name = "Last_lead_lowerCI_3ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_3ev")
Last_leadday_skill_5ev_var <- ncvar_def(name = "Last_lead_lowerCI_5ev", units = "th_day", dim = list(londim,latdim), missval = fillvalue,
                                        longname = "Last_lead_day_median_lower_5thperc_CI_5ev")
nb_days_1ev_obs_var <- ncvar_def(name = "nb_days_1ev_obs", units = "nb_days", dim = list(londim,latdim,leaddim), missval = fillvalue,
                                 longname = "nb_days_with1evormore_obs")
nb_days_2ev_obs_var <- ncvar_def(name = "nb_days_2ev_obs", units = "nb_days", dim = list(londim,latdim,leaddim), missval = fillvalue,
                                 longname = "nb_days_with2evormore_obs")
nb_days_1ev_hind_var <- ncvar_def(name = "nb_days_1ev_hind", units = "nb_days", dim = list(londim,latdim,leaddim), missval = fillvalue,
                                  longname = "nb_days_with1evormore_hind")
nb_days_2ev_hind_var <- ncvar_def(name = "nb_days_2ev_hind", units = "nb_days", dim = list(londim,latdim,leaddim), missval = fillvalue,
                                  longname = "nb_days_with2evormore_hind")

ncout <- nc_create(filename = paste0(last_skillful_day_folder,"/BLF_WR",ref_WR,"_spatiACC_seas_Last_lead_day_median_lower_5thperc_CI_",perc,"th_perc.nc"),
                   vars = list(Last_leadday_skill_1ev_var,Last_leadday_skill_2ev_var,Last_leadday_skill_3ev_var,Last_leadday_skill_5ev_var,
                               nb_days_1ev_obs_var,nb_days_2ev_obs_var,nb_days_1ev_hind_var,nb_days_2ev_hind_var), force_v4=TRUE)

ncvar_put(nc = ncout, varid = Last_leadday_skill_1ev_var, vals = Last_leadday_skill_1ev)
ncvar_put(nc = ncout, varid = Last_leadday_skill_2ev_var, vals = Last_leadday_skill_2ev)
ncvar_put(nc = ncout, varid = Last_leadday_skill_3ev_var, vals = Last_leadday_skill_3ev)
ncvar_put(nc = ncout, varid = Last_leadday_skill_5ev_var, vals = Last_leadday_skill_5ev)
ncvar_put(nc = ncout, varid = nb_days_1ev_obs_var, vals = nb_1ev_obs_mat)
ncvar_put(nc = ncout, varid = nb_days_2ev_obs_var, vals = nb_2ev_obs_mat)
ncvar_put(nc = ncout, varid = nb_days_1ev_hind_var, vals = nb_1ev_hind_mat)
ncvar_put(nc = ncout, varid = nb_days_2ev_hind_var, vals = nb_2ev_hind_mat)

nc_close(ncout)

