library(ncdf4) 
# to adapt accordingly:
Brier_folder = ""
last_skillful_day_folder = ""

# Seasonal Local ----------------------------------------------------------------

perc = 95
season = "MJJASO"

nc_brier = nc_open(paste0(Brier_folder,"Brier_score_perlead_",perc,"th_perc_",season,".nc"))

BS = ncvar_get(nc_brier, varid = "Brier")
lon_dim = ncvar_get(nc_brier, varid = "lon")
lat_dim = ncvar_get(nc_brier, varid = "lat")
unit_lon = nc_brier$dim$lon$units
unit_lat = nc_brier$dim$lat$units
leadtime = ncvar_get(nc_brier, varid = "leadtime")
nc_close(nc_brier)

nc_climbrier = nc_open(paste0(Brier_folder,"Brier_score_climato_perlead_",perc,"th_perc_",season,"_CONSTCLIMATO.nc"))
BS_clim = ncvar_get(nc_climbrier, varid = "Brier")
nc_close(nc_climbrier)

BSS <- array(dim=dim(BS))
for (lead in 1:46) {
  BSS[,,lead] <- 1-(BS[,,lead]/BS_clim)
}#end for lead


Last_pos_BSS = array(dim=dim(BSS)[c(1,2)])

for (LON in 1:length(lon_dim)) {
  for (LAT in 1:length(lat_dim)) {
    Last_pos_BSS[LON,LAT] = min(which(BSS[LON, LAT,]<0), na.rm = T) - 1
  }#end for LAT
}#end for LON


# Save in ncdf file -------------------------------------------------------

fillvalue <- -9999
latdim <- ncdim_def("lat",unit_lat,lat_dim)
londim <- ncdim_def("lon",unit_lon,lon_dim)
leaddim <- ncdim_def("leadtime","days",1:46)

Brier_var <- ncvar_def(name = "BSS", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Brier_skill_score")

Lastday_var <- ncvar_def(name = "Last_week", units = "no", dim = list(londim,latdim), missval = fillvalue,
                         longname = "Last_skillful_week")

ncout <- nc_create(filename = paste0(last_skillful_day_folder,"BSS_lead_",perc,"th_perc_",season,".nc"),
                   vars = list(Brier_var, Lastday_var), force_v4=TRUE)

ncvar_put(nc = ncout, varid = Brier_var, vals = BSS)
ncvar_put(nc = ncout, varid = Lastday_var, vals = Last_pos_BSS)

nc_close(ncout)