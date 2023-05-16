library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)

# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2.3)) #27 cores
perc=95

ref_WR = 6
# 0=no regime, 
# 1=Atlantic Trough (AT/AT), 
# 2=Atlantic Ridge (AR/AR), 
# 3=Greenland Blocking (GL/NAO-), 
# 4=European Blocking (EUBL/ZOWE),
# 5=Scandinavian Blocking (ScBL/BL), 
# 6=Zonal regime(ZO/ZO=NAO+),
# 7=Scandinavian Trough(ScTr/ZOEA)


# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""
WR_folder = ""


# Get the lon lat chunks --------------------------------------------------------------
size_chunks_lat = rep(3,28)
size_1_chunk_lat = 3
lat_start_ref = c(1,cumsum(size_chunks_lat[-length(size_chunks_lat)])+1)

dX <- 3*55.4327 # window size, in km, to get 3 GP in delta latitudes
dLat <- geosphere::destPoint(c(0,30),0,d=dX*1E3)[2]-30
# dLat %/% 0.5

latitudess=seq(30,72,0.5)[lat_start_ref+1]
dLon <- 0*latitudess
for (j in 1:length(latitudess)){
  dLon[j] <- geosphere::destPoint(c(0,latitudess[j]),90,d=dX*1E3)[1]
}
size_chunks_lon = dLon %/% 0.5

lon_start_ref_list = list()
for (lat_inter in 1:length(size_chunks_lon)) {
  sum_cumu = cumsum(rep(size_chunks_lon[lat_inter], (130 %/% size_chunks_lon[lat_inter])))
  lon_start_ref_list[[lat_inter]]=c(1,(sum_cumu[-length(sum_cumu)]+1))
}#end for lat_inter

source("Binary_loss_index_functions.R")

# Load weather regime, observations and details hindcast --------------------------------------------------------
A=read.csv("WR_folder/WR.csv")

all_dates_WR = character()
for (dd in 1:length(A$date.in.yyyymmdd_hh)) {
  all_dates_WR[dd] <- paste0(substr(A$date.in.yyyymmdd_hh[dd], start = 1, stop = 4),"-",
                             substr(A$date.in.yyyymmdd_hh[dd], start = 5, stop = 6),"-",
                             substr(A$date.in.yyyymmdd_hh[dd], start = 7, stop = 8), " ",
                             substr(A$date.in.yyyymmdd_hh[dd], start = 10, stop = 11), ":00")
}


data_WR = A$LC.attribution.based.on.WR.index[(all_dates_WR>="2001-01-04 00:00" & all_dates_WR<="2021-02-14 00:00" & hour(all_dates_WR) == 0)]
hours_since_19790101_00_WR = A$hour.since.19790101_00[(all_dates_WR>="2001-01-04 00:00" & all_dates_WR<="2021-02-14 00:00" & hour(all_dates_WR) == 0)]
dates_WRs = as.Date(as.POSIXct(hours_since_19790101_00_WR*3600,origin='1979-01-01 00:00'))
rm(A, all_dates_WR)

dates_1WR = dates_WRs[data_WR==ref_WR]

print(paste("number of days with WR",ref_WR,"=",length(dates_1WR)))

nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_", perc,"thperc_season_dependant.nc"))
Exceed_obs = ncvar_get(nc_obs, varid = paste0("exceed",perc))
time_obs =  ncvar_get(nc_obs, varid = "time")
nc_close(nc_obs)

time_obs_days = as.Date(time_obs, origin = "1970-01-01")

File_names_exceed <- list.files(hind_median_exceed_folder) # hindcast file names, one file per init
init_dates <- substr(File_names_exceed, 42, 51)
init_days = as.Date(init_dates)

gc()


# Compute BLF --------------------------------------------------------

BLF_spa_neighb_1ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spa_neighb_2ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spa_neighb_3ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spa_neighb_5ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
nb_geq1ev_WR <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
nb_geq2ev_WR <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
nb_geq1ev_WR_hind <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
nb_geq2ev_WR_hind <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))

for (lat in lat_start_ref) {
  ref_lat_chunck = which(lat_start_ref == lat) #have an index
  
  for (lon in lon_start_ref_list[[ref_lat_chunck]]) {
    size_1_chunk_lon = size_chunks_lon[ref_lat_chunck] #size of the chunk here
    
    sel_neighb_obs = Exceed_obs[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),]
    
    Exceed_list <- foreach (INI = 1:length(init_days)) %dopar% {
      Init = init_dates[INI]
      nc_hind = nc_open(paste0(hind_median_exceed_folder,"Exceedences_Hincast_",perc ,"thperc_init_",Init,".nc"))
      A = ncvar_get(nc_hind, varid = paste0("exceed_perc",perc), start = c(lon,lat,1,1), count = c(size_1_chunk_lon,size_1_chunk_lat,-1,-1))
      nc_close(nc_hind)
      gc()
      return(A)
    }#end for INI
    
    sel_hind_neighb <- array(dim = c(size_1_chunk_lon, size_1_chunk_lat, 11, 46, length(init_days)))
    for (INI in 1:length(init_days)) {
      sel_hind_neighb[,,,,INI] <- Exceed_list[[INI]]
    }#end for INI
    
    
    RES <- foreach(lead = 1:46) %dopar% {
      
      OUT_1 = Spati_acc_WR_Loss_func_lead(lead_time = lead, nb_ev = 1, dates_WR = dates_1WR, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                          square_hind_neighb = sel_hind_neighb, date_init = init_days)
      A_1=OUT_1$BLS
      
      nb_ev_1=OUT_1$nb_ev_obs
      nb_ev_1_hind=OUT_1$nb_ev_hind
      
      OUT_2 = Spati_acc_WR_Loss_func_lead(lead_time = lead, nb_ev = 2, dates_WR = dates_1WR, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                          square_hind_neighb = sel_hind_neighb, date_init = init_days)
      A_2=OUT_2$BLS
      
      nb_ev_2=OUT_2$nb_ev_obs
      nb_ev_2_hind=OUT_2$nb_ev_hind
      
      A_3=Spati_acc_WR_Loss_func_lead(lead_time = lead, nb_ev = 3, dates_WR = dates_1WR, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                      square_hind_neighb = sel_hind_neighb, date_init = init_days)$BLS
      A_5=Spati_acc_WR_Loss_func_lead(lead_time = lead, nb_ev = 5, dates_WR = dates_1WR, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                      square_hind_neighb = sel_hind_neighb, date_init = init_days)$BLS
      
      
      return(list(BLS_1=A_1, BLS_2=A_2, BLS_3=A_3, BLS_5=A_5,
                  nb_ev_1=nb_ev_1, nb_ev_2=nb_ev_2, nb_ev_1_hind=nb_ev_1_hind, nb_ev_2_hind=nb_ev_2_hind))
    }#end foreach lead
    
    for (lead in 1:46) {
      BLF_spa_neighb_1ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_1
      BLF_spa_neighb_2ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_2
      BLF_spa_neighb_3ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_3
      BLF_spa_neighb_5ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_5
      nb_geq1ev_WR[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$nb_ev_1
      nb_geq2ev_WR[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$nb_ev_2
      nb_geq1ev_WR_hind[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$nb_ev_1_hind
      nb_geq2ev_WR_hind[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$nb_ev_2_hind
    }#end for lead
    
  }#end for lon
}#end for lat




# Save in netcdf file ---------------------------------------------------

nc_hind = nc_open(paste0("/scratch3/pauline/S2S_precip_clustering/Data_05_Europe/Exceedances/perc",
                         perc,"/Exceedences_Hincast_",perc,"thperc_init_2001-01-04.nc"))
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

BLF_var_1 <- ncvar_def(name = "BLF_1ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_1event_or_more")
BLF_var_2 <- ncvar_def(name = "BLF_2ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_2events_or_more")
BLF_var_3 <- ncvar_def(name = "BLF_3ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_3events_or_more")
BLF_var_5 <- ncvar_def(name = "BLF_5ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_5events_or_more")
nb_geq1ev_WR_var <- ncvar_def(name = "nb_obs_geq1ev_WR", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                              longname = "nb_obs_times_nb_ev_geq1_WR")
nb_geq2ev_WR_var <- ncvar_def(name = "nb_obs_geq2ev_WR", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                              longname = "nb_obs_times_nb_ev_geq2_WR")
nb_geq1ev_WR_var_hind <- ncvar_def(name = "nb_hind_geq1ev_WR", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                                   longname = "nb_hind_times_nb_ev_geq1_WR")
nb_geq2ev_WR_var_hind <- ncvar_def(name = "nb_hind_geq2ev_WR", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                                   longname = "nb_hind_times_nb_ev_geq2_WR")

ncout <- nc_create(filename = paste0(BLI_folder,"/WR",ref_WR,"_spatiACC_seas_BLF_perlead_",perc,"th_perc.nc"),
                   vars = list(BLF_var_1,BLF_var_2,BLF_var_3,BLF_var_5,
                               nb_geq1ev_WR_var, nb_geq2ev_WR_var, nb_geq1ev_WR_var_hind, nb_geq2ev_WR_var_hind), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLF_var_1, vals = BLF_spa_neighb_1ev)
ncvar_put(nc = ncout, varid = BLF_var_2, vals = BLF_spa_neighb_2ev)
ncvar_put(nc = ncout, varid = BLF_var_3, vals = BLF_spa_neighb_3ev)
ncvar_put(nc = ncout, varid = BLF_var_5, vals = BLF_spa_neighb_5ev)
ncvar_put(nc = ncout, varid = nb_geq1ev_WR_var, vals = nb_geq1ev_WR)
ncvar_put(nc = ncout, varid = nb_geq2ev_WR_var, vals = nb_geq2ev_WR)
ncvar_put(nc = ncout, varid = nb_geq1ev_WR_var_hind, vals = nb_geq1ev_WR_hind)
ncvar_put(nc = ncout, varid = nb_geq2ev_WR_var_hind, vals = nb_geq2ev_WR_hind)

nc_close(ncout)
