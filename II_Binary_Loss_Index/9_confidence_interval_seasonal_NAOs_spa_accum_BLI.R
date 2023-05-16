library(ncdf4);library(lubridate);library(geosphere)
library(foreach);library(iterators);library(parallel);library(doParallel)



# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""
WR_folder = ""


# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2.3)) #27 cores climcal 3
perc=95

alpha = 0.05
boot_size=500
ref_WR = 6 # 0:7 NAO-=3, NAO+=6 
# 0=no regime, 
# 1=Atlantic Trough (AT/AT), 
# 2=Atlantic Ridge (AR/AR), 
# 3=Greenland Blocking (GL/NAO-), 
# 4=European Blocking (EUBL/ZOWE),
# 5=Scandinavian Blocking (ScBL/BL), 
# 6=Zonal regime(ZO/ZO=NAO+),
# 7=Scandinavian Trough(ScTr/ZOEA)


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



# Useful function --------------------------------------------------------

BOOT_Spati_acc_WR_CI_BLF_lead <- function(nb_boot,lead_time, nb_ev, dates_WR, square_obs_neighb, date_obs, date_init){
  
  Nb_events_obs_square = numeric()
  
  sub_init = which((date_init+lead_time-1) %in% dates_WR)
  
  for (INIT in 1:length(sub_init)) {
    Nb_events_obs_square[INIT] <- sum(square_obs_neighb[,,which(date_obs==(date_init[sub_init[INIT]]+lead_time-1))])
  }#end for INIT
  
  Obs_sel = Nb_events_obs_square>=nb_ev
  
  Boostrapped_spa_BLF = numeric()
  
  for (boot_m in 1:nb_boot) {
    
    Nb_events_random_hind_square = numeric()
    sample_date_obs = sample(which(date_obs %in% dates_WR), size = length(sub_init), replace = F)
    
    for (INIT in 1:length(sub_init)) {
      Nb_events_random_hind_square[INIT] <- sum(square_obs_neighb[,,sample_date_obs[INIT]])
    }#end for INIT
    
    Hind_rand_sel = Nb_events_random_hind_square>=nb_ev
    
    Boostrapped_spa_BLF[boot_m] = mean((Obs_sel + Hind_rand_sel) == 1) / mean((Obs_sel + Hind_rand_sel) >= 1)
  }#end for boot_m
  
  return(Boostrapped_spa_BLF)
}# end spatial accumulation WR loss func




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

nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_", perc,"thperc_season_dependant.nc"))
Exceed_obs = ncvar_get(nc_obs, varid = paste0("exceed",perc))
time_obs =  ncvar_get(nc_obs, varid = "time")
nc_close(nc_obs)

time_obs_days = as.Date(time_obs, origin = "1970-01-01")

File_names_exceed <- list.files(hind_exceed_folder) # hindcast file names, one file per init
init_dates <- substr(File_names_exceed, 35, 44)
init_days = as.Date(init_dates)

gc()


# Compute BLF --------------------------------------------------------

BLF_spat_climato95th_1ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spat_climato95th_2ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spat_climato95th_3ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spat_climato95th_5ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))


for (lat in lat_start_ref) {
  ref_lat_chunck = which(lat_start_ref == lat) #have an index
  
  for (lon in lon_start_ref_list[[ref_lat_chunck]]) {
    size_1_chunk_lon = size_chunks_lon[ref_lat_chunck] #size of the chunk here
    
    sel_neighb_obs = Exceed_obs[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),]
    
    
    RES <- foreach(lead = 1:46) %dopar% {
      boot_BLF_1 <- BOOT_Spati_acc_WR_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 1, dates_WR = dates_1WR,
                                                  square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_1))){
        lowerbound_1ev <- NA
      } else {
        lowerbound_1ev <- quantile(boot_BLF_1, alpha)
      }
      
      boot_BLF_2 <- BOOT_Spati_acc_WR_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 2, dates_WR = dates_1WR,
                                                  square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_2))){
        lowerbound_2ev <- NA
      } else {
        lowerbound_2ev <- quantile(boot_BLF_2, alpha)
      }
      
      boot_BLF_3 <- BOOT_Spati_acc_WR_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 3, dates_WR = dates_1WR,
                                                  square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_3))){
        lowerbound_3ev <- NA
      } else {
        lowerbound_3ev <- quantile(boot_BLF_3, alpha)
      }
      
      boot_BLF_5 <- BOOT_Spati_acc_WR_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 5, dates_WR = dates_1WR,
                                                  square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_5))){
        lowerbound_5ev <- NA
      } else {
        lowerbound_5ev <- quantile(boot_BLF_5, alpha)
      }
      
      return(list(lowerbound_1ev=lowerbound_1ev, lowerbound_2ev=lowerbound_2ev, lowerbound_3ev=lowerbound_3ev, lowerbound_5ev=lowerbound_5ev))
    }#end foreach lead
    
    for (lead in 1:46) {
      BLF_spat_climato95th_1ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_1ev
      BLF_spat_climato95th_2ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_2ev
      BLF_spat_climato95th_3ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_3ev
      BLF_spat_climato95th_5ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_5ev
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
leaddim <- ncdim_def("leadtime","days",1:46)

BLF_var_1 <- ncvar_def(name = "BLF_climato_1ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_climato_1event_or_more")
BLF_var_2 <- ncvar_def(name = "BLF_climato_2ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_climato_2events_or_more")
BLF_var_3 <- ncvar_def(name = "BLF_climato_3ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_climato_3events_or_more")
BLF_var_5 <- ncvar_def(name = "BLF_climato_5ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                       longname = "Binary_Loss_Function_climato_5events_or_more")

ncout <- nc_create(filename = paste0(BLI_folder,"WR", ef_WR,"_CI",100*(alpha),"thlowerbound_spatiACC_BLF_perlead_andmember_",perc,"th_perc.nc"),
                   vars = list(BLF_var_1,BLF_var_2,BLF_var_3,BLF_var_5), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLF_var_1, vals = BLF_spat_climato95th_1ev)
ncvar_put(nc = ncout, varid = BLF_var_2, vals = BLF_spat_climato95th_2ev)
ncvar_put(nc = ncout, varid = BLF_var_3, vals = BLF_spat_climato95th_3ev)
ncvar_put(nc = ncout, varid = BLF_var_5, vals = BLF_spat_climato95th_5ev)

nc_close(ncout)
