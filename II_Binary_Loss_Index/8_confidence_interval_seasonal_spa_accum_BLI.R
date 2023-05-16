library(ncdf4);library(lubridate);library(geosphere)
library(foreach);library(iterators);library(parallel);library(doParallel)

# Parameters --------------------------------------------------------------

# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""


registerDoParallel(cores=floor(detectCores()/2.3)) 
perc=95
# ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO" #16h 27 cores
ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA" #16.5h 27 cores
alpha = 0.05
boot_size=500


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

BOOT_Spati_acc_seas_CI_BLF_lead <- function(nb_boot,lead_time, nb_ev, month_in_seas, square_obs_neighb, date_obs, date_init){
  
  Nb_events_obs_square = numeric()
  
  sub_init = which(month(date_init+lead_time-1) %in% month_in_seas)
  
  for (INIT in 1:length(sub_init)) {
    Nb_events_obs_square[INIT] <- sum(square_obs_neighb[,,which(date_obs==(date_init[sub_init[INIT]]+lead_time-1))])
  }#end for INIT
  
  Obs_sel = Nb_events_obs_square>=nb_ev
  
  Boostrapped_spa_BLF = numeric()
  
  for (boot_m in 1:nb_boot) {
    
    Nb_events_random_hind_square = numeric()
    sample_date_obs = sample(which(month(date_obs) %in% month_in_seas), size = length(sub_init), replace = F)
    
    for (INIT in 1:length(sub_init)) {
      Nb_events_random_hind_square[INIT] <- sum(square_obs_neighb[,,sample_date_obs[INIT]])
    }#end for INIT
    
    Hind_rand_sel = Nb_events_random_hind_square>=nb_ev
    
    Boostrapped_spa_BLF[boot_m] = mean((Obs_sel + Hind_rand_sel) == 1) / mean((Obs_sel + Hind_rand_sel) >= 1)
  }#end for boot_m
  
  return(Boostrapped_spa_BLF)
}# end spatial accumulation seas loss func




# Load observation and details hindcast --------------------------------------------------------


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
BLF_spat_climato95th_10ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spat_climato95th_15ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))
BLF_spat_climato95th_20ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))


for (lat in lat_start_ref) {
  ref_lat_chunck = which(lat_start_ref == lat) #have an index
  
  for (lon in lon_start_ref_list[[ref_lat_chunck]]) {
    size_1_chunk_lon = size_chunks_lon[ref_lat_chunck] #size of the chunk here
    
    sel_neighb_obs = Exceed_obs[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),]
    
    
    RES <- foreach(lead = 1:46) %dopar% {
      boot_BLF_1 <- BOOT_Spati_acc_seas_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 1, month_in_seas = ext_seas,
                                                    square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_1))){
        lowerbound_1ev <- NA
      } else {
        lowerbound_1ev <- quantile(boot_BLF_1, alpha)
      }
      
      boot_BLF_2 <- BOOT_Spati_acc_seas_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 2, month_in_seas = ext_seas,
                                                    square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_2))){
        lowerbound_2ev <- NA
      } else {
        lowerbound_2ev <- quantile(boot_BLF_2, alpha)
      }
      
      boot_BLF_3 <- BOOT_Spati_acc_seas_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 3, month_in_seas = ext_seas,
                                                    square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_3))){
        lowerbound_3ev <- NA
      } else {
        lowerbound_3ev <- quantile(boot_BLF_3, alpha)
      }
      
      boot_BLF_5 <- BOOT_Spati_acc_seas_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 5, month_in_seas = ext_seas,
                                                    square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_5))){
        lowerbound_5ev <- NA
      } else {
        lowerbound_5ev <- quantile(boot_BLF_5, alpha)
      }
      
      boot_BLF_10 <- BOOT_Spati_acc_seas_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 10, month_in_seas = ext_seas,
                                                     square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_10))){
        lowerbound_10ev <- NA
      } else {
        lowerbound_10ev <- quantile(boot_BLF_10, alpha)
      }
      
      boot_BLF_15 <- BOOT_Spati_acc_seas_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 15, month_in_seas = ext_seas,
                                                     square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_15))){
        lowerbound_15ev <- NA
      } else {
        lowerbound_15ev <- quantile(boot_BLF_15, alpha)
      }
      
      boot_BLF_20 <- BOOT_Spati_acc_seas_CI_BLF_lead(nb_boot = boot_size, lead_time = lead, nb_ev = 20, month_in_seas = ext_seas,
                                                     square_obs_neighb = sel_neighb_obs, date_init = init_days, date_obs = time_obs_days)
      if(is.na(sum(boot_BLF_20))){
        lowerbound_20ev <- NA
      } else {
        lowerbound_20ev <- quantile(boot_BLF_20, alpha)
      }
      
      
      return(list(lowerbound_1ev=lowerbound_1ev, lowerbound_2ev=lowerbound_2ev, lowerbound_3ev=lowerbound_3ev, lowerbound_5ev=lowerbound_5ev,
                  lowerbound_10ev=lowerbound_10ev, lowerbound_15ev=lowerbound_15ev, lowerbound_20ev=lowerbound_20ev))
    }#end foreach lead
    
    for (lead in 1:46) {
      BLF_spat_climato95th_1ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_1ev
      BLF_spat_climato95th_2ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_2ev
      BLF_spat_climato95th_3ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_3ev
      BLF_spat_climato95th_5ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_5ev
      BLF_spat_climato95th_10ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_10ev
      BLF_spat_climato95th_15ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_15ev
      BLF_spat_climato95th_20ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$lowerbound_20ev
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
BLF_var_10 <- ncvar_def(name = "BLF_climato_10ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                        longname = "Binary_Loss_Function_climato_10events_or_more")
BLF_var_15 <- ncvar_def(name = "BLF_climato_15ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                        longname = "Binary_Loss_Function_climato_15events_or_more")
BLF_var_20 <- ncvar_def(name = "BLF_climato_20ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                        longname = "Binary_Loss_Function_climato_20events_or_more")

ncout <- nc_create(filename = paste0(BLI_folder,"CI_", 100*(alpha),"thlowerbound_spatiACC_seasonal_BLF_perlead_andmember_",
                                     perc,"th_perc_",name_seas,".nc"),
                   vars = list(BLF_var_1,BLF_var_2,BLF_var_3,BLF_var_5,BLF_var_10,BLF_var_15,BLF_var_20), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLF_var_1, vals = BLF_spat_climato95th_1ev)
ncvar_put(nc = ncout, varid = BLF_var_2, vals = BLF_spat_climato95th_2ev)
ncvar_put(nc = ncout, varid = BLF_var_3, vals = BLF_spat_climato95th_3ev)
ncvar_put(nc = ncout, varid = BLF_var_5, vals = BLF_spat_climato95th_5ev)
ncvar_put(nc = ncout, varid = BLF_var_10, vals = BLF_spat_climato95th_10ev)
ncvar_put(nc = ncout, varid = BLF_var_15, vals = BLF_spat_climato95th_15ev)
ncvar_put(nc = ncout, varid = BLF_var_20, vals = BLF_spat_climato95th_20ev)

nc_close(ncout)
