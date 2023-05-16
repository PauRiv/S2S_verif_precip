library(ncdf4);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)



# to adapt accordingly:
obs_data_path = "" 
hind_data_path = "" 
exceed_folder = ""
hind_exceed_folder = ""
hind_median_exceed_folder = ""
BLI_folder = ""

# Parameters --------------------------------------------------------------

registerDoParallel(cores=floor(detectCores()/2.8)) 
perc=95
# ext_seas = c(5,6,7,8,9,10) ; name_seas="MJJASO"
ext_seas = c(11,12,1,2,3,4) ; name_seas="NDJFMA"


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

# Load observation and details hindcast --------------------------------------------------------
nc_obs=nc_open(paste0(exceed_folder,"Exceed_Observ_",perc,"thperc_season_dependant.nc"))
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
BLF_spa_neighb_10ev <- array(dim=c(dim(Exceed_obs)[c(1,2)], 46))

for (lat in lat_start_ref) {
  ref_lat_chunck = which(lat_start_ref == lat) #have an index
  
  for (lon in lon_start_ref_list[[ref_lat_chunck]]) {
    size_1_chunk_lon = size_chunks_lon[ref_lat_chunck] #size of the chunk here
    
    sel_neighb_obs = Exceed_obs[lon:(lon+size_1_chunk_lon-1),lat:(lat+size_1_chunk_lat-1),]
    
    Exceed_list <- foreach (INI = 1:length(init_days)) %dopar% {
      Init = init_dates[INI]
      nc_hind = nc_open(paste0(hind_exceed_folder,"Exceedences_Hincast_",perc ,"thperc_init_",Init,".nc"))
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
      
      A_1=Spati_acc_seas_Loss_func_lead(lead_time = lead, nb_ev = 1, month_in_seas = ext_seas, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                        square_hind_neighb = sel_hind_neighb, date_init = init_days)
      A_2=Spati_acc_seas_Loss_func_lead(lead_time = lead, nb_ev = 2, month_in_seas = ext_seas, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                        square_hind_neighb = sel_hind_neighb, date_init = init_days)
      A_3=Spati_acc_seas_Loss_func_lead(lead_time = lead, nb_ev = 3, month_in_seas = ext_seas, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                        square_hind_neighb = sel_hind_neighb, date_init = init_days)
      A_5=Spati_acc_seas_Loss_func_lead(lead_time = lead, nb_ev = 5, month_in_seas = ext_seas, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                        square_hind_neighb = sel_hind_neighb, date_init = init_days)
      A_10=Spati_acc_seas_Loss_func_lead(lead_time = lead, nb_ev = 10, month_in_seas = ext_seas, square_obs_neighb = sel_neighb_obs, date_obs = time_obs_days,
                                         square_hind_neighb = sel_hind_neighb, date_init = init_days)
      
      return(list(BLS_1=A_1, BLS_2=A_2, BLS_3=A_3, BLS_5=A_5, BLS_10=A_10))
    }#end foreach lead
    
    for (lead in 1:46) {
      BLF_spa_neighb_1ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_1
      BLF_spa_neighb_2ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_2
      BLF_spa_neighb_3ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_3
      BLF_spa_neighb_5ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_5
      BLF_spa_neighb_10ev[lon+(floor(size_1_chunk_lon/2)),lat+(floor(size_1_chunk_lat/2)),lead] <- RES[[lead]]$BLS_10
    }#end for lead
    
  }#end for lon
  print(paste("end lat", ref_lat_chunck, "out of", length(lat_start_ref)))
}#end for lat




# Save in netcdf file ---------------------------------------------------

nc_hind = nc_open(paste0(hind_exceed_folder,"/Exceedences_Hincast_",perc,"thperc_init_2001-01-04.nc"))
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
BLF_var_10 <- ncvar_def(name = "BLF_10ev_spa", units = "no", dim = list(londim,latdim,leaddim), missval = fillvalue,
                        longname = "Binary_Loss_Function_10events_or_more")

ncout <- nc_create(filename = paste0(BLI_folder,"spatiACC_seas_BLF_perlead_",  perc,"th_perc_",name_seas,".nc"),
                   vars = list(BLF_var_1,BLF_var_2,BLF_var_3,BLF_var_5,BLF_var_10), force_v4=TRUE)

ncvar_put(nc = ncout, varid = BLF_var_1, vals = BLF_spa_neighb_1ev)
ncvar_put(nc = ncout, varid = BLF_var_2, vals = BLF_spa_neighb_2ev)
ncvar_put(nc = ncout, varid = BLF_var_3, vals = BLF_spa_neighb_3ev)
ncvar_put(nc = ncout, varid = BLF_var_5, vals = BLF_spa_neighb_5ev)
ncvar_put(nc = ncout, varid = BLF_var_10, vals = BLF_spa_neighb_10ev)

nc_close(ncout)
