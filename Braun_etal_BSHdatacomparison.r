
# prep --------------------------------------------------------------

library(reshape2)
library(raster)
library(dismo)
library(ggplot2)
library(dplyr)

# get env --------------------------------------------------------------

## see ./pipelines/cmems/get_cmems_daily.r

## see ./pipelines/cmems/get_cmems_monthly.r


# prep iccat --------------------------------------------------------------
## make sure our downloaded, qc'd iccat data is up to date (all species)
#system('dvc repro iccat_combine')
sp <- c('ALB', 'BET', 'BFT', 'BSH', 'BUM', 'POR', 'SAI', 'SKJ', 'SMA', 'SWO', 'WHM', 'YFT')
for (i in sp){
  system(paste0('docker-compose run iccat ./R/iccat_qc.r /data/iccat/download/_tag', i, '.xlsx /data/bathy/global_bathy_0.01.nc /data/iccat/qc/_tag', i, '.csv'))
}

## combine all qc'd into master iccat
system('docker-compose run py_tools ./combine.py /data/iccat/qc/ /data/iccat/combined_tags.csv')

## subset that master data to our study area
xl <- c(-100, -5); yl <- c(10, 55)
df <- data.table::fread('./data/iccat/combined_tags.csv')
df <- df %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2] & date >= as.POSIXct('1993-01-01') & date <= as.POSIXct('2019-12-31'))
## mask to atlantic only (i.e. remove any weird stuff in great lakes or in ETP)
atl <- rgdal::readOGR(dsn='./data/shapefiles/',layer='atlantic') ## prevents great lakes and others that bathymetry mask doesn't catch
tmp <- df
coordinates(tmp) <- ~lon + lat
proj4string(tmp) <- proj4string(atl)
df <- df[which(!is.na(over(tmp, as(atl, "SpatialPolygons")))),]
data.table::fwrite(df, './data/iccat/combined_tags_sub.csv')

## generate pseudoabs using background sampling
#system('dvc repro iccat_pseudoabs')
system('docker-compose run etag ./R/generate_pseudoabs.r /data/iccat/combined_tags_sub.csv --index_var SpeciesCode --bathy_file /data/bathy/global_bathy_0.01.nc --abs_ratio 5 --force_180 /data/iccat/with_pseudoabs.csv')

## split daily
system('docker compose run enhance ./split_by_date.py /data/iccat/with_pseudoabs.csv /data/enhance/iccat/split_daily/')

## split monthly
system('chmod +x ./pipelines/split_by_month.r')
system('./split_by_month.r ../data/iccat/with_pseudoabs.csv ../data/enhance/iccat/split_monthly/')


# enhance iccat --------------------------------------------------------------
## enhance iccat w/ cmems daily
fList <- list.files('./data/enhance/iccat/split_daily/', recursive = TRUE, full.names = TRUE)
cmems_files <- list.files('~/work/EnvData/glorys_daily/', recursive = TRUE, full.names = TRUE)
for (i in 1:length(fList)){
  df <- data.table::fread(fList[i])
  raster_name <- cmems_files[grep(substr(fList[i], 35, nchar(fList[i]) - 4), cmems_files)]
  if (length(raster_name) == 0) next
  raster_name <- raster_name[grep('.grd', raster_name)]
  extr <- data.frame(raster::extract(raster::brick(raster_name), cbind(df$lon, df$lat)))
  names(extr) <- c('sst','sss','ssh','mld', 'log_eke','sst_sd','ssh_sd','sss_sd','bathy','rugosity')
  df <- cbind(df, extr)
  out_name <- paste0('./data/enhance/iccat/enhanced_daily/', substr(fList[i], 35, nchar(fList[i])))
  data.table::fwrite(df, out_name)
  rm(df); rm(raster_name); rm(out_name)
  print(substr(fList[i], 35, nchar(fList[i])))
}

## recombine all the daily enhanced data to master csv
#system('dvc repro enhance_iccat_combine')
#system('docker compose run enhance ./combine.py /data/enhance/iccat/enhanced/ /data/enhance/iccat/iccat-enhanced.csv')
input_dir <- '~/work/RCode/NASA-FaCeT/data/enhance/iccat/enhanced_daily/'
output_csv <- '~/work/RCode/NASA-FaCeT/data/enhance/iccat/iccat-enhanced.csv'
fList <- list.files(input_dir, full.names = TRUE)
fList <- fList[grep('.csv', fList)]
for (i in 1:length(fList)){
  df <- data.table::fread(fList[i])
  if (file.exists(output_csv)){
    data.table::fwrite(df, output_csv, append=TRUE)
  } else{
    data.table::fwrite(df, output_csv, append=FALSE)
  }
}

## sample the data to a 1:1 ratio before model fit
df <- data.table::fread('./data/enhance/iccat/iccat-enhanced.csv')
df <- data.frame(na.omit(df, cols = c(which(names(df) == 'sst'):ncol(df))))
#nms <- names(df)
#nms <- nms[which(nms != 'bathy')]
#df <- df %>% dplyr::select(all_of(nms))
#names(df)[which(names(df) == 'bathy.1')] = 'bathy'
df.split <- split(df, df$SpeciesCode)
df <- lapply(df.split, FUN = function(x){
  df.1 <- x %>% filter(pres == 1)
  df.0 <- x %>% filter(pres == 0)
  set.seed(311)
  df.0 <- df.0[sample(1:nrow(df.0), size = nrow(df.1)),]
  df <- rbind(df.1, df.0)
}) %>% do.call(rbind, .)
df <- df[order(df$SpeciesCode, df$date),]
data.table::fwrite(df, file = './data/enhance/iccat/iccat-enhanced_1to1.csv')


# fit ICCAT BRTs --------------------------------------------------------------
#system('dvc repro iccat_fit_brt')
#for (i in sp){
#  system(paste0('docker compose run brt ./R/model_fit_brt.r /data/enhance/iccat/iccat-enhanced_1to1.csv /data/model-brt/iccat-config/', i, '/', i, '_config_brt.csv data/model-brt/iccat-fit/', i, '_fit_brt.RDS data/model-brt/iccat-fit/', i, '_fit_brt_eval.csv'))
#}
df <- data.table::fread("./data/enhance/iccat/iccat-enhanced_1to1.csv")
for (bb in sp){
  df.bb <- df %>% filter(SpeciesCode == bb) 
  data.table::fwrite(df.bb, paste0("./data/enhance/iccat/enhanced_sp/iccat-enhanced_1to1_", bb, '.csv'))
}

home_dir <- getwd()
brt_dir <- './pipelines/model-brt/R/'
setwd(brt_dir)
args <- list()
for (bb in sp){
  args$input_csv <- paste0("../../../data/enhance/iccat/enhanced_sp/iccat-enhanced_1to1_", bb, '.csv')
  args$config_file <- paste0("../../../data/model-brt/iccat-config/", bb, "/", bb, "_config_brt.csv")
  args$output_model <- paste0("../../../data/model-brt/iccat-fit/", bb, "_fit_brt.RDS")
  args$output_eval <- paste0("../../../data/model-brt/iccat-fit/", bb, "_fit_brt_eval.csv")
  system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))
}
setwd(home_dir)

# BSH specific iccat models

## read in data
df <- data.table::fread(paste0("./data/enhance/iccat/enhanced_sp/iccat-enhanced_1to1_BSH.csv"), sep=',', header=T)

## randomly sample to decrease N
df.1 <- df %>% filter(pres == 1)
df.0 <- df %>% filter(pres == 0)
set.seed(311)
df.1 <- df.1[sample(1:nrow(df.1), size = 5115),] ## sample size to match smallest e-tag dataset
df.0 <- df.0[sample(1:nrow(df.0), size = nrow(df.1)),]
df_1b <- rbind(df.1, df.0)
rm(df.1); rm(df.0)
data.table::fwrite(df_1b, "./data/enhance/iccat/enhanced_sp/iccat-enhanced_1to1_BSH_1b.csv")

## sample down to BOTH temporal limits and N
df_1c <- df %>% filter(as.Date(date) >= as.Date("2009-01-01") & as.Date(date) <= as.Date("2017-12-31")) 
data.table::fwrite(df_1c, "./data/enhance/iccat/enhanced_sp/iccat-enhanced_1to1_BSH_1c.csv")

bb='BSH'
setwd(brt_dir)
args$input_csv <- paste0("../../../data/enhance/iccat/enhanced_sp/iccat-enhanced_1to1_BSH_1b.csv")
args$config_file <- paste0("../../../data/model-brt/iccat-config/", bb, "/", bb, "_config_brt.csv")
args$output_model <- paste0("../../../data/model-brt/iccat-fit/", bb, "_fit_brt_1b.RDS")
args$output_eval <- paste0("../../../data/model-brt/iccat-fit/", bb, "_fit_brt_eval_1b.csv")
system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))

args$input_csv <- paste0("../../../data/enhance/iccat/enhanced_sp/iccat-enhanced_1to1_BSH_1c.csv")
args$config_file <- paste0("../../../data/model-brt/iccat-config/", bb, "/", bb, "_config_brt.csv")
args$output_model <- paste0("../../../data/model-brt/iccat-fit/", bb, "_fit_brt_1c.RDS")
args$output_eval <- paste0("../../../data/model-brt/iccat-fit/", bb, "_fit_brt_eval_1c.csv")
system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))
setwd(home_dir)

# eTUFF -> tracks --------------------------------------------------------------

## psat_track:
#  docker compose run etag ./R/etag_track.r /data/etag/etuff/psat/160424_2009_91068_eTUFF.txt /data/etag/track/psat/160424_2009_91068_eTUFF_track.csv

## sat_track:
#  docker compose run etag ./R/etag_track.r /data/etag/etuff/sat/160424_2006_gsmp00099_eTUFF.txt /data/etag/track/sat/160424_2006_gsmp00099_eTUFF_track.csv

## etag_combine:
system('docker-compose run py_tools ./combine.py /data/etag/track/psat/ /data/etag/combined_tags_psat.csv')
system('docker-compose run py_tools ./combine.py /data/etag/track/sat/ /data/etag/combined_tags_sat.csv')

# pseudoabs for etags --------------------------------------------------------------

## etag_crw:
system('docker compose run etag ./R/generate_pseudoabs.r /data/etag/combined_tags_psat.csv --index_var instrument_name --bathy_file /data/bathy/global_bathy_0.01.nc --n_sims 5 --force_180 --pseudo crw /data/etag/pseudoabs/with_pseudoabs_crw_psat.csv')
system('docker compose run etag ./R/generate_pseudoabs.r /data/etag/combined_tags_sat.csv --index_var instrument_name --bathy_file /data/bathy/global_bathy_0.01.nc --n_sims 5 --force_180 --pseudo crw /data/etag/pseudoabs/with_pseudoabs_crw_sat.csv')

## etag_bg:
system('docker compose run etag ./R/generate_pseudoabs.r /data/etag/combined_tags_psat.csv --index_var instrument_name --bathy_file /data/bathy/global_bathy_0.01.nc --abs_ratio 5 --force_180 --pseudo bg /data/etag/pseudoabs/with_pseudoabs_bg_psat.csv')
system('docker compose run etag ./R/generate_pseudoabs.r /data/etag/combined_tags_sat.csv --index_var instrument_name --bathy_file /data/bathy/global_bathy_0.01.nc --abs_ratio 5 --force_180 --pseudo bg /data/etag/pseudoabs/with_pseudoabs_bg_sat.csv')

## sample bg for sat tags using an extent rather than the spatial limits/extent of a given tag as this leads to not enough environmental niche separation between pres and abs
df <- data.table::fread('./data/etag/combined_tags_sat.csv')
system(paste0('docker compose run etag ./R/generate_pseudoabs.r /data/etag/combined_tags_sat.csv --index_var instrument_name --bathy_file /data/bathy/global_bathy_0.01.nc --abs_ratio 5 --force_180 --pseudo bg --extent ',
              '--xmin ', min(df$longitude),
              ' --xmax ', -5,
              ' --ymin ', 10,
              ' --ymax ', max(df$latitude),
              ' /data/etag/pseudoabs/with_pseudoabs_bg_sat_extent2.csv'))

df <- data.table::fread('./data/etag/combined_tags_psat.csv')
system(paste0('docker compose run etag ./R/generate_pseudoabs.r /data/etag/combined_tags_psat.csv --index_var instrument_name --bathy_file /data/bathy/global_bathy_0.01.nc --abs_ratio 5 --force_180 --pseudo bg --extent ',
              '--xmin ', min(df$longitude),
              ' --xmax ', -5,
              ' --ymin ', 10,
              ' --ymax ', max(df$latitude),
              ' /data/etag/pseudoabs/with_pseudoabs_bg_psat_extent2.csv'))

## split
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_bg_sat.csv /data/enhance/etag/split/bg_sat/')
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_bg_psat.csv /data/enhance/etag/split/bg_psat/')
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_crw_sat.csv /data/enhance/etag/split/crw_sat/')
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_crw_psat.csv /data/enhance/etag/split/crw_psat/')
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_bg_sat_extent.csv /data/enhance/etag/split/bg_sat_extent/')
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_bg_psat_extent.csv /data/enhance/etag/split/bg_psat_extent/')
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_bg_sat_extent2.csv /data/enhance/etag/split/bg_sat_extent2/')
system('docker compose run enhance ./split_by_date.py /data/etag/pseudoabs/with_pseudoabs_bg_psat_extent2.csv /data/enhance/etag/split/bg_psat_extent2/')

# enhance etag --------------------------------------------------------------
## enhance etag w/ cmems daily
etag_vec <- c('bg_sat','bg_psat','crw_sat','crw_psat', 'bg_sat_extent', 'bg_psat_extent', 'bg_psat_extent_error')
for (bb in etag_vec){
  fList <- list.files(paste0('./data/enhance/etag/split/', bb, '/'), recursive = TRUE, full.names = TRUE)
  #cmems_files <- list.files('./data/glorys_daily/', recursive = TRUE, full.names = TRUE)
  for (i in 1:length(fList)){
    df <- data.table::fread(fList[i])
    raster_name <- cmems_files[grep(substr(fList[i], nchar(fList[i]) - 13, nchar(fList[i]) - 4), cmems_files)]
    if (length(raster_name) == 0) next
    raster_name <- raster_name[grep('.grd', raster_name)]
    if (length(grep('_psat', bb)) == 0) extr <- data.frame(raster::extract(raster::brick(raster_name), cbind(df$lon, df$lat))) ## pointwise extract
    if (length(grep('_psat', bb)) > 0) extr <- data.frame(mean(raster::extract(raster::brick(raster_name), extent(df$lon - df$longitudeError,
                                                                                                                  df$lon + df$longitudeError,
                                                                                                                  df$lat - df$latitudeError,
                                                                                                                  df$lat + df$latitudeError)), na.rm = TRUE))
    names(extr) <- c('sst','sss','ssh','mld', 'log_eke','sst_sd','ssh_sd','sss_sd','bathy','rugosity')
    df <- cbind(df, extr)
    out_name <- paste0('~/Documents/work/RCode/NASA-FaCeT/data/enhance/etag/enhanced/', bb, '/', substr(fList[i], nchar(fList[i]) - 13, nchar(fList[i])))
    data.table::fwrite(df, out_name)
    print(out_name)
  }
}

## recombine all the daily enhanced data to master csv
#for (bb in etag_vec) system(paste0('docker compose run enhance ./combine.py /data/enhance/etag/enhanced/', bb, '/ /data/enhance/etag/enhanced/', bb, '-enhanced.csv'))
for (bb in etag_vec){
  fList <- list.files(paste0('~/Documents/work/RCode/NASA-FaCeT/data/enhance/etag/enhanced/', bb, '/'), recursive = TRUE, full.names = TRUE)
  fList <- fList[grep('.csv', fList)]
  output_csv <- paste0('~/Documents/work/RCode/NASA-FaCeT/data/enhance/etag/etag-enhanced-', bb, '.csv')
  for (i in 1:length(fList)){
    df <- data.table::fread(fList[i])
    if (file.exists(output_csv)){
      data.table::fwrite(df, output_csv, append=TRUE)
    } else{
      data.table::fwrite(df, output_csv, append=FALSE)
    }
  }
}

## sample the data to a 1:1 ratio before model fit
etag_vec <- c('bg_sat_extent2', 'bg_psat_extent2')
for (bb in etag_vec){
  df <- data.table::fread(paste0('./data/enhance/etag/enhanced/etag-enhanced-', bb, '.csv'))
  df <- data.frame(na.omit(df, cols = c(which(names(df) == 'sst'):ncol(df))))
  df <- df %>% filter(as.Date(date) > as.Date('1993-01-01') & as.Date(date) <= as.Date('2019-12-31') & 
                        lon >= xl[1] & lon <= xl[2] & lat >= yl[1] & lat <= yl[2])
  df.split <- split(df, df$platform)
  df <- lapply(df.split, FUN = function(x){
    df.1 <- x %>% filter(pres == 1)
    df.0 <- x %>% filter(pres == 0)
    set.seed(311)
    df.0 <- df.0[sample(1:nrow(df.0), size = nrow(df.1)),]
    df <- rbind(df.1, df.0)
  }) %>% do.call(rbind, .)
  df <- df[order(df$date),]
  data.table::fwrite(df, file = paste0('./data/enhance/etag/enhanced/', bb, '-enhanced_1to1.csv'))
}

## filter SAT data to remove error prone positions, then sample to 1:1 ratio before model fit
df <- data.table::fread(paste0('./data/enhance/etag/enhanced/bg_sat_extent-enhanced_fixedMedianError.csv'))
df <- data.frame(na.omit(df, cols = c(which(names(df) == 'sst'):ncol(df))))
df <- df %>% filter(as.Date(date) > as.Date('1993-01-01') & as.Date(date) <= as.Date('2019-12-31') & lon >= xl[1] & lon <= xl[2] & lat >= yl[1] & lat <= yl[2])
df.split <- split(df, df$platform)
df <- lapply(df.split, FUN = function(x){
  df.1 <- x %>% filter(pres == 1)# & longitudeError < 0.25 & latitudeError < 0.25)
  df.0 <- x %>% filter(pres == 0)
  set.seed(311)
  df.0 <- df.0[sample(1:nrow(df.0), size = nrow(df.1)),]
  df <- rbind(df.1, df.0)
}) %>% do.call(rbind, .)
df <- df[order(df$date),]
data.table::fwrite(df, file = paste0('./data/enhance/etag/enhanced/bg_sat_extent_fixedMedianError-enhanced_1to1.csv'))


# ETAG - fit BRTs --------------------------------------------------------------
#for (bb in etag_vec){
#  system(paste0('docker compose run brt ./R/model_fit_brt.r /data/enhance/etag/enhanced/', bb, '-enhanced_1to1.csv /data/model-brt/etag-config/', bb, '/', bb, '_config_brt.csv data/model-brt/etag-fit/', bb, '_fit_brt.RDS data/model-brt/etag-fit/', bb, '_fit_brt_eval.csv'))
#}
home_dir <- getwd()
brt_dir <- './pipelines/model-brt/R/'
#system('chmod +x ./pipelines/model-brt/R/model_fit_brt.r')
#bg_sat_extent-enhanced_fixedMedianError
setwd(brt_dir)
args <- list()
for (bb in etag_vec){
  args$input_csv <- paste0("../../../data/enhance/etag/enhanced/", bb, "-enhanced_1to1.csv")
  args$config_file <- paste0("../../../data/model-brt/etag-config/etag_config_brt.csv")
  args$output_model <- paste0("../../../data/model-brt/etag-fit/", bb, "_fit_brt_error_1to1.RDS")
  args$output_eval <- paste0("../../../data/model-brt/etag-fit/", bb, "_fit_brt_eval_error_1to1.csv")
  system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))
}
setwd(home_dir)



# OBSERVER - split & enhance --------------------------------------------------------------
obs <- data.table::fread('./data/observer/bsh_catch_new_1992_2019_braun.csv')
obs$date <- obs$END_HAUL_DT_UTC
obs$lon <- obs$HAUL_LONGITUDE
obs$lat <- obs$HAUL_LATITUDE
obs <- obs %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2] & date >= as.POSIXct('1993-01-01') & date <= as.POSIXct('2019-12-31'))
data.table::fwrite(obs, './data/observer/bsh_pop.csv')

## enhance
#setwd('~/Google Drive/Shared drives/MPG_WHOI/data/facet/')
fList <- list.files('./data/enhance/observer/split/', recursive = TRUE, full.names = TRUE)
cmems_files <- list.files('~/work/EnvData/glorys_daily/', recursive = TRUE, full.names = TRUE)
drop_names <- c('BATHYMETRY','RUGOSITY','SST','SSH','MLD','UO','VO','SSTSD')
for (i in 1:length(fList)){
  df <- data.table::fread(fList[i])
  df <- df %>% dplyr::select(!(all_of(drop_names)))
  raster_name <- cmems_files[grep(substr(fList[i], 32, nchar(fList[i]) - 4), cmems_files)]
  if (length(raster_name) == 0) next
  raster_name <- raster_name[grep('.grd', raster_name)]
  extr <- data.frame(raster::extract(raster::brick(raster_name), cbind(df$lon, df$lat)))
  names(extr) <- c('sst','sss','ssh','mld', 'log_eke','sst_sd','ssh_sd','sss_sd','bathy','rugosity')
  df <- cbind(df, extr)
  out_name <- paste0('./data/enhance/observer/enhanced_daily/', substr(fList[i], 32, nchar(fList[i])))
  data.table::fwrite(df, out_name)
  rm(df); rm(raster_name); rm(out_name)
  print(substr(fList[i], 32, nchar(fList[i])))
}

## recombine
input_dir <- '~/Documents/work/RCode/NASA-FaCeT/data/enhance/observer/enhanced_daily/'
output_csv <- '~/Documents/work/RCode/NASA-FaCeT/data/enhance/observer/observer-enhanced.csv'
fList <- list.files(input_dir, full.names = TRUE)
fList <- fList[grep('.csv', fList)]
for (i in 1:length(fList)){
  
  df <- data.table::fread(fList[i])
  if (file.exists(output_csv)){
    data.table::fwrite(df, output_csv, append=TRUE)
  } else{
    data.table::fwrite(df, output_csv, append=FALSE)
  }
}

#-------------------
## same as above but with generating bg pseudoabsences

obs <- obs %>% filter(CATCH != 0)
data.table::fwrite(obs, './data/observer/bsh_pop_presOnly.csv')
dir.create('./data/observer/pseudoabs/', recursive = TRUE)

system(paste0('docker compose run etag ./R/generate_pseudoabs.r /data/observer/bsh_pop_presOnly.csv --bathy_file /data/bathy/global_bathy_0.01.nc --abs_ratio 5 --force_180 --pseudo bg --extent ',
              '--xmin ', xl[1],
              ' --xmax ', -39,
              ' --ymin ', yl[1],
              ' --ymax ', 51,
              ' /data/observer/pseudoabs/with_pseudoabs_bg.csv'))

dir.create('./data/enhance/observer/split_pseudo/', recursive = TRUE)
system('docker compose run enhance ./split_by_date.py /data/observer/pseudoabs/with_pseudoabs_bg.csv /data/enhance/observer/split_pseudo/')

## enhance
#setwd('~/Google Drive/Shared drives/MPG_WHOI/data/facet/')
fList <- list.files('./data/enhance/observer/split_pseudo/', recursive = TRUE, full.names = TRUE)
cmems_files <- list.files('~/work/EnvData/glorys_daily/', recursive = TRUE, full.names = TRUE)
drop_names <- c('BATHYMETRY','RUGOSITY','SST','SSH','MLD','UO','VO','SSTSD')
for (i in 1:length(fList)){
  df <- data.table::fread(fList[i])
  df <- df %>% dplyr::select(!(all_of(drop_names)))
  raster_name <- cmems_files[grep(substr(fList[i], 32, nchar(fList[i]) - 4), cmems_files)]
  if (length(raster_name) == 0) next
  raster_name <- raster_name[grep('.grd', raster_name)]
  extr <- data.frame(raster::extract(raster::brick(raster_name), cbind(df$lon, df$lat)))
  names(extr) <- c('sst','sss','ssh','mld', 'log_eke','sst_sd','ssh_sd','sss_sd','bathy','rugosity')
  df <- cbind(df, extr)
  out_name <- paste0('./data/enhance/observer/enhanced_daily_pseudo/', substr(fList[i], 32, nchar(fList[i])))
  data.table::fwrite(df, out_name)
  rm(df); rm(raster_name); rm(out_name)
  print(substr(fList[i], 32, nchar(fList[i])))
}

## recombine
input_dir <- '~/Documents/work/RCode/NASA-FaCeT/data/enhance/observer/enhanced_daily_pseudo/'
output_csv <- '~/Documents/work/RCode/NASA-FaCeT/data/enhance/observer/observer-enhanced_pseudo.csv'
fList <- list.files(input_dir, full.names = TRUE)
fList <- fList[grep('.csv', fList)]
for (i in 1:length(fList)){
  df <- data.table::fread(fList[i])
  if (file.exists(output_csv)){
    data.table::fwrite(df, output_csv, append=TRUE)
  } else{
    data.table::fwrite(df, output_csv, append=FALSE)
  }
}


# OBSERVER - subsample to 1:1 --------------------------------------------------------------

## sample the data to a 1:1 ratio before model fit
obs <- data.table::fread('./data/enhance/observer/observer-enhanced.csv')
obs <- data.frame(na.omit(obs, cols = c(which(names(obs) == 'sst'):ncol(obs))))
obs$pres <- ifelse(obs$CATCH == 0, 0, 1)
data.table::fwrite(obs, file = './data/enhance/observer/observer-enhanced_allAbs.csv')
obs.1 <- obs %>% filter(pres == 1)
obs.0 <- obs %>% filter(pres == 0)
set.seed(311)
obs.0 <- obs.0[sample(1:nrow(obs.0), size = nrow(obs.1)),]
obs <- rbind(obs.1, obs.0)
obs <- obs[order(obs$date),]
data.table::fwrite(obs, file = './data/enhance/observer/observer-enhanced_1to1.csv')

## sample the data to a 1:1 ratio before model fit
obs <- data.table::fread('./data/enhance/observer/observer-enhanced_pseudo.csv')
obs <- data.frame(na.omit(obs, cols = c(which(names(obs) == 'sst'):ncol(obs))))
#obs$pres <- ifelse(obs$CATCH == 0, 0, 1)
#data.table::fwrite(obs, file = './data/enhance/observer/observer-enhanced_allAbs.csv')
obs.1 <- obs %>% filter(pres == 1)
obs.0 <- obs %>% filter(pres == 0)
set.seed(311)
obs.0 <- obs.0[sample(1:nrow(obs.0), size = nrow(obs.1)),]
obs <- rbind(obs.1, obs.0)
obs <- obs[order(obs$date),]
data.table::fwrite(obs, file = './data/enhance/observer/observer-enhanced_pseudo_1to1.csv')

# OBSERVER - fit BRTs --------------------------------------------------------------
#for (bb in etag_vec){
#  system(paste0('docker compose run brt ./R/model_fit_brt.r /data/enhance/etag/enhanced/', bb, '-enhanced_1to1.csv /data/model-brt/etag-config/', bb, '/', bb, '_config_brt.csv data/model-brt/etag-fit/', bb, '_fit_brt.RDS data/model-brt/etag-fit/', bb, '_fit_brt_eval.csv'))
#}
home_dir <- getwd()
brt_dir <- './pipelines/model-brt/R/'
#system('chmod +x ./pipelines/model-brt/R/model_fit_brt.r')
setwd(brt_dir)
args <- list()
args$input_csv <- paste0("../../../data/enhance/observer/observer-enhanced_allAbs.csv")
args$config_file <- paste0("../../../data/model-brt/observer-config/observer_config_brt.csv")
args$output_model <- paste0("../../../data/model-brt/observer-fit/observer_fit_brt_allAbs.RDS")
args$output_eval <- paste0("../../../data/model-brt/observer-fit/observer_fit_brt_eval_allAbs.csv")
system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))

args$input_csv <- paste0("../../../data/enhance/observer/observer-enhanced_1to1.csv")
args$config_file <- paste0("../../../data/model-brt/observer-config/observer_config_brt.csv")
args$output_model <- paste0("../../../data/model-brt/observer-fit/observer_fit_brt_1to1.RDS")
args$output_eval <- paste0("../../../data/model-brt/observer-fit/observer_fit_brt_eval_1to1.csv")
system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))

args$input_csv <- paste0("../../../data/enhance/observer/observer-enhanced_pseudo_1to1.csv")
args$config_file <- paste0("../../../data/model-brt/observer-config/observer_config_brt.csv")
args$output_model <- paste0("../../../data/model-brt/observer-fit/observer_fit_brt_pseudo.RDS")
args$output_eval <- paste0("../../../data/model-brt/observer-fit/observer_fit_brt_eval_pseudo.csv")
system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))


# MODEL SELECTION for indiv data types --------------------------------------------------------------

fList <- list.files('./data/model-brt/', recursive = TRUE, full.names = TRUE)
fList <- fList[grep('fit', fList)]
fList <- fList[grep('.csv', fList)]

for (i in 1:length(fList)){
  eval <- data.table::fread(fList[i])
  eval$name <- fList[i]
  data.table::fwrite(eval, './data/model-brt/combine_eval.csv', append=TRUE)
}

# ENSEMBLE --------------------------------------------------------------

## combine all the pres/abs data and fit one overall model
## AND
## ensemble models fit individually

## predict to each presence dataset and the observer absence data with:
# - brt_5b
# - rowMeans(brt_1x, brt_2x, ...)

## calc metrics from predictions:
# - deviance expl
# - sensitivity
# - specificity
# - TSS
# - AUC


## FIT MODELS TO ALL COMBINED DATA
## combine the 1:1 datasets across data types
var_names <- c('sst','sss','ssh','mld', 'log_eke','sst_sd','ssh_sd','sss_sd','bathy','rugosity')            
df1 <- data.table::fread('./data/enhance/iccat/iccat-enhanced_1to1.csv')
df1 <- df1 %>% filter(SpeciesCode == 'BSH') %>% dplyr::select(lon, lat, date, pres, all_of(var_names))
df1$type <- 'iccat'
df2d <- data.table::fread("./data/enhance/etag/enhanced/bg_sat_extent_error-enhanced_1to1.csv")
df2d <- df2d %>% dplyr::select(lon, lat, date, pres, all_of(var_names))
df2d$type <- 'sat'
df3d <- data.table::fread('./data/enhance/etag/enhanced/bg_psat_extent_ERROR-enhanced_1to1.csv')
df3d <- df3d %>% dplyr::select(lon, lat, date, pres, all_of(var_names))
df3d$type <- 'psat'
obs <- data.table::fread(file = './data/enhance/observer/observer-enhanced_pseudo_1to1.csv')
obs <- obs %>% dplyr::select(lon, lat, date, pres, all_of(var_names))
obs$type <- 'pop'
obs$date <- as.POSIXct(obs$date, tz='UTC')
all_df <- rbind(df1, df2d, df3d, obs)
data.table::fwrite(all_df, './data/enhance/all_data_pseudo_1to1.csv')


## subsample so N matches across data types (i.e. "equal" treatment of each data type)
df.split <- split(all_df, all_df$type)
all_df_sub <- lapply(df.split, FUN = function(x){
  df.1 <- x %>% filter(pres == 1)
  if (nrow(df.1) > 4913) df.1 <- df.1[sample(1:nrow(df.1), size = 4913),]
  df.0 <- x %>% filter(pres == 0)
  set.seed(311)
  df.0 <- df.0[sample(1:nrow(df.0), size = nrow(df.1)),]
  df <- rbind(df.1, df.0)
}) %>% do.call(rbind, .)
data.table::fwrite(all_df_sub, './data/enhance/all_data_equalN_1to1.csv')

setwd(brt_dir)
args <- list()
args$input_csv <- paste0("../../../data/enhance/all_data_pseudo_1to1.csv")
args$config_file <- paste0("../../../data/model-brt/observer-config/observer_config_brt.csv")
args$output_model <- paste0("../../../data/model-brt/all_data_fit_brt_pseudo.RDS")
args$output_eval <- paste0("../../../data/model-brt/all_data_fit_brt_pseudo_eval.csv")
system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))

args$input_csv <- paste0("../../../data/enhance/all_data_equalN_1to1.csv")
args$config_file <- paste0("../../../data/model-brt/observer-config/observer_config_brt.csv")
args$output_model <- paste0("../../../data/model-brt/all_data_equalN_fit_brt.RDS")
args$output_eval <- paste0("../../../data/model-brt/all_data_equalN_fit_brt_eval.csv")
system(paste0('./model_fit_brt.r ', args$input_csv, ' ', args$config_file, ' ', args$output_model, ' ', args$output_eval))
setwd(home_dir)


## ENSEMBLE ACROSS SELECTED MODELS

brt_1a <- readRDS('./data/model-brt/iccat-fit/BSH_fit_brt.RDS')
#brt_1b <- readRDS('./data/model-brt/iccat-fit/BSH_fit_brt_1b.RDS')
#brt_1c <- readRDS('./data/model-brt/iccat-fit/BSH_fit_brt_1c.RDS')
#brt_2a <- readRDS('./data/model-brt/etag-fit/crw_sat_fit_brt_1to1.RDS')
#brt_2b <- readRDS('./data/model-brt/etag-fit/bg_sat_fit_brt_1to1.RDS')
#brt_2c <- readRDS('./data/model-brt/etag-fit/bg_sat_extent_fit_brt_error_1to1.RDS')
brt_2d <- readRDS('./data/model-brt/etag-fit/bg_sat_extent_error_fit_brt_error_1to1.RDS')
#brt_3a <- readRDS('./data/model-brt/etag-fit/crw_psat_fit_brt_1to1.RDS')
#brt_3b <- readRDS('./data/model-brt/etag-fit/bg_psat_fit_brt_1to1.RDS') ##
#brt_3c <- readRDS('./data/model-brt/etag-fit/bg_psat_extent_fit_brt_1to1.RDS') ##
brt_3d <- readRDS('./data/model-brt/etag-fit/bg_psat_extent_fit_brt_error_1to1.RDS') ##
#brt_4a <- readRDS('./data/model-brt/observer-fit/observer_fit_brt_eval_allAbs.RDS')
#brt_4b <- readRDS('./data/model-brt/observer-fit/observer_fit_brt_1to1.RDS')
brt_4b <- readRDS('./data/model-brt/observer-fit/observer_fit_brt_pseudo.RDS')
brt_5a <- readRDS('./data/model-brt/all_data_fit_brt_pseudo.RDS')
#brt_5b <- readRDS('./data/model-brt/all_data_equalN_fit_brt.RDS')

models <- list(brt_1a, brt_2d, brt_3d, brt_4b, brt_5a)#, 'brt_ens')
models[[6]] <- NA
names(models) <- c('brt_1a', 'brt_2d', 'brt_3d', 'brt_4b', 'brt_5a', 'brt_ens')
all_df <- data.table::fread('./data/enhance/all_data_pseudo_1to1.csv')
datasets <- split(all_df, all_df$type)
datasets[[5]] <- all_df
names(datasets)[5] <- 'all_df'

res_mat <- list()
#res_mat <- as.data.frame(matrix(NA, ncol=6, nrow=length(datasets) * 3))
#colnames(res_mat) <- c('model','data','auc','tss','sens','spec')
counter <- 1
for (i in 1:length(models)){
  for (b in 1:length(datasets)){
    
    if (i == length(models)){
      ## ensemble
      ens_preds <- list()
      for (tt in 1:(length(models) - 1)){
        ens_preds[[tt]] <- gbm::predict.gbm(models[[tt]], datasets[[b]],
                                           n.trees = 2000,
                                           type = 'response',
                                           na.rm = FALSE)
        
      }
      pred <- rowMeans(as.data.frame(ens_preds, col.names=seq(1:4)))
      
      ## presence only
      ens_preds <- list()
      for (tt in 1:(length(models) - 1)){
        ens_preds[[tt]] <- gbm::predict.gbm(models[[tt]], datasets[[b]][which(datasets[[b]]$pres == 1),],
                                            n.trees = 2000,
                                            type = 'response',
                                            na.rm = FALSE)
        
      }
      pred.1 <- rowMeans(as.data.frame(ens_preds, col.names=seq(1:4)))
      
      ## presence only
      ens_preds <- list()
      for (tt in 1:(length(models) - 1)){
        ens_preds[[tt]] <- gbm::predict.gbm(models[[tt]], datasets[[b]][which(datasets[[b]]$pres == 0),],
                                            n.trees = 2000,
                                            type = 'response',
                                            na.rm = FALSE)
        
      }
      pred.0 <- rowMeans(as.data.frame(ens_preds, col.names=seq(1:4)))
      
      
    } else{
        pred <- gbm::predict.gbm(models[[i]], datasets[[b]],
                             n.trees = 2000,
                             type = 'response',
                             na.rm = FALSE)
    }
    
    #evaluate
    pred_cm <- pred
    pred_cm[which(pred_cm <= 0.25)] <- 0
    pred_cm[which(pred_cm >= 0.75)] <- 1
    pred_cm[which(pred_cm > 0.25 & pred_cm < 0.75)] <- NA
    idx <- which(is.na(pred_cm))
    pred_cm <- pred_cm[-idx]
    cm <- caret::confusionMatrix(factor(pred_cm),
                                 factor(datasets[[b]]$pres[-idx]), positive = '1')
    d <- cbind(datasets[[b]]$pres, pred)
    pres <- as.numeric(d[d[,1] == 1,2])
    abs <- as.numeric(d[d[,1] == 0,2])
    e <- dismo::evaluate(p=pres, a=abs)
    auc <- round(e@auc, 2)
    tss <- round(max(e@TPR + e@TNR-1), 2)
    sens <- round(cm$byClass['Sensitivity'], 2)
    spec <- round(cm$byClass['Specificity'], 2)
    acc <- cm$overall['Accuracy']
    acc_p <- cm$overall['AccuracyPValue']
    
    pred_q <- pred
    pred_q[which(pred_q <= quantile(pred)[2])] <- 0
    pred_q[which(pred_q >= quantile(pred)[4])] <- 1
    pred_q[which(pred_q > quantile(pred)[2] & pred_q < quantile(pred)[4])] <- NA
    idx <- which(is.na(pred_q))
    pred_q <- pred_q[-idx]
    cm <- caret::confusionMatrix(factor(pred_q),
                                 factor(datasets[[b]]$pres[-idx]), positive = '1')
    d <- cbind(datasets[[b]]$pres, pred)
    pres <- as.numeric(d[d[,1] == 1,2])
    abs <- as.numeric(d[d[,1] == 0,2])
    e <- dismo::evaluate(p=pres, a=abs)
    sens_q <- round(cm$byClass['Sensitivity'], 2)
    spec_q <- round(cm$byClass['Specificity'], 2)
    accq <- cm$overall['Accuracy']
    accq_p <- cm$overall['AccuracyPValue']
    
    if (i == length(models)){
      median_pred_at_pres <- median(pred.1)
      median_pred_at_abs <- median(pred.0)
    } else {
      median_pred_at_pres <- median(gbm::predict.gbm(models[[i]], datasets[[b]][which(datasets[[b]]$pres == 1),], n.trees = models[[i]]$gbm.call$best.trees, type="response"))
      median_pred_at_abs <- median(gbm::predict.gbm(models[[i]], datasets[[b]][which(datasets[[b]]$pres == 0),], n.trees = models[[i]]$gbm.call$best.trees, type="response"))
    }
    
    res_mat[[counter]] <- list(names(models)[i], names(datasets)[b], auc, tss, sens, spec, acc, acc_p, q25=quantile(pred)[2], q75=quantile(pred)[4], sens_q, spec_q, accq, accq_p)
    counter <- counter + 1
    rm(auc); rm(pred); rm(tss); rm(sens); rm(spec)
  }
}

res_mat <- do.call(rbind.data.frame, res_mat)
colnames(res_mat) <- c('model','data','auc','tss','sens','spec', 'accuracy', 'accuracy_pvalue', 'q25','q75','sens_quantile','spec_quantile', 'accuracy_quantile', 'accuracy_pvalue_quantile')
data.table::fwrite(res_mat, file='./data/bsh_resmat_20221223.csv')


for (i in 1:length(models)){
  for (b in 1:length(datasets)){
    
    median_pred_at_pres <- median(gbm::predict.gbm(models[[i]], datasets[[b]][which(datasets[[b]]$pres == 1),], n.trees = models[[i]]$gbm.call$best.trees, type="response"))
    median_pred_at_abs <- median(gbm::predict.gbm(models[[i]], datasets[[b]][which(datasets[[b]]$pres == 0),], n.trees = models[[i]]$gbm.call$best.trees, type="response"))
    
    res_mat[[counter]] <- list(names(models)[i], names(datasets)[b], median_pred_at_pres, median_pred_at_abs)
    counter <- counter + 1
    
  }
}


# METRICS & FIGURES --------------------------------------------------------------


# FIGURE 1 - plot distribution of data --------------------------------------------------------------

## get world map data
world <- map_data('world')

df1 <- data.table::fread('./data/enhance/iccat/iccat-enhanced_1to1.csv')
df1 <- df1 %>% filter(SpeciesCode == 'BSH')
p1 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill='grey') +
  coord_fixed(xlim=xl, ylim=yl, ratio=1.3) + 
  xlab('') + ylab('') +# ggtitle(paste('ICCAT', ctag.new$SpeciesCode[1], 'conventional tags N =', nrow(ctag.new))) +
  geom_point(data=df1[which(df1$pres == 0),], aes(x=lon, y=lat), colour=rgb(239/255,138/255,98/255, alpha=0.25)) +
  geom_point(data=df1[which(df1$pres == 1),], aes(x=lon, y=lat), colour=rgb(103/255, 169/255, 207/255, alpha=0.75)) +
  theme_bw(base_size = 10) + theme(panel.grid=element_blank(), panel.border = element_blank())  + theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"))
#p1

#df2 <- data.table::fread("./data/enhance/etag/enhanced/crw_sat-enhanced_1to1.csv")
#df2b <- data.table::fread("./data/enhance/etag/enhanced/bg_sat-enhanced_1to1.csv")
bb=etag_vec[7]
#df2c <- data.table::fread(paste0("./data/enhance/etag/enhanced/", bb, "-enhanced_1to1.csv"))
#df2c <- data.table::fread("./data/enhance/etag/enhanced/bg_sat_extent-enhanced_1to1.csv")
sat_start <- data.table::fread('./data/etag/combined_tags_sat.csv') %>% group_by(instrument_name) %>% slice(1) #%>% filter(lon >= xl[1] & lon <= xl[2] & lat >= yl[1] & lat <= yl[2])
p2 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill='grey') +
  coord_fixed(xlim=xl, ylim=yl, ratio=1.3) + 
  xlab('') + ylab('') +# ggtitle(paste('ICCAT', ctag.new$SpeciesCode[1], 'conventional tags N =', nrow(ctag.new))) +
  geom_point(data=df2c[which(df2c$pres == 0),], aes(x=lon, y=lat), colour=rgb(239/255,138/255,98/255, alpha=0.25)) +
  geom_point(data=df2c[which(df2c$pres == 1),], aes(x=lon, y=lat), colour=rgb(103/255, 169/255, 207/255, alpha=0.75)) +
  geom_point(data=sat_start, aes(x=longitude, y=latitude), shape=24, fill='green') +
  theme_bw(base_size = 10) + theme(panel.grid=element_blank(), panel.border = element_blank())  + theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"))
p2

#df3 <- data.table::fread("./data/enhance/etag/enhanced/crw_psat-enhanced_1to1.csv")
#df3b <- data.table::fread("./data/enhance/etag/enhanced/bg_psat-enhanced_1to1.csv")
#df3c <- data.table::fread("./data/enhance/etag/enhanced/bg_psat_extent-enhanced_1to1.csv")
df3d <- data.table::fread('./data/enhance/etag/enhanced/bg_psat_extent_error-enhanced_1to1.csv')
psat_start <- data.table::fread('./data/etag/combined_tags_psat.csv') %>% group_by(instrument_name) %>% slice(1) #%>% filter(lon >= xl[1] & lon <= xl[2] & lat >= yl[1] & lat <= yl[2])
p3 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill='grey') +
  coord_fixed(xlim=xl, ylim=yl, ratio=1.3) + 
  xlab('') + ylab('') +# ggtitle(paste('ICCAT', ctag.new$SpeciesCode[1], 'conventional tags N =', nrow(ctag.new))) +
  geom_point(data=df3d[which(df3d$pres == 0),], aes(x=lon, y=lat), colour=rgb(239/255,138/255,98/255, alpha=0.25)) +
  geom_point(data=df3d[which(df3d$pres == 1),], aes(x=lon, y=lat), colour=rgb(103/255, 169/255, 207/255, alpha=0.75)) +
  geom_point(data=psat_start, aes(x=longitude, y=latitude), shape=24, fill='green') +
  theme_bw(base_size = 10) + theme(panel.grid=element_blank(), panel.border = element_blank()) + theme(plot.margin = unit(c(.5,.5,.5,.5), "cm")) ## top, right, bottom, left
#p3

lay <- rbind(c(1),
             c(2),
             c(3))
g <- gridExtra::arrangeGrob(grobs = list(p1, p2, p3), heights = c(4,4,4),
                            width = c(6), layout_matrix = lay)
ggsave(file = './figures/biocomp_fig1_v3.pdf', width=6, height=10, units = 'in', g)



## just to see what looks like if plotted like observer data as counts per cell
df1a.1 <- df1 %>% filter(pres == 1)
df1a.0 <- df1 %>% filter(pres == 0)
df1a_pres_ras <- rasterize(cbind(df1a.1$lon, df1a.1$lat), r, fun='count')
df1a_abs_ras <- rasterize(cbind(df1a.0$lon, df1a.0$lat), r, fun='count')
pp1a <- plot_observer(df1a_pres_ras)
pp1a_abs <- plot_observer(df1a_abs_ras)

df2d <- data.table::fread("./data/enhance/etag/enhanced/bg_sat_extent2-enhanced_1to1.csv")
df2d.1 <- df2d %>% filter(pres == 1)
df2d.0 <- df2d %>% filter(pres == 0)
df2d_pres_ras <- rasterize(cbind(df2d.1$lon, df2d.1$lat), r, fun='count')
df2d_abs_ras <- rasterize(cbind(df2d.0$lon, df2d.0$lat), r, fun='count')
pp2d <- plot_observer(df2d_pres_ras)
pp2d <- pp2d + geom_point(data=sat_start, aes(x=longitude, y=latitude), shape=24, fill='green')
pp2d_abs <- plot_observer(df2d_abs_ras)

df3d <- data.table::fread("./data/enhance/etag/enhanced/bg_psat_extent2-enhanced_1to1.csv")
df3d.1 <- df3d %>% filter(pres == 1)
df3d.0 <- df3d %>% filter(pres == 0)
df3d_pres_ras <- rasterize(cbind(df3d.1$lon, df3d.1$lat), r, fun='count')
df3d_abs_ras <- rasterize(cbind(df3d.0$lon, df3d.0$lat), r, fun='count')
pp3d <- plot_observer(df3d_pres_ras)
pp3d <- pp3d + geom_point(data=psat_start, aes(x=longitude, y=latitude), shape=24, fill='green')
pp3d_abs <- plot_observer(df3d_abs_ras)


lay <- rbind(c(1, 2),
             c(3, 4),
             c(5, 6),
             c(7, 8))
g <- gridExtra::arrangeGrob(grobs = list(pp1a, pp1a_abs, pp2d, pp2d_abs, pp3d, pp3d_abs, po.1, po.0), heights = c(5, 5, 5, 5),
                            width = c(6, 6), layout_matrix = lay)
ggsave(file = './figures/all_data_fig2.pdf', width=14, height=22, units = 'in', g)


## observer data figure

## CANNOT plot observer data pointwise as below
#p4 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill='grey') +
#  coord_fixed(xlim=xl, ylim=yl, ratio=1.3) + 
#  xlab('') + ylab('') +# ggtitle(paste('ICCAT', ctag.new$SpeciesCode[1], 'conventional tags N =', nrow(ctag.new))) +
#  geom_point(data=obs[which(obs$pres == 0),], aes(x=lon, y=lat), colour=rgb(239/255,138/255,98/255, alpha=0.25)) +
#  geom_point(data=obs[which(obs$pres == 1),], aes(x=lon, y=lat), colour=rgb(103/255, 169/255, 207/255, alpha=0.75)) +
#  theme_bw(base_size = 0) + theme(panel.grid=element_blank(), panel.border = element_blank()) + theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"))
#p4


obs <- data.table::fread('./data/enhance/observer/observer-enhanced_allAbs.csv')
r <- raster(xmn=xl[1], ymn=yl[1], xmx=xl[2], ymx=yl[2], res=.1)
obs.1 <- obs %>% filter(pres == 1)
obs.0 <- obs %>% filter(pres == 0)
obs_pres_mask <- rasterize(cbind(obs.1$lon, obs.1$lat), r, field = obs.1$POP_VESSEL_CODE, fun=function(x,...) length(unique(x)))
obs_pres_mask[obs_pres_mask <= 3] <- NA
obs_abs_mask <- rasterize(cbind(obs.0$lon, obs.0$lat), r, field = obs.0$POP_VESSEL_CODE, fun=function(x,...) length(unique(x)))
obs_abs_mask[obs_abs_mask <= 3] <- NA
obs_pres_ras <- rasterize(cbind(obs.1$lon, obs.1$lat), r, fun='count')
obs_abs_ras <- rasterize(cbind(obs.0$lon, obs.0$lat), r, fun='count')
#par(mfrow=c(2,1))
#plot(obs_pres_ras)
#plot(mask(obs_pres_ras, obs_pres_mask))

plot_observer <- function(r){
  
  df_map = rasterToPoints(r) %>% as.data.frame()
  colnames(df_map)= c("rows","cols","value")    
  
  map.world = map_data(map="world")
  testt=map.world %>% filter(long<=180)
  
  a=df_map %>% 
    ggplot() + 
    # geom_tile(aes(x = rows, y = cols, fill = ntile(value,100))) + coord_equal()+
    geom_tile(aes(x = rows, y = cols, fill = value), color=NA) + #coord_equal()+
    #scale_fill_gradientn("Habitat suitability",colours = pals::parula(100),limits = c(0,1),breaks=c(0,.5,1),labels=c("0",".5","1"))+
    scale_fill_gradientn("# of observations", colours = cmocean::cmocean('haline')(100)) + 
                           #breaks = c(0, 200, 400, 600, 800), limits = c(1, 811),#labels = format(bzz), 
                           #guide = guide_colorbar(barwidth = 1.5, barheight = 10)) +
    coord_fixed(xlim=xl, ylim=yl, ratio=1.3) + xlab('') + ylab('') +
    #scale_x_continuous(expand=c(0,0),limits = xl) +
    #scale_y_continuous(expand=c(0,0),limits = yl)+
    geom_map(data=testt, map=testt, aes(map_id=region), fill="darkgrey", color="darkgrey")+
    #geom_sf(data = eez, color = "black",fill=NA,size=1)+
    theme_bw(base_size = 10) + theme(panel.grid=element_blank(), panel.border = element_blank()) #+ theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"))
    #theme(#legend.position = "none",
          #legend.text = element_text(size=10),
          #legend.title = element_text(size=12),
    #      plot.margin = unit(c(0,0,3,0), "mm")) + theme(axis.ticks=element_blank(), 
    #                                                    panel.background=element_rect(fill = "white"), 
    #                                                    axis.text.x=element_blank(), axis.text.y=element_blank(),  
    #                                                    panel.grid.major = element_blank(),
    #                                                    panel.grid.minor = element_blank(),
    #                                                    axis.title.x=element_blank(), axis.title.y=element_blank())
  
}

po.1 <- plot_observer(mask(obs_pres_ras, obs_pres_mask))
po.1
po.0 <- plot_observer(mask(obs_abs_ras, obs_abs_mask))
po.0

lay <- rbind(c(1, 2))
g <- gridExtra::arrangeGrob(grobs = list(po.0, po.1), heights = c(5),
                            width = c(6, 6), layout_matrix = lay)
ggsave(file = './figures/observer_data_fig.pdf', width=12, height=5, units = 'in', g)

#p4 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill='grey') +
#  coord_fixed(xlim=xl, ylim=yl, ratio=1.3) + 
#  xlab('') + ylab('') +# ggtitle(paste('ICCAT', ctag.new$SpeciesCode[1], 'conventional tags N =', nrow(ctag.new))) +
#  geom_bin2d(data = obs, aes(x=lon, y=lat), binwidth=1) + 
#  facet_wrap(.~pres) + 
#  scale_fill_gradientn(" Frequency of\n occurence\n", colours = cmocean::cmocean('haline')(100), 
#                       breaks = c(5, 100, 200, 300, 400, 500, 600), #labels = format(bzz), 
#                       guide = guide_colorbar(barwidth = 1.5, barheight = 10)) +
#  theme_bw(base_size = 10) + theme(panel.grid=element_blank(), panel.border = element_blank()) #+ theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"))
#p4
#ggsave(file = './figures/biocomp_fig2_observer.pdf', width=8, height=4, units = 'in', p4)


# FIGURE X - Sensitivity / Specificity --------------------------------------------------------------

## look at sensitivity/specificity as models are predicted against other datasets
## esp interested in how the two fishery-dependent models predict to fishery-independent locations from satellite tags

sens <- ggplot(data = res_mat, aes(data, model, fill = accuracy_quantile))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours = cmocean::cmocean('balance')(100),
                       limits = c(0,1), space = "Lab",
                       name="Sensitivity") +
  geom_text(aes(data, model, label = round(accuracy_quantile, 2)), color = "black", size = 4) +
  theme_minimal()+ 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

res_mat$fishery_model <- NA
res_mat$fishery_model[which(res_mat$model %in% c('brt_1a','brt_4b'))] <- 'dependent'
res_mat$fishery_model[which(res_mat$model %in% c('brt_2d','brt_3d'))] <- 'independent'
res_mat$fishery_model[which(res_mat$model %in% c('brt_ens'))] <- 'ensemble'

p1 <- res_mat %>% filter(!is.na(fishery_model)) %>% 
  ggplot(aes(x=factor(fishery_model, levels=c('dependent','independent','ensemble')), y=sens)) + 
  geom_violin(fill='#728DC2') +
  geom_point() + 
  stat_summary(fun.y="median",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), width=1, show.legend = FALSE) +
  xlab('model type') +
  ylab('Proportion of correct predictions') +
  theme_minimal() #+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())
  
    
ggsave('uncertainty_boxplot.pdf', height=8, width=8, p1)  
res_mat <- res_mat %>% filter(!is.na(fishery_model))
data.table::fwrite(res_mat, file='./data/res_mat_odp.csv')

res_mat_true <- list()
counter <- 1
obs_abs <- data.table::fread('./data/enhance/observer/observer-enhanced_allAbs.csv')
for (i in 1:length(models)){
  if (i == length(models)){
    ## ensemble
    ens_preds <- list()
    for (tt in 1:(length(models) - 1)){
      ens_preds[[tt]] <- gbm::predict.gbm(models[[tt]], obs_abs,
                                          n.trees = 2000,
                                          type = 'response',
                                          na.rm = FALSE)
      
    }
    pred <- rowMeans(as.data.frame(ens_preds, col.names=seq(1:4)))
    
  } else{
    pred <- gbm::predict.gbm(models[[i]], obs_abs,
                             n.trees = 2000,
                             type = 'response',
                             na.rm = FALSE)
  }
  
  #evaluate
  pred_q <- pred
  pred_q[which(pred_q <= quantile(pred)[2])] <- 0
  pred_q[which(pred_q >= quantile(pred)[4])] <- 1
  pred_q[which(pred_q > quantile(pred)[2] & pred_q < quantile(pred)[4])] <- NA
  idx <- which(is.na(pred_q))
  pred_q <- pred_q[-idx]
  cm <- caret::confusionMatrix(factor(pred_q),
                               factor(obs_abs$pres[-idx]), positive = '1')
  d <- cbind(obs_abs$pres, pred)
  pres <- as.numeric(d[d[,1] == 1,2])
  abs <- as.numeric(d[d[,1] == 0,2])
  e <- dismo::evaluate(p=pres, a=abs)
  sens_q <- round(cm$byClass['Sensitivity'], 2)
  spec_q <- round(cm$byClass['Specificity'], 2)
  
  res_mat_true[[counter]] <- list(names(models)[i], 'obs_trueabs', q25=quantile(pred)[2], q75=quantile(pred)[4], sens_q, spec_q)
  counter <- counter + 1
  #rm(auc); rm(pred); rm(tss); rm(sens); rm(spec)
}

res_mat_true <- do.call(rbind.data.frame, res_mat_true)
names(res_mat_true) <- c('model','data','q25','q75','sens_quantile','spec_quantile')

spec_true <- ggplot(data = res_mat_true, aes(data, model, fill = spec_quantile))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours = cmocean::cmocean('balance')(100),
                       limits = c(0,1), space = "Lab",
                       name="Specificity") +
  geom_text(aes(data, model, label = round(spec_quantile, 2)), color = "black", size = 4) +
  theme_minimal()+ 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

spec <- ggplot(data = res_mat, aes(data, model, fill = spec_quantile))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours = cmocean::cmocean('balance')(100),
                       limits = c(0,1), space = "Lab",
                       name="Specificity") +
  geom_text(aes(data, model, label = round(spec_quantile, 2)), color = "black", size = 4) +
  theme_minimal()+ 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

ggsave(file = './figures/sensitivity_v6.pdf', width=6, height=6, units = 'in', sens)
ggsave(file = './figures/specificity_v6.pdf', width=6, height=6, units = 'in', spec_true)


# FIGURE x - EXAMPLE DAILY PREDICTIONS --------------------------------------------------------------
## no reason why I picked these dates, totally random in two v different seasons
jan <- raster::stack('/Volumes/Elements/glorys_daily/2019/cmems_mod_glo_phy_my_0.083_P1D-m_2019-01-01.grd')
jul <- raster::stack('/Volumes/Elements/glorys_daily/2019/cmems_mod_glo_phy_my_0.083_P1D-m_2019-07-01.grd')
models <- list(brt_1a, brt_2d, brt_3d, brt_4b, brt_5a)#, brt_ens)
names(models) <- c('iccat','sat','psat','pop','all')#,'ens')
preds <- list()

for (i in 1:length(models)){
#  pred_jan <- predict(jan, models[[i]], type="response", n.trees = 2000, na.rm=F)
  pred_jul <- predict(jul, models[[i]], type="response", n.trees = 2000, na.rm=F)
  #preds[[i]] <- list(pred_jan, pred_jul)
  preds[[i]] <- list(pred_jul)
}

# function
predict_models <- function(modtype, patype, rasStack, pred, date, out_type = 'png', data=NULL, extrap=NA, mask=FALSE, xl=NA, yl=NA){
  require(glue)
  
  df_map = rasterToPoints(pred) %>% as.data.frame()
  colnames(df_map)= c("rows","cols","value")    
  
  if (is.na(xl)) xl <- c(raster::extent(pred)@xmin, raster::extent(pred)@xmax)
  if (is.na(yl)) yl <- c(raster::extent(pred)@ymin, raster::extent(pred)@ymax)
  
  map.world = map_data(map="world")
  testt=map.world %>% filter(long<=180)
  
  a=df_map %>% 
    ggplot() + 
    # geom_tile(aes(x = rows, y = cols, fill = ntile(value,100))) + coord_equal()+
    geom_tile(aes(x = rows, y = cols, fill = value),color=NA) + coord_equal()+
    #scale_fill_gradientn("Habitat suitability",colours = pals::parula(100),limits = c(0,1),breaks=c(0,.5,1),labels=c("0",".5","1"))+
    scale_fill_gradientn("Habitat suitability", colours = cmocean::cmocean('haline')(100),limits = c(0,1), breaks=c(0,.5,1), labels=c("0",".5","1")) +
    scale_x_continuous(expand=c(0,0),limits = xl) +
    scale_y_continuous(expand=c(0,0),limits = yl)
  
  if (class(extrap) != 'logical' & !mask){
    ## Extract polygons
    #pp <- rasterToPolygons(mask, dissolve=TRUE)
     
    ## Convert SpatialPolygons to a format usable by ggplot2
    #outline <- fortify(pp)
    #a = a + geom_tile(data = outline, aes(x = long, y = lat, group = group),  
    #                  size=0.25, fill=rgb(1,1,1, alpha=0.25))
  } else if(class(extrap) != 'logical' & mask){
    
    extrap <- raster::resample(extrap, pred)
    extrap[!is.na(extrap)] <- 1
    #pred <- raster::mask(pred, extrap, maskvalue=1)
    df_map = rasterToPoints(pred) %>% as.data.frame()
    colnames(df_map)= c("rows","cols","value")    
    
    map.world = map_data(map="world")
    testt=map.world %>% filter(long<=180)
    
    a=df_map %>% 
      ggplot() + 
      # geom_tile(aes(x = rows, y = cols, fill = ntile(value,100))) + coord_equal()+
      geom_tile(aes(x = rows, y = cols, fill = value),color=NA) + coord_equal()+
      #scale_fill_gradientn("Habitat suitability",colours = pals::parula(100),limits = c(0,1),breaks=c(0,.5,1),labels=c("0",".5","1"))+
      scale_fill_gradientn("Habitat suitability", colours = cmocean::cmocean('haline')(100),limits = c(0,1), breaks=c(0,.5,1), labels=c("0",".5","1")) +
      scale_x_continuous(expand=c(0,0),limits = xl) +
      scale_y_continuous(expand=c(0,0),limits = yl)
    
    pp <- rasterToPoints(extrap) %>% as.data.frame()
    colnames(pp)= c("x","y","mic")    
    a = a + geom_tile(data=pp, aes(x = x, y = y), fill = rgb(1,1,1, alpha=0.55), color=NA)
    
  }
  
  if (!is.null(data)) a = a + geom_bin2d(data = data[which(lubridate::month(data$date) == 7 & data$pres == 1),], aes(x=lon, y=lat), binwidth=1, colour = 'black', fill = rgb(0,0,0, alpha=0), size=0.5) 
  
    a=a + geom_map(data=testt, map=testt, aes(map_id=region), fill="darkgrey", color="darkgrey")+
    #geom_sf(data = eez, color = "black",fill=NA,size=1)+
    theme(legend.position = "none",
          #legend.text = element_text(size=10),
          #legend.title = element_text(size=12),
          plot.margin = unit(c(0,0,3,0), "mm")) + theme(axis.ticks=element_blank(), 
                                                        panel.background=element_rect(fill = "white"), 
                                                        axis.text.x=element_blank(), axis.text.y=element_blank(),  
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        axis.title.x=element_blank(), axis.title.y=element_blank()) #+
  #ggtitle(glue("{modtype} {patype} {date}"))
  
  
  #if (out_type == 'pdf') pdf(file=glue("./figures/predictions/{modtype}_{patype}_{date}.pdf"), height = 5, width = 6)
  if (out_type == 'png') png(file=glue("./figures/predictions/{modtype}_{patype}_{date}.png"), height = 5, width = 6, units='in', res=600)
  par(mar=c(0,0,0,0))
  #cowplot::plot_grid(backplot,buffplot,crwplot,revplot, ncol=2, align="h",axis = "bt")
  gridExtra::grid.arrange(a, ncol=1)
  dev.off()
  #return(a)
}


for (i in 1:length(models)){
  require(dsmextra)
  # https://densitymodelling.github.io/dsmextra/articles/dsmextra.html
  
  covars <- models[[i]]$gbm.call$gbm.x
  
  if (names(models)[i] != 'all'){
    segs <- all_df %>% dplyr::filter(pres == 1 & type == names(models)[i]) %>% dplyr::select(lon, lat, all_of(covars))
    datasets_plot <- all_df %>% dplyr::filter(pres == 1 & type == names(models)[i])
  } else{
    segs <- all_df %>% dplyr::filter(pres == 1) %>% dplyr::select(lon, lat, all_of(covars))
    datasets_plot <- all_df %>% dplyr::filter(pres == 1)
  }
  names(segs)[1:2] <- c('x','y')
  
  # Target system - JULY
  target <- raster::mask(jul, mask = atl) %>%
    raster::as.data.frame(., xy = TRUE, na.rm = TRUE)

  preds_plot <- raster::mask(preds[[i]][[1]], atl)
  
  predict_models(modtype = "BRT", 
                 patype = names(models)[i], 
                 rasStack = jul, 
                 pred=preds_plot, 
                 date="2019.07.01_grid_v2", 
                 out_type = 'png', 
                 data=datasets_plot)
  
}

## now that we have the spatial predictions, we can add an ensemble "model"
#pred_ens <- mean(stack(predict(jul, brt_1a, type="response", n.trees = 2000, na.rm=F), ## using ICCAT equal N model
#                  stack(lapply(preds[2:4], FUN=function(x) stack(x)))))
pred_ens <- mean(stack(lapply(preds[1:4], FUN=function(x) stack(x))))
pred_ens
datasets_plot <- datasets[c(names(models)[1:4], 'all_df')]
for (i in 1:length(models)){
  #predict_models(modtype = "BRT", patype = names(models)[i], rasStack = jan, pred=preds[[i]][[1]], date="2019.01.01", out_type = 'png')
  predict_models(modtype = "BRT", patype = names(models)[i], rasStack = jul, pred=mask(preds[[i]][[1]], atl), date="2019.07.01_grid_v2", out_type = 'png', data=datasets_plot[[i]])
}
predict_models(modtype = "BRT", patype = 'ens', rasStack = jul, pred=mask(pred_ens, atl), date="2019.07.01_v2", out_type = 'png')

## can also try plotting these like this following ELH - https://github.com/elhazen/PA-paper/blob/main/6%20-%20FinalFigs.R


