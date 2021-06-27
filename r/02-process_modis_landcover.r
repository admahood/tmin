#make yearly landcover mosaics
library(tidyverse)
library(gdalUtils)
library(raster)
dir.create("data/MCD12Q1_tifs/")
# getting the tifs

years <- 2016:2019

for(y in years){
  system(paste0("aws s3 sync s3://earthlab-natem/modis-burned-area/input/landcover/",y, " data/MCD12Q1/", y))
}

for(y in years){
lc_files <- list.files(file.path("data/MCD12Q1/",y), full.names = TRUE)
dir.create(paste0("data/MCD12Q1_tifs/",y))

for(f in lc_files){
  tile <- str_extract(f, "h\\d{2}v\\d{2}")
  out_file <- paste0("data/MCD12Q1_tifs/",y,"/",tile,".tif")
  sds <- f %>%
    gdalUtils::get_subdatasets()
  gdalUtils::gdal_translate(sds[1], dst_dataset = out_file)
 print(paste(tile))
}}


dir.create("data/MCD12Q1_mosaics")
library(raster)
for(y in years){
  ff <- list.files(file.path("data/MCD12Q1_tifs", y), full.names = TRUE)
  rr <- list()
  for(fff in ff){
    rr[[fff]] <- raster(fff)
  }
  names(rr) <- NULL
  rr$fun = mean
  rr$na.rm = TRUE
  xx <- do.call(mosaic,rr)
  writeRaster(xx,paste0("data/MCD12Q1_mosaics/","lc_mosaic_",y,".tif"))
}

system(paste0("aws s3 sync data/MCD12Q1_mosaics/ s3://earthlab-natem/modis-burned-area/input/landcover/MCD12Q1_mosaics"))
