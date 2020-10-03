#make yearly landcover mosaics
library(tidyverse)
library(gdalUtils)
library(raster)
dir.create("data/MCD12Q1_tifs/")
# getting the tifs


lc_files <- list.files(file.path("data/hdf"), full.names = TRUE)
for(f in lc_files){
  tile <- str_extract(f, "h\\d{2}v\\d{2}")
  out_file <- paste0("data/MCD12Q1_tifs/",tile,".tif")
  sds <- f %>%
    gdalUtils::get_subdatasets()
  gdalUtils::gdal_translate(sds[1], dst_dataset = out_file)
 print(paste(tile))
}
