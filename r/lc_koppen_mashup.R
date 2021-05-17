# script to create landcover layer for use in paper. For the revised analysis,
# we're using both the Koppen climate classifications and the MODIS MCDQ1
# landcover to generate vpd thresholds for the burnable landscape. This script
# is going to combine those two layers for visualization and use with climate 
# data

# setup =============
libs <- c("tidyverse", "raster", "ncdf4", "stars")
sapply(libs, library, character.only=T)

# get landcover
lc_path_s3 <- "s3://earthlab-amahood/night_fires/landcover/MCD12Q1_mosaics"
lc_path_local <- "data/lc"
fn <- "lc_mosaic_2010.tif"

system(paste("aws s3 cp",
             file.path(lc_path_s3,fn ),
             file.path(lc_path_local, fn)))

# get koppen
download.file("https://ndownloader.figshare.com/files/12407516",
              "data/koppen.zip")
unzip(zipfile = "data/koppen.zip", exdir = "data/koppen")

# what the koppen values mean - we only need the first letter
lut_kop <- c("Af",  "Am", "Aw",   
             "BWh" ,"BWk","BSh" ,"BSk" ,
             "Csa", "Csb" ,"Csc" ,"Cwa", "Cwb", "Cwc","Cfa" , "Cfb", "Cfc" ,
             "Dsa" ,"Dsb" ,"Dsc" ,"Dsd","Dwa" ,"Dwb" ,"Dwc" ,"Dwd" ,
             "Dfa", "Dfb","Dfc", "Dfd" , 
             "ET",   "EF")

lut_kop <- c(rep(1,3), rep(2,4), rep(3,9), rep(4,12), rep(5,2))
names(lut_kop)<-c(1:30)

# get template
template_path_s3 <- "s3://earthlab-amahood/night_fires/test1.nc"
template_path_local <- "data/test1.nc"
system(paste("aws s3 cp",
             template_path_s3,
             template_path_local))

# data ingest =================
template <-raster("data/test1.nc") %>%
  st_as_stars()

lc_path_s3 <- "s3://earthlab-amahood/night_fires/lc.Rda"
lc_path_local <- "data/lc.Rda"
system(paste("aws s3 cp",
             lc_path_s3,
             lc_path_local))


if(!file.exists("data/lc.Rda")){
  lc <- read_stars(file.path(lc_path_local, fn)) %>%
    st_warp(dest = template, use_gdal=TRUE,
            method = "mode")
  save(lc, file="data/lc.Rda")
  system("aws s3 cp data/lc.Rda s3://earthlab-amahood/night_fires/lc.Rda")
}else{load("data/lc.Rda")}


kop_path_s3 <- "s3://earthlab-amahood/night_fires/kop.Rda"
kop_path_local <- "data/kop.Rda"
system(paste("aws s3 cp",
             kop_path_s3,
             kop_path_local))

if(!file.exists("data/kop.Rda")){
  kop <- raster('data/koppen/Beck_KG_V1_present_0p083.tif')  
  # reclassifying the kop the old-fashioned way
  kop[kop<4] <- 1
  kop[kop>3 & kop<8] <- 2
  kop[kop>7 & kop<17] <- 3
  kop[kop>16 & kop < 29] <- 4
  kop[kop > 28] <- 5
  
  kop <- kop %>%
    st_as_stars %>%
    st_warp(dest = template, use_gdal=TRUE,
            method = "mode")
  save(kop, file="data/kop.Rda")
  system("aws s3 cp data/kop.Rda s3://earthlab-amahood/night_fires/kop.Rda")
}else{load("data/kop.Rda")}


# putting them together
kop_lc2010 <- (kop*100) + lc
stars::write_stars(kop_lc2010, "lc_koppen_2010_mode.tif")


raster::raster("lc_koppen_2010_mode.tif")->xx

xx[xx==117] <- NA
xx[xx==517] <- NA
xx[xx==417] <- NA
xx[xx==317] <- NA
xx[xx==217] <- NA
plot(xx)

writeRaster(xx, "data/lc_koppen_2010_mode.tif", overwrite=TRUE)
system("aws s3 cp data/lc_koppen_2010_mode.tif s3://earthlab-amahood/night_fires/lc_koppen_2010_mode.tif")
