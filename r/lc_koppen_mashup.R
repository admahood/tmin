# script to create landcover layer for use in paper. For the revised analysis,
# we're using both the Koppen climate classificaitons and the MODIS MCDQ1
# landcover to generate vpd thresholds for the burnable landscape. This script
# is going to combine those two layers for visualization and use with climate 
# data

# setup =============
libs <- c("tidyverse", "raster", "ncdf4", "exactextractr", "stars")
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

# if(!exists("kop")){
lc <- read_stars(file.path(lc_path_local, fn)) %>%
  st_warp(dest = template)

kop <- raster('data/koppen/Beck_KG_V1_present_0p0083.tif') 
# reclassifying the kop the old-fashioned way
kop[kop<4] <- 1
kop[kop>3 & kop<8] <- 2
kop[kop>7 & kop<17] <- 3
kop[kop>16 & kop < 29] <- 4
kop[kop > 28] <- 5
# take a while so we'll save it
# save(lc,kop, file="data/lc_kop.Rda")
# }else{load("data/lc_kop.Rda")}


# # SANITY CHECK
# freq(kop)
# freq(lc)
