# getting landcover for fired events

# Setup and data prep ==================
system("aws s3 sync s3://earthlab-natem/modis-burned-area/input/landcover/MCD12Q1_mosaics data/MCD12Q1_mosaics")

library(sf);library(tidyverse); library(raster)
install.packages("exactextractr");library(exactextractr)
install.packages("fasterize");library(fasterize)

download.file("https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/VEG/GLDASp4_domveg_025d.nc4",
              "data/gldas.nc4")
gldas <- raster("data/gldas.nc4")

# koppen climate classes:
# http://koeppen-geiger.vu-wien.ac.at/shifts.htm
download.file("http://koeppen-geiger.vu-wien.ac.at/data/1976-2000_GIS.zip",
              "data/koppen.zip")
unzip(zipfile = "data/koppen.zip", exdir = "data/koppen")

lut_kop <- c("Af",  "Am",  "As" , "Aw",  "BWk", "BWh" ,"BSk" ,"BSh" ,"Cfa" ,
             "Cfb", "Cfc" ,"Csa", "Csb" ,"Csc" ,"Cwa", "Cwb", "Cwc", "Dfa",
             "Dfb","Dfc", "Dfd" ,"Dsa" ,"Dsb" ,"Dsc" ,"Dsd" ,"Dwa" ,"Dwb" ,
             "Dwc" ,"Dwd" ,"EF",  "ET" )

names(lut_kop)<-c(11:14,  21,  22,  26,  27,  31:39,  41:52,  61,62)

lut_lc<-c( "Evergreen Needleleaf Forests",
 "Evergreen Broadleaf Forests",
 "Deciduous Needleleaf Forests",
 "Deciduous Broadleaf Forests",
 "Mixed Forests",
 "Closed Shrublands",
 "Open Shrublands",
 "Woody Savannas",
 "Savannas",
 "Grasslands",
 "Permanent Wetlands",
 "Croplands",
 "Urban and Built-up Lands",
 "Cropland/Natural  Vegetation  Mosaics",
 "Permanent Snow and Ice",
 "Barren",
 "Water Bodies")
names(lut_lc) <- 1:17
# lut_kop <- c("Af"=11,
#               "Am"=12,
#               "As"=13,
#               "Aw" = 14,
#               "BWk" =21,
#               "BWh" = 22,
#               "BSk" = 26,
#               "BSh" = 27,
#               "Cfa" = 31,
#              "Cfb" = 32 ,
#               "Cfc" = 33,
#               "Csa" = 34,
#               "Csb" = 35,
#               "Csc" = 36,
#               "Cwa" = 37,
#               "Cwb" = 38,
#               "Cwc" = 39,
#               "Dfa" = 41,
#               "Dfb" = 42,
#               "Dfc" = 43,
#               "Dfd" = 44,
#               "Dsa" = 45,
#               "Dsb" = 46,
#               "Dsc" = 47,
#               "Dsd" = 48,
#               "Dwa" = 49,
#               "Dwb" = 50,
#               "Dwc" = 51,
#               "Dwd" = 52,
#               "EF" = 61,
#               "ET" = 62)
kop <- st_read("data/koppen")
kop_r<- fasterize(sf=kop, raster = gldas, field="GRIDCODE", fun="max")

NAvalue(gldas) <- 0
NAvalue(kop_r) <- 0
years <- 2016:2019


# The business ============
na <- st_read("data/fired_na_2017-nids.gpkg") %>%
  mutate(lc_year = as.numeric(str_sub(first_date_7,1,4))-1)

res<-list()
for(y in 1:length(years)){
  lcfiles <- list.files("data/MCD12Q1_mosaics/", full.names = TRUE)
  r <- raster(lcfiles[[y]])
  NAvalue(r) <- 0
  res[[y]] <- na %>%
    filter(lc_year == years[y]) %>%
    mutate(lc = exact_extract(x=r,y= ., 'mode'))
}

na_lc<- do.call('rbind',res) %>%
  dplyr::select(-lc_year) %>%
  mutate(gldas = exact_extract(x=gldas, y=., 'mode'),
         koppen = exact_extract(x=kop_r, y=., 'mode'))
st_write(na_lc, "data/fired_na_2017-nids_lc.gpkg", delete_dsn = TRUE)
system("aws s3 cp data/fired_na_2017-nids_lc.gpkg s3://earthlab-amahood/night_fires/fired_na_2017-nids_lc.gpkg")

# sa 
sa <- st_read("data/fired_sa_2017-nids.gpkg") %>%
  mutate(lc_year = as.numeric(str_sub(first_date_7,1,4))-1)

res_s<-list()
for(y in 1:length(years)){
  lcfiles <- list.files("data/MCD12Q1_mosaics/", full.names = TRUE)
  r <- raster(lcfiles[[y]])
  NAvalue(r) <- 0
  res_s[[y]] <- sa %>%
    filter(lc_year == years[y]) %>%
    mutate(lc = exact_extract(x=r,y= ., 'mode'))
}

sa_lc<- do.call('rbind',res_s) %>%
  dplyr::select(-lc_year)%>%
  mutate(gldas = exact_extract(x=gldas, y=., 'mode'),
         koppen = exact_extract(x=kop_r, y=., 'mode'))
st_write(sa_lc, "data/fired_sa_2017-nids_lc.gpkg", delete_dsn = TRUE)
system("aws s3 cp data/fired_sa_2017-nids_lc.gpkg s3://earthlab-amahood/night_fires/fired_sa_2017-nids_lc.gpkg")

# trying to figure out gldas ===== moral of the story - don't bother ===========
na_tundra <- na_lc %>%
  filter(gldas > 17)

ggplot(na_tundra, aes(x=as.factor(gldas), fill = as.factor(lc))) +
  geom_bar(stat="count", position="dodge", color="black") +
  scale_fill_brewer(palette = "Paired") +
  geom_text(x=as.factor(19), y=750, label=nrow(na_lc)%>%formatC(big.mark = ","))

# splitting up events by landcover and latitude ================================
lut_clim<- c("Equatorial", "Arid", "Temperate", "Boreal", "Polar")
names(lut_clim)<-c("A", "B", "C", "D", "E")
test<-na_lc %>%
  rbind(sa_lc)%>%
  dplyr::select(-gldas) %>%
  na.omit() %>%
  mutate(kop_c=lut_kop[as.character(koppen)],
         lc_c = lut_lc[as.character(lc)],
         main_clim=lut_clim[str_sub(kop_c,1,1)])

ggplot(test, aes(x=main_clim, fill=lc_c)) +
  geom_bar(stat="count", position="dodge") +
  scale_y_continuous(labels=scales::label_comma())

dir.create("data/lc_splits")
for (i in test$main_clim){
  for( j in test$lc){
    test %>%
      filter(main_clim == i & lc == j) %>%
      st_write(paste0("data/lc_splits/clim_",i,"_lc_",j,".gpkg"))
  }
}

