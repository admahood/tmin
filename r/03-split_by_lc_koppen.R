# getting landcover and koppen-gieger climate classifications for fired events

# Setup and data prep ==================
system("aws s3 sync s3://earthlab-natem/modis-burned-area/input/landcover/MCD12Q1_mosaics data/MCD12Q1_mosaics")

library(sf);library(tidyverse); library(raster)
install.packages("exactextractr");library(exactextractr)
install.packages("fasterize");library(fasterize)

download.file("https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/VEG/GLDASp4_domveg_025d.nc4",
              "data/gldas.nc4")
gldas <- raster("data/gldas.nc4")

download.file("https://ndownloader.figshare.com/files/12407516",
              "data/koppen.zip")
unzip(zipfile = "data/koppen.zip", exdir = "data/koppen")

# koppen from https://www.nature.com/articles/sdata2018214#Sec8
lut_kop <- c("Af",  "Am", "Aw",   
             "BWh" ,"BWk","BSh" ,"BSk" ,
             "Csa", "Csb" ,"Csc" ,"Cwa", "Cwb", "Cwc","Cfa" , "Cfb", "Cfc" ,
             "Dsa" ,"Dsb" ,"Dsc" ,"Dsd","Dwa" ,"Dwb" ,"Dwc" ,"Dwd" ,
             "Dfa", "Dfb","Dfc", "Dfd" , 
             "ET",   "EF")
names(lut_kop)<-c(1:30)

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
 "Cropland/Natural Vegetation Mosaics",
 "Permanent Snow and Ice",
 "Barren",
 "Water Bodies")
names(lut_lc) <- 1:17

kop_r <- raster('data/koppen/Beck_KG_V1_present_0p0083.tif')

NAvalue(gldas) <- 0
NAvalue(kop_r) <- -Inf
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

# splitting up events by landcover and koppen class ============================
lut_clim<- c("Equatorial", "Arid", "Temperate", "Boreal", "Polar")
names(lut_clim)<-c("A", "B", "C", "D", "E")
test<-na_lc %>%
  rbind(sa_lc)%>%
  dplyr::select(-gldas) %>%
  na.omit() %>%
  mutate(kop_c=lut_kop[as.character(koppen)],
         lc_c = lut_lc[as.character(lc)],
         main_clim=lut_clim[str_sub(kop_c,1,1)],
         lc_c = str_replace_all(lc_c," ","_"),
         lc_c = str_replace_all(lc_c,"/","_"))

ggplot(test, aes(x=main_clim, fill=lc_c)) +
  geom_bar(stat="count", position="dodge") +
  scale_y_continuous(labels=scales::label_comma())

dir.create("data/tallys")
test%>%
  st_set_geometry(NULL) %>%
  group_by(lc_c) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  write_csv("data/tallys/event_tallys_lc.csv")

test%>%
  st_set_geometry(NULL) %>%
  group_by(main_clim,lc_c) %>%
  summarise(n=n())%>%
  write_csv("data/tallys/event_tallys_kop_lc.csv")

test%>%
  st_set_geometry(NULL) %>%
  group_by(kop_c,lc_c) %>%
  summarise(n=n())%>%
  write_csv("data/tallys/event_tallys_kopfull_lc.csv")

dir.create("data/lc_splits")
for (i in unique(test$main_clim)){
  for( j in unique(test$lc_c)){
    outfile<-paste0(i,"_",j,".gpkg")
    out<- test %>%
      filter(main_clim == i & lc_c == j) 
    if(nrow(out)>0) st_write(out, paste0("data/", "lc_splits/", outfile), delete_dsn=TRUE)
  }
}
system("aws s3 sync data/lc_splits s3://earthlab-amahood/night_fires/lc_splits")
system("aws s3 sync data/tallys s3://earthlab-amahood/night_fires/tallys")
