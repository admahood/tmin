# getting landcover for fired events

system("aws s3 sync s3://earthlab-natem/modis-burned-area/input/landcover/MCD12Q1_mosaics data/MCD12Q1_mosaics")

library(sf);library(tidyverse); library(raster)
install.packages("exactextractr");library(exactextractr)

download.file("https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/VEG/GLDASp4_domveg_025d.nc4",
              "data/gldas.nc4")
gldas <- raster("data/gldas.nc4")
years <- 2016:2019

na <- st_read("data/fired_na_2017-nids.gpkg") %>%
  mutate(lc_year = as.numeric(str_sub(first_date_7,1,4))-1)

res<-list()
for(y in 1:length(years)){
  lcfiles <- list.files("data/MCD12Q1_mosaics/", full.names = TRUE)
  r <- raster(lcfiles[[y]])
  r[r==0] <- NA
  res[[y]] <- na %>%
    filter(lc_year == years[y]) %>%
    mutate(lc = exact_extract(x=r,y= ., 'mode'))
}

na_lc<- do.call('rbind',res) %>%
  dplyr::select(-lc_year) %>%
  mutate(gldas = exact_extract(x=gldas, y=., 'mode'))
st_write(na_lc, "data/fired_na_2017-nids_lc.gpkg", delete_dsn = TRUE)
system("aws s3 cp data/fired_na_2017-nids_lc.gpkg s3://earthlab-amahood/night_fires/fired_na_2017-nids_lc.gpkg")

# sa 
sa <- st_read("data/fired_sa_2017-nids.gpkg") %>%
  mutate(lc_year = as.numeric(str_sub(first_date_7,1,4))-1)

res<-list()
for(y in 1:length(years)){
  lcfiles <- list.files("data/MCD12Q1_mosaics/", full.names = TRUE)
  r <- raster(lcfiles[[y]])
  r[r == 0] <- NA
  res[[y]] <- sa %>%
    filter(lc_year == years[y]) %>%
    mutate(lc = exact_extract(x=r,y= ., 'mode'))
}

sa_lc<- do.call('rbind',res) %>%
  dplyr::select(-lc_year)%>%
  mutate(gldas = exact_extract(x=gldas, y=., 'mode'))
st_write(sa_lc, "data/fired_sa_2017-nids_lc.gpkg", delete_dsn = TRUE)
system("aws s3 cp data/fired_sa_2017-nids_lc.gpkg s3://earthlab-amahood/night_fires/fired_sa_2017-nids_lc.gpkg")
