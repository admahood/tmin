# streamline ids for NA and SA
library(tidyverse) ; library(sf); library(vroom)

# sa first

sa<- st_read("data/fired_sa_2017-.gpkg") %>%
  dplyr::select(-is_sa) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(lat = st_coordinates(st_transform(st_centroid(.), 4326))[,2])


na <- st_read("data/fired_na_2017-.gpkg") %>%
  dplyr::select(-is_na)

lut_naids <- 1:nrow(na) + nrow(sa)
names(lut_naids) <- na$id

na_nid <- na %>%
  mutate(id=as.character(id),
         nid = lut_naids[id]) %>%
  mutate(lat = st_coordinates(st_transform(st_centroid(.), 4326))[,2])

na_vpd <- vroom("data/na_vpd_long_2017-.csv") %>%
  mutate(fireID=as.character(fireID),
         nid = lut_naids[fireID])

system("aws s3 cp data/fired_na_2017-.gpkg s3://earthlab-amahood/night_fires/fired_na_2017-.gpkg")
system("aws s3 cp data/fired_sa_2017-.gpkg s3://earthlab-amahood/night_fires/fired_sa_2017-.gpkg")
