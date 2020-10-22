# streamline ids for NA and SA
library(tidyverse) ; library(sf); library(vroom)

system("aws s3 sync s3://earthlab-amahood/night_fires data")

# sa first

sa<- st_read("data/fired_sa_2017-.gpkg") %>%
  dplyr::select(-is_sa)

sa_rows <- nrow(sa)

lut_saids <- 1:sa_rows
names(lut_saids) <- sa$id

sa %>%
  mutate(id=as.character(id),
         nid = lut_saids[id] %>% as.numeric(),
         id = as.numeric(id)) %>%
  mutate(lat = st_coordinates(st_transform(st_centroid(.), 4326))[,2])%>%
  st_write("data/fired_sa_2017-nids.gpkg", delete_dsn=TRUE)

system("aws s3 cp data/fired_sa_2017-nids.gpkg s3://earthlab-amahood/night_fires/fired_sa_2017-nids.gpkg")

rm(sa)
# north america

na <- st_read("data/fired_na_2017-.gpkg") %>%
  dplyr::select(-is_na)

lut_naids <- 1:nrow(na) + sa_rows
names(lut_naids) <- na$id

save(lut_naids, lut_saids, file = "data/luts.Rda")
system("aws s3 cp data/luts.Rda s3://earthlab-amahood/night_fires/luts.Rda")

na %>%
  mutate(id=as.character(id),
         nid = lut_naids[id] %>% as.numeric(),
         id = as.numeric(id)) %>%
  mutate(lat = st_coordinates(st_transform(st_centroid(.), 4326))[,2]) %>%
  st_write("data/fired_na_2017-nids.gpkg")
system("aws s3 cp data/fired_na_2017-nids.gpkg s3://earthlab-amahood/night_fires/fired_na_2017-nids.gpkg")

rm(na)
# now vpd monsters
sa_vpd <- vroom("data/sa_vpd_long_2017-test.csv") %>%
  mutate(fireID=as.character(fireID),
         nid = lut_saids[fireID],
         fireID=as.numeric(fireID),
         nid=as.numeric(nid))
write_csv(sa_vpd, "data/sa_vpd_long_2017-nids.csv")
system("aws s3 cp data/sa_vpd_long_2017-nids.csv s3://earthlab-amahood/night_fires/sa_vpd_long_2017-nids.csv")
rm(sa_vpd)
gc()

na_vpd <- vroom("data/na_vpd_long_2017-.csv") %>%
  mutate(fireID=as.character(fireID),
         nid = lut_saids[fireID],
         fireID=as.numeric(fireID),
         nid=as.numeric(nid))
write_csv(na_vpd, "data/na_vpd_long_2017-nids.csv")

system("aws s3 cp data/na_vpd_long_2017-nids.csv s3://earthlab-amahood/night_fires/na_vpd_long_2017-nids.csv")
rm(na_vpd)


