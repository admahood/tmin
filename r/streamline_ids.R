# streamline ids for NA and SA
library(tidyverse) ; library(sf); library(vroom); library(lubridate); library(data.table)

system("aws s3 sync s3://earthlab-amahood/night_fires/fired_polys data")

# effort (number of goes scenes per hour)
system("aws s3 cp s3://earthlab-mkoontz/goes16meta/sampling-effort-goes16.csv data/sampling-effort-goes16.csv")

# polygons =====================================================================
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

system("aws s3 cp data/fired_sa_2017-nids.gpkg s3://earthlab-amahood/night_fires/fired_polys/fired_sa_2017-nids.gpkg")

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
  st_write("data/fired_na_2017-nids.gpkg", delete_dsn=T)
system("aws s3 cp data/fired_na_2017-nids.gpkg s3://earthlab-amahood/night_fires/fired_polys/fired_na_2017-nids.gpkg")

rm(na)
# now vpd monsters =============================================================
sa_vpd <- vroom("data/sa_vpd_long_2017-test.csv") %>%
  mutate(fireID=as.character(fireID),
         nid = lut_saids[fireID],
         fireID=as.numeric(fireID),
         nid=as.numeric(nid))
write_csv(sa_vpd, "data/sa_vpd_long_2017-nids.csv")
system("aws s3 cp data/sa_vpd_long_2017-nids.csv s3://earthlab-amahood/night_fires/sa_vpd_long_2017-nids.csv")
rm(sa_vpd)
gc()

system("aws s3 cp s3://earthlab-amahood/night_fires/na_vpd_long_2017-.csv data/na_vpd_long_2017-.csv")
system("aws s3 cp s3://earthlab-amahood/night_fires/luts.Rda data/luts.Rda")
load("data/luts.Rda")

na_vpd <- vroom("data/na_vpd_long_2017-.csv") %>%
  mutate(fireID=as.character(fireID),
         nid = lut_naids[fireID],
         fireID=as.numeric(fireID),
         nid=as.numeric(nid))
vroom_write(na_vpd, "data/na_vpd_long_2017-nids.csv")
system("aws s3 cp data/na_vpd_long_2017-nids.csv s3://earthlab-amahood/night_fires/na_vpd_long_2017-nids.csv")
rm(na_vpd)

# joining vpd csvs, then splitting by landcover and koppen =====================

system("aws s3 cp s3://earthlab-amahood/night_fires/na_vpd_long_2017-nids.csv data/na.csv")
system("aws s3 cp s3://earthlab-amahood/night_fires/sa_vpd_long_2017-nids.csv data/sa.csv")
system("aws s3 sync s3://earthlab-amahood/night_fires/lc_splits data/lc_splits")

files <- list.files("data", full.names = TRUE, pattern = ".csv")[c(2,4)]
wh <- vroom(files) 

# wh <- wh %>% 
#   separate(hour, into = c("day", "hour"),sep="_",) # maybe do this in smaller chunks

vroom_write(wh, "data/wh_vpd.csv")
system("aws s3 cp data/wh_vpd.csv s3://earthlab-amahood/night_fires/wh_vpd.csv")


dir.create("data/vpd_lc")
fired_files <- list.files("data/lc_splits", pattern = ".gpkg", full.names = TRUE)
out_files<- list.files("data/vpd_lc", pattern = "csv")
wh<- vroom("data/wh_vpd.csv")

for(f in fired_files){
  out_fn <- str_replace(f, "lc_splits", "") %>%
    str_replace(".gpkg", "_vpds.csv") %>%
    str_replace_all("/", "") %>%
    str_replace("data", "")
  if(!file.exists(file.path("data","vpd_lc",out_fn))){
    firez<- st_read(f)
    
    lut_dates <- firez$first_date_7
    names(lut_dates) <- firez$nid
    
    ids <- firez %>%
      pull(nid)
    
    
    
    
    subsettt <- filter(wh, nid %in% ids) %>%
      separate(hour, into = c("day", "hour"), sep="_", convert=TRUE) %>%
      mutate(first_date = as.Date(lut_dates[as.character(nid)]),
             date = first_date + day)
    
    vroom_write(subsettt, file.path("data", "vpd_lc", out_fn))
    system(paste("aws s3 cp",
                 file.path("data", "vpd_lc", out_fn),
                 file.path("s3://earthlab-amahood", "night_fires","vpd_lc",out_fn)))
    
    rm(subsettt); rm(firez);gc()
  }
}
