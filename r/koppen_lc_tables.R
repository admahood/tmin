library(tidyverse)
library(stars)
library(s2)

# koppen look up tables
lut_kop<- c("Equatorial", "Arid", "Temperate", "Boreal", "Polar")
names(lut_kop) <- c(1,2,3,4,5)

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
           "Cropland Natural Vegetation Mosaics",
           "Permanent Snow and Ice",
           "Barren",
           "Water Bodies")
names(lut_lc) <- str_pad(1:17, width = 2, side = "left",pad = "0")

lck<-read_stars("data/lc_koppen_2010_mode.tif")
lck_poly<-lck %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>%
  as_tibble() %>% 
  select(lck = lc_koppen_2010_mode.tif, geometry)  %>% 
  mutate(geometry = as_s2_geography(geometry)) %>%
  mutate(area_km2 = s2_area(geometry)/1000000)

lck_tab <- lck_poly %>%
  st_as_sf %>%
  st_set_geometry(NULL) %>%
  group_by(lck) %>%
  summarise(area_km2 = sum(area_km2)) %>%
  ungroup() %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         lc_kop = paste(kop, lc)) %>%
  na.omit

thresholds <- read_csv("data/updated-goes-af-vpd-thresholds.csv") %>%
  dplyr::select(lc_kop = lc_name, vpd_thresh_hpa, sd) %>%
  left_join(lck_tab %>% dplyr::select(lc_kop, area_km2)) %>%
  mutate(millions_of_km2 = area_km2/1000000) %>%
  # dplyr::select(-area_km2) %>%
  arrange(lc_kop)

write_csv(thresholds %>% dplyr::select(-area_km2),"data/lck_thresh_area.csv")

#148940000 km2 (from wikipedia)
total_land_area = lck_tab %>% pull(area_km2) %>% sum
burnable_land_area = thresholds %>% pull(area_km2) %>% sum()
burnable_land_area/(total_land_area)

xxx<-st_crs(read_stars("data/adjusted_counts/op-adjusted_D_200301.tif"))
yyy<-st_crs(lck)
xxx==yyy

# getting extended table 1 stats right here baby (michael's script has a lot of 
# adjusting that is no longer necessary)
