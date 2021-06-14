library(tidyverse)
library(stars)
library(raster)
library(s2)
library(foreach)
library(doParallel)

# koppen look up tables=========================================================
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

# reading lck in ===============================================================
lck<-read_stars("in/lc_koppen_2010_mode.tif")
lck_s <- read_stars("out/aggregations_2003-2020/lck_shifted.tif")
d_ovp <- read_stars("out/aggregations_2003-2020/2003-2020_day_overpass-count.tif")
n_ovp <- read_stars("out/aggregations_2003-2020/2003-2020_night_overpass-count.tif")

# calculating area =============================================================
lck_poly<-lck_s %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>%
  as_tibble() %>% 
  dplyr::select(lck = lck_shifted.tif, geometry)  %>% 
  mutate(geometry = as_s2_geography(geometry)) %>%
  mutate(area_km2 = s2_area(geometry)/1000000)


lck_sf<-lck_s %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>%
  mutate(area_km2 = st_area(.)/1000000)

d_ovp_area<-d_ovp %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>% 
  mutate(area_km2 =(st_area(.) %>% as.numeric)/1000000)

n_ovp_area<-n_ovp %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>% 
  mutate(area_km2 =(st_area(.) %>% as.numeric)/1000000)

st_join(lck_sf, d_ovp_area)


  
  

# summarising area by landcover/koppen =========================================
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

# making a table with the thresholds ===========================================
thresholds <- read_csv("in/updated-goes-af-vpd-thresholds.csv") %>%
  dplyr::select(lc_kop = lc_name, vpd_thresh_hpa, sd) %>%
  left_join(lck_tab %>% dplyr::select(lc_kop, area_km2)) %>%
  mutate(millions_of_km2 = area_km2/1000000) %>%
  # dplyr::select(-area_km2) %>%
  arrange(lc_kop)

write_csv(thresholds %>% dplyr::select(-area_km2),"out/lck_thresh_area.csv")

# making a table of burnable koppen area =======================================

burnable_koppen <-
  thresholds %>%
  separate(lc_kop, into = "koppen") %>%
  group_by(koppen) %>%
  summarise(area_Mkm2 = sum(millions_of_km2)) %>%
  ungroup()

write.csv(burnable_koppen, "out/burnable_koppen.csv")

# making a burnable globe mask

burnable <- lck_poly %>%
  st_as_sf %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         lc_kop = paste(kop, lc)) %>%
  filter(lc_kop %in% pull(thresholds, lc_kop))

# grabbing some quick numbers for the paper
# earth has 148940000 km2 of land surface area (from wikipedia)
total_land_area = lck_tab %>% pull(area_km2) %>% sum #145M km2 -- close enough
burnable_land_area = thresholds %>% pull(area_km2) %>% sum() # 90M km2
burnable_land_area/(total_land_area) # this is where I got the 62% for the paper

# getting extended table 1 stats right here baby (michael's script has a lot of 
# adjusting that is no longer necessary)

# doing the overpass adjustment for 0.25 degrees================================
# dir.create("data/adjusted_counts_025", recursive=T)
# dir.create("data/adjusted_frp_025", recursive=T)

# getting files ================================================================
system(paste("aws s3 sync",
             "s3://earthlab-jmcglinchy/night_fire/gridded/vars_refresh_may2021/CSV_nocorn_grid_0_25_degree_vars",
             "data/gridded_mod14_025",
             "--only-show-errors"))

system(paste("aws s3 sync",
             "s3://earthlab-mkoontz/MODIS-overpass-counts_0.25_analysis-ready/2003-2020",
             "out/aggregations_2003-2020",
             "--only-show-errors"))

# aggregation ==================================================================
dir.create("out/aggregations_2003-2020")

mod14_day_counts_025 <- list.files("data/gridded_mod14_025/AFC_num", 
                                  full.names = TRUE,
                                  pattern = "_D_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002)%>%
  pull(value) %>%
  raster::stack()
beginCluster()
raw_day_counts <- clusterR(mod14_day_counts_025, function(x)sum(x, na.rm=T), verbose=T)
endCluster()

writeRaster(raw_day_counts, "out/aggregations_2003-2020/raw_day_counts.tif",
            overwrite=TRUE)

mod14_night_counts_025 <- list.files("data/gridded_mod14_025/AFC_num", 
                                   full.names = TRUE,
                                   pattern = "_N_") %>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  pull(value) %>%
  raster::stack()

beginCluster()
raw_night_counts <- clusterR(mod14_night_counts_025, 
                             function(x)sum(x, na.rm=T), verbose=T)
endCluster()

writeRaster(raw_night_counts, "out/aggregations_2003-2020/raw_night_counts.tif",
            overwrite=TRUE)

# day and night total frp
mod14_day_frp_025<- list.files("data/gridded_mod14_025/FRP_total", 
                                  full.names = TRUE,
                                  pattern = "_D_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  pull(value) %>%
  raster::stack()

beginCluster()
raw_day_frp_sum <- clusterR(mod14_day_frp_025, 
                            function(x)sum(x, na.rm=T), verbose=T)
endCluster()

writeRaster(raw_day_frp_sum, "out/aggregations_2003-2020/raw_day_frp_sum.tif",
            overwrite=TRUE)

mod14_night_frp_025 <- list.files("data/gridded_mod14_025/FRP_total", 
                                     full.names = TRUE,
                                     pattern = "_N_") %>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  pull(value) %>%
  raster::stack()

beginCluster()
raw_night_frp_sum <- clusterR(mod14_night_frp_025, 
                              function(x)sum(x, na.rm=T), verbose=T)
endCluster()

writeRaster(raw_night_frp_sum, "out/aggregations_2003-2020/raw_night_frp_sum.tif",
            overwrite=TRUE)

# shifting the lck layer =======================================================

raster::shift(x = raster::raster("in/lc_koppen_2010_mode.tif"), dx = -179.75, dy = -0.25) %>%
  writeRaster("out/aggregations_2003-2020/lck_shifted.tif")

# making the table from the aggregated rasters =================================

lck_sf<-c(lck_s, d_ovp, n_ovp, 
          read_stars("out/aggregations_2003-2020/raw_day_counts.tif"),
          read_stars("out/aggregations_2003-2020/raw_night_counts.tif"),
                     along = "band") %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>%
  mutate(area_km2 = (st_area(.) %>% as.numeric)/1e6) %>%
  st_set_geometry(NULL)%>%
  dplyr::rename(lck = lck_shifted.tif.V1, 
                d_ovp = lck_shifted.tif.V2, 
                n_ovp = lck_shifted.tif.V3,
                dc = lck_shifted.tif.V4,
                nc = lck_shifted.tif.V5) %>%
  filter(!is.na(lck))%>%
  mutate(day_afd_per_ovp = 1e6*(dc/d_ovp)/area_km2,
         night_afd_per_ovp = 1e6*(nc/n_ovp)/area_km2,
         kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         lc_kop = paste(kop, lc)) %>%
  left_join(thresholds %>% dplyr::select(-area_km2),
            ., by = "lc_kop")%>%
  group_by(lc_kop) %>%
  summarise(land_area_Mkm2 = sum(area_km2)/1e6,
            total_day_detections_per_ovp_per_Mkm2 = mean(day_afd_per_ovp, na.rm=T),
            total_night_detections_per_ovp_per_Mkm2 = mean(night_afd_per_ovp, na.rm=T),
            pct_night_afd = 100*total_night_detections_per_ovp_per_Mkm2/(total_night_detections_per_ovp_per_Mkm2+total_day_detections_per_ovp_per_Mkm2)) %>%
  ungroup
            
            

final_table <- list.files("out/aggregations_2003-2020",
                                 full.names = TRUE) %>%
  raster::stack() %>%
  as.data.frame() %>%
  dplyr::select(lck = lck_shifted,
                raw_day_counts,
                raw_night_counts,
                raw_day_frp_sum,
                raw_night_frp_sum)%>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         lc_kop = paste(kop, lc)) %>%
  left_join(thresholds %>% dplyr::select(-area_km2),
                         ., by = "lc_kop")%>%
  group_by(lc_kop) %>%
  summarise(sum_day_counts = sum(raw_day_counts, na.rm=T), 
            sum_night_counts = sum(raw_night_counts, na.rm=T),
            sum_day_frp = sum(raw_day_frp_sum, na.rm=T),
            sum_night_frp = sum(raw_night_frp_sum, na.rm=T),
            millions_of_km2 = first(millions_of_km2),
            sum_day_ovp = sum(day_ovp),
            sum_night_ovp= sum(night_ovp)) %>%
  ungroup  %>%
  mutate(day_counts_per_ovp_perMkm2 = (sum_day_counts/sum_day_ovp)*millions_of_km2,
         night_counts_per_ovp_perMkm2 = (sum_night_counts/sum_night_ovp)*millions_of_km2,
         day_frp_per_ovp_perMkm2 = sum_day_frp/sum_day_ovp*millions_of_km2,
         night_frp_per_ovp_perMkm2 = sum_night_frp/sum_night_ovp*millions_of_km2,
         night_count_pct = night_counts_per_ovp_perMkm2/(day_counts_per_ovp_perMkm2+night_counts_per_ovp_perMkm2)*100,
         night_frp_pct = night_frp_per_ovp_perMkm2/(day_frp_per_ovp_perMkm2+night_frp_per_ovp_perMkm2)*100) %>%
  dplyr::select(lc_kop, millions_of_km2, 
                day_counts_per_ovp_perMkm2, night_counts_per_ovp_perMkm2, night_count_pct,
                day_frp_per_ovp_perMkm2, night_frp_per_ovp_perMkm2, night_frp_pct)


# final_final_table<- rbind(final_table, total_row) %>% 
#   dplyr::mutate_if(is.numeric, sprintf, fmt = "%.1f")

# write_csv(final_final_table,"out/table_s1.csv")
