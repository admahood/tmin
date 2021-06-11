# modeling overpass count-adjusted active fire data

# setup ==========

library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(stars)
library(mblm)
library(terra)
library(doParallel)
library(foreach)
library(mgcv)

# thresholds ===================================================================
thresholds <- read_csv("data/updated-goes-af-vpd-thresholds.csv")
burnable_lcs<- pull(thresholds, lc_name)
# functions ====================================================================
parallel_theilsen <- function(stack_df, zero_to_na=FALSE, pb =TRUE,
                              minimum_sample = 10, workers = 1){
  require(doParallel)
  require(foreach)
  cells<-unique(stack_df$cell)
  
  
  registerDoParallel(workers)
  result<-foreach(i = cells, .combine = bind_rows)%dopar%{
    
    dd<-filter(stack_df, cell==i) 
    
    if(zero_to_na) dd<-dd %>% filter(value>0)
    
    if(nrow(dd)>=minimum_sample & sum(dd$value)>0){
      mod<-mblm(value~timestep, data = dd, repeated =TRUE)
      sum<-summary(mod)
      df<-tibble(cell = dd$cell[1],
                 x = dd$x[1], 
                 y = dd$y[1],
                 p = sum$coefficients[2,4],
                 beta = sum$coefficients[2,1],
                 n = nrow(dd))
      if(pb) system(paste("echo", df[1,2], "p=", round(df[1,4],2),"b=", round(df[1,5],2)))
      return(df)
    }
  }
}

parallel_theilsen_lc <- function(stack_df, zero_to_na=FALSE, pb =TRUE,
                              minimum_sample = 10, workers = 1){
  require(doParallel)
  require(foreach)
  cells<-unique(stack_df$cell)
  
  
  registerDoParallel(workers)
  result<-foreach(i = cells, .combine = bind_rows)%dopar%{
    
    dd<-filter(stack_df, cell==i) %>%
      replace_na(list(value = 0))
    
    if(zero_to_na) dd<-dd %>% filter(value>0)
    
    if(nrow(dd)>=minimum_sample & sum(dd$value)>0){
      mod<-mblm(value~timestep, data = dd, repeated =TRUE)
      sum<-summary(mod)
      df<-tibble(cell = dd$cell[1],
                 p = sum$coefficients[2,4],
                 beta = sum$coefficients[2,1],
                 n = nrow(dd))
      if(pb) system(paste("echo", df[1,1], "p=", round(df[1,2],2),"b=", round(df[1,3],2)))
      return(df)
    }
  }
}

# file getting ============

system(paste("aws s3 sync",
             "s3://earthlab-mkoontz/MODIS-overpass-counts_1_analysis-ready",
             "data/overpass_counts"))
system(paste("aws s3 sync",
             "s3://earthlab-mkoontz/MODIS-overpass-counts_2.5_analysis-ready",
             "data/overpass_counts_2_5"))



system(paste("aws s3 sync",
             "s3://earthlab-jmcglinchy/night_fire/gridded/vars_refresh_may2021/CSV_nocorn_grid_1_0_degree_vars",
             "data/gridded_mod14"))

system(paste("aws s3 sync",
             "s3://earthlab-jmcglinchy/night_fire/gridded/vars_refresh_may2021/CSV_nocorn_grid_2_5_degree_vars",
             "data/gridded_mod14_2_5"))

op_days <- list.files("data/overpass_counts/year-month", 
                      full.names = TRUE, pattern = "*day*") %>%
  as_tibble() %>%
  mutate(ym = str_extract(value,"\\d{4}-\\d{2}"),
         year = str_sub(ym, 1,4)%>% as.numeric,
         month = str_sub(ym,6,7)%>% as.numeric)
op_nights <- list.files("data/overpass_counts/year-month", 
                        full.names = TRUE, pattern = "*night*") %>%
  as_tibble() %>%
  mutate(ym = str_extract(value,"\\d{4}-\\d{2}"),
         year = str_sub(ym, 1,4)%>% as.numeric,
         month = str_sub(ym,6,7) %>% as.numeric)

op_days_annual <- list.files("data/overpass_counts/annual", 
                      full.names = TRUE, pattern = "day") %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}") %>% as.numeric)
op_nights_annual <- list.files("data/overpass_counts/annual", 
                        full.names = TRUE, pattern = "night") %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}") %>% as.numeric)

mod14_day_counts<- list.files("data/gridded_mod14/AFC_num", 
                              full.names = TRUE,
                              pattern = "_D_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002)%>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g4","g5","g7","month","g6"), remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date))

mod14_night_counts<- list.files("data/gridded_mod14/AFC_num", 
                              full.names = TRUE,
                              pattern = "_N_") %>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g4","g5","g7","month","g6"), remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date)) 


# do the adjusting =================
dir.create("data/adjusted_counts")

# annual =====================
raw_n <- list.files("data/gridded_mod14/AFC_num", full.names = TRUE, pattern = "_N_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

raw_d <- list.files("data/gridded_mod14/AFC_num", full.names = TRUE, pattern = "_D_")%>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

raw_frp_mean <- list.files("data/gridded_mod14/FRP_mean", pattern = "_N_", full.names = TRUE) %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

raw_frp_total <- list.files("data/gridded_mod14/FRP_total", pattern = "_N_", full.names = TRUE) %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

dir.create("data/annual_adjusted_counts")
years <- 2003:2020
for(y in years){
  print(y)
  d<-filter(raw_d, year == y) %>%
    pull(value) %>%
    terra::rast()%>% 
    sum
  op_d<-filter(op_days_annual, year == y) %>%
    pull(value) %>%
    terra::rast()
  adjusted_d <-d/op_d
  terra::writeRaster(adjusted_d, 
              filename = paste0("data/annual_adjusted_counts/day_adj_annual_counts_",y, ".tif"),
              overwrite = TRUE)
  n<-filter(raw_n, year == y) %>%
    pull(value) %>%
    terra::rast() %>%
    sum
  op_n<-filter(op_nights_annual, year == y) %>%
    pull(value) %>%
    terra::rast()
  adjusted_n <- n/op_n
  terra::writeRaster(adjusted_n, 
              filename = paste0("data/annual_adjusted_counts/night_adj_annual_counts_",y, ".tif"),
              overwrite = TRUE)
  nf<- adjusted_n/(adjusted_n+adjusted_d)
  terra::writeRaster(nf, 
              filename = paste0("data/annual_adjusted_counts/night_fraction_adj_annual_counts_",y, ".tif"),
              overwrite = TRUE)
  
  frp_mean <- filter(raw_frp_mean, year == y) %>%
    pull(value) %>%
    terra::rast() %>%
    mean
  terra::writeRaster(frp_mean, 
              filename = paste0("data/annual_adjusted_counts/annual_night_frp_mean_MW_per_detection_",y, ".tif"),
              overwrite = TRUE)
  frp_total <- filter(raw_frp_total, year == y) %>%
    pull(value) %>%
    terra::rast() %>%
    mean
  adjusted_frp_total <- frp_total/op_n
  terra::writeRaster(adjusted_frp_total, 
                     filename = paste0("data/annual_adjusted_counts/annual_night_frp_MW_per_ovp_",y, ".tif"),
                     overwrite = TRUE)
  
}


# TIME SERIES ANALYSIS =========================================================
# ANNUAL =======================================================================
## annual day counts============================================================
day_counts <- list.files("data/annual_adjusted_counts", 
                         pattern = "day",
                         full.names = TRUE) %>%
  terra::rast()

dir.create("out")

dc_area <- read_stars("data/annual_adjusted_counts/day_adj_annual_counts_2003.tif") %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>%
  mutate(area_km2 = (st_area(.)/1000000)%>% as.numeric) %>%
  st_centroid() %>%
  mutate(lat = st_coordinates(.)[,2])%>%
  st_set_geometry(NULL) %>%
  dplyr::select(area_km2, lat) %>%
  unique()

# 1 degree
day_counts_df <- day_counts %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>% 
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1))

day_counts_trends<-parallel_theilsen(day_counts_df,
                                     zero_to_na = FALSE, 
                                     pb=TRUE,
                                     workers=4,
                                     minimum_sample = 10) %>%
  left_join(dc_area,by=c("y"="lat"))

save(day_counts_trends, file = "data/day_counts_trends.Rda")

p_dc<-ggplot(day_counts_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle("Annual sums of daytime active fire detections", 
          "1 degree, 2003-2018, p<0.05")+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/day_count_trend_1_deg_2003-2018.png")

## annual night counts =========================================================
night_counts <- list.files("data/annual_adjusted_counts", 
                           pattern = "night_adj_annual_counts", 
                           full.names = TRUE) %>%
  terra::rast()

# 1 degree
night_counts_df<-night_counts %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],names_to = "layer", values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1))

night_counts_trends<-parallel_theilsen(night_counts_df,
                                       zero_to_na = FALSE, 
                                       workers=4,
                                       minimum_sample = 10)%>%
  left_join(dc_area,by=c("y"="lat"))

p_nc<-ggplot(night_counts_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle("Annual sums of nighttime active fire detections",
          "1 degree, 2003-2018, p<0.05")+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_count_trend_1_deg_2003-2018.png")

save(night_counts_trends, file = "night_counts_trends.Rda")

## annual night fraction ===========
night_fractions <- list.files("data/annual_adjusted_counts", pattern = "night_fraction", full.names = TRUE) %>%
  terra::rast()

# 1 degree
night_fraction_df<-night_fractions%>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],names_to = "layer", values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1))

night_fraction_trends<-parallel_theilsen(night_fraction_df,
                                       zero_to_na = FALSE, 
                                       workers=4,
                                       minimum_sample = 10)%>%
  left_join(dc_area,by=c("y"="lat"))

p_nf<-ggplot(night_fraction_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle("Fraction of detections that occur at night", "1 degree, 2003-2018, p<0.05")+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_fraction_trend_1_deg_2003-2018.png")
save(night_fraction_trends,file="night_fraction_trends.Rda")
## annual night frp mean =======================================================

night_frp <- list.files("data/annual_adjusted_counts", 
                        pattern = "annual_night_frp_mean_MW_per_detection_", 
                        full.names = TRUE) %>%
  terra::rast()


night_frp_df<-night_frp %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,6,7) %>% as.numeric)) %>%
  replace_na(list(timestep = 1))

night_frp_trends<-parallel_theilsen(night_frp_df,
                                    zero_to_na = TRUE,
                                    workers=4, 
                                    minimum_sample = 10)%>%
  left_join(dc_area,by=c("y"="lat"))

p_frp <- ggplot(night_frp_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle("Annual mean nighttime FRP, MW per active fire detection", 
          "1 degree, 2003-2018, p<0.05")+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_frp_trend_1_deg_2003-2018.png")

save(night_frp_trends, file = "night_frp_mean_trends.Rda")

## annual night frp total ======================================================

night_frp_mw_per_ovp <- list.files("data/annual_adjusted_counts", 
                        pattern = "annual_night_frp_MW_per_ovp", 
                        full.names = TRUE) %>%
  terra::rast()


adj_night_frp_df<-night_frp_mw_per_ovp %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,6,7) %>% as.numeric)) %>%
  replace_na(list(timestep = 1))

adj_night_frp_trends<-parallel_theilsen(adj_night_frp_df,
                                    zero_to_na = FALSE,
                                    workers=4, minimum_sample = 10)%>%
  left_join(dc_area,by=c("y"="lat"))

p_adj_frp <- ggplot(adj_night_frp_trends %>% filter(p<0.05) %>% 
                  mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle("Annual Mean Nighttime FRP MW per Overpass", 
          "1 degree, 2003-2018, p<0.05")+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/adjusted_total_night_frp_trend_1_deg_2003-2018.png")

save(adj_night_frp_trends, file = "night_frp_total_trends.Rda")

# all together
library(ggpubr)
ggarrange(p_dc, p_nc, p_nf, p_adj_frp, nrow = 2, ncol=2) +
  ggsave(height = 7, width=12, filename = "out/annual_4pan.png")

# table of trends per area

bind_rows("daytime_afd"=day_counts_trends, 
          "nighttime_afd"= night_counts_trends, 
          "percent_nighttime_afd"= night_fraction_trends,
          "frp_MW_per_afd"=night_frp_trends, 
          "frp_MW_per_overpass"=adj_night_frp_trends, .id="id")%>%
  filter(p<0.05) %>% 
  mutate(Trend = ifelse(beta > 0, "positive", "negative")) %>%
  group_by(Trend, id) %>%
  summarise(area_Mkm2 = sum(area_km2)/1e6) %>%
  ungroup %>%
  pivot_wider(names_from = "id", values_from = area_Mkm2,
              names_glue = "{id}_Mkm2") %>%
  dplyr::select(Trend, daytime_afd_Mkm2, nighttime_afd_Mkm2, 
                percent_nighttime_afd_Mkm2,
                frp_MW_per_afd_Mkm2,frp_MW_per_overpass_Mkm2) %>%
  write_csv("out/land_area_trends.csv")
  

# TIME SERIES BY LANDCOVER TYPE ================================================
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

# landcover
lck_shifted <- terra::rast("in/lc_koppen_2010_mode.tif") %>%
  terra::aggregate(4, fun = "modal") %>% 
  terra::shift(dx = -179.5, dy = -0.5) %>%
  terra::crop(day_counts) # getting rid of an extra polar row
terra::writeRaster(lck_shifted,
                   "out/aggregations_2003-2020/lck_shifted_1_deg.tif",
                   overwrite=TRUE)
# terra::ext(day_counts);terra::ext(lck_shifted)

# day counts ===================================================================
day_counts <- list.files("data/annual_adjusted_counts", 
                         pattern = "day",
                         full.names = TRUE) %>%
  terra::rast()

day_counts_lc_df <- c(lck_shifted, day_counts) %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(rsum = rowSums(.[5:ncol(.)]),
         lck = as.character(lck)) %>% 
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[5:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc)) %>% # calling it cell for the function
  na.omit %>%
  group_by(cell, timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  filter(cell %in% burnable_lcs)

day_counts_trends_lc<-parallel_theilsen_lc(day_counts_lc_df,
                                     zero_to_na = FALSE, 
                                     pb=TRUE,
                                     workers=4,
                                     minimum_sample = 10)

day_counts_global_df <- c(lck_shifted, day_counts) %>%
  as.data.frame() %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(lck = as.character(lck)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc))  %>%
  filter(cell %in% burnable_lcs)%>% 
  dplyr::select(-lc,-kop, -cell) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>% # calling it cell for the function
  na.omit %>%
  group_by(timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

global_d<-mblm(value~timestep, day_counts_global_df) %>% summary()
dctlc<- bind_rows(day_counts_trends_lc,
          tibble(cell = "Global",
                 p = global_d$coefficients[2,4],
                 beta = global_d$coefficients[2,1],
                 n = nrow(day_counts_global_df)))

save(dctlc, file = "data/day_counts_trends_lc.Rda")

# night counts =================================================================
night_counts <- list.files("data/annual_adjusted_counts", 
                           pattern = "night_adj_annual_counts",
                           full.names = TRUE) %>%
  terra::rast()

night_counts_lc_df <- c(lck_shifted, night_counts) %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(rsum = rowSums(.[5:ncol(.)]),
         lck = as.character(lck)) %>% 
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[5:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc)) %>% # calling it cell for the function
  na.omit %>%
  filter(cell %in% burnable_lcs)%>%
  group_by(cell, timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

night_counts_trends_lc<-parallel_theilsen_lc(night_counts_lc_df,
                                        zero_to_na = FALSE, 
                                        pb=TRUE,
                                        workers=4,
                                        minimum_sample = 10)

night_counts_global_df <- c(lck_shifted, night_counts) %>%
  as.data.frame() %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(lck = as.character(lck)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc))  %>%
  filter(cell %in% burnable_lcs)%>% 
  dplyr::select(-lc,-kop, -cell) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>% # calling it cell for the function
  na.omit %>%
  group_by(timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

global_n<-mblm(value~timestep, night_counts_global_df) %>% summary()
nctlc<- bind_rows(night_counts_trends_lc,
                  tibble(cell = "Global",
                         p = global_n$coefficients[2,4],
                         beta = global_n$coefficients[2,1],
                         n = nrow(night_counts_global_df)))

save(nctlc, file = "data/night_counts_trends_lc.Rda")

# night_fraction ===============================================================

night_fraction <- list.files("data/annual_adjusted_counts", 
                           pattern = "night_fraction",
                           full.names = TRUE) %>%
  terra::rast()

percent_night_afd_lc_df <- c(lck_shifted, night_fraction) %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(rsum = rowSums(.[5:ncol(.)]),
         lck = as.character(lck)) %>% 
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[5:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc)) %>% # calling it cell for the function
  na.omit%>%
  filter(cell %in% burnable_lcs)%>%
  group_by(cell, timestep) %>%
  summarise(value = mean(value)*100) %>%
  ungroup()

percent_night_afd_trends_lc<-parallel_theilsen_lc(percent_night_afd_lc_df,
                                          zero_to_na = FALSE, 
                                          pb=TRUE,
                                          workers=4,
                                          minimum_sample = 10)
percent_night_afd_global_df <- c(lck_shifted, night_fraction) %>%
  as.data.frame() %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(lck = as.character(lck)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc))  %>%
  filter(cell %in% burnable_lcs)%>% 
  dplyr::select(-lc,-kop, -cell) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>% # calling it cell for the function
  na.omit %>%
  group_by(timestep) %>%
  summarise(value = mean(value)*100) %>%
  ungroup()

global_nf<-mblm(value~timestep, percent_night_afd_global_df) %>% summary()
nfctlc<- bind_rows(percent_night_afd_trends_lc,
                  tibble(cell = "Global",
                         p = global_nf$coefficients[2,4],
                         beta = global_nf$coefficients[2,1],
                         n = nrow(percent_night_afd_global_df)))

save(nfctlc, file = "data/percent_night_afd_trends_lc.Rda")

# night frp ====================================================================
night_frp_mean <- list.files("data/annual_adjusted_counts", 
                             pattern = "night_frp_mean",
                             full.names = TRUE) %>%
  terra::rast()

night_frp_mean_lc_df <- c(lck_shifted, night_frp_mean) %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(rsum = rowSums(.[5:ncol(.)]),
         lck = as.character(lck)) %>% 
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[5:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,6,7) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc)) %>% # calling it cell for the function
  na.omit%>%
  filter(cell %in% burnable_lcs)%>%
  group_by(cell, timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

night_frp_mean_trends_lc<-parallel_theilsen_lc(night_frp_mean_lc_df,
                                            zero_to_na = TRUE, 
                                            pb=TRUE,
                                            workers=4,
                                            minimum_sample = 10)

night_frp_mean_global_df <- c(lck_shifted, night_frp_mean) %>%
  as.data.frame() %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(lck = as.character(lck)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc))  %>%
  filter(cell %in% burnable_lcs)%>% 
  dplyr::select(-lc,-kop, -cell) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,6,7) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>% # calling it cell for the function
  na.omit %>%
  group_by(timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

global_nfrpm<-mblm(value~timestep, night_frp_mean_global_df) %>% summary()
nfrpmtlc<- bind_rows(night_frp_mean_trends_lc,
                     tibble(cell = "Global",
                            p = global_nfrpm$coefficients[2,4],
                            beta = global_nfrpm$coefficients[2,1],
                            n = nrow(night_frp_mean_global_df)))

save(nfrpmtlc, file = "data/night_frp_mean_trends_lc.Rda")
night_frp_mean_global_df %>% 
  mutate(year = timestep+2003) %>%
  write_csv("data/global_yearly_means_night_frp_per_afd.csv")


# night frp adjusted ===========================================================
night_frp_total <- list.files("data/annual_adjusted_counts", 
                             pattern = "night_frp_MW_per_ovp",
                             full.names = TRUE) %>%
  terra::rast()

night_frp_total_lc_df <- c(lck_shifted, night_frp_total) %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(rsum = rowSums(.[5:ncol(.)]),
         lck = as.character(lck)) %>% 
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[5:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,6,7) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc)) %>% # calling it cell for the function
  na.omit%>%
  filter(cell %in% burnable_lcs)%>%
  group_by(cell, timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

night_frp_total_trends_lc<-parallel_theilsen_lc(night_frp_total_lc_df,
                                                zero_to_na = FALSE, 
                                                pb=TRUE,
                                                workers=4,
                                                minimum_sample = 10)

night_frp_total_global_df <- c(lck_shifted, night_frp_total) %>%
  as.data.frame() %>%
  dplyr::rename(lck = lc_koppen_2010_mode) %>%
  mutate(lck = as.character(lck)) %>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         cell = paste(kop, lc))  %>%
  filter(cell %in% burnable_lcs)%>% 
  dplyr::select(-lc,-kop, -cell) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "layer", 
               values_to = "value") %>%
  mutate(timestep = 1+(str_sub(layer,6,7) %>% as.numeric)) %>%
  replace_na(list(timestep = 1)) %>% # calling it cell for the function
  na.omit %>%
  group_by(timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

global_nfrpt<-mblm(value~timestep, night_frp_total_global_df) %>% summary()
nfrpttlc<- bind_rows(night_frp_total_trends_lc,
                  tibble(cell = "Global",
                         p = global_nfrpt$coefficients[2,4],
                         beta = global_nfrpt$coefficients[2,1],
                         n = nrow(night_frp_total_global_df)))

save(nfrpttlc, file = "data/night_frp_total_trends_lc.Rda")
night_frp_total_global_df %>% 
  mutate(year = timestep+2003) %>%
  write_csv("data/global_yearly_means_night_frp_per_overpass.csv")
# making a table of lc trends ==================================================

list.files("data", pattern = "lc.Rda", full.names = TRUE) %>%
  lapply(load)

trends_table <- bind_rows(list(dctlc %>% 
                                 mutate(variable = "day_afd_per_ovp"),
                               nctlc %>% 
                                 mutate(variable = "night_afd_per_ovp"),
                               nfctlc %>% 
                                 mutate(variable = "percent_night_afd"),
                               nfrpttlc %>%
                                 mutate(variable = "night_frp_MW_per_ovp"),
                               nfrpmtlc %>% 
                                 mutate(variable = "night_frp_MW_per_afd"))) %>%
  dplyr::select(-n) %>%
  dplyr::mutate(p = ifelse(p<0.05, "*", ""),
                beta = str_c(
                  signif(beta, 3),
                  " ",
                  p)%>% trimws) %>%
  dplyr::select(-p) %>%
  pivot_wider(names_from = "variable",
              values_from = "beta") %>%
  na.omit() %>%
  dplyr::rename(landcover = cell) %>%
  mutate(landcover = replace(landcover, landcover == "Global", "zGlobal")) %>%
  arrange(landcover)
  
write_csv(trends_table, "out/annual_trends_by_lc.csv")

# table of yearly means ========================================================

yearly_means<-bind_rows("day_afd_per_ovp"=day_counts_global_df,
          "night_afd_per_ovp"=night_counts_global_df,
          "percent_night_afd"=percent_night_afd_global_df,
          'night_frp_MW_per_afd'=night_frp_mean_global_df,
          'night_frp_MW_per_ovp'=night_frp_total_global_df,
          .id = "id") %>%
  mutate(timestep = timestep + 2002) %>%
  pivot_wider(names_from = "id", values_from = value)
write_csv(yearly_means, "out/yearly_means.csv")
# MONTHLY GLOBAL ===============================================================
wide_df<-read_csv("data/mcd14ml-global-trend-by-month_wide.csv") %>%
  mutate(percent_n_night = prop_n_night*100)%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, time2 = min(year_month), units = "days")))

long_df <- read_csv("data/mcd14ml-global-trend-by-month.csv") 


day_afd_ts <- long_df %>%
  filter(dn_detect == "day")%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days"))) %>%
  mblm(n_per_op_per_Mkm2~time, data=.);summary(day_afd_ts)


night_afd_ts<-long_df %>%
  filter(dn_detect == "night")%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days"))) %>%
  mblm(n_per_op_per_Mkm2~time, data=.);summary(night_afd_ts)

nf_ts <-mblm(percent_n_night~time, data=wide_df);summary(nf_ts)

frp_ts <- long_df %>%
  filter(dn_detect == "night") %>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days"))) %>%
  mblm(mean_frp_per_detection~time, data=.);summary(frp_ts)

bind_rows(day_afd_ts, night_afd_ts, nf_ts, frp_ts)-> global_trends

# global predictions
global_preds <- tibble("day_afd_per_ovp"=predict(day_afd_ts), 
       "night_afd_per_ovp"=   predict(night_afd_ts), 
       "percent_n_afd"=  predict(nf_ts), 
       "night_frp_MW_per_afd"= predict(frp_ts)) %>%
  slice(c(1,216)) %>%
  mutate(date = c("January 2003", "December 2020"))

write_csv(global_preds, "out/global_prediction_bookends.csv")
predict(nf_ts)

# https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
nf_gam <- gamm(percent_n_night ~ 
                 s(acq_month, bs = "cc", k = 12) + 
                 s(time),
               correlation = corAR1(form = ~ time),
               data = wide_df)
plot(nf_gam$lme)
plot(nf_gam$gam)
summary(nf_gam$lme)
summary(nf_gam$gam)

day_gam <- long_df %>%
  filter(dn_detect == "day")%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days")))%>%
  gamm(n_per_op_per_Mkm2 ~
                  s(acq_month, bs = "cc", k = 12) + 
                  s(time),
       correlation = corCAR1(form = ~ time), 
                data = .)
plot(day_gam$lme)
plot(day_gam$gam)
summary(day_gam$lme)
summary(day_gam$gam)

night_gam <- long_df %>%
  filter(dn_detect == "night")%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days")))%>%
  gamm(n_per_op_per_Mkm2 ~
         s(acq_month, bs = "cc", k = 12) + 
         s(time),
       correlation = corAR1(form = ~ time), 
       data = .)
plot(night_gam$lme)
plot(night_gam$gam)
summary(night_gam$lme)
summary(night_gam$gam)

night_frp_gam <- long_df %>%
  filter(dn_detect == "night") %>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days")))%>%
  gamm(mean_frp_per_detection ~ 
                        s(acq_month, bs = "cc", k = 12) + 
                        s(time),
       correlation = corAR1(form = ~ time), 
                      data = .)

plot(night_frp_gam$lme)
plot(night_frp_gam$gam)
summary(night_frp_gam$lme)
summary(night_frp_gam$gam)

layout(matrix(1:2, ncol = 2))
acf(resid(nf_gam$lme), lag.max = 36, main = "ACF")
pacf(resid(nf_gam$lme), lag.max = 36, main = "pACF")
layout(1)

# MONTHLY by Koppen=============================================================

lut_kop<- c("Equatorial", "Arid", "Temperate", "Boreal", "Polar")
names(lut_kop) <- c(1,2,3,4,5)

bind_rows(
koppen_percent_n_night <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen_wide.csv") %>%
  mutate(value = prop_n_night*100)%>%
  dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days")),
                cell = lut_kop[koppen])%>%
  parallel_theilsen_lc() %>%
  mutate(variable = "percent_afd_night")
,
koppen_day_afd <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv")%>%
  filter(dn_detect == "day")%>%
  dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days")),
                cell = lut_kop[koppen]) %>%
  dplyr::rename(value = n_per_op_per_Mkm2)%>%
  parallel_theilsen_lc() %>%
  mutate(variable = "day_afd_per_op_per_Mkm2")
,
koppen_night_afd <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv")%>%
  filter(dn_detect == "night")%>%
  dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                               time2 = min(year_month), 
                                               units = "days")),
                cell = lut_kop[koppen]) %>%
  dplyr::rename(value = n_per_op_per_Mkm2)%>%
  parallel_theilsen_lc() %>%
  mutate(variable = "night_afd_per_op_per_Mkm2")
,
koppen_night_frp <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv")%>%
  filter(dn_detect == "night")%>%
  dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                               time2 = min(year_month), 
                                               units = "days")),
                cell = lut_kop[koppen]) %>%
  dplyr::rename(value = mean_frp_per_detection)%>%
  parallel_theilsen_lc() %>%
  mutate(variable = "mean_frp_per_detection")
)-> koppen_trends

write_csv(koppen_trends, "out/koppen_trends_raw.csv")

koppen_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         sign = ifelse(beta>0, "+", "-"),
         sign_sig = ifelse(p<0.05, sign, sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "sign_sig") %>%
  write_csv("out/koppen_trends_pretty.csv")

koppen_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta*365.25,2),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  write_csv("out/koppen_trends_pretty_year.csv")

koppen_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta,5),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  write_csv("out/koppen_trends_pretty_day.csv")

# MONTHLY by landcover =========================================================
lck_tab <- read_csv("in/csvs_from_michael/koppen-modis-landcover-lookup-table.csv")
lut_lc <- pull(lck_tab, koppen_modis_name)
names(lut_lc)<- pull(lck_tab, koppen_modis_code)

bind_rows(
  landcover_percent_n_night <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-landcover_wide.csv") %>%
    mutate(value = prop_n_night*100)%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_lc[as.character(lc)])%>%
    parallel_theilsen_lc() %>%
    mutate(variable = "percent_afd_night")
  ,
  landcover_day_afd <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-landcover.csv")%>%
    filter(dn_detect == "day")%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_lc[as.character(lc)]) %>%
    dplyr::rename(value = n_per_op_per_Mkm2)%>%
    parallel_theilsen_lc() %>%
    mutate(variable = "day_afd_per_op_per_Mkm2")
  ,
  landcover_night_afd <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-landcover.csv")%>%
    filter(dn_detect == "night")%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_lc[as.character(lc)]) %>%
    dplyr::rename(value = n_per_op_per_Mkm2)%>%
    parallel_theilsen_lc() %>%
    mutate(variable = "night_afd_per_op_per_Mkm2")
  ,
  landcover_night_frp <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-landcover.csv")%>%
    filter(dn_detect == "night")%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_lc[as.character(lc)]) %>%
    dplyr::rename(value = mean_frp_per_detection)%>%
    parallel_theilsen_lc() %>%
    mutate(variable = "mean_frp_per_detection")
)-> landcover_trends

write_csv(landcover_trends, "out/landcover_trends_raw.csv")

landcover_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         sign = ifelse(beta>0, "+", "-"),
         sign_sig = ifelse(p<0.05, sign, sig)) %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "sign_sig") %>%
  write_csv("out/landcover_trends_pretty.csv")

landcover_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta*365.25,2),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  write_csv("out/landcover_trends_pretty_year.csv")

landcover_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta,5),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  write_csv("out/landcover_trends_pretty_day.csv")
