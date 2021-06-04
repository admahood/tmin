# modeling overpass count-adjusted active fire data

# setup ==========

library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(mblm)
library(terra)
library(doParallel)
library(foreach)

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
    
    dd<-filter(stack_df, cell==i) 
    
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

# monthly ==========================
registerDoParallel(detectCores()/2)
foreach(i = 1:nrow(mod14_day_counts))%dopar%{
  system(paste("echo",i))
  monthi <- mod14_day_counts$month_n[i]
  yeari <- mod14_day_counts$year[i]
  
  counts <- filter(mod14_day_counts, year == yeari, month_n == monthi) %>%
    pull(value) %>%
    raster
  overpasses<- filter(op_days, year == yeari, month == monthi) %>%
    pull(value) %>%
    raster
  
  adjusted <- counts/overpasses
  outfile <- paste0("data/adjusted_counts/op-adjusted_D_", 
                   yeari,
                   str_pad(monthi, width=2, side="left", pad = "0"),
                   ".tif"
                   )
  writeRaster(adjusted,outfile, overwrite=TRUE)  

  counts_n <- filter(mod14_night_counts, year == yeari, month_n == monthi) %>%
    pull(value) %>%
    raster
  overpasses_n<- filter(op_nights, year == yeari, month == monthi) %>%
    pull(value) %>%
    raster
  
  adjusted_n <- counts_n/overpasses_n
  outfile_n <- paste0("data/adjusted_counts/op-adjusted_N_", 
                    yeari,
                    str_pad(monthi, width=2, side="left", pad = "0"),
                    ".tif"
  )
  writeRaster(adjusted_n,outfile_n, overwrite=TRUE)
  
  if(str_extract(outfile, "\\d{6}")==str_extract(outfile_n, "\\d{6}")){
  
    nf<- adjusted_n/(adjusted+adjusted_n)
    outfile_nf <- paste0("data/adjusted_counts/op-adjusted_NF_", 
                        yeari,
                        str_pad(monthi, width=2, side="left", pad = "0"),
                        ".tif"
    )
    writeRaster(nf,outfile_nf, overwrite=TRUE)
  }
}

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
                                     minimum_sample = 10)
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
                                       minimum_sample = 10)

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
                                       minimum_sample = 10)

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
                                    minimum_sample = 10)

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
                                    workers=4, minimum_sample = 10)

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

# all together
library(ggpubr)
ggarrange(p_dc, p_nc, p_nf, p_adj_frp, nrow = 2, ncol=2) +
  ggsave(height = 7, width=12, filename = "out/annual_4pan.png")


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
  filter(cell %in% burnable_lcs) %>%
  group_by(cell, timestep) %>%
  summarise(value = mean(value)) %>%
  ungroup()

day_counts_trends_lc<-parallel_theilsen_lc(day_counts_lc_df,
                                     zero_to_na = FALSE, 
                                     pb=TRUE,
                                     workers=4,
                                     minimum_sample = 10)
save(day_counts_trends_lc, file = "data/day_counts_trends_lc.Rda")

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
save(night_counts_trends_lc,file= "data/night_counts_trends_lc.Rda")

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
save(percent_night_afd_trends_lc,file= "data/night_fraction_trends_lc.Rda")

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
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
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

save(night_frp_mean_trends_lc,file= "data/night_frp_mean_trends_lc.Rda")

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
  mutate(timestep = 1+(str_sub(layer,5,6) %>% as.numeric)) %>%
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

save(night_frp_total_trends_lc,file= "data/night_frp_total_trends_lc.Rda")

# making a table of lc trends ==================================================

list.files("data", pattern = "lc.Rda", full.names = TRUE) %>%
  lapply(load)

trends_table <- bind_rows(list(day_counts_trends_lc %>% 
                                 mutate(variable = "day_afd_per_ovp"),
                               night_counts_trends_lc %>% 
                                 mutate(variable = "night_afd_per_ovp"),
                               percent_night_afd_trends_lc %>% 
                                 mutate(variable = "percent_night_afd"),
                               night_frp_total_trends_lc %>%
                                 mutate(variable = "night_frp_MW_per_ovp"),
                               night_frp_mean_trends_lc %>% 
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
  arrange(landcover)
  
write_csv(trends_table, "out/annual_trends_by_lc.csv")
# MONTHLY ======================================================================

workers<- detectCores()/2
## day counts ===========
dc_m <- list.files("data/adjusted_counts",
                          pattern = "_D_", full.names = TRUE)%>%
  as_tibble  %>%
  mutate(ym = str_extract(value, "\\d{6}") %>% as.numeric,
         year = str_sub(ym, 1,4),
         month = str_sub(ym, 5,6)) %>%
  filter(year > 2002) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%m-%d"),
         month_n = lubridate::month(date)) %>%
  dplyr::arrange(date)

dc_stack_m<-terra::rast(pull(dc_m, value))

dc_df_m <- dc_stack_m  %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],names_to = "layer", values_to = "value")  %>%
  mutate(ym = str_extract(layer, "\\d{6}") %>% as.numeric,
         year = str_sub(ym, 1,4),
         month = str_sub(ym, 5,6)) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%m-%d")) %>%
  dplyr::arrange(cell,date) %>%
  group_by(cell) %>%
  mutate(timestep = row_number()) %>%
  ungroup()


workers=4
dc_trends_m<-parallel_theilsen(dc_df_m, pb=TRUE,
                                      zero_to_na = FALSE,
                                      workers=workers, 
                                      minimum_sample = 30)

p_dc_m <- ggplot(dc_trends_m %>% 
                  filter(p<0.05->alpha) %>% 
                  mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Monthly daytime detections, 1 degree, 2003-2020, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/monthly_day_count_trend_1_deg_2003-2020.png")

## night counts ===========
nc_m <- list.files("data/adjusted_counts",
                   pattern = "_N_", full.names = TRUE)%>%
  as_tibble  %>%
  mutate(ym = str_extract(value, "\\d{6}") %>% as.numeric,
         year = str_sub(ym, 1,4),
         month = str_sub(ym, 5,6)) %>%
  filter(year > 2002) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%m-%d"),
         month_n = lubridate::month(date)) %>%
  dplyr::arrange(date)

nc_stack_m<- terra::rast(pull(nc_m, value))

nc_df_m <- nc_stack_m  %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],names_to = "layer", values_to = "value")  %>%
  mutate(ym = str_extract(layer, "\\d{6}") %>% as.numeric,
         year = str_sub(ym, 1,4),
         month = str_sub(ym, 5,6)) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%m-%d")) %>%
  dplyr::arrange(cell,date) %>%
  group_by(cell) %>%
  mutate(timestep = row_number()) %>%
  ungroup()

nc_trends_m<-parallel_theilsen(nc_df_m, pb=TRUE,
                               zero_to_na = FALSE,
                               workers=workers, 
                               minimum_sample = 30)

p_nc_m <- ggplot(nc_trends_m %>% 
                   filter(p<0.05->alpha) %>% 
                   mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Monthly nighttime detections, 1 degree, 2003-2020, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/monthly_night_count_trend_1_deg_2003-2020.png")

## night fraction ===========
nf_df_m <- list.files("data/adjusted_counts",
                   pattern = "_NF_", full.names = TRUE)%>%
  terra::rast() %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  # mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  # filter(rsum>0) %>%
  # dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],names_to = "layer", values_to = "value")  %>%
  mutate(ym = str_extract(layer, "\\d{6}") %>% as.numeric,
         year = str_sub(ym, 1,4),
         month = str_sub(ym, 5,6)) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%m-%d")) %>%
  dplyr::arrange(cell,date) %>%
  group_by(cell) %>%
  mutate(timestep = row_number()) %>%
  ungroup()

nf_trends_m<-parallel_theilsen(nf_df_m, pb=TRUE,
                               zero_to_na = FALSE,
                               workers=workers, 
                               minimum_sample = 30)

p_nf_m <- ggplot(nf_trends_m %>% 
                   filter(p<0.05->alpha) %>% 
                   mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Monthly fraction of nighttime detections, 1 degree, 2003-2020, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/monthly_night_fraction_trend_1_deg_2003-2020.png")

## night frp ===========
# need to get the files in the correct order
night_frp_m <- list.files("data/gridded_mod14/FRP_mean",
                          pattern = "_N_", full.names = TRUE)%>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  tidyr::separate(value, sep = "_",
           into = c("g1", "g2","g3","g4","g5","g7","month","g6"), 
           remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date)) %>%
  dplyr::arrange(date)

night_frp_stack_m <- terra::rast(pull(night_frp_m, value))

night_frp_df_m <- night_frp_stack_m %>%
  as.data.frame(xy=TRUE, cells=TRUE) %>%
  mutate(rsum = rowSums(.[4:ncol(.)])) %>%
  filter(rsum>0) %>%
  dplyr::select(-rsum) %>%
  pivot_longer(cols = names(.)[4:ncol(.)],names_to = "layer", values_to = "value")%>%
  tidyr::separate(layer, sep = "_",
                  into = c("g1", "g2","g3","g4","month","g6"), 
                  remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date))%>%
  dplyr::arrange(cell,date) %>%
  group_by(cell) %>%
  mutate(timestep = row_number()) %>%
  ungroup() %>%
  filter(value>0)

night_frp_trends_m <- parallel_theilsen(night_frp_df_m,
                                      zero_to_na = TRUE,
                                      workers=workers, 
                                      minimum_sample = 30)

p_frp_m <- ggplot(night_frp_trends_m %>% 
                  filter(p<0.05->alpha) %>% 
                  mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Monthly nighttime FRP, 1 degree, 2003-2020, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_frp_trend_1_deg_2003-2020.png")


# 4pan
library(ggpubr)
ggarrange(p_dc_m, p_nc_m, p_nf_m, p_frp_m) +
  ggsave(filename = "out/monthly_4pan.png", height = 7, width = 12)
