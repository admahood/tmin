# modeling overpass count-adjusted active fire data

# setup ==========

library(sf)
library(cowplot)
library(raster)
library(tidyverse)
library(lubridate)
library(stars)
library(mblm)
library(terra)
library(doParallel)
library(foreach)
library(mgcv)
library(ggthemes)
library(viridis)
library(ggpubr)
library(broom)

posneg_cols <- c("Positive"="#CF6630", "Negative"="#07484D") # colours.cafe palette 582
mask_col<- #DDCEBF
daynight_cols <- c("#B2182B","#2166AC") # red is #B2182B
dir.create("out")

# thresholds ===================================================================
thresholds <- read_csv("in/csvs_from_michael/zero-goes-af-vpd-thresholds-with-landcover-codes.csv")
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

parallel_gamm_lc <- function(stack_df, zero_to_na=FALSE, pb =TRUE,
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
      mod<- gamm(value ~ 
               s(acq_month, bs = "cc", k = 12) + 
               s(timestep),
             correlation = corAR1(form = ~ timestep), 
             data = dd)
      sum<-summary(mod$lme)
      df<-tibble(cell = dd$cell[1],
                 p = sum$tTable[2,5],
                 beta = sum$coefficients$fixed[2],
                 n = nrow(dd))
      #if(pb) system(paste("echo", df[1,1], "p=", round(df[1,2],2),"b=", round(df[1,3],2)))
      return(df)
    }
  }
}
# lc koppen setup ==============================================================
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
template <- list.files("data/annual_adjusted_counts", 
                         pattern = "day",
                         full.names = TRUE) %>%
  terra::rast()

lck_shifted <- terra::rast("in/lc_koppen_2010_mode.tif") %>%
  terra::aggregate(4, fun = "modal") %>% 
  terra::shift(dx = -179.5, dy = -0.5) %>%
  terra::crop(template) # getting rid of an extra polar row

# creating the burnable land mask ==============================================
lck_s <- terra::rast("out/aggregations_2003-2020/lck_shifted.tif")
lck_s[lck_s==100] <- NA
lck_s[lck_s==200] <- NA
lck_s[lck_s==300] <- NA
lck_s[lck_s==400] <- NA
lck_s[lck_s==500] <- NA
brnbl<-thresholds %>% 
  dplyr::select(lc_kop=lc_name) %>% 
  mutate(Burnable = "yes")

lck_df<-as.data.frame(lck_shifted, xy=TRUE)%>%
  mutate(Koppen = lut_kop[str_sub(lc_koppen_2010_mode,1,1)],
         Landcover = lut_lc[str_sub(lc_koppen_2010_mode,2,3)],
         lc_kop = paste(Koppen, Landcover)) %>%
  left_join(brnbl) %>%
  replace_na(list(Burnable = "no")) %>%
  mutate(x = ifelse(x< 0, x+360, x)-181)

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

# pointing to the overpass counts ==============================================
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

# pointing to the monthly grids with month =====================================
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


# pointing to the raw monthly grids ============================================
raw_n <- list.files("data/gridded_mod14/AFC_num", 
                    full.names = TRUE, pattern = "_N_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

raw_d <- list.files("data/gridded_mod14/AFC_num", 
                    full.names = TRUE, pattern = "_D_")%>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

raw_frp_mean <- list.files("data/gridded_mod14/FRP_mean",
                           pattern = "_N_", full.names = TRUE) %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

raw_frp_total <- list.files("data/gridded_mod14/FRP_total", 
                            pattern = "_N_", full.names = TRUE) %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

# monthly adjust ===============================================================
dir.create("data/monthly_adjusted_counts")
ym <- op_days %>%
  filter(year>2002) %>%
  pull(ym)
for(y in ym){
  print(y)
  year_i <- str_extract(y, "\\d{4}")
  month_i <- str_extract(y, "-\\d{2}") %>%
    str_remove("-")
  
  d<-filter(mod14_day_counts, 
            year == year_i &
            month_n == as.numeric(month_i)) %>%
    pull(value) %>%
    terra::rast()
  
  op_d<-filter(op_days, 
               year == year_i &
               month == as.numeric(month_i)) %>%
    pull(value) %>%
    terra::rast()
  
  adjusted_d <-d/op_d
  terra::writeRaster(adjusted_d, 
                     filename = paste0("data/monthly_adjusted_counts/day_afd_per_overpass_",
                                       year_i,"-",month_i, ".tif"),
                     overwrite = TRUE)
  n<-filter(mod14_night_counts, 
            year == year_i &
              month_n == as.numeric(month_i)) %>%
    pull(value) %>%
    terra::rast()
  
  op_n<-filter(op_nights, 
               year == year_i &
                 month == as.numeric(month_i)) %>%
    pull(value) %>%
    terra::rast()
  
  adjusted_n <-n/op_n
  terra::writeRaster(adjusted_n, 
                     filename = paste0("data/monthly_adjusted_counts/night_afd_per_overpass_",
                                       year_i,"-",month_i, ".tif"),
                     overwrite = TRUE)
  nf<- adjusted_n/(adjusted_n+adjusted_d) * 100
  terra::writeRaster(nf, 
                     filename = paste0("data/monthly_adjusted_counts/percent_night_afd_",
                                       year_i,"-",month_i, ".tif"),
                     overwrite = TRUE)
}


# annual sum and adjust ========================================================
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
## annual day counts ===========================================================
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


if(!file.exists("data/day_counts_trends.Rda")){
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
  }else{
  load("data/day_counts_trends.Rda")}

p_dc<-ggplot(day_counts_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "Positive", "Negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25, fill="white")+
  geom_raster(data = lck_df%>%
                filter(Burnable == "yes"), 
              aes(x=x,y=y), fill = "grey90")+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = posneg_cols)+
  theme_void()+
  labs(caption="Daytime Active Fire Detections")+
  ylim(c(-52.625, 75.125)) + 
  theme(text = element_text(family="DejaVuSans"),
        plot.caption = element_text(hjust=0.5, size=rel(1.2)),
        legend.position = "none",
        legend.justification = c(0,0),
        panel.background = element_rect(color = "black", size=1)) +
  ggsave("out/day_count_trend_1_deg_2003-2020.png")

## annual night counts =========================================================
if(!file.exists("data/night_counts_trends.Rda")){
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
                                         minimum_sample = 10) %>%
    left_join(dc_area,by=c("y"="lat"))
  
  save(night_counts_trends, file="data/night_counts_trends.Rda")
  }else{
  load("data/night_counts_trends.Rda")}

p_nc<-ggplot(night_counts_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "Positive", "Negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25, fill="white")+
  geom_raster(data = lck_df%>%
                filter(Burnable == "yes"), 
              aes(x=x,y=y), fill = "grey90")+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = posneg_cols)+
  theme_void()+
  labs(caption="Nighttime Active Fire Detections") +
  ylim(c(-52.625, 75.125)) + 
  theme(text = element_text(family="DejaVuSans"),
        plot.caption = element_text(hjust=0.5, size=rel(1.2)),
        legend.position = "none",
        legend.justification = c(0,0),
        panel.background = element_rect(color = "black", size=1)) +
  ggsave("out/night_count_trend_1_deg_2003-2020.png")

save(night_counts_trends, file = "night_counts_trends.Rda")

## annual night fraction ===========
if(!file.exists("data/night_fraction_trends.Rda")){
  
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
  
  save(night_fraction_trends, file= "data/night_fraction_trends.Rda")
}else{
  load("data/night_fraction_trends.Rda")}

p_nf<-ggplot(night_fraction_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "Positive", "Negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25, fill="white")+
  geom_raster(data = lck_df%>%
                filter(Burnable == "yes"), 
              aes(x=x,y=y), fill = "grey90")+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = posneg_cols)+
  theme_void()+
  labs(caption="Percent Nighttime Active Fire Detections")+
  ylim(c(-52.625, 75.125)) + 
  theme(text = element_text(family="DejaVuSans"),
        plot.caption = element_text(hjust=0.5, size=rel(1.2)),
        legend.position = "none",
        legend.justification = c(0,0),
        panel.background = element_rect(color = "black", size=1)) +
  ggsave("out/night_fraction_trend_1_deg_2003-2020.png")
save(night_fraction_trends,file="night_fraction_trends.Rda")
## annual night frp mean =======================================================
if(!file.exists("data/night_frp_trends.Rda")){
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
  
  save(night_frp_trends,file= "data/night_frp_trends.Rda")
}else{
  load("data/night_frp_trends.Rda")
}
p_frp <- ggplot(night_frp_trends %>% filter(p<0.05) %>% 
         mutate(Trend = ifelse(beta > 0, "Positive", "Negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25, fill="white")+
  geom_raster(data = lck_df %>%
                filter(Burnable == "yes"), 
              aes(x=x,y=y), fill = "grey90")+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = posneg_cols)+
  theme_void()+
  labs(caption="Nighttime Fire Radiative Power")+
  ylim(c(-52.625, 75.125)) + 
  theme(text = element_text(family="DejaVuSans"),
        plot.caption = element_text(hjust=0.5, size=rel(1.2)),
        legend.position = "none",
        legend.justification = c(0,0),
        panel.background = element_rect(color = "black", size=1)) +
  ggsave("out/night_frp_trend_1_deg_2003-2020.png")

## 4pan plot ======================================================

df_legend <- tibble(x=c(1,2,3), y=c(1,2,3), 
                    Trend=c("Negative","Positive","Burnable, but\nNot Significant"))  %>%
  ggplot(aes(x=x,y=y,fill=Trend)) +
  geom_raster() +
  theme_transparent()+
  # theme(legend.title = element_blank())+
  scale_fill_manual(values = c("grey90",posneg_cols[2], posneg_cols[1])) 

trend_leg<-get_legend(df_legend)

# all together
p_trends<-ggarrange(p_dc, p_nc, p_nf, p_frp,
          nrow = 2, ncol=2,
          labels = c("a", "b", "c", "d"), 
          font.label = "bold",
          label.y = 0.99)

ggdraw()+
  draw_plot(p_trends) +
  draw_plot(df_legend, .015,.12,.1,.1) +
ggsave(height = 5, width=12, filename = "out/annual_4pan.png")

# table of trends per area =====================================================

bind_rows("daytime_afd"=day_counts_trends, 
          "nighttime_afd"= night_counts_trends, 
          "percent_nighttime_afd"= night_fraction_trends,
          "frp_MW_per_afd"=night_frp_trends,.id = "id")%>%
  filter(p<0.05) %>% 
  mutate(Trend = ifelse(beta > 0, "positive", "negative")) %>%
  group_by(Trend, id) %>%
  summarise(area_Mkm2 = sum(area_km2)/1e6) %>%
  ungroup %>%
  pivot_wider(names_from = "id", values_from = area_Mkm2,
              names_glue = "{id}_Mkm2") %>%
  dplyr::select(Trend, daytime_afd_Mkm2, 
                nighttime_afd_Mkm2, 
                percent_nighttime_afd_Mkm2,
                frp_MW_per_afd=frp_MW_per_afd_Mkm2) %>%
  write_csv("out/land_area_trends.csv")





# MONTHLY GLOBAL ===============================================================
wide_df<-read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-month_wide.csv") %>%
  mutate(percent_n_night = prop_n_night*100)%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, time2 = min(year_month), units = "days")))

long_df <- read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-month.csv") 


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

day_frp_ts <- long_df %>%
  filter(dn_detect == "day") %>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, 
                                           time2 = min(year_month), 
                                           units = "days"))) %>%
  mblm(mean_frp_per_detection~time, data=.);summary(day_frp_ts)


global_row<-bind_rows(tidy(day_afd_ts)%>% mutate(variable = "day_afd_per_op_per_Mkm2"), 
          tidy(night_afd_ts)%>% mutate(variable = "night_afd_per_op_per_Mkm2"), 
          tidy(nf_ts)%>% mutate(variable = "percent_afd_night"), 
          tidy(frp_ts)%>% mutate(variable = "mean_frp_per_night_detection"),
          tidy(day_frp_ts)%>% mutate(variable = "mean_frp_per_day_detection")) %>%
  filter(term == "time")%>%
  mutate(sig = ifelse(p.value<0.05, "*", ""),
         bs = paste(round(estimate*365.25,2),sig))  %>%
  pivot_wider(id_cols = term, names_from = "variable", values_from = "bs") %>%
  dplyr::rename(cell=term)

global_row %>%
  write_csv("out/global_trends_pretty_year.csv")

global_row_simple <- bind_rows(tidy(day_afd_ts)%>% mutate(variable = "day_afd_per_op_per_Mkm2"), 
                      tidy(night_afd_ts)%>% mutate(variable = "night_afd_per_op_per_Mkm2"), 
                      tidy(nf_ts)%>% mutate(variable = "percent_afd_night"), 
                      tidy(frp_ts)%>% mutate(variable = "mean_frp_per_night_detection"),
                      tidy(day_frp_ts)%>% mutate(variable = "mean_frp_per_day_detection")) %>%
  filter(term == "time")%>%
  mutate(sig = ifelse(p.value<0.05, "*", ""),
         sign = ifelse(estimate>0, "+", "-"),
         sign_sig = ifelse(p.value<0.05, sign, sig))  %>%
  dplyr::select(-sign, -sig,-p.value, -estimate) %>%
  pivot_wider(id_cols = term, names_from = "variable", values_from = "sign_sig") %>%
  dplyr::rename(cell=term)

global_row_simple %>%
  write_csv("out/global_trends_simple_year.csv")

# global predictions
global_preds <- tibble("day_afd_per_ovp"=predict(day_afd_ts), 
       "night_afd_per_ovp"=   predict(night_afd_ts), 
       "percent_n_afd"=  predict(nf_ts), 
       "night_frp_MW_per_afd"= predict(frp_ts)) %>%
  slice(c(1,216)) %>%
  mutate(date = c("January 2003", "December 2020"))

write_csv(global_preds, "out/global_prediction_bookends.csv")
xxxx<-predict(nf_ts, interval = "confidence")
(xxxx[216,1]-xxxx[12,1])/18


# global gams =========
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

day_frp_gam <- long_df %>%
  filter(dn_detect == "day") %>%
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

# MONTHLY by Koppen ============================================================

lut_kop<- c("Equatorial", "Arid", "Temperate", "Boreal", "Polar")
names(lut_kop) <- c(1,2,3,4,5)

# ts =====
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
  mutate(variable = "mean_frp_per_night_detection")
,
koppen_day_frp <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv")%>%
  filter(dn_detect == "day")%>%
  dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                               time2 = min(year_month), 
                                               units = "days")),
                cell = lut_kop[koppen]) %>%
  dplyr::rename(value = mean_frp_per_detection)%>%
  parallel_theilsen_lc() %>%
  mutate(variable = "mean_frp_per_day_detection")
)-> koppen_trends

write_csv(koppen_trends, "out/koppen_trends_raw.csv")

koppen_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         sign = ifelse(beta>0, "+", "-"),
         sign_sig = ifelse(p<0.05, sign, sig))  %>%
  dplyr::select(-sign, -sig,-p, -beta) %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "sign_sig") %>%
  bind_rows(global_row_simple)%>%
  mutate(cell=replace(cell, cell=="time", "Global")) %>%
  dplyr::rename(Koppen = cell) %>%
  write_csv("out/koppen_trends_pretty.csv")

koppen_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta*365.25,2),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  bind_rows(global_row)%>%
  mutate(cell=replace(cell, cell=="time", "Global")) %>%
  dplyr::rename(Koppen = cell) %>%
  write_csv("out/koppen_trends_pretty_year.csv")

koppen_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta,5),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  write_csv("out/koppen_trends_pretty_day.csv")




# gam ====
bind_rows(
  koppen_percent_n_night <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen_wide.csv") %>%
    mutate(value = prop_n_night*100)%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_kop[koppen])%>%
    parallel_gamm_lc() %>%
    mutate(variable = "percent_afd_night")
  ,
  koppen_day_afd <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv")%>%
    filter(dn_detect == "day")%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_kop[koppen]) %>%
    dplyr::rename(value = n_per_op_per_Mkm2)%>%
    parallel_gamm_lc() %>%
    mutate(variable = "day_afd_per_op_per_Mkm2")
  ,
  koppen_night_afd <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv")%>%
    filter(dn_detect == "night")%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_kop[koppen]) %>%
    dplyr::rename(value = n_per_op_per_Mkm2)%>%
    parallel_gamm_lc() %>%
    mutate(variable = "night_afd_per_op_per_Mkm2")
  ,
  koppen_night_frp <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv")%>%
    filter(dn_detect == "night")%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_kop[koppen]) %>%
    dplyr::rename(value = mean_frp_per_detection)%>%
    parallel_gamm_lc() %>%
    mutate(variable = "mean_frp_per_detection")
)-> koppen_trends_gamm

koppen_trends_gamm %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         sign = ifelse(beta>0, "+", "-"),
         sign_sig = ifelse(p<0.05, sign, sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "sign_sig") %>%
  write_csv("out/koppen_trends_pretty_gamm.csv")

koppen_trends_gamm %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta*365.25,2),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  write_csv("out/koppen_trends_pretty_year_gamm.csv")
##

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
    mutate(variable = "mean_frp_per_night_detection")
  ,
  landcover_night_frp <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-landcover.csv")%>%
    filter(dn_detect == "day")%>%
    dplyr::mutate(timestep = as.numeric(difftime(time1 = year_month, 
                                                 time2 = min(year_month), 
                                                 units = "days")),
                  cell = lut_lc[as.character(lc)]) %>%
    dplyr::rename(value = mean_frp_per_detection)%>%
    parallel_theilsen_lc() %>%
    mutate(variable = "mean_frp_per_day_detection")
)-> landcover_trends

write_csv(landcover_trends, "out/landcover_trends_raw.csv")

landcover_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         sign = ifelse(beta>0, "+", "-"),
         sign_sig = ifelse(p<0.05, sign, sig)) %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "sign_sig")%>%
  bind_rows(global_row_simple)%>%
  mutate(cell=replace(cell, cell=="time", "Global")) %>%
  dplyr::rename(Koppen = cell) %>%
  write_csv("out/landcover_trends_pretty.csv")

landcover_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta*365.25,2),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs")%>%
  bind_rows(global_row)%>%
  mutate(cell=replace(cell, cell=="time", "Global")) %>%
  dplyr::rename(Koppen = cell) %>%
  write_csv("out/landcover_trends_pretty_year.csv")

landcover_trends %>%
  dplyr::select(-n) %>%
  mutate(sig = ifelse(p<0.05, "*", ""),
         bs = paste(round(beta,5),sig))  %>%
  pivot_wider(id_cols = cell, names_from = "variable", values_from = "bs") %>%
  write_csv("out/landcover_trends_pretty_day.csv")

# MONTHLY global line plots ====================================================

wide_df<-read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-month_wide.csv") %>%
  mutate(percent_n_night = prop_n_night*100)%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, time2 = min(year_month), units = "days")))

long_df <- read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-month.csv") 


preds_afd <- predict(day_afd_ts, interval = "confidence") %>%
  as.data.frame() %>%
  mutate(dn_detect = "day") %>%
  rbind(.,predict(night_afd_ts,  interval = "confidence") %>%
          as.data.frame() %>%
          mutate(dn_detect = "night")
        ) %>%
  mutate(year_month = pull(long_df%>% arrange(dn_detect), year_month))

preds_nf <- predict(nf_ts, interval = "confidence") %>%
  as.data.frame() %>%
  mutate(year_month = pull(wide_df, year_month))

preds_frp <- predict(day_frp_ts, interval = "confidence") %>%
  as.data.frame() %>%
  mutate(dn_detect = "day") %>%
  rbind(.,predict(frp_ts,  interval = "confidence") %>%
          as.data.frame() %>%
          mutate(dn_detect = "night")
  ) %>%
  mutate(year_month = pull(long_df%>% arrange(dn_detect), year_month))

p_afd <- ggplot(long_df, aes(x = year_month, y = n_per_op_per_Mkm2, color = dn_detect)) +
  geom_line(alpha=0.75) +
  geom_line(data = preds_afd, aes(y=fit), lwd=1)+
  geom_line(data = preds_afd, aes(y=upr), lty=2)+
  geom_line(data = preds_afd, aes(y=lwr), lty=2)+
  scale_color_manual(values=daynight_cols)+
  theme_clean() +
  ggtitle("Active Fire Detections per Overpass per Mkm2")+
  theme(legend.position = "none",
        axis.title = element_blank())

p_np <- ggplot(wide_df, aes(x=year_month, y=percent_n_night)) +
  geom_line(alpha=0.75) +
  geom_line(data = preds_nf, aes(y=fit), lwd=1)+
  geom_line(data = preds_nf, aes(y=upr), lty=2)+
  geom_line(data = preds_nf, aes(y=lwr), lty=2)+
  theme_clean() +
  ggtitle("Percent of Detections that Occurred at Night") +
  theme(axis.title = element_blank())

p_frp <- ggplot(long_df, aes(x = year_month, y = mean_frp_per_detection, color = dn_detect)) +
  geom_line(alpha=0.75) +
  geom_line(data = preds_frp, aes(y=fit), lwd=1)+
  geom_line(data = preds_frp, aes(y=upr), lty=2)+
  geom_line(data = preds_frp, aes(y=lwr), lty=2)+
  scale_color_manual(values=daynight_cols)+
  theme_clean() +
  ggtitle("Fire Radiative Power per Detection (MW)")+
  xlab("Date (Monthly Increments)")+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        axis.title.y = element_blank())

ggarrange(p_afd, p_np, p_frp, nrow=3) +
  ggsave(filename = "out/global_trends_line_plots.png",
         height =12, width = 5)

# MONTHLY koppen line plots ====================================================


wide_df_k<-read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen_wide.csv") %>%
  mutate(percent_n_night = prop_n_night*100)%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, time2 = min(year_month), units = "days"))) %>%
  mutate(koppen = lut_kop[koppen])

long_df_k <- read_csv("in/csvs_from_michael/mcd14ml-trend-by-month-koppen.csv") %>%
  mutate(koppen = lut_kop[koppen])

p_k_afd <- ggplot(long_df_k, aes(x = year_month, y = n_per_op_per_Mkm2, color = dn_detect)) +
  geom_line(alpha=0.75) +
  geom_smooth(method="lm")+
  scale_color_manual(values=daynight_cols)+
  facet_wrap(~koppen, nrow=1, ncol=4, scales = "free")+
  theme_clean() +
  ggtitle("Active Fire Detections per Overpass per Mkm2")+
  theme(legend.position = "none",
        axis.title = element_blank())

p_k_np <- ggplot(wide_df_k, aes(x=year_month, y=percent_n_night)) +
  geom_line(alpha=0.75) +
  theme_clean() +
  facet_wrap(~koppen, nrow=1, ncol=4, scales = "free")+
  geom_smooth(method="lm")+
  ggtitle("Percent of Detections that Occurred at Night") +
  theme(axis.title = element_blank())

p_k_frp <- ggplot(long_df_k, aes(x = year_month, y = mean_frp_per_detection, color = dn_detect)) +
  geom_line(alpha=0.75) + 
  facet_wrap(~koppen, nrow=1, ncol=4, scales = "free")+
  geom_smooth(method="lm")+
  scale_color_manual(values=daynight_cols)+
  theme_clean() +
  ggtitle("Fire Radiative Power per Detection (MW)")+
  xlab("Date (Monthly Increments)")+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        axis.title.y = element_blank())

ggarrange(p_k_afd, p_k_np, p_k_frp, nrow=3) +
  ggsave(filename = "out/global_trends_by_koppen_line_plots.png",
         height =12, width =12)

# time series decomposition ====================================================
library(forecast)
# devtools:: install_github("brisneve/ggplottimeseries")
library(ggplottimeseries)

wide_df<-read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-month_wide.csv") %>%
  mutate(percent_n_night = prop_n_night*100)%>%
  dplyr::mutate(time = as.numeric(difftime(time1 = year_month, time2 = min(year_month), units = "days")))

long_df <- read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-month.csv") 

df_nf <- dts1(pull(wide_df, year_month),
              pull(wide_df, percent_n_night),
              12, type = "additive")
deco_nf <- ggdecompose(df_nf)+
  xlab("Date")+
  ylab("Percent of Active Fire Detctions at Night")

df_afd_day <- dts1(pull(long_df%>% filter(dn_detect == "day"), year_month),
                   pull(long_df%>% filter(dn_detect == "day"), n_per_op_per_Mkm2),
                   12, type = "additive")

deco_day_afd <- ggdecompose(df_afd_day)+
  xlab("Date")+
  ylab("Daytime Active Fire Detctions per Overpass per Mkm2")

df_afd_night <- dts1(pull(long_df%>% filter(dn_detect == "night"), year_month),
                   pull(long_df%>% filter(dn_detect == "night"), n_per_op_per_Mkm2),
                   12, type = "additive")

deco_night_afd <- ggdecompose(df_afd_night)+
  xlab("Date")+
  ylab("Nighttime Active Fire Detctions per Overpass per Mkm2")

df_frp_night <- dts1(pull(long_df%>% filter(dn_detect == "night"), year_month),
                     pull(long_df%>% filter(dn_detect == "night"), mean_frp_per_detection),
                     12, type = "additive")

deco_night_frp <- ggdecompose(df_frp_night)+
  xlab("Date")+
  ylab("Mean FRP per Nighttime Detection")

df_frp_day <- dts1(pull(long_df%>% filter(dn_detect == "day"), year_month),
                     pull(long_df%>% filter(dn_detect == "day"), mean_frp_per_detection),
                     12, type = "additive")

deco_day_frp <- ggdecompose(df_frp_day)+
  xlab("Date")+
  ylab("Mean FRP per Daytime Detection")


ggarrange(deco_day_afd, deco_night_afd, deco_day_frp, deco_night_frp, deco_nf, nrow=3, ncol=2) +
  ggsave("out/ts_decompose.png", width = 12, height =12)





# night fraction with inset frp trends =========================================
load("data/night_fraction_trends.Rda")

wide_df<-read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-year_wide.csv") %>%
  mutate(percent_n_night = prop_n_night*100)
long_df<-read_csv("in/csvs_from_michael/mcd14ml-global-trend-by-year.csv") 

p_nf_single <-ggplot(night_fraction_trends %>% filter(p<0.05) %>% 
                       mutate(`Trend in Percent Night Active Fire Detections` = ifelse(beta > 0,
                                                         "Positive", 
                                                         "Negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25, fill="white")+
  geom_raster(data = lck_df%>%
                filter(Burnable == "yes"), 
              aes(x=x,y=y), fill = "grey90")+
  geom_raster(aes(x=x,y=y,fill=`Trend in Percent Night Active Fire Detections`)) +
  scale_fill_manual(values = posneg_cols)+
  theme_void()+
  ylim(c(-52.625, 75.125)) + 
  theme(text = element_text(family="DejaVuSans"),
        plot.caption = element_text(hjust=0.5, size=rel(1.2)),
        legend.position = "bottom",
        legend.text = element_text(size = 30),
        legend.title = element_text(size=30),
        panel.background = element_rect(color = "black", size=2))

diff_val<-filter(long_df, acq_year==2003, dn_detect=="day")%>% 
  pull(mean_frp_per_detection)-
    filter(long_df, acq_year==2003, dn_detect=="night")%>% 
      pull(mean_frp_per_detection)

dayminmax <- filter(long_df, dn_detect =="day")%>% 
  pull(mean_frp_per_detection) %>%
  range

nightminmax <- filter(long_df, dn_detect =="night")%>% 
  pull(mean_frp_per_detection) %>%
  range

p_inset2 <- long_df %>%
  mutate(adj_frp = ifelse(dn_detect == "day", 
                          mean_frp_per_detection-diff_val,
                          mean_frp_per_detection),
         dn_detect = str_to_title(dn_detect)%>% str_sub(1,1)) %>%
  ggplot(aes(x=acq_year, y=adj_frp, color = dn_detect))+
  geom_line(lwd=1) +
  scale_color_manual(values=daynight_cols)+
  theme_void() +
  scale_y_continuous(name = "", breaks = c(28.5, 40),labels = c(41,49), 
                     sec.axis = dup_axis(name = "",breaks = c(28.5,40), 
                                         labels = c(28.5,40)))+ 
  ggtitle("Fire Radiative Power\n(MW per detection)")+
  theme(axis.line=element_blank(),
        axis.text.y.left = element_text(color = daynight_cols[1], size=20),
        axis.text.y.right = element_text(color = daynight_cols[2], size=20),
        axis.text.x.bottom = element_text(size=20), 
        axis.ticks = element_blank(),
        legend.position = c(0.01,.95),
        legend.justification = c(0,1),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=20),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(color = "black", 
                                        size=1, 
                                        fill="transparent"))+
  scale_x_continuous(breaks = c(2003, 2020), 
                     labels = c("   2003", "2020   "),
                     name=""); p_inset2



ggdraw() +
  draw_plot(p_nf_single) +
  draw_plot(p_inset2, 0.02,0.2, 0.25,0.29) +
  ggsave("out/percent_night_w_inset.png", height=8, width=20)

