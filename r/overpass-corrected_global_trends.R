# modeling overpass count-adjusted active fire data

# setup ==========

library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(mblm)
library(doParallel)
library(foreach)

# install_github("leandroroser/EcoGenetics-devel")
# library(EcoGenetics)

# functions ========
parallel_theilsen <- function(stack_df, zero_to_na=FALSE, pb =TRUE,
                              minimum_sample = 10, workers = 1){
  require(doParallel)
  require(foreach)
  cells<-unique(stack_df$xy)
  
  
  registerDoParallel(workers)
  result<-foreach(i = cells, .combine = bind_rows)%dopar%{
    
    dd<-filter(stack_df, xy==i) %>%
      mutate(timestep = 1:nrow(.))
    
    if(zero_to_na) dd<-dd %>% filter(value>0)
    
    if(nrow(dd)>=minimum_sample & sum(dd$value)>0){
      mod<-mblm(value~timestep, data = dd, repeated =TRUE)
      sum<-summary(mod)
      df<-tibble(x = dd$x[1], y=dd$y[1], 
                 p = sum$coefficients[2,4],
                 beta = sum$coefficients[2,1],
                 n = nrow(dd))
      if(pb) system(paste("echo", df[1,1], df[1,2],"p=", df[1,3],"b=", df[1,4]))
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

op_days <- list.files("data/overpass_counts/month_2003-2020", 
                     full.names = TRUE, pattern = "*day*") 
op_nights <- list.files("data/overpass_counts/month_2003-2020", 
                     full.names = TRUE, pattern = "*night*")

op_days25 <- list.files("data/overpass_counts_2_5/month_2003-2020", 
                      full.names = TRUE, pattern = "*day*") 
op_nights25 <- list.files("data/overpass_counts_2_5/month_2003-2020", 
                        full.names = TRUE, pattern = "*night*")

system(paste("aws s3 sync",
             "s3://earthlab-jmcglinchy/night_fire/gridded/vars_refresh_may2021/CSV_nocorn_grid_1_0_degree_vars",
             "data/gridded_mod14"))

system(paste("aws s3 sync",
             "s3://earthlab-jmcglinchy/night_fire/gridded/vars_refresh_may2021/CSV_nocorn_grid_2_5_degree_vars",
             "data/gridded_mod14_2_5"))

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

mod14_day_counts_25<- list.files("data/gridded_mod14_2_5/AFC_num", 
                              full.names = TRUE,
                              pattern = "_D_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002)%>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g4","g5","g7","g8","g9","month","g6"), remove = FALSE) %>%
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

mod14_night_counts_25<- list.files("data/gridded_mod14_2_5/AFC_num", 
                                full.names = TRUE,
                                pattern = "_N_") %>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g4","g5","g7","g8","g9","month","g6"), remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date)) 

# do the adjusting =================
dir.create("data/adjusted_counts")
dir.create("data/adjusted_counts_25")

registerDoParallel(detectCores()-1)
foreach(i = 1:nrow(mod14_day_counts))%dopar%{
  system(paste("echo",i))
  month <- mod14_day_counts$month_n[i]
  
  counts <- raster::raster(mod14_day_counts$value[i])
  overpasses<- raster::raster(op_days[month])
  
  adjusted <- counts/overpasses
  outfile <- paste0("data/adjusted_counts/op-adjusted_D_", 
                   mod14_day_counts$year[i],
                   str_pad(month, width=2, side="left", pad = "0"),
                   ".tif"
                   )
  writeRaster(adjusted,outfile, overwrite=TRUE)  

  month_n <- mod14_night_counts$month_n[i]
  
  counts_n <- raster::raster(mod14_night_counts$value[i])
  overpasses_n <- raster::raster(op_nights[month])
  
  adjusted_n <- counts_n/overpasses_n
  outfile_n <- paste0("data/adjusted_counts/op-adjusted_N_", 
                    mod14_night_counts$year[i],
                    str_pad(month, width=2, side="left", pad = "0"),
                    ".tif"
  )
  writeRaster(adjusted_n,outfile_n, overwrite=TRUE)
  
  
  if(str_extract(outfile, "\\d{6}")==str_extract(outfile_n, "\\d{6}")){
  
  nf<- adjusted_n/(adjusted+adjusted_n)
  outfile_nf <- paste0("data/adjusted_counts/op-adjusted_NF_", 
                      mod14_night_counts$year[i],
                      str_pad(month, width=2, side="left", pad = "0"),
                      ".tif"
  )
  writeRaster(nf,outfile_nf, overwrite=TRUE)}
}

# making things annual =====================
adjusted_n <- list.files("data/adjusted_counts_25", full.names = TRUE, pattern = "*_N_*") %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

adjusted_d <- list.files("data/adjusted_counts_25", full.names = TRUE, pattern = "*_D_*")%>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

frp <- list.files("data/gridded_mod14_2_5/FRP_mean", pattern = "*_N_*", full.names = TRUE) %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

dir.create("data/annual_adjusted_counts_25")
years <- 2003:2020
for(y in years){
  print(y)
  d<-filter(adjusted_d, year == y) %>%
    pull(value) %>%
    raster::stack() %>%
    sum
  writeRaster(d, 
              filename = paste0("data/annual_adjusted_counts_25/day_adj_annual_counts_",y, ".tif"),
              overwrite = TRUE)
  n<-filter(adjusted_n, year == y) %>%
    pull(value) %>%
    raster::stack() %>%
    sum
  writeRaster(n, 
              filename = paste0("data/annual_adjusted_counts_25/night_adj_annual_counts_",y, ".tif"),
              overwrite = TRUE)
  nf<- n/(n+d)
  writeRaster(nf, 
              filename = paste0("data/annual_adjusted_counts_25/night_fraction_adj_annual_counts_",y, ".tif"),
              overwrite = TRUE)
  
  ff <- filter(frp, year == y) %>%
    pull(value) %>%
    raster::stack() %>%
    mean
  writeRaster(ff, 
              filename = paste0("data/annual_adjusted_counts_25/night_frp_annual_mean_",y, ".tif"),
              overwrite = TRUE)
  
}


# time series analysis=========

# annual =================

# day counts===========
day_counts <- list.files("data/annual_adjusted_counts", pattern = "day*", full.names = TRUE) %>%
  raster::stack()
day_counts_25 <- list.files("data/annual_adjusted_counts_25", pattern = "day*", full.names = TRUE) %>%
  raster::stack()

# removing cells with all 0s
sum_day_counts <- raster::calc(day_counts,function(x)sum(x))
sum_day_counts_25 <- raster::calc(day_counts_25,function(x)sum(x))

day_counts[sum_day_counts == 0] <- NA
day_counts_25[sum_day_counts_25 == 0] <- NA

# day_trends<- eco.theilsen(day_counts, dates = 2003:2020)
dir.create("out")
# system("mv slope.tif out/ts_estimate_daycount_annual.tif")
# system("mv pvalue.tif out/ts_p_daycount_annual.tif")

# 1 degree
day_counts_df<-day_counts %>%
  as.data.frame(xy=TRUE) %>%
  na.omit %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4", "year"), 
           remove = TRUE, sep = "_") %>%
  mutate(xy = str_c(x,y)) %>%
  arrange(year)

day_counts_trends<-parallel_theilsen(day_counts_df,
                                     zero_to_na = FALSE, 
                                     workers=4,
                                     minimum_sample = 10)

p_dc<-ggplot(day_counts_trends %>% filter(p<0.05->alpha) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in day counts, 1 degree, 2003-2018, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/day_count_trend_1_deg_2003-2018.png")

# 2.5 degree 

day_counts_df_25<-day_counts_25 %>%
  as.data.frame(xy=TRUE) %>%
  na.omit %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4", "year"), 
           remove = TRUE, sep = "_") %>%
  mutate(xy = str_c(x,y)) %>%
  arrange(year)

day_counts_trends_25<-parallel_theilsen(day_counts_df_25,workers=4,
                                        zero_to_na = FALSE,
                                        minimum_sample = 10)

ggplot(day_counts_trends_25 %>% filter(p<0.01->alpha) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in day counts, 2.5 degrees, 2003-2018, p<0.01"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/day_count_trend_2_5_deg_2003-2018.png")

# night counts ===========================
night_counts <- list.files("data/annual_adjusted_counts", pattern = "night_adj*", full.names = TRUE) %>%
  raster::stack()
night_counts_25 <- list.files("data/annual_adjusted_counts_25", pattern = "night_adj*", full.names = TRUE) %>%
  raster::stack()

sum_night_counts <- night_counts %>% raster::calc(function(x)sum(x))
sum_night_counts_25 <- night_counts_25 %>% raster::calc(function(x)sum(x))

night_counts[sum_night_counts == 0] <- NA
night_counts_25[sum_night_counts_25 == 0] <- NA


# night_trends<- eco.theilsen(night_counts, dates = 2003:2020)

# system("mv slope.tif out/ts_estimate_nightcount_annual.tif")
# system("mv pvalue.tif out/ts_p_nightcount_annual.tif")


# 1 degree
night_counts_df<-night_counts %>%
  as.data.frame(xy=TRUE) %>%
  na.omit %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4", "year"), 
           remove = TRUE, sep = "_") %>%
  mutate(xy = str_c(x,y)) %>%
  arrange(year)

night_counts_trends<-parallel_theilsen(night_counts_df,
                                     zero_to_na = FALSE, 
                                     workers=4,
                                     minimum_sample = 10)

p_nc<-ggplot(night_counts_trends %>% filter(p<0.05->alpha) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in night counts, 1 degree, 2003-2018, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_count_trend_1_deg_2003-2018.png")

# 2.5 degree 

night_counts_df_25<-night_counts_25 %>%
  as.data.frame(xy=TRUE) %>%
  na.omit %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4", "year"), 
           remove = TRUE, sep = "_") %>%
  mutate(xy = str_c(x,y)) %>%
  arrange(year)

night_counts_trends_25<-parallel_theilsen(night_counts_df_25,workers=4,
                                        zero_to_na = FALSE,
                                        minimum_sample = 10)

ggplot(night_counts_trends_25 %>% filter(p<0.05->alpha) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in night counts, 2.5 degrees, 2003-2018, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_count_trend_2_5_deg_2003-2018.png")

# night fraction ===========
night_fractions <- list.files("data/annual_adjusted_counts", pattern = "night_fraction*", full.names = TRUE) %>%
  raster::stack()
night_fractions_25 <- list.files("data/annual_adjusted_counts_25", pattern = "night_fraction*", full.names = TRUE) %>%
  raster::stack()

night_fractions[sum_night_counts+sum_day_counts == 0] <- NA
night_fractions_25[sum_night_counts_25+sum_day_counts_25 == 0] <- NA

# night_fraction_trends<- eco.theilsen(night_fractions, dates = 1:nlayers(night_fractions))
# 
# system("mv slope.tif out/ts_estimate_nightfraction_annual.tif")
# system("mv pvalue.tif out/ts_p_nightfraction_annual.tif")

# 1 degree
night_fraction_df<-night_fractions %>%
  as.data.frame(xy=TRUE) %>%
  na.omit %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4","x5", "year"), 
           remove = TRUE, sep = "_") %>%
  mutate(xy = str_c(x,y)) %>%
  arrange(year)

night_fraction_trends<-parallel_theilsen(night_fraction_df,
                                       zero_to_na = FALSE, 
                                       workers=4,
                                       minimum_sample = 10)

p_nf<-ggplot(night_fraction_trends %>% filter(p<0.05->alpha) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in night fractions, 1 degree, 2003-2018, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_fraction_trend_1_deg_2003-2018.png")

# 2.5 degree 

night_fraction_df_25<-night_fractions_25 %>%
  as.data.frame(xy=TRUE) %>%
  na.omit %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4","x5", "year"), 
           remove = TRUE, sep = "_") %>%
  mutate(xy = str_c(x,y)) %>%
  arrange(year)

night_fraction_trends_25<-parallel_theilsen(night_fraction_df_25,workers=4,
                                          zero_to_na = FALSE,
                                          minimum_sample = 10)

ggplot(night_fraction_trends_25 %>% filter(p<0.05->alpha) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in night fractions, 2.5 degrees, 2003-2018, p < 0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_fraction_trend_2_5_deg_2003-2018.png")

# night frp ===============

night_frp <- list.files("data/annual_adjusted_counts", pattern = "night_frp_annual*", full.names = TRUE) %>%
  raster::stack()

night_frp[night_frp == 0] <- NA # 0 should be NA

night_frp_df<-night_frp %>%
  as.data.frame(xy=TRUE) %>%
  na.omit %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4", "year"), 
           remove = FALSE, sep = "_") %>%
  mutate(xy = str_c(x,y))

night_frp_trends<-parallel_theilsen(night_frp_df,zero_to_na = TRUE,workers=4, minimum_sample = 10)

p_frp <- ggplot(night_frp_trends %>% filter(p<0.05->alpha) %>% 
         mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in night FRP, 1 degree, 2003-2018, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_frp_trend_1_deg_2003-2018.png")



# all together
library(ggpubr)
ggarrange(p_dc, p_nc, p_nf, p_frp, nrow = 2, ncol=2) +
  ggsave(height = 7, width=12, filename = "out/annual_4pan.png")

# monthly ============

workers<- detectCores()/2
# day counts ===========
dc_m <- list.files("data/adjusted_counts",
                          pattern = "*_D_*", full.names = TRUE)%>%
  as_tibble  %>%
  mutate(ym = str_extract(value, "\\d{6}") %>% as.numeric,
         year = str_sub(ym, 1,4),
         month = str_sub(ym, 5,6)) %>%
  filter(year > 2002) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%m-%d"),
         month_n = lubridate::month(date)) %>%
  dplyr::arrange(date)

dc_stack_m<-raster::stack(pull(dc_m, value))
sum_dc_m <- dc_stack_m %>% raster::calc(function(x)sum(x))
dc_stack_m[sum_dc_m== 0] <- NA # 5 minutes

dc_df_m <- dc_stack_m %>%
  as.data.frame(xy=TRUE) %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(ym = str_extract(value, "\\d{6}") %>% as.numeric,
         year = str_sub(ym, 1,4),
         month = str_sub(ym, 5,6)) %>%
  filter(year > 2002) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%m-%d"),
         month_n = lubridate::month(date)) %>%
  dplyr::arrange(date) %>%
  mutate(xy = str_c(x,y))

dc_trends_m<-parallel_theilsen(dc_df_m,
                                      zero_to_na = TRUE,
                                      workers=workers, 
                                      minimum_sample = 30)

p_dc_m <- ggplot(dc_trends_m %>% 
                  filter(p<0.05->alpha) %>% 
                  mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Monthly day detections, 1 degree, 2003-2018, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/monthly_day_count_trend_1_deg_2003-2018.png")

# night frp ===========
# need to get the files in the correct order
night_frp_m <- list.files("data/gridded_mod14/FRP_mean",
                          pattern = "*_N_*", full.names = TRUE)%>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  separate(value, sep = "_",
           into = c("g1", "g2","g3","g4","g5","g7","month","g6"), 
           remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date)) %>%
  dplyr::arrange(date)

night_frp_stack_m<-raster::stack(pull(night_frp_m, value))

# takes too long, just use the one from the annual calculations
# sum_night_frp_m <- night_frp_m %>% raster::calc(sum)

# night_frp_stack_m[night_frp_stack_m== 0] <- NA # 5 minutes

night_frp_df_m <- night_frp_stack_m %>%
  as.data.frame(xy=TRUE) %>%
  pivot_longer(cols = names(.)[3:ncol(.)],names_to = "filename", values_to = "value") %>%
  separate(filename, into = c("x1","x2", "x3","x4","month", "year"), 
           remove = FALSE, sep = "_") %>%
  mutate(xy = str_c(x,y)) %>%
  filter(value>0)

night_frp_trends_m<-parallel_theilsen(night_frp_df_m,
                                      zero_to_na = TRUE,
                                      workers=workers, 
                                      minimum_sample = 30)

p_frp <- ggplot(night_frp_trends_m %>% 
                  filter(p<0.05->alpha) %>% 
                  mutate(Trend = ifelse(beta > 0, "positive", "negative"))) +
  geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
  geom_raster(aes(x=x,y=y,fill=Trend)) +
  scale_fill_manual(values = c("blue","red"))+
  theme_void()+
  ggtitle(paste("Trends in night FRP, 1 degree, 2003-2018, p<0.05"))+
  theme(legend.position = c(0.1,0.2),
        legend.justification = c(0,0)) +
  ggsave("out/night_frp_trend_1_deg_2003-2018.png")

# plots ========
plot_sig_ts <- function(nf, pv, title, binary=TRUE){

  nf<- raster(nf)
  pv<- raster(pv)
  
  x<-pv
  x[pv<0.05] <- 1
  x[pv>0.05] <- 0
  
  z<-nf*x
  
  z[z==0]<-NA
  if(binary ==T){
  z[z>0]<- 1
  z[z<0]<- -1}
  
  
  z %>%
    as.data.frame(xy=TRUE) %>%
    ggplot()+
    geom_sf(data = st_read("world.gpkg"), lwd=0.25)+
    geom_raster(aes(x=x,y=y,fill=layer))+
    scale_fill_gradient2(na.value = "transparent", low = "blue", high = "red") +
    theme_void()+
    theme(legend.position = c(0.1,0.2),
          legend.justification = c(0,0))+
    ggtitle(title)
}

plot_sig_ts(nf="slope.tif",
            pv="pvalue.tif",
            title = "night frp",
            binary = TRUE)

plot_sig_ts(nf="out/ts_estimate_nightfraction_annual.tif",
            pv="out/ts_p_nightfraction_annual.tif",
            title = "Night Fraction")

plot_sig_ts(nf="out/ts_estimate_nightcount_annual.tif",
            pv="out/ts_p_nightcount_annual.tif",
            title = "Night Count")

plot_sig_ts(nf="out/ts_estimate_daycount_annual.tif",
            pv="out/ts_p_daycount_annual.tif",
            title = "Day Count")
