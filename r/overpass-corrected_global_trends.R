# modeling overpass counts
library(devtools)

# install_github("leandroroser/EcoGenetics-devel")

library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(EcoGenetics)
# file management

system(paste("aws s3 sync",
             "s3://earthlab-mkoontz/MODIS-overpass-counts_1_analysis-ready",
             "data/overpass_counts"))

op_days <- list.files("data/overpass_counts/month_2003-2020", 
                     full.names = TRUE, pattern = "*day*") 
op_nights <- list.files("data/overpass_counts/month_2003-2020", 
                     full.names = TRUE, pattern = "*night*")

system(paste("aws s3 sync",
             "s3://earthlab-jmcglinchy/night_fire/gridded/vars_refresh_may2021/CSV_nocorn_grid_1_0_degree_vars",
             "data/gridded_mod14"))

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
for(i in 1:nrow(mod14_day_counts)){
  print(i)
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
}

for(i in 1:nrow(mod14_night_counts)){
  print(i)
  month <- mod14_night_counts$month_n[i]
  
  counts <- raster::raster(mod14_night_counts$value[i])
  overpasses<- raster::raster(op_nights[month])
  
  adjusted <- counts/overpasses
  outfile <- paste0("data/adjusted_counts/op-adjusted_N_", 
                    mod14_night_counts$year[i],
                    str_pad(month, width=2, side="left", pad = "0"),
                    ".tif"
  )
  writeRaster(adjusted,outfile, overwrite=TRUE)
}

# making things annual
adjusted_n <- list.files("data/adjusted_counts", full.names = TRUE, pattern = "*_N_*") %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

adjusted_d <- list.files("data/adjusted_counts", full.names = TRUE, pattern = "*_D_*")%>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))
dir.create("data/annual_adjusted_counts")
years <- 2003:2020
for(y in years){
  print(y)
  d<-filter(adjusted_d, year == y) %>%
    pull(value) %>%
    raster::stack() %>%
    sum
  writeRaster(d, filename = paste0("data/annual_adjusted_counts/day_adj_annual_counts_",y, ".tif"))
  n<-filter(adjusted_n, year == y) %>%
    pull(value) %>%
    raster::stack() %>%
    sum
  writeRaster(n, filename = paste0("data/annual_adjusted_counts/night_adj_annual_counts_",y, ".tif"))
  nf<- n/(n+d)
  writeRaster(nf, filename = paste0("data/annual_adjusted_counts/night_fraction_adj_annual_counts_",y, ".tif"))
  
}


# time series analysis
# day counts
day_counts <- list.files("data/annual_adjusted_counts", pattern = "day*", full.names = TRUE) %>%
  raster::stack()

sum_day_counts <- day_counts %>% raster::calc(sum)

day_counts[sum_day_counts == 0] <- NA

day_trends<- eco.theilsen(day_counts, dates = 2003:2020)
dir.create("out")
system("mv slope.tif out/ts_estimate_daycount_annual.tif")
system("mv pvalue.tif out/ts_p_daycount_annual.tif")

# night counts
night_counts <- list.files("data/annual_adjusted_counts", pattern = "night_adj*", full.names = TRUE) %>%
  raster::stack()

sum_night_counts <- night_counts %>% raster::calc(sum)

night_counts[sum_night_counts == 0] <- NA

night_trends<- eco.theilsen(night_counts, dates = 2003:2020)

system("mv slope.tif out/ts_estimate_nightcount_annual.tif")
system("mv pvalue.tif out/ts_p_nightcount_annual.tif")

# night fraction
night_fractions <- list.files("data/annual_adjusted_counts", pattern = "night_fraction*", full.names = TRUE) %>%
  raster::stack()

night_fractions[sum_night_counts+sum_day_counts == 0] <- NA

night_fraction_trends<- eco.theilsen(night_fractions, dates = 1:nlayers(night_fractions))

system("mv slope.tif out/ts_estimate_nightfraction_annual.tif")
system("mv pvalue.tif out/ts_p_nightfraction_annual.tif")

# night frp



plot_sig_ts <- function(nf, pv, title){

  nf<- raster(nf)
  pv<- raster(pv)
  
  x<-pv
  x[pv<0.05] <- 1
  x[pv>0.05] <- 0
  
  z<-nf*x
  
  z[z==0]<-NA
  z[z>0]<- 1
  z[z<0]<- -1
  plot(z)
  
  z %>%
    as.data.frame(xy=TRUE) %>%
    ggplot()+
    geom_raster(aes(x=x,y=y,fill=layer))+
    scale_fill_gradient2(na.value = "transparent", low = "blue", high = "red") +
    theme_minimal() +
    ggtitle(title)
}
nf <- raster("out/ts_estimate_nightfraction_annual.tif")
pv <- raster("out/ts_p_nightfraction_annual.tif")

plot_sig_ts(nf="out/ts_estimate_nightfraction_annual.tif",
            pv="out/ts_p_nightfraction_annual.tif",
            title = "Night Fraction")

plot_sig_ts(nf="out/ts_estimate_nightcount_annual.tif",
            pv="out/ts_p_nightcount_annual.tif",
            title = "Night Count")

plot_sig_ts(nf="out/ts_estimate_daycount_annual.tif",
            pv="out/ts_p_daycount_annual.tif",
            title = "Day Count")
