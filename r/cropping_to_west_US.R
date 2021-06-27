# getting data for john for western US
library(tidyverse)
library(terra)

westUSext <- terra::ext(c(-125,-103,31,49))

rasts<-list.files("data/adjusted_counts", pattern = "N_", full.names = TRUE) %>%
  terra::rast() %>% 
  terra::crop(westUSext)

terra::writeRaster(rasts, filename = "west_us_night_aft_per_overpass_monthly_1deg_2003-2020.tif")


mean_FRP <- list.files("data/gridded_mod14/FRP_mean", full.names=TRUE,
                       pattern = "N_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002)%>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g4","g5","g7","month","g6"), remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date)) %>%
  arrange(date) %>%
  pull(value) %>%
  terra::rast() %>% 
  terra::crop(westUSext)

terra::writeRaster(mean_FRP, filename = "west_us_night_frp_per_afd_1deg_2003-2020.tif", overwrite=TRUE)
