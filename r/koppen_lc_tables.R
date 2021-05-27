library(tidyverse)
library(stars)
library(raster)
library(s2)
library(foreach)
library(doParallel)

# koppen look up tables=========
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

# reading stuff in ===============
lck<-read_stars("in/lc_koppen_2010_mode.tif")

# calculating area ==========
lck_poly<-lck %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>%
  as_tibble() %>% 
  dplyr::select(lck = lc_koppen_2010_mode.tif, geometry)  %>% 
  mutate(geometry = as_s2_geography(geometry)) %>%
  mutate(area_km2 = s2_area(geometry)/1000000)

# summarising by landcover/koppen =======
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

# making a table with the thresholds ==========
thresholds <- read_csv("in/updated-goes-af-vpd-thresholds.csv") %>%
  dplyr::select(lc_kop = lc_name, vpd_thresh_hpa, sd) %>%
  left_join(lck_tab %>% dplyr::select(lc_kop, area_km2)) %>%
  mutate(millions_of_km2 = area_km2/1000000) %>%
  # dplyr::select(-area_km2) %>%
  arrange(lc_kop)

write_csv(thresholds %>% dplyr::select(-area_km2),"out/lck_thresh_area.csv")

#148940000 km2 (from wikipedia)
total_land_area = lck_tab %>% pull(area_km2) %>% sum
burnable_land_area = thresholds %>% pull(area_km2) %>% sum() # 90M km2
burnable_land_area/(total_land_area) # this is where I got the 62% for the paper

# getting extended table 1 stats right here baby (michael's script has a lot of 
# adjusting that is no longer necessary)

# doing the overpass adjustment for 0.25 degrees================================
dir.create("data/adjusted_counts_025", recursive=T)
dir.create("data/adjusted_frp_025", recursive=T)

# getting files
system(paste("aws s3 sync",
             "s3://earthlab-jmcglinchy/night_fire/gridded/vars_refresh_may2021/CSV_nocorn_grid_0_25_degree_vars",
             "data/gridded_mod14_025",
             "--only-show-errors"))

system(paste("aws s3 sync",
             "s3://earthlab-mkoontz/MODIS-overpass-counts_0.25_analysis-ready/month_2003-2020/",
             "data/overpass_counts_025",
             "--only-show-errors"))

# creating the lists of files ==================================================
# overpass counts
op_days_025 <- list.files("data/overpass_counts_025", 
                      full.names = TRUE, pattern = "*day*") 
op_nights_025 <- list.files("data/overpass_counts_025", 
                        full.names = TRUE, pattern = "*night*")

# day and night counts
mod14_day_counts_025<- list.files("data/gridded_mod14_025/AFC_num", 
                              full.names = TRUE,
                              pattern = "_D_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002)%>%
  separate(value, sep = "_", 
           into = c("g1", "g2","g3","g9","g4","g5","g7","month","g6"),
           remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date))

mod14_night_counts_025 <- list.files("data/gridded_mod14_025/AFC_num", 
                                   full.names = TRUE,
                                   pattern = "_N_") %>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g4","g5","g6","g7","month", "g8"), remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date)) 

# day and night total frp
mod14_day_frp_025<- list.files("data/gridded_mod14_025/FRP_total", 
                                  full.names = TRUE,
                                  pattern = "_D_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002)%>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g9","g4","g5","g7","month","g6"), remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date))

mod14_night_frp_025 <- list.files("data/gridded_mod14_025/FRP_total", 
                                     full.names = TRUE,
                                     pattern = "_N_") %>%
  as_tibble  %>%
  mutate(year = str_extract(value, "\\d{4}") %>% as.numeric) %>%
  filter(year > 2002) %>%
  separate(value, sep = "_", into = c("g1", "g2","g3","g4","g5","g6","g7","month", "g8"), remove = FALSE) %>%
  dplyr::select(-starts_with("g")) %>%
  mutate(date = as.Date(paste(year, month, "01", sep="-"), "%Y-%B-%d"),
         month_n = lubridate::month(date))

# adjusting the counts in parallel ===========
registerDoParallel(detectCores()-1)
foreach(i = 1:nrow(mod14_day_counts_025))%dopar%{
  # counts ==========
  # day counts
  system(paste("echo",i))
  month <- mod14_day_counts_025$month_n[i]
  
  counts <- raster::raster(mod14_day_counts_025$value[i])
  overpasses<- raster::raster(op_days_025[month])
  
  adjusted <- counts/overpasses
  outfile <- paste0("data/adjusted_counts_025/op-adjusted_D_", 
                    mod14_day_counts_025$year[i],
                    str_pad(month, width=2, side="left", pad = "0"),
                    ".tif"
  )
  writeRaster(adjusted,outfile, overwrite=TRUE)  
  
  # night counts
  month_n <- mod14_night_counts_025$month_n[i]
  
  counts_n <- raster::raster(mod14_night_counts_025$value[i])
  overpasses_n <- raster::raster(op_nights_025[month_n])
  
  adjusted_n <- counts_n/overpasses_n
  outfile_n <- paste0("data/adjusted_counts_025/op-adjusted_N_", 
                      mod14_night_counts_025$year[i],
                      str_pad(month, width=2, side="left", pad = "0"),
                      ".tif"
  )
  writeRaster(adjusted_n,outfile_n, overwrite=TRUE)
  
  # night fraction
  if(str_extract(outfile, "\\d{6}")==str_extract(outfile_n, "\\d{6}")){
    
    nf<- adjusted_n/(adjusted+adjusted_n)
    outfile_nf <- paste0("data/adjusted_counts_025/op-adjusted_NF_", 
                         mod14_night_counts_025$year[i],
                         str_pad(month, width=2, side="left", pad = "0"),
                         ".tif"
    )
    writeRaster(nf,outfile_nf, overwrite=TRUE)}
  rm(month, counts, overpasses, adjusted, outfile, 
     month_n, counts_n, overpasses_n, adjusted_n, outfile_n)
  ## FRP ==========
  # day frp
  month <- mod14_day_frp_025$month_n[i]
  
  counts <- raster::raster(mod14_day_frp_025$value[i])
  overpasses<- raster::raster(op_days_025[month])
  
  adjusted <- counts/overpasses
  outfile <- paste0("data/adjusted_frp_025/op-adjusted_D_", 
                    mod14_day_frp_025$year[i],
                    str_pad(month, width=2, side="left", pad = "0"),
                    ".tif"
  )
  writeRaster(adjusted,outfile, overwrite=TRUE)  
  
  # night frp
  month_n <- mod14_night_frp_025$month_n[i]
  
  counts_n <- raster::raster(mod14_night_frp_025$value[i])
  overpasses_n <- raster::raster(op_nights_025[month_n])
  
  adjusted_n <- counts_n/overpasses_n
  outfile_n <- paste0("data/adjusted_frp_025/op-adjusted_N_", 
                      mod14_night_frp_025$year[i],
                      str_pad(month, width=2, side="left", pad = "0"),
                      ".tif"
  )
  writeRaster(adjusted_n,outfile_n, overwrite=TRUE)
  
  # night fraction
  if(str_extract(outfile, "\\d{6}")==str_extract(outfile_n, "\\d{6}")){
    
    nf<- adjusted_n/(adjusted+adjusted_n)
    outfile_nf <- paste0("data/adjusted_frp_025/op-adjusted_NF_", 
                         mod14_night_frp_025$year[i],
                         str_pad(month, width=2, side="left", pad = "0"),
                         ".tif"
    )
    writeRaster(nf,outfile_nf, overwrite=TRUE)}
}

# aggregation ==============

adjusted_n <- list.files("data/adjusted_counts_025", full.names = TRUE, pattern = "_N_") %>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

adjusted_d <- list.files("data/adjusted_counts_025", full.names = TRUE, pattern = "_D_")%>%
  as_tibble() %>%
  mutate(year = str_extract(value,"\\d{4}"))

# frp <- list.files("data/gridded_mod14_2_5/FRP_mean", pattern = "*_N_*", full.names = TRUE) %>%
#   as_tibble() %>%
#   mutate(year = str_extract(value,"\\d{4}"))

dir.create("out/aggregations_2003-2020")


# dc<-adjusted_d %>%
#   pull(value) %>%
#   read_stars(along = "band") %>%
#   st_apply(MARGIN = 1:2, FUN = function(x)sum(x, na.rm = TRUE))
# 
# write_stars(dc, 
#             dsn = paste0("data/aggregations_2003-2020/day_counts.tif"))


dcr<-adjusted_d %>%
  pull(value) %>%
  raster::stack()

beginCluster()
dcrxx <- clusterR(dcr, function(x)sum(x, na.rm=T), verbose=T)
endCluster()

writeRaster(dcrxx, 
            filename = "out/aggregations_2003-2020/day_counts_rast.tif",
            overwrite=T)


ncr<-adjusted_n %>%
  pull(value) %>%
  raster::stack()

beginCluster()  
ncrxx <- clusterR(ncr, function(x)sum(x, na.rm=T), verbose=T)
endCluster()

writeRaster(ncrxx, 
            filename = "out/aggregations_2003-2020/night_counts_rast.tif", 
            overwrite=T)

ncfr<- ncrxx/(ncrxx+dcrxx)
writeRaster(ncfr, 
            filename = "out/aggregations_2003-2020/night_fraction_rast.tif", 
            overwrite=T)

# 
# nc<-adjusted_n %>%
#   pull(value) %>%
#   read_stars(along = "band") %>%
#   st_apply(MARGIN = 1:2, FUN = function(x)sum(x, na.rm = TRUE))
# write_stars(nc, 
#             dsn = paste0("data/aggregations_2003-2020/night_counts.tif"))
# ncf<- nc/(nc+dc)
# write_stars(ncf, 
#             dsn = paste0("data/aggregations_2003-2020/night_fraction.tif"))

# then match the lck layer, stack the aggregations, make a table =====

lckrastsh<-raster::shift(x = lckrast, dx = -179.75, dy = -0.25) %>%
  writeRaster("out/aggregations_2003-2020/aaalck_shifted.tif")

aggregation_tables <- list.files("out/aggregations_2003-2020", full.names = TRUE) %>%
  sort()%>%
  raster::stack() %>%
  as.data.frame() %>%
  filter(!is.na(aaalck_shifted)) %>%
  group_by(lck = aaalck_shifted) %>%
  summarise(day_counts = sum(day_counts_rast, na.rm=TRUE),
            night_counts = sum(night_counts_rast, na.rm=TRUE),
            night_fraction = mean(night_fraction_rast, na.rm=TRUE)) %>%
  ungroup%>%
  mutate(kop = lut_kop[str_sub(lck,1,1)],
         lc = lut_lc[str_sub(lck,2,3)],
         lc_kop = paste(kop, lc)) %>%
  dplyr::select(lc_kop, day_counts, night_counts, night_fraction) %>%
  na.omit()
