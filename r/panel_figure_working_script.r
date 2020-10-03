# setup ========================================================================
libs<- c("tidyverse", "sf", "raster", "gganimate", "gifski", "stars","viridis",
         "gdalUtils")
lapply(libs, library, character.only = TRUE)

s3path <- "s3://earthlab-natem/modis-burned-area/MCD64A1/C6/tif_converted_alltiles"
local_data <- "/home/a/data"

mylabels = label_bquote(
  cols = .(daynight) ~ .(event_day)
)
options(stringsAsFactors = F)
prism_crs <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
meat <- "PRISM_tmin_stable_4kmD1_"
bits <- "_bil.bil"


# read in AF data for CUS

viirs_af <- st_read(file.path(local_data,"fire", "modis_viirs_active_fire", 
                              "CUS_V1","fire_archive_V1_54276.shp"))
modis_af <- st_read(file.path(local_data,"fire", "modis_viirs_active_fire", 
                              "CUS_M6","fire_archive_M6_54275.shp"))


classification_table <- data.frame(
  name = c("Evergreen Needleleaf Forests", "Evergreen Broadleaf Forests",
           "Deciduous Needleleaf Forests","Deciduous Broadleaf Forests",
           "Mixed Forests","Closed Shrublands",
           "Open Shrublands","Woody Savannas",
           "Savannas","Grasslands",
           "Permanent Wetlands","Croplands",
           "Urban and Built-up Lands","Cropland/Natural  Vegetation  Mosaics",
           "Permanent Snow and Ice","Barren",
           "Water Bodies", "Unclassified"),
  value = c(1:17,255)
)

lut_lc <-  c("1"="Evergreen Needleleaf Forests", "2"= "Evergreen Broadleaf Forests",
          "3"= "Deciduous Needleleaf Forests","4"= "Deciduous Broadleaf Forests",
          "5"= "Mixed Forests","6"= "Closed Shrublands",
          "7"= "Open Shrublands","8"= "Woody Savannas",
          "9"=  "Savannas","10"= "Grasslands",
          "11"= "Permanent Wetlands","12"= "Croplands",
          "13"=  "Urban and Built-up Lands","14"= "Cropland/Natural  Vegetation  Mosaics",
          "15"=  "Permanent Snow and Ice","16"= "Barren",
          "17"=  "Water Bodies", "255"= "Unclassified")

# hot pot fire (july 2016) ==========================================================
hpextent <- c(-10874329, -10842484,3815485, 3832978) %>%
  extent

hpfile <- "bd_numeric_2016_183_h08v05.tif"

if(!file.exists(file.path("data/hotpot",hpfile))){
  system(paste("aws s3 cp",
               file.path(s3path, hpfile),
               file.path("data/hotpot",hpfile)))
}

hpm <- read_stars(file.path("data/hotpot",hpfile)) 

hp_bbox<-st_bbox(c(xmin=-10874329, xmax=-10842484,ymin=3815485, ymax=3832978), 
                 crs = hpm) %>%
  st_as_sfc()

hp_prsm <- hp_bbox %>%
  st_transform(prism_crs)%>%
  st_bbox()%>%
  st_as_sfc()

hpmc <- st_crop(hpm, hp_bbox)
hpmc[hpmc<1]<-NA
names(hpmc) <- "burn_date"

# getting land cover from the year before
hp_lc_file <- "/home/a/projects/tmin/data/lc/lc_h08v05_2015.tif"
sds <- '/home/a/data/MCD12Q1/MCD12Q1.A2012001.h08v05.006.2018054212715.hdf' %>%
  gdalUtils::get_subdatasets()
gdalUtils::gdal_translate(sds[1], dst_dataset = hp_lc_file)

hp_lc<- read_stars(hp_lc_file) %>%
  st_crop(hp_bbox)
names(hp_lc) <- "lc"
hp_lc[is.na(hpmc)==TRUE] <- NA
hp_lc<- st_as_sf(hp_lc, as_points = TRUE, merge = FALSE) %>%
  mutate(classes = lut_lc[lc],
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL)



hp_tmin_path <- "/home/a/data/PRISM/daily/tmin/y2016m7"
hp_tmin_files <- list.files(hp_tmin_path, full.names = TRUE, pattern = "\\.bil$")

if(!file.exists("data/hotpot/hp_afv.gpkg")){
  hp_afm <- modis_af %>%
    filter(substr(ACQ_DATE,1,4)%>% as.numeric > 2015 &
             substr(ACQ_DATE,1,4)%>% as.numeric < 2017)  %>%
    st_transform(crs=st_crs(hpm))%>%
    st_intersection(hp_bbox) %>%
    filter(ACQ_DATE<as.Date("2016-08-01"), ACQ_DATE>as.Date("2016-07-22")) %>%
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
           daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
           daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
           dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
           daynight = as.character(daynight)) %>%
    dplyr::select(-DAYNIGHT)%>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
  group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  hp_day0 <- min(hp_afm$ACQ_DATE -1)
  
  ll <- list()
  for(d in 1:length(unique(hp_afm$event_day))){
    ed <- as.character(hp_day0+unique(hp_afm$event_day)[d])
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(hp_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(hp_prsm, "Spatial")) %>%
      projectRaster(as(hpmc, "Raster")) 
    
    ll[[d]] <- filter(hp_afm, event_day == unique(hp_afm$event_day)[d])%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  
  hp_afm<- do.call("rbind",ll)
  
  hp_afv <- viirs_af %>%
    filter(substr(ACQ_DATE,1,4)%>% as.numeric > 2015 &
             substr(ACQ_DATE,1,4)%>% as.numeric < 2017)  %>%
    st_transform(crs=st_crs(hpm))%>%
    st_intersection(hp_bbox)%>%
    filter(ACQ_DATE<as.Date("2016-08-01"), ACQ_DATE>as.Date("2016-07-22")) %>%
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
           daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
           daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
           dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
           daynight = as.character(daynight))%>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
    group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  hp_day0 <- min(hp_afv$ACQ_DATE -1)
  
  ll <- list()
  for(d in unique(hp_afv$event_day)){
    ed <- as.character(hp_day0+d)
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(hp_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(hp_prsm, "Spatial")) %>%
      projectRaster(as(hpmc, "Raster")) 
    
    ll[[d]] <- filter(hp_afv, event_day == d)%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  
  hp_afv<- do.call("rbind",ll)
  
  st_write(hp_afv, "data/hotpot/hp_afv.gpkg", delete_dsn = TRUE)
  st_write(hp_afm, "data/hotpot/hp_afm.gpkg", delete_dsn = TRUE)
}else{
  hp_afv <- st_read("data/hotpot/hp_afv.gpkg") 
  hp_afm <- st_read("data/hotpot/hp_afm.gpkg")
  }

ggplot() +
  #geom_stars(data=hpmc, aes(x=x,y=y,fill = as.Date(burn_date, origin = "1970-01-01")))+
  geom_raster(data=hp_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = c("darkgreen","turquoise4", "lightyellow", "olivedrab",
                               "lightgreen", "grey70", "lightblue"),
                    name = "Land Cover")+
  #scale_fill_brewer(palette = "Pastel2")+
  theme_void()+
  theme(#legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
        ) +
  geom_sf(data=hp_afv %>% filter(ACQ_DATE > as.Date("2016-07-22"),
                                 ACQ_DATE < as.Date("2016-07-27")),
          aes(color = tmin), shape = 1, alpha = 0.95) +
  coord_sf(datum=NA)+
  geom_text(data = hp_afv %>% filter(ACQ_DATE > as.Date("2016-07-22"),
                                     ACQ_DATE < as.Date("2016-07-27")),
            x=-10852000, y=3816500, aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  scale_color_viridis(option="B", direction = -1, "Minimum Temperature",end = .85, begin=0  )+
  facet_wrap(~event_day + daynight, ncol=2, labeller = mylabels) +
  ggtitle("Hot Pot Fire: VIIRS Active Fire") +
  ggsave("images/tmin_hotpot_viirs_lc.png")


ggplot() +
  # geom_stars(data=hpmc,
  #            aes(x=x,y=y, fill = as.Date(burn_date, origin = "1970-01-01")))+
 #scale_fill_continuous(na.value = "white")+
  geom_raster(data=hp_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = c("darkgreen","turquoise4", "lightyellow", "olivedrab",
                               "lightgreen", "grey70", "lightblue"),
                    name = "Land Cover")+
  
  theme_void()+
  theme(#legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  ) +
  geom_sf(data=hp_afm, aes(color = tmin), shape = 19, alpha = 0.95) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = -1, "Minimum Temperature",end = .85, begin=0  )+
  #scale_colour_gradient(low = "white", high = "red")+
  geom_text(data = hp_afm,
            x=-10852000, y=3816500, 
            aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  facet_wrap(~event_day+daynight, ncol=2, labeller = mylabels) +
  ggtitle("Hot Pot Fire: MODIS Active Fire") +
  ggsave("images/tmin_hotpot_modis_lc.png")

p1 <- ggplot() +
  geom_raster(data = hpmdf, aes(x=x,y=y,fill = burn_doy))+
  ggtitle("Modis Burned Area")
p2 <-ggplot()+
  geom_sf(data = hp_afv, alpha=0.5, aes(color = daynight)) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Viirs Active Fire")
p3 <- ggplot() +
  geom_sf(data = hp_afm, alpha=0.5, aes(color = daynight)) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Modis Active Fire")
ggarrange(p1,p2,p3)

both_afs <- rbind(dplyr::select(hp_afv, LONGITUDE, LATITUDE, daynight, dt),
                  dplyr::select(hp_afm, LONGITUDE, LATITUDE, daynight, dt))

anim<-ggplot(both_afs, aes(x=LONGITUDE,y=LATITUDE))+ 
  geom_point(aes(color = daynight)) + # the business
  theme_void() + # this is a theme without any annoying lines  
  scale_color_manual(name = "", values=c("orange","darkblue")) + # setting the color scheme and naming the legend
  labs(title = 'Day: {frame_time}') + # this is setting the label on top
  transition_time(dt) +
  shadow_mark()

aa<-gganimate::animate(anim, fps=2, nframes = length(unique(both_afs$dt)));aa
anim_save(aa, filename="data/hotpot/af_both.gif")

ggplot()+
  geom_density(data = hp_afm, aes(x=FRP)) +
  geom_density(data = hp_afv, aes(x=FRP), col = "red")

# martin fire (2018) ===========================================================

martin_p <- st_read("data/martin/Martin_perimeter.shp")
martin_extent <- martin_p %>%
  st_transform(crs=st_crs(hpmc))%>%
  st_buffer(dist = 1000) %>%
  st_bbox() %>%
  st_as_sfc()

martin_prsm <- martin_extent %>%
  st_transform(prism_crs)%>%
  st_bbox()%>%
  st_as_sfc()

martin_mcdf <- "bd_numeric_2018_182_h09v04.tif" 

if(!file.exists(file.path("data/martin",martin_mcdf))){
   system(paste("aws s3 cp",
             file.path(s3path, martin_mcd),
             file.path("data/martin",martin_mcd)))
}
martin_mcd <-read_stars(file.path("data/martin",martin_mcdf)) %>%
  st_crop(martin_extent)
martin_mcd[martin_mcd<1]<-NA
names(martin_mcd) <- "burn_date"

martin_lc_file <- "/home/a/projects/tmin/data/lc/lc_h09v04_2016.tif"
sds <- '/home/a/data/MCD12Q1/MCD12Q1.A2016001.h09v04.006.2018055070404.hdf' %>%
  gdalUtils::get_subdatasets()
gdalUtils::gdal_translate(sds[1], dst_dataset = martin_lc_file)

martin_lc<- read_stars(martin_lc_file) %>%
  st_crop(martin_extent)
names(martin_lc) <- "lc"
martin_lc[is.na(martin_mcd)==TRUE] <- NA
martin_lc<- st_as_sf(martin_lc, as_points = TRUE, merge = FALSE) %>%
  mutate(classes = lut_lc[lc],
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL)

martin_tmin_path <-"/home/a/data/PRISM/daily/tmin/y2018m7"
meat = "PRISM_tmin_stable_4kmD1_"
if(!file.exists("data/martin/martin_afv.gpkg")){
  martin_afm <- modis_af %>%
    filter(substr(ACQ_DATE,1,4)%>% as.numeric > 2017 &
             substr(ACQ_DATE,1,4)%>% as.numeric < 2019)  %>%
    st_transform(crs=st_crs(martin_extent))%>%
    st_intersection(martin_extent) %>%
    filter(ACQ_DATE>as.Date("2018-07-01"), ACQ_DATE<as.Date("2018-08-01")) %>%
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
           daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
           daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
           dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
           daynight = as.character(daynight)) %>%
    dplyr::select(-DAYNIGHT)%>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
    group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  day0 <- min(martin_afm$ACQ_DATE -1)
  
  ll <- list()
  for(d in 1:length(unique(martin_afm$event_day))){
    ed <- as.character(day0+unique(martin_afm$event_day)[d])
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(martin_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(martin_prsm, "Spatial")) %>%
      projectRaster(as(martin_mcd, "Raster")) 
    
    ll[[d]] <- filter(martin_afm, event_day == unique(martin_afm$event_day)[d])%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  martin_afm <- do.call("rbind",ll)
  
  martin_afv <- viirs_af %>%
    filter(substr(ACQ_DATE,1,4)%>% as.numeric > 2017 &
             substr(ACQ_DATE,1,4)%>% as.numeric < 2019)  %>%
    st_transform(crs=st_crs(martin_extent))%>%
    st_intersection(martin_extent)%>%
    filter(ACQ_DATE>as.Date("2018-07-03"), ACQ_DATE<as.Date("2018-08-01")) %>%
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
           daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
           daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
           dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
           daynight = as.character(daynight)) %>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
    group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  day0 <- min(martin_afv$ACQ_DATE -1)
  
  ll <- list()
  for(d in 1:length(unique(martin_afv$event_day))){
    ed <- as.character(day0+unique(martin_afv$event_day)[d])
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(martin_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(martin_prsm, "Spatial")) %>%
      projectRaster(as(martin_mcd, "Raster")) 
    
    ll[[d]] <- filter(martin_afv, event_day == unique(martin_afv$event_day)[d])%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  martin_afv <- do.call("rbind",ll)
  
  
  st_write(martin_afv, "data/martin/martin_afv.gpkg", delete_dsn = TRUE)
  st_write(martin_afm, "data/martin/martin_afm.gpkg", delete_dsn = TRUE)
}else{
  martin_afv <- st_read("data/martin/martin_afv.gpkg") 
    
  martin_afm <- st_read("data/martin/martin_afm.gpkg")
}

ggplot() +
  geom_raster(data=martin_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = "lightgreen",
                    name = "Land Cover")+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()
  ) +
  geom_sf(data=martin_afm %>% filter(ACQ_DATE < as.Date("2018-07-09")),
          aes(color = tmin), shape = 19, alpha = 0.95) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = -1 
                      ,end = .85, begin=0, name = "Minimum Temperature")+
  geom_text(data = martin_afm%>% filter(ACQ_DATE < as.Date("2018-07-09")),
           x=-9698000, y=4610000, aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  facet_wrap(~event_day+daynight,ncol=2, labeller = mylabels) +
  ggtitle("Martin Fire: MODIS Active Fire") +
  ggsave("images/tmin_martin_modis_lc.png")


ggplot() +
  geom_raster(data=martin_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values =  "lightgreen",
                    name = "Land Cover")+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_sf(data=martin_afv %>% filter(ACQ_DATE < as.Date("2018-07-09")),
          aes(color = tmin), shape = 19, alpha = 0.65) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = -1, "Minimum Temperature",end = .85, begin=0  )+
  geom_text(data = martin_afv %>% filter(ACQ_DATE < as.Date("2018-07-09")),
            x=-9698000, y=4610000, aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  facet_wrap(~event_day+daynight,ncol=2, labeller = mylabels) +
  ggtitle("Martin Fire: VIIRS Active Fire")+
  ggsave("images/tmin_martin_viirs_lc.png")


p1 <- ggplot() +
  geom_raster(data = hpmdf, aes(x=x,y=y,fill = burn_doy))+
  ggtitle("Modis Burned Area")
p2 <-ggplot()+
  geom_sf(data = hp_afv, alpha=0.5, aes(color = daynight)) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Viirs Active Fire")
p3 <- ggplot() +
  geom_sf(data = hp_afm, alpha=0.5, aes(color = daynight)) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Modis Active Fire")
ggarrange(p1,p2,p3)

both_afs <- rbind(dplyr::select(martin_afv, LONGITUDE, LATITUDE, daynight, dt),
                  dplyr::select(martin_afm, LONGITUDE, LATITUDE, daynight, dt))

anim<-ggplot(both_afs, aes(x=LONGITUDE,y=LATITUDE))+ 
  geom_point(aes(color = daynight)) + # the business
  theme_void() + # this is a theme without any annoying lines  
  scale_color_manual(name = "", values=c("orange","darkblue")) + # setting the color scheme and naming the legend
  labs(title = 'Day: {frame_time}') + # this is setting the label on top
  transition_time(dt) +
  shadow_mark()

aa<-gganimate::animate(anim, fps=2, nframes = length(unique(both_afs$dt)));aa
anim_save(aa, filename="data/martin/af_both.gif")

# rim fire (2013) ==============================================================

rim_mcdf <- "bd_numeric_2018_213_h08v05.tif" 
if(!file.exists(file.path("data/rim",rim_mcdf))){
  system(paste("aws s3 cp",
               file.path(s3path, rim_mcdf),
               file.path("data/rim",rim_mcdf)))
}
rim_mcd <-read_stars(file.path("data/rim",rim_mcdf))
rim_p <- st_read(file.path(local_data, "fire", "mtbs", "mtbs_perims_DD.shp"))
  filter(Year == 2013, Acres >200000) 
rim_extent <- rim_p %>%
  st_transform(crs = st_crs(rim_mcd))
  st_buffer(dist = 10000) %>%
  st_bbox() %>%
  st_as_sfc()

rim_prsm <- rim_extent %>%
  st_transform(prism_crs)%>%
  st_bbox()%>%
  st_as_sfc()
  
rim_mcdc <- rim_mcd %>%
  st_crop(rim_extent)
# rim_mcdc[rim_mcdc<1]<-NA
names(rim_mcdc) <- "burn_date"

rim_lc_file <- "/home/a/projects/tmin/data/lc/lc_h08v05_2012.tif"
sds <- '/home/a/data/MCD12Q1/MCD12Q1.A2012001.h08v05.006.2018054212715.hdf' %>%
  gdalUtils::get_subdatasets()
gdalUtils::gdal_translate(sds[1], dst_dataset =rim_lc_file)

rim_lc<- read_stars(rim_lc_file) %>%
  st_crop(rim_extent)
rim_lc[is.na(rim_mcdc)==TRUE] <- NA
rim_lc<- st_as_sf(rim_lc, as_points = TRUE, merge = FALSE) %>%
  mutate(classes = lut_lc[lc_h08v05_2012.tif],
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL)

ggplot() +
  geom_raster(data=rim_lc, aes(x=x,y=y,fill=classes))

rim_tmin_path <-"/home/a/data/PRISM/daily/tmin/y2013m8"

if(!file.exists("data/rim/rim_afv.gpkg")){
  rim_afm <- modis_af %>%
    filter(substr(ACQ_DATE,1,4)%>% as.numeric > 2012 &
             substr(ACQ_DATE,1,4)%>% as.numeric < 2014)  %>%
    st_transform(crs=st_crs(rim_extent))%>%
    st_intersection(rim_extent) %>%
    filter(ACQ_DATE>as.Date("2013-08-16"), ACQ_DATE < as.Date("2013-09-01")) %>% 
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
            daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
            daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
            dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
            daynight = as.character(daynight)) %>%
    dplyr::select(-DAYNIGHT)%>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
    group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  day0 <- min(rim_afm$ACQ_DATE -1)
  
  ll <- list()
  for(d in 1:length(unique(rim_afm$event_day))){
    ed <- as.character(day0+unique(rim_afm$event_day)[d])
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(rim_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(rim_prsm, "Spatial")) %>%
      projectRaster(as(rim_mcd, "Raster")) 
    
    ll[[d]] <- filter(rim_afm, event_day == unique(rim_afm$event_day)[d])%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  rim_afm <- do.call("rbind",ll)

## afv 
  rim_afv <- viirs_af %>%
    filter(substr(ACQ_DATE,1,4)%>% as.numeric > 2012 &
             substr(ACQ_DATE,1,4)%>% as.numeric < 2014)  %>%
    st_transform(crs=st_crs(rim_extent))%>%
    st_intersection(rim_extent)%>%
    filter(ACQ_DATE>as.Date("2013-08-16"), ACQ_DATE<as.Date("2013-09-01")) %>%
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
           daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
           daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
           dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
           daynight = as.character(daynight)) %>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
    group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  day0 <- min(rim_afv$ACQ_DATE -1)
  
  ll <- list()
  for(d in 1:length(unique(rim_afv$event_day))){
    ed <- as.character(day0+unique(rim_afv$event_day)[d])
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(rim_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(rim_prsm, "Spatial")) %>%
      projectRaster(as(rim_mcd, "Raster")) 
    
    ll[[d]] <- filter(rim_afv, event_day == unique(rim_afv$event_day)[d])%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  rim_afv <- do.call("rbind",ll)
  
  st_write(rim_afv, "data/rim/rim_afv.gpkg", delete_dsn = TRUE)
  st_write(rim_afm, "data/rim/rim_afm.gpkg", delete_dsn = TRUE)
}else{
  rim_afv <- st_read("data/rim/rim_afv.gpkg")
  rim_afm <- st_read("data/rim/rim_afm.gpkg")
}

ggplot() +
  geom_raster(data=rim_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = c("gold","darkgreen", "lightyellow", 
                               "lightgreen", "grey70", "lightblue",
                               "olivedrab"),
                    name = "Land Cover")+
  theme_void()+
  theme(#legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  geom_sf(data=rim_afm, aes(color = tmin), shape = 19, alpha = 0.95) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = -1, "Minimum Temperature",end = .85, begin=0  )+
  #scale_colour_gradient(low = "white", high = "red")+
  geom_text(data = rim_afm,
            x=-10510000, y=4200000, 
            aes(label=paste("FRP:",round(mean_FRP,0))))+
  facet_wrap(~event_day+daynight, ncol=4, labeller = mylabels) +
  ggtitle("Rim Fire: MODIS Active Fire")

ggplot() +
  geom_raster(data=rim_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = c("gold","darkgreen", "lightyellow", 
                               "lightgreen", "grey70", "lightblue",
                               "olivedrab"),
                    name = "Land Cover")+
  theme_void()+
  theme(#legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  geom_sf(data=rim_afv, aes(color = tmin), shape = 19, alpha = 0.95) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = -1, "Minimum Temperature",end = .85, begin=0  )+
  #scale_colour_gradient(low = "white", high = "red")+
  geom_text(data = rim_afv,
            x=-10510000, y=4200000, 
            aes(label=paste("FRP:",round(mean_FRP,0))))+
  facet_wrap(~event_day+daynight, ncol=4, labeller = mylabels) +
  ggtitle("Rim Fire: VIIRS Active Fire")

p1 <- ggplot() +
  geom_raster(data = hpmdf, aes(x=x,y=y,fill = burn_doy))+
  ggtitle("Modis Burned Area")
p2 <-ggplot()+
  geom_sf(data = hp_afv, alpha=0.5, aes(color = daynight)) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Viirs Active Fire")
p3 <- ggplot() +
  geom_sf(data = hp_afm, alpha=0.5, aes(color = daynight)) +
  scale_color_viridis_d(option = "C") +
  ggtitle("Modis Active Fire")
ggarrange(p1,p2,p3)

both_afs <- rbind(dplyr::select(rim_afv, LONGITUDE, LATITUDE, daynight, dt),
                  dplyr::select(rim_afm, LONGITUDE, LATITUDE, daynight, dt))

anim<-ggplot(both_afs, aes(x=LONGITUDE,y=LATITUDE))+ 
  geom_point(aes(color = daynight)) + # the business
  theme_void() + # this is a theme without any annoying lines  
  scale_color_manual(name = "", values=c("orange","darkblue")) + # setting the color scheme and naming the legend
  labs(title = 'Day: {frame_time}') + # this is setting the label on top
  transition_time(dt) +
  shadow_mark()

aa<-gganimate::animate(anim, fps=2, nframes = length(unique(both_afs$dt)));aa
anim_save(aa, filename="data/rim/af_both.gif")


# pioneer fire (idaho 2016) ============================================================
pio_p<- st_read(file.path(local_data, "fire", "mtbs", "mtbs_perims_DD.shp")) %>%
  filter(Fire_Name == "PIONEER", Year == 2016) %>%
  st_transform(crs = st_crs(hpm))

pio_mcdf <- "bd_numeric_2016_245_h09v04.tif" 
if(!file.exists(file.path("data/pio",pio_mcdf))){
  system(paste("aws s3 cp",
               file.path(s3path, pio_mcdf),
               file.path("data/pio",pio_mcdf)))
} #manually messed with this a bit (there are three files)
pio_mcd <-read_stars(file.path("data/pio",pio_mcdf))

pio_extent <- pio_p %>%
  st_transform(crs = st_crs(pio_mcd))%>%
  st_buffer(dist = 1000) %>%
  st_bbox() %>%
  st_as_sfc()

pio_prsm <- pio_extent %>%
  st_transform(prism_crs)%>%
  st_bbox()%>%
  st_as_sfc()

pio_mcdc <- pio_mcd %>%
  st_crop(pio_extent)
pio_mcdc[pio_mcdc<1]<-NA
names(pio_mcdc) <- "burn_date"

pio_lc_file <- "/home/a/projects/tmin/data/lc/lc_h09v04_2015.tif"
sds <- '/home/a/data/MCD12Q1/MCD12Q1.A2015001.h09v04.006.2018055062819.hdf' %>%
  gdalUtils::get_subdatasets()
gdalUtils::gdal_translate(sds[1], dst_dataset = pio_lc_file)

pio_lc<- read_stars(pio_lc_file)[pio_p] 

pio_lc <- pio_lc%>%
  st_warp(crs=st_crs(modis_af),) %>%
  st_as_sf(as_points = TRUE, merge = FALSE) %>%
  mutate(classes = lut_lc[lc_h09v04_2015.tif],
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL)

ggplot()+
  geom_raster(data = pio_lc, aes(x=x,y=y, fill = classes))


pio_tmin_path <- "/home/a/data/PRISM/daily/tmin/y2016m7"
if(!file.exists("data/pio/pio_afv.gpkg")){
  pio_afm <- modis_af %>%
    filter(ACQ_DATE>as.Date("2016-07-16"), ACQ_DATE < as.Date("2016-09-01")) %>% 
    st_transform(crs=st_crs(pio_extent))%>%
    st_intersection(pio_extent) %>%
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
           daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
           daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
           dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
           daynight = as.character(daynight)) %>%
    dplyr::select(-DAYNIGHT)%>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
    group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  day0 <- min(pio_afm$ACQ_DATE -1)
  
  ll <- list()
  for(d in 1:length(unique(pio_afm$event_day))){
    ed <- as.character(day0+unique(pio_afm$event_day)[d])
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(pio_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(pio_prsm, "Spatial")) %>%
      projectRaster(as(pio_mcd, "Raster")) 
    
    ll[[d]] <- filter(pio_afm, event_day == unique(pio_afm$event_day)[d])%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  pio_afm <- do.call("rbind",ll)
  
  ## afv 
  pio_afv <- viirs_af %>%
    filter(ACQ_DATE>as.Date("2016-07-16"), ACQ_DATE<as.Date("2016-09-01")) %>%
    st_transform(crs=st_crs(pio_extent))%>%
    st_intersection(pio_extent)%>%
    mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),
           daynight = ifelse(as.numeric(ACQ_TIME)<700, "_night", "day"),
           daynight = ifelse(as.numeric(ACQ_TIME)>2100, "_night", daynight),
           dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
           daynight = as.character(daynight)) %>%
    mutate(event_day = (ACQ_DATE - min(ACQ_DATE))+1,
           event_day = ifelse(as.numeric(ACQ_TIME)>2100, event_day+1, event_day))%>%
    group_by(event_day, daynight) %>%
    mutate(mean_FRP = mean(FRP)) %>%
    ungroup()
  
  day0 <- min(pio_afv$ACQ_DATE -1)
  
  ll <- list()
  for(d in 1:length(unique(pio_afv$event_day))){
    ed <- as.character(day0+unique(pio_afv$event_day)[d])
    year <- substr(ed,1,4)
    month<- substr(ed,6,7)
    day <- substr(ed, 9,10)
    
    rf <- file.path(pio_tmin_path, paste0(meat, year, month, day, bits)) 
    r <- raster(rf) %>%
      crop(as(pio_prsm, "Spatial")) %>%
      projectRaster(as(pio_mcd, "Raster")) 
    
    ll[[d]] <- filter(pio_afv, event_day == unique(pio_afv$event_day)[d])%>%
      mutate(tmin = raster::extract(r,.)%>% as.numeric)
  }
  pio_afv <- do.call("rbind",ll)
  
  st_write(pio_afv, "data/pio/pio_afv.gpkg", delete_dsn = TRUE)
  st_write(pio_afm, "data/pio/pio_afm.gpkg", delete_dsn = TRUE)
}else{
  pio_afv <- st_read("data/pio/pio_afv.gpkg")
  pio_afm <- st_read("data/pio/pio_afm.gpkg")
}

ggplot()+
  geom_sf(data=st_transform(pio_p, 4326))+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_sf(data=st_transform(pio_afv,4326), 
          aes(color = tmin#, 
              #shape = daynight
              ), shape =16,
          alpha = 0.50) +
  scale_color_viridis(option = "B")+
  coord_sf(datum=NA)+
  facet_wrap(~event_day+daynight, nrow = 6, labeller = mylabels) +
  ggtitle("Pioneer Fire: VIIRS Active Fire")+
  ggsave("images/tmin_pioneer_viirs.png", width = 11, height = 6)
  
  
ggplot()+
  geom_sf(data=st_transform(pio_p, 4326))+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_sf(data=st_transform(pio_afm,4326), 
          aes(color = tmin), 
          shape =16,
          alpha = 0.50) +
  scale_color_viridis(option = "B")+
  coord_sf(datum=NA)+
  facet_wrap(~event_day+daynight, nrow = 6, labeller = mylabels) +
  ggtitle("Pioneer Fire: MODIS Active Fire")+
  ggsave("images/tmin_pioneer_modis.png", width = 11, height = 6)
  
# day vs night time series

modis_dn <- modis_af %>%
  st_set_geometry(NULL) %>%
  mutate(year = substr(ACQ_DATE, 1,4)) %>%
  group_by(DAYNIGHT, year) %>%
  summarise(count = n()) %>%
  ungroup()

spr <- modis_dn %>%
  filter(year>2002, year < 2019) %>%
  spread(key = DAYNIGHT, value = count) %>%
  mutate(pct_n = N/D*100)

ggplot(spr, aes(x=year, y = pct_n)) +
  geom_point()
