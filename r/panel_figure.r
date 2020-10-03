# panel figure

# setup ========================================================================
libs<- c("tidyverse", "sf", "raster", "stars","viridis",
         "gdalUtils")
lapply(libs, library, character.only = TRUE)
source("r/tmin_functions.R")
# parameters ===================================================================
s3path <- "s3://earthlab-natem/modis-burned-area/MCD64A1/C6/tif_converted_alltiles"
local_data <- "/home/a/data"

mylabels = label_bquote(
  cols = .(daynight) ~ .(event_day)
)

tz <- st_read("/home/a/data/background/timezones/world_timezones.shp") %>%
  mutate(direction = ifelse(substr(UTC_OFFSET,4,4) =="-", -1,1),
         magnitude_hrs = as.numeric(substr(UTC_OFFSET, 5, 6))*60*60,
         magnitude_mns = as.numeric(substr(UTC_OFFSET, 8, 9))*60,
         total_seconds_offset = direction*(magnitude_hrs+magnitude_mns)) %>%
  mutate(total_seconds_offset = replace(total_seconds_offset,is.na(total_seconds_offset), 0)) %>%
  dplyr::select(total_seconds_offset)

options(stringsAsFactors = F)
prism_crs <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
meat <- "PRISM_tmin_stable_4kmD1_"
bits <- "_bil.bil"

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

# af data import ===============================================================

viirs_af <- st_read(file.path(local_data,"fire", "modis_viirs_active_fire", 
                              "CUS_V1","fire_archive_V1_54276.shp"))
modis_af <- st_read(file.path(local_data,"fire", "modis_viirs_active_fire", 
                              "CUS_M6","fire_archive_M6_54275.shp"))

# hotpot data wrangling ========================================================
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

hpmc <- st_crop(hpm, hp_bbox)
hpmc[hpmc<1]<-NA
names(hpmc) <- "burn_date"

ggplot() +
  geom_stars(data=hpmc, aes(fill=burn_date, x=x,y=y)) +
  theme_void()+
  coord_fixed()+
  ggsave("~/projects/gis_talk/hp_mcd64.png")

hp_prsm <- hp_bbox %>%
  st_transform(prism_crs)%>%
  st_bbox()%>%
  st_as_sfc()

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

hp_afm <- get_af_event(polygon = hp_bbox, p_crs = st_crs(hp_bbox),
                       af_product =  modis_af, 
                       tmin_path = "/home/a/data/PRISM/daily/tmin/y2016m7", 
                       r_prism = hp_prsm,r_native = hpmc, start = "2016-07-23",
                       end = "2016-07-26")

hp_afv <- get_af_event(polygon = hp_bbox, p_crs = st_crs(hp_bbox), 
                       af_product =  viirs_af, 
                       tmin_path = "/home/a/data/PRISM/daily/tmin/y2016m7", 
                       r_native = hpmc, start = "2016-07-23",
                       end = "2016-07-26")

# hotpot plots =================================================================
ggplot() +
  geom_raster(data=hp_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = c("darkgreen","turquoise4", "lightyellow", "olivedrab",
                               "lightgreen", "grey70", "lightblue"),
                    name = "Land Cover")+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_sf(data=hp_afv,
          aes(color = tmin, size = FRP), shape = 1, alpha = 0.95, stroke = 1) +
  coord_sf(datum=NA)+
  geom_text(data = hp_afv,
            x=-10852000, y=3816500, aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  scale_color_viridis(option="B", direction = 1, "Minimum Temperature",end = .85, begin=0  )+
  facet_wrap(~event_day + daynight, ncol=2, labeller = mylabels) +
  ggtitle("Hot Pot Fire: VIIRS Active Fire") +
  ggsave("images/tmin_hotpot_viirs_lc.png")


ggplot() +
  geom_raster(data=hp_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = c("darkgreen","turquoise4", "lightyellow", "olivedrab",
                               "lightgreen", "grey70", "lightblue"),
                    name = "Land Cover")+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_sf(data=hp_afm, aes(color = tmin, size = FRP), shape = 1, alpha = 0.75, stroke =1) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = 1, "Minimum Temperature",
                      end = .85, begin=0  )+
  geom_text(data = hp_afm,
            x=-10852000, y=3816500, 
            aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  facet_wrap(~event_day+daynight, ncol=2, labeller = mylabels) +
  ggtitle("Hot Pot Fire: MODIS Active Fire") +
  ggsave("images/tmin_hotpot_modis_lc.png")

# plots for jan 2020 figure ============================================

hp_viirs <- viirs_af %>%
  filter(ACQ_DATE>=as.Date("2016-07-23"), ACQ_DATE<=as.Date("2016-07-26")) %>%
  st_transform(st_crs(hp_bbox))%>%
  st_intersection(hp_bbox) %>%
  st_intersection(st_transform(tz, crs=st_crs(hp_bbox))) %>%
  mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),#,
         dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
         lt = dt+total_seconds_offset,
         acq_date_l = as.Date(substr(lt, 1,10)),
         acq_hour_l = as.numeric(substr(lt,12,13)),
         daynight = ifelse(acq_hour_l<7, "night", "day"),
         daynight = ifelse(acq_hour_l>=21, "night", daynight),
         daynight = as.character(daynight)
  ) %>%
  mutate(event_day = (acq_date_l - min(acq_date_l)),
         event_day = ifelse(acq_hour_l>=21, event_day+1, event_day),
         daynight = factor(daynight, levels = c("day", "night")),
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2],
         dnd = paste(daynight, event_day+1) %>%
           factor(levels=c("day 1", "night 1", "day 2", "night 2", "day 3", "night 3",
                           "day 4", "night 4"))
  )

hp_lc <- mutate(hp_lc, fire_perimeter = "MCD64 fire perimeter")

ggplot() +
  theme_void()+
  geom_raster(data=hp_lc, aes(x=x,y=y, fill=fire_perimeter))+
  geom_point(data=hp_viirs, aes(x=x,y=y,color = dnd, shape = dnd), 
          # shape = 18, 
          size =6,alpha = 0.75, stroke =1) +
  scale_shape_manual(values = c(18,18,1,1,17,17,3,3)) +
  scale_color_manual(values = rep(c("red","blue"),4))+
  # facet_wrap(~dnd)+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank())+
  scale_fill_manual(values="grey85")+
  ggtitle("2016 Hot Pot Fire: VIIRS Active Fire Detections") +
  ggsave("images/jan_2020_draft_figure.png")


# camp fire data wrangling =====================================================
camp_p <- st_read("data/camp/ca_camp_20181125_1822_dd83.shp")
camp_tmin_path <- path.expand("~/data/PRISM/daily/tmin/y2018m11")

camp_hdf_file <- "/home/a/projects/tmin/data/camp/MCD64A1.A2018244.h08v05.006.2018316192508.hdf"
camp_mcd_file<- path.expand("~/projects/tmin/data/camp/mcd_camp.tif")
sds <- camp_hdf_file %>%
  gdalUtils::get_subdatasets()
gdalUtils::gdal_translate(sds[1], dst_dataset = camp_mcd_file)

camp_mcd <- read_stars(camp_mcd_file)
camp_mcdp <- st_transform(camp_p, st_crs(camp_mcd))
camp_mcdc <- camp_mcd[camp_mcdp]

camp_prism <- camp_p%>%
  st_transform(prism_crs)%>%
  st_bbox()%>%
  st_as_sfc()

camp_afm <- get_af_event(polygon = camp_p,
                         p_crs = st_crs(camp_p),
                         af_product = modis_af,
                         tmin_path = camp_tmin_path,
                         r_prism = camp_prism,r_native = camp_mcdc,
                         start = "2018-11-08",end = "2018-11-12")
camp_afv <- get_af_event(polygon = camp_p,
                         p_crs = st_crs(camp_p),
                         af_product = viirs_af,
                         tmin_path = camp_tmin_path,
                         r_prism = camp_prism,r_native = camp_mcdc,
                         start = "2018-11-08",end = "2018-11-12")

camp_lc_file <- "/home/a/projects/tmin/data/lc/lc_h08v05_2016.tif"

camp_lc<- read_stars(camp_lc_file) %>%
  st_crop(camp_mcdp)
names(camp_lc) <- "lc"
camp_lc[is.na(camp_mcdc)==TRUE] <- NA
camp_lc<- st_as_sf(camp_lc, as_points = TRUE, merge = FALSE) %>%
  mutate(classes = lut_lc[lc],
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL)

# camp plots ===================================================================
ggplot() +
  geom_raster(data=camp_lc, aes(x=x,y=y, fill=classes))+
  scale_fill_manual(values = c( "turquoise3","orange","turquoise4", "turquoise4",
                                "darkgreen","lightyellow", "olivedrab","olivedrab",
                                "lightblue", "lightgreen", "grey70", "lightgreen"),
                    name = "2016 Land Cover")+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_sf(data=camp_afm %>% filter(event_day < 5) %>% st_transform(crs = st_crs(camp_mcd))
          , aes(color = tmin, size=FRP), shape = 1, alpha = 0.95, stroke =1) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = 1, "Minimum Temperature",end = .85, begin=0 )+
  # geom_text(data = camp_afm %>% filter(event_day < 5),
  #           x=-10852000, y=3816500, 
  #           aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  facet_wrap(~event_day+daynight, ncol=2, labeller = mylabels) +
  ggtitle("Camp Fire: MODIS Active Fire") +
  ggsave("images/tmin_camp_modis_lc.png")

ggplot() +
  geom_raster(data=camp_lc, aes(x=x,y=y, fill=classes))+
   scale_fill_manual(values = c( "turquoise3","orange","turquoise4", "turquoise4",
                               "darkgreen","lightyellow", "olivedrab","olivedrab",
                               "lightblue", "lightgreen", "grey70", "lightgreen"),
                     name = "2016 Land Cover")+
  theme_void()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_sf(data=camp_afv %>% filter(event_day < 5) %>% st_transform(crs = st_crs(camp_mcd))
          , aes(color = tmin, size=FRP), shape = 1, stroke =1, alpha = 0.75) +
  coord_sf(datum=NA)+
  scale_color_viridis(option="B", direction = 1, "Minimum Temperature",end = .85, begin=0 )+
  # geom_text(data = camp_afm %>% filter(event_day < 5),
  #           x=-10852000, y=3816500, 
  #           aes(label=paste("mean FRP:",round(mean_FRP,2))))+
  facet_wrap(~event_day+daynight, ncol=2, labeller = mylabels) +
  ggtitle("Camp Fire: VIIRS Active Fire") +
  ggsave("images/tmin_camp_viirs_lc.png")
