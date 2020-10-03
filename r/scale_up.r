# tring to use only mtbs polygon and lc
source("r/tmin_functions.R")
libs<- c("tidyverse", "sf", "raster", "stars","viridis")
lapply(libs, library, character.only = TRUE)

local_data <- "/home/a/data"

tz <- st_read("/home/a/data/timezones/world_timezones.shp") %>%
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


lc_e<- read_stars("/home/a/data/MCD12Q1_mosaics/usa_lc_mosaic_2001.tif") %>%
  st_bbox() %>%
  st_as_sfc()


mtbs <- st_read("/home/a/data/fire/mtbs/mtbs_perims_DD.shp") %>%
  filter(Year > 2001) %>%
  mutate(start_date = as.Date(paste(Year,
                                    str_pad(StartMonth,2, pad="0"),
                                    str_pad(StartDay,2, pad = "0"),
                                    sep = "-"))) %>%
  st_intersection(st_transform(lc_e, crs=st_crs(.)))

years = 2002:2016

res_list <- list()
counter <- 1
for(year in years){
  mtbs_s <- filter(mtbs,Year == year)
  r_native <- read_stars(paste0("/home/a/data/MCD12Q1_mosaics/usa_lc_mosaic_",
                                as.numeric(year)-1,".tif"))
  tmin_path <- file.path(local_data, "PRISM/daily/tmin", paste0("y", year))
  for(i in 1:nrow(mtbs_s)){
    
    dd <- get_af_event(polygon = mtbs_s[i,],
               p_crs = st_crs(mtbs_s),
               af_product = modis_af,
               tmin_path = tmin_path,
               r_native = r_native,
               start = mtbs_s[i,]$start_date,
               end = mtbs_s[i,]$start_date +4,
               timezone = tz)
    if(!is.null(dd)){
      res_list[[counter]] <- dd
      counter <- counter +1
    }
  print(round(i/nrow(mtbs)*100,1))
  }
}
res <- do.call("rbind", res_list) %>%
  mutate(lc_class = lut_lc[lc])
st_write(res, "data/res1_june6.gpkg", delete_dsn = T)
ggplot(res, aes(x=event_day, y = tmin)) +
  geom_point(aes(color = lc_class)) +
  facet_wrap(~lc_class)

ggplot(res, aes(x=as.factor(event_day), y = tmin)) +
  geom_boxplot(aes(color = lc_class)) +
  facet_wrap(~lc_class)

ggplot(res, aes(x=tmin, y =FRP, color = lc_class)) +
  geom_point(alpha = 0.25) +
  geom_smooth() +
  facet_wrap(~lc_class)

res_events<- res %>%
  group_by(Fire_Name) %>%
  summarise(tmin_mean = mean(tmin, na.rm=T),
            frp_mean = mean(FRP),
            pixels = n()) 

ggplot(res_events%>%filter(pixels>10), aes(x=tmin_mean, y = pixels)) +
  geom_point() +
  geom_smooth(method="lm")
