# camp fire figure jan 2020

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

# af data import ===============================================================
# this takes a loooooong time (5 minutes per file)
viirs_af <- st_read(file.path(local_data,"fire", "modis_viirs_active_fire", 
                              "CUS_V1","fire_archive_V1_54276.shp"))
modis_af <- st_read(file.path(local_data,"fire", "modis_viirs_active_fire", 
                              "CUS_M6","fire_archive_M6_54275.shp"))

camp_p <- st_read("data/camp/ca_camp_20181125_1822_dd83.shp")
start = "2018-11-08"
end = "2018-11-14"

library(lubridate)
camp_viirs <- viirs_af %>%
  filter(ACQ_DATE>=as.Date(start), ACQ_DATE<=as.Date(end)) %>%
  st_transform(st_crs(camp_p))%>%
  st_intersection(camp_p) %>%
  st_intersection(st_transform(tz, crs=st_crs(camp_p))) %>%
  dplyr::rename_all(tolower) %>%
  dplyr::mutate(acq_hour = as.numeric(substr(acq_time, start = 1, stop = 2)),
                acq_min = as.numeric(substr(acq_time, start = 3, stop = 4)),
                acq_datetime = ymd_hm(paste0(acq_date, " ", acq_hour, ":", acq_min)),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night"),
                x = st_coordinates(.)[,1],
                y = st_coordinates(.)[,2])
  

camp_modis <- modis_af %>%
  filter(ACQ_DATE>=as.Date(start), ACQ_DATE<=as.Date(end)) %>%
  st_transform(st_crs(camp_p))%>%
  st_intersection(camp_p) %>%
  st_intersection(st_transform(tz, crs=st_crs(camp_p))) %>%
  mutate(daynight = ifelse(DAYNIGHT == "D", "day", "night"),
         daynight = factor(daynight, levels = c("day", "night")),
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2])


camp_m_and_v <- rbind(dplyr::select(camp_viirs, daynight, x, y),
      dplyr::select(camp_modis, daynight, x,y))

paradise <- data.frame(y = 39.767380, x=-121.633759, label = "Paradise, CA")

ggplot() +
  theme_void()+
  geom_sf(data=camp_p)+
  geom_point(data=camp_viirs, aes(x=x,y=y,color = daynight), 
             shape = 18, 
             size = 4,alpha = 0.75) +
  geom_sf(data=camp_p, fill="transparent")+
  geom_point(data = paradise, aes(x=x,y=y), stroke =1, 
             color = "black", shape = 8, size = 8)+
  geom_text(data=paradise, aes(x=x,y=y,label=label), 
            hjust="left", 
            nudge_x = .012,
            fontface="bold", 
            size=6)+
  scale_color_discrete(name="2018 Camp Fire:\nVIIRS Active Fire Detections Nov 8-14")+
  facet_wrap(~daynight)+
  theme(legend.position = c(0,0.9),
        legend.justification = c(0,1),
        legend.title = element_text(size=15,face = "bold"))+
  scale_fill_manual(values="grey85")+
  # ggtitle() +
  ggsave("images/jan_2020_draft_figure.png")
