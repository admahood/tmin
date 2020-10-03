# australian wildfires fig

library(tidyverse)
library(sf)
library(raster)
library(ggsn)
library(stars)
library(ggpubr)
library(lubridate)
library(cowplot)

theme_set(theme_void())

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# add a scale bar!
# get landcover type thresholds

# landcover classification =====================================================
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
             "13"=  "Urban and Built-up Lands",
             "14"= "Cropland/Natural  Vegetation  Mosaics",
             "15"=  "Permanent Snow and Ice","16"= "Barren",
             "17"=  "Water Bodies", "255"= "Unclassified")

lut_colors <- c("Evergreen Needleleaf Forests" = "#006400",
                "Evergreen Broadleaf Forests" = "#228B22",
                "Deciduous Needleleaf Forests"= "#458B00",
                "Deciduous Broadleaf Forests" = "#008B45", 
                "Mixed Forests"= "#3CB371",
                "Closed Shrublands" = "#6E8B3D",
                "Open Shrublands" = "#9ACD32", 
                "Woody Savannas" = "#6B8E23",
                "Savannas" = "#8B8B00",
                "Grasslands" = "#CDC673", 
                "Permanent Wetlands" = "#00868B", 
                "Croplands" = "#EE9572",
                "Urban and Built-up Lands" = "#B3B3B3", 
                "Cropland/Natural  Vegetation  Mosaics" = "#EE8262",
                "Permanent Snow and Ice" = "#FFFFFF",
                "Barren" = "#DEB887",
                "Water Bodies" = "#87CEEB",
                "Unclassified" = "#BEBEBE"
)

snowy_lc_file <- "data/MCD12Q1_tifs/h29v12.tif"
braz_lc_file <- "data/MCD12Q1_tifs/h12v09.tif"
tubbs_lc_file <- "data/MCD12Q1_tifs/h08v05.tif"

modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

# snowy complex ================================================================

snowy <- st_read("data/snowy_complex.gpkg") %>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))%>%
  mutate(acq_date = as.Date(acq_date)) %>%
  filter(acq_date > as.Date("2019-12-28"),
         acq_date < as.Date("2020-01-6"))

snowy_n <- snowy %>%
  mutate(acq_datetime = round_date(ymd_hms(acq_datetime), unit = "hour")) %>%
  group_by(acq_datetime) %>%
  summarise(value = n(),
            daynight = first(daynight)) %>%
  ungroup() %>%
  mutate(variable = "Detections",
         lty = 0) %>%
  mutate(label_y = max(value)*0.85) %>%
  dplyr::select(datetime = acq_datetime, variable, value,lty, label_y, daynight) %>%
  st_set_geometry(NULL)

snowy_clim <- read_csv("data/aussie2020.csv",col_names = F) %>%
  dplyr::rename(day=X1, hour=X2, value = X3)%>% 
  mutate(year = ifelse(day>300,2019,2020),
         month = ifelse(day>300,"12","01"),
         value = value/10) %>%
  mutate(day = replace(day, day > 300, day-334),
         day = str_pad(day, 2, "left", "0"),
         date = paste0(year, "-",month,"-",day) %>% as.Date(),
         datetime = lubridate::ymd_hms(paste0(date," ", hour, ":00:00")),
         variable = "VPD (kPa)",
         lty=2) %>%
  mutate(label_y = max(value)*0.85) %>%
  filter(date > as.Date("2019-12-28"),
         date < as.Date("2020-01-6"))%>%
  dplyr::mutate(latitude = first(snowy$latitude),
                longitude = first(snowy$longitude),
                acq_datetime = datetime,
                acq_hour = hour(acq_datetime),
                acq_min = minute(acq_datetime),
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
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  dplyr::select(datetime, variable, value,lty, label_y,daynight) %>%
  rbind(snowy_n)

lc <- snowy %>%
  mutate(lc=raster::extract(x=raster(snowy_lc_file), y=.)) %>%
  mutate(lc_name = lut_lc[lc]) %>%
  group_by(lc_name) %>%
  summarise(percent_lc = round((n()/nrow(.))*100)) %>%
  ungroup() %>%
  st_set_geometry(NULL)

stats <- snowy %>%
  st_set_geometry(NULL) %>%
  group_by(daynight) %>%
  summarise(percent_dn = round((n()/nrow(.))*100,2)) %>%
  ungroup() 


daterange <- snowy$acq_date %>% as.Date() %>% range()

snowy_area <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Australia") %>%
  st_crop(snowy %>% st_buffer(1))

australia <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Australia") 

bbox <- st_bbox(snowy)%>% as.numeric()

if(!file.exists("data/fishnet.RDS")){
  fishnet <- st_bbox(snowy)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01)%>%
    as("Spatial")%>%
    raster::extract(x=raster(snowy_lc_file), y=., 
                    fun = function(x,...) getmode(x), 
                    method="simple") 
  saveRDS(fishnet, "data/fishnet.RDS")
}else{
  fishnet <- readRDS("data/fishnet.RDS")
}

if(!file.exists("data/fishnet_lc.RDS")){
  fishnet_lc <- st_bbox(snowy) %>%
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01) %>%
    st_as_sf() %>%
    mutate(lc= fishnet,
           classes = lut_lc[lc])
  saveRDS(fishnet_lc, "data/fishnet_lc.RDS")
}else{
  fishnet_lc<-readRDS("data/fishnet_lc.RDS")
}


locator_box <- st_bbox(snowy) %>% st_as_sfc()

main_plot_s <- ggplot() +
  # geom_sf(data = fishnet_lc, 
  #         aes(fill=classes),
  #         color = "transparent", alpha = 0.5)+
  geom_sf(data = snowy_area, fill = "transparent")+
  geom_sf(data = snowy, 
          aes(color=daynight, alpha = daynight), show.legend = "point") +
  scale_fill_manual(values = lut_colors,
                    name = "Land Cover")+
  scale_alpha_manual(values = c(0.5,0.15))+
  scale_color_manual(values = c("blue", "red"))+
  ggsn::scalebar(data = snowy, location = "topleft", model = "WGS84",
                 dist = 30, dist_unit = "km",transform = TRUE,st.dist = 0.05) +
  guides(fill = FALSE)+
  xlim(c(bbox[c(1,3)])) +
  ylim(c(bbox[c(2,4)])) +
  ggtitle(paste("A. Snowy Complex. October 2, 2019 - January 19, 2020")) +
  theme(legend.position ="none",
        # legend.justification = c(0.95,0),
        # legend.title = element_blank(),
        panel.border = element_rect(color="black", fill=NA))

locator_plot_s <- ggplot(australia) +  
  theme(panel.border = element_rect(color="black", fill=NA),
        plot.background = element_rect(fill="white"))+
  geom_sf(fill="white") +
  geom_sf(data=st_centroid(locator_box), color = "red", size=2)+
  ylim(c(-43, -11))+
  xlim(c(114,153))

inset_s <- ggplot(snowy_clim, aes(x=datetime, y=value)) +
  geom_line(aes(color=daynight, group = "none")) +
  geom_abline(aes(slope = 0,intercept = 1, lty=as.factor(lty)))+
  scale_linetype_manual(values = c(0,2))+
  scale_color_manual(values = c("red", "blue"))+
  facet_wrap(~variable, scales = "free_y", 
             nrow = 2, strip.position = "left") +
  xlab("Date") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

snowy_cow <- ggdraw() +
  draw_plot(main_plot_s) +
  draw_plot(locator_plot_s, 0.8,.7,.25,.25)

# some fires in brazil ==========================================================
# pull in global forest cover and/or mod17 to make sure it's not an ag fire

braz <- st_read("data/brazil_fire.gpkg")%>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))%>%
  mutate(acq_date = as.Date(acq_date)) %>%
  filter(acq_date > as.Date("2019-07-28"),
         acq_date < as.Date("2019-08-18"))

braz_n <- braz %>%
  mutate(acq_datetime = round_date(ymd_hms(acq_datetime), unit = "hour")) %>%
  group_by(acq_datetime) %>%
  summarise(value = n()) %>%
  ungroup() %>%
  mutate(variable = "Detections",
         lty = 0,
         label_y = max(value)*0.85) %>%
  dplyr::select(datetime = acq_datetime, variable, value,lty, label_y) %>%
  st_set_geometry(NULL)

braz_clim <- read_csv("data/brazil2019vpd.csv") %>%
  dplyr::rename(value = `vpd(hPa)`)%>% 
  mutate(value = value/10,
         year = 2019,
         month = ifelse(day>212,"08","07")) %>%
  mutate(day = ifelse(day > 212, day-212, day-181))%>%
  mutate(day = str_pad(day, 2, "left", "0"),
         date = paste0(year, "-",month,"-",day) %>% as.Date(),
         datetime = lubridate::ymd_hms(paste0(date," ", hour, ":00:00")),
         variable = "VPD (kPa)",
         lty = 2,
         label_y = max(value)*0.85) %>%
  filter(date > as.Date("2019-07-28"),
         date < as.Date("2019-08-18"))%>%
  dplyr::select(datetime, variable, value,lty, label_y) %>%
  rbind(braz_n)



lc <- braz %>%
  mutate(lc=raster::extract(x=raster(braz_lc_file), y=.)) %>%
  mutate(lc_name = lut_lc[lc]) %>%
  group_by(lc_name) %>%
  summarise(percent_lc = round((n()/nrow(.))*100)) %>%
  ungroup() %>%
  st_set_geometry(NULL)

stats <- braz %>%
  st_set_geometry(NULL) %>%
  group_by(daynight) %>%
  summarise(percent_dn = round((n()/nrow(.))*100,2)) %>%
  ungroup() 


daterange <- braz$acq_date %>% as.Date() %>% range()

braz_area <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Brazil") %>%
  st_crop(braz %>% st_buffer(1))

brazil <- st_read("~/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  filter(NAME_EN == "Brazil") 


bbox <- st_bbox(braz)%>% as.numeric()

locator_box <- st_bbox(braz) %>% st_as_sfc()
if(!file.exists("data/braz_fishnet.RDS")){
  fishnet <- st_bbox(braz)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.005)%>%
    as("Spatial")%>%
    raster::extract(x=raster(braz_lc_file), y=., 
                    fun = function(x,...) getmode(x), 
                    method="simple") 
  saveRDS(fishnet, "data/braz_fishnet.RDS")
}else{
  fishnet <- readRDS("data/braz_fishnet.RDS")
}

if(!file.exists("data/braz_fishnet_lc.RDS")){
  fishnet_lc <- st_bbox(braz) %>%
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.005) %>%
    st_as_sf() %>%
    mutate(lc= fishnet,
           classes = lut_lc[lc])
  saveRDS(fishnet_lc, "data/braz_fishnet_lc.RDS")
}else{
  fishnet_lc<-readRDS("data/braz_fishnet_lc.RDS")
}

main_plot_b <- ggplot() +
  # geom_sf(data = fishnet_lc, 
  #         aes(fill=classes),
  #         color = "transparent", alpha = 0.5)+
  geom_sf(data = braz_area, fill = "transparent")+
  geom_sf(data = braz, 
          aes(color=daynight, alpha = daynight), show.legend = "point") +
  scale_fill_manual(values = lut_colors,
                    name = "Land Cover")+
  scale_alpha_manual(values = c(0.75,0.5))+
  scale_color_manual(values = c("blue", "red"))+
  guides(fill = FALSE, color=FALSE)+
  ggsn::scalebar(data = braz, location = "topleft", model = "WGS84",st.dist = 0.03,
                 dist = 5, dist_unit = "km",transform = TRUE) +
  xlim(c(bbox[c(1,3)])) +
  ylim(c(bbox[c(2,4)])) +
  ggtitle(paste("C. Fires in Brazil. July 14 - September 30, 2019")) +
  theme(legend.position ="none",
        legend.justification = c(0,0),
        legend.title = element_blank(),
        panel.border = element_rect(color="black", fill=NA))


locator_plot_b <- ggplot(brazil) +  
  theme(panel.border = element_rect(color="black", fill=NA),
        plot.background = element_rect(fill="white"))+
  geom_sf(fill="white") +
  geom_sf(data=locator_box, fill="transparent", color = "red", lwd=2)

inset_b <- ggplot(braz_clim,aes(x=datetime, y=value)) +
  geom_line() +
  geom_abline(aes(slope = 0,intercept = 1, lty=as.factor(lty)))+
  scale_linetype_manual(values = c(0,2))+
  facet_wrap(~variable, scales = "free_y", nrow = 2) +
  theme_bw() +
  xlab("Date")+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

braz_cow <- ggdraw() +
  draw_plot(main_plot_b) +
  draw_plot(locator_plot_b, .72,.65,.25,.25) 

# tubbs ==========================================================

tubbs <- st_read("data/tubbs.gpkg") %>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))

tubbs_p <- st_read("~/data/fire/mtbs/mtbs_perimeter_data/") %>%
  filter(Fire_Name == "TUBBS")

lct <- tubbs %>%
  mutate(lc=raster::extract(x=raster(tubbs_lc_file), y=.)) %>%
  mutate(lc_name = lut_lc[lc]) %>%
  group_by(lc_name) %>%
  summarise(percent_lc = round((n()/nrow(.))*100)) %>%
  ungroup() %>%
  st_set_geometry(NULL)

statst <- tubbs %>%
  st_set_geometry(NULL) %>%
  group_by(daynight) %>%
  summarise(percent_dn = round((n()/nrow(.))*100,2)) %>%
  ungroup() 

daterange <- tubbs$acq_date %>% as.Date() %>% range()

tubbs_area <- st_read("~/data/background/CUS/CUS.shp") %>%
  st_transform(crs=st_crs(tubbs)) %>%
  st_crop(tubbs %>% st_buffer(1))

usa <- st_read("~/data/background/CUS/CUS.shp") %>%
  st_transform(crs=st_crs(tubbs)) %>%
  filter(STUSPS == "CA")

bbox <- st_bbox(tubbs)%>% as.numeric()

locator_box <- st_bbox(tubbs_area) %>% st_as_sfc()

if(!file.exists("data/fishnett.RDS")){
  fishnet <- st_bbox(tubbs)%>% 
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01)%>%
    as("Spatial")%>%
    raster::extract(x=raster(tubbs_lc_file), y=., 
                    fun = function(x,...) getmode(x), 
                    method="simple") 
  saveRDS(fishnet, "data/fishnett.RDS")
}else{
  fishnet <- readRDS("data/fishnett.RDS")
}

if(!file.exists("data/fishnet_lct.RDS")){
  fishnet_lc <- st_bbox(tubbs) %>%
    st_as_sfc() %>%
    st_buffer(0.5) %>%
    st_make_grid(cellsize = 0.01) %>%
    st_as_sf() %>%
    mutate(lc= fishnet,
           classes = lut_lc[lc])
  saveRDS(fishnet_lc, "data/fishnet_lct.RDS")
}else{
  fishnet_lc<-readRDS("data/fishnet_lct.RDS")
}

main_plot_t <- ggplot() +
  # geom_sf(data = fishnet_lc, 
  #         aes(fill=classes),
  #         color = "transparent", alpha = 0.5)+
  geom_sf(data = tubbs_area, fill = "transparent")+
  geom_sf(data = tubbs_p, fill = "transparent",lwd=1) +
  geom_sf(data = tubbs, aes(color=daynight, alpha = daynight), size=3,show.legend = "point") +
  scale_alpha_manual(values = c(0.5,0.5))+
  scale_color_manual(values = c("blue","red"))+
  scale_fill_manual(values = lut_colors)+
  xlim(c(bbox[c(1,3)])) +
  ylim(c(bbox[c(2,4)])) +
  ggsn::scalebar(data = tubbs, location = "topleft", model = "WGS84",st.dist = 0.03,
                 dist = 2, dist_unit = "km",transform = TRUE) +
  guides(fill = FALSE) +
  ggtitle(paste("B. Tubbs Fire. Oct 9-15, 2017")) +
  theme(legend.position = "none",
        legend.justification = c(1,0),
        legend.title = element_blank(),
        panel.border = element_rect(color="black", fill=NA));main_plot_t

inset_t <- read_csv("data/tubbs-hourly.csv") %>%
  mutate(lty = ifelse(name == "VPD (kPa)",2,0)) %>%
  group_by(name) %>%
  mutate(label_y = max(value)*0.85) %>%
  ungroup()%>%
  ggplot(aes(x=datetime_utc, y=value)) +
  geom_line() +
  xlab("Date") +
  geom_hline(aes(yintercept = 1, lty=as.factor(lty)))+
  scale_linetype_manual(values=c(0,2)) +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

locator_plot_t <- ggplot(usa) +  
  theme(panel.border = element_rect(color="black", fill=NA),
        plot.background = element_rect(fill="white"))+
  geom_sf(fill="white") +
  geom_sf(data=st_centroid(locator_box), color = "red", size=2)

tubbs_cow<-ggdraw() +
  draw_plot(main_plot_t) +
  draw_plot(locator_plot_t, .75,.15,.2,.2)

# all together =================================================================
leg <- get_legend(snowy%>%
                    mutate(daynight = ifelse(daynight=="day", "Day", "Night")) %>%
                    ggplot() +
                    geom_point(aes(x=latitude, y=longitude, color=daynight))+
                    scale_color_manual(values= c("blue", "red"))+
                    theme(legend.background = element_rect(color="black"),
                          legend.margin = margin(1,1,1,1, "mm"),
                          legend.title = element_blank()))


insets <- ggarrange(inset_s, inset_t, inset_b, 
                    nrow=1, ncol=3, labels = c("D. Snowy", "E. Tubbs", "F. Brazil"),
                    label.x = .95, label.y = .95, hjust="right", vjust="top");insets

finalfig <- ggdraw(xlim = c(0,7.5), ylim = c(0, 11.75)) +
  draw_plot(snowy_cow, x = 0, y = 7.75, width = 7.5, height = 4) +
  draw_plot(leg, x=6.14,y=7.1, width=2,height=2) +
  draw_plot(tubbs_cow, x = 0, y = 3, width = 3, height = 4.75) +
  draw_plot(braz_cow, x = 3, y=3, width = 4.5, height = 4.75) +
  draw_plot(insets, x=0,y=0, width=7.5, height=3) +
  draw_label(label = "1 kPa", x = 2.18, y=0.7, size=12,fontface = "bold") +
  ggsave("images/panel_fig_no_lc.png", height = 11.75, width = 7.5)

