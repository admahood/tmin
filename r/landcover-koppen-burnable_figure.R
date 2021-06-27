# setup ========================================================================
library(tidyverse)
library(stars)
library(raster)
library(s2)
library(foreach)
library(ggpubr)
library(doParallel)

# koppen look up tables=========================================================
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
           "Urban and Built-up",
           "Cropland Natural Vegetation Mosaics",
           "Permanent Snow and Ice",
           "Barren",
           "Water Bodies")
names(lut_lc) <- str_pad(1:17, width = 2, side = "left",pad = "0")

# color palette ================================================================
df_colors <- read_csv("data/landcover-colors-2021.csv") %>%
  mutate(lc_name= replace(lc_name, lc_name == "Urban and Built-up Lands",
                          "Urban and Built-up"))
lc_cols <- pull(df_colors, color)
names(lc_cols) <- pull(df_colors, lc_name)

# reading lck in ===============================================================
lck<-read_stars("in/lc_koppen_2010_mode.tif")
lck_s <- read_stars("out/aggregations_2003-2020/lck_shifted.tif")
d_ovp <- read_stars("out/aggregations_2003-2020/2003-2020_day_overpass-count.tif")
n_ovp <- read_stars("out/aggregations_2003-2020/2003-2020_night_overpass-count.tif")

# calculating area =============================================================
lck_poly<-lck_s %>% 
  stars::st_xy2sfc(as_points = FALSE) %>%
  st_as_sf() %>%
  as_tibble() %>% 
  dplyr::select(lck = lck_shifted.tif, geometry)  %>% 
  mutate(geometry = as_s2_geography(geometry)) %>%
  mutate(area_km2 = s2_area(geometry)/1000000)

# summarising area by landcover/koppen =========================================
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

# making a table with the thresholds ===========================================
thresholds <- read_csv("in/csvs_from_michael/zero-goes-af-vpd-thresholds-with-landcover-codes.csv") %>%
  dplyr::select(lc_kop = lc_name, vpd_thresh_hpa, sd) %>%
  left_join(lck_tab %>% dplyr::select(lc_kop, area_km2)) %>%
  mutate(millions_of_km2 = area_km2/1000000) %>%
  arrange(lc_kop)

write_csv(thresholds %>% dplyr::select(-area_km2),"out/lck_thresh_area.csv")

# making a table of burnable koppen area =======================================

burnable_koppen <-
  thresholds %>%
  separate(lc_kop, into = "koppen") %>%
  group_by(koppen) %>%
  summarise(area_Mkm2 = sum(millions_of_km2)) %>%
  ungroup()

write.csv(burnable_koppen, "out/burnable_koppen.csv")

# making a 3 panel figure - lc, kop, burnable globe mask =====================


lck_s <- terra::rast("out/aggregations_2003-2020/lck_shifted.tif")
lck_s[lck_s==100] <- NA
lck_s[lck_s==200] <- NA
lck_s[lck_s==300] <- NA
lck_s[lck_s==400] <- NA
lck_s[lck_s==500] <- NA

brnbl<-thresholds %>% 
  dplyr::select(lc_kop) %>% 
  mutate(Burnable = "Burnable")

lck_df<-as.data.frame(lck_s, xy=TRUE)%>%
  mutate(Koppen = lut_kop[str_sub(lck_shifted,1,1)],
         Landcover = lut_lc[str_sub(lck_shifted,2,3)],
         lc_kop = paste(Koppen, Landcover)) %>%
  left_join(brnbl) %>%
  replace_na(list(Burnable = "Not Burnable")) %>%
  mutate(x = ifelse(x< 13, x+360, x))
  
# theme clean? or plot.border or something like that
p1 <- ggplot(lck_df) +
  geom_raster(aes(x=x,y=y,fill=Koppen)) +
  coord_equal() +
  theme_void()+
  labs(caption="Köppen-Geiger Climate Classification") + 
  theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
  # ggtitle("Koppen-Gieger Climate Classification") +
  theme(legend.justification = c(0,0),
        legend.position = c(0.05,0.2),
        text = element_text(family="DejaVuSans"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", size=1)) +
  ylim(c(-52.625, 75.125))

p2 <- ggplot(lck_df%>%
               mutate(Landcover, replace(Landcover, 
                                       Landcover == "Urban and Built-up Lands",
                                       "Urban and Built-up"))) +
  geom_raster(aes(x=x,y=y,fill=Landcover))  +
  coord_equal() +
  theme_void()+
  labs(caption="MOD12Q1 Landcover Classification") + 
  theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
  # ggtitle("MOD12Q1 Landcover Classification") +
  theme(legend.title = element_blank(),
        text = element_text(family="DejaVuSans"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", size=1))+
  scale_fill_manual(values = lc_cols)+
  ylim(c(-52.625, 75.125))

p3 <- ggplot(lck_df) +
  geom_raster(aes(x=x,y=y,fill=Burnable)) +
  coord_equal() +
  theme_void()+
  # ggtitle("Burnable Land Area") +
  labs(caption="Burnable Land Area (> 100 fires 2017 - 2020)") + 
  theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))+
  theme(legend.justification = c(0,0),
        legend.position = c(0.05,0.2),
        text = element_text(family="DejaVuSans"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = "black", size=1)) +
  scale_fill_manual(values = c("firebrick", "grey"))+
  ylim(c(-52.625, 75.125))

leg2 <- get_legend(p2)


plots <- ggarrange(p1,
                   p2 + theme(legend.position = "none"), leg2, p3,
                   ncol=1, nrow=4,
                   heights = c(2,2,0.75, 2),
                   labels = c("a", "b", "","c"),label.y = 0.95,
                   font.label = "bold") +
# ggarrange(plots, leg2, ncol=2, widths = c(5,2))+
  ggsave(filename = "out/three_panel_lc_map.png", width = 8.5, height = 12)


# grabbing some quick numbers for the paper
# earth has 148940000 km2 of land surface area (from wikipedia)
total_land_area = lck_tab %>% pull(area_km2) %>% sum #145M km2 -- close enough
burnable_land_area = thresholds %>% pull(area_km2) %>% sum() # 90M km2
burnable_land_area/(total_land_area) # this is where I got the 62% for the paper

# getting extended table 1 stats right here baby (michael's script has a lot of 
# adjusting that is no longer necessary)

