# script to join fired burned area perimeters, af product, 
# and john's analysis product (vpd trends thing)
# then do some basic stats

library(sf)
library(tidyverse)
library(lubridate)
library(vroom)

s3_path <- "s3://earthlab-mkoontz/mcd14ml_analysis-ready"
flnm_base <- "mcd14ml_c006_v03_"
years<-2000:2020 # maybe go 2001 to 2019 for official analysis (or redo fired)

# na fired
fired_na <- st_read("data/fired_NorthAmerica_s1_t5.gpkg") %>%
  mutate(year = lubridate::year(first_date_7+7))

result <- list()
counter<-1
for(y in years){
  source <- file.path(s3_path, paste0(flnm_base, y, ".csv"))
  target <- file.path("data", "mcd14", paste0(flnm_base, y, ".csv"))
  if(!file.exists(target)){
    system(paste("aws s3 cp",
                 source,
                 target
    ))
    }
  
  result[[counter]]<-vroom(target) %>%
    filter(longitude < 0) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs=4326) %>%
    st_transform(crs = st_crs(fired_na)) %>%
    st_join(x=fired_na %>% filter(year == y),
            y=.) %>%
    filter(acq_date >= first_date_7 & acq_date <= last_date_7) %>%
    st_set_geometry(NULL) %>%
    group_by(id) %>%
    mutate(dn_01 = ifelse(dn_detect == "day", 1,0)) %>%
    summarise(year = first(year),
              lc = min(na.rm=T),
              n = n(),
              day_detections = sum(dn_01),
              night_detections = n - day_detections,
              night_fraction = night_detections/n) %>%
    ungroup()
  
  counter <- counter+1
}

afna <- do.call(bind_rows, result) %>%
  mutate(kop = str_sub(lc, 1,1),
         lc_raw = str_sub(lc,2,3))
# active fire western hemisphere

save(afna, file = "data/afna.Rda")

# explore ========================
ggplot(afna, aes(x=n)) +
  geom_histogram()+
  scale_x_log10()

ggplot(afna %>% filter(n>5, year >2002, year < 2020), aes(x=as.factor(year), y=night_fraction,
       color = as.factor(kop))) +
  geom_boxplot() +
  facet_wrap(~kop)

ggplot(afna %>% filter(n>5, year >2002, year < 2020), aes(x=year,
                                                           y=night_detections, 
                                                           color = as.factor(kop))) +
  geom_smooth(method = "lm")

# model - binomial glm, trials = n   ====================================

