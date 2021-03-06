# script to join fired burned area perimeters, af product, 
# and john's analysis product (vpd trends thing)
# then do some basic stats

library(sf)
library(tidyverse)
library(lubridate)
library(vroom)

s3_path_af <- "s3://earthlab-mkoontz/mcd14ml_analysis-ready"
flnm_base_af <- "mcd14ml_c006_v03_"
years<-2000:2020 # maybe go 2001 to 2019 for official analysis (or redo fired)

s3_path_fired <- "s3://earthlab-amahood"
flnm_na<- "fired_NorthAmerica_s1_t5.gpkg"
flnm_sa <- "fired_SouthAmerica_s1_t5.gpkg"

system(paste("aws s3 cp",
             file.path(s3_path_fired, flnm_na),
             file.path("data", flnm_na)))
system(paste("aws s3 cp",
             file.path(s3_path_fired, flnm_sa),
             file.path("data", flnm_sa)))

# fired
fired_na <- st_read("data/fired_NorthAmerica_s1_t5.gpkg") %>%
  mutate(year = lubridate::year(first_date_7+7))

fired_sa <- st_read("data/fired_SouthAmerica_s1_t5.gpkg") %>%
  mutate(year = lubridate::year(first_date_7+7)) 

fired_wh <- bind_rows(fired_na, fired_sa)

result <- list()
counter<-1
for(y in years){
  outfile<-paste0("af_fired_wh_", y, ".gpkg")
  if(!file.exists(file.path("data",outfile))){
  source <- file.path(s3_path_af, paste0(flnm_base_af, y, ".csv"))
  target <- file.path("data", "mcd14", paste0(flnm_base_af, y, ".csv"))
  if(!file.exists(target)){
    system(paste("aws s3 cp",
                 source,
                 target
    ))
    }
  
  af_joined <-vroom(target) %>%
    filter(longitude < 0) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs=4326) %>%
    st_transform(crs = st_crs(fired_wh)) %>%
    st_join(x=fired_wh %>% filter(year == y),
            y=.) %>%
    filter(acq_date >= first_date_7 & acq_date <= last_date_7) %>%
    # st_set_geometry(NULL) %>%
    group_by(id) %>%
    mutate(dn_01 = ifelse(dn_detect == "day", 1,0)) %>%
    summarise(year = first(year),
              lc = first(na.omit(lc)),
              n = n(),
              day_detections = sum(dn_01),
              night_detections = n - day_detections,
              night_fraction = night_detections/n) %>%
    ungroup()%>%
    mutate(kop = str_sub(lc, 1,1),
           lc_raw = str_sub(lc,2,3))
  
  st_write(af_joined, file.path("data", outfile),delete_dsn = TRUE)
  
  system(paste(
    "aws s3 cp",
    file.path("data",outfile),
    file.path("s3://earthlab-mkoontz", "mcd14ml_analysis-ready", "joined_with_fired",outfile)
  ))
  }
  # result[[counter]] <- af_joined
  counter <- counter+1
}

afwh <- do.call(bind_rows, result) %>%
  mutate(kop = str_sub(lc, 1,1),
         lc_raw = str_sub(lc,2,3))
# active fire western hemisphere

save(afna, file = "data/afna.Rda")

# explore ========================
system(paste(
  "aws s3 sync",
  file.path("s3://earthlab-mkoontz", "mcd14ml_analysis-ready", "joined_with_fired"),
  file.path("data", "af_fired_joined")
))

joined<- list.files(file.path("data", "af_fired_joined"), full.names = T) %>%
  lapply(st_read) %>%
  bind_rows %>%
  st_set_geometry(NULL)

ggplot(joined, aes(x=n)) +
  geom_histogram()+
  scale_x_log10()

ggplot(joined %>% filter(n>10, year >2002, year < 2020), 
       aes(x=as.factor(year), y=night_fraction,
       color = as.factor(kop))) +
  geom_boxplot() +
  facet_wrap(~kop)

ggplot(joined %>% filter(n>5, year >2002, year < 2020) %>% na.omit, 
       aes(x=year,
           y=night_detections, 
           color = as.factor(lc_raw))) +
  geom_smooth() +
  facet_grid(lc_raw~kop, scales = "free")


ggplot(joined %>% filter(year >2002, year < 2020) %>% na.omit, 
       aes(x=year,
           y=log(night_detections+1))) +
  geom_smooth(method = "glm", 
              method.args=list(family = "quasipoisson"))

ggplot(joined %>% filter(year >2002, year < 2020) %>% na.omit, 
       aes(x=year,
           y=night_fraction)) +
  geom_smooth()

ggplot(joined %>% filter(n>10,year >2002, year < 2020) %>% na.omit, 
       aes(x=year %>% as.factor,
           y=night_fraction)) +
  geom_boxplot(outlier.shape=NA) +
  facet_grid(kop~lc_raw)

# model - binomial glm, trials = n   ====================================
library(effects)
library(lme4)
mod1<-glm(night_fraction ~ year*as.factor(lc), 
          weights = n,
          data = joined, family = "binomial")
summary(mod1)

plot(allEffects(mod1))


mod2<-glm(night_fraction ~ year*as.factor(kop), 
            weights = n,
            data = joined, 
            family = "binomial")
summary(mod2)

mod3<-glm(night_fraction ~ year*as.factor(lc_raw), 
          weights = n,
          data = joined %>% filter(kop ==3), 
          family = "binomial")
summary(mod3)

# theil sen =======
library(mblm)
library(foreach)
library(doParallel)



lcs <- joined$lc %>% na.omit %>% unique

registerDoParallel(cores=detectCores()-1)
tslist <- foreach(i = 1:length(lcs), .combine=bind_rows)%dopar%{
# for(i in 1:length(lcs)){
  system(paste("echo", i))
  d <- filter(joined, lc == lcs[i], n>10)
  # if(nrow(d) > 100) d <- slice(d, 1:1000) # for testing
  
  if(nrow(d)>10){
    tsmod <- mblm(night_fraction ~ year,
                  repeated = TRUE,
                  dataframe = d)
  }
  xx<- summary(tsmod)
  tibble(lckop = lcs[i], 
         intercept = tsmod$coefficients[1], 
         slope = tsmod$coefficients[2],
         p = xx$coefficients[2,4])
}
tslist
