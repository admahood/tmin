# event day vs tmin
library(lmerTest)
library(lme4)
xx <- st_set_geometry(dplyr::select(hp_afv, event_day, tmin),NULL) %>% mutate(product = "viirs", fire = "hot_pot")%>%
  rbind(st_set_geometry(dplyr::select(hp_afm, event_day, tmin),NULL) %>% mutate(product = "modis", fire = 'hot_pot'))%>%
   rbind(st_set_geometry(dplyr::select(camp_afm, event_day, tmin),NULL) %>% mutate(product = "modis", fire = "camp"))%>%
  rbind(st_set_geometry(dplyr::select(camp_afv, event_day, tmin),NULL) %>% mutate(product = "viirs", fire = "camp"))

ggplot(xx, aes(x=as.factor(event_day), y=tmin)) +
  geom_boxplot() +
  facet_wrap(~product+fire,dir = "v")

m=lme4::lmer(tmin~event_day + (1|fire), xx)
summary(m)
