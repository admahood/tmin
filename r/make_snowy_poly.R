library(sf)
library(tidyverse)

snowy <- st_read("data/snowy_complex.gpkg") %>%
  mutate(daynight = factor(daynight, levels = c("night", "day")))%>%
  mutate(acq_date = as.Date(acq_date)) %>%
  filter(acq_date > as.Date("2019-12-28"),
         acq_date < as.Date("2020-01-6"))%>%
  st_buffer(dist = 0.1)%>%
  group_by(version) %>%
  summarise(max_date = max(acq_date),
            min_date = min(acq_date))%>%
  ungroup()

plot(snowy[0])

st_write(snowy, "data/snowy_poly_buff.gpkg")
