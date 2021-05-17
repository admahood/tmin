# script to join fired burned area perimeters, af product, 
# and john's analysis product (vpd trends thing)
# then do some basic stats

library(sf)
library(tidyverse)

# na first

fired_na <- st_read("data/fired_NorthAmerica_s1_t5.gpkg")


# fired_with_afd <-
#   sf::st_join(fired_ca_event, modis_afd) %>%
#   filter(ACQ_DATE >= ignition_date & ACQ_DATE <= last_date)