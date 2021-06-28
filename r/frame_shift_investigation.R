library(tidyverse)

read_csv("data/Temperate_Evergreen_Needleleaf_Forests.csv") %>%
  filter(nid == "442353")
read_csv("data/Temperate_Evergreen_Needleleaf_Forests_gamready.csv") %>%
  filter(nid == "442353")

read_csv("data/Temperate_Woody_Savannas.csv") %>%
  filter(nid == "474994") %>%
  ggplot(aes(x=rounded_datetime, y=n)) +
  geom_line() +
  scale_x_datetime(date_breaks = "12 hours") +
  theme(axis.text.x = element_text(angle = 90))

read_csv("data/Temperate_Woody_Savannas_gamready.csv") %>%
  filter(nid == "474994")
