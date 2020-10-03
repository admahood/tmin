library(sf);library(tidyverse)

sa <- st_read("data/fired_SouthAmerica_s1_t5.gpkg")
na <- st_read("data/fired_NorthAmerica_s1_t5.gpkg")

conts <- st_read("/home/a/data/background/world_borders/ne_50m_admin_0_countries.shp") %>%
  group_by(CONTINENT) %>%
  summarise()
  
sa_cont <- filter(conts, CONTINENT == "South America") %>%
  st_transform(crs = st_crs(sa))
na_cont <- filter(conts, CONTINENT == "North America") %>%
  st_transform(crs = st_crs(na))

# south america =======
sa_goes <- sa %>%
  filter(first_date_7 > as.Date("2016-12-31"))

x1<-st_intersects(sa_goes, sa_cont,sparse = F) %>%
  rowSums()

sa_goes <- sa_goes %>%
  mutate(is_sa = x1)%>%
  filter(is_sa > 0)
st_write(sa_goes, "data/fired_sa_2017-.gpkg")


# north america ================
na_goes <- na %>%
  filter(first_date_7 > as.Date("2016-12-31"))

x1<-st_intersects(na_goes, na_cont,sparse = F) %>%
  rowSums()

na_goes <- na_goes %>%
  mutate(is_na = x1)%>%
  filter(is_na > 0)
st_write(na_goes, "data/fired_na_2017-.gpkg", delete_dsn = TRUE)

