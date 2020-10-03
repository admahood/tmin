
library(sf); library(tidyverse); library(purrr)

system("aws s3 sync s3://earthlab-amahood/night_fires data")

path<- "data/sa_vpd.csv"
output_fn <- "data/sa_vpd_long_2017-test.csv"
# x is a large character string, lets subset it to be more usable (but still big)

good_nums <- st_read("data/fired_sa_2017-.gpkg") %>%
  st_set_geometry(NULL) %>%
  dplyr::select(id) %>%
  pull() %>%
  unique()


ragged_reader<- function(x){
  
  x <-read_lines(path, progress = TRUE) %>%
    purrr::map(.x = ., .f = function(x) {
      char_list <- stringr::str_split(string = x, pattern = ",")
      char_vec <- unlist(char_list)
      num_vec <- as.numeric(char_vec)
      return(num_vec[!is.na(num_vec)])
    })
  
  max_col_nums <- purrr::map_int(x, .f = length) %>% max()
  
  # total number
  max_days <- (max_col_nums - 1) %/% 24
  
  # new column names
  datetime_combos <- as.vector(t(base::outer(X = 0:max_days, Y = 0:23, FUN = paste, sep = "_")))
  datetime_combos <- datetime_combos[1:(max_col_nums - 1)]
  new_col_names <- c("fireID", datetime_combos)
  
  
  d<-do.call(rbind, x) %>%
    as_tibble() %>%
    `colnames<-`(new_col_names)
  return(d)
}
memory.limit()
memory.limit(size=56000)

d<-ragged_reader(path) #%>%
  
d %>% 
  filter(fireID %in% good_nums) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "hour",
               values_to = "VPD_hPa")%>%
  #separate(day_hour, c("day", "hour"), sep="_") # this adds some serious time
  write_csv(output_fn)
