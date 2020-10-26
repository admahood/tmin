# extract GOES16 data to fired polygons

libs<- c("tidyverse", "exactextractr", "sf", "lubridate", "data.table", 
         "doParallel", "foreach")
invisible(sapply(libs, library, character.only=TRUE, quietly=TRUE))

# data import ==================================================================

# goes
system("aws s3 sync s3://earthlab-mkoontz/goes16 data/goes16 --only-show-errors")

# fired polygons
system("aws s3 sync s3://earthlab-amahood/night_fires/lc_splits data/fired --only-show-errors")

goes_files<- list.files("data/goes16", pattern = ".csv") %>%
  as_tibble() %>%
  dplyr::rename(filename = value) %>%
  separate(filename, c("datetime", "datetime_mid"), sep = "_", remove=FALSE) %>%
  mutate(date = str_sub(datetime,1,8) %>% as.Date(date, format="%Y%m%d"),
         hour = str_sub(datetime, 9, 10))
fired_files <- list.files("data/fired", pattern = ".gpkg", full.names = TRUE)
dir.create("data/out")

corz <- detectCores()-1
registerDoParallel(corz)

foreach(i = fired_files){
  
  fired<- st_read(fired_files[i])
  
  out_file <- fired_files[i] %>%
    str_replace(".gpkg", ".csv") %>%
    str_replace("data/fired/", "")
  
  fire_counts <- list()
  for(f in 1:nrow(fired)){
    goes <- goes_files %>%
      dplyr::filter(date >= fired[f,]$first_date_7)%>%
      dplyr::filter(date <= fired[f,]$last_date_7)
    
    glist<- list()
    for(g in 1:nrow(goes)) {
      
      glist[[g]] <- read_csv(file.path("data","goes16",goes$filename[g]), col_types = c("TTdddddddddd"))%>%
          dplyr::select(-cellindex, -x,-y,-Area,-Temp,-Power, -DQF,-scan_center)# %>%
          # mutate(date = goes[g,]$date,
          #        hour = goes[g,]$hour)
        
    }
    big_thing<-bind_rows(glist) %>%
      st_as_sf(coords=c("sinu_x", "sinu_y"), crs=st_crs(fired))
    
    ecos3 <- big_thing %>%
      mutate(is_fire = st_intersects(big_thing, fired[f,],sparse = F) %>%
               rowSums())%>%
      filter(is_fire > 0) 
    
    fire_counts[[f]]%>%
      st_set_geometry(NULL) %>%
      group_by(rounded_datetime) %>%
      dplyr::summarise(n=n()) %>%
      mutate(nid = fired[f,]$nid)
  }
  bind_rows(fire_counts) %>%
    write_csv(file.path("data","out", out_file))
  system(paste0("aws s3 cp ", file.path("data","out", out_file), 
                " s3://earthlab-amahood/night_fires/goes_counts/", out_file
                ))
}