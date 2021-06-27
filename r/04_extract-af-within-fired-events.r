# extract GOES16 data to fired polygons

libs<- c("tidyverse", "sf", "lubridate", "data.table", "vroom",
         "doParallel", "foreach")
invisible(sapply(libs, library, character.only=TRUE, quietly=TRUE))

# data import ==================================================================

# goes
system("aws s3 sync s3://earthlab-mkoontz/goes16 data/goes16 --only-show-errors")
# system("aws s3 sync s3://earthlab-amahood/night_fires/gamready data/gamready --only-show-errors")

# fired polygons
system("aws s3 sync s3://earthlab-amahood/night_fires/lc_splits data/fired --only-show-errors")

# effort (number of goes scenes per hour)
system("aws s3 cp s3://earthlab-mkoontz/goes16meta/sampling-effort-goes16.csv data/segoes.csv")


# The Business =================================================================

# extracting the needed info from the goes extract file names
goes_files<- list.files("data/goes16", pattern = ".csv", full.names = TRUE) %>%
  as_tibble() %>%
  dplyr::rename(filename = value) %>%
  separate(filename, c("datetime", "datetime_mid"), sep = "_", remove=FALSE) %>%
  mutate(date = str_sub(datetime,13,20) %>% as.Date(date, format="%Y%m%d"))
fired_files <- list.files("data/fired", pattern = ".gpkg", full.names = TRUE)
dir.create("data/out")

effort <- vroom("data/segoes.csv") %>%
  dplyr::select(n_scenes = n_scenes_per_hour, rounded_datetime)

# extracting the detection counts to each fire perimeter in a nested for loop 
# (the interior loop is parallel)


for(i in 1:length(fired_files)){
  
  t0 <- Sys.time()
  fired <- st_read(fired_files[i])
  fired_crs <- st_crs(fired)
  out_file <- fired_files[i] %>%
    str_replace(".gpkg", ".csv") %>%
    str_replace("data/fired/", "")
  
  corz<- detectCores()-1
  registerDoParallel(corz)
  fc<-foreach(f = 1:nrow(fired), .combine= bind_rows)%dopar%{
    # fc<-foreach(f = 1:150, .combine= bind_rows)%dopar%{
      
    
    goes <- goes_files %>%
      dplyr::filter(date >= fired[f,]$first_date_7)%>%
      dplyr::filter(date <= fired[f,]$last_date_7) %>%
      dplyr::select(filename) %>%
      pull()
    
    if(length(goes)>0){
      tbl <- goes %>% 
        map_df(~read_csv(., col_types = c("TTdddddddddd")))%>%
        dplyr::select(rounded_datetime, Mask, sinu_x, sinu_y, 
                      exact_time = scan_center,
                      goes_cellindex=cellindex) %>%
        na.omit()%>%
        st_as_sf(coords=c("sinu_x", "sinu_y"), crs=fired_crs)
      
      ints<- tbl%>%
        mutate(is_fire = st_intersects(., fired[f,],sparse = F) %>%
                 rowSums()) %>%
        filter(is_fire > 0) 
      
      if(nrow(ints)>0){
        system(paste("echo", f, round(f/nrow(fired)*100,2), "%", out_file))
        fire_counts <- ints %>%
          st_set_geometry(NULL) %>%
          group_by(rounded_datetime, exact_time) %>%
          dplyr::summarise(n=n()) %>%
          mutate(nid = fired[f,]$nid)%>%
          ungroup() %>%
          left_join(effort)
        return(fire_counts)
        
      }
    }
  }
  print(Sys.time()-t0)
  
  if(!is.null(fc)){  
  write_csv(fc,file.path("data","out", out_file))
  system(paste0("aws s3 cp ", file.path("data","out", out_file), 
                " s3://earthlab-amahood/night_fires/goes_counts/", out_file))}
}
