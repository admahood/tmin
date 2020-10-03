# functions for tmin

get_af_event<- function(polygon, p_crs, af_product, tmin_path, r_native,
                        meat = "PRISM_tmin_stable_4kmD1_", bits ="_bil.bil",
                        start, end, timezone = tz){
  
  xx <- af_product %>%
    filter(ACQ_DATE>=as.Date(start), ACQ_DATE<=as.Date(end)) 
  if(nrow(xx)>0){
    xx <- xx %>%
      st_transform(p_crs)%>%
      st_intersection(polygon) %>%
      st_intersection(st_transform(tz, crs=p_crs)) %>%
      mutate(ACQ_TIME = as.numeric(as.character(ACQ_TIME))%>% str_pad(4, pad="0"),#,
             dt = as.POSIXct(strptime(paste(ACQ_DATE, ACQ_TIME), "%Y-%m-%d %H%M")),
             lt = dt+total_seconds_offset,
             acq_date_l = as.Date(substr(lt, 1,10)),
             acq_hour_l = as.numeric(substr(lt,12,13)),
             daynight = ifelse(acq_hour_l<7, "_night", "day"),
             daynight = ifelse(acq_hour_l>=21, "_night", daynight),
             daynight = as.character(daynight)
            ) %>%
      mutate(event_day = (acq_date_l - min(acq_date_l)),
             event_day = ifelse(acq_hour_l>=21, event_day+1, event_day)
             )%>%
      group_by(event_day, daynight) %>%
      mutate(mean_FRP = mean(FRP)) %>%
      ungroup()
    
    if(nrow(xx)>0){
      if(min(xx$event_day) == 0) xx$event_day <- xx$event_day +1
      if(any(names(xx) == "DAYNIGHT"))  xx <- xx %>% dplyr::select(-DAYNIGHT)
      
      # making sure that the day that tmin is extracted from is correct
      hp_day0 <- max(xx[xx$event_day ==1,]$acq_date_l)-1
      
      r_native <- r_native[st_transform(polygon, crs=st_crs(r_native))]
      r_prism <- r_native %>%
        st_bbox %>%
        st_as_sfc() %>%
        st_transform(prism_crs)%>%
        st_bbox()%>%
        st_as_sfc()
      
      ll <- list()
      for(d in 1:length(unique(xx$event_day))){
        ed <- as.character(hp_day0+unique(xx$event_day)[d])
        year <- substr(ed,1,4)
        month<- substr(ed,6,7)
        day <- substr(ed, 9,10)
        
        rf <- file.path(tmin_path, paste0(meat, year, month, day, bits)) 
        r <- raster(rf) %>%
          crop(as(r_prism, "Spatial")) %>%
          projectRaster(as(r_native, "Raster")) 
        
        ll[[d]] <- filter(xx, event_day == unique(xx$event_day)[d])%>%
          mutate(tmin = raster::extract(r,.)%>% as.numeric,
                 tmin_date = paste0(year,month,day),
                 lc = raster::extract(as(r_native,"Raster"),.)) # just to doublecheck
      }
      
      return(do.call("rbind",ll))
    }
  }
}
