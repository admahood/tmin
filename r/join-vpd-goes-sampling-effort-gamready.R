library(data.table)
library(stringr)
library(lubridate)

# joining vpd and goes counts ===================================
dir.create("data/out/gamready", recursive = TRUE, showWarnings = FALSE)

# # test using Arid Deciduous Broadleaf (~24 events)
# system2(command = "aws", args = "s3 cp s3://earthlab-amahood/night_fires/vpd_lc/Arid_Deciduous_Broadleaf_Forests_vpds.csv data/vpd_lc/Arid_Deciduous_Broadleaf_Forests_vpds.csv")
# system2(command = "aws", args = "s3 cp s3://earthlab-amahood/night_fires/goes_counts/Arid_Deciduous_Broadleaf_Forests.csv data/goes_counts/Arid_Deciduous_Broadleaf_Forests.csv")
# system2(command = "aws", args = "s3 cp s3://earthlab-amahood/night_fires/gamready/Arid_Deciduous_Broadleaf_Forests_gamready.csv data/out/Arid_Deciduous_Broadleaf_Forests_gamready.csv")
# 
# # bigger test using Boreal_Evergreen_Needleleaf_Forests (~176 events)
# system2(command = "aws", args = "s3 cp s3://earthlab-amahood/night_fires/vpd_lc/Boreal_Evergreen_Needleleaf_Forests_vpds.csv data/vpd_lc/Boreal_Evergreen_Needleleaf_Forests_vpds.csv")
# system2(command = "aws", args = "s3 cp s3://earthlab-amahood/night_fires/goes_counts/Boreal_Evergreen_Needleleaf_Forests.csv data/goes_counts/Boreal_Evergreen_Needleleaf_Forests.csv")
# system2(command = "aws", args = "s3 cp s3://earthlab-amahood/night_fires/gamready/Boreal_Evergreen_Needleleaf_Forests_gamready.csv data/out/Boreal_Evergreen_Needleleaf_Forests_gamready.csv")

# These are pretty big file syncs (50 GB) so don't run them if you're not doing the whole
# analysis. Instead, get the necessary vpd_lc and goes data using the `aws s3 cp` commands
# that are commented out above. (Two separate examples; one big and one small)
system("aws s3 sync s3://earthlab-amahood/night_fires/vpd_lc data/vpd_lc")
system("aws s3 sync s3://earthlab-amahood/night_fires/goes_counts data/goes_counts")
system(str_c("aws s3 sync ",
             file.path("s3://earthlab-amahood","night_fires","gamready")," ",
             file.path("data", "out", "gamready")))

vpd_files <- list.files("data/vpd_lc", pattern=".csv")
goes_files <- list.files("data/goes_counts", pattern = ".csv")

for(f in vpd_files){
  if(!file.exists(file.path("data", "out", "gamready", str_replace(f, "vpds", "gamready")))){
    if(file.exists(file.path("data", "goes_counts", str_replace(f, "_vpds", "")))){
      print(f)
      gc()
      
      # read GOES data first, because we only need VPD data for rows of VPD data frame
      # that are part of a FIRED event that has GOES detections
      goes_raw <- data.table::fread(file.path("data", "goes_counts", str_replace(f, "_vpds", "")))
      
      # If there were no GOES fire detections for this landcover type, then there
      # will be 0 rows in the dataframe, and the loop should iterate to the next file
      if(nrow(goes_raw) == 0) next() 
      
      # FIRST KEY STEP HERE for getting dataframe in analysis-ready format
      # Collapse the goes detections data.table so that we count the number of "at least 1 detection"
      # GOES images per fire event, per rounded datetime
      goes <- goes_raw[, .(fire_scenes_per_hour = .N), by = .(rounded_datetime, nid)]
      
      data.table::setkey(goes, nid, rounded_datetime)
      
      # Read VPD data and set the key
      vpds <- data.table::fread(file.path("data", "vpd_lc", f), key = c("nid"))
      
      # Subset VPD dataframe to just rows where there are GOES detections associated
      # with the FIRED event; this should speed up the creation of the rounded_datetime
      # column and might speed up the table join
      vpds <- vpds[.(unique(goes$nid)), nomatch = NULL]
      vpds[, `:=`(rounded_datetime = lubridate::ymd_h(paste(date, hour)),
                  fireID = NULL,
                  day = NULL,
                  hour = NULL,
                  first_date = NULL,
                  date = NULL)]
      
      # Set another key that is the rounded datetime for joining to GOES data
      data.table::setkey(vpds, nid, rounded_datetime)
      
      # left join the VPD dataframe with the GOES dataframe, then (using 
      # data.table chaining), assign all the GOES counts (column 'n') to 0 when
      # the VPD values for that 'nid' and 'rounded_datetime' don't have 
      # GOES detections
      vpds_goes <- goes[vpds][is.na(fire_scenes_per_hour), fire_scenes_per_hour := 0]
      
      # SECOND KEY STEP HERE for getting dataframe into analysis-ready format
      # We need to remove all the leading and trailing hours of each fire event that
      # have 0 GOES detections. Sometimes there are lots of rows like this, because
      # we extracted ERA5 climate data based on the start and end dates for each 
      # FIRED event which are based on the MODIS burned area product. So often the
      # MODIS burned area product has a burn duration that is much longer than
      # when GOES was detecting active fire.
      # This code originally lived in the count_ts.R script and only was run just
      # before modeling. We need to do this step earlier now, because we are going
      # to lose the rounded_datetime information (necessary for determining which 0's
      # occur before or after the string of GOES counts) when we collapse the binomial
      # distributions into sufficient statistics based on unique combinations of
      # nid and VPD_hPa
      
      vpds_goes_trimmed <-
        vpds_goes[, `:=`(first_detection = min(rounded_datetime[fire_scenes_per_hour > 0]),
                         hour_after_last_detection = max(rounded_datetime[fire_scenes_per_hour > 0]) + 60*60),
                  by = .(nid)][rounded_datetime >= first_detection & rounded_datetime <= hour_after_last_detection]
      
      data.table::setkey(vpds_goes_trimmed, rounded_datetime)
      
      # THIRD KEY STEP HERE for getting dataframe into analysis-ready format is to join with the sampling effort
      # dataframe (sampling_effort)
      
      sampling_effort <- data.table::fread("data/sampling-effort-goes16.csv", key = "rounded_datetime")
      
      vpds_goes_sampling_effort <- sampling_effort[vpds_goes_trimmed]
      
      # Sometimes there were no GOES scenes in a particular hour at all. These show up as NA after
      # joining the sampling effort dataframe to the VPD dataframe (because the VPD dataframe has
      # an entry for *every* hour, but there aren't rows in the sampling effort dataframe for
      # hours during which there were no GOES detections). We set these to 0 (they'll need to be
      # filtered out prior to modeling, but the model ought to do that automatically)
      vpds_goes_sampling_effort <- vpds_goes_sampling_effort[is.na(n_scenes_per_hour), n_scenes_per_hour := 0]
      
      # FOURTH KEY STEP HERE for getting dataframe into analysis-ready format is to further aggregate the binomial
      # fire_scenes_per_hour/n_scenes_per_hour by summing each component (successes and trials) across unique
      # combinations of nid and VPD_hPa
      
      vpds_goes_gam_ready <- vpds_goes_sampling_effort[, .(fire_scenes = sum(fire_scenes_per_hour),
                                                           n_scenes = sum(n_scenes_per_hour)),
                                                       by = .(nid, VPD_hPa)]
      
      # Write to disk
      data.table::fwrite(vpds_goes_gam_ready, file = file.path("data", "out", "gamready", str_replace(f, "vpds", "gamready")))
      
      system(str_c("aws s3 cp ",
                   file.path("data", "out", "gamready", str_replace(f, "vpds", "gamready")), " ",
                   file.path("s3://earthlab-amahood","night_fires","gamready",
                             str_replace(f, "vpds", "gamready"))))
    }
  }
}
