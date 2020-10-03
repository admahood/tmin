# preprocessing modis/viirs active fire

# preprocessing ================================================================

read_csv("data/fire_nrt_V1_101893.csv")%>%
  dplyr::mutate(acq_hour = as.numeric(substr(acq_time, start = 1, stop = 2)),
                acq_min = as.numeric(substr(acq_time, start = 3, stop = 4)),
                acq_datetime = ymd_hm(paste0(acq_date, " ", acq_hour, ":", acq_min)),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  write_csv("data/fire_nrt_australia.csv")

# preprocessing brazil =========================================================

read_csv("data/fire_nrt_V1_104769.csv")%>%
  dplyr::mutate(acq_hour = as.numeric(substr(acq_time, start = 1, stop = 2)),
                acq_min = as.numeric(substr(acq_time, start = 3, stop = 4)),
                acq_datetime = ymd_hm(paste0(acq_date, " ", acq_hour, ":", acq_min)),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  write_csv("data/fire_nrt_brazil.csv")

# preprocessing tubbs ==============


read_csv("data/fire_archive_V1_104769.csv")%>%
  dplyr::mutate(acq_hour = as.numeric(substr(acq_time, start = 1, stop = 2)),
                acq_min = as.numeric(substr(acq_time, start = 3, stop = 4)),
                acq_datetime = ymd_hm(paste0(acq_date, " ", acq_hour, ":", acq_min)),
                acq_year = year(acq_datetime),
                acq_month = month(acq_datetime),
                acq_day = day(acq_datetime),
                solar_offset = longitude / 15,
                hemisphere = ifelse(latitude >= 0, yes = "Northern hemisphere", no = "Southern hemisphere"),
                acq_datetime_local = acq_datetime + as.duration(solar_offset * 60 * 60),
                local_doy = lubridate::yday(acq_datetime_local),
                local_hour_decmin = ((acq_hour) + (acq_min / 60) + solar_offset + 24) %% 24,
                local_solar_hour_decmin_round = round(local_hour_decmin),
                local_solar_hour_decmin_round0.5 = round(local_hour_decmin * 2) / 2,
                h = (local_hour_decmin - 12) * 15 * pi / 180,
                phi = latitude * pi / 180,
                delta = -asin(0.39779 * cos(pi / 180 * (0.98565 * (local_doy + 10) + 360 / pi * 0.0167 * sin(pi / 180 * (0.98565 * (local_doy - 2)))))),
                solar_elev_ang = (asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(h))) * 180 / pi,
                daynight = ifelse(solar_elev_ang > 0, yes = "day", no = "night")) %>%
  write_csv("data/fire_archive_brazil.csv")
