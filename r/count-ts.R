library(tidyverse)
library(mgcv)
library(vroom)
library(sf)
library(parallel)
library(lubridate)
library(pbapply)
library(lutz)
library(units)
library(Matrix)

theme_set(theme_minimal() + 
            theme(panel.grid.minor = element_blank()))

system("aws s3 sync s3://earthlab-amahood/night_fires/gamready data/gamready")
system("aws s3 cp s3://earthlab-amahood/night_fires/lut_ba.Rda data/lut_ba.Rda")
system("aws s3 cp s3://earthlab-mkoontz/goes16meta/sampling-effort-goes16.csv data/s_effort.csv")

# fired <- st_read("data/fired/events_w_attributes.gpkg") %>%
#   st_transform(4326)
# 
# fire_events <- fired %>%
#   filter(ignition_date > as.Date('2014-01-07')) %>%
#   mutate(duration = last_date - ignition_date)
# 
# fire_areas <- fire_events %>%
#   as.data.frame %>%
#   select(-geom) %>%
#   distinct(id, total_area_km2, lc_name, lc) %>%
#   as_tibble
# 
# fire_tz <- fired %>%
#   st_transform(crs = 4326) %>%
#   st_centroid() %>%
#   tz_lookup(method = 'accurate') %>%
#   tibble %>%
#   mutate(id = fired$id) %>%
#   rename(timezone = ".")
# 
# # system("aws s3 cp s3://earthlab-mbjoseph/goes-survival/out/goes-af-era5.csv out/goes-af-era5.csv")
# 
# d <- vroom('out/goes-af-era5.csv') %>%
#   left_join(fire_areas) %>%
#   left_join(fire_tz)
# 
# event_summaries <- d %>%
#   filter(n > 0) %>%
#   group_by(id, total_area_km2) %>%
#   summarize(n_days = length(unique(day_of_year)), 
#             total_n = sum(n)) %>%
#   arrange(total_n)

gamready_files <- list.files("data/gamready", full.names = TRUE, pattern = ".csv$")
load("data/lut_ba.Rda")


lut_ba <- drop_units(lut_ba)
hourdf<- vroom("data/s_effort.csv") %>%
  mutate(rounded_datetime = ymd_hm(rounded_hour)) %>%
  dplyr::rename(n_scenes=n)

# lut_n<-hourdf$rounded_hour
# names(lut_n)<-hourdf$n
# lut_ba
for(f in gamready_files){}
# subset data to include time of first detection to first terminal zero

out_fn_base <- str_replace(f, "data/gamready/","") %>%
  str_replace("_gamready.csv","")

# d<- vroom(f) %>%
#   mutate(rounded_hour = str_replace_all(str_sub(rounded_datetime, 1,13), "-", "") %>%
#            str_replace_all(" ", "")%>% as.numeric)

events <- vroom(f) %>%
  mutate(nid = factor(as.character(nid))) %>%
  group_by(nid) %>%
  mutate(first_detection = min(rounded_datetime[n > 0]), 
         hour_after_last_detection = max(rounded_datetime[n > 0]) + 60*60) %>%
  filter(rounded_datetime >= first_detection, 
         rounded_datetime <= hour_after_last_detection) %>%
  ungroup %>%
  dplyr::select(-hour_after_last_detection, first_detection) %>%
  mutate(ba = lut_ba[nid]) %>%
  left_join(hourdf, by="rounded_datetime") %>%
  mutate(effort = log(ba) * n_scenes)

# write_csv(events, "out/events.csv")

# length(unique(events$id))
# nrow(events)
# 
# events %>%
#   group_by(lc_name) %>%
#   summarize(n_events = length(unique(fid))) %>%
#   arrange(n_events)



# HGAM fitting
# split_events <- events #%>%
  #split(.$lc_name)

# lc_counts <- events %>%
#   group_by(lc_name) %>%
#   distinct(id) %>%
#   ungroup %>%
#   count(lc_name)%>%
#   arrange(n) %>%
#   filter(n >= 15)
# 
# events %>%
#   filter(lc_name %in% lc_counts$lc_name) %>%
#   nrow

# for (i in seq_along(split_events)) {
#   print(i)
#   df <- split_events[[i]] 
  if (nrow(events) == 0) next
  outname <- paste0("data/mods/",out_fn_base, 
                    "-gam.rds")
  if (file.exists(outname)) next
  events <- droplevels(events)
  if (length(unique(events$nid)) < 15) next
  levels(events$nid) <- c(levels(events$nid), "NewFire")
  m <- bam(n  ~ s(VPD_hPa) + s(VPD_hPa, nid, bs = "fs"), 
           data = events, 
           nthreads = parallel::detectCores(), 
           offset = effort, 
           family = nb(), 
           discrete = TRUE,
           drop.unused.levels = FALSE)
  write_rds(m, outname)
# }


# For each land cover type, simulate from the predictive distribution
cover_types <- str_replace_all(gamready_files, "data/gamready/","") %>%
  str_replace("_gamready.csv","")

predictive_sim <- function(cover_type) {
  
  # subd <- vroom(str_c("data/gamready/",cover_type))
  model_path <- paste0(out_fn_base, 
                       "-gam.rds")
  if (!file.exists(model_path)) {
    return(NA)
  }
  prediction_csv_path <- gsub("-gam.rds", "-preds.csv", model_path)
  if (file.exists(prediction_csv_path)) {
    return(NA)
  }
  m <- model_path %>%
    read_rds
  model_summary <- summary(m)
  write_rds(model_summary, gsub("-gam", "-gam-summary", model_path))
  subd %>%
    mutate(n_pred = fitted(m)) %>%
    distinct(vpd_kPa, n_pred, total_area_km2, fid, lc_name) %>%
    ggplot(aes(vpd_kPa, n_pred / total_area_km2, group = fid)) + 
    geom_path(alpha = .05) +
    scale_y_log10() + 
    facet_wrap(~lc_name)
  
  # Simulate from predictive distribution
  pred_df <- subd %>%
    distinct(lc_name) %>%
    mutate(total_area_km2 = events %>% 
             distinct(fid, total_area_km2) %>% 
             summarize(median(total_area_km2)) %>% 
             unlist %>% 
             c,
           vpd_kPa = list(seq(0, max(events$vpd_kPa), by = .05)), 
           fid = "NewFire") %>%
    unnest(cols = c(vpd_kPa)) %>%
    mutate(idx = 1:n())
  Xp <- predict(m, pred_df, type = "lpmatrix")
  nonzero_cols <- colSums(Xp) != 0
  beta <- coef(m)[nonzero_cols]
  Vb <- vcov(m)[nonzero_cols, nonzero_cols] ## posterior mean and cov of coefs
  n <- 1000
  br <- MASS::mvrnorm(n, beta, Vb) 
  dfs <- list()
  pb <- txtProgressBar(max = n, style = 3)
  stopifnot(length(unique(pred_df$total_area_km2)) == 1)
  for (j in 1:n) {
    dfs[[j]] <- tibble(Xb = Xp[, nonzero_cols] %*% br[j, ] %>% 
                         as.vector, 
                       j = j) %>%
      mutate(idx = 1:n(), 
             log_mu = Xb + log(pred_df$total_area_km2[1]))
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  posterior_predictions <- dfs %>%
    bind_rows %>%
    left_join(pred_df) %>%
    arrange(lc_name, j, vpd_kPa)  %>%
    mutate(y_pred = MASS::rnegbin(n(), 
                                  exp(log_mu), 
                                  theta = exp(m$family$getTheta())), 
           pr_zero = dnbinom(0, size = exp(m$family$getTheta()), mu = exp(log_mu)))
  posterior_predictions %>%
    write_csv(prediction_csv_path)
}

sims <- pblapply(cover_types, predictive_sim)