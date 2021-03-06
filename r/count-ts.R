library(tidyverse)
library(mgcv)
library(vroom)
library(sf)
library(brms)
library(parallel)
library(lubridate)
library(pbapply)
library(lutz)
library(units)
library(Matrix)

theme_set(theme_minimal() + 
            theme(panel.grid.minor = element_blank()))

system("aws s3 sync s3://earthlab-amahood/night_fires/gamready data/out/gamready")

### MJK notes: 
# 1. we don't need to add in the sampling effort here; it is already in the gamready data
# 2. We also don't need to trim off the rows of the data that represented leading or 
# trailing 0's for GOES fire detections. That had to be done in the previous script 
# prior to aggregating across unique combinations of nid/VPD_hPa (because we lose the 
# rounded datetime info when we do that aggregation, so lose the ability to see which 0's
# come before or after the sequence of GOES fire detections)
# 3. We don't need to add in the total burned area for each fire event anymore, as that
# no longer needs to be an offset in the model given the response variable that we've 
# chosen to proceed with.
# 4. The distribution of the response has changed from nb() to binomial(), so the syntax
# of the model is a little different
# 5. If the GOES-16 satellite was offline, there might have been times where there were
# no GOES images taken at all. We filter these out in the dataframe prior to modeling using
# filter(n_scenes != 0)
# 6. The analysis ready data have 4 columns: nid (the unique identifier per FIRED event), 
# VPD_hPa (the VPD), fire_scenes (the count of all GOES scenes that had at least one fire
# detection within the nid event perimeter when the VPD was observed to be the value in the
# VPD_hPa column), and n_scenes (the count of *all* GOES scenes for a given nid/VPD_hPa
# combination)


# system("aws s3 sync s3://earthlab-amahood/night_fires/gam_progress data/gam_progress")
# system("aws s3 cp s3://earthlab-amahood/night_fires/lut_ba.Rda data/lut_ba.Rda")
# system("aws s3 cp s3://earthlab-mkoontz/goes16meta/sampling-effort-goes16.csv data/s_effort.csv")

gamready_files <- list.files("data/out/gamready", full.names = TRUE, pattern = ".csv$") %>%
  file.info() %>%
  as_tibble(rownames = "file") %>%
  arrange(size)

# load("data/lut_ba.Rda")
# lut_ba <- drop_units(lut_ba)

# hourdf<- vroom("data/s_effort.csv") %>%
#   mutate(rounded_datetime = ymd_hm(rounded_hour)) %>%
#   dplyr::rename(n_scenes=n)
# dir.create("data/gam_progress")

dir.create("data/mods")

# build models

for(i in 1:nrow(gamready_files)){

t0<- Sys.time()  

f<- gamready_files[i, "file"] %>% pull

out_fn_base <- 
  str_replace(f, "data/out/gamready/", "") %>%
  str_replace("_gamready.csv", "")

outname <- paste0("data/mods/", out_fn_base, 
                  "-gam.rds")

if (file.exists(outname)) next

print(out_fn_base)

events <- 
  vroom(f) %>%
  mutate(nid = factor(as.character(nid))) %>%
  # group_by(nid) %>%
  # mutate(first_detection = min(rounded_datetime[n > 0]), 
  #        hour_after_last_detection = max(rounded_datetime[n > 0]) + 60*60) %>%
  # filter(rounded_datetime >= first_detection, 
  #        rounded_datetime <= hour_after_last_detection) %>%
  # ungroup %>%
  # dplyr::select(-hour_after_last_detection, first_detection) %>%
  # mutate(ba = lut_ba[nid]) %>%
  # left_join(hourdf, by="rounded_datetime") %>%
  # mutate(effort = log(ba) * n_scenes) %>%
  filter(n_scenes != 0) %>% # drop any VPD_hPa/nid combinations that had 0 GOES images associated with it (due to orbital maneuvering and GOES imager being offline)
  filter(VPD_hPa < 30) %>% # there were some vpds of 70000 or so that were messing up the model... most dont go above 20 (mjk note: this isn't quite true for all landcovers; the landcover with the most events/rows is equatorial savannas and only 65% of data are below 30; we should consider a higher cutoff; 97.6% of data are below 100 and 97.3% of data are below 60. The same amount of data are below 100 as are below 1000 for this landcover type)
  droplevels()

  if (nrow(events) == 0) next
  if (length(unique(events$nid)) < 15) next

  levels(events$nid) <- c(levels(events$nid), "NewFire")
  m <- bam(cbind(fire_scenes, n_scenes - fire_scenes) ~ s(VPD_hPa) + s(VPD_hPa, nid, bs = "fs"), 
           data = events, 
           nthreads = parallel::detectCores()-1, 
           family = binomial(), 
           discrete = TRUE,
           drop.unused.levels = FALSE)
  
  # m <- brm(fire_scenes | trials(n_scenes) ~ s(VPD_hPa) + s(VPD_hPa, nid, bs = "fs"), 
  #          data = events, 
  #          chains = 4,
  #          cores = 4,
  #          family = binomial())
  
  write_rds(m, outname)
  
  t1 <- Sys.time() - t0
  
  blank_fn <- paste(out_fn_base,
                    round(t1), units(t1), 
        length(unique(events$nid)), "events", ".csv", sep="_")
  
  data.frame(x=1) %>%
    write_csv(file.path("data","gam_progress",blank_fn))
  system(paste("aws s3 cp", 
                file.path("data", "gam_progress", blank_fn),
                file.path("s3://earthlab-amahood", "night_fires", "gam_progress", blank_fn)))
  system(paste("aws s3 cp", 
               file.path("data", "mods", paste0(out_fn_base, "-gam.rds")),
               file.path("s3://earthlab-amahood", "night_fires", "gam_mods", paste0(out_fn_base, "-gam.rds"))))
}

# gotta fix everything below here:
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