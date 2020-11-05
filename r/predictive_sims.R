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
system("aws s3 sync s3://earthlab-amahood/night_fires/gam_progress data/gam_progress")
system("aws s3 sync s3://earthlab-amahood/night_fires/gam_mods data/gam_mods")

system("aws s3 cp s3://earthlab-amahood/night_fires/lut_ba.Rda data/lut_ba.Rda")
system("aws s3 cp s3://earthlab-mkoontz/goes16meta/sampling-effort-goes16.csv data/s_effort.csv")


gamready_files <- list.files("data/gamready", full.names = TRUE, pattern = ".csv$") %>%
  file.info() %>%
  as_tibble(rownames = "file") %>%
  arrange(size)
load("data/lut_ba.Rda")
lut_ba <- drop_units(lut_ba)

hourdf<- vroom("data/s_effort.csv") %>%
  mutate(rounded_datetime = ymd_hm(rounded_hour)) %>%
  dplyr::rename(n_scenes=n)
dir.create("data/gam_progress")
dir.create("data/mods")
# build models


# gotta fix everything below here:
# For each land cover type, simulate from the predictive distribution
cover_types <- list.files("data/gam_mods", pattern="rds$")%>%
  str_replace_all("-gam.rds","")

predictive_sim <- function(cover_type) {
  
  subd <- vroom(str_c("data/gamready/",cover_type, "_gamready.csv")) %>%
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
    mutate(effort = log(ba) * n_scenes) %>%
    na.omit() %>%
    filter(VPD_hPa < 30) %>% # there were some vpds of 70000 or so that were messing up the model... most dont go above 20
    droplevels()
  levels(subd$nid) <- c(levels(subd$nid), "NewFire")
  
  model_path <- paste0("data/gam_mods/",cover_type, 
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
    distinct(VPD_hPa, n_pred, ba, nid) %>%
    ggplot(aes(VPD_hPa, n_pred / ba, group = nid)) + 
    geom_path(alpha = .05) +
    scale_y_log10() 
  
  # Simulate from predictive distribution
  pred_df <- subd %>%
    # distinct(lc_name) %>%
    mutate(ba = subd %>% 
             distinct(nid, effort) %>% 
             summarize(median(effort)) %>% 
             unlist %>% 
             c,
           VPD_hPa = list(seq(0, max(subd$VPD_hPa), by = .05)), 
           nid = "NewFire") %>%
    unnest(cols = c(VPD_hPa)) %>%
    mutate(idx = 1:n())
  Xp <- predict(m, pred_df, type = "lpmatrix")
  nonzero_cols <- colSums(Xp) != 0
  beta <- coef(m)[nonzero_cols]
  Vb <- vcov(m)[nonzero_cols, nonzero_cols] ## posterior mean and cov of coefs
  n <- 1000
  br <- MASS::mvrnorm(n, beta, Vb) 
  dfs <- list()
  pb <- txtProgressBar(max = n, style = 3)
  stopifnot(length(unique(pred_df$effort)) == 1)
  for (j in 1:n) {
    dfs[[j]] <- tibble(Xb = Xp[, nonzero_cols] %*% br[j, ] %>% 
                         as.vector, 
                       j = j) %>%
      mutate(idx = 1:n(), 
             log_mu = Xb + log(pred_df$effort[1]))
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  posterior_predictions <- dfs %>%
    bind_rows %>%
    left_join(pred_df) %>%
    arrange(j, VPD_hPa)  %>%
    mutate(y_pred = MASS::rnegbin(n(), 
                                  exp(log_mu), 
                                  theta = exp(m$family$getTheta())), 
           pr_zero = dnbinom(0, size = exp(m$family$getTheta()), mu = exp(log_mu)))
  posterior_predictions %>%
    write_csv(prediction_csv_path)
}

sims <- pblapply(cover_types, predictive_sim)