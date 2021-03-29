library(tidyverse)
library(mgcv)
library(vroom)
library(sf)
library(parallel)
library(lubridate)
library(pbapply)
library(units)
library(Matrix)

set.seed(1234)

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


gamready_files <- list.files("data/out/gamready", full.names = TRUE, pattern = ".csv$") %>%
  file.info() %>%
  as_tibble(rownames = "file") %>%
  arrange(-size)

dir.create("data/mods", showWarnings = FALSE, recursive = TRUE)

# Fit models for each ecoregion
for(i in 1:nrow(gamready_files)){
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
    filter(n_scenes != 0) %>% # drop any VPD_hPa/nid combinations that had 0 GOES images associated with it (due to orbital maneuvering and GOES imager being offline)
    filter(VPD_hPa < 100) %>% 
    droplevels()

  if (nrow(events) == 0) next
  if (length(unique(events$nid)) < 200) next
  
  # subsample for ecoregions with lots of rows
  max_nrow <- 100000
  if (nrow(events) > max_nrow) {
    print(paste("subsampling for", out_fn_base))
    print(paste("initial nrow =", nrow(events)))
    while (nrow(events) > max_nrow) {
      current_nrow <- nrow(events)
      frac_to_keep <- max_nrow / current_nrow * .99
      events_to_keep <- sample(unique(events$nid), 
                               size = round(frac_to_keep * 
                                              length(unique(events$nid))), 
                               replace = FALSE)
      events <- filter(events, nid %in% events_to_keep) %>%
        droplevels()
    }
    print(paste("final nrow =", nrow(events)))
  }

  levels(events$nid) <- c(levels(events$nid), "NewFire")
  m <- bam(cbind(fire_scenes, n_scenes - fire_scenes) ~ VPD_hPa + s(nid, bs = "re"),
           data = events,
           nthreads = parallel::detectCores(),
           family = binomial(),
           discrete = TRUE,
           drop.unused.levels = FALSE)
  m2 <- bam(cbind(fire_scenes, n_scenes - fire_scenes) ~ s(VPD_hPa, k=3) + s(nid, bs = "re"),
            data = events,
            nthreads = parallel::detectCores(),
            family = binomial(),
            discrete = TRUE,
            drop.unused.levels = FALSE)
  write_rds(m, outname)
  
  # simulate from the posterior predictive distribution
  pred_df <- tibble(lc_name = out_fn_base,
                    VPD_hPa = 0:92, 
                    nid = "NewFire", 
                    n_scenes = 1) %>%
    mutate(idx = 1:n())
  
  # generate approximate draws from the posterior 
  # (assuming that we can approximate the posterior with a multivariate normal)
  # note: treats the random effect standard deviation as "known"
  # in other words, if our "NewFire" factor level adjustment is z, this is:
  # z ~ Normal(0, sd_z), where sd_z is a point estimate from the model
  Xp <- predict(m, pred_df, type = "lpmatrix")
  nonzero_cols <- colSums(Xp) != 0
  beta <- coef(m)[nonzero_cols]
  Vb <- vcov(m)[nonzero_cols, nonzero_cols] ## posterior mean and cov of coefs
  n <- 1000
  br <- MASS::mvrnorm(n, beta, Vb) 
  
  # for each posterior draw, compute logit pr of detection
  dfs <- list()
  pb <- txtProgressBar(max = n, style = 3)
  for (j in 1:n) {
    dfs[[j]] <- tibble(logit_p = Xp[, nonzero_cols] %*% br[j, ] %>% 
                         as.vector, 
                       j = j) %>%
      mutate(idx = 1:n())
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  # join all draws into big df and simulate values/compute probabilities
  # this data frame should have (n_draws x n_vpd_values) rows
  posterior_predictions <- dfs %>%
    bind_rows %>%
    left_join(pred_df) %>%
    arrange(lc_name, j, VPD_hPa)  %>%
    # simulating values from the approximate posterior
    mutate(y_pred = rbinom(n(), size = 1, prob = plogis(logit_p)), 
           pr_zero = dbinom(0, size = 1, prob = plogis(logit_p)), 
           n_events = length(unique(events$nid)))
  posterior_predictions %>%
    write_csv(gsub("-gam.rds", "-gam-predictions.csv", outname))
}
