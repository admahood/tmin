# Plot the VPD/AF time series for the case study fire
library(tidyverse)
library(vroom)
library(patchwork)
library(here)

long_pts <- read_csv(here("data", "case_study_data.csv"))

last_date_df <- long_pts %>%
  group_by(id) %>%
  filter(n > 0) %>%
  summarize(last_date = max(datetime_utc) + 60 * 60 * 24 * 1, 
            first_date = min(datetime_utc) - 60 * 60 * 24 * 1)

split_events <- long_pts %>%
  select(datetime_utc, n, id, vpd_kPa, Event,
         solar_radiation_w_per_m2) %>%
  left_join(last_date_df) %>%
  filter(datetime_utc <= last_date, 
         datetime_utc >= first_date) %>%
  pivot_longer(cols = c("n", "vpd_kPa")) %>%
  mutate(value = ifelse(!(name == "n" & value == 0), value, NA)) %>% 
  mutate(name = ifelse(name == "n", "AF counts", "VPD (kPa)")) %>%
  split(.$id)

# TODO: extract new land cover classes for focal events (using ignition lat/lon)
# and add a dashed line for the estimated VPD threshold
plot_list <- split_events %>%
  lapply(function(x) {
    event_name <- unique(x$Event)
    
    pct_night <- x %>%
      filter(name == "AF counts", !is.na(value)) %>%
      mutate(day = solar_radiation_w_per_m2 > .001) %>%
      group_by(day) %>%
      summarize(n_af = sum(value)) %>%
      ungroup()
    
    pct_night <- pct_night$n_af[!pct_night$day] / sum(pct_night$n_af)
    title <- paste0(event_name, ": ", 100 * round(pct_night, 2), 
                    "% nighttime fire detections")
    
    ggplot(x, aes(datetime_utc, value, color = solar_radiation_w_per_m2 > .001)) + 
      geom_point(size = .6, alpha= .8) + 
      geom_line(alpha = 0.3, color = "black") + 
      theme_minimal() + 
      facet_wrap(~name, scales = "free_y", ncol = 1, strip.position = "left") + 
      xlab("") + 
      ylab("") + 
      theme(panel.grid.minor = element_blank(), 
            legend.position = "none") +
      scale_color_manual(values = c("#2166ac", "red")) + 
      ggtitle(title) + 
      theme(plot.title = element_text(size=12))
  })

p <- wrap_plots(plot_list, ncol = 1)
p
ggsave("fig/case-study-ts.pdf", plot = p, width = 8, height = 6.5)
ggsave("fig/case-study-ts.png", plot = p, width = 8, height = 6.5)

