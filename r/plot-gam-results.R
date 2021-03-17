
# Visualize GAM predictions -----------------------------------------------

library(tidyverse)
library(vroom)
library(ggrepel)
library(ggthemes)
library(patchwork)
library(sf)

predictions <- list.files(pattern = "-predictions.csv$", 
                          path = "data/mods", 
                          full.names = TRUE) %>%
  lapply(vroom) %>%
  bind_rows

predictions

partial_df <- predictions %>%
  group_by(lc_name, VPD_hPa) %>%
  summarize(lo = quantile(logit_p, .025), 
            hi = quantile(logit_p, .975),
            mu = mean(logit_p), 
            n = n()) %>%
  ungroup


# Visualize partial effects -----------------------------------------------

plot_df <- partial_df %>%
  left_join(thresholds) %>%
  mutate(lc_name = paste0(lc_name, " (", round(vpd_thresh_hpa, 1), ')'))

partial_plot <- plot_df %>%
  ggplot(aes(VPD_hPa, plogis(mu), color = lc_name)) +
  geom_path() +
  geom_text_repel(aes(label = lc_name, x = VPD_hPa), 
                  data = partial_df %>%
                    group_by(lc_name) %>%
                    filter(VPD_hPa == max(VPD_hPa)) %>%
                    ungroup, 
                  nudge_x      = 0.4,
                  direction    = "y",
                  hjust        = 0,
                  segment.size = 0.1, size = 3, 
                  segment.alpha = .5,
                  segment.color = "black") + 
  scale_x_continuous(breaks = c(0:4) * 10, limits = c(0, 7) * 10) +
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(hjust = .24)) + 
  xlab("Vapor pressure deficit (hPa)") + 
  ylab("Active fire detection probability")
partial_plot

dir.create("fig", showWarnings = FALSE)
ggsave(plot = partial_plot, 
       filename = "fig/vpd-partial-effects.pdf", width = 7.5, height = 3.5)
ggsave(plot = partial_plot, 
       filename = "fig/vpd-partial-effects.png", width = 7.5, height = 3.5)





# Partial effects plots, with uncertainty
bayes_df <- predictions %>%
  mutate(new_mu = exp(Xb + log(mean(fire_areas$total_area_km2))))

bayes_plot <- bayes_df %>%
  droplevels() %>%
  left_join(thresholds) %>%
  ggplot(aes(VPD_hPa, new_mu, color = lc_name)) +
  geom_path(aes(group = j), alpha = .01) +
  scale_y_log10(breaks = c(0.01, 1, 100), 
                labels = c(0.01, 1, 100)) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank()) + 
  xlab("Vapor pressure deficit (kPa)") + 
  ylab("Active fire detections per hour") + 
  scale_alpha_manual(values = c(.1, 1)) + 
  scale_color_ptol() + 
  facet_wrap(~fct_reorder(lc_name, vpd_thresh_kpa), 
             nrow = 4, , labeller = label_wrap_gen()) + 
  geom_path(data = partial_df %>%
              mutate(lc_name = trimws(gsub("\\(.*", "", lc_name))) %>%
              left_join(thresholds), 
            aes(y = exp(mu)), size = 1) + 
  coord_cartesian(ylim = c(1e-4, 1e4)) + 
  geom_rug(data = filter(events, lc_name %in% bayes_df$lc_name) %>%
             filter(VPD_hPa <= quantile(events$VPD_hPa, 0.99)) %>%
             left_join(thresholds), 
           aes(x = VPD_hPa, color = lc_name), 
           inherit.aes = FALSE, alpha = .008)
bayes_plot
ggsave("fig/bayes_plot.png", bayes_plot, width = 5, height = 6)
ggsave("fig/bayes_plot.pdf", bayes_plot, width = 5, height = 6)
