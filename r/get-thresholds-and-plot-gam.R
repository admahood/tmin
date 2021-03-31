library(tidyverse)
library(vroom)
library(ggrepel)
library(ggthemes)
library(patchwork)
library(sf)
library(stringr)

predictions <- list.files(pattern = "-predictions.csv$", 
                          path = "data/mods", 
                          full.names = TRUE) %>%
  lapply(vroom) %>%
  bind_rows %>%
  mutate(lc_name = gsub("_", " ", lc_name))

# this data frame contains posterior predictive draws for every land cover type
predictions



# Compute vpd thresholds --------------------------------------------------
# These are posterior means of VPD values for which the probability of a zero
# (not detecting a fire) is 0.95

thresholds <- predictions %>%
  mutate(abs_diff = abs(pr_zero - 0.95)) %>%
  group_by(j, lc_name) %>%
  filter(abs_diff == min(abs_diff)) %>%
  ungroup %>%
  group_by(lc_name) %>%
  summarize(vpd_thresh_hpa = mean(VPD_hPa), 
            sd = sd(VPD_hPa)) %>%
  mutate(probability_of_zero = 0.95) %>%
  arrange(vpd_thresh_hpa)

thresholds %>%
  write_csv("data/mods/zero-goes-af-vpd-thresholds.csv")



# Visualize partial effects -----------------------------------------------

partial_df <- predictions %>%
  mutate(p = plogis(logit_p)) %>%
  group_by(lc_name, VPD_hPa) %>%
  summarize(lo = quantile(p, .05), 
            hi = quantile(p, .95),
            mu = median(p), 
            n = n()) %>%
  ungroup %>%
  mutate(color_group = gsub( " .*$", "", lc_name), 
         facet_label = paste(color_group, "ecoregions"))

plot_df <- partial_df %>%
  left_join(thresholds) %>%
  mutate(lc_name = paste0(lc_name, " (", round(vpd_thresh_hpa, 1), ')'))

partial_plot <- plot_df %>%
  ggplot(aes(VPD_hPa, mu, color = facet_label)) +
  geom_path(alpha = .5, aes(group = lc_name)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = lc_name), 
              color = NA, alpha = .02) +
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(hjust = .24)) + 
  xlab("Vapor pressure deficit (hPa)") + 
  ylab("Active fire detection probability") + 
  facet_wrap(~reorder(facet_label, vpd_thresh_hpa))
partial_plot

ggsave(plot = partial_plot, 
       filename = "images/vpd-partial-effects.pdf", width = 7.5, height = 3.5)
ggsave(plot = partial_plot, 
       filename = "images/vpd-partial-effects.png", width = 7.5, height = 3.5)


# Supplementary partial effects plots, with uncertainty ----------------------

bayes_plot <- predictions %>%
  left_join(thresholds) %>%
  ggplot(aes(VPD_hPa, 1-pr_zero)) +
  geom_path(aes(group = j), alpha = .01) +
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank()) + 
  xlab("Vapor pressure deficit (hPa)") + 
  ylab("Active fire detection probability") + 
  facet_wrap(~fct_reorder(lc_name, vpd_thresh_hpa), 
             nrow = 5, , labeller = label_wrap_gen()) + 
  geom_path(data = partial_df %>%
              mutate(lc_name = trimws(gsub("\\(.*", "", lc_name))) %>%
              left_join(thresholds), 
            aes(y = mu), size = 1)
bayes_plot
ggsave("images/bayes_plot.png", bayes_plot, width = 9, height = 6)
ggsave("images/bayes_plot.pdf", bayes_plot, width = 5, height = 6)
