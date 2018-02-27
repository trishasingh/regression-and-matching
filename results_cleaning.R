# Setup
library(tidyverse)
library(stringr)
setwd("~/Dropbox/Matching Paper/Matching Git")

results <- read_csv("results_for_figure.csv") %>%
  mutate(`Method (line)` = str_sub(`Method (line)`, 7, 7)) %>% 
  group_by(`PropHTE (facet)`, `Strength of confounding (x axis)`, `Method (line)`, PropModel, HTE) %>% 
  summarize(`Bias (y axis)` = mean(`Bias (y axis)`), `SE (y axis)` = mean(`SE (y axis)`),
            `RMSE (y axis)` = sqrt(mean(`RMSE (y axis)`)))

# Relative Bias Plot

baseplot <- ggplot(results, aes(x = `Strength of confounding (x axis)`, col = `Method (line)`, shape = `Method (line)`)) +
  facet_grid(HTE ~ PropModel) +
  labs(x = "Strength of Confounding", col = "Method")

baseplot +
  geom_line(aes(y = `Bias (y axis)`)) +
  geom_point(aes(y = `Bias (y axis)`), show.legend = F) +
  labs(y = "Bias", title = "Bias")
ggsave("figures/bias.png", width=8, height=4.5)

# SE Plot

baseplot +
  geom_line(aes(y = `SE (y axis)`)) +
  geom_point(aes(y = `SE (y axis)`), show.legend = F) +
  labs(y = "Standard Error", title = "Standard Error")
ggsave("figures/se.png", width=8, height=4.5)

# RMSE Plot

baseplot +
  geom_line(aes(y = `RMSE (y axis)`)) +
  geom_point(aes(y = `RMSE (y axis)`), show.legend = F) +
  labs(y = "RMSE", title = "RMSE")
ggsave("figures/rmse.png", width=8, height=4.5)