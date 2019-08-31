library(tidyverse)
library(stringr)
results <- read_csv("results_for_figure_0224.csv") %>%
  mutate(`Method (line)` = str_sub(`Method (line)`, 7, 7)) %>% 
  # Number strength of counfounding 1 thru 3, not 0 thru 2:
  mutate(`Strength of confounding (x axis)` = `Strength of confounding (x axis)` + 1) %>% 
  # Explicitly label heterogeneous treatment effects:
  mutate(HTE = ifelse(HTE == 0, "HTE Scenario 1", "HTE Scenario 2")) %>% 
  # Explicitly label propensity score method:
  mutate(
    PropModel = case_when(
      PropModel == "prop1" ~ "Prop. Score 1",
      PropModel == "prop2" ~ "Prop. Score 2",
      PropModel == "prop3" ~ "Prop. Score 3",
      PropModel == "prop4" ~ "Prop. Score 4",
      PropModel == "prop5" ~ "Prop. Score 5"
    )
  ) %>% 
  # Explictly label method
  mutate(
    `Method (line)` = case_when(
      `Method (line)` == 1 ~ "Regression only",
      `Method (line)` == 2 ~ "Regression w/ \nbalancing",
      `Method (line)` == 3 ~ "Regression w/ \nprop. scores"
    ),
    `Method (line)` = 
      factor(`Method (line)`, 
             levels = c("Regression only", "Regression w/ \nbalancing", "Regression w/ \nprop. scores")
      )
  ) %>% 
  group_by(`PropHTE (facet)`, `Strength of confounding (x axis)`, `Method (line)`, PropModel, HTE) %>% 
  summarize(
    `Bias (y axis)` = mean(`Bias (y axis)`), 
    `Variance (y axis)` = mean(`Variance (y axis)`)
  )

# Two options for faceting
# baseplot_1 <- 
#   ggplot(results, aes(x = `Strength of confounding (x axis)`, col = `Method (line)`, shape = `Method (line)`)) +
#   facet_wrap(~`PropHTE (facet)`, nrow = 2) +
#   labs(x = "Strength of Confounding Case", col = "Method") +
#   # Changed to more colorblind friendly color palette:
#   scale_color_brewer(palette = "Set2") +
#   # Change x-axis to only have 1 thru 3 as tick marks:
#   scale_x_continuous(breaks=c(1, 2, 3))

baseplot_2 <- 
  ggplot(results, aes(x = `Strength of confounding (x axis)`, col = `Method (line)`, shape = `Method (line)`)) +
  facet_grid(HTE ~ PropModel) +
  labs(x = "Strength of Confounding Case", col = "Method") +
  # Changed to more colorblind friendly color palette:
  scale_color_brewer(palette = "Set2") +
  # Change x-axis to only have 1 thru 3 as tick marks:
  scale_x_continuous(breaks=c(1, 2, 3))


# Bias
# baseplot_1 +
#   geom_line(aes(y = `Bias (y axis)`, linetype = `Method (line)`)) +
#   geom_point(aes(y = `Bias (y axis)`), show.legend = F) +
#   labs(y = "Bias", title = "Bias", linetype = "Method")
# ggsave("figures/bias_1_new.png", width=8, height=4.5)

baseplot_2 +
  geom_line(aes(y = `Bias (y axis)`, linetype = `Method (line)`)) +
  geom_point(aes(y = `Bias (y axis)`), show.legend = F) +
  labs(y = "Bias", title = "Bias", linetype = "Method")
ggsave("figures/bias_2_new.png", width=8, height=4.5)


# Variance
# baseplot_1 +
#   geom_line(aes(y = `Variance (y axis)`)) +
#   geom_point(aes(y = `Variance (y axis)`), show.legend = F) +
#   labs(y = "Variance", title = "Variance")
# ggsave("figures/variance_1_new.png", width=8, height=4.5)

baseplot_2 +
  geom_line(aes(y = `Variance (y axis)`)) +
  geom_point(aes(y = `Variance (y axis)`), show.legend = F) +
  labs(y = "Variance", title = "Variance")
ggsave("figures/variance_2_new.png", width=8, height=4.5)


