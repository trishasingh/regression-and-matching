library(tidyverse)
library(stringr)
results <- read_csv("results_for_figure_0224.csv") %>%
  mutate(`Method (line)` = str_sub(`Method (line)`, 7, 7)) %>% 
  group_by(`PropHTE (facet)`, `Strength of confounding (x axis)`, `Method (line)`, PropModel) %>% 
  summarize(`Bias (y axis)` = mean(`Bias (y axis)`), `Variance (y axis)` = mean(`Variance (y axis)`))

baseplot_1 <- ggplot(results, aes(x = `Strength of confounding (x axis)`, col = `Method (line)`, shape = `Method (line)`)) +
  facet_wrap(~`PropHTE (facet)`, nrow = 2) +
  labs(x = "Strength of Confounding", col = "Method")

baseplot_1 +
  geom_line(aes(y = `Bias (y axis)`)) +
  geom_point(aes(y = `Bias (y axis)`), show.legend = F) +
  labs(y = "Bias", title = "Bias")
ggsave("figures/bias_1.png", width=8, height=4.5)

baseplot_2 <- ggplot(results, aes(x = `Strength of confounding (x axis)`, col = `Method (line)`, shape = `Method (line)`)) +
  facet_grid(HTE ~ PropModel) +
  labs(x = "Strength of Confounding", col = "Method")

baseplot_2 +
  geom_line(aes(y = `Bias (y axis)`)) +
  geom_point(aes(y = `Bias (y axis)`), show.legend = F) +
  labs(y = "Bias", title = "Bias")
ggsave("figures/bias_2.png", width=8, height=4.5)