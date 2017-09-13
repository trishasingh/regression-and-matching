# Bert's Simulations

# Final Simulations

# Setup -------------------------------------------------------------------
# setwd("~/Documents/Stats Indep Study/Simulations")

library(Matching)
library(MatchIt)
library(stringr)
library(broom)
library(forcats)
library(knitr)
# tidyverse includes ggplot2, tibble, tidyr, readr, purrr, and dplyr.
# load this last as there are command conflicts with filter() and select() from
# MASS package, which gets loaded by Matching.
library(tidyverse)


# Set random number generator seed value
set.seed(76)

# Create data frame of all possible covariate combinations
A_values <- seq(from = .01, to = 1, by = .01)
B_values <- seq(from = .01, to = 1, by = .01)
covariates <- data_frame(
  A = rep(A_values, times=length(B_values)),
  B = rep(B_values, each=length(A_values))
) %>%
  mutate(ID=1:n()) %>%
  select(ID, A, B)

# Number of simulations
n_sim <- 1000

# Names of estimation models
model_names <- c("1_0", "1_1", "1_2", "2_1", "2_2", "3_1", "3_2")





# Simulation parameters -------------------------------------------------

# Functions to compute propensity scores based on covariate values
compute_prop_score_1 <- function(a, b) {
  x <- a + a^2 + b
  prop_score <- 1/(1+exp(-x))
  return(prop_score)
}
compute_prop_score_2 <- function(a, b) {
  x <- a + b
  prop_score <- 1/(1+exp(-x))
  return(prop_score)
}

# Compute propensity score for each combination of independent variables and
# then generate potential outcomes
P_1 <- covariates %>%
  mutate(
    prop_score = compute_prop_score_1(A, B),
    y0 = 100 + 3*A + 2*B + rnorm(n = n(), mean = 0, sd = 5),
    y1 = 102 + 6*A + 4*B + rnorm(n = n(), mean = 0, sd = 5),
    delta = y1 - y0
  )
P_2 <- covariates %>%
  mutate(
    prop_score = compute_prop_score_2(A, B),
    y0 = 100 + 3*A + 2*B + rnorm(n = n(), mean = 0, sd = 5),
    y1 = 102 + 6*A + 4*B + 5*prop_score + rnorm(n = n(), mean = 0, sd = 5),
    delta = y1 - y0
  )

# Compare plots of propensity scores
bind_rows(
  mutate(P_1, Simulation = "Simulation 1"),
  mutate(P_2, Simulation = "Simulation 2")
) %>%
  ggplot(aes(x=A, y=B)) +
  geom_tile(aes(fill=prop_score)) +
  geom_contour(aes(z=prop_score)) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(
    x = "Covariate A",
    y = "Covariate B",
    fill = "Propensity Score",
    title = "Figure 1: Propensity Score Function for Each Simulation"
  ) +
  facet_wrap(~Simulation) +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# Only one of two propensity scores
ggplot(P_1, aes(x=A, y=B)) +
  geom_tile(aes(fill=prop_score)) +
  geom_contour(aes(z=prop_score)) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(
    x = "Covariate A",
    y = "Covariate B",
    fill = "Propensity Score",
    title="Figure 1: Propensity Score Function for Simulation 1"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))





# Main Simulation Loop ----------------------------------------------------
# Save results here
simulation_results <- NULL
saved_P <- NULL

for(sim in 1:2){
  # Set propensity scores, potential outcomes & delta based on simulation number
  if(sim == 1){
    P <- P_1
  } else if (sim == 2) {
    P <- P_2
  }
  
  # Main for loop
  for (k in 1:n_sim) {
    # Randomly assign treatment with probability = propensity score, and then
    # accordingly set which of potential outcomes to set as observed outcome
    P <- P %>%
      mutate(
        di = rbinom(n = n(), size = 1, prob = prop_score),
        yi = ifelse(di == 1, y1, y0)
      )
    # Save treated/untreated variable and resulting observation
    saved_P <- P %>%
      mutate(
        simulation = sim,
        iteration = k
      ) %>%
      # select(simulation, iteration, ID, di) %>%
      select(di) %>%
      bind_rows(saved_P, .)
    
    # Estimation 1.0
    model_1_0 <- P %>%
      lm(yi ~ di, data = .) %>%
      tidy()
    
    # Estimation 1.1
    model_1_1 <- P %>%
      lm(yi ~ di + B, data = .) %>%
      tidy()
    
    # Estimation 1.2
    model_1_2 <- P %>%
      lm(yi ~ di + A + B, data = .) %>%
      tidy()
    
    # Estimation 2.1
    model_2_1 <- P %>%
      # Run matching
      matchit(di ~ B, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi ~ di + B, weights = weights, data = .) %>%
      tidy()
    
    #Estimation 2.2
    model_2_2 <- P %>%
      # Run matching
      matchit(di ~ A + B, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi ~ di + A + B, weights = weights, data = .) %>%
      tidy()
    
    # Estimation 3.1
    model_3_1 <- P %>%
      # Compute fitted probability of
      glm(di ~ B, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di + (pi*(1-di))/(1-pi)) %>%
      lm(yi ~ di + B, weights = wi, data = .) %>%
      tidy()
    
    # Estimation 3.2
    model_3_2 <- P %>%
      # Compute fitted probability of
      glm(di ~ A + B, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di + (pi*(1-di))/(1-pi)) %>%
      lm(yi ~ di + A + B, weights = wi, data = .) %>%
      tidy()
    
    # Save simulation results for this k
    ID_variables <- data_frame(
      simulation = sim,
      iteration = k,
      model = model_names
    )
    beta_di <-
      bind_rows(
        model_1_0, model_1_1, model_1_2,
        model_2_1, model_2_2, model_3_1, model_3_2
      ) %>%
      filter(term == "di") %>%
      select(estimate, std.error)
    
    # Append to overall results
    simulation_results <- simulation_results %>%
      bind_rows(bind_cols(ID_variables, beta_di))
    
    # Show simulation counter and save results
    if(k %% 10 == 1){
      str_c("Sim", sim, "Iter", k, "of", n_sim, "done at", Sys.time(), sep=" ") %>%
        print()
      save(simulation_results, file="simulation_results.RData")
      save(saved_P, file="saved_P.RData")
    }
  }
}
save(simulation_results, file="simulation_results.RData")
save(saved_P, file="saved_P.RData")

# Above (for n_sim=1000) took X time to run on Amherst College RStudio server






# Analyze Results ---------------------------------------------------------
# Compute average treatment effect for both simulations
ATE <- bind_rows(
  P_1 %>%
    summarise(ATE = mean(delta)) %>%
    mutate(simulation = 1),
  P_2 %>%
    summarise(ATE = mean(delta)) %>%
    mutate(simulation = 2)
)

# Load saved simulation results and compute bias
load("simulation_results.RData")
# Delete later:
# n_sim <- simulation_results$iteration %>% max

simulation_results <- simulation_results %>%
  mutate(
    method = str_sub(model, start=1, end=1),
    method = factor(method),
    method = fct_recode(method,
                        "Just Regression" = "1",
                        "Matching" = "2",
                        "Inverse Prop Score Weights" = "3"
    )
  ) %>%
  left_join(ATE, by="simulation") %>%
  mutate(bias = estimate - ATE) %>%
  select(model, method, simulation, bias, std.error)


# Plot
simulation_results %>%
  gather(type, value, -c(model, method, simulation)) %>%
  ggplot() +
  geom_boxplot(aes(x=model, y=value, fill=method)) +
  facet_grid(type~simulation, scales = "free") +
  labs(title = str_c(n_sim, " iterations"), y = "") +
  scale_fill_brewer(palette = "Set2")

# Table
simulation_results %>%
  group_by(simulation, model) %>%
  summarise(bias=mean(bias), std.error = mean(std.error)) %>%
  kable(digits=4)
