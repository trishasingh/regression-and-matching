# Final Simulations

# Setup -------------------------------------------------------------------
# setwd("~/Documents/Stats Indep Study/Simulations")

# tidyverse includes ggplot2, tibble, tidyr, readr, purrr, and dplyr:
library(tidyverse)
library(stringr)
library(broom)
library(Matching)
library(MatchIt)

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
  dplyr::select(ID, A, B)

# Number of simulations
n_sim <- 10

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

for(i in 1:2){
  # Set propensity scores, potential outcomes & delta based on simulation number
  if(i == 1){
    P <- P_1
  } else if (i == 2) {
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
        simulation = i,
        iteration = k
      ) %>% 
      # dplyr::select(simulation, iteration, ID, di) %>% 
      dplyr::select(di) %>% 
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
      match.data() %>% 
      lm(yi ~ di + B, data = .) %>% 
      tidy()
    
    #Estimation 2.2
    model_2_2 <- P %>% 
      # Run matching
      matchit(di ~ A + B, method = "nearest", replace = TRUE, data = .) %>% 
      match.data() %>% 
      lm(yi ~ di + A + B, data = .) %>% 
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
      mutate(wi = di/pi + (1-di)/(1-pi)) %>% 
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
      mutate(wi = di/pi + (1-di)/(1-pi)) %>% 
      lm(yi ~ di + A + B, weights = wi, data = .) %>% 
      tidy()
    
    # Save simulation results for this k
    ID_variables <- data_frame(
      simulation = i,
      iteration = k,
      model = model_names
    )
    beta_di <- 
      bind_rows(
        model_1_0, model_1_1, model_1_2, 
        model_2_1, model_2_2, model_3_1, model_3_2
      ) %>% 
      filter(term == "di") %>% 
      dplyr::select(estimate, std.error)
    
    # Append to overall results
    simulation_results <- simulation_results %>% 
      bind_rows(bind_cols(ID_variables, beta_di))
    
    # Show simulation counter and save results
    if(k %% 10 == 1){
      str_c("Sim", i, "Iter", k, "of", n_sim, "done at", Sys.time(), sep=" ") %>%
        print()
      save(simulation_results, file="simulation_results.RData")
      save(saved_P, file="saved_P.RData")
    }
  }
}
save(simulation_results, file="simulation_results.RData")
save(saved_P, file="saved_P.RData")




load("simulation_results.RData")

ATE <- bind_rows(
  P_1 %>%
    summarise(ATE = mean(delta)) %>% 
    mutate(simulation = 1),
  P_2 %>%
    summarise(ATE = mean(delta)) %>% 
    mutate(simulation = 2)
)

simulation_results_1 <- simulation_results_1 %>% 
  mutate(simulation = 1)
simulation_results_2 <- simulation_results_2 %>% 
  mutate(simulation = 2)

simulation_results_2 %>% 
  #bind_rows(simulation_results_2) %>% 
  mutate(
    model = str_sub(model, start=7),
    method = str_sub(model, start=1, end=1),
    method = as.factor(method)
    ) %>% 
  left_join(ATE, by="simulation") %>% 
  mutate(bias = estimate - ATE) %>% 
  dplyr::select(model, method, simulation, bias, std.error) %>% 
  gather(type, value, -c(model, method, simulation)) %>% 
  ggplot() +
  geom_boxplot(aes(x=model, y=value, fill=method)) +
  facet_wrap(~type, scales = "free") +
  labs(title = "Simulation 2", y = "") +
  scale_fill_brewer(palette = "Set2")
  # facet_wrap(~simulation) + 
  # geom_hline(yintercept = 0, linetype = "dashed") +


simulation_results %>% 
  ggplot() +
  geom_boxplot(aes(x=model, y=bias)) +
  # facet_wrap(~simulation) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title="Bias of Estimator")

simulation_results %>% 
  ggplot() +
  geom_boxplot(aes(x=model, y=std.error)) +
  # facet_wrap(~simulation) +
  labs(title="Standard Error of Estimator")




