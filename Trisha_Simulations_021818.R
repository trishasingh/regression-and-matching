# Trisha's Simulations

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
library(MASS)


# Set random number generator seed value
set.seed(76)

# Create data frame of all possible covariate combinations
n_rows <- 1000
# Generate covariance matrix
M = matrix(c(1, 0, 0, 0, 0.2, 0, 0, 0, 0, 0,
             0, 1, 0, 0, 0, 0.9, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0, 0, 0.2, 0, 0,
             0, 0, 0, 1, 0, 0, 0, 0, 0.9, 0,
             0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow = 10, ncol = 10)
# Generate values from normal random variables using the covariance matrix
L = chol(M)
nvars = dim(L)[1]
r = t(L) %*% matrix(rnorm(nvars*n_rows, mean=0, sd=1), nrow=nvars, ncol=n_rows)
r = t(r)
covariates <- data.frame(r) %>% 
  mutate(ID=1:n()) %>% 
  select(ID, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)

# Number of simulations
n_sim <- 1

# Names of estimation models
model_names <- c("1_prop1_conf1_hte0", "1_prop2_conf1_hte0", "1_prop3_conf1_hte0", "1_prop4_conf1_hte0", "1_prop5_conf1_hte0",
                 "1_prop1_conf2_hte0", "1_prop2_conf2_hte0", "1_prop3_conf2_hte0", "1_prop4_conf2_hte0", "1_prop5_conf2_hte0",
                 "1_prop1_conf3_hte0", "1_prop2_conf3_hte0", "1_prop3_conf3_hte0", "1_prop4_conf3_hte0", "1_prop5_conf3_hte0",
                 "1_prop1_conf1_hte1", "1_prop2_conf1_hte1", "1_prop3_conf1_hte1", "1_prop4_conf1_hte1", "1_prop5_conf1_hte1",
                 "1_prop1_conf2_hte1", "1_prop2_conf2_hte1", "1_prop3_conf2_hte1", "1_prop4_conf2_hte1", "1_prop5_conf2_hte1",
                 "1_prop1_conf3_hte1", "1_prop2_conf3_hte1", "1_prop3_conf3_hte1", "1_prop4_conf3_hte1", "1_prop5_conf3_hte1",
                 
                 "2_prop1_conf1_hte0", "2_prop2_conf1_hte0", "2_prop3_conf1_hte0", "2_prop4_conf1_hte0", "2_prop5_conf1_hte0",
                 "2_prop1_conf2_hte0", "2_prop2_conf2_hte0", "2_prop3_conf2_hte0", "2_prop4_conf2_hte0", "2_prop5_conf2_hte0",
                 "2_prop1_conf3_hte0", "2_prop2_conf3_hte0", "2_prop3_conf3_hte0", "2_prop4_conf3_hte0", "2_prop5_conf3_hte0",
                 "2_prop1_conf1_hte1", "2_prop2_conf1_hte1", "2_prop3_conf1_hte1", "2_prop4_conf1_hte1", "2_prop5_conf1_hte1",
                 "2_prop1_conf2_hte1", "2_prop2_conf2_hte1", "2_prop3_conf2_hte1", "2_prop4_conf2_hte1", "2_prop5_conf2_hte1",
                 "2_prop1_conf3_hte1", "2_prop2_conf3_hte1", "2_prop3_conf3_hte1", "2_prop4_conf3_hte1", "2_prop5_conf3_hte1",
                 
                 "3_prop1_conf1_hte0", "3_prop2_conf1_hte0", "3_prop3_conf1_hte0", "3_prop4_conf1_hte0", "3_prop5_conf1_hte0",
                 "3_prop1_conf2_hte0", "3_prop2_conf2_hte0", "3_prop3_conf2_hte0", "3_prop4_conf2_hte0", "3_prop5_conf2_hte0",
                 "3_prop1_conf3_hte0", "3_prop2_conf3_hte0", "3_prop3_conf3_hte0", "3_prop4_conf3_hte0", "3_prop5_conf3_hte0",
                 "3_prop1_conf1_hte1", "3_prop2_conf1_hte1", "3_prop3_conf1_hte1", "3_prop4_conf1_hte1", "3_prop5_conf1_hte1",
                 "3_prop1_conf2_hte1", "3_prop2_conf2_hte1", "3_prop3_conf2_hte1", "3_prop4_conf2_hte1", "3_prop5_conf2_hte1",
                 "3_prop1_conf3_hte1", "3_prop2_conf3_hte1", "3_prop3_conf3_hte1", "3_prop4_conf3_hte1", "3_prop5_conf3_hte1")




# Simulation parameters -------------------------------------------------

# Functions to compute propensity scores based on covariate values (from Leacy and Stuart)
a0 <- -1.897
a1 <- 0.8
a2 <- -0.25
a3 <- 0.6
a4 <- -0.4
a5 <- -0.8
a6 <- -0.5
a7 <- 0.7

#Additivity and linearity
compute_prop_score_1 <- function(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) {
  x <- a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 + a5*X5 + a6*X6 + a7*X7
  prop_score <- 1/(1+exp(-x))
  return(prop_score)
}
#Moderate non-additivity
compute_prop_score_2 <- function(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) {
  x <- a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 + a5*X5 + a6*X6 + a7*X7 +
    0.5*a1*X1*X3 + 0.7*a2*X2*X4 + 0.5*a3*X3*X5 +
    0.7*a4*X4*X6 + 0.5*a5*X5*X7 + 0.5*a1*X1*X6 +
    0.7*a2*X2*X3 + 0.5*a3*X3*X4 + 0.5*a4*X4*X5 +
    0.5*a5*X5*X6
  prop_score <- 1/(1+exp(-x))
  return(prop_score)
}
#Mild non-additivity and non-linearlity
compute_prop_score_3 <- function(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) {
  x <- a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 + a5*X5 + a6*X6 + a7*X7 +
    a2*X2*X2 + 0.5*a1*X1*X3 + 0.7*X2*X4 + 0.5*a4*X4*X5 + 0.5*a5*X5*X6 
  prop_score <- 1/(1+exp(-x))
  return(prop_score)
}
#Moderate non-linearity
compute_prop_score_4 <- function(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) {
  x <- a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 + a5*X5 + a6*X6 + a7*X7 +
    a2*X2*X2 + a4*X4*X4 + a7*X7*X7
  prop_score <- 1/(1+exp(-x))
  return(prop_score)
}
#Moderate non-additivity and non-linearity
compute_prop_score_5 <- function(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10) {
  x <- a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 + a5*X5 + a6*X6 + a7*X7 +
    a2*X2*X2 + a4*X4*X4 + a7*X7*X7 + 0.5*a1*X1*X3 + 0.7*a2*X2*X4 + 
    0.5*a3*X3*X5 + 0.7*a4*X4*X6 + 0.5*a5*X5*X7 + 0.5*a1*X1*X6 +
    0.7*a2*X2*X3 + 0.5*a3*X3*X4 + 0.5*a4*X4*X5 + 0.5*a5*X5*X6
  prop_score <- 1/(1+exp(-x))
  return(prop_score)
}

# Compute propensity score for each combination of independent variables and
# then generate potential outcomes
P <- covariates %>% 
  mutate(prop1 = compute_prop_score_1(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10),
         prop2 = compute_prop_score_2(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10),
         prop3 = compute_prop_score_3(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10),
         prop4 = compute_prop_score_4(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10),
         prop5 = compute_prop_score_5(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10))

b0 = -1.386
b1 = 0.3
b2 = -0.36
b3 = -0.73
b4 = -0.2
b5 = 0.71
b6 = -0.19
b7 = 0.26
gmma = 0.4 #this is -ve in stuart, but i make it +ve for easier understanding. Treatment increases outcome.

P <- P %>% 
  mutate(y0_hte0 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + rnorm(n = n(), mean = 0, sd = .5),
         y1_hte0 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + gmma + rnorm(n = n(), mean = 0, sd = .5),
         y0_prop1_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop1 + rnorm(n = n(), mean = 0, sd = .5),
         y0_prop2_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop2 + rnorm(n = n(), mean = 0, sd = .5),
         y0_prop3_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop3 + rnorm(n = n(), mean = 0, sd = .5),
         y0_prop4_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop4 + rnorm(n = n(), mean = 0, sd = .5),
         y0_prop5_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop5 + rnorm(n = n(), mean = 0, sd = .5),
         y1_prop1_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop1 + gmma + rnorm(n = n(), mean = 0, sd = .5),
         y1_prop2_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop2 + gmma + rnorm(n = n(), mean = 0, sd = .5),
         y1_prop3_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop3 + gmma + rnorm(n = n(), mean = 0, sd = .5),
         y1_prop4_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop4 + gmma + rnorm(n = n(), mean = 0, sd = .5),
         y1_prop5_hte1 = b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X8 + b6*X9 + b7*X10 + prop5 + gmma + rnorm(n = n(), mean = 0, sd = .5),
         delta_hte0 = y1_hte0 - y0_hte0,
         delta_prop1_hte1 = y1_prop1_hte1 - y0_prop1_hte1,
         delta_prop2_hte1 = y1_prop2_hte1 - y0_prop2_hte1,
         delta_prop3_hte1 = y1_prop3_hte1 - y0_prop3_hte1,
         delta_prop4_hte1 = y1_prop4_hte1 - y0_prop4_hte1,
         delta_prop5_hte1 = y1_prop5_hte1 - y0_prop5_hte1)


# Main Simulation Loop ----------------------------------------------------
# Save results here
simulation_results <- NULL
saved_P <- NULL
  
  # Main for loop
  for (k in 1:n_sim) {
    # Randomly assign treatment with probability = propensity score, and then
    # accordingly set which of potential outcomes to set as observed outcome
    P <- P %>%
      mutate(
        di_prop1 = rbinom(n = n(), size = 1, prob = prop1),
        di_prop2 = rbinom(n = n(), size = 1, prob = prop2),
        di_prop3 = rbinom(n = n(), size = 1, prob = prop3),
        di_prop4 = rbinom(n = n(), size = 1, prob = prop4),
        di_prop5 = rbinom(n = n(), size = 1, prob = prop5),
        yi_prop1_hte0 = ifelse(di_prop1 == 1, y1_hte0, y0_hte0),
        yi_prop2_hte0 = ifelse(di_prop2 == 1, y1_hte0, y0_hte0),
        yi_prop3_hte0 = ifelse(di_prop3 == 1, y1_hte0, y0_hte0),
        yi_prop4_hte0 = ifelse(di_prop4 == 1, y1_hte0, y0_hte0),
        yi_prop5_hte0 = ifelse(di_prop5 == 1, y1_hte0, y0_hte0),
        yi_prop1_hte1 = ifelse(di_prop1 == 1, y1_prop1_hte1, y0_prop1_hte1),
        yi_prop2_hte1 = ifelse(di_prop2 == 1, y1_prop2_hte1, y0_prop2_hte1),
        yi_prop3_hte1 = ifelse(di_prop3 == 1, y1_prop3_hte1, y0_prop3_hte1),
        yi_prop4_hte1 = ifelse(di_prop4 == 1, y1_prop4_hte1, y0_prop4_hte1),
        yi_prop5_hte1 = ifelse(di_prop5 == 1, y1_prop5_hte1, y0_prop5_hte1))
    
    # Save Average Treatment Effect on Treated for analysis later
    scratch1 <- P %>% group_by(di_prop1) %>% summarise(att = mean(delta_hte0))
    scratch2 <- P %>% group_by(di_prop2) %>% summarise(att = mean(delta_hte0))
    scratch3 <- P %>% group_by(di_prop3) %>% summarise(att = mean(delta_hte0))
    scratch4 <- P %>% group_by(di_prop4) %>% summarise(att = mean(delta_hte0))
    scratch5 <- P %>% group_by(di_prop5) %>% summarise(att = mean(delta_hte0))
    scratch6 <- P %>% group_by(di_prop1) %>% summarise(att = mean(delta_prop1_hte1))
    scratch7 <- P %>% group_by(di_prop2) %>% summarise(att = mean(delta_prop2_hte1))
    scratch8 <- P %>% group_by(di_prop3) %>% summarise(att = mean(delta_prop3_hte1))
    scratch9 <- P %>% group_by(di_prop4) %>% summarise(att = mean(delta_prop4_hte1))
    scratch10 <- P %>% group_by(di_prop5) %>% summarise(att = mean(delta_prop5_hte1))
    
    att_prop1_hte0 <- scratch1$att[2]
    att_prop2_hte0 <- scratch2$att[2]
    att_prop3_hte0 <- scratch3$att[2]
    att_prop4_hte0 <- scratch4$att[2]
    att_prop5_hte0 <- scratch5$att[2]
    att_prop1_hte1 <- scratch6$att[2]
    att_prop2_hte1 <- scratch7$att[2]
    att_prop3_hte1 <- scratch8$att[2]
    att_prop4_hte1 <- scratch9$att[2]
    att_prop5_hte1 <- scratch10$att[2]
    
    # Different types of att in different scenarios
    list_att_hte0 <- c(att_prop1_hte0, att_prop2_hte0, att_prop3_hte0, att_prop4_hte0, att_prop5_hte0)
    list_att_hte1 <- c(att_prop1_hte1, att_prop2_hte1, att_prop3_hte1, att_prop4_hte1, att_prop5_hte1)
    list_att <- rep(c(rep(list_att_hte0, times=3), rep(list_att_hte1, times=3)), times=3)
    
    
    
  
    # Just Regression
    model_1_prop1_conf1_hte0 <- P %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop2_conf1_hte0 <- P %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop3_conf1_hte0 <- P %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop4_conf1_hte0 <- P %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop5_conf1_hte0 <- P %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    
    model_1_prop1_conf2_hte0 <- P %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>% #removed X1 and X4
      tidy()
    model_1_prop2_conf2_hte0 <- P %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop3_conf2_hte0 <- P %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop4_conf2_hte0 <- P %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop5_conf2_hte0 <- P %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    
    model_1_prop1_conf3_hte0 <- P %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>% #removed all the confounders X1,X2,X3,X4
      tidy()
    model_1_prop2_conf3_hte0 <- P %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop3_conf3_hte0 <- P %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop4_conf3_hte0 <- P %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop5_conf3_hte0 <- P %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    
    
    model_1_prop1_conf1_hte1 <- P %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop2_conf1_hte1 <- P %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop3_conf1_hte1 <- P %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop4_conf1_hte1 <- P %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop5_conf1_hte1 <- P %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X1 + X2 + X3 + X4 + +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    
    model_1_prop1_conf2_hte1 <- P %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>% #removed X1 and X4
      tidy()
    model_1_prop2_conf2_hte1 <- P %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop3_conf2_hte1 <- P %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop4_conf2_hte1 <- P %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop5_conf2_hte1 <- P %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X2 + X3 +X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    
    model_1_prop1_conf3_hte1 <- P %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>% #removed all the confounders X1,X2,X3,X4
      tidy()
    model_1_prop2_conf3_hte1 <- P %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop3_conf3_hte1 <- P %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop4_conf3_hte1 <- P %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()
    model_1_prop5_conf3_hte1 <- P %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X5 + X6 + X7 + X8 + X9 + X10, data = .) %>%
      tidy()

    # Regression on a balanced sample
    model_2_prop1_conf1_hte0 <- P %>%
      # Run matching
      matchit(di_prop1 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop2_conf1_hte0 <- P %>%
      # Run matching
      matchit(di_prop2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop3_conf1_hte0 <- P %>%
      # Run matching
      matchit(di_prop3 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop4_conf1_hte0 <- P %>%
      # Run matching
      matchit(di_prop4 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop5_conf1_hte0 <- P %>%
      # Run matching
      matchit(di_prop5 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    
    model_2_prop1_conf2_hte0 <- P %>%
      # Run matching
      matchit(di_prop1 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop2_conf2_hte0 <- P %>%
      # Run matching
      matchit(di_prop2 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop3_conf2_hte0 <- P %>%
      # Run matching
      matchit(di_prop3 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop4_conf2_hte0 <- P %>%
      # Run matching
      matchit(di_prop4 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop5_conf2_hte0 <- P %>%
      # Run matching
      matchit(di_prop5 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    
    model_2_prop1_conf3_hte0 <- P %>%
      # Run matching
      matchit(di_prop1 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop2_conf3_hte0 <- P %>%
      # Run matching
      matchit(di_prop2 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop3_conf3_hte0 <- P %>%
      # Run matching
      matchit(di_prop3 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop4_conf3_hte0 <- P %>%
      # Run matching
      matchit(di_prop4 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop5_conf3_hte0 <- P %>%
      # Run matching
      matchit(di_prop5 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    
    
    model_2_prop1_conf1_hte1 <- P %>%
      # Run matching
      matchit(di_prop1 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop2_conf1_hte1 <- P %>%
      # Run matching
      matchit(di_prop2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop3_conf1_hte1 <- P %>%
      # Run matching
      matchit(di_prop3 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop4_conf1_hte1 <- P %>%
      # Run matching
      matchit(di_prop4 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop5_conf1_hte1 <- P %>%
      # Run matching
      matchit(di_prop5 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    
    model_2_prop1_conf2_hte1 <- P %>%
      # Run matching
      matchit(di_prop1 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop2_conf2_hte1 <- P %>%
      # Run matching
      matchit(di_prop2 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop3_conf2_hte1 <- P %>%
      # Run matching
      matchit(di_prop3 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop4_conf2_hte1 <- P %>%
      # Run matching
      matchit(di_prop4 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop5_conf2_hte1 <- P %>%
      # Run matching
      matchit(di_prop5 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    
    model_2_prop1_conf3_hte1 <- P %>%
      # Run matching
      matchit(di_prop1 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop2_conf3_hte1 <- P %>%
      # Run matching
      matchit(di_prop2 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop3_conf3_hte1 <- P %>%
      # Run matching
      matchit(di_prop3 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop4_conf3_hte1 <- P %>%
      # Run matching
      matchit(di_prop4 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    model_2_prop5_conf3_hte1 <- P %>%
      # Run matching
      matchit(di_prop5 ~ X5 + X6 + X7 + X8 + X9 + X10, method = "nearest", replace = TRUE, data = .) %>%
      match.data(., weight = "weights") %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X5 + X6 + X7 + X8 + X9 + X10, weights = weights, data = .) %>%
      tidy()
    
    # Estimation 3
    model_3_prop1_conf1_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop1 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop1))/(1-pi)) %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop2_conf1_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop2))/(1-pi)) %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop3_conf1_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop3 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop3 + (pi*(1-di_prop3))/(1-pi)) %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop4_conf1_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop4 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop4 + (pi*(1-di_prop4))/(1-pi)) %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop5_conf1_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop5 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop5 + (pi*(1-di_prop5))/(1-pi)) %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    
    model_3_prop1_conf2_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop1 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop1))/(1-pi)) %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop2_conf2_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop2 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop2))/(1-pi)) %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop3_conf2_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop3 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop3 + (pi*(1-di_prop3))/(1-pi)) %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop4_conf2_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop4 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop4 + (pi*(1-di_prop4))/(1-pi)) %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop5_conf2_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop5 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop5 + (pi*(1-di_prop5))/(1-pi)) %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    
    model_3_prop1_conf3_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop1 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop1))/(1-pi)) %>%
      lm(yi_prop1_hte0 ~ di_prop1 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop2_conf3_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop2 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop2))/(1-pi)) %>%
      lm(yi_prop2_hte0 ~ di_prop2 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop3_conf3_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop3 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop3 + (pi*(1-di_prop3))/(1-pi)) %>%
      lm(yi_prop3_hte0 ~ di_prop3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop4_conf3_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop4 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop4 + (pi*(1-di_prop4))/(1-pi)) %>%
      lm(yi_prop4_hte0 ~ di_prop4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop5_conf3_hte0 <- P %>%
      # Compute fitted probability of
      glm(di_prop5 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop5 + (pi*(1-di_prop5))/(1-pi)) %>%
      lm(yi_prop5_hte0 ~ di_prop5 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    
    
    model_3_prop1_conf1_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop1 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop1))/(1-pi)) %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop2_conf1_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop2))/(1-pi)) %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop3_conf1_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop3 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop3 + (pi*(1-di_prop3))/(1-pi)) %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop4_conf1_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop4 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop4 + (pi*(1-di_prop4))/(1-pi)) %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop5_conf1_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop5 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop5 + (pi*(1-di_prop5))/(1-pi)) %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    
    model_3_prop1_conf2_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop1 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop1))/(1-pi)) %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop2_conf2_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop2 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop2))/(1-pi)) %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop3_conf2_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop3 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop3 + (pi*(1-di_prop3))/(1-pi)) %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop4_conf2_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop4 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop4 + (pi*(1-di_prop4))/(1-pi)) %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop5_conf2_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop5 ~ X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop5 + (pi*(1-di_prop5))/(1-pi)) %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    
    model_3_prop1_conf3_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop1 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop1))/(1-pi)) %>%
      lm(yi_prop1_hte1 ~ di_prop1 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop2_conf3_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop2 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop1 + (pi*(1-di_prop2))/(1-pi)) %>%
      lm(yi_prop2_hte1 ~ di_prop2 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop3_conf3_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop3 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop3 + (pi*(1-di_prop3))/(1-pi)) %>%
      lm(yi_prop3_hte1 ~ di_prop3 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop4_conf3_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop4 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop4 + (pi*(1-di_prop4))/(1-pi)) %>%
      lm(yi_prop4_hte1 ~ di_prop4 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    model_3_prop5_conf3_hte1 <- P %>%
      # Compute fitted probability of
      glm(di_prop5 ~ X5 + X6 + X7 + X8 + X9 + X10, family = binomial, data = .) %>%
      predict(type = "response") %>%
      tbl_df() %>%
      rename(pi = value) %>%
      bind_cols(P) %>%
      # Compute weights
      mutate(wi = di_prop5 + (pi*(1-di_prop5))/(1-pi)) %>%
      lm(yi_prop5_hte1 ~ di_prop5 + X5 + X6 + X7 + X8 + X9 + X10, weights = wi, data = .) %>%
      tidy()
    
    
    # Save simulation results for this k
    ID_variables <- data_frame(
      iteration = k,
      model = model_names
#      ATT = att
    )
    beta_di <-
      bind_rows(
        model_1_prop1_conf1_hte0, model_1_prop2_conf1_hte0, model_1_prop3_conf1_hte0, model_1_prop4_conf1_hte0, model_1_prop5_conf1_hte0,
        model_1_prop1_conf2_hte0, model_1_prop2_conf2_hte0, model_1_prop3_conf2_hte0, model_1_prop4_conf2_hte0, model_1_prop5_conf2_hte0,
        model_1_prop1_conf3_hte0, model_1_prop2_conf3_hte0, model_1_prop3_conf3_hte0, model_1_prop4_conf3_hte0, model_1_prop5_conf3_hte0,
        model_1_prop1_conf1_hte1, model_1_prop2_conf1_hte1, model_1_prop3_conf1_hte1, model_1_prop4_conf1_hte1, model_1_prop5_conf1_hte1,
        model_1_prop1_conf2_hte1, model_1_prop2_conf2_hte1, model_1_prop3_conf2_hte1, model_1_prop4_conf2_hte1, model_1_prop5_conf2_hte1,
        model_1_prop1_conf3_hte1, model_1_prop2_conf3_hte1, model_1_prop3_conf3_hte1, model_1_prop4_conf3_hte1, model_1_prop5_conf3_hte1,
        
        model_2_prop1_conf1_hte0, model_2_prop2_conf1_hte0, model_2_prop3_conf1_hte0, model_2_prop4_conf1_hte0, model_2_prop5_conf1_hte0,
        model_2_prop1_conf2_hte0, model_2_prop2_conf2_hte0, model_2_prop3_conf2_hte0, model_2_prop4_conf2_hte0, model_2_prop5_conf2_hte0,
        model_2_prop1_conf3_hte0, model_2_prop2_conf3_hte0, model_2_prop3_conf3_hte0, model_2_prop4_conf3_hte0, model_2_prop5_conf3_hte0,
        model_2_prop1_conf1_hte1, model_2_prop2_conf1_hte1, model_2_prop3_conf1_hte1, model_2_prop4_conf1_hte1, model_2_prop5_conf1_hte1,
        model_2_prop1_conf2_hte1, model_2_prop2_conf2_hte1, model_2_prop3_conf2_hte1, model_2_prop4_conf2_hte1, model_2_prop5_conf2_hte1,
        model_2_prop1_conf3_hte1, model_2_prop2_conf3_hte1, model_2_prop3_conf3_hte1, model_2_prop4_conf3_hte1, model_2_prop5_conf3_hte1,
        
        model_3_prop1_conf1_hte0, model_3_prop2_conf1_hte0, model_3_prop3_conf1_hte0, model_3_prop4_conf1_hte0, model_3_prop5_conf1_hte0,
        model_3_prop1_conf2_hte0, model_3_prop2_conf2_hte0, model_3_prop3_conf2_hte0, model_3_prop4_conf2_hte0, model_3_prop5_conf2_hte0,
        model_3_prop1_conf3_hte0, model_3_prop2_conf3_hte0, model_3_prop3_conf3_hte0, model_3_prop4_conf3_hte0, model_3_prop5_conf3_hte0,
        model_3_prop1_conf1_hte1, model_3_prop2_conf1_hte1, model_3_prop3_conf1_hte1, model_3_prop4_conf1_hte1, model_3_prop5_conf1_hte1,
        model_3_prop1_conf2_hte1, model_3_prop2_conf2_hte1, model_3_prop3_conf2_hte1, model_3_prop4_conf2_hte1, model_3_prop5_conf2_hte1,
        model_3_prop1_conf3_hte1, model_3_prop2_conf3_hte1, model_3_prop3_conf3_hte1, model_3_prop4_conf3_hte1, model_3_prop5_conf3_hte1) %>%
#      filter(term == "di") %>%
      filter(grepl("di",term)) %>% 
      select(estimate, std.error)
    
    list_att <- data.frame(list_att)
    # Append to overall results
    simulation_results <- simulation_results %>%
      bind_rows(bind_cols(ID_variables, beta_di, list_att))
    
    # Show simulation counter and save results
    if(k %% 10 == 1){
      str_c("Iter", k, "of", n_sim, "done at", Sys.time(), sep=" ") %>%
        print()
      save(simulation_results, file="simulation_results.RData")
      save(saved_P, file="saved_P.RData")
    }
  }

save(simulation_results, file="simulation_results.RData")
save(saved_P, file="saved_P.RData")

# Above (for n_sim=1000) took X time to run on Amherst College RStudio server






# Analyze Results ---------------------------------------------------------
# Compute average treatment effect for both simulations
# ATE <- bind_rows(
#   P_1 %>%
#     summarise(ATE = mean(delta)) %>%
#     mutate(simulation = 1),
#   P_2 %>%
#     summarise(ATE = mean(delta)) %>%
#     mutate(simulation = 2)
# )

# Load saved simulation results and compute bias
load("simulation_results.RData")
# Delete later:
# n_sim <- simulation_results$iteration %>% max



# Bert's code to write simulation_results data frame to results_for_figure.csv
# so that format matches Placeholder_Results.xlsx
simulation_results_new <- simulation_results %>% 
  # Separate model variable
  separate(model, c("method", "PropModel", "confounding", "HTE"), sep = "_") %>% 
  # Reformat data
  mutate(
    prop_hte_facet = str_c(PropModel, HTE), 
    method = str_c("method", method),
    HTE = str_sub(HTE, 4, 4) %>% as.integer(),
    confounding = str_sub(confounding, 5, 5) %>% as.integer() %>% `-`(1) %>% as.integer(),
    #
    # Not sure about these:
    #
    bias = estimate - 76,
    variance = std.error^2
  ) %>% 
  # Match ordering of Placeholder_Results.xlsx
  select(prop_hte_facet, PropModel, HTE, confounding, method, bias, variance) %>% 
  arrange(PropModel, desc(HTE), method) %>% 
  # Make column names match
  rename(
    `PropHTE (facet)` = prop_hte_facet,
    `Strength of confounding (x axis)` = confounding,
    `Method (line)` = method,
    `Bias (y axis)` = bias,
    `Variance (y axis)` = variance
  )

# Write to CSV
write_csv(simulation_results_new, path = "results_for_figure.csv")




# ORIGINAL CODE
simulation_results <- simulation_results %>%
  mutate(
    method = str_sub(grepl("model",.), start=1, end=1),
    method = factor(method),
    method = fct_recode(method,
                        "Just Regression" = "1",
                        "Matching" = "2",
                        "Inverse Prop Score Weights" = "3"
    )
  ) %>%
  #left_join(ATE, by="simulation") %>%
  mutate(bias = estimate - ATT) %>%
  select(model, method, bias, std.error)


# Plot
# simulation_results %>%
#   gather(type, value, -c(model, method, simulation)) %>%
#   ggplot() +
#   geom_boxplot(aes(x=model, y=value, fill=method)) +
#   facet_grid(type~simulation, scales = "free") +
#   labs(title = str_c(n_sim, " iterations"), y = "") +
#   scale_fill_brewer(palette = "Set2")

# Table
simulation_results %>%
  group_by(simulation, model) %>%
  summarise(bias=mean(bias), std.error = mean(std.error)) %>%
  kable(digits=4)
