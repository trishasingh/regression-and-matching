# Final Simulations with 1000 reps

# Setting up
setwd("~/Documents/Stats Indep Study/Simulations")
library(tidyr)
library(stringr)
library(dplyr)
library(broom)
library(Matching)
library(MatchIt)
library(ggplot2)

# Creating a and b, independent variables
A <- seq(.01,1,.01)
B <- seq(.01,1,.01)

# Simulation 1.1:
p1 <- function(a,b) {
  x <- a + a^2 + b
  return(1/(1+exp(-x)))
}

#-------------------------------------------------------------------------------------------------------------------------------
# Propensity score matrix for 10000 combinations of A and B
Pmat <- matrix(nrow=100, ncol=100)
for (i in 1:100) {
  for (j in 1:100) {
    Pmat[i,j] <- p1(i/100, j/100)
  }
}
# Converting prop score matrix into a dataframe with A and B variables for each person
rownames(Pmat) <- A
colnames(Pmat) <- B
P <- data.frame(Pmat)
P <- gather(P, "B", "PScore", 1:100) %>% 
  mutate(B = as.numeric(str_sub(B, 2, 5)))
P <- cbind(A, P)

# Generate potential outcomes
P <- P %>% 
  mutate(y0 = 100 + 3*A + 2*B + rnorm(1,0,5)) %>% 
  mutate(y1 = 102 + 6*A + 4*B + rnorm(1,0,5)) %>% 
  mutate(te = y1-y0)

ate1 <- mean(P$te)

# Treatment assignment function (treatment is bernoulli distributed wp propensity score)
treatment <- function(prob) {
  return(rbinom(1,1,prob))
}

pscore_sim1 <- ggplot(P, aes(x=A, y=B)) +
  geom_tile(aes(fill=PScore)) + 
  labs(fill="Propensity Score") +
  ggtitle("Figure 1: Propensity Score Function for Simulation 1") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# Start loop here

b_1_1 <- numeric(1000)
b_1_2 <- numeric(1000)
b_1_3 <- numeric(1000)
se_1_1 <- numeric(1000)
se_1_2 <- numeric(1000)
se_1_3 <- numeric(1000)

b_2_1 <- numeric(1000)
b_2_2 <- numeric(1000)
se_2_1 <- numeric(1000)
se_2_2 <- numeric(1000)

b_3_1 <- numeric(1000)
b_3_2 <- numeric(1000)
se_3_1 <- numeric(1000)
se_3_2 <- numeric(1000)


for (k in 1:1000) {
  # Creating matrix with binomial treatment assignment
  Dmat <- matrix(nrow=100, ncol=100)
  for (i in 1:100) {
    for (j in 1:100) {
      Dmat[i,j] <- rbinom(1,1,Pmat[i,j])
    }
  }
  Dmat <- data.frame(Dmat)
  Dmat <- gather(Dmat, "B", "D", 1:100)
  
  P$D <- NULL
  P$y <- NULL
  
  P <- P %>% 
    arrange(A,B) %>% 
    cbind2(Dmat$D, P) %>%
    rename(D = y) %>% 
    tbl_df()
  
  
  #-----------DATASET P GENERATED------------------------------------------------------------------------------------------------------------
  
  # Generate potential outcomes
  P <- P %>% 
    mutate(yi = ifelse(D==1, y1, y0))
  
  #Estimation 1.1
  model_1_1 <- lm(yi~D, data=P)
  b_1_1[k] <- coef(model_1_1)["D"]
  se_1_1[k] <- coef(summary(model_1_1))[, 2]["D"]
  # Estimation 1.2
  model_1_2 <- lm(yi~D+B, data=P)
  b_1_2[k] <- coef(model_1_2)["D"]
  se_1_2[k] <- coef(summary(model_1_2))[, 2]["D"]
  # Estimation 1.3
  model_1_3 <- lm(yi~D+A+B, data=P)
  tidy(model_1_3)
  b_1_3[k] <- coef(model_1_3)["D"]
  se_1_3[k] <- coef(summary(model_1_3))[, 2]["D"]
  
  #Estimation 2.1
  match_2_1 <- matchit(D ~ B, method = "nearest", replace = TRUE, data = P)
  dta_2_1 <- match.data(match_2_1)
  model_2_1 <- lm(yi ~ D+B, data = dta_2_1)
  b_2_1[k] <- coef(model_2_1)["D"]
  se_2_1[k] <- coef(summary(model_2_1))[, 2]["D"]
  #Estimation 2.2
  match_2_2 <- matchit(D ~ A+B, method = "nearest", replace = TRUE, data = P)
  dta_2_2 <- match.data(match_2_2)
  model_2_2 <- lm(yi ~ D+A+B, data = dta_2_2)
  b_2_2[k] <- coef(model_2_2)["D"]
  se_2_2[k] <- coef(summary(model_2_2))[, 2]["D"]
  
  #Estimation 3.1
  prop1_3_1 <- glm(D ~ B, family = binomial(link="logit"), data = P)
  P <- P %>% 
    mutate(pi=predict(prop1_3_1, type="response")) %>% 
    mutate(ip = D/pi + (1-D)/(1-pi))
  model_3_1 <- lm(yi ~ D+B, data = P, weights = ip)
  b_3_1[k] <- coef(model_3_1)["D"]
  se_3_1[k] <- coef(summary(model_3_1))[, 2]["D"]
  #Estimation 3.2
  prop1_3_2 <- glm(D ~ A+B, family = binomial(link="logit"), data = P)
  P <- P %>% 
    mutate(pi=predict(prop1_3_2, type="response")) %>% 
    mutate(ip = D/pi + (1-D)/(1-pi))
  model_3_2 <- lm(yi ~ D+A+B, data = P, weights = ip)
  b_3_2[k] <- coef(model_3_2)["D"]
  se_3_2[k] <- coef(summary(model_3_2))[, 2]["D"]
}


################
# SIMULATION 2 #
################

p2 <- function(a,b) {
  x <- a + b
  return(1/(1+exp(-x)))
}

#-------------------------------------------------------------------------------------------------------------------------------
# Propensity score matrix for 10000 combinations of A and B
Pmat2 <- matrix(nrow=100, ncol=100)
for (i in 1:100) {
  for (j in 1:100) {
    Pmat2[i,j] <- p2(i/100, j/100)
  }
}
# Converting prop score matrix into a dataframe with A and B variables for each person
rownames(Pmat2) <- A
colnames(Pmat2) <- B
P2 <- data.frame(Pmat2)
P2 <- gather(P2, "B", "PScore", 1:100) %>% 
  mutate(B = as.numeric(str_sub(B, 2, 5)))
P2 <- cbind(A, P2)

# Generate potential outcomes
P2 <- P2 %>% 
  mutate(y0 = 100 + 3*A + 2*B + rnorm(1,0,5)) %>% 
  mutate(y1 = 102 + 6*A + 4*B + 5*PScore + rnorm(1,0,5)) %>% 
  mutate(te = y1-y0)

ate2 <- mean(P2$te)

# Treatment assignment function (treatment is bernoulli distributed wp propensity score)
treatment <- function(prob) {
  return(rbinom(1,1,prob))
}

pscore_sim2 <- ggplot(P2, aes(x=A, y=B)) +
  geom_tile(aes(fill=PScore)) + 
  labs(fill="Propensity Score") +
  ggtitle("Figure 2: Propensity Score Function for Simulation 2") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))



# Start loop here

b2_1_1 <- numeric(1000)
b2_1_2 <- numeric(1000)
b2_1_3 <- numeric(1000)
se2_1_1 <- numeric(1000)
se2_1_2 <- numeric(1000)
se2_1_3 <- numeric(1000)

b2_2_1 <- numeric(1000)
b2_2_2 <- numeric(1000)
se2_2_1 <- numeric(1000)
se2_2_2 <- numeric(1000)

b2_3_1 <- numeric(1000)
b2_3_2 <- numeric(1000)
se2_3_1 <- numeric(1000)
se2_3_2 <- numeric(1000)


for (k in 1:1000) {
  # Creating matrix with binomial treatment assignment
  Dmat2 <- matrix(nrow=100, ncol=100)
  for (i in 1:100) {
    for (j in 1:100) {
      Dmat2[i,j] <- rbinom(1,1,Pmat2[i,j])
    }
  }
  Dmat2 <- data.frame(Dmat2)
  Dmat2 <- gather(Dmat2, "B", "D", 1:100)
  
  P2$D <- NULL
  P2$y <- NULL
  
  P2 <- P2 %>% 
    arrange(A,B) %>% 
    cbind2(Dmat2$D, P2) %>%
    rename(D = y) %>% 
    tbl_df()
  
  # Generate potential outcomes
  P2 <- P2 %>% 
    mutate(yi = ifelse(D==1, y1, y0))
  
  #Estimation 1.1
  model2_1_1 <- lm(yi~D, data=P2)
  b2_1_1[k] <- coef(model2_1_1)["D"]
  se2_1_1[k] <- coef(summary(model2_1_1))[, 2]["D"]
  # Estimation 1.2
  model2_1_2 <- lm(yi~D+B, data=P2)
  b2_1_2[k] <- coef(model2_1_2)["D"]
  se2_1_2[k] <- coef(summary(model2_1_2))[, 2]["D"]
  # Estimation 1.3
  model2_1_3 <- lm(yi~D+A+B, data=P2)
  tidy(model2_1_3)
  b2_1_3[k] <- coef(model2_1_3)["D"]
  se2_1_3[k] <- coef(summary(model2_1_3))[, 2]["D"]
  
  #Estimation 2.1
  match2_2_1 <- matchit(D ~ B, method = "nearest", replace = TRUE, data = P2)
  dta2_2_1 <- match.data(match2_2_1)
  model2_2_1 <- lm(yi ~ D+B, data = dta2_2_1)
  b2_2_1[k] <- coef(model2_2_1)["D"]
  se2_2_1[k] <- coef(summary(model2_2_1))[, 2]["D"]
  #Estimation 2.2
  match2_2_2 <- matchit(D ~ A+B, method = "nearest", replace = TRUE, data = P2)
  dta2_2_2 <- match.data(match2_2_2)
  model2_2_2 <- lm(yi ~ D+A+B, data = dta2_2_2)
  b2_2_2[k] <- coef(model2_2_2)["D"]
  se2_2_2[k] <- coef(summary(model2_2_2))[, 2]["D"]
  
  #Estimation 3.1
  prop2_3_1 <- glm(D ~ B, family = binomial(link="logit"), data = P2)
  P2 <- P2 %>% 
    mutate(pi=predict(prop2_3_1, type="response")) %>% 
    mutate(ip = D/pi + (1-D)/(1-pi))
  model2_3_1 <- lm(yi ~ D+B, data = P2, weights = ip)
  b2_3_1[k] <- coef(model2_3_1)["D"]
  se2_3_1[k] <- coef(summary(model2_3_1))[, 2]["D"]
  #Estimation 3.2
  prop2_3_2 <- glm(D ~ A+B, family = binomial(link="logit"), data = P2)
  P2 <- P2 %>% 
    mutate(pi=predict(prop2_3_1, type="response")) %>% 
    mutate(ip = D/pi + (1-D)/(1-pi))
  model2_3_2 <- lm(yi ~ D+A+B, data = P2, weights = ip)
  b2_3_2[k] <- coef(model2_3_2)["D"]
  se2_3_2[k] <- coef(summary(model2_3_2))[, 2]["D"]
}