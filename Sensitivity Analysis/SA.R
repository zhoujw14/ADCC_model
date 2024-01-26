# project: ADCC
# Author: Jiawei Zhou
# Objectives: perform sensitivity analysis to evaluate model uncertainty
# Figure S6

rm(list=ls())

library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(deSolve)

#simulate contiuoue dosing, 3 days a cycle
#cell dynamics under treatment
trt.model <- function(t,x,params){
  # first extract state variables
  S <- x[1]
  R <- x[2]
  
  ## now extract parameters
  g <- params['g'] #cell growth rate
  mu_SR <- params['mu_SR'] #S to R transition rate
  ds <- params['ds'] # S death rate under drug
  dr <- params['dr'] #R death rate under drug
  w  <- params['w'] #R death rate decay rate
  k <- params['k'] #drug efficacy decay rate
  n <- params['n'] # dosing cycle
  
  n_time <- (n-1)*3+t
  
  ## now code the model
  dSdt <- g*S - mu_SR*S -ds*exp(-k*t)*S
  dRdt <- g*R + mu_SR*S -dr*exp(-w*n_time)*exp(-k*t)*R
  
  ## combine results into a vector
  dxdt <- c(dSdt, dRdt)
  
  ## return result as a list
  list(dxdt)
}


# Set up initial conditions
N <- 5000 # initial total cell number
f <- 0.239 # initial S fraction

# Set up time points
times <- seq(from = 0, to = 3, by = 1) #sim 3 days

# Set up empty data frame to store results
results <- data.frame()

# Sensitivity step
sa.step <- c(0.8, 0.9, 1, 1.1, 1.2)

# mu_SR sensitivity
for (i in sa.step) {
  #set up parameters
  parms <- c(g = 0.4, mu_SR = 0.347*i, ds = 5, dr = 4.26, k=2.67, w=0.028, n=i)
  # the first cycle of treatment
    xstart <- c(S = f*N, R = (1-f)*N)
  
  # Run ODE model
  ode(
    func = trt.model,
    y = xstart,
    times = times,
    parms = parms
  ) %>% as.data.frame() -> out 
  
  # Add simulation results to the data frame
  out2 <- out %>% mutate(Total = S + R, SA = i)
  results <- rbind(results, out2)
}

### Plot mu_SR parameter sensitivity
p1 <- results %>% mutate(sensitivity = factor(SA)) %>%
  ggplot(aes(x=time,y=Total,color=sensitivity))+
  geom_line(size=1.5)+
  theme_classic()+
  scale_y_continuous(trans='log10')+
  labs(x='Time (day)',y='Total Cell number')+ggtitle(expression("Sensitivity on "~mu*"_SR"))


# ds sensitivity
for (i in sa.step) {
  #set up parameters
  parms <- c(g = 0.4, mu_SR = 0.347, ds = 5*i, dr = 4.26, k=2.67, w=0.028, n=i)
  # the first cycle of treatment
  xstart <- c(S = f*N, R = (1-f)*N)
  
  # Run ODE model
  ode(
    func = trt.model,
    y = xstart,
    times = times,
    parms = parms
  ) %>% as.data.frame() -> out 
  
  # Add simulation results to the data frame
  out2 <- out %>% mutate(Total = S + R, SA = i)
  results <- rbind(results, out2)
}

### Plot ds parameter sensitivity
p2 <- results %>% mutate(sensitivity = factor(SA)) %>%
  ggplot(aes(x=time,y=Total,color=sensitivity))+
  geom_line(size=1.5)+
  theme_classic()+
  scale_y_continuous(trans='log10')+
  labs(x='Time (day)',y='Total Cell number')+ggtitle("Sensitivity on ds")


# merge two plots
library(ggpubr)
p3 <- ggarrange(plotlist=list(p1, p2))
ggsave(p3, file='sensitivity.tiff', width=8, height=4)
