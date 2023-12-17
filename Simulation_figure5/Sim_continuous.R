# project: ADCC
# Author: Jiawei Zhou
# Objectives: simulate R and S and R+S under different dosing regimen
# Updated Model

rm(list=ls())


library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(deSolve)

#write a loop for 10 cycles of treatment
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

# Loop through simulation 10 times
for (i in 1:40) {
  #set up parameters
  parms <- c(g = 0.4, mu_SR = 0.347, ds = 5, dr = 4.26, k=2.67, w=0.028, n=i)
  # Set up initial conditions
  if (i == 1) {
    xstart <- c(S = f*N, R = (1-f)*N)
  } else {
    xstart <- c(S = out$S[length(out$S)], R = out$R[length(out$R)])
  }
  
  # Run ODE model
  ode(
    func = trt.model,
    y = xstart,
    times = times,
    parms = parms
  ) %>% as.data.frame() -> out 
  
  # Add simulation results to the data frame
  out2 <- out %>% mutate(Total = S + R, iteration = i)
  results <- rbind(results, out2)
}

results2 <- results %>%
  mutate(time2 = 3*(iteration-1)+time)%>%
  mutate(S =  ifelse(S<=1, 1, S)) 


#plot
p1 <- results2 %>%gather(variable,value,-time2, -iteration, -time, -Total) %>%
  ggplot(aes(x=time2,y=value,color=variable))+
  geom_line(size=1.5)+
  theme_classic()+
  scale_y_continuous(trans='log10',limits=c(1e0, 1e20))+
  labs(x='Time (day)',y='Cell number')+ggtitle('Continuous Dosing')


ggsave(p1, file='result/continuous.tiff', width=4, height=2)

write.csv(results2, file='result/continuous.csv', quote=FALSE, row.names=FALSE)
