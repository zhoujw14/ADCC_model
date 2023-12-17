
# project: ADCC
# Author: Jiawei Zhou
# Objectives: figure 6
# Updated model

rm(list=ls())



library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(deSolve)

# the fold to evaluate
fold <- c(0.0001, 0.001, 0.01, 0.1, 0.5,1, 5,10, 100, 1000, 10000, 100000)

# continuous treatment model
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


total.con <- c()
for (kk in 1:length(fold)){
  # Set up initial conditions
  N <- 5000 # initial total cell number
  f <- 0.239 # initial S fraction
  tmp.mu_RS <- 0.0014*fold[kk]
  # Set up time points
  times <- seq(from = 0, to = 3, by = 1) #sim 3 days
  
  # Set up empty data frame to store results
  results <- data.frame()
  
  # Loop through simulation 40 times
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
    mutate(time2 = 3*(iteration-1)+time)
  
  total.con[kk] <- results2[results2$time2 == 100,]$Total
}

con <- data.frame(
  fold_w =  fold,
  total = total.con
)


# generate on 3 days off 3 days

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
  
  ## now code the model
  dSdt <- g*S - mu_SR*S -ds*exp(-k*t)*S
  dRdt <- g*R + mu_SR*S -dr*exp(-w*t)*exp(-k*t)*R
  
  ## combine results into a vector
  dxdt <- c(dSdt, dRdt)
  
  ## return result as a list
  list(dxdt)
}


rest.model <- function(t,x,params){
  # first extract state variables
  S <- x[1]
  R <- x[2]
  
  ## now extract parameters
  g <- params['g'] #cell growth rate
  mu_RS <- params['mu_RS'] #S to R transition rate
  
  ## now code the model
  dSdt <- g*S + mu_RS*R
  dRdt <- g*R - mu_RS*R
  
  ## combine results into a vector
  dxdt <- c(dSdt, dRdt)
  
  ## return result as a list
  list(dxdt)
}

# params
total.on3off3 <- c()

for(kk in 1:length(fold)){
  
  tmp.muRS <- fold[kk]*0.0014
  parms <- c(g = 0.4, mu_SR = 0.347, mu_RS = tmp.muRS, ds = 5, dr = 4.26, k=2.67, w=0.028)
  N <- 5000
  f <- 0.239
  # Set up empty data frame to store results
  out <- data.frame()
  # 20 cycles of treatment
  for (i in 1:20){
    # Set up initial conditions
    if (i == 1) {
      xstart <- c(S = f*N, R = (1-f)*N)
    } else {
      xstart <- c(S = out$S[length(out$S)], R = out$R[length(out$R)])
    }
    
    # Simulate first 3 days with trt.model
    times1 <- seq(from = 0, to = 3, by = 1)
    out1 <- ode(func = trt.model, y = xstart, times = times1, parms = parms) %>% as.data.frame()
    # Get initial values for rest.model from last day of trt.model simulation
    xstart2 <-  c(S = out1$S[length(out1$S)], R = out1$R[length(out1$R)])
    # Simulate last 3 days with rest.model
    times2 <- seq(from = 3, to = 6, by = 1)
    out2 <- ode(func = rest.model, y = xstart2, times = times2, parms = parms) %>% as.data.frame()
    
    # Combine results from both simulations and add to output data frame
    out_cycle <- rbind(out1, out2)
    out_cycle$cycle <- i
    out <- rbind(out, out_cycle)
    
    
  }
  
  results2 <- out %>%
    mutate(time2 = 6*(cycle-1)+time,
           N = S + R)
  
  total.on3off3[kk] <- results2[results2$time2==100,]$N
}

on3off3 <- data.frame(
  fold_w =  fold,
  total = total.on3off3
)
#write.csv(on3off3, file='on3off3_muRS.csv', quote=FALSE, row.names=FALSE)
p1 <- ggplot(con)+
  geom_line(aes(x=fold_w, y=total, col='Continuous'), size=1.5)+
  geom_line(data=on3off3, aes(x=fold_w, y=total, col='ON 3 days + OFF 3 days'), size=1.5)+
  scale_colour_manual("Dosing Strategy", breaks = c("Continuous", 
                                                    "ON 3 days + OFF 3 days"),
                      values = c("Continuous" = "red", 
                                 "ON 3 days + OFF 3 days" = "dodgerblue4"))+
  theme_classic()+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  labs(x=expression("Fold on "~mu*"_RS"),y='Total cell number', title=expression("Reversibility "~mu*"_RS"))
ggsave(p1, file='result/muRS_on3off3_con.tiff', width=5, height=2)


# on 3 off 10

#simulate 3 days on, 10 days off
# Set up empty data frame to store results
total.on3off10 <- c()

for(kk in 1:length(fold)){
  
  tmp.muRS <- fold[kk]*0.0014
  parms <- c(g = 0.4, mu_SR = 0.347, mu_RS = tmp.muRS, ds = 5, dr = 4.26, k=2.67, w=0.028)
  N <- 5000
  f <- 0.239
  out <- data.frame()
  # 5 cycles of treatment
  for (i in 1:10){
    # Set up initial conditions
    if (i == 1) {
      xstart <- c(S = f*N, R = (1-f)*N)
    } else {
      xstart <- c(S = out$S[length(out$S)], R = out$R[length(out$R)])
    }
    
    # Simulate first 3 days with trt.model
    times1 <- seq(from = 0, to = 3, by = 1)
    out1 <- ode(func = trt.model, y = xstart, times = times1, parms = parms) %>% as.data.frame()
    # Get initial values for rest.model from last day of trt.model simulation
    xstart2 <-  c(S = out1$S[length(out1$S)], R = out1$R[length(out1$R)])
    # Simulate last 5 days with rest.model
    times2 <- seq(from = 3, to = 13, by = 1)
    out2 <- ode(func = rest.model, y = xstart2, times = times2, parms = parms) %>% as.data.frame()
    
    # Combine results from both simulations and add to output data frame
    out_cycle <- rbind(out1, out2)
    out_cycle$cycle <- i
    out <- rbind(out, out_cycle)
    
    
  }
  
  results2 <- out %>%
    mutate(time2 = 13*(cycle-1)+time, N=S+R)
  total.on3off10[kk] <- results2[results2$time2==100,]$N
}

on3off10 <- data.frame(
  fold_w =  fold,
  total = total.on3off10
)


p2 <- ggplot(on3off3)+
  geom_line(aes(x=fold_w, y=total, col='ON 3 days + OFF 3 days'), size=1.5)+
  geom_line(data=on3off10, aes(x=fold_w, y=total, col='ON 3 days + OFF 10 days'), size=1.5)+
  scale_colour_manual("Dosing Strategy", breaks = c("ON 3 days + OFF 3 days", 
                                                    "ON 3 days + OFF 10 days"),
                      values = c("ON 3 days + OFF 3 days" = "dodgerblue4", 
                                 "ON 3 days + OFF 10 days" = "deepskyblue1"))+
  theme_classic()+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  labs(x=expression("Fold on "~mu*"_RS"),y='Total cell number', title=expression("Reversibility "~mu*"_RS"))
ggsave(p2, file='result/muRS_on3off10_on3off3.tiff', width=5, height=2)

p3 <- ggplot(con)+
  geom_line(aes(x=fold_w, y=total, col='Continuous'), size=1.5)+
  geom_line(data=on3off10, aes(x=fold_w, y=total, col='ON 3 days + OFF 10 days'), size=1.5)+
  scale_colour_manual("Dosing Strategy", breaks = c("Continuous", 
                                                    "ON 3 days + OFF 10 days"),
                      values = c("Continuous" = "red", 
                                 "ON 3 days + OFF 10 days" = "deepskyblue1"))+
  theme_classic()+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  labs(x=expression("Fold on "~mu*"_RS"),y='Total cell number', title=expression("Reversibility "~mu*"_RS"))
ggsave(p3, file='result/muRS_on3off10_con.tiff', width=5, height=2)


total.on3off5 <- c()

for(kk in 1:length(fold)){
  
  tmp.muRS <- fold[kk]*0.0014
  parms <- c(g = 0.4, mu_SR = 0.347, mu_RS = tmp.muRS, ds = 5, dr = 4.26, k=2.67, w=0.028)
  N <- 5000
  f <- 0.239
  
  #simulate 3 days on, 5 days off
  # Set up empty data frame to store results
  out <- data.frame()
  # 10 cycles of treatment
  for (i in 1:15){
    # Set up initial conditions
    if (i == 1) {
      xstart <- c(S = f*N, R = (1-f)*N)
    } else {
      xstart <- c(S = out$S[length(out$S)], R = out$R[length(out$R)])
    }
    
    # Simulate first 3 days with trt.model
    times1 <- seq(from = 0, to = 3, by = 1)
    out1 <- ode(func = trt.model, y = xstart, times = times1, parms = parms) %>% as.data.frame()
    # Get initial values for rest.model from last day of trt.model simulation
    xstart2 <-  c(S = out1$S[length(out1$S)], R = out1$R[length(out1$R)])
    # Simulate last 5 days with rest.model
    times2 <- seq(from = 3, to = 8, by = 1)
    out2 <- ode(func = rest.model, y = xstart2, times = times2, parms = parms) %>% as.data.frame()
    
    # Combine results from both simulations and add to output data frame
    out_cycle <- rbind(out1, out2)
    out_cycle$cycle <- i
    out <- rbind(out, out_cycle)
    
    
  }
  
  results2 <- out %>%
    mutate(time2 = 8*(cycle-1)+time, N= S+R)
  total.on3off5[kk] <- results2[results2$time2==100,]$N
}

on3off5 <- data.frame(
  fold_w =  fold,
  total = total.on3off5
)



p4 <- ggplot(on3off3)+
  geom_line(aes(x=fold_w, y=total, col='ON 3 days + OFF 3 days'), size=1.5)+
  geom_line(data=on3off5, aes(x=fold_w, y=total, col='ON 3 days + OFF 5 days'), size=1.5)+
  scale_colour_manual("Dosing Strategy", breaks = c("ON 3 days + OFF 3 days", 
                                                    "ON 3 days + OFF 5 days"),
                      values = c("ON 3 days + OFF 3 days" = "dodgerblue4", 
                                 "ON 3 days + OFF 5 days" = "dodgerblue2"))+
  theme_classic()+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  labs(x=expression("Fold on "~mu*"_RS"),y='Total cell number', title=expression("Reversibility "~mu*"_RS"))
ggsave(p4, file='result/muRS_on3off5_on3off3.tiff', width=5, height=2)

p5 <- ggplot(con)+
  geom_line(aes(x=fold_w, y=total, col='Continuous'), size=1.5)+
  geom_line(data=on3off5, aes(x=fold_w, y=total, col='ON 3 days + OFF 5 days'), size=1.5)+
  scale_colour_manual("Dosing Strategy", breaks = c("Continuous", 
                                                    "ON 3 days + OFF 5 days"),
                      values = c("Continuous" = "red", 
                                 "ON 3 days + OFF 5 days" = "dodgerblue2"))+
  theme_classic()+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  labs(x=expression("Fold on "~mu*"_RS"),y='Total cell number', title=expression("Reversibility "~mu*"_RS"))
ggsave(p5, file='result/muRS_on3off5_con.tiff', width=5, height=2)



p6 <- ggplot(on3off5)+
  geom_line(aes(x=fold_w, y=total, col='ON 3 days + OFF 5 days'), size=1.5)+
  geom_line(data=on3off10, aes(x=fold_w, y=total, col='ON 3 days + OFF 10 days'), size=1.5)+
  scale_colour_manual("Dosing Strategy", breaks = c("ON 3 days + OFF 5 days", 
                                                    "ON 3 days + OFF 10 days"),
                      values = c("ON 3 days + OFF 5 days" = "dodgerblue2", 
                                 "ON 3 days + OFF 10 days" = "deepskyblue1"))+
  theme_classic()+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  labs(x=expression("Fold on "~mu*"_RS"),y='Total cell number', title=expression("Reversibility "~mu*"_RS"))
ggsave(p6, file='result/muRS_on3off5_on3off10.tiff', width=5, height=2)
