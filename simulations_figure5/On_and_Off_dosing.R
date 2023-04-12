# project: ADCC
# Author: Jiawei Zhou
# Objectives: simulate R and S and R+S under on/off treatment

rm(list=ls())

library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(deSolve)


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
parms <- c(g = 0.4, mu_SR = 4.3e-5, mu_RS = 0.0014, ds = 6.03, dr = 3.51, k=2.48, w=0.025)
N <- 5000
f <- 0.239
# Set up empty data frame to store results
out <- data.frame()
# 10 cycles of treatment
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
write.csv(results2, file='3on3off.csv', quote=FALSE, row.names=FALSE)
# for plot
#plot
p1 <- results2 %>%gather(variable,value,-time2, -cycle, -time, -N) %>%
  ggplot(aes(x=time2,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  scale_y_continuous(trans='log10')+
  labs(x='Time (day)',y='Cell number')+ggtitle('ON 3 days + OFF 3 days')
ggsave(p1, file='on3off3.tiff', width=4, height=2)

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
write.csv(results2, file='3on5off.csv', quote=FALSE, row.names=FALSE)
#plot
p2 <- results2 %>%gather(variable,value,-time2, -cycle, -time, -N) %>%
  ggplot(aes(x=time2,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  scale_y_continuous(trans='log10')+
 labs(x='Time (day)',y='Cell number')+ggtitle('ON 3 days + OFF 5 days')
ggsave(p2, file='on3off5.tiff', width=4, height=2)


#simulate 3 days on, 10 days off
# Set up empty data frame to store results
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
write.csv(results2, file='3on10off.csv', quote=FALSE, row.names=FALSE)
#plot
p3 <- results2 %>%gather(variable,value,-time2, -cycle, -time, -N) %>%
  ggplot(aes(x=time2,y=value,color=variable))+
  geom_line(size=1)+
  theme_classic()+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(breaks = c(0,25,50,75,100,125),labels = c(0,25,50,75,100,125),lim = c(0,125))+
  labs(x='Time (day)',y='Cell number')+ggtitle('ON 3 days + OFF 10 days')
ggsave(p3, file='on3off10.tiff', width=4, height=2)


#plot all scenarios
con <- read.csv(file='continuous.csv',header=T)
on3off3 <- read.csv(file='3on3off.csv',header=T)
on3off5 <- read.csv(file='3on5off.csv',header=T)
on3off10 <- read.csv(file='3on10off.csv',header=T)

p5 <- ggplot(con)+
  geom_line(aes(x=time2, y=Total, col='Continuous'), size=1.5)+
  geom_line(data=on3off3, aes(x=time2, y=N, col='ON 3 days + OFF 3 days'), size=1)+
  geom_line(data=on3off5, aes(x=time2, y=N, col='ON 3 days + OFF 5 days'), size=1)+
  geom_line(data=on3off10, aes(x=time2, y=N, col='ON 3 days + OFF 10 days'), size=1)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(breaks = c(0,25,50,75,100,125),labels = c(0,25,50,75,100,125),lim = c(0,125))+
  labs(x='Time (day)', y='Total cell number', title='Total = R + S')+
  scale_colour_manual("Dosing Strategy", breaks = c("Continuous", 
                                              "ON 3 days + OFF 3 days",
                                              "ON 3 days + OFF 5 days" ,
                                              "ON 3 days + OFF 10 days"),
                      values = c("Continuous" = "red", 
                                              "ON 3 days + OFF 3 days" = "dodgerblue4",
                                              "ON 3 days + OFF 5 days" = "dodgerblue2",
                                              "ON 3 days + OFF 10 days" = "deepskyblue1"))+
  theme_classic()
ggsave(p5, file='total.tiff', width=8, height=3)
