# deterministic mumps model
# impulsive exclusion

library(deSolve)
library(tidyverse)
library(viridis)

mumps_model <- function(time, y, parms) {
  with(as.list(c(y, parms)), {
    
    contacts <- beta*(A + Y + A_v + Y_v + c_*(A_e + Y_e)) # transmission * contacts
    
    
    # unvaccinated in school
    dS <- - S*contacts
    dE <- S*contacts - mu*E
    dA <- epsilon*mu*E - gam*A
    dY <- (1-epsilon)*mu*E - gam*Y
    dR <- gam*(A + Y)
    
    # vaccinated in school
    dS_v <- - v*S_v*contacts # scaled by vaccine failure rate (v)
    dE_v <- v*S_v*contacts - mu*E_v
    dA_v <- epsilon*mu*E_v - gam*A_v
    dY_v <- (1-epsilon)*mu*E_v - gam*Y_v
    dR_v <- gam*(A_v + Y_v)
    
    dS_e <- - c_*S_e*contacts # scaled by exclusion induced contact reduction (c)
    dE_e <- c_*S_e*contacts - mu*E_e
    dA_e <- epsilon*mu*E_e - gam*A_e
    dY_e <- (1-epsilon)*mu*E_e - gam*Y_e
    dR_e <- gam*(A_e + Y_e)
    
    flag <- 0 # to control when exclusion happens
    t_detected <- 0 # time of outbreak detection
    
    return(list(c(dS, dE, dA, dY, dR,
                  dS_v, dE_v, dA_v, dY_v, dR_v,
                  dS_e, dE_e, dA_e, dY_e, dR_e,
                  flag, t_detected)))
  }) # with
} # function


exclusion <- function(time, y, parms) {
  out <- with(as.list(c(y, parms)),{
    
    if (((Y + Y_v) >= detection_threshold) & (flag == 0)) {
      t_detected <- time # exclusion will happen at time == t_detected + t_delay
    # print(paste("detected at t = ", t_detected))
      flag <- 1
    }
    
    if ((time >= (t_detected + t_delay)) & (flag == 1)){
      
      # report that exclusion happened
    #  print(paste("exclusion at t = ", time))
      
      # calculated number to move in each compartment
      S_move <- efficacy*S
      E_move <- efficacy*E
      A_move <- efficacy*A
      Y_move <- efficacy*Y
      R_move <- efficacy*R
      
      # remove excluded students from school
      S <- S - S_move
      E <- E - E_move
      A <- A - A_move
      Y <- Y - Y_move
      R <- R - R_move
      
      # add excluded students to exclusion category
      S_e <- S_e + S_move
      E_e <- E_e + E_move
      A_e <- A_e + A_move
      Y_e <- Y_e + Y_move
      R_e <- R_e + R_move
      
      # prevent exclusion from happening more than once
      flag <- 2
    }
    
    new_states <- c(S, E, A, Y, R,
                    S_v, E_v, A_v, Y_v, R_v,
                    S_e, E_e, A_e, Y_e, R_e,
                    flag, t_detected)
    
    names(new_states) <- c("S", "E", "A", "Y", 'R',
                           "S_v", "E_v", "A_v", "Y_v", "R_v",
                           "S_e", "E_e", "A_e", "Y_e", "R_e",
                           "flag", "t_detected")
    new_states 
  })
  return(out)
}

# init
y0 <- c(S = 99, E = 0, A = 1, Y = 0, R = 0, 
        S_v = 400, E_v = 0, A_v = 0, Y_v = 0, R_v = 0,
        S_e = 0, E_e = 0, A_e = 0, Y_e = 0, R_e = 0,
        flag = 0, t_detected = 0)

tmax <- 1*30 # t and event_times must be same length
t <- seq(0, tmax, 0.01)
event_times <- seq(0, tmax, 1)

pars <- c(beta =  0.001683,  # 7/6,       # transmission rate
          gam = 1/6,        # recovery: 1/infection duration
          epsilon = 1/3,    # proportion asymptomatic
          mu = 1/13,        # becoming infectious: 1/latent period
          v = 0.0, c_ = 0,  # vax failure rate, contact rate of excluded students 
          detection_threshold = 0.1, efficacy = 0.9, t_delay = 0) # control pars

print(paste("R0 = ",unname(pars["beta"])*unname(y0["S"])/unname(pars["gam"]),""))

out <- ode(y = y0, times = t, parms = pars, func = mumps_model,
           method = "rk4")


out_event <- ode(y = y0, times = t, parms = pars, func = mumps_model,
                 events = list(func = exclusion, times = event_times), 
                 method = "rk4")

out_df <- 
  out %>% 
  as.data.frame() %>% 
  gather(compartment, number, S:R_e) %>%
  select(-flag,-t_detected) %>%
  separate(compartment, into = c("class", "status"), sep = "\\_") %>%
  mutate(status = case_when(is.na(status) == TRUE ~ "unvaccinated",
                            status =="v" ~ "vaccinated",
                            status == "e" ~ "excluded"))

out_event_df <- 
  out_event %>% 
  as.data.frame() %>% 
  gather(compartment, number, S:R_e) %>%
  select(-flag,-t_detected) %>%
  separate(compartment, into = c("class", "status"), sep = "\\_") %>%
  mutate(status = case_when(is.na(status) == TRUE ~ "unvaccinated",
                            status =="v" ~ "vaccinated",
                            status == "e" ~ "excluded"))
  
# plot without events
ggplot(out_df, aes(x = time, y = number, color = class)) + 
  geom_line() +
  facet_wrap(~status, nrow = 3)

# plot with exclusion
ggplot(out_event_df, aes(x = time, y = number, color = class)) + 
  geom_line() +
  facet_wrap(~status, nrow = 3)

# summarizing
out_df %>% group_by(status) %>%
  filter(class %in% c("S","R") & time == tmax) %>% mutate(number = round(number, digits = 1))

# summarizing
out_event_df %>% group_by(status) %>%
  filter(class %in% c("S","R") & time == tmax) %>% mutate(number = round(number, digits = 1))


# Vaccine protection is perfect
# dS_V is the fraction of vaccinated individuals that did not develop immune response (primary vaccine failure)
# i.e., S_v = (1-(fail rate))(total vaccinated) 

mumps_model_2 <- function(time, y, parms) {
  with(as.list(c(y, parms)), {
    
    contacts <- beta*(A + Y + A_v + Y_v + c_*(A_e + Y_e)) # transmission * contacts
    
    
    # unvaccinated in school
    dS <- - S*contacts
    dE <- S*contacts - mu*E
    dA <- epsilon*mu*E - gam*A
    dY <- (1-epsilon)*mu*E - gam*Y
    dR <- gam*(A + Y)
    
    # vaccinated in school
    dS_v <- - S_v*contacts # scaled by vaccine failure rate (v)
    dE_v <- S_v*contacts - mu*E_v
    dA_v <- epsilon_v*mu*E_v - gam*A_v
    dY_v <- (1-epsilon_v)*mu*E_v - gam*Y_v
    dR_v <- gam*(A_v + Y_v)
    
    dS_e <- - c_*S_e*contacts # scaled by exclusion induced contact reduction (c)
    dE_e <- c_*S_e*contacts - mu*E_e
    dA_e <- epsilon*mu*E_e - gam*A_e
    dY_e <- (1-epsilon)*mu*E_e - gam*Y_e
    dR_e <- gam*(A_e + Y_e)
    
    flag <- 0 # to control when exclusion happens
    t_detected <- 0 # time of outbreak detection
    
    return(list(c(dS, dE, dA, dY, dR,
                  dS_v, dE_v, dA_v, dY_v, dR_v,
                  dS_e, dE_e, dA_e, dY_e, dR_e,
                  flag, t_detected)))
  }) # with
} # function

fail2 <- 0.05
y02 <- c(S = 99, E = 0, A = 1, Y = 0, R = 0, 
         S_v = 400*fail2, E_v = 0, A_v = 0, Y_v = 0, R_v = 0,
         S_e = 0, E_e = 0, A_e = 0, Y_e = 0, R_e = 0,
         flag = 0, t_detected = 0)
# R0 is around 4 - 7 (Paul Fine. Herd Immunity: history, theory, practice. Epid. Rev. 1993)
# so, if primary vax failure (fail2) is 0.05, (estimate for protection at 6 months from Lewnard & Grad 2018; in the future, should use diff values to account for waning)
# [4,7] * (400*fail2 + unvax)/N_total = [0.952,1.66]
# assume Rv is always > 0 for these simulations and that gamma = 1/6, 
# beta = [1/(119/(1/6)), 1.66/(119/(1/6))] = [1/714, 5/(3*714)]
# beta = [4/(119/(1/6)), 7/(119/(1/6))] = [4/714, 7/714)]

pars2 <- c(beta =  7/714,  # transmission rate
          gam = 1/6,        # recovery: 1/infection duration
          epsilon = 1/3,    # proportion asymptomatic
          epsilon_v = 0.672, # vax asympmtomatic rate
          mu = 1/13,        # becoming infectious: 1/latent period
          v = NA, c_ = 0,  # vax failure rate, contact rate of excluded students 
          detection_threshold = 1, efficacy = 1, t_delay = 21) # control pars

tmax <- 6*30 # t and event_times must be same length
t <- seq(0, tmax, 0.01)
event_times2 <- seq(0, tmax, 0.1)

out_2 <-  ode(y = y02, times = t, parms = pars2, func = mumps_model_2,
              events = list(func = exclusion, times = event_times2), 
              method = "rk4")
out_2_df <- 
  out_2 %>% 
  as.data.frame() %>% 
  gather(compartment, number, S:R_e) %>%
  select(-flag,-t_detected) %>%
  separate(compartment, into = c("class", "status"), sep = "\\_") %>%
  mutate(status = case_when(is.na(status) == TRUE ~ "unvaccinated",
                            status =="v" ~ "vaccinated",
                            status == "e" ~ "excluded"))
# look at epi curve
out_2 %>% 
  as.data.frame() %>% 
  filter(time %in% seq(0, 60, by = 1)) %>%
  mutate(cases_u = lead(Y) - Y,
         cases_v = lead(Y_v) - Y_v) %>% 
  select(time, Y, Y_v, cases_u, cases_v) -> epicurve_df

plot(cases_u ~ time, data = epicurve_df)

# plot with exclusion
ggplot(out_2_df, aes(x = time, y = number, color = class)) + 
  geom_line() +
  facet_wrap(~status, nrow = 3)

ggplot(filter(out_2_df, status %in% c("unvaccinated","excluded") ), aes(x = time, y = number, color = class)) + 
  geom_point()


# summarizing
out_2_df %>% group_by(status) %>%
  filter(class %in% c("S","R") & time == tmax) %>% mutate(number = round(number, digits = 3)) -> summary_2
summary_2
sum(summary_2$number)

# Simulate across parameter space of control pars
effRange <- seq(0,1, by = 0.1) # range of exclusion efficacy
delRange <- seq(0,28,3.5) # range of delay
cRange <- seq(0,1, by = 0.1) # proportional reduction in contacts when excluded

controlPars <- expand.grid(effRange,delRange,cRange)
names(controlPars) <- c("efficacy", "delay", "contact reduction")
store_wide <- vector("list", length = nrow(controlPars))
store_long <- store_wide

ptm <- proc.time()
for(i in 1:nrow(controlPars)) { # 
  
  # pars with control pars from row i
  pars2 <- c(beta =  4/714,  # 7/6,       # transmission rate
            gam = 1/6,        # recovery: 1/infection duration
            epsilon = 0.3,    # unvax proportion asymptomatic
            epsilon_v = 0.672, # vax asympmtomatic rate
            mu = 1/13,        # becoming infectious: 1/latent period
            v = NA, c_ = controlPars[i,3],  # vax failure rate, contact rate of excluded students 
            detection_threshold = 1, efficacy = controlPars[i,1], t_delay = controlPars[i,2]) # control pars
  # simulate for controlPars[i,]
  out <-  ode(y = y02, times = t, parms = pars2, func = mumps_model_2,
                events = list(func = exclusion, times = event_times2), 
                method = "rk4")
  # reformat
  # out_2_df <- 
  #   out %>% 
  #   as.data.frame() %>% 
  #   gather(compartment, number, S:R_e) %>%
  #   select(-flag,-t_detected) %>%
  #   separate(compartment, into = c("class", "status"), sep = "\\_") %>%
  #   mutate(status = case_when(is.na(status) == TRUE ~ "unvaccinated",
  #                             status =="v" ~ "vaccinated",
  #                             status == "e" ~ "excluded"))
  # store df
  store_wide[[i]] <- out
  # store_long[[i]] <- out_2_df
  
  # print
  if(i %% 50 == 0) print(i)
}
runtime1 <- proc.time() - ptm
runtime1

# moke df to store the output for each parcombo
result_df <- cbind(controlPars, matrix(0, nrow = nrow(controlPars), ncol = 20))
cNames <- names(round(store_wide[[1]][length(t),], digits = 3)[c(2:16,18)]) # get vec of names.
names(result_df)[4:19] <- cNames
names(result_df)[20:23] <- c("tPeakExp","hPeakExp","tPeakInf","hPeakInf")
names(result_df)

# fill result_df 
for (i in 1:length(store_wide)){ 
  # pop status at t_final, rounded to thousandths
  result_df[i,4:19] <- round(store_wide[[i]][length(t),], digits = 3)[c(2:16,18)]
  
  # peak of exposed
  ePeakPos <- which.max(store_wide[[i]][,3]+store_wide[[i]][,8]+store_wide[[i]][,13]) # index of peak
  result_df[i,20] <- store_wide[[i]][ePeakPos,1] # tPeakExp
  result_df[i,21] <- max(store_wide[[i]][,3]+store_wide[[i]][,8]+store_wide[[i]][,13]) # heigth of peak exp
  
  # peak of infected (A + Y)
  ePeakPos <- which.max(store_wide[[i]][,4]+store_wide[[i]][,9]+store_wide[[i]][,14]+
                        store_wide[[i]][,5]+store_wide[[i]][,10]+store_wide[[i]][,15]) # index of peak
  result_df[i,22] <- store_wide[[i]][ePeakPos,1] # time of peak inf
  result_df[i,23] <- max(store_wide[[i]][,4]+store_wide[[i]][,9]+store_wide[[i]][,14]+ # height of peak inf
                         store_wide[[i]][,5]+store_wide[[i]][,10]+store_wide[[i]][,15]) 
}

saveRDS(result_df, file = "data/deterministic_heat_low_r0.RDS") # need to increase t_final, epidemic is not over w/i 90 days

# high r0, beta = 7/714


ptm <- proc.time()
for(i in 1:nrow(controlPars)) { # 
  
  # pars with control pars from row i
  pars2 <- c(beta =  7/714,  # 7/6,       # transmission rate
             gam = 1/6,        # recovery: 1/infection duration
             epsilon = 0.3,    # unvax proportion asymptomatic
             epsilon_v = 0.672, # vax asympmtomatic rate
             mu = 1/13,        # becoming infectious: 1/latent period
             v = NA, c_ = controlPars[i,3],  # vax failure rate, contact rate of excluded students 
             detection_threshold = 1, efficacy = controlPars[i,1], t_delay = controlPars[i,2]) # control pars
  # simulate for controlPars[i,]
  out <-  ode(y = y02, times = t, parms = pars2, func = mumps_model_2,
              events = list(func = exclusion, times = event_times2), 
              method = "rk4")
  # reformat
  # out_2_df <- 
  #   out %>% 
  #   as.data.frame() %>% 
  #   gather(compartment, number, S:R_e) %>%
  #   select(-flag,-t_detected) %>%
  #   separate(compartment, into = c("class", "status"), sep = "\\_") %>%
  #   mutate(status = case_when(is.na(status) == TRUE ~ "unvaccinated",
  #                             status =="v" ~ "vaccinated",
  #                             status == "e" ~ "excluded"))
  # store df
  store_wide[[i]] <- out
  # store_long[[i]] <- out_2_df
  
  # print
  if(i %% 50 == 0) print(i)
}
runtime2 <- proc.time() - ptm
runtime2

# moke df to store the output for each parcombo
result_df <- cbind(controlPars, matrix(0, nrow = nrow(controlPars), ncol = 20))
cNames <- names(round(store_wide[[1]][length(t),], digits = 3)[c(2:16,18)]) # get vec of names.
names(result_df)[4:19] <- cNames
names(result_df)[20:23] <- c("tPeakExp","hPeakExp","tPeakInf","hPeakInf")
names(result_df)

# fill result_df 
for (i in 1:length(store_wide)){ 
  # pop status at t_final, rounded to thousandths
  result_df[i,4:19] <- round(store_wide[[i]][length(t),], digits = 3)[c(2:16,18)]
  
  # peak of exposed
  ePeakPos <- which.max(store_wide[[i]][,3]+store_wide[[i]][,8]+store_wide[[i]][,13]) # index of peak
  result_df[i,20] <- store_wide[[i]][ePeakPos,1] # tPeakExp
  result_df[i,21] <- max(store_wide[[i]][,3]+store_wide[[i]][,8]+store_wide[[i]][,13]) # heigth of peak exp
  
  # peak of infected (A + Y)
  ePeakPos <- which.max(store_wide[[i]][,4]+store_wide[[i]][,9]+store_wide[[i]][,14]+
                          store_wide[[i]][,5]+store_wide[[i]][,10]+store_wide[[i]][,15]) # index of peak
  result_df[i,22] <- store_wide[[i]][ePeakPos,1] # time of peak inf
  result_df[i,23] <- max(store_wide[[i]][,4]+store_wide[[i]][,9]+store_wide[[i]][,14]+ # height of peak inf
                           store_wide[[i]][,5]+store_wide[[i]][,10]+store_wide[[i]][,15]) 
}

saveRDS(result_df, file = "data/deterministic_heat_high_r0.RDS") # also need to increase t_final

