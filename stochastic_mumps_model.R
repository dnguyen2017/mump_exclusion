# compile the stochastic simulation algorithm (is called inside MumpSim function)
library(Rcpp)
sourceCpp(file = "ssa.cpp") 

# for simulating an outbreak. Returns a list containing a matrix of all times and states.
# can run N simulations by setting argument iter = N.
MumpSim <- function(Init = c(0,   # 1. time
                             100, # 2. susc
                             0,   # 3. exp
                             1,   # 4. asy
                             0,   # 5. symp
                             0,   # 6. rec
                             0,   # 7. s exc
                             0,   # 8. e exc
                             0,   # 9. a exc
                             0,   # 10. y exc
                             0,   # 11. r exc
                             50,  # 12. s vax
                             0,   # 13. e vax
                             0,   # 14. a vax
                             0,   # 15. y vax
                             0,   # 16. r vax
                             0),  # 17. y vax exc, this is excluded under the assumption that symptomatic individuals will be isolated
                    param = c((7/6),        # transmission
                              (1/6),       # recovery
                              (1/3),       # prop. asymp
                              (1/13),      # latent
                              0,           # exclusion rate
                              0,           # return rate      
                              0.05,        # vax failure rate  
                              0),          # scales contacts of excluded students
                    twait = 1.0, eff = 0.9, iter = 1, Tfinal = 3*365){
  runs <- vector("list", iter ) # init list to contain each simulations
  for (k in seq_along(runs)){
    old <- matrix(ncol=length(Init),nrow=1)
    
    old[1,] <- Init
    #print(Init)
    i <- 1; flag <- 0; Pop <- sum(Init); tex = Tfinal
    
    # set seed so that I can compare different simulations
    # set.seed(k)
    #print(k)
    #in loop
    while ((old[i,1] < Tfinal) && 
           ((old[i,2]+old[i,6]+old[i,7]+old[i,11]+old[i,12]+old[i,16]) < Pop)) { #ensures RateTotal is never 0
      
      # check if outbreak is detected (symptomatic > 1)
      if ((old[i,10]+old[i,17]) == 1 & flag == 0){
        flag <- 1
        tex <- (old[i,1] + twait)
      }
      
      if (old[i,1] >= tex & flag == 1){
        flag <- 2
        old[i,2] <- old[i-1,2] - floor(eff*old[i-1,2]) #flow out of susc. state
        old[i,3] <- old[i-1,3] - floor(eff*old[i-1,3]) #flow out of expd. state
        old[i,4] <- old[i-1,4] - floor(eff*old[i-1,4]) #flow out of asym. state
        old[i,5] <- old[i-1,5] - floor(eff*old[i-1,5]) #flow out of symp. state
        old[i,6] <- old[i-1,6] - floor(eff*old[i-1,6]) #flow out of recv. state
        old[i,7] <- old[i-1,7] + floor(eff*old[i-1,2]) #flow into s.exc state
        old[i,8] <- old[i-1,8] + floor(eff*old[i-1,3]) #flow into e.exc state
        old[i,9] <- old[i-1,9] + floor(eff*old[i-1,4]) #flow into a.exc state
        old[i,10]<- old[i-1,10] + floor(eff*old[i-1,5]) #flow into y.exc state
        old[i,11]<- old[i-1,11] + floor(eff*old[i-1,6]) #flow into r.exc state
      }  
      
      if (flag!=2)
      {
        new <- ssa(old[i,], param)
        old <- rbind(old,new)
        i <- i + 1
        #print(i);print(old[i-1,])
      }
      
      if (flag==2) 
      {
        flag <- 3
        next
      }
      
    } # while  
    
    row.names(old) <- NULL # remove row index numbers
    runs[[k]] <- old #data.frame(old, k)# store each simulation in list named "runs"
  } # for k in Iter
  
  # collect each run into a single df
  runs <- do.call(rbind, lapply(seq_along(runs), function(x) data.frame(runs[[x]], x)))
  # name cols
  names(runs) <- c("time",
                 "s","e","a","y","r",
                 "s_e", "e_e", "a_e", "y_e", "r_e",
                 "s_v", "e_v", "a_v",
                 "y_v", # should always be zero, since structural assumption in model that they are always immediately isolated
                 "r_v",
                 "y_v_e", # vaccinated students that were immediatly isolated b/c they were symptomatic
                 "iteration")
  # add parameter values to meta-data of output
  # attr(runs, "parameters") <- param # not very useful
  return(runs)
} # end function

# calculate expected proportion of unprotected individuals as function of age since last MMR dose
# these age specific probabilities can be used to estimate the expected proportion of vaccinated students in each age class
# that are susceptible to mumps

mumps_protection <- function(years_post_vax = seq(1, 20, by = 1), 
                             waning_rate = (1/27.4), 
                             init_protected = 0.964, 
                             age_last_dose = 4) {
  # parameters based on annual waning rate of protection after last dose using an exponential model (Lewnard and Grad, 2018)
  # dose scheduling https://www.cdc.gov/vaccines/vpd/mmr/public/index.html accessed on 3/14/2020
  
  prob_waned <- pexp(years_post_vax, rate = waning_rate)
  return(init_protected * (1 - prob_waned))
}