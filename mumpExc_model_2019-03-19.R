library(Rcpp)

# set number of samples to take (iter)
iter <- 1000000
# set how many samples should have return rate = 0
enforced <- round(0.05*iter)

### not totally sure if runif() is the correct sampling for each parm. I think maybe taking
### a sample(seq(min,max,by), iter) may be better for some things, e.g., twait, efficacy

# create param combos

set.seed(123)

twait <- sample(1:7, iter, replace = TRUE) # 1:7 integers
eff_vec <- rep_len(0, length.out = iter) #runif(iter, 0, 1)  # make 1 run all 0; another run at 0.9
ret_vec <- sample(c(rep(0,times=enforced), 
                    runif(iter-enforced,10^(-3.2),10^(-0.8))), size=iter)
beta_vec <- runif(iter,1.2, 2.3) # R0[7,14] and gamma = 1/6 -> beta[1.66,2.33]
vax_fail <- runif(iter,0.05,0.69)
asymp_vec <- runif(iter,0.15,0.3)
exempt_vec <- runif(iter,0.15,0.5) # this range ensures that r0 > 1 for all combos

# concatenate all parm vectors
par_vec <- c(twait,eff_vec,ret_vec,beta_vec,vax_fail,asymp_vec,exempt_vec)

# make par_vec into a matrix. Each col is one of the parameters
parms <- matrix(par_vec,
                nrow = iter, ncol=(length(par_vec)/iter))
rm(twait,eff_vec,ret_vec,beta_vec,vax_fail,asymp_vec,exempt_vec)

# compile the stochastic simulation algorithm (is called inside MumpSim function)
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
                             0),  # 17. y vax exc
                    param = c((7/6),        # transmission
                              (1/6),       # recovery
                              (1/3),       # prop. asymp
                              (1/13),      # latent
                              0,           # exclusion rate
                              0,           # return rate      
                              0.05),       # vax failure rate  
                    twait = 1.0, eff = 0.9, iter = 1, Tfinal = 3*365){
  runs <- vector("list", iter ) # init list to contain each simulations
  for (k in 1:iter){
    old <- matrix(ncol=length(Init),nrow=1)
    
    old[1,] <- Init
    #print(Init)
    i <- 1; flag <- 0; Pop <- sum(Init); tex = Tfinal
    
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
    runs[[k]] <- old # store each simulation in list named "runs"
  } # for k in Iter
  
  return(runs)
} # end function

# demo of MumpSim function
# out <- MumpSim(Init = c(0, 100, 0, 1, 0, 0, 0, 0, 0, 0, 0, 400, 0, 0, 0, 0, 0))
# df <- as.data.frame(out[[1]])

# Takes in matrix of parameter combinations. Output is a matrix of the param combinations
# with the total duration and mumps cases. The unvax and vax cases are tracked seperately.
# this function also automatically saves the output according to the argument "fileName"
sensitivityMumps <- function(nstudent,parameters, reps, 
                             # use saveRDS instead of write.csv
                             
                             fileName = paste0("data/",format(Sys.time(), "%Y%m%d_%H%M_"), "epidata.rds"),
                             Tfinal = 3*365,
                             saveFreq = 100) 
                               #paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"), "epidata.csv"))
  {
    # initialize data matrices and output list 
    width1 <- ncol(parameters)+1
    width2 <- ncol(parameters)+2
    width3 <- ncol(parameters)+3
    
    epidata_mat     <- matrix(data = NA, nrow = nrow(parameters)*reps, ncol = width3)
    epidata_mat[,1] <- rep(parameters[,1], each = reps)
    epidata_mat[,2] <- rep(parameters[,2], each = reps)
    epidata_mat[,3] <- rep(parameters[,3], each = reps)
    epidata_mat[,4] <- rep(parameters[,4], each = reps)
    epidata_mat[,5] <- rep(parameters[,5], each = reps)
    epidata_mat[,6] <- rep(parameters[,6], each = reps)
    epidata_mat[,7] <- rep(parameters[,7], each = reps)
    
    #output <- vector("list", length = 2) # to return epitime_mat and episize_mat from sensitivityExp
    
    i <- 1 # counter for each param combo
    simCount <- 0 # count for simulation number; final should = reps*nrow(parameters)
    
    for (i in 1:nrow(parameters)){
      
      # epidemic simulations stored in list 'result'
      result <- MumpSim(Init = c(0, round(nstudent*parameters[i,7]), 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 (nstudent-round(nstudent*parameters[i,7])), 0, 0, 0, 0,0),
                        param = c(parameters[i,4], # beta
                                  (1/6.0),         # rec
                                  parameters[i,6], # asymp  
                                  (1/13.0),        # latent
                                  0,              # exclusion rate
                                  parameters[i, 3], # return rate
                                  parameters[i,5]), # vax eff rate
                        twait = parameters[i, 1], 
                        eff = parameters[i, 2], 
                        iter = reps,
                        Tfinal = Tfinal)
      countResult <- 1 # counter for extracting correct element in result
      
      for (simCount in (simCount+1):(simCount+reps)){
        # store total duration
        epidata_mat[simCount, width1] <- tail(result[[countResult]], 1)[1]     # time
        
        # store total cases
        epidata_mat[simCount, width2] <- tail(result[[countResult]], 1)[6] +  # rec
          tail(result[[countResult]], 1)[11]    # r.exc
        
        epidata_mat[simCount, width3] <- tail(result[[countResult]], 1)[16]  # v.rec
        
        countResult <- countResult + 1
      } # FOR inputting data into epitime_mat and episize_mat
      
      rm(result)   
      gc(verbose = FALSE)   # clear memory after rm() & dont print memory stats
      
      # periodically save results
      if((i %% saveFreq)==0){
        print(i)
        # write.csv(epidata_mat, file = fileName, row.names = FALSE)
        saveRDS(epidata_mat, file = fileName)
        print("saved")
        
      } # IF
    } # FOR generating data with MumpSim
    
    # if not already saved, save final datasets  
    if ( ((nrow(parameters)) %% saveFreq) != 0){
      #write.csv(epidata_mat, file = fileName, row.names = FALSE)
      saveRDS(epidata_mat, file = fileName)
    } # if  
    
    return(epidata_mat)
  } # function

ptm <- proc.time()
out1 <- sensitivityMumps(nstudent = 500, parameters = parms, 
                         reps = 1, fileName = "data/mumps_run_2019-03-19.rds", saveFreq = 10000)
runtime <- proc.time() - ptm
runtime

# ptm <- proc.time()
# out1 <- sensitivityMumps(nstudent = 500, parameters = parms, reps = 1,
#                          fileName = "data/testrun.csv")
# runtime <- proc.time() - ptm

df1 <- as.data.frame(out1)
names(df1) <- c("delay", "efficacy", "return","beta","vaccine failure",
                "asympt","exempt.rate","duration","unvax","vax")
write.csv(df1, file="data/mumps_run_2019_-3_19.csv", row.names = FALSE)
# need to have cols for init vax and unvax pop sizes
