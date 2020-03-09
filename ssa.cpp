#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ssa(NumericVector old, NumericVector param) {
              
              // gillespie's direct algorithm (pg. 201 in Keeling and Rohani)
              
              // step 1. label all possible events - each rate[] is an event
              
              NumericVector rate(18);                      // initialize rate vec
              NumericVector store(4);                      // init vec to store random nums, time, and P_
              NumericVector cumsum(18);                    //= {R::cumsum(rate)};
              
              double asymp = param[2]*param[3];             // calc some rates 
              double symp  = (1-param[2])*param[3];
              
              // step 2. determine rate each event occurs
              
              rate[0]  = param[0]*old[1]*(old[3]+old[13]); //exposure
              rate[1]  = asymp*old[2];                     //asymptomatic
              rate[2]  = symp*old[2];                      //symptomatic isolation  
              rate[3]  = param[1]*old[3];                  //asymp rec
              rate[4]  = param[1]*old[4];                  //symp rec (this never happens b/c symp individuals are immediately excluded)
              rate[5]  = param[4]*old[1];                  //s exc
              rate[6]  = param[5]*old[6];                  //s reentry
              rate[7]  = asymp*old[7];                     //asymp exc
              rate[8]  = symp*old[7];                      //sympt exc
              rate[9]  = param[1]*old[8];                  //asymp rec  
              rate[10] = param[1]*old[9];                  //symp recov
              rate[11] = param[5]*old[7];                  //e reentry
              rate[12] = param[5]*old[8];                  //a reentry
              rate[13] = param[6]*param[0]*old[11]*(old[3]+old[13]);//v exposure
              rate[14] = asymp*old[12];                    //v asymptomatic
              rate[15] = symp*old[12];                     //v symptomatic isolation
              rate[16] = param[1]*old[13];                 //v asymp rec
              rate[17] = param[1]*old[16];                 //v symp rec
              
              // step 3. the rate any event occurs is r_total = cumsum[n-1] 
              // note, n-1 is r_total because c++ index starts at 0
              
              int n = rate.size();
              double sum = 0;
              for(int i = 0; i < n; i++){  // calculate cumsum
                sum += rate[i];
                cumsum[i] = sum;} // end for
              
              
              store[0] = R::runif(0,1); // rand1
              store[1] = R::runif(0,1); // rand2
              
              // step 4. the time until next event occurs
              
              store[2] = (-1.0/cumsum[n-1])*log(store[0]); // size of time step
              
              // step 5. set p = rand2 * R_total
              
              store[3] = store[1]*cumsum[n-1];             // p = Rand2 * R_total
              
              // step 6. event p occurs if the cumsum[p - 1] < p <= cumsum[p]
              // use if else structure to sequentially test all values of p
 
              if(store[3] <= cumsum[0])       {old[0] += store[2]; --old[1]; ++old[2];  // s -> e
              } else if(store[3] <= cumsum[1]){old[0] += store[2]; --old[2]; ++old[3];  // e -> a
              } else if(store[3] <= cumsum[2]){old[0] += store[2]; --old[2]; ++old[9];  // e -> y_excluded
              } else if(store[3] <= cumsum[3]){old[0] += store[2]; --old[3]; ++old[5];  // a -> r
              } else if(store[3] <= cumsum[4]){old[0] += store[2]; --old[4]; ++old[5];  // y -> r (doesn't happen if we assume immediate isolation of symptomatic individuals, since there are no y (only y_exlcuded))
              } else if(store[3] <= cumsum[5]){old[0] += store[2]; --old[1]; ++old[6];  // s -> s_excluded (doesn't happen if we assume all excluded instantaneously, since exclusion rate == 0)
              } else if(store[3] <= cumsum[6]){old[0] += store[2]; ++old[1]; --old[6];  // s_exc -> s (re-entry; disobey exclusion)
              } else if(store[3] <= cumsum[7]){old[0] += store[2]; --old[7]; ++old[8];  // e_exc -> a_exc
              } else if(store[3] <= cumsum[8]){old[0] += store[2]; --old[7]; ++old[9];  // e_exc -> y_exc
              } else if(store[3] <= cumsum[9]){old[0] += store[2]; --old[8]; ++old[10]; // a_exc -> r_exc
              } else if(store[3] <= cumsum[10]){old[0] += store[2]; --old[9]; ++old[10];// y_exc -> r_exc 
              } else if(store[3] <= cumsum[11]){old[0] += store[2]; ++old[2]; --old[7]; // e_exc -> e (re-entry; disobey exclusion)
              } else if(store[3] <= cumsum[12]){old[0] += store[2]; ++old[3]; --old[8]; // y_exc -> y (re-entry; discobey exclusion)
              } else if(store[3] <= cumsum[13]){old[0] += store[2]; --old[11]; ++old[12]; // s_vax -> e_vax
              } else if(store[3] <= cumsum[14]){old[0] += store[2]; --old[12]; ++old[13]; // e_vax -> a_vax
              } else if(store[3] <= cumsum[15]){old[0] += store[2]; --old[12]; ++old[16]; // e_vax -> y_vax_exc
              } else if(store[3] <= cumsum[16]){old[0] += store[2]; --old[13]; ++old[15]; // a_vax -> r_vax
              } else if(store[3] <= cumsum[17]){old[0] += store[2]; --old[16]; ++old[15]; // y_vax_exc -> r_vax (assume ppl return to school after recovery)
              }  
              
              return old;
              } //end of ssa
  