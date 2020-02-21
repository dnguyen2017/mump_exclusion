#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ssa(NumericVector old, NumericVector param) {
              
              NumericVector rate(18);                      // initialize rate vec
              NumericVector store(4);                      // init vec to store random nums, time, and P_
              NumericVector cumsum(18);                    //= {R::cumsum(rate)};
              
              double asymp = param[2]*param[3];             // calc some rates 
              double symp  = (1-param[2])*param[3];
              
              rate[0]  = param[0]*old[1]*(old[3]+old[13]); //exposure
              rate[1]  = asymp*old[2];                     //asymptomatic
              rate[2]  = symp*old[2];                      //symptomatic isolation  
              rate[3]  = param[1]*old[3];                  //asymp rec
              rate[4]  = param[1]*old[4];                  //symp rec
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
              
              int n = rate.size();
              double sum = 0;
              for(int i = 0; i < n; i++){
                sum += rate[i];
                cumsum[i] = sum;} // end for
  
              store[0] = R::runif(0,1);
              store[1] = R::runif(0,1);
              store[2] = (-1.0/cumsum[n-1])*log(store[0]);
              store[3] = store[1]*cumsum[n-1];
 
              if(store[3] <= cumsum[0])       {old[0] += store[2]; --old[1]; ++old[2];
              } else if(store[3] <= cumsum[1]){old[0] += store[2]; --old[2]; ++old[3];
              } else if(store[3] <= cumsum[2]){old[0] += store[2]; --old[2]; ++old[9];
              } else if(store[3] <= cumsum[3]){old[0] += store[2]; --old[3]; ++old[5];
              } else if(store[3] <= cumsum[4]){old[0] += store[2]; --old[4]; ++old[5];
              } else if(store[3] <= cumsum[5]){old[0] += store[2]; --old[1]; ++old[6];
              } else if(store[3] <= cumsum[6]){old[0] += store[2]; ++old[1]; --old[6];
              } else if(store[3] <= cumsum[7]){old[0] += store[2]; --old[7]; ++old[8];
              } else if(store[3] <= cumsum[8]){old[0] += store[2]; --old[7]; ++old[9];
              } else if(store[3] <= cumsum[9]){old[0] += store[2]; --old[8]; ++old[10];
              } else if(store[3] <= cumsum[10]){old[0] += store[2]; --old[9]; ++old[10];
              } else if(store[3] <= cumsum[11]){old[0] += store[2]; ++old[2]; --old[7];
              } else if(store[3] <= cumsum[12]){old[0] += store[2]; ++old[3]; --old[8];
              } else if(store[3] <= cumsum[13]){old[0] += store[2]; --old[11]; ++old[12];
              } else if(store[3] <= cumsum[14]){old[0] += store[2]; --old[12]; ++old[13];
              } else if(store[3] <= cumsum[15]){old[0] += store[2]; --old[12]; ++old[16];
              } else if(store[3] <= cumsum[16]){old[0] += store[2]; --old[13]; ++old[15];
              } else if(store[3] <= cumsum[17]){old[0] += store[2]; --old[16]; ++old[15];
              }  
              
              return old;
              } //end of ssa
  