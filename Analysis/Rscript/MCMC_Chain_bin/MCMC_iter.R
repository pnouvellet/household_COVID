#' MCMC iterate
#'
#' run the MCMC to sample posterior of R and initial coniditions at each location 
#' FYI: this is called internally by adapt_tuning
#' 
#' @param incidence the incidence for the time window during which we assume Rt to be constant.  
#'           I is a dataframe, first column are dates then incidence for all locations
#'           nb of row is the size of time widows, dates must be sequential
#' 
#' @param N_geo integer of  numbers of locations
#'                   
#' @param iter integer, the number of iteration for the MCMC
#'
#' @param theta0 vector of inital parameters, here taken from the last MCMC iteration after tuning (save some burn-in)
#'
#' @param s variance of proposal distributions (log-normal) - tuned previously
#' 
#' @param SI Serial interval distribution (see SI_gamma_dist_EpiEstim)
#' 
#' @param mu0: initial conidtions to guaranty that if R=1, then we predict the number of cases in the future will stablise at the mean number of cases observed in the time window
#'              mu0 is also used as the mean of the (exponential) prior for intial conditions estimated
#' 
#' @details  res a list containing 2 matrices: theta: matrix of posterior samples
#'                      and logL: matrix of associated log-likelihood
#' @export
#' 

MCMC_iter <- function(n_size, N,iter,theta0,s,prior){
  # options(warn=2)
   #############################################################################
  # for MCMC
  thetas <- matrix(NA,iter,length(theta0))  # Rt's and initial conditions for each location
  L <- thetas                         # store likelihoods
  Ps <- array(NA,list(length(n_size)+1,length(n_size)+1,iter))
  #############################################################################
  # get the daily 'forces of infections'
 
  P <- lambda_fct(n_s = n_size, theta = theta0)
  L1 <- Like1(n = N, p = P)

  L[1,] <- L1
  thetas[1,] <- theta0    
  Ps[,,1] <- P
  
  #############################################################################
  # sampling
  for (i in 2:iter){   
    #print(i)
    for (j in 1:length(theta0)){
      Ts <- theta0
      Ts[j] <- Ts[j]*exp(s[j]*rnorm(1,0,1))
      
      if ( ((Ts[j] < prior[j,1]) | (Ts[j] > prior[j,2])) ){
        r <- 0 
      }else{
        
        Pint <- lambda_fct(n_s = n_size, theta = Ts)
        Lint <- Like1(n = N, p = Pint)
        
        #get the ratio with previous value of parameter and correct for porposal (and, only for initial conditions, prior distribution)
        r <- exp(Lint-L1)*Ts[j]/theta0[j]
      }
      
      
      # accept or reject
    
        if (runif(1,0,1) <= r){
          theta0 <- Ts  # if accept, keep new parameter value
          L1 <- Lint         # if accept, keep new lieklihood
          P <- Pint
        } 
      
      L[i,j] <- L1    # store final likelihood
    }
    
    thetas[i,] <- theta0  # store final parameter values for this iteration
    Ps[,,i] <- P
  }
  
  ####
  res <- list(theta = thetas, logL = L, Ps = Ps)
  return(res)
}
