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

MCMC_iter <- function(death,iter,theta0,s,over_disp = NA){
  
   #############################################################################
  # for MCMC
  thetas <- matrix(NA,iter,length(theta0))  # Rt's and initial conditions for each location
  L <- thetas                         # store likelihoods
  Traj <- matrix(NA,iter,nb_day)
  Rt <-  matrix(NA,iter,length(times))
  
  #############################################################################
  # get the daily 'forces of infections'
  input <- list(   beta_w = theta0[n_change+1],
                   delta = 1/4,  # rate for exposed period
                   sigma  = 1/ 4.5,  # rate for infectiousness individual
                   mu = 0.01,  # IFR
                   hs = 1/(18.5 - 4.5 - 4),
                   N_houshold = x0,
                   I0 = theta0[n_change+2]*6,
                   beta_b = theta0[1:n_change],
                   t_Btw_change = c(1,seq(0,length(times), 
                                          length.out = n_change)[-1]) )
  
  # running within and between
  output <- run.syst(input,times,x0)
  Traj[1,] <- output$D_expect
  Rt[1,] <- c(theta0[1],diff(output$out_b[,2]))
  # lambda <- lambda_fct(param = theta0, death, output)

  L1 <- Like1(output$D_expect,death,theta0)

  L[1,] <- L1
  thetas[1,] <- theta0       
  
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
        # get the new 'death'
        input$beta_w = Ts[n_change+1]
        input$beta_b = Ts[1:n_change]
        input$I0 <- Ts[n_change+2]*6
        # running within and between
        output_T <- run.syst(input,times,x0)
        # lambda <- lambda_fct(param = Ts, death, output)
        
        # get the likelihood for proposae value
        
        Lint <- Like1(output_T$D_expect,death, theta = Ts)
        
        #get the ratio with previous value of parameter and correct for porposal (and, only for initial conditions, prior distribution)
        r <- exp(Lint-L1)*Ts[j]/theta0[j]
      }
      
      
      # accept or reject
    
        if (runif(1,0,1) <= r){
          theta0 <- Ts  # if accept, keep new parameter value
          output <- output_T
          L1 <- Lint         # if accept, keep new lieklihood
        } 
      
      L[i,j] <- L1    # store final likelihood
    }
    
    thetas[i,] <- theta0  # store final parameter values for this iteration
    Traj[i,] <- output$D_expect
    Rt[i,] <- c(theta0[1],diff(output$out_b[,2]))
  }
  
  ####
  res <- list(theta = thetas, logL = L, Traj = Traj, Rt = Rt)
  return(res)
}
