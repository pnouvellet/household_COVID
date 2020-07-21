#' 'force of infection' 
#'
#' return incidence weighted by serial interval for the time window of interest 
#' internal to the MCMC
#' 
#' @param param: 'force of infection' matrix (incidence weighted by serial interval),
#'                  column number of days, row: number of locations   
#'                  
#' @param I matrix of observed incidence, same dimension as lambda
#' 
#' @param N_l integer of  numbers of locations
#'
#' @param ws vector reversed serial interval distribution (output from SI_gamma_dist_EpiEstim reverse in the MCMC function) )
#'
#' @param SItrunc integer, threshold of serial interval distribution (see SI_gamma_dist_EpiEstim)
#'
#' 
#' @details lambda incidence weighted by serial interval
#' @export
#' 

lambda_fct<-function(n_s, theta){
  
  P <- matrix(NA,length(n_s)+1,length(n_s)+1)
  B <- theta[1]
  Q <- theta[2]
  # size 1:
  P[1,1] <- 1
  P[1:2,1+1] <- c(B,1-B) 
  for (k in 2:length(n_s)){
    for (j in 0:(k-1)){
      P[j+1,k+1] <- choose(n = k, k = j) * P[ 1+j, 1+j] * B^(k-j) *Q^(j*(k-j))
    }
    P[k+1,k+1] <- 1-sum(P[1:k,k+1])
  }
  return(P)
}



