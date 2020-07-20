#' get likelihood  
#'
#' get Poisson likelihood  
#' internal to the MCMC
#' 
#' @param lambda: 'force of infection' matrix (incidence weighted by serial interval),
#'                  column number of days, row: number of locations   
#'                  
#' @param I matrix of observed incidence, same dimension as lambda
#' 
#' @param R0 vector of reproduction numbers per locations
#'
#' @details  L log likelihood
#' @export
#' 

Like1<-function(lambda,death,theta){
  
  # L <- sum(dpois(x = death,lambda = lambda,log = TRUE))
  L <- sum(dnbinom(x = death,mu = lambda, size = theta[n_change+3], log = TRUE))
  
  return(L)
}

