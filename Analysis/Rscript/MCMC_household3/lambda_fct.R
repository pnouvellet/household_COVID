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

lambda_fct<-function(param, death, output){
  
  # peak <- param[n_change+2]
  # f <- which(output$D_expect == max(output$D_expect))
  # 
  # lambda <- output$D_expect[(f-peak+1):(f-peak+length(death))]
  # 
  # return(lambda)
}



