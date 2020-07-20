# functions


#########################
############  within
within_m <- odin::odin ({
  
  
  # household dyn
  deriv(S[])   <-  - beta * S[i] * I[i] / N[i] 
  deriv(E[])   <-  beta * S[i] * I[i]  /N[i] - delta * E[i]
  deriv(I[])   <-  (1-mu) * delta * E[i]  - sigma * I[i]
  deriv(Id[])   <-  (mu) * delta * E[i]  - sigma * Id[i]
  deriv(H[])   <-   sigma * Id[i] - hs * H[i]
  deriv(R[])   <-   sigma * I[i]  
  deriv(D[])   <-  hs * H[i]
  
  # initial conditions
  initial(S[]) <- S0[i]
  initial(E[]) <- E0[i]
  initial(I[]) <- I0[i]
  initial(Id[]) <- Id0[i]
  initial(H[]) <- H0[i]
  initial(R[]) <- R0[i]
  initial(D[]) <- D0[i]
  
  # input
  S0[] <- user()
  E0[] <- user()
  I0[] <- user()
  Id0[] <- user()
  H0[] <- user()
  R0[] <- user()
  D0[] <- user()
  
  # paramters
  beta <- user()
  delta <- user()
  mu <- user()
  sigma <- user()
  hs <- user()
  N[] <- user()
  
  
  #dimension
  dim(N) <- user()
  n_size <- length(N)
  
  dim(S) <- n_size
  dim(E) <- n_size
  dim(I) <- n_size
  dim(Id) <- n_size
  dim(H) <- n_size
  dim(R) <- n_size
  dim(D) <- n_size
  
  dim(S0) <- n_size
  dim(E0) <- n_size
  dim(I0) <- n_size
  dim(Id0) <- n_size
  dim(H0) <- n_size
  dim(R0) <- n_size
  dim(D0) <- n_size
  
  config(base) <- "lv4"
  
})

#########################
############ between

Between_m <- odin::odin ({
  
  
  # between household dyn
  I_tot <- sum(I)
  
  deriv(S[])   <-  - beta * S[i] * I_tot / sum(N_h) 
  deriv(E[])   <-  beta * S[i] * I_tot / sum(N_h)  - delta * E[i]
  deriv(I[])   <-  delta * E[i] - sigma[i] * I[i]
  deriv(R[])   <-  sigma[i] * I[i]  
  deriv(O[])   <-  beta * S[i] * I_tot / sum(N_h)
  deriv(Bt)    <- beta
  
  # change in beta
  beta <- interpolate(beta_b_time,beta_b,'linear')
  
  # initial conditions
  initial(S[]) <- S0[i]
  initial(E[]) <- E0[i]
  initial(I[]) <- I0[i]
  initial(R[]) <- R0[i]
  initial(O[]) <- O0[i]
  initial(Bt) <- 0
  
  # input
  S0[] <- user()
  E0[] <- user()
  I0[] <- user()
  R0[] <- user()
  O0[] <- user()
  
  # paramters
  beta_b[] <- user()
  beta_b_time[] <- user()
  delta <- user()
  sigma[] <- user()
  N_h[] <- user()
  
  #dimension
  dim(beta_b) <-user()
  dim(beta_b_time) <- user()
  
  dim(N_h) <- user()
  n_size <- length(N_h)
  
  dim(sigma) <- n_size
  dim(S) <- n_size
  dim(E) <- n_size
  dim(I) <- n_size
  dim(R) <- n_size
  dim(O) <- n_size
  
  dim(S0) <- n_size
  dim(E0) <- n_size
  dim(I0) <- n_size
  dim(R0) <- n_size
  dim(O0) <- n_size
  
  config(base) <- "lv4"
 
  
})

#########################
############ solve within
run.syst <- function(input,times,x0){

  n_size <- length(input$N_houshold)
 
  # Specifying parameters
  pars_w <- list(
    # parameters
    delta    = input$delta,  # latency
    sigma  = input$sigma,  # infection
    beta    = input$beta_w,
    mu = input$mu, #mortality
    hs = input$hs,
    N = seq(1,n_size, by = 1),
    
    # initial conditions
    S0 = seq(1,n_size, by = 1) - 1,
    E0 = rep(1,n_size),
    I0 = rep(0,n_size),
    Id0 = rep(0,n_size),
    H0 = rep(0,n_size),
    R0 = rep(0,n_size),
    D0 = rep(0,n_size)
  )
  # for (i in 1:1e2){
  mod <- within_m(user = pars_w)

  out <- mod$run(times)
  # }
  

  # inverse of mean infectious time per households
  # integral(t*(R+D)/sum(R+D)) with R+D newly recovered or death ~ stop being infectious
  I_period <- rep(NA,n_size)
  for (i in 1:n_size){
    I_period[i] <- sum( (out[1:(length(times)-1),1]-0.5)*
                          diff(out[,1+5*n_size+i]+out[,1+6*n_size+i])/
                          ( out[length(times),1+5*n_size+i] + 
                              out[length(times),1+6*n_size+i]) ) - (1/input$delta)
                       
  # c(1/sum((out[,1]*diff(out[,]+out[,]) ),
  #            1/sum((out$time[1:(length(times)-1)]-0.5)*diff(out$R2)/(input$N2*(1-0.01))),
  #            1/sum((out$time[1:(length(times)-1)]-0.5)*diff(out$R3)/(input$N3*(1-0.01))),
  #            1/sum((out$time[1:(length(times)-1)]-0.5)*diff(out$R4)/(input$N4*(1-0.01))),
  #            1/sum((out$time[1:(length(times)-1)]-0.5)*diff(out$R5)/(input$N5*(1-0.01))),
  #            1/sum((out$time[1:(length(times)-1)]-0.5)*diff(out$R6)/(input$N6*(1-0.01))) )
  }
  
  
  #########################
  ########### solve between
  
  # Specifying parameters
  pars_b <- list(
    # parameters
    delta    = input$delta,  # latency
    sigma  = 1/I_period,  # infection
    beta_b    = input$beta_b,
    beta_b_time = input$t_Btw_change,
    N_h = input$N_houshold,
    
    # initial conditions
    S0 = input$N_houshold - input$I0/n_size,
    E0 = rep(input$I0/n_size,n_size),
    I0 = rep(0,n_size),
    R0 = rep(0,n_size),
    O0 = rep(0,n_size)
    )
  # for (i in 1:1e2){
  mod2 <- Between_m(user = pars_b)

  out_b <- mod2$run(times)
  # }
  
  # full dynamics
  out_f <- out
  out_f[,2:ncol(out)] <- 0 
  
  H_I <- round(out_b[,2+4*n_size+seq(1,n_size,1)])
  H_I <- rbind(rep(0,6),diff(as.matrix(H_I)))
  H_I[1,] <- input$I0/6
  # n_d <- nrow(H_I)
  # for(i in 1:6){
  #   f <- which(H_I[,i] > 0)
  #   for (j in f){
  #     out_f[j:n_d,1+(i-1)*6+c(1:6)] <- out_f[j:n_d,1+(i-1)*6+c(1:6)] + 
  #       H_I[j,i]*out[1:(n_d-j+1),1+(i-1)*6+c(1:6)]
  #   }
  #   
  #   out_f[,1+(i-1)*6+c(1)] <- (x0[i] - cumsum(H_I[,i]))*i
  # }
  #############################
  for(i in 1:n_size){
    f <- i+seq(1,n_size*7,by=6)
    a <- out[,f]
    b <- H_I[,i]
    c <- out_f[,f]
    # f <- which(b > 0)
    
    for (j in 2:7){ 
      temp <- a[,j]
      ac <- matrix(temp[use_temp],(nb_day+add_days),(nb_day+add_days))
      ac[is.na(ac)] <- 0
      c[,j] <- ac %*% b
    }
    c[,1] <- ( input$N_houshold[i] *i - rowSums(c[,2:7]) ) 
    out_f[,f] <- c
  }
  #########################################
  out_f[,-1] <- out_f[,-1]
  
  D_expect <- c(0,diff(rowSums(out_f[,1+n_size*6+seq(1,n_size,1)])))
  
  res <- list(D_expect = D_expect[-c(1:add_days)], 
              out_t = out_f, out_w = out, out_b = out_b, H_I = H_I,
              pars_w = pars_w, pars_b = pars_b)
  
  return(res)
  
}