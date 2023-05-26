# used packages
require("somebm") # package for generating brownian motion
library("writexl") # package for writing Excel files
library("pracma") # package for a function minimizing 

# function to generate process X
generate.Xkn <- function(n, h, H, sigma, kappa) {
  tkn <- (0:(n+3)) * h # time partition
  Wkn <- bm(x0 = 0, t0 = 0, t = (n+3)*h, n = (n+3)) # path of Bm on [0, 2^m] with step 1/2^n
  Bkn <- fbm(hurst = H, n = (n+3)) # path of fBm on [0, 1]
  Bkn <- ts(((n+3)*h)^(H)*Bkn, start = tkn[1], end = tkn[length(tkn)], frequency = 1/h) # rescale fBm to [o, 2^m]
  Xkn <- sigma * Wkn + kappa * Bkn # considered process X
  return(Xkn)
}

# function \log_{2+}
log2.nan <- function(x) {
  if (x > 0) {
    log2(x)
  } else {
    NaN
  }
}

# function \sqrt{+}
sqrt.nan <- function(x) {
  if (is.nan(x)) {
    NaN
  } else if (x > 0) {
    sqrt(x)
  } else {
    NaN
  }
}

# covariance matrix
Gamma <- function(Y, sigma, kappa, H, h)
{ N <- length(Y)
  G <- matrix(0, N,N)
  hs2 <- h*sigma^2
  hk2 <- h^(2*H)*kappa^2/2
  H2 <- 2*H
  for (i in 1:N) {
    for (j in 1:N) {
      G[i,j] <- hs2*(min(i,j)) + hk2*(i^H2+j^H2-abs(i-j)^H2)
    }
  }
  return(G)
}

# negative log likelihood
mLF <- function(val, hs = h, Y = Xkn){
  N <- length(Y)
  G <- Gamma(Y,sqrt(val[2]), sqrt(val[3]), val[1], hs)
  if (det(G)!= 0) {
    like <- (1/sqrt(det(G)*(2*pi)^N))*exp((-0.5)*t(Y) %*% solve(G) %*% Y)
  }
  else {like <- NaN}
  return(-log(like))
}

# function to generate process X
generate.Xkn <- function(n, h, H, sigma, kappa) {
  tkn <- (0:(n+3)) * h # time partition
  Wkn <- bm(x0 = 0, t0 = 0, t = (n+3)*h, n = (n+3)) # path of Bm on [0, 2^m] with step 1/2^n
  Bkn <- fbm(hurst = H, n = (n+3)) # path of fBm on [0, 1]
  Bkn <- ts(((n+3)*h)^(H)*Bkn, start = tkn[1], end = tkn[length(tkn)], frequency = 1/h) # rescale fBm to [o, 2^m]
  Xkn <- sigma * Wkn + kappa * Bkn # considered process X
  return(Xkn)
}

# estimators of all parameters
estimators <- function(ts, k, n, h) {
  if (k > n) {
    print('There are not enough number of points! Use less k.')
    return(NaN)
  } else {
    ksi_n <- sum((ts[2:(k+1)] - ts[1:(k)])^2)/k
    eta_n <- sum((ts[2:(k+1)] - ts[1:(k)])*(ts[3:(k+2)] - ts[2:(k+1)]))/k
    zeta_n <- sum((ts[3:(k+2)] - ts[1:(k)])*(ts[5:(k+4)] - ts[3:(k+2)]))/k
    
    H_n <- log2.nan(zeta_n/eta_n)/2
    kappa_n <- sqrt.nan(eta_n/h^(2*H_n)/(2^(2*H_n-1)-1))
    sigma_n <- sqrt.nan((ksi_n-kappa_n^2*h^(2*H_n))/h)
    return(c(H_n,  kappa_n^2,  sigma_n^2))
  }
}

# parameters of the model
#H <- 0.7
h <- 1
#sigma <- 1.5
#kappa <- 2.5

n <- 2^7  # horizon of observations

n_iter <- 1000 # number of simulations for one set of parameters

columns <- c('n_trajectories', 'h', 'sigma', 'kappa', 'H', 'type',
             'estimator', 'measure', 'n_2pow7')

for (sigma in c(1.5)) {
  for (kappa in c(2.5)) {
    for (H in c(0.1,0.2,0.3,0.4,0.6, 0.7)) {
      print(paste('sigma =', sigma, ', kappa =', kappa, ', H =', H))
      results <- c()
      
      # matrices for all estimators to write estimates for each simulation with various partitions
      matHn_er <- c()
      matsigman_er <- c()
      matkappan_er <- c()
      mattime_er <- c()
      
      matHn_ml <- c()
      matsigman_ml <- c()
      matkappan_ml <- c()
      mattime_ml <- c()
      
      for (ii in 1:n_iter) {
        print(paste('Running simulation number', ii, '...'))
        st <- Sys.time()
        Xkn <- generate.Xkn(n, h, H, sigma, kappa) # generate process X
        st_er <- Sys.time()
        # vectors for all estimators to write estimates for one particular simulation with various partitions
        vecest <- c()
        for (j in c(2^(7))) {
          vecest <- rbind(vecest, estimators(Xkn, j, n, h))
        }
        et_er <-  Sys.time()
        
        matHn_er     <- rbind(matHn_er, vecest[,1])
        matkappan_er <- rbind(matkappan_er, vecest[,2])
        matsigman_er <- rbind(matsigman_er, vecest[,3])
        mattime_er   <- rbind(mattime_er, et_er-st_er)
        
        st_ml <- Sys.time()
        min <- fminsearch(mLF, c(0.4,2,2), method = "Hooke-Jeeves", maxiter = 10000, 
                          lower = c(10^(-10),10^(-10),10^(-10)), upper = c(1, Inf, Inf))
        et_ml <-  Sys.time()
         
        matHn_ml     <- rbind(matHn_ml, min$xmin[1])
        matkappan_ml <- rbind(matkappan_ml, min$xmin[3])
        matsigman_ml <- rbind(matsigman_ml, min$xmin[2])
        mattime_ml   <- rbind(mattime_ml, et_ml-st_ml)
       
        et <- Sys.time()
        print(et-st)
      }
      
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'ERG', 'H_n', 'mean', format(apply(matHn_er, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'ERG', 'H_n', 'std dev', format(apply(matHn_er, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'ERG', 'sigma_n', 'mean', format(apply(matsigman_er, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'ERG', 'sigma_n', 'std dev', format(apply(matsigman_er, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'ERG', 'kappa_n', 'mean', format(apply(matkappan_er, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'ERG', 'kappa_n', 'std dev', format(apply(matkappan_er, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'ERG', 'time', 'mean', format(apply(mattime_er, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'MLE', 'H_n', 'mean', format(apply(matHn_ml, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'MLE', 'H_n', 'std dev', format(apply(matHn_ml, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'MLE', 'sigma_n', 'mean', format(apply(matsigman_ml, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'MLE', 'sigma_n', 'std dev', format(apply(matsigman_ml, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'MLE', 'kappa_n', 'mean', format(apply(matkappan_ml, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'MLE', 'kappa_n', 'std dev', format(apply(matkappan_ml, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'MLE', 'time', 'mean', format(apply(mattime_ml, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
      results <- rbind(results, rep('==========', dim(results)[2]))
      
      name <- paste("\\Estim. H=",H, ", sigma=", sigma,
                    ", kappa=",kappa, ".xlsx", sep="")
      colnames(results) <- columns
      
      file <- paste("Z:\\University\\Ральченко\\2022-06-10 Асимптотична нормальність\\Modeling\\Update\\mle", name, sep="")
      res <- data.frame(results)
      write_xlsx(res, file)
    }
  }
}