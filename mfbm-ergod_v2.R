# used packages
require("somebm") # package for generating brownian motion
library("writexl") # package for writing Excel files

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

# function to generate values of wrho
wrho <- function(j, par){
  res <- par$k^2*par$h^(2*par$H)*0.5*(abs(j+1)^(2*par$H) - 2*abs(j)^(2*par$H) + abs(j-1)^(2*par$H)) + par$s^2*par$h*(j==0)
  return(res)
}

#series convergency
sum_wrho <- function(d1, d2, par){
  n_it <- 0
  delta <- 10
  s <- 0
  s_prv <- 0
  r_cur <- 10
  rho <- c()
  while (TRUE) {
    if (n_it == 100000001) {
      print("Series convergent issue")
      break
    }
    s_prv <- s
    r_cur <- wrho(n_it+d1, par)*wrho(n_it+d2,par)
    s <- s + r_cur
    if ((n_it >2)&(abs(s - s_prv)<10^(-4))) {
      break
    }
    n_it <- n_it+1
  }
  return(s)
}

#part of asymptotic matrix wrho
Z_11 <- function(par){
  res <- 2*sum_wrho(0,0,par) + 2*sum_wrho(1,1,par)
  return(res)
}
#part of asymptotic matrix wrho
Z_12 <- function(par){
  res <- 4*sum_wrho(0,1,par)
  return(res)
}
#part of asymptotic matrix wrho
Z_22 <- function(par){
  res <- sum_wrho(0,0,par) + sum_wrho(1,1,par) + sum_wrho(1,-1,par)+sum_wrho(2-1,1-(-1),par)
  return(res)
}
#part of asymptotic matrix wrho
Z_33 <- function(par){
  res <- 0
  for (a in c(0,1)) {
    for (b in c(2,3)) {
      for (c in c(0,1)) {
        for (d in c(2,3)) {
          for (gamma in c(0,1)) {
            res <-res + sum_wrho(c-a-gamma*(b-a),d-b-gamma*(a-b),par) + sum_wrho(1+(a-c)-gamma*(a-b),1+(b-d)-gamma*(b-a),par)
          }
        }
      }
    }
  }
  return(res)
}
#part of asymptotic matrix wrho
Z_13 <- function(par){
  s <- 0
  for (a in c(0,1)) {
    for (b in c(2,3)) {
      s <- s + sum_wrho(a,b,par) + sum_wrho(1-a,1-b,par)
    }
  }
  res <- 2*s
  return(res)
}
#part of asymptotic matrix wrho
Z_23 <- function(par){
  res <- 0
  for (a in c(0,1)) {
    for (b in c(2,3)) {
      for (gamma in c(0,1)) {
        res <-res + sum_wrho(a-gamma,b-1+gamma,par) + sum_wrho(1 - (a-gamma),1 - (b-1+gamma),par)
      }
    }
  }
  return(res)
}
#part of asymptotic matrix wrho
Z_rho <- function(par){
  Z <- t(matrix(c(Z_11(par),Z_12(par),Z_13(par),
                  Z_12(par),Z_22(par),Z_23(par),
                  Z_13(par),Z_23(par),Z_33(par)),
                nrow = 3, ncol = 3))
  return(Z)
}
#matrix of partial derivatives
g_der <- function(par){
  g11 <- 0
  g12 <- (-1)/(2*par$k^2*par$h^(2*par$H)*(2^(2*par$H-1)-1)*log(2))
  g13 <- (1) /(2*par$k^2*par$h^(2*par$H)*2^(2*par$H)*(2^(2*par$H-1)-1)*log(2))
  g21 <- 0
  g22 <-(2/par$h^(2*par$H))*(((2+log2(par$h))*2^(2*par$H) 
                          - 2*(log2(par$h)+1))/(2^(2*par$H)-2)^2)  
  g23 <-(2/par$h^(2*par$H))*((2*log2(par$h)-(log2(par$h)+1) * 
                            (2^(2*par$H)))/(2^(2*par$H)*(2^(2*par$H)-2)^2)) 
  g31 <- 1/(par$h)
  g32 <- (4*(1-2^(2*par$H)))/(par$h*(2^(2*par$H)-2)^2)
  g33 <- 2/(par$h*(2^(2*par$H)-2)^2)
  g_der <- t(matrix(c(g11,g12,g13,
                      g21,g22,g23,
                      g31,g32,g33),
                    nrow = 3, ncol = 3))
  return(g_der)
}
#main asymptotic covariance matrix
Z_asm <- function(par){
  g   <- g_der(par)
  Z_r <- Z_rho(par)
  Z <- g %*% Z_r %*% t(g)
  return(Z)
}


# estimators of all parameters
estimators <- function(ts, k, n, h, H, sigma, kappa, alpha = 0.05) {
  if (k > n) {
    print('There are not enough number of points! Use less k.')
    return(NaN)
  } else {
    ksi_n <- sum((ts[2:(k+1)] - ts[1:(k)])^2)/k
    eta_n <- sum((ts[2:(k+1)] - ts[1:(k)])*(ts[3:(k+2)] - ts[2:(k+1)]))/k
    zeta_n <- sum((ts[3:(k+2)] - ts[1:(k)])*(ts[5:(k+4)] - ts[3:(k+2)]))/k
    ksi  <- sigma^2*h + kappa^2*h^(2*H)
    eta  <- kappa^2 * h^(2*H)*(2^(2*H-1)-1)
    zeta <- kappa^2 * h^(2*H)*(2^(2*H-1)-1)*2^(2*H)
    
    
    H_n <- log2.nan(zeta_n/eta_n)/2
    kappa_n <- sqrt.nan(eta_n/h^(2*H_n)/(2^(2*H_n-1)-1))
    sigma_n <- sqrt.nan((ksi_n-kappa_n^2*h^(2*H_n))/h)
  
    if ((is.nan(H_n))|(is.nan(kappa_n))|(is.nan(sigma_n))) {
      return(c(H_n,  kappa_n^2,  sigma_n^2,  ksi_n, eta_n,  zeta_n,
               NaN, NaN, NaN, NaN, NaN, NaN,
               NaN, NaN, NaN, NaN, NaN, NaN, 
               1, 0))
    } 
    else if ((H_n<=0)|(H_n>=0.75)){
      return(c(H_n,  kappa_n^2,  sigma_n^2,  ksi_n, eta_n,  zeta_n,
               NaN, NaN, NaN, NaN, NaN, NaN,
               NaN, NaN, NaN, NaN, NaN, NaN, 
               0, 1))
      }
    else{{
      par <- list(H = H_n, k = kappa_n,
                  s = sigma_n, h = h)
      g   <- g_der(par)
      Zrho_n <- Z_rho(par)
      Zasm_n <- g %*% Zrho_n %*% t(g)

      zr_xi    <- sqrt.nan(Zrho_n[1,1]/k)
      zr_eta   <- sqrt.nan(Zrho_n[2,2]/k)
      zr_zeta  <- sqrt.nan(Zrho_n[3,3]/k)
      za_H     <- sqrt.nan(Zasm_n[1,1]/k)
      za_kappa <- sqrt.nan(Zasm_n[2,2]/k)
      za_sigma <- sqrt.nan(Zasm_n[3,3]/k)
      
      q = qnorm(1-alpha/2)
      cp_H    = abs(H-H_n) < q*za_H
      
      cp_kappa= (kappa^2-kappa_n^2< q*za_kappa)&
        (kappa^2-kappa_n^2>-q*za_kappa)
      
      cp_sigma= (sigma^2-sigma_n^2< q*za_sigma)&
        (sigma^2-sigma_n^2>-q*za_sigma)
      
      cp_xi   = (ksi-ksi_n<q*zr_xi)    &(ksi-ksi_n  >-q*zr_xi)
      cp_eta  = (eta-eta_n<q*zr_eta)   &(eta-eta_n  >-q*zr_eta)
      cp_zeta = (zeta-zeta_n<q*zr_zeta)&(zeta-zeta_n>-q*zr_zeta)

      return(c(H_n,  kappa_n^2,  sigma_n^2,  ksi_n, eta_n,  zeta_n,
               za_H, za_kappa, za_sigma, zr_xi, zr_eta, zr_zeta,
               cp_H, cp_kappa, cp_sigma, cp_xi, cp_eta, cp_zeta,
               0, 0))
    }
    }
  }
}

# parameters of the model
#H <- 0.35
h <- 1
#sigma <- 1.5
#kappa <- 2.5

n <- 2^20  # horizon of observations

n_iter <- 1000 # number of simulations for one set of parameters


columns <- c('n_trajectories', 'h', 'sigma', 'kappa', 'H',
             'estimator', 'measure',
             paste0('n_2pow', c(8,10,12,14,16,18,20)))

for (sigma in c(1.5)) {
  for (kappa in c(2.5)) {
    for (H in c(0.1,0.2,0.3,0.4,0.6,0.7)) {
      print(paste('sigma =', sigma, ', kappa =', kappa, ', H =', H))
      results <- c()
      
      # matrices for all estimators to write estimates for each simulation with various partitions
      matHn <- c()
      matsigman <- c()
      matkappan <- c()
      matxin <- c()
      matetan <- c()
      matzetan <- c()
      # variances
      matHn_z <- c()
      matsigman_z <- c()
      matkappan_z <- c()
      matxin_z <- c()
      matetan_z <- c()
      matzetan_z <- c()
      #cover probability
      matHn_cp <- c()
      matsigman_cp <- c()
      matkappan_cp <- c()
      matxin_cp <- c()
      matetan_cp <- c()
      matzetan_cp <- c()
      
      #Inappropriate values  
      matH_in <- c()
      matNoNan <- c()

        for (ii in 1:n_iter) {
          print(paste('Running simulation number', ii, '...'))
          st <- Sys.time()
          Xkn <- generate.Xkn(n, h, H, sigma, kappa) # generate process X
          
          # vectors for all estimators to write estimates for one particular simulation with various partitions
          vecest <- c()
          for (j in c(2^(8),2^(10),2^(12),2^(14),2^(16),2^(18),2^(20))) {
            vecest <- rbind(vecest, estimators(Xkn, j, n, h,
                                               H, sigma, kappa))
          }
          
          matHn     <- rbind(matHn, vecest[,1])
          matkappan <- rbind(matkappan, vecest[,2])
          matsigman <- rbind(matsigman, vecest[,3])
          matxin    <- rbind(matxin, vecest[,4])
          matetan   <- rbind(matetan, vecest[,5])
          matzetan  <- rbind(matzetan, vecest[,6])
          
          matHn_z     <- rbind(matHn_z, vecest[,7])
          matkappan_z <- rbind(matkappan_z, vecest[,8])
          matsigman_z <- rbind(matsigman_z, vecest[,9])
          matxin_z    <- rbind(matxin_z, vecest[,10])
          matetan_z   <- rbind(matetan_z, vecest[,11])
          matzetan_z  <- rbind(matzetan_z, vecest[,12])
          
          matHn_cp     <- rbind(matHn_cp, vecest[,13])
          matkappan_cp <- rbind(matkappan_cp, vecest[,14])
          matsigman_cp <- rbind(matsigman_cp, vecest[,15])
          matxin_cp    <- rbind(matxin_cp, vecest[,16])
          matetan_cp   <- rbind(matetan_cp, vecest[,17])
          matzetan_cp  <- rbind(matzetan_cp, vecest[,18])
          
          matNoNan  <- rbind(matNoNan, vecest[,19])
          matH_in   <- rbind(matH_in, vecest[,20])
          et <- Sys.time()
          print(et-st)
        }
        
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'H_n', 'mean', format(apply(matHn, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'H_n', 'std dev', format(apply(matHn, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'sigma_n', 'mean', format(apply(matsigman, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'sigma_n', 'std dev', format(apply(matsigman, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'kappa_n', 'mean', format(apply(matkappan, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'kappa_n', 'std dev', format(apply(matkappan, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
        
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'xi_n', 'mean', format(apply(matxin, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'xi_n', 'std dev', format(apply(matxin, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'etan_n', 'mean', format(apply(matetan, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'etan_n', 'std dev', format(apply(matetan, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'zeta_n', 'mean', format(apply(matzetan, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'zeta_n', 'std dev', format(apply(matzetan, 2, sd, na.rm = TRUE), digits = 4, nsmall = 4)))
        
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'H_n asigma', 'mean', format(apply(matHn_z, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'sigma_n asigma', 'mean', format(apply(matsigman_z, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'kappa_n asigma', 'mean', format(apply(matkappan_z, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'xi_n asigma', 'mean', format(apply(matxin_z, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'etan_n asigma', 'mean', format(apply(matetan_z, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'zeta_n asigma', 'mean', format(apply(matzetan_z, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'H_n CP', 'mean', format(100*apply(matHn_cp, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'sigma_n CP', 'mean', format(100*apply(matsigman_cp, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'kappa_n CP', 'mean', format(100*apply(matkappan_cp, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'xi_n CP', 'mean', format(100*apply(matxin_cp, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'etan_n CP', 'mean', format(100*apply(matetan_cp, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, 'zeta_n CP', 'mean', format(100*apply(matzetan_cp, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, '% of NaN', 'mean', format(100*apply(matNoNan, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, c(n_iter, h, sigma, kappa, H, '% H out (0,1)', 'mean', format(100*apply(matH_in, 2, mean, na.rm = TRUE), digits = 4, nsmall = 4)))
        results <- rbind(results, rep('==========', dim(results)[2]))
        
        name <- paste("\\Estim. H=",H, ", sigma=", sigma,
                      ", kappa=",kappa, ".xlsx", sep="")
        colnames(results) <- columns
        
        file <- paste("Z:\\University\\Ральченко\\2022-06-10 Асимптотична нормальність\\Modeling\\Update\\results", name, sep="")
        res <- data.frame(results)
        write_xlsx(res, file)
      }
  }
}