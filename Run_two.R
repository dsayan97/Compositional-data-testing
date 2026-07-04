library(MASS)
library(zCompositions)
source("EBC.R")

set.seed(1)
## ------------------------------------------------------------
## User controls: change these for each table
## ------------------------------------------------------------
n1 <- 50 # sample sizes
n2 <- 60 # sample sizes
p <- 200 # dimension
d <- 1 # delta
dist <- 1  # 1 = Normal, 2 = t5
COV <- 1 # Setups M.1, M.2, M.3
s <- 0.1 #sparcity

alpha <- 0.05 # level
zero_extra <- 0.30 # eta - excess zeroes


ITER <- 30 # Monte Carlo iterations

impute_method <- "CZM"
## ------------------------------------------------------------
## helper functions
## ------------------------------------------------------------
make_L_matrices <- function(p, COV) {
  pd <- floor(p / 4)
  
  if (COV %in% c(1, 3)) {
    r1 <- 0.999; r2 <- 0.9; r3 <- 0.7; r4 <- 0.5
    
    A1 <- (1 - r1) * diag(pd) + matrix(r1, pd, pd)
    A2 <- (1 - r2) * diag(pd) + matrix(r2, pd, pd)
    A3 <- (1 - r3) * diag(pd) + matrix(r3, pd, pd)
    A4 <- (1 - r4) * diag(pd) + matrix(r4, pd, pd)
    
    D2 <- matrix(0, pd, pd)
    D1 <- matrix(1 / pd, pd, pd)
    D3 <- matrix(-1 / pd, pd, pd)
    
    EC1 <- rbind(
      cbind(A1, D1, D2, D3),
      cbind(D1, A2, D1, D2),
      cbind(D2, D1, A3, D1),
      cbind(D3, D2, D1, A4)
    )
    
    eig1 <- eigen(EC1, symmetric = TRUE)
    L1 <- eig1$vectors %*% diag(sqrt(pmax(eig1$values, 0)))
    
    if (COV == 1) {
      L2 <- L1
    } else {
      EC2 <- rbind(
        cbind(A4, D1, D2, D3),
        cbind(D1, A3, D1, D2),
        cbind(D2, D1, A2, D1),
        cbind(D3, D2, D1, A1)
      )
      
      eig2 <- eigen(EC2, symmetric = TRUE)
      L2 <- eig2$vectors %*% diag(sqrt(pmax(eig2$values, 0)))
    }
  }
  
  if (COV == 2) {
    r1 <- 0.5; r2 <- 0.3; r3 <- 0.1
    
    A1 <- toeplitz(r1^(0:(pd - 1)))
    A2 <- toeplitz(r2^(0:(pd - 1)))
    A3 <- toeplitz(r3^(0:(pd - 1)))
    A4 <- diag(pd)
    
    D2 <- matrix(0, pd, pd)
    D1 <- matrix(1 / pd, pd, pd)
    D3 <- matrix(-1 / pd, pd, pd)
    
    EC1 <- rbind(
      cbind(A1, D1, D2, D3),
      cbind(D1, A2, D1, D2),
      cbind(D2, D1, A3, D1),
      cbind(D3, D2, D1, A4)
    )
    
    eig1 <- eigen(EC1, symmetric = TRUE)
    L1 <- eig1$vectors %*% diag(sqrt(pmax(eig1$values, 0)))
    
    EC2 <- rbind(
      cbind(A4, D1, D2, D3),
      cbind(D1, A1, D1, D2),
      cbind(D2, D1, A2, D1),
      cbind(D3, D2, D1, A3)
    )
    
    eig2 <- eigen(EC2, symmetric = TRUE)
    L2 <- eig2$vectors %*% diag(sqrt(pmax(eig2$values, 0)))
  }
  
  list(L1 = L1, L2 = L2)
}

get_NB_params <- function(p) {
  if (p == 200) return(list(mu = 3.0 * p, theta = 0.05 * p))
  if (p == 500) return(list(mu = 1.5 * p, theta = 0.02 * p))
  if (p == 1000) return(list(mu = 1.5 * p, theta = 0.02 * p))
}

one_rep <- function(task, mu) {
  n1 <- task[1]
  n2 <- task[2]
  p <- task[3]
  d <- task[4]
  dist <- task[5]
  COV <- task[6]
  s <- task[7]
  
  Ls <- make_L_matrices(p, COV)
  L1 <- Ls$L1
  L2 <- Ls$L2
  
  G <- diag(p) - matrix(1, p, p) / p
  
  Cc1 <- matrix(0, n1, p)
  Cc2 <- matrix(0, n2, p)
  
  m1 <- rep(0, p)
  
  m2 <- rep(0, p)
  m2[1:floor(p * s)] <- mu
  
  nb <- get_NB_params(p)
  N1 <- pmax(p, rnegbin(n1, mu = nb$mu, theta = nb$theta))
  N2 <- pmax(p, rnegbin(n2, mu = nb$mu, theta = nb$theta))
  
  drop1 <- matrix(rbinom(n1 * p, 1, zero_extra), n1, p)
  drop2 <- matrix(rbinom(n2 * p, 1, zero_extra), n2, p)
  
  gen_noise <- function() {
    if (dist == 1) {
      rnorm(p)
    } else {
      rt(p, 5) / sqrt(5 / 3)
    }
  }
  
  for (i in seq_len(n1)) {
    Z <- as.numeric(m1 + L1 %*% gen_noise())
    Pi <- exp(Z - max(Z))
    Pi <- Pi / sum(Pi)
    Cc1[i, ] <- rmultinom(1, N1[i], Pi)
  }
  
  for (i in seq_len(n2)) {
    Z <- as.numeric(m2 + L2 %*% gen_noise())
    Pi <- exp(Z - max(Z))
    Pi <- Pi / sum(Pi)
    Cc2[i, ] <- rmultinom(1, N2[i], Pi)
  }
  
  zero_before_1 <- mean(Cc1 == 0)
  zero_before_2 <- mean(Cc2 == 0)
  
  Cc1[drop1 == 1] <- 0
  Cc2[drop2 == 1] <- 0
  
  zero_after_1 <- mean(Cc1 == 0)
  zero_after_2 <- mean(Cc2 == 0)
  
  X1 <- cmultRepl(
    Cc1,
    method = impute_method,
    output = "prop",
    z.warning = 0.999,
    z.delete = FALSE,
    suppress.print = TRUE
  )
  
  X2 <- cmultRepl(
    Cc2,
    method = impute_method,
    output = "prop",
    z.warning = 0.999,
    z.delete = FALSE,
    suppress.print = TRUE
  )
  
  Y1 <- t(G %*% t(log(X1)))
  Y2 <- t(G %*% t(log(X2)))
  
  Y <- list(Y1, Y2)
  
  res_EBC <- as.integer(EBC(Y) < alpha)
  
  res_EBC
}

## ------------------------------------------------------------
## Run
## ------------------------------------------------------------
task <- c(n1, n2, p, d, dist, COV, s)

res <- numeric(ITER)
mu <- runif(floor(p * s), -d / 2, d / 2) 
for(i in 1:ITER){
  res[i] <- one_rep(task, mu)
}

mean(res)
