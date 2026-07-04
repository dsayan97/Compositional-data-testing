library(MASS)
library(zCompositions)
source("EBC.R")

set.seed(1)
## ------------------------------------------------------------
## User controls: change these for each table
## ------------------------------------------------------------
n1 <- 50 # sample sizes
n2 <- 60 # sample sizes
n3 <- 70 # sample sizes
n4 <- 80 # sample sizes
n5 <- 90 # sample sizes
p <- 200 # dimension
d <- 1 # delta
dist <- 1  # 1 = Normal, 2 = t5
s <- 0.1 #sparcity

alpha <- 0.05 # level
zero_extra <- 0.30 # eta - excess zeroes


ITER <- 30 # Monte Carlo iterations

impute_method <- "CZM"
## ------------------------------------------------------------
## helper functions
## ------------------------------------------------------------
make_L_matrices_K <- function(p) {
  pd <- floor(p / 4)
  
  D2 <- matrix(0, pd, pd)
  D1 <- matrix(1 / pd, pd, pd)
  D3 <- matrix(-1 / pd, pd, pd)
  
  make_CS <- function(A1, A2, A3, A4) {
    rbind(
      cbind(A1, D1, D2, D3),
      cbind(D1, A2, D1, D2),
      cbind(D2, D1, A3, D1),
      cbind(D3, D2, D1, A4)
    )
  }
  
  sqrt_mat <- function(Sigma) {
    eig <- eigen(Sigma, symmetric = TRUE)
    eig$vectors %*% diag(sqrt(pmax(eig$values, 0)))
  }
  
  EQ <- function(rho) (1 - rho) * diag(pd) + rho * matrix(1, pd, pd)
  
  A999 <- EQ(0.999)
  A9   <- EQ(0.9)
  A7   <- EQ(0.7)
  A5   <- EQ(0.5)
  
  Sigma1 <- make_CS(A999, A9,   A7,   A5)
  Sigma2 <- make_CS(A5,   A999, A9,   A7)
  Sigma3 <- make_CS(A7,   A5,   A999, A9)
  Sigma4 <- make_CS(A9,   A7,   A5,   A999)
  Sigma5 <- Sigma1
  
  list(
    L1 = sqrt_mat(Sigma1),
    L2 = sqrt_mat(Sigma2),
    L3 = sqrt_mat(Sigma3),
    L4 = sqrt_mat(Sigma4),
    L5 = sqrt_mat(Sigma5)
  )
}

get_NB_params <- function(p) {
  if (p == 200)  return(list(mu = 4.0 * p, theta = 0.10 * p))
  if (p == 500)  return(list(mu = 3.0 * p, theta = 0.08 * p))
  if (p == 1000) return(list(mu = 2.5 * p, theta = 0.06 * p))
  stop("Unsupported p.")
}

make_m <- function(a, b) {
  out <- rep(0, p)
  sidx <- sample(seq_len(p), floor(p * s))
  out[sidx] <- runif(length(sidx), a, b)
  out
}

one_rep <- function(task, mu_list) {
  nvec <- task[1:5]
  p <- task[6]
  d <- task[7]
  dist <- task[8]
  s <- task[9]
  
  Ls <- make_L_matrices_K(p)
  Llist <- list(Ls$L1, Ls$L2, Ls$L3, Ls$L4, Ls$L5)
  
  G <- diag(p) - matrix(1, p, p) / p
  
  mlist <- mu_list
  
  nb <- get_NB_params(p)
  
  gen_noise <- function() {
    if (dist == 1) {
      rnorm(p)
    } else {
      rt(p, 5) / sqrt(5 / 3)
    }
  }
  
  C_list <- vector("list", 5)
  zero_before <- numeric(5)
  zero_after <- numeric(5)
  
  for (k in 1:5) {
    nk <- nvec[k]
    Ck <- matrix(0, nk, p)
    
    Ndepth <- pmax(p, rnegbin(nk, mu = nb$mu, theta = nb$theta))
    dropk <- matrix(rbinom(nk * p, 1, zero_extra), nk, p)
    
    for (i in seq_len(nk)) {
      Z <- as.numeric(mlist[[k]] + Llist[[k]] %*% gen_noise())
      Pi <- exp(Z - max(Z))
      Pi <- Pi / sum(Pi)
      Ck[i, ] <- rmultinom(1, Ndepth[i], Pi)
    }
    
    zero_before[k] <- mean(Ck == 0)
    
    Ck[dropk == 1] <- 0
    
    zero_after[k] <- mean(Ck == 0)
    
    C_list[[k]] <- Ck
  }
  
  ## Zero replacement + CLR
  Y_list <- lapply(C_list, function(Ck) {
    Xk <- cmultRepl(
      Ck,
      method = impute_method,
      output = "prop",
      z.warning = 0.999,
      z.delete = FALSE,
      suppress.print = TRUE
    )
    
    t(G %*% t(log(Xk)))
  })
  
  ## K = 3, 4, 5 tests
  res_K3 <- as.integer(EBC(Y_list[1:3]) < alpha)
  res_K4 <- as.integer(EBC(Y_list[1:4]) < alpha)
  res_K5 <- as.integer(EBC(Y_list[1:5]) < alpha)
  
  c(res_K3, res_K4, res_K5)
}

## ------------------------------------------------------------
## Run simulation
## ------------------------------------------------------------
task <- c(n1, n2, n3, n4, n5, p, d, dist, s)

mlist <- list(
  rep(0, p),
  make_m(-d / 4,  d / 4),
  make_m(-d / 2,  d / 2),
  make_m(-3 * d / 4, 3 * d / 4),
  make_m(-d , d)
)

res <- matrix(0, ITER, 3)
for(i in 1:ITER){
  res[i,] <- one_rep(task, mlist)
}

colnames(res) <- c("K_1", "K_2", "K_3")
colMeans(res)
