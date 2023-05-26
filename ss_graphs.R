library(tidyverse)
library(Matrix)

sl_evect <- function(n, lam = 4 * pi^2, truncation_factor = 10) {
  # Computes the eigenvectors for the schrodinger laplacian
  m <- truncation_factor * n
  quad_diag <- 4 * (pi^2) * seq(-m, m)^2 + 2 * lam
  off_diag <- -lam * rep(1, 2 * m)
  
  td_matrix <- diag(quad_diag)
  td_matrix[row(td_matrix) - col(td_matrix) == 1] <- td_matrix[row(td_matrix) - col(td_matrix) == -1] <- off_diag
  
  eigen_result <- eigen(td_matrix)
  v <- eigen_result$vectors[, ncol(eigen_result$vectors)+1-seq_len(n)]
  return(v)
}

dtft_inv <- function(n, res = 64) {
  # Returns an inversion matrix
  c <- 2*pi*outer(seq(0,res+1)/res-1/2, seq(-n,n))
  inv <- cos(c) + 1i*sin(c)
  return(inv)
}

ssmooth <- function(d, lam, res=64, truncation_factor=10){
  h <- sl_evect(d, lam, truncation_factor=truncation_factor)
  ss <- dtft_inv(d*truncation_factor, res=res) %*% h
  return(ss)
}

optsssAndrews <- function(x, y, cs, labels, alpha=0.1, res=256, truncation_factor=10){
  
}
  
