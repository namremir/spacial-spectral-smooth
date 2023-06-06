library(tidyverse)
library(Matrix)
library(reshape2)

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
  c <- 2*pi*outer(seq(0,res)/res-1/2, seq(-n,n))
  inv <- cos(c) + 1i*sin(c)
  return(inv)
}

sssmooth <- function(d, lam, res=64, truncation_factor=10){
  h <- sl_evect(d, lam, truncation_factor=truncation_factor)
  ss <- dtft_inv(d*truncation_factor, res=res) %*% h
  return(ss)
}

optsssAndrews <- function(x, y, cs, labels, alpha=0.1, res=256, truncation_factor=10){
  n <- nrow(x)
  d <- ncol(x)
  
  t <- seq(0,1,length.out=res+1)
  
  svd <- svd(x)
  f <- t((sssmooth(d, 1/alpha, res=res, truncation_factor=truncation_factor) * 
            matrix(rep(svd$d,each=res+1),nrow=res+1)) %*% t(svd$u))
  f <- Re(f) + Im(f)
  
  f_long <- melt(t(f))
  f_long$t <- rep(t,times=n)
  f_long$cat <- rep(y,each=res+1)
  
  #plot 1
  plot1 <- ggplot(f_long)+
    geom_line(aes(x=t,y=value,group=factor(Var2),color=cat),alpha=.5)+
    scale_color_manual(values=cs)+
    labs(color="")+
    theme_minimal()
  
  minmax = data.frame(cat = rep(labels,each=res+1))
  minmax$t = rep(t,times=length(labels))
  
  min <- c()
  for (label in labels){
    for (n in 1:(res+1)){
      min <- append(min, min(f_long$value[f_long$Var1==n & f_long$cat==label]))
    }
  }  
  
  max <- c()
  for (label in labels){
    for (n in 1:(res+1)){
      max <- append(max, max(f_long$value[f_long$Var1==n & f_long$cat==label]))
    }
  }
  
  minmax$min <- min
  minmax$max <- max
  
  #plot 2
  plot2 <- ggplot(minmax)+
    geom_ribbon(aes(x=t,ymin=min,ymax=max,group=cat,color=cat,fill=cat), alpha=.5)+
    scale_fill_manual(values=cs)+
    scale_color_manual(values=cs)+
    labs(fill="",color="")+
    theme_minimal()
  
  list(plot1,plot2)
}




data(iris)
x <- iris[,seq(1,4)]  
y <- iris[,5]
labels <- levels(y)
cs <- c('#e41a1c', '#377eb8', '#4daf4a')
optsssAndrews(x,y,cs,labels,alpha=100)
optsssAndrews(x,y,cs,labels,alpha=.01)


bc <- read.table("wdbc.data",sep=",")
x <- bc[,seq(3,32)]
y <- bc[,2]
y <- as.factor(y)
labels <- levels(y)
cs = c( '#377eb8', '#e41a1c')
optsssAndrews(x, y, cs, labels, alpha=100) # Andrews
optsssAndrews(x, y, cs, labels, alpha=0.01) # spatial-spectral smooth


