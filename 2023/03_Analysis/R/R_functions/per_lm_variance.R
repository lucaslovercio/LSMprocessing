# From Nat Comms 2021 Goswami paper. Made a couple of changes.
library(viridis)
per_lm_variance <- function(shape.data){
  
  variances<-rowSums(apply(shape.data ,c(1,2),var))
  
  cols<-viridis(100)
  
  #calculate log rates:
  x=(log10(variances))
  xlims<-NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin=100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  variancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  
  
  
  variance_table <- tibble("Per_Lm_Variance" = variances, "Log_Variance" = x, "Variance_Colors" = variancecolors)
  return(variance_table)
}