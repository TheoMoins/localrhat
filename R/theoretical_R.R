############ THEORETICAL VERSION #############

local_r_dist <- function(x, dists){
  m <- length(dists)
  W <- 0
  B <- 0
  for (j in 1:m){
    W <- W + dists[[j]](x)* (1- dists[[j]](x))
    if (j < m){
      for (k in (j+1):m){
        B <- B + (dists[[j]](x) - dists[[k]](x))**2
      }
    }
  }
  B <- B * (1/m**2)
  W <- W * (1/m)
  return (sqrt((W+B)/W))
}

r_dist_values <- function(npoints, xlim, dists){
  grid <- seq(xlim[1], xlim[2], length.out = npoints)
  theoretical_r <- unlist(lapply(grid, (function(x) local_r_dist(x, dists))))
  return(array(c(grid, theoretical_r), c(npoints,2)))
}

max_r_dist_bivariate <- function(npoints, xlim, dists){
  grid1 <- seq(xlim[1], xlim[2], length.out = npoints)
  grid2 <- seq(xlim[1], xlim[2], length.out = npoints)
  max_r = 0
  for (x1 in grid1){
    for (x2 in grid2){
      max_r = max(max_r, local_r_dist(c(x1,x2), dists))
    }
  }
  return(max_r)
}
