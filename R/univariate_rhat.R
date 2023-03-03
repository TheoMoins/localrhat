#### Rhat functions

#' Compute the traditional Gelman-Rubin diagnostic.
#'
#' @param chains an array of size \eqn{n \times m} where \eqn{n} is the length of
#' the chains and \eqn{m \geq 2} is the number of chains.
#'
#' @return \eqn{\hat{R}} computed on the \eqn{m} chains.
#'
#' @details The function return the tradition Gelman-Rubin diagnostic \eqn{\hat{R}}.
#'
#' @references TO INCLUDE
#'
#'
trad_rhat <- function(chains) {
    dimm <- dim(chains)
    N <- dimm[1]
    M <- dimm[2]
    B <- 0
    W <- 0
    for (i in 1:M) {
        B <- B + (mean(chains[, i]) - mean(chains))^2
        W <- W + mean(chains[, i]^2) - mean(chains[, i])^2
    }
    # B <- (1/(M - 1)) * B
    # W <- (1/M) * W
    return(sqrt(1 + B/W))
}

############ UNIVARIATE VERSION ############


#' Local version of the Gelman Rubin diagnostic \eqn{\hat{R}}
#'
#' Compute \eqn{\hat{R}(x)}, a version of univariate \eqn{\hat{R}} computed
#' on indicator variables for a given quantile x.
#'
#' @param x a float number corresponding to the quantile used for the computation
#' of \eqn{\hat{R}(x)}.
#' @param chains an array of size \eqn{n \times m} where \eqn{n} is the length of
#' the chains and \eqn{m \geq 2} is the number of chains.
#' @param verbose a boolean indicating if a suggested threshold is printed with
#' the value of \eqn{\hat{R}(x)}. This threshold is associated to a type I error
#' of 0.05.
#'
#' @return The local-\eqn{\hat{R}(x)} computed on the \eqn{m} chains.
#'
#' @details The function return \eqn{\hat{R}} computed on \eqn{I(\theta^{(i,j)} \leq x)}
#' for a given value of \eqn{x}:
#' \deqn{\hat{R}(x) = \sqrt{\frac{\frac{n-1}{n} \hat{W}(x) + \hat{B}(x)}{\hat{W}(x)}},}
#' with \eqn{\hat{W}(x)} and \eqn{\hat{B}(x)} the estimated local within-chain variance
#' and between-chain variance:
#' \deqn{
#' \hat{W}(x) = \frac{1}{m} \sum_{j=1}^m (F_j(x) - F^2_j(x))
#' }
#' And
#' \deqn{
#' \hat{B}(x) = \frac{1}{m^2}\sum_{j<k} \left(F_j(x)-F_{k}(x)\right)^2.
#' }
#'
#' @references TO INCLUDE
#'
#' @examples
#' library(localrhat)
#'
#' N <- 500 # length of chains
#' M <- 4 # number of chains
#'
#' # Toy example with 3 i.i.d chains uniform in [-0.5, 0.5] and 1 in [-1, 1]:
#' chains <- array(c(runif((M-1)*N, -1/2, 1/2), runif(N, -1,1)), c(N,M))
#'
#' # Quantile to evaluate:
#' x <- 0.5
#' # local_rhat(x, chains)
#'
local_rhat <- function(x, chains, verbose=F) {
    res <- trad_rhat(chains <= x)
    if (is.na(res)){
        res <- 1
    }
    if (verbose){
        rlim <- rlim_x(m = dim(chains)[2])
        cat(paste("Threshold at confidence level 5%: ", rlim))
        cat(paste("\nLocal R-hat obtained: ", round(res, digits = 4)))
        pval <- p_value_r_x(res, dim(chains)[2])
        if (pval == 0){
            cat(paste("\np-value: < 0.001"))
        } else {
            cat(paste("\np-value: ", pval))
        }
        if (res > rlim){
            cat("\nWARNING: A convergence issue has been diagnosed")
        } else {
            cat("\nAt 5%, no convergence issues have been diagnosed\n")
        }
    }
    return(res)
}


#' Quantile values used for the computation of \eqn{\hat{R}_\infty}
#'
#' Return the set of points that will be used to estimate the supremum of
#' \eqn{\hat{R}(x)} over \eqn{x}. The number of different values taken by
#' \eqn{\hat{R}(x)} can not exceed the number of samples, so the number of points
#' used for the computation is the minimum between \eqn{n \times m} and a threshold
#' value that can be specified in the argument \code{max_nb_points}
#'
#' @param chains an array of size \eqn{n \times m} where \eqn{n} is the length of
#' the chains and \eqn{m \geq 2} is the number of chains.
#' @param max_nb_points the maximal length of the grid in the case where the total
#' number of samples is larger. By default, \code{max_nb_points = 500}.
#'
#' @return A list that contains the different \eqn{x} used for the computation of
#' \eqn{\hat{R}(x)}.
#'
#' @details As the computation of local-\eqn{\hat{R}} are based on indicator
#' values \eqn{I(\theta^{(i,j)} \leq x)} and because the chains are finite, it is
#' enough to evaluate these indicators in \eqn{x = \theta^{(i',j')}} to obtain all
#' the possible values of \eqn{\hat{R}(x)}. However, to keep the computation time
#' reasonable, if \eqn{n \times m} is greater than a given \code{max_nb_points},
#' the values are sorted and uniform thinning is apply (which means picking off every k
#' iterations) so that only \code{max_nb_points} is kept at the end.
#'
#'
#' @references TO INCLUDE
#'
#' @examples
#' library(localrhat)
#'
#' N <- 500 # length of chains
#' M <- 4 # number of chains
#'
#' # Toy example with 3 i.i.d chains uniform in [-0.5, 0.5] and 1 in [-1, 1]:
#' chains <- array(c(runif((M-1)*N, -1/2, 1/2), runif(N, -1,1)), c(N,M))
#'
#' # plot_local_r(grid_for_R(chains), all_local_rhat(chains),
#' #             xlim = c(-1, 1), ylim=c(0.999,1.1), title ="Gaussian distributions")
#'
grid_for_R <- function(chains, max_nb_points = 500) {
    if (length(chains) < max_nb_points | max_nb_points == "ALL") {
        grid <- sort(chains)
    } else {
        ind <- as.integer(seq(1, length(chains), length.out = max_nb_points))
        grid <- sort(chains)[ind]
    }
    return(grid)
}

#' Computation of \eqn{\hat{R}(x)} on a given grid
#'
#' Return the vector of \eqn{\hat{R}(x)} values that will be used to compute
#' the supremum \eqn{\hat{R}_\infty}.
#' The function use \code{grid_for_R} for the grid of values of \eqn{x}.
#'
#' @param chains an array of size \eqn{n \times m} where \eqn{n} is the length of
#' the chains and \eqn{m \geq 2} is the number of chains.
#' @param max_nb_points the maximal length of the grid in the case where the total
#' number of samples is larger. By default, \code{max_nb_points = 500}.
#'
#' @return A list that contains the different \eqn{\hat{R}(x)} values used for
#' the computation of \eqn{\hat{R}_\infty}.
#'
#' @details The computation of \eqn{\hat{R}_\infty} require to compute a supremum
#' of \eqn{\hat{R}(x)} values over the quantiles \eqn{x}. See the documentation
#' of \code{grid_for_R} to see how the choice of the different \eqn{x} is done.
#' Thus, this function combine \code{grid_for_R} for the computation of the grid
#' and \code{local_rhat} for the computation of \eqn{\hat{R}(x)} on a given
#' \eqn{x} to obtain all the values of \eqn{\hat{R}(x)} on the generated grid.
#'
#'
#' @references TO INCLUDE
#'
#' @examples
#' library(localrhat)
#'
#' N <- 500 # length of chains
#' M <- 4 # number of chains
#'
#' # Toy example with 3 i.i.d chains uniform in [-0.5, 0.5] and 1 in [-1, 1]:
#' chains <- array(c(runif((M-1)*N, -1/2, 1/2), runif(N, -1,1)), c(N,M))
#'
#' # plot_local_r(grid_for_R(chains), all_local_rhat(chains),
#' #              xlim = c(-1, 1), ylim=c(0.999,1.1), title ="Gaussian distributions")
#'
all_local_rhat <- function(chains, max_nb_points = 500) {
    grid <- grid_for_R(chains, max_nb_points)
    tab_rhat <- sapply(grid, function(x) local_rhat(x, chains))
    return(tab_rhat)
}
