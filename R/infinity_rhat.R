############ FUNCTIONS R-hat INFINITY ############


#' Scalar version of the local Gelman Rubin diagnostic \eqn{\hat{R}(x)}.
#'
#' Compute \eqn{\hat{R}_\infty}, a scalar summary of the function \eqn{\hat{R}(x)}
#' corresponding to the supremum over the quantiles \eqn{x}.
#'
#' @param chains an array of size \eqn{n \times m} in the univariate case and of
#' size \eqn{n \times d \times m} in the d-variate one, where \eqn{n} is the length of
#' the chains, \eqn{m \geq 2} is the number of chains and \eqn{d} is the dimension.
#' @param dir a vector specifying which indicator to use for the multivariate case.
#' See the function \code{multivariate_directed_local_rhat} for more details. If no
#' direction is given, the computation is done on the indicator variable with the
#' "\eqn{\leq}" sign on all dimension.
#' @param max_nb_points the maximal length of the grid in the case where the total
#' number of samples is larger. By default, \code{max_nb_points = 500}.
#'
#' @return The value of \eqn{\hat{R}_\infty} computed on the \eqn{m} chains.
#'
#' @details \eqn{\hat{R}_\infty} is based on \eqn{\hat{R}(x)}, the computation of
#' \eqn{\hat{R}} on \eqn{I(\theta^{(i,j)} \leq x)} (in the univariate case):
#' \deqn{R_\infty = \sup_{x \in R} R(x).}
#' This require to compute \eqn{\hat{R}(x)} on different values of \eqn{x} to
#' estimate this supremum. See the function \code{all_local_rhat} for more details.
#' In the multivariate case, \eqn{\hat{R}} can be computed on
#' \eqn{I\{\theta_1^{(j)} \leq x_1, \ldots, \theta_d^{(j)} \leq x_d\}},
#' but also with indicator where some "\eqn{\leq} are replaced by "\eqn{\geq}",
#' which leads to \eqn{2^{d-1}} possibilities and potentially as many different
#' results of \eqn{\hat{R}_\infty}. Here, a 0-1 vector can be given in the function
#' to specify which indicator to use.
#'
#' @references TO INCLUDE
#'
#' @examples
#' library(localrhat)
#' library("MASS") # For Multivariate Normal Distribution
#' rho <- 0.9
#' m <- 4
#' n <- 100
#' reps <- 50
#'
#' # Function to generate bivariate normal chains, with one that differs by its sigma:
#' gen_bvnormal_chains <- function(M, N, rho){
#'     sig_matrix <- (1-rho) * diag(2) + matrix(rho, nrow=2, ncol=2)
#'     return (array(c(mvrnorm((M-1)*N, mu = rep(0, 2), Sigma = diag(2)),
#'                     mvrnorm(N, mu = rep(0, 2), Sigma = sig_matrix)), c(N,2,M)))
#' }
#'
#' R_values <- c()
#' for (e in 1:reps){
#'     chains <- gen_bvnormal_chains(m, n, rho)
#'     # R_values <- c(R_values, rhat_infinity_max_directions(chains))
#' }
#' # R_mat <- matrix(data = R_values, ncol = 1)
#' # colnames(R_mat) <- "R-hat-inf_max_dir"
#'
#' # plot_hist(R_mat, bin_size = 0.004,
#' #           lim_y_axis = 15, vaxis_pos = 1.015,
#' #           plot_threshold = F)
#'
rhat_infinity <- function(chains, dir = NULL, max_nb_points = 500) {
    # Check the dimension to use the univariate or multivariate R-hat :
    if (length(dim(chains)) == 2) {
        val_rhat = all_local_rhat(chains, max_nb_points)
    } else {
        val_rhat = multivariate_all_local_rhat(chains, dir, max_nb_points)
    }
    return(max(val_rhat[is.finite(val_rhat)]))
}



#' Symmetrical version of \eqn{\hat{R}_\infty} in the multivariate case.
#'
#' Compute the multivariate version of \eqn{\hat{R}_\infty} in all possible
#' dependence direction and return the maximum values.
#'
#' @param chains an array of size \eqn{n \times d \times m}, where \eqn{n} is
#' the length of the chains, \eqn{m \geq 2} is the number of chains and \eqn{d}
#' is the dimension.
#' @param max_nb_points the maximal length of the grid in the case where the total
#' number of samples is larger. By default, \code{max_nb_points = 500}.
#' @param nb_directions an integer giving the number of directions to compute
#' \eqn{\hat{R}_\infty}: if given, it has to be between 1 and \eqn{2^{d-1}} and a
#' a random set of directions of this size will be used.
#'
#' @return The maximum value of all \eqn{\hat{R}_\infty} computed.
#'
#' @details In the multivariate case, the sensitivity of \eqn{R_\infty} strongly
#' depends on the sign of dependence of the variables. This function implement
#' the naive way to consider all direction, by computing on the \eqn{2^{d-1}}
#' possible one. \strong{Warning}: this method becomes extremely costly when the
#' dimension grow, and we recommend to avoid it as soon as \eqn{d > 5} (pending future
#' work on multidimensional improvements).
#' Alternatively, a number of directions can be specified, to avoid computational
#' issues and compute \eqn{R_\infty} on a subset of the possibilities: then,
#' random directions will be generated to diagnostic the convergence on the
#' corresponding dependence direction.
#'
#' @examples
#' library(localrhat)
#' library("MASS") # For Multivariate Normal Distribution
#' rho <- 0.9
#' m <- 4
#' n <- 100
#' reps <- 50
#'
#' # Function to generate bivariate normal chains, with one that differs by its sigma:
#' gen_bvnormal_chains <- function(M, N, rho){
#'     sig_matrix <- (1-rho) * diag(2) + matrix(rho, nrow=2, ncol=2)
#'     return (array(c(mvrnorm((M-1)*N, mu = rep(0, 2), Sigma = diag(2)),
#'                     mvrnorm(N, mu = rep(0, 2), Sigma = sig_matrix)), c(N,2,M)))
#' }
#'
#' R_values <- c()
#' for (e in 1:reps){
#'     chains <- gen_bvnormal_chains(m, n, rho)
#'     # R_values <- c(R_values, rhat_infinity_max_directions(chains))
#' }
#' # R_mat <- matrix(data = R_values, ncol = 1)
#' # colnames(R_mat) <- "R-hat-inf_max_dir"
#'
#' # plot_hist(R_mat, bin_size = 0.005,
#' #           lim_y_axis = 15, vaxis_pos = 1.025,
#' #           plot_threshold = F)
#'
rhat_infinity_max_directions <- function(chains,
                                         max_nb_points = 500,
                                         nb_directions = NULL) {

    dim_params = dim(chains)[2]

    # Vector that will contains the different R-hat infinity:
    directed_rhat_inf = c()

    # Vector with the ids of the different directions:
    directions_idx = 0:(2^(dim_params - 1) - 1)

    # If a number is given, sample this number of random directions:
    if (!is.null(nb_directions)) {
        directions_idx = sample(directions_idx, nb_directions, replace = F)
    }
    for (i in directions_idx) {

        directions = as.integer(intToBits(i))

        directed_rhat_inf = c(directed_rhat_inf, rhat_infinity(chains, directions, max_nb_points))

    }
    return(max(directed_rhat_inf))
}


