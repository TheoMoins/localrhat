############ MULTIVARIATE VERSION ############

#' Multivariate version of the local-\eqn{\hat{R}}: \eqn{\hat{R}(x)}
#'
#' Compute \eqn{\hat{R}(x)} in the multivariate case, which correspond to \eqn{\hat{R}}
#' computed on a multivariate indicator variable for a given quantile x.
#'
#' @param x a vector of size \eqn{d} corresponding to the quantile used for the computation
#' of \eqn{\hat{R}(x)}.
#' @param chains an array of size \eqn{n \times d \times m}, where \eqn{n} is
#' the length of the chains, \eqn{m \geq 2} is the number of chains and \eqn{d}
#' is the dimension.
#'
#' @return The multivariate \eqn{\hat{R}(x)} computed on the \eqn{m} chains.
#'
#' @details The function return \eqn{\hat{R}} computed on
#' \deqn{I\{\theta_1^{(j)} \leq x_1, \ldots, \theta_d^{(j)} \leq x_d\}},
#' for a given value of \eqn{x \in \mathbb{R}^d}.
#'
#' @references TO INCLUDE
#'
#' @examples
#' library("localrhat")
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
#' chains <- gen_bvnormal_chains(m, n, rho)
#'
#' # Quantile to evaluate:
#' x <- c(0.5, 0.5)
#' # multivariate_local_rhat(x, chains)
#'
multivariate_local_rhat <- function(x, chains) {
    dim_ch <- dim(chains)
    bool_ch <- rep(NA, dim_ch[3])
    for (ch_id in 1:dim_ch[3]) {
        bool_ch[ch_id] <- apply(t(t(chains[, , ch_id]) <= x), 1, all)
    }
    return (trad_rhat(array(bool_ch, c(dim_ch[1], dim_ch[3]))))
}

#' Quantile values used for the computation of  multivariate \eqn{\hat{R}_\infty}
#'
#' Return the set of points that will be used to estimate the supremum of
#' \eqn{\hat{R}(x)} over \eqn{x \in \mathbb{R}^d}. The function works in the
#' same way as \code{grid_for_R} in the univariate case, see the corresponding
#' documentation for more details.
#'
#' @param chains an array of size \eqn{n \times d \times m}, where \eqn{n} is
#' the length of the chains, \eqn{m \geq 2} is the number of chains and \eqn{d}
#' is the dimension.
#' @param max_nb_points the maximal length of the grid in the case where the total
#' number of samples is larger. By default, \code{max_nb_points = 500}.
#'
#' @return A list that contains the different \eqn{x} used for the computation of
#' \eqn{\hat{R}(x)}.
#'
#'
#' @references TO INCLUDE
#'
multivariate_grid_for_R <- function(chains, max_nb_points = 500) {
    grid = chains[, , 1]
    for (ch in 1:dim(chains)[3]) {
        grid = rbind(grid, chains[, , ch])
    }
    if (dim(grid)[1] > max_nb_points & max_nb_points != "ALL") {
        ind = as.integer(seq(1, dim(grid)[1], length.out = max_nb_points))
        grid = grid[ind, ]
    }
    return(grid)
}

#' Multivariate \eqn{\hat{R}(x)} with a specified dependence direction
#'
#' Compute \eqn{\hat{R}(x)} in the multivariate case, for a given quantile x and
#' a given sense for the \eqn{d} signs in the indicator variable.
#'
#' @param x a vector of size \eqn{d} corresponding to the quantile used for the computation
#' of \eqn{\hat{R}(x)}.
#' @param chains an array of size \eqn{n \times d \times m}, where \eqn{n} is
#' the length of the chains, \eqn{m \geq 2} is the number of chains and \eqn{d}
#' is the dimension.
#' @param dir a binary vector of size \eqn{d} indicating the signs in the indicator
#' variable. For example, \code{dir = c(0, ..., 0)} correspond to the computation
#' of \eqn{\hat{R}} on \eqn{I\{\theta_1^{(j)} \leq x_1, \ldots, \theta_d^{(j)} \leq x_d\}},
#' and if \code{dir[i] == 1}, then the indicator with \eqn{\theta_d^{(j)} \geq x_d}
#' will be used.
#'
#' @return The corresponding multivariate \eqn{\hat{R}(x)} computed on the
#' given dependence direction.
#'
#' @details In the multivariate case, \eqn{\hat{R}} can be computed on
#' \eqn{I\{\theta_1^{(j)} \leq x_1, \ldots, \theta_d^{(j)} \leq x_d\}},
#' but also with indicator where some "\eqn{\leq}" are replaced by "\eqn{\geq}",
#' which leads to \eqn{2^{d-1}} possibilities and potentially as many different
#' results of \eqn{\hat{R}_\infty}. Here, a 0-1 vector has to be mentioned in
#' the function to specify which indicator to use, indicating which size on each
#' dimension.
#'
#' @references TO INCLUDE
#'
#' @examples
#' library("localrhat")
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
#' chains <- gen_bvnormal_chains(m, n, rho)
#'
#' # Quantile to evaluate:
#' x <- c(0.5, 0.5)
#' # multivariate_local_rhat(x, chains)
#' # multivariate_directed_local_rhat(x, chains)
#'
multivariate_directed_local_rhat <- function(x, chains, dir) {
    dim_ch <- dim(chains)
    bool_ch = rep(NA, dim_ch[3])
    for (ch_id in 1:dim_ch[3]) {
        tmp_ind = rep(NA, dim_ch[2])
        for (dim_id in 1:dim_ch[2]) {
            if (dir[dim_id] == 1) {
                tmp_ind[dim_id] = chains[, dim_id, ch_id] >= x[dim_id]
            } else {
                tmp_ind[dim_id] = chains[, dim_id, ch_id] <= x[dim_id]
            }
        }
        tmp_ch = array(tmp_ind, c(dim_ch[1], dim_ch[2]))

        bool_ch[ch_id] = apply(tmp_ch, 1, all)
    }
    return (trad_rhat(array(bool_ch, c(dim_ch[1], dim_ch[3]))))
}

#' Computation of multivariate \eqn{\hat{R}(x)} on a given grid
#'
#' Return the vector of \eqn{\hat{R}(x)} values that will be used to compute
#' the supremum \eqn{\hat{R}_\infty} in the multivariate case.
#' See \code{multivariate_grid_for_R} for the grid computation for \eqn{x}.
#'
#' @param chains an array of size \eqn{n \times d \times m}, where \eqn{n} is
#' the length of the chains, \eqn{m \geq 2} is the number of chains and \eqn{d}
#' is the dimension.
#' @param max_nb_points the maximal length of the grid in the case where the total
#' number of samples is larger. By default, \code{max_nb_points = 500}.
#' @param dir a binary vector of size \eqn{d} indicating the signs in the indicator
#' variable. For example, \code{dir = c(0, ..., 0)} correspond to the computation
#' of \eqn{\hat{R}} on \eqn{I\{\theta_1^{(j)} \leq x_1, \ldots, \theta_d^{(j)} \leq x_d\}}.
#' This direction will be used if no arguments is given.
#'
#' @return A list that contains the different \eqn{\hat{R}(x)} values used for
#' the computation of \eqn{\hat{R}_\infty}.
#'
#' @details The computation of \eqn{\hat{R}_\infty} require to compute a supremum
#' of \eqn{\hat{R}(x)} values over the quantiles \eqn{x}. See the documentation
#' of \code{multivariate_grid_for_R} and \code{grid_for_R} to see how the choice
#' of the different \eqn{x} is done.
#' Thus, this function combine \code{multivariate_grid_for_R} for the computation
#' of the grid and \code{multivariate_directed_local_rhat} or
#' \code{multivariate_directed_local_rhat}, depending on if a direction is specified.
#'
#' @references TO INCLUDE
#'
#'
multivariate_all_local_rhat <- function(chains, dir = NULL, max_nb_points = 500) {
    grid = multivariate_grid_for_R(chains, max_nb_points)
    tab_rhat = rep(NA, dim(grid)[1])
    for (i in 1:dim(grid)[1]) {
        if (is.null(dir)) {
            tab_rhat[i] = multivariate_local_rhat(grid[i, ], chains)
        } else {
            tab_rhat[i] = multivariate_directed_local_rhat(grid[i, ], chains, dir)
        }
    }
    return(tab_rhat)
}
