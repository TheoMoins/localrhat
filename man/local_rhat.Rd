% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_rhat.R
\name{local_rhat}
\alias{local_rhat}
\title{Local version of the Gelman Rubin diagnostic \eqn{\hat{R}}}
\usage{
local_rhat(x, chains)
}
\arguments{
\item{x}{a float number corresponding to the quantile used for the computation
of \eqn{\hat{R}(x)}.}

\item{chains}{an array of size \eqn{n \times m} where \eqn{n} is the length of
the chains and \eqn{m \geq 2} is the number of chains.}
}
\value{
The local-\eqn{\hat{R}(x)} computed on the \eqn{m} chains.
}
\description{
Compute \eqn{\hat{R}(x)}, a version of univariate \eqn{\hat{R}} computed
on indicator variables for a given quantile x.
}
\details{
The function return \eqn{\hat{R}} computed on \eqn{I(\theta^{(i,j)} \leq x)}
for a given value of \eqn{x}:
\deqn{\hat{R}(x) = \sqrt{\frac{\frac{n-1}{n} \hat{W}(x) + \hat{B}(x)}{\hat{W}(x)}},}
with \eqn{\hat{W}(x)} and \eqn{\hat{B}(x)} the estimated local within-chain variance
and between-chain variance:
\deqn{
\hat{W}(x) = \frac{1}{m} \sum_{j=1}^m (F_j(x) - F^2_j(x))
}
And
\deqn{
\hat{B}(x) = \frac{1}{m^2}\sum_{j<k} \left(F_j(x)-F_{k}(x)\right)^2.
}
}
\examples{
library(localrhat)

N <- 500 # length of chains
M <- 4 # number of chains

# Toy example with 3 i.i.d chains uniform in [-0.5, 0.5] and 1 in [-1, 1]:
chains <- array(c(runif((M-1)*N, -1/2, 1/2), runif(N, -1,1)), c(N,M))

# Quantile to evaluate:
x <- 0.5
# local_rhat(x, chains)

}
\references{
TO INCLUDE
}
