\name{moments}
\alias{tnorm.moment}
\alias{unif.moment}
\alias{emp.moment}
\title{Moments}

\description{
  If data (vector) is provided use \code{emp.moment()} function, otherwise for truncated normal distribution use \code{tnorm.moment()}, and for uniform distribution use \code{unif.moment()}.
}

\usage{
  tnorm.moment(mu = 0, sigma = 1, k1 = -Inf, k2 = Inf, order = 1, central = FALSE)
  unif.moment(k1 = 0, k2 = 1, order = 1, central = FALSE)
  emp.moment(x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE)
}

\arguments{
  \item{mu}{mean of truncated normal - applies to \code{tnorm.moment()}.}
  \item{sigma}{standard deviation of truncated normal - applies to \code{tnorm.moment()}.}
  \item{k1}{left truncation point for truncated normal distribution or lower bound for uniform distribution.}
  \item{k2}{right truncation point for truncated normal distribution or upper bound for uniform distribution.}
  \item{order}{+ int; order of moment}
  \item{x}{a vector of values - applies to \code{emp.moment()}.}
  \item{central}{logical; if \code{TRUE} produces central moments.}
  \item{absolute}{logical; if \code{TRUE} produces absolute moments - applies to \code{emp.moment()}.}
  \item{na.rm}{logical; if \code{TRUE} removes missing values - applies to \code{emp.moment()}.}
}

\examples{
tnorm.moment(k1 = -20, k2 = 20, order = 4, central = FALSE)
emp.moment(rnorm(10000), order = 4, central = FALSE)
unif.moment(k1 = 0, k2 = 1, order = 4, central = FALSE)
emp.moment(runif(10000), order = 4, central = FALSE)
}
