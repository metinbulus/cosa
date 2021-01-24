\name{inspect.score}
\alias{inspect.score}
\title{Computes Regression Discontinuity Design Effects}

\description{
  Computes Regression Discontinuity Design Effects (RDDE) either based on analytic deviations (up to second order with interactions), an empirical score variable, or simulation.
}

\usage{
  inspect.score(score = NULL, p = NULL, cutoff = NULL,
              treat.lower = FALSE, order = 1, interaction = FALSE,
              mu = 0, sigma = 1, k1 = -Inf, k2 =  Inf,
              dists = "normal", sim = FALSE, ndraw = 1000, nsim = 1000)
}

\arguments{
  \item{sim}{logical; if \code{TRUE} results are based on simulation.}
  \item{score}{vector; score variable.}
  \item{p}{proportion of units in the treatment condition.}
  \item{cutoff}{decision threshold.}
  \item{treat.lower}{logical; if \code{TRUE} units below cutoff are treated.}
  \item{order}{integer >= 0; order of polynomial functional form specification for the score variable.}
  \item{interaction}{logical; if \code{TRUE} polynomial specification interacts with the treatment variable.}
  \item{mu}{mean of (uncentered) truncated normal - applies when \code{score = NULL} and \code{dists = "normal"}.}
  \item{sigma}{standard deviation of (uncentered) truncated normal - applies when \code{score = NULL} and \code{dists = "normal"}.}
  \item{k1}{left truncation point for (uncentered) empirical, truncated normal, or uniform distribution.}
  \item{k2}{right truncation point for (uncentered) empirical, truncated normal, or uniform distribution.}
  \item{dists}{char; type of distribution, \code{"normal"} or \code{"uniform"}.}
  \item{ndraw}{number of draws - applies when \code{sim = TRUE}.}
  \item{nsim}{number of simulations - applies when \code{sim = TRUE}.}
}

\value{
  \item{parms}{list; list of parameters used in the computation.}
  \item{cutoff}{decision threshold (computed if \code{p} is provided).}
  \item{treat.lower}{if \code{TRUE} units below cutoff are treated.}
  \item{p}{proportion of subjects treated (computed if \code{cutoff} is provided).}
  \item{order}{order of polynomial specification for the score variable.}
  \item{interaction}{if \code{TRUE} polynomial specification interacts with the treatment variable.}
  \item{center}{if \code{TRUE} the score variable is centered on the cutoff}
  \item{rdde}{regression discontinuity design effect.}
}

\examples{
  # based on an empirical score variable
  inspect.score(score = rnorm(10000), cutoff = 0)
  # based on analytic derivation
  inspect.score(cutoff = 0)
  # based on simulation
  inspect.score(sim = TRUE, cutoff = 0)
}
