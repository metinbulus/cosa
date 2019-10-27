\name{inspect.score}
\alias{inspect.score}
\title{Inspects Relations between Treatment, Score and Score^2 Triad}

\description{
  Inpects relations between Treatment, Score and Score^2 triad, outputs correlations and design effects for linear and linear + quadratic functional forms for the score variable.
}

\usage{
  inspect.score(score = NULL, sim = FALSE, p = NULL, cutoff = NULL, treat.lower = FALSE,
                          mu = 0, sigma = 1, k1 = -1e+10, k2 =  1e+10,
                          dists = "normal", ndraw = 1000, nsim = 1000)
}

\arguments{
  \item{sim}{logical; if \code{TRUE} results are based on simulation.}
  \item{score}{vector; score variable.}
  \item{p}{proportion of units in treatment condition.}
  \item{cutoff}{cutoff.}
  \item{treat.lower}{logical; if \code{TRUE} subjects below cutoff are treated.}
  \item{mu}{mean of truncated normal - applies when \code{score = NULL} and \code{dists = "normal"}.}
  \item{sigma}{standard deviation of truncated normal - applies when \code{score = NULL} and \code{dists = "normal"}.}
  \item{k1}{left truncation point for empirical, truncated normal, or uniform distribution.}
  \item{k2}{right truncation point for empirical, truncated normal, or uniform distribution.}
  \item{dists}{char; type of distribution, "normal" or "uniform".}
  \item{ndraw}{number of draws - applies when \code{sim = TRUE}.}
  \item{nsim}{number of simulations - applies when \code{sim = TRUE}.}
}

\value{
  \item{parms}{list; list of parameters used in the computation.}
  \item{cutoff}{cutoff score (computed if \code{p} is provided).}
  \item{treat.lower}{logical; not used very much but may be utilized in the future.}
  \item{p}{proportion of subjects treated (computed if \code{cutoff} is provided).}
  \item{rhots}{correlation between Treatment and Score.}
  \item{rhots2}{correlation between Treatment and Score^2.}
  \item{rhoss2}{{correlation between Score and Score^2}.}
  \item{d1}{design effect for linear functional form.}
  \item{d2}{design effect for linear + quadratic functional form.}
}

\examples{
  inspect.score(score = rnorm(10000), p = .50)

  # default based on ~ N(0,1)
  inspect.score(p = .50)
  inspect.score(sim = TRUE, p = .50)
}
