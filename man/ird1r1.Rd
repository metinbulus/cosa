\name{ird}
\alias{ird1r1}
\alias{ird}
\alias{bird2f1}
\alias{cosa.bird2f1}
\alias{mdes.ird1r1}
\alias{mdes.ird}
\alias{mdes.bird2f1}
\alias{power.ird1r1}
\alias{power.ird}
\alias{power.bird2f1}
\title{Simple Individual-level Regression Discontinuity (w/ or w/o Strata or Fixed Blocks)}

\description{
 Use \code{mdes.ird()} to calculate minimum detectable effect size and \code{power.ird()} to calculate statistical power. If higher level strata or fixed blocks exist, use \code{mdes.bird2f1()} to calculate minimum detectable effect size, \code{power.bird2f1()} to calculate statistical power, and \code{cosa.bird2f1()} for bound constrained optimal sample size allocation (BCOSSA).
}

\usage{
mdes.ird(score = NULL, dists = "normal", k1 = -6, k2 = 6,
         order = 1, interaction = FALSE,
         treat.lower = TRUE, cutoff = 0, p = NULL,
         power = .80, alpha = .05, two.tailed = TRUE,
         df = n1 - g1 - order * (1 + interaction) - 2,
         r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, n1)

power.ird(score = NULL, dists = "normal", k1 = -6, k2 = 6,
          order = 1, interaction = FALSE,
          treat.lower = TRUE, cutoff = 0, p = NULL,
          es = .25, alpha = .05, two.tailed = TRUE,
          df = n1 - g1 - order * (1 + interaction) - 2,
          r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, n1)

mdes.bird2f1(score = NULL, dists = "normal", k1 = -6, k2 = 6,
             order = 1, interaction = FALSE,
             treat.lower = TRUE, cutoff = 0, p = NULL,
             power = .80, alpha = .05, two.tailed = TRUE,
             df = n2 * (n1 - 2) - g1 - order * (1 + interaction),
             r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, n1, n2 = 1)

power.bird2f1(score = NULL, dists = "normal", k1 = -6, k2 = 6,
              order = 1, interaction = FALSE,
              treat.lower = TRUE, cutoff = 0, p = NULL,
              es = .25, alpha = .05, two.tailed = TRUE,
              df = n2 * (n1 - 2) - g1 - order * (1 + interaction),
              r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, n1, n2 = 1)

cosa.bird2f1(score = NULL, dists = "normal", k1 = -6, k2 = 6, rhots = NULL,
             order = 1, interaction = FALSE,
             treat.lower = TRUE, cutoff = 0, p = NULL,
             cn1 = 0, cn2 = 0, cost = NULL,
             n1 = NULL, n2 = NULL,
             n0 = c(400, 5), p0 = .499,
             constrain = "power", round = TRUE, max.power = FALSE,
             local.solver = c("LBFGS", "SLSQP"),
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             g1 = 0, r21 = 0)
}

\arguments{
  \item{score}{vector or list; an empirical score variable or an object with class 'score' returned from the \code{inspect.score()} function.}
  \item{dists}{character; distribution of the score variable, \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k1 = -6} and \code{k2 = 6}.}
  \item{k1}{left truncation point for (uncentered) empirical, truncated normal, or uniform distribution. Ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{k2}{right truncation point for (uncentered) empirical, truncated normal, or uniform distribution. Ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{order}{integer >= 0; order of polynomial functional form specification for the score variable.}
  \item{interaction}{logical; if \code{TRUE} polynomial specification interacts with the treatment variable.}
  \item{rhots}{obsolote; use \code{order = 0} to obtain results equivalent to random assignment designs.}
  \item{treat.lower}{logical; if \code{TRUE} units below cutoff the are treated.}
  \item{cutoff}{decision threshold.}
  \item{p}{proportion of units in the treatment condition.}

  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{numeric > 0; effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{df}{degrees of freedom.}
  \item{g1}{number of covariates.}
  \item{r21}{proportion of variance in the outcome explained by covariates.}
  \item{rate.tp}{treatment group participation rate.}
  \item{rate.cc}{control group crossover rate.}
  \item{n1}{sample size (per stratum or block, if exists).}
  \item{n2}{number of stratum or fixed blocks.}

  \item{cn1}{marginal cost per unit in treatment and control conditions, e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal cost per stratum or fixed block.}
  \item{cost}{total cost or budget. Ignored when \code{constrain = "power"} or \code{constrain = "es"}. }
  \item{constrain}{character; constrains one of the \code{"cost"}, \code{"power"}, or \code{"es"} at the specified value.}
  \item{n0}{starting value for \code{n1} or \code{n1, n2}. Starting value is replaced with the average when sample size is constrained by bounds.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{round}{logical; \code{TRUE} for rounded BCOSSA solution.}
  \item{max.power}{logical; \code{TRUE} for maximizing power instead of minimizing variance, applies when \code{constrain = "cost"}}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP")}}
}

\value{
  \item{parms}{list of parameters used in the function.}
  \item{df}{degrees of freedom.}
  \item{sse}{standardized standard error.}
  \item{cosa}{BCOSSA solution.}
  \item{mdes}{minimum detectable effect size and (1 - \eqn{\alpha})\% confidence limits.}
  \item{power}{statistical power (1 - \eqn{\beta})}
}

\examples{
score.obj <- inspect.score(rnorm(1000),
                           order = 1, interaction = FALSE,
                           cutoff = 0, k1 = -1, k2 = 1)

# single site (no blocks)
power.ird(score.obj, g1 = 0, r21 = 0,
          es = 0.25, n = 100)

# with 5 blocks (note that r21 is modified but g1 remains the same)
power.bird2f1(score.obj,  g1 = 0, r21 = .30,
              es = 0.25, n1 = 100, n2 = 5)

# minimum required sample size for each block
cosa.bird2f1(score.obj,  g1 = 0, r21 = .30,
             n1 = NULL, n2 = 5)
}
