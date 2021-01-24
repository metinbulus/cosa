\name{crd2}
\alias{crd2r2}
\alias{crd2}
\alias{bcrd3f2}
\alias{cosa.crd2r2}
\alias{cosa.crd2}
\alias{cosa.bcrd3f2}
\alias{power.crd2r2}
\alias{power.crd2}
\alias{power.bcrd3f2}
\alias{mdes.crd2r2}
\alias{mdes.crd2}
\alias{mdes.bcrd3f2}
\title{Cluster-level Regression Discontinuity (Two-level Design, Discontinuity at Level 2, w/ or w/o Strata or Fixed Blocks)}

\description{
  Use \code{mdes.crd2()} to calculate minimum detectable effect size, \code{power.crd2()} to calculate statistical power, and \code{cosa.crd2()} for bound constrained optimal sample size allocation (BCOSSA). If higher level strata or fixed blocks exist, use \code{mdes.bcrd3f2()} to calculate minimum detectable effect size, \code{power.bcrd3f2()} to calculate statistical power, and \code{cosa.bcrd3f2()} for BCOSSA.
}

\usage{
mdes.crd2(score = NULL, dists = "normal", k1 = -6, k2 = 6,
          order = 1, interaction = FALSE,
          treat.lower = TRUE, cutoff = 0, p = NULL,
          power = .80, alpha = .05, two.tailed = TRUE,
          df = n2 - g2 - order * (1 + interaction) - 2,
          rho2, r21 = 0, r22 = 0, g2 = 0, rate.tp = 1, rate.cc = 0, n1, n2)

power.crd2(score = NULL, dists = "normal", k1 = -6, k2 = 6,
           order = 1, interaction = FALSE,
           treat.lower = TRUE, cutoff = 0, p = NULL,
           es = .25, alpha = .05, two.tailed = TRUE,
           df = n2 - g2 - order * (1 + interaction) - 2,
           rho2, r21 = 0, r22 = 0, g2 = 0, rate.tp = 1, rate.cc = 0, n1, n2)

cosa.crd2(score = NULL, dists = "normal", k1 = -6, k2 = 6, rhots = NULL,
          order = 1, interaction = FALSE,
          treat.lower = TRUE, cutoff = 0, p = NULL,
          cn1 = 0, cn2 = 0, cost = NULL,
          n1 = NULL, n2 = NULL, n0 = c(10, 100), p0 = .499,
          constrain = "power", round = TRUE,
          max.power = FALSE, local.solver = c("LBFGS", "SLSQP"),
          power = .80, es = .25, alpha = .05, two.tailed = TRUE,
          rho2, g2 = 0, r21 = 0, r22 = 0)

mdes.bcrd3f2(score = NULL, dists = "normal", k1 = -6, k2 = 6,
             order = 1, interaction = FALSE,
             treat.lower = TRUE, cutoff = 0, p = NULL,
             power = .80, alpha = .05, two.tailed = TRUE,
             df = n3 * (n2 - 2) - g2 - order * (1 + interaction),
             rho2, r21 = 0, r22 = 0, g2 = 0,
             rate.tp = 1, rate.cc = 0, n1, n2, n3)

power.bcrd3f2(score = NULL, dists = "normal", k1 = -6, k2 = 6,
              order = 1, interaction = FALSE,
              treat.lower = TRUE, cutoff = 0, p = NULL,
              es = .25, alpha = .05, two.tailed = TRUE,
              df = n3 * (n2 - 2) - g2 - order * (1 + interaction),
              rho2, r21 = 0, r22 = 0, g2 = 0,
              rate.tp = 1, rate.cc = 0, n1, n2, n3)

cosa.bcrd3f2(score = NULL, dists = "normal", k1 = -6, k2 = 6, rhots = NULL,
             order = 1, interaction = FALSE,
             treat.lower = TRUE, cutoff = 0, p = NULL,
             cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
             n1 = NULL, n2 = NULL, n3 = NULL,
             n0 = c(10, 100, 5), p0 = .499,
             constrain = "power", round = TRUE, max.power = FALSE,
             local.solver = c("LBFGS", "SLSQP"),
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             rho2, g2 = 0, r21 = 0, r22 = 0)
}

\arguments{
  \item{score}{vector or list; an empirical score variable or an object with class 'score' returned from the \code{inspect.score()} function.}
  \item{dists}{character; distribution of the score variable, \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k1 = -6} and \code{k2 = 6}.}
  \item{k1}{left truncation point for (uncentered) empirical, truncated normal, or uniform distribution. Ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{k2}{right truncation point for (uncentered) empirical, truncated normal, or uniform distribution. Ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{order}{integer >= 0; order of polynomial functional form specification for the score variable.}
  \item{interaction}{logical; if \code{TRUE} polynomial specification interacts with the treatment variable.}
  \item{rhots}{obsolote; use \code{order = 0} to obtain results equivalent to random assignment designs.}
  \item{treat.lower}{logical; if \code{TRUE} units below the cutoff are treated.}
  \item{cutoff}{decision threshold.}
  \item{p}{proportion of level 2 units in the treatment condition.}

  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{df}{degrees of freedom.}
  \item{rho2}{proportion of variance in the outcome between level 2 units (unconditional ICC2).}
  \item{g2}{number of covariates at level 2.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{rate.tp}{treatment group participation rate.}
  \item{rate.cc}{control group crossover rate.}
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{number of level 2 units (per stratum or block, if exists).}
  \item{n3}{number of stratum or fixed blocks.}

  \item{cn1}{marginal costs per level 1 unit in treatment and control conditions (positional), e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal costs per level 2 unit in treatment and control conditions (positional), e.g. \code{c(50, 30)}.}
  \item{cn3}{marginal cost per stratum or fixed block.}
  \item{cost}{total cost or budget. Ignored when \code{constrain = "power"} or \code{constrain = "es"}. }
  \item{n0}{vector of starting values for \code{n1, n2} or \code{n1, n2, n3} (positional). Starting values are replaced with the averages when sample sizes are constrained by bounds.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} or \code{order = 0}, and \code{p = NULL}. Starting value is replaced with the average when \code{p} is constrained by bounds.}
  \item{constrain}{character; constrains one of the \code{"cost"}, \code{"power"}, or \code{"es"} at the specified value.}
  \item{round}{logical; \code{TRUE} for rounded BCOSSA solution.}
  \item{max.power}{logical; \code{TRUE} for maximizing the power rate instead of minimizing the variance. Applies when \code{constrain = "cost"}.}
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
power.crd2(score.obj,
           es = .25, rho2 = .20, g2 = 0, r22 = 0,
           n1 = 50, n2 = 30)

# with 5 blocks (note that r22 is modified but g2 remains the same)
power.bcrd3f2(score.obj,
              es = .25, rho2 = .20, g2 = 0, r22 = .30,
              n1 = 50, n2 = 30, n3 = 5)

# minimum required number of level 2 units for each block
cosa.bcrd3f2(score.obj,
             es = .25, rho2 = .20, g2 = 0, r22 = .30,
             n1 = 50, n2 = NULL, n3 = 5)

}
