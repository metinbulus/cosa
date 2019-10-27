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
  Use \code{mdes.crd2()} to calculate minimum detectable effect size, \code{power.crd2()} to calculate statistical power, and \code{cosa.crd2()} for constrained optimal sample allocation. If higher level strata or fixed blocks exist, use \code{mdes.bcrd3f2()} to calculate minimum detectable effect size, \code{power.bcrd3f2()} to calculate statistical power, and \code{cosa.bcrd3f2()} for constrained optimal sample allocation. Alternatively modify degrees of freedom in \code{<output>.crd2()} functions as \code{n2 - 2*nb - g2 - order} where \code{n2} is total number of level 2 units across blocks, and \code{nb} is number of blocks. Keep in mind that \code{r22} now includes information about blocks, but this fact will not be reflected in \code{g2}. See examples below.
}

\usage{
mdes.crd2(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
          power = .80, alpha = .05, two.tailed = TRUE, df = n2 - g2 - order - 2,
          rho2, r21 = 0, r22 = 0, g2 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1, n2)

power.crd2(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
           es = .25, alpha = .05, two.tailed = TRUE, df = n2 - g2 - order - 2,
           rho2, r21 = 0, r22 = 0, g2 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1, n2)

cosa.crd2(score = NULL, order = 2, rhots = NULL,
          k1 = -6, k2 = 6, dists = "normal",
          cn1 = 0, cn2 = 0, cost = NULL,
          n1 = NULL, n2 = NULL, p = NULL,
          n0 = c(10, 100 + g2 + order), p0 = .499,
          constrain = "power", round = TRUE,
          max.power = FALSE, local.solver = c("LBFGS", "SLSQP"),
          power = .80, es = .25, alpha = .05, two.tailed = TRUE,
          rho2, g2 = 0, r21 = 0, r22 = 0)

mdes.bcrd3f2(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             power = .80, alpha = .05, two.tailed = TRUE, df = n3 * (n2 - 2) - g2 - order,
             rho2, r21 = 0, r22 = 0, g2 = 0,
             rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3)

power.bcrd3f2(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
              es = .25, alpha = .05, two.tailed = TRUE, df = n3 * (n2 - 2) - g2 - order,
              rho2, r21 = 0, r22 = 0, g2 = 0,
              rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3)

cosa.bcrd3f2(score = NULL, order = 2, rhots = NULL,
             k1 = -6, k2 = 6, dists = "normal",
             cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
             n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
             n0 = c(10, 100 + g2, 5), p0 = .499,
             constrain = "power", round = TRUE, max.power = FALSE,
             local.solver = c("LBFGS", "SLSQP"),
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             rho2, g2 = 0, r21 = 0, r22 = 0)
}

\arguments{
  \item{score}{list; an object with class 'score' returned from \code{inspect.score()} function.}
  \item{order}{integer; order of functional form for the score variable, 0 for corresponding random assignment designs, 1 for RD design with linear score variable, 2 for RD design with linear + quadratic score variable}
  \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} or \code{order = 0} to obtain results equivalent to random assignment designs.}
  \item{k1}{numeric; left truncation point for truncated normal dist., or lower bound for uniform dist., ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{k2}{numeric; right truncation point for truncated normal dist., or upper bound for uniform dist., ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{dists}{character; distribution of the score variable, \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k1 = -6} and \code{k2 = 6}.}

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
  \item{p}{proportion of level 2 units in treatment condition.}
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{number of level 2 units (per stratum or block, if exists).}
  \item{n3}{number of stratum or fixed blocks.}

  \item{cn1}{marginal costs per level 1 unit in treatment and control conditions, e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal costs per level 2 unit in treatment and control conditions, e.g. \code{c(50, 30)}.}
  \item{cn3}{marginal costs per stratum or fixed block.}
  \item{cost}{total cost or budget.}
  \item{n0}{vector of starting values for \code{n1, n2} or \code{n1, n2, n3} (positional). Starting values are replaced with averages when sample sizes are constrained by bounds.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} or \code{order = 0}, and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"es"}.}
  \item{round}{logical; \code{TRUE} for rounded COSA solution.}
  \item{max.power}{logical; \code{TRUE} for maximizing power instead of minimizing variance, applies when \code{constrain = "cost"}}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP")}}
}

\value{
  \item{parms}{list of parameters used in the function.}
  \item{df}{degrees of freedom.}
  \item{sse}{standardized standard error.}
  \item{cosa}{constrained optimal sample allocation.}
  \item{mdes}{minimum detectable effect size and (1 - \eqn{\alpha})\% confidence limits.}
  \item{power}{statistical power (1 - \eqn{\beta})}
}

\examples{
score.obj <- inspect.score(rnorm(10000), cutoff = 0)
power.crd2(score.obj, order = 2,
           es = .25, rho2 = .20, g2 = 0, r22 = 0,
           n1 = 50, n2 = 100)

# with 5 blocks df = n2 - 2*(n blocks) - order - g2
# n2: number of level 2 units across five blocks
power.crd2(score.obj, order = 2, df = 100 - 2*5 - 2 - 0,
           es = .25, rho2 = .20, g2 = 0, r22 = .30,
           n1 = 50, n2 = 100)

# compare
# n2: number of level 2 units per block, n3: number of blocks
power.bcrd3f2(score.obj, order = 2,
              es = .25, rho2 = .20, g2 = 0, r22 = .30,
              n1 = 50, n2 = 20, n3 = 5)

# optimal combination of sample sizes for level 1 and level 2
# that produce power = .80 (given range restriction for level 1 sample size)
cosa.bcrd3f2(score.obj, order = 2,
             constrain = "power", power = .80,
             es = .25, rho2 = .20, g2 = 0, r22 = .30,
             n1 = c(20, 60), n2 = NULL, n3 = 5)
}
