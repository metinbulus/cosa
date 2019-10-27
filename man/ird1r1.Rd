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
 Use \code{mdes.ird()} to calculate minimum detectable effect size and \code{power.ird()} to calculate statistical power. If higher level strata or fixed blocks exist, use \code{mdes.bird2f1()} to calculate minimum detectable effect size, \code{power.bird2f1()} to calculate statistical power, and \code{cosa.bird2f1()} for constrained optimal sample allocation. Alternatively modify degrees of freedom in \code{<output>.ird()} functions as \code{n1 - 2*nb - g1 - order} where \code{n1} is total number of subjects across blocks, and \code{nb} is number of blocks. Keep in mind that \code{r21} now includes information about blocks, but this fact will not be reflected in \code{g1}. See examples below.
}

\usage{
mdes.ird(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
         power = .80, alpha = .05, two.tailed = TRUE, df = n1 - g1 - order - 2,
         r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1)

power.ird(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
          es = .25, alpha = .05, two.tailed = TRUE, df = n1 - g1 - order - 2,
          r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1)

mdes.bird2f1(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             power = .80, alpha = .05, two.tailed = TRUE, df = n2 * (n1 - 2) - g1 - order,
             r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1, n2 = 1)

power.bird2f1(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
              es = .25, alpha = .05, two.tailed = TRUE, df = n2 * (n1 - 2) - g1 - order,
              r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1, n2 = 1)

cosa.bird2f1(score = NULL, order = 2, rhots = NULL,
             k1 = -6, k2 = 6, dists = "normal",
             cn1 = 0, cn2 = 0, cost = NULL,
             n1 = NULL, n2 = NULL, p = NULL,
             n0 = c(400 + g1, 5), p0 = .499,
             constrain = "power", round = TRUE, max.power = FALSE,
             local.solver = c("LBFGS", "SLSQP"),
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             g1 = 0, r21 = 0)
}

\arguments{
  \item{score}{list; an object with class 'score' returned from \code{inspect.score()} function.}
  \item{order}{integer; order of functional form for the score variable, 0 for corresponding random assignment designs, 1 for RD design with linear score variable, 2 for RD design with linear + quadratic score variable}
  \item{rhots}{correlation between the treatment and the scoring variable. Specify \code{rhots = 0} or \code{order = 0} to obtain results equivalent to random assignment designs.}
  \item{k1}{numeric; left truncation point for truncated normal dist., or lower bound for uniform dist., ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{k2}{numeric; right truncation point for truncated normal dist., or upper bound for uniform dist., ignored when \code{rhots = 0} or \code{order = 0}.}
  \item{dists}{character; distribution of the score variable, \code{"normal"} or \code{"uniform"}. By default, \code{dists = "normal"} specification implies a truncated normal distribution with \code{k1 = -6} and \code{k2 = 6}.}

  \item{power}{statistical power (1 - \eqn{\beta}).}
  \item{es}{numeric > 0; effect size (Cohen's d).}
  \item{alpha}{probability of type I error (\eqn{\alpha}).}
  \item{two.tailed}{logical; \code{TRUE} for two-tailed hypothesis testing.}
  \item{df}{degrees of freedom.}
  \item{g1}{number of covariates.}
  \item{r21}{proportion of variance in the outcome explained by covariates.}
  \item{p}{proportion of units in treatment condition.}
  \item{rate.tp}{treatment group participation rate.}
  \item{rate.cc}{control group crossover rate.}
  \item{n1}{sample size (per stratum or block, if exists).}
  \item{n2}{number of stratum or fixed blocks.}

  \item{cn1}{marginal cost per unit in treatment and control conditions, e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal cost per stratum or fixed block.}
  \item{cost}{total cost or budget.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"es"}.}
  \item{n0}{starting value for \code{n1} or \code{n1, n2}. Starting value is replaced with average when sample size is constrained by bounds.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
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
power.ird(score.obj, order = 2,
          es = 0.25, g1 = 0, r21 = 0, n = 400)

# with 5 blocks df = n1 - 2*(n blocks) - order - g1
# n1: number of subjects across five blocks
power.ird(score.obj, order = 2, df = 400 - 2*5 - 2 - 0,
          es = 0.25, g1 = 0, r21 = .30, n = 400)

# compare
# n1: number of subjects per block, n2: number of blocks
power.bird2f1(score.obj, order = 2,
              es = 0.25, g1 = 0, r21 = .30,
              n1 = 80, n2 = 5)

# optimal combination of sample sizes for subjects and blocks
# that produce power = .80 (given range restrictions)
cosa.bird2f1(score.obj, order = 2,
             constrain = "power", power = .80,
             es = 0.25, g1 = 0, r21 = .30,
             n1 = c(100, 200), n2 = c(5, 10))
}
