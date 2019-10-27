\name{crd3}
\alias{crd3r3}
\alias{crd3}
\alias{bcrd4f3}
\alias{cosa.crd3r3}
\alias{cosa.crd3}
\alias{cosa.bcrd4f3}
\alias{power.crd3r3}
\alias{power.crd3}
\alias{power.bcrd4f3}
\alias{mdes.crd3r3}
\alias{mdes.crd3}
\alias{mdes.bcrd4f3}
\title{Cluster-level Regression Discontinuity (Three-level Design, Discontinuity at Level 3, w/ or w/o Strata or Fixed Blocks)}

\description{
  Use \code{mdes.crd3()} to calculate minimum detectable effect size, \code{power.crd3()} to calculate statistical power, and \code{cosa.crd3()} for constrained optimal sample allocation. If higher level strata or fixed blocks exist, use \code{mdes.bcrd4f3()} to calculate minimum detectable effect size, \code{power.bcrd4f3()} to calculate statistical power, and \code{cosa.bcrd4f3()} for constrained optimal sample allocation. Alternatively modify degrees of freedom in \code{<output>.crd3()} functions as \code{n3 - 2*nb - g3 - order} where \code{n3} is total number of level 3 units across blocks, and \code{nb} is number of blocks. Keep in mind that \code{r23} now includes information about blocks, but this fact will not be reflected in \code{g3}. See examples below.
}

\usage{
mdes.crd3(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
          power = .80, alpha = .05, two.tailed = TRUE, df = n3 - g3 - order - 2,
          rho2, rho3, r21 = 0, r22 = 0, r23 = 0,
          g3 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3)

power.crd3(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
           es = .25, alpha = .05, two.tailed = TRUE, df = n3 - g3 - order - 2,
           rho2, rho3, r21 = 0, r22 = 0, r23 = 0,
           g3 = 0, rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3)

cosa.crd3r3(score = NULL, order = 2, rhots = NULL,
            k1 = -6, k2 = 6, dists = "normal",
            cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
            n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
            n0 = c(10, 3, 100 + g3 + order), p0 = .499,
            constrain = "power", round = TRUE, max.power = FALSE,
            local.solver = c("LBFGS", "SLSQP"),
            power = .80, es = .25, alpha = .05, two.tailed = TRUE,
            rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0)

mdes.bcrd4f3(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
             power = .80, alpha = .05, two.tailed = TRUE, df = n4 * (n3 - 2) - g3 - order,
             rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
             rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3, n4)

power.bcrd4f3(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
              es = .25, alpha = .05, two.tailed = TRUE, df = n4 * (n3 - 2) - g3 - order,
              rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
              rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3, n4)

cosa.bcrd4f3(score = NULL, order = 2, rhots = NULL,
             k1 = -6, k2 = 6, dists = "normal",
             cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
             n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL,
             p = NULL, n0 = c(10, 3, 100 + g3, 5), p0 = .499,
             constrain = "power", round = TRUE, max.power = FALSE,
             local.solver = c("LBFGS", "SLSQP"),
             power = .80, es = .25, alpha = .05, two.tailed = TRUE,
             rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0)
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
  \item{rho3}{proportion of variance in the outcome between level 3 units (unconditional ICC3).}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r22}{proportion of level 2 variance in the outcome explained by level 2 covariates.}
  \item{r23}{proportion of level 3 variance in the outcome explained by level 3 covariates.}
  \item{rate.tp}{treatment group participation rate.}
  \item{rate.cc}{control group crossover rate.}
  \item{p}{proportion of level 3 units in treatment condition.}
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{average number of level 2 units per level 3 unit.}
  \item{n3}{number of level 3 units(per stratum or block, if exists).}
  \item{n4}{number of stratum or fixed blocks.}

  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions, e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal cost per level 2 unit in treatment and control conditions, e.g. \code{c(50, 30)}.}
  \item{cn3}{marginal cost per level 3 unit in treatment and control conditions, e.g. \code{c(80, 50)}.}
  \item{cn4}{marginal cost per stratum or fixed block.}
  \item{cost}{total cost or budget.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{n0}{vector of starting values for \code{n1, n2, n3} or \code{n1, n2, n3, n4} (positional). Starting values are replaced with averages when sample sizes are constrained by bounds.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{round}{logical; \code{TRUE} for rounded COSA solution.}
  \item{max.power}{logical; \code{TRUE} for maximizing power instead of minimizing variance.}
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
power.crd3(score.obj, order = 2,
           es = .25, rho2 = .20, rho3 = .10,
           g3 = 0, r23 = 0, n1 = 20, n2 = 3, n3 = 100)

# with 5 blocks df = n3 - 2*(n blocks) - order - g3
# n3: number of level 3 units across five blocks
power.crd3(score.obj, order = 2, df = 100 - 2*5 - 2 - 0,
           es = .25, rho2 = .20, rho3 = .10,
           g3 = 0, r23 = .30, n1 = 20, n2 = 3, n3 = 100)

# compare
# n3: number of level 3 units per block, n4: number of blocks
power.bcrd4f3(score.obj, order = 2,
              es = .25, rho2 = .20, rho3 = .10,
              g3 = 0, r23 = .30,
              n1 = 20, n2 = 3, n3 = 20, n4 = 5)

# optimal combination of sample sizes for level 1 and level 3
# that produce power = .80 (given range restriction for level 1 sample size)
cosa.bcrd4f3(score.obj, order = 2,
             constrain = "power", power = .80,
             es = .25, rho2 = .20, rho3 = .10, g3 = 0, r23 = .30,
             n1 = c(20, 60), n2 = 2, n3 = NULL, n4 = 5)
}
