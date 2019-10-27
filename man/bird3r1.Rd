\name{bird3}
\alias{bird3r1}
\alias{bird3}
\alias{cosa.bird3r1}
\alias{cosa.bird3}
\alias{mdes.bird3r1}
\alias{mdes.bird3}
\alias{power.bird3r1}
\alias{power.bird3}
\title{Blocked Individual-level Regression Discontinuity (Three-level Design, Discontinuity at Level 1)}

\description{
 Use \code{mdes.bird3()} to calculate minimum detectable effect size, \code{power.bird3()} to calculate statistical power, and \code{cosa.bird3()} for constrained optimal sample allocation. To consider fixed block effects, modify degrees of freedom in \code{<output>.bird3()} functions as \code{n3 - 2*nb - g3} where \code{n3} is total number of level 3 units across blocks, and \code{nb} is number of blocks. Keep in mind that \code{r2t3} now includes information about blocks, but this fact will not be reflected in \code{g3}. See examples below.
}

\usage{
mdes.bird3(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
           power = .80, alpha = .05, two.tailed = TRUE, df = n3 - g3 - 1,
           rho2, rho3, omega2, omega3, r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
           rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3)

power.bird3(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            es = .25, alpha = .05, two.tailed = TRUE, df = n3 - g3 - 1,
            rho2, rho3, omega2, omega3, r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
            rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3)

cosa.bird3(score = NULL, order = 2, rhots = NULL,
           k1 = -6, k2 = 6, dists = "normal",
           cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
           n1 = NULL, n2 = NULL, n3 = NULL, p = NULL,
           n0 = c(10, 3, 100 + g3), p0 = .499,
           constrain = "power", round = TRUE, max.power = FALSE,
           local.solver = c("LBFGS", "SLSQP"),
           power = .80, es = .25, alpha = .05, two.tailed = TRUE,
           rho2, rho3, omega2, omega3,
           g3 = 0, r21 = 0, r2t2 = 0, r2t3 = 0)
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
  \item{omega2}{ratio of the treatment effect variance between level 2 units to the variance in the outcome between level 2 units.}
  \item{omega3}{ratio of the treatment effect variance between level 3 units to the variance in the outcome between level 3 units.}
  \item{g3}{number of covariates at level 3.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r2t2}{proportion of treatment effect variance between level 2 units explained by level 2 covariates.}
  \item{r2t3}{proportion of treatment effect variance between level 3 units explained by level 3 covariates.}
  \item{rate.tp}{treatment group participation rate.}
  \item{rate.cc}{control group crossover rate.}
  \item{p}{proportion of level 1 units in treatment condition.}
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{average number of level 2 units per level 3 unit.}
  \item{n3}{number of level 3 units.}


  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions, e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal cost per level 2 unit.}
  \item{cn3}{marginal cost per level 3 unit.}
  \item{cost}{total cost or budget.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{n0}{vector of starting values for \code{n1, n2, n3} (positional). Starting values are replaced with averages when sample sizes are constrained by bounds.}
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
power.bird3(score.obj, order = 2,
            es = 0.25, rho2 = .20, rho3 = .10,
            omega2 = .30, omega3 = .30,
            g3 = 0, r2t3 = 0, n1 = 50, n2 = 3, n3 = 15)


# with 5 blocks df = n3- 2*(n blocks) - g3
# n3: number of level 3 units across five blocks
# increase in r2t3 does not make up for reduction in df
power.bird3(score.obj, order = 2, df = 15 - 2*5 - 0,
            es = 0.25, rho2 = .20, rho3 = .10,
            omega2 = .30, omega3 = .30,
            g3 = 0, r2t3 = .30, n1 = 50, n2 = 3, n3 = 15)


# optimal combination of sample sizes for level 1, level 2 and level 3
# that produce power = .80 (given range restrictions)
cosa.bird3(score.obj, order = 2,
           constrain = "power", power = .80,
           es = 0.25, rho2 = .20, rho3 = .10,
           omega2 = .30, omega3 = .30,
           g3 = 0, r2t3 = 0,
           n1 = c(15,30), n2 = c(3, 5), n3 = c(10,30))
}
