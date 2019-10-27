\name{bird4}
\alias{bird4r1}
\alias{bird4}
\alias{cosa.bird4r1}
\alias{cosa.bird4}
\alias{mdes.bird4r1}
\alias{mdes.bird4}
\alias{power.bird4r1}
\alias{power.bird4}
\title{Blocked Individual-level Regression Discontinuity (Four-level Design, Discontinuity at Level 1)}

\description{
 Use \code{mdes.bird4()} to calculate minimum detectable effect size, \code{power.bird4()} to calculate statistical power, and \code{cosa.bird4()} for constrained optimal sample allocation.
}

\usage{
mdes.bird4(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
           power = .80, alpha = .05, two.tailed = TRUE, df = n4 - g4 - 1,
           rho2, rho3, rho4, omega2, omega3, omega4,
           r21 = 0, r2t2 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
           rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3, n4)

power.bird4(score = NULL, order = 2, rhots = NULL, k1 = -6, k2 = 6, dists = "normal",
            es = .25, alpha = .05, two.tailed = TRUE, df = n4 - g4 - 1,
            rho2, rho3, rho4, omega2, omega3, omega4,
            r21 = 0, r2t2 = 0, r2t3 = 0, r2t4 = 0, g4 = 0,
            rate.tp = 1, rate.cc = 0, p = .50, n1, n2, n3, n4)

cosa.bird4(score = NULL, order = 2, rhots = NULL,
           k1 = -6, k2 = 6, dists = "normal",
           cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
           n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL,
           n0 = c(10, 3, 100, 5 + g4), p0 = .499,
           constrain = "power", round = TRUE, max.power = FALSE,
           local.solver = c("LBFGS", "SLSQP"),
           power = .80, es = .25, alpha = .05, two.tailed = TRUE,
           rho2, rho3, rho4, omega2, omega3, omega4,
           g4 = 0, r21 = 0, r2t2 = 0, r2t3 = 0, r2t4 = 0)
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
  \item{rho4}{proportion of variance in the outcome between level 4 units (unconditional ICC4).}
  \item{omega2}{ratio of the treatment effect variance between level 2 units to the variance in the outcome between level 2 units.}
  \item{omega3}{ratio of the treatment effect variance between level 3 units to the variance in the outcome between level 3 units.}
  \item{omega4}{ratio of the treatment effect variance between level 4 units to the variance in the outcome between level 4 units.}
  \item{g4}{number of covariates at level 4.}
  \item{r21}{proportion of level 1 variance in the outcome explained by level 1 covariates.}
  \item{r2t2}{proportion of treatment effect variance between level 2 units explained by level 2 covariates.}
  \item{r2t3}{proportion of treatment effect variance between level 3 units explained by level 3 covariates.}
  \item{r2t4}{proportion of treatment effect variance between level 4 units explained by level 4 covariates.}
  \item{rate.tp}{treatment group participation rate.}
  \item{rate.cc}{control group crossover rate.}
  \item{p}{proportion of level 1 units in treatment condition.}
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{average number of level 2 units per level 3 unit.}
  \item{n3}{average number of level 3 units per level 4 unit.}
  \item{n4}{number of level 4 units.}

  \item{cn1}{marginal cost per level 1 unit in treatment and control conditions.}
  \item{cn2}{marginal cost per level 2 unit.}
  \item{cn3}{marginal cost per level 3 unit.}
  \item{cn4}{marginal cost per level 4 unit.}
  \item{cost}{total cost or budget.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with average when \code{p} is constrained by bounds.}
  \item{n0}{vector of starting values for \code{n1, n2, n3, n4} (positional). Starting values are replaced with averages when sample sizes are constrained by bounds.}
  \item{constrain}{character; \code{"cost"}, \code{"power"}, or \code{"mdes"}.}
  \item{round}{logical; \code{TRUE} for rounded COSA solution.}
  \item{max.power}{logical; \code{TRUE} for maximizing power instead of minimizing variance.}
  \item{local.solver}{subset of \code{c("LBFGS", "SLSQP")}.}
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
mdes.bird4(score.obj, order = 2,
           power = .80, rho2 = .20, rho3 = .10, rho4 = .05,
           omega2 = .30, omega3 = .30, omega4 = .30,
           g4 = 0, r2t4 = 0, n1 = 20, n2 = 3, n3 = 20, n4 = 10)

power.bird4(score.obj, order = 2,
            es = .152, rho2 = .20, rho3 = .10, rho4 = .05,
            omega2 = .30, omega3 = .30, omega4 = .30,
            g4 = 0, r2t4 = 0, n1 = 20, n2 = 3, n3 = 20, n4 = 10)

# optimal combination of sample sizes for level 1, level 2, level 3 and level 4
# that produce power = .80 (given range restrictions)
cosa.bird4(score.obj, order = 2,
           constrain = "power", power = .80,
           es = .25, rho2 = .20, rho3 = .10, rho4 = .05,
           omega2 = .30, omega3 = .30, omega4 = .30,
           g4 = 0, r2t4 = 0,
           n1 = c(15, 30), n2 = c(2, 5),
           n3 = c(10, 30), n4 = c(5, 20))
}
