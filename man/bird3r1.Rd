\name{bird3}
\alias{bird3r1}
\alias{bird3}
\alias{cosa.bird3r1}
\alias{cosa.bird3}
\alias{mdes.bird3r1}
\alias{mdes.bird3}
\alias{power.bird3r1}
\alias{power.bird3}
\title{Blocked (Random) Individual-level Regression Discontinuity (Three-level Design, Discontinuity at Level 1)}

\description{
 Use \code{mdes.bird3()} to calculate minimum detectable effect size, \code{power.bird3()} to calculate statistical power, and \code{cosa.bird3()} for bound constrained optimal sample size allocation (BCOSSA).
}

\usage{
mdes.bird3(score = NULL, dists = "normal", k1 = -6, k2 = 6,
           order = 1, interaction = FALSE,
           treat.lower = TRUE, cutoff = 0, p = NULL,
           power = .80, alpha = .05, two.tailed = TRUE, df = n3 - g3 - 1,
           rho2, rho3, omega2, omega3, r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
           rate.tp = 1, rate.cc = 0, n1, n2, n3)

power.bird3(score = NULL, dists = "normal", k1 = -6, k2 = 6,
            order = 1, interaction = FALSE,
            treat.lower = TRUE, cutoff = 0, p = NULL,
            es = .25, alpha = .05, two.tailed = TRUE, df = n3 - g3 - 1,
            rho2, rho3, omega2, omega3, r21 = 0, r2t2 = 0, r2t3 = 0, g3 = 0,
            rate.tp = 1, rate.cc = 0, n1, n2, n3)

cosa.bird3(score = NULL, dists = "normal", k1 = -6, k2 = 6, rhots = NULL,
           order = 1, interaction = FALSE,
           treat.lower = TRUE, cutoff = 0, p = NULL,
           cn1 = 0, cn2 = 0, cn3 = 0, cost = NULL,
           n1 = NULL, n2 = NULL, n3 = NULL,
           n0 = c(10, 3, 100), p0 = .499,
           constrain = "power", round = TRUE, max.power = FALSE,
           local.solver = c("LBFGS", "SLSQP"),
           power = .80, es = .25, alpha = .05, two.tailed = TRUE,
           rho2, rho3, omega2, omega3,
           g3 = 0, r21 = 0, r2t2 = 0, r2t3 = 0)
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
  \item{p}{proportion of level 1 units in the treatment condition.}

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
  \item{n1}{average number of level 1 units per level 2 unit.}
  \item{n2}{average number of level 2 units (blocks) per level 3 unit.}
  \item{n3}{number of level 3 units (blocks).}


  \item{cn1}{marginal costs per level 1 unit in treatment and control conditions (positional), e.g. \code{c(10, 5)}.}
  \item{cn2}{marginal cost per level 2 unit.}
  \item{cn3}{marginal cost per level 3 unit.}
  \item{cost}{total cost or budget. Ignored when \code{constrain = "power"} or \code{constrain = "es"}.}
  \item{p0}{starting value for \code{p} when \code{rhots = 0} and \code{p = NULL}. Starting value is replaced with the average when \code{p} is constrained by bounds.}
  \item{n0}{vector of starting values for \code{n1, n2, n3} (positional). Starting values are replaced with the averages when sample sizes are constrained by bounds.}
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

power.bird3(score.obj,
            es = 0.25, rho2 = .20, rho3 = .10,
            omega2 = .30, omega3 = .30,
            g3 = 0, r2t3 = 0,
            n1 = 20, n2 = 3, n3 = 20)

# minimum required number of level 1 units for each one of the level 2 block
cosa.bird3(score.obj,
           es = 0.25, rho2 = .20, rho3 = .10,
           omega2 = .30, omega3 = .30,
           g3 = 0, r2t3 = 0,
           n1 = NULL, n2 = 3, n3 = 20)
}
