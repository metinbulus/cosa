\docType{package}
\name{cosa-package}
\alias{cosa-package}
\title{Bound Constrained Optimal Design of MRDDs and MRTs}

\description{
Bound Constrained Optimal Sample Size Allocation (BCOSSA) functions are designed to optimize sample sizes at one or more levels subject to budget, statistical power, or effect size constraints. BCOSSA can be found in the following forms; (i) under budgetary constraints given marginal costs per unit while minimizing sampling variance of the treatment effect (or, alternatively, while maximizing power rate), (ii) under statistical power or effect size (ES) constraints while minimizing the total cost, and (iii) under sample size constraints for one or more levels along with (i) or (ii). Specifying \code{rhots = 0} or \code{order = 0} produces results equivalent to corresponding random assignment designs, which means there is no relationship between the treatment [random] and the score variable. Therefore, BCOSSA functions also allow optimization of treatment group sampling rate (\code{p}) under unequal marginal costs when primary constraint is placed on the total cost. Different starting values and algorithms may produce different results when marginal cost information is not provided and sample sizes at two or more levels and \code{p} are optimized. In such cases, experimenting different starting values and/or comparing several algorithms may faciliate decisions regarding sample sizes and \code{p}.

Designs available in the \pkg{cosa} package:
\tabular{rccc}{
    \emph{Design}  \tab \emph{Total Levels} \tab \emph{Treatment Level}  \tab \emph{Top Level}\cr
    \code{ird}  \tab 1 \tab 1 \tab random  \cr
    \code{bird2} \tab 2 \tab 1 \tab random  \cr
    \code{bird2f1} \tab 2 \tab 1 \tab fixed   \cr
    \code{bird3} \tab 3 \tab 1 \tab random  \cr
    \code{bird4} \tab 4 \tab 1 \tab random  \cr
    \code{crd2}  \tab 2 \tab 2 \tab random  \cr
    \code{bcrd3f2} \tab 3 \tab 2 \tab fixed   \cr
    \code{bcrd3r2} \tab 3 \tab 2 \tab random  \cr
    \code{bcrd4r2} \tab 4 \tab 2 \tab random  \cr
    \code{crd3}  \tab 3 \tab 3 \tab random  \cr
    \code{bcrd4f3} \tab 4 \tab 3 \tab fixed   \cr
    \code{bcrd4r3} \tab 4 \tab 3 \tab random  \cr
    \code{crd4}  \tab 4 \tab 4 \tab random
}
  \code{ird}: individual-level regression discontinuity.
  \code{bird}: blocked individual-level regression discontinuity.
  \code{crd}: cluster-level regression discontinuity.
  \code{bcrd}: blocked cluster-level regression discontinuity.

  Design parameters follow a sequential order. Numbers at the end of a sequential parameter refers to the corresponding level. For example \code{rho2} is the proportion of variance in the outcome between level 2 units, \code{rho3} is the proportion of variance in the outcome between level 3 units. Similiarly, \code{r21} is the proportion of the variance in the outcome explained by level 1 covariates, \code{r22} is the proportion of the variance in the outcome explained by level 2 covariates and so on. Similiar naming conventions applies to other design parameters.
}
