\name{Vectorize BCOSSA}
\alias{vectorize.cosa}
\title{Vectorizes BCOSSA Solutions}

\description{
  Vectorizes bound constrained optimal sample size allocation (BCOSSA) solutions based on multiple sets of parameter values. This is particularly useful when multiple values of design parameters are to be considered.
}

\usage{
  vectorize.cosa(design, args.grid, args.names = NULL, ordered = TRUE, ncase = 10L)
}

\arguments{
  \item{design}{an object returned from one of the \code{cosa.<design>()} functions.}
  \item{args.grid}{vector or matrix: arguments' grid consisting of sets of parameter values. A vector of values (for a single parameter) or a matrix (for multiple parameters).}
  \item{args.names}{character list; arguments' names. Default option \code{args.names = NULL} uses column names from \code{args.grid}}
  \item{ordered}{logical: whether results should be ordered (cases with worst power rate or highest total cost are on top).}
  \item{ncase}{integer: number of cases to be subsetted, ignored if \code{ordered = FALSE}.}
}

\examples{

design <- cosa.crd2(order = 0, round = FALSE,
                    constrain = "power", power = .80,
                    cn1 = c(20, 10), cn2 = c(200, 50),
                    es = .25, rho2 = .10,
                    g2 = 3, r22 = .30,
                    n1 = NULL, n2 = NULL, p = NULL)

args.grid <- expand.grid(
  rho2 = seq(.15, .25, .05)
)

vectorize.cosa(design, args.grid = args.grid, ordered = FALSE)

}
