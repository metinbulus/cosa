\name{plot}
\alias{plot.power}
\alias{plot.mdes}
\alias{plot.cosa}
\title{Power and MDES Curves}

\description{
  Plots statistical power or minimum detectable effect size curves with (1-\eqn{\alpha})x100 \% confidence interval for the design of interest.
}

\usage{
\method{plot}{power}(x, ypar = "mdes",  xpar = NULL,
          xlim = NULL, ylim = NULL,
          xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL,
          locate = FALSE, benchmark = NULL, ...)

\method{plot}{mdes}(x, ypar = "mdes",  xpar = NULL,
          xlim = NULL, ylim = NULL,
          xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL,
          locate = FALSE, benchmark = NULL, ...)

\method{plot}{cosa}(x, ypar = "mdes",  xpar = NULL,
          xlim = NULL, ylim = NULL,
          xlab = NULL, ylab = NULL,
          main = NULL, sub = NULL,
          locate = FALSE, benchmark = NULL, ...)
}

\arguments{
  \item{x}{an object returned from functions in \pkg{cosa} package.}
  \item{ypar}{character; \code{"mdes"} or  \code{"power"} on \code{y} axis.}
  \item{xpar}{character; one of the sample sizes on \code{x} axis.}
  \item{xlim}{limits for \code{xpar}.}
  \item{ylim}{limits for \code{ypar}.}
  \item{xlab}{\code{x} axis label.}
  \item{ylab}{\code{y} axis label.}
  \item{main}{title for the plot.}
  \item{sub}{subtitle for the plot.}
  \item{locate}{logical; \code{TRUE} locates parameter values for design \code{x} on the plot.}
  \item{benchmark}{benchmark line.}
  \item{...}{other graphical parameters to pass to \code{plot.new()}.}
}


\examples{
d1 <- mdes.bcrd3r2(rho2 = .10, rho3 = .20, omega3 = .30,
                   n1 = 20, n2 = 44, n3 = 50)
plot(d1, xpar = "n3", xlim = c(30, 100))
}
