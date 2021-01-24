plot.power <- plot.cosa <- plot.mdes <- function(x, score = NULL, ypar = "mdes",  xpar = NULL,
                                                 xlim = NULL, ylim = NULL,
                                                 xlab = NULL, ylab = NULL,
                                                 main = NULL, sub = NULL,
                                                 locate = FALSE, benchmark = NULL, ...) {

  if(!inherits(x, "mdes") & !inherits(x, "power") & !inherits(x, "cosa")) {
    stop("'x' should be an object returned from functions in the 'cosa' library")
  }
  if(inherits(x, "main") | inherits(x, "mod") | inherits(x, "med")) {
    stop("'x' should be an object returned from functions in the 'cosa' library")
  }
  if(!is.logical(locate)) {
    stop("Non-logical value for argument 'locate'", call. = FALSE)
  }

  ypar <- tolower(ypar)
  if(ypar == "mdes") {
    if(inherits(x, "power")) {
      capture.output(x <- .power2mdes(x, score))
    } else if(inherits(x, "cosa")) {
      capture.output(x <- .cosa2mdes(x, score))
    }
  } else if(ypar == "power") {
    if(inherits(x, "mdes")) {
      capture.output(x <- .mdes2power(x, score))
    } else if(inherits(x, "cosa")) {
      capture.output(x <- .cosa2power(x, score))
    }
  } else {
    stop("Incorrect value for 'ypar'")
  }

  design <- class(x)[2]
  block <- substr(design, nchar(design) - 1, nchar(design) - 1)
  nlevels <- substr(design, nchar(design) - 2, nchar(design) - 2)
  order <- x$parms$order
  inter <- x$parms$interaction

  N <- c("n1", "n2", "n3", "n4")
  if(is.null(xpar)) {
    if(block == "r") {
      xpar <- switch(nlevels,
                     "1"= "n1",
                     "2"= "n2",
                     "3"= "n3",
                     "4"= "n4")
    } else {
      xpar <- switch(nlevels,
                     "2"= "n1",
                     "3"= "n2",
                     "4"= "n3")
    }

    if(is.null(xlim)) {
      if(block == "r") {
        xseq <- switch(nlevels,
                       "1"= seq(round(x$parms$g1) + order * (1 + inter) + 3, 1.5 * round(x$parms$n1), by = .5),
                       "2"= seq(round(x$parms$g2) + order * (1 + inter) + 3, 1.5 * round(x$parms$n2), by = .5),
                       "3"= seq(round(x$parms$g3) + order * (1 + inter) + 3, 1.5 * round(x$parms$n3), by = .5),
                       "4"= seq(round(x$parms$g4) + order * (1 + inter) + 3, 1.5 * round(x$parms$n4), by = .5))
      } else {
        xseq <- switch(nlevels,
                       "2"= seq(round(x$parms$g1) + order * (1 + inter) + 3, 1.5 * round(x$parms$n1), by = .5),
                       "3"= seq(round(x$parms$g2) + order * (1 + inter) + 3, 1.5 * round(x$parms$n2), by = .5),
                       "4"= seq(round(x$parms$g3) + order * (1 + inter) + 3, 1.5 * round(x$parms$n3), by = .5))
      }
    } else {
      if(!is.numeric(xlim)) {
        stop("Incorrect values for argument 'xlim'")
      }
      if(length(xlim) > 2) {
        xlim <- range(xlim)
      }
      xseq <- seq(xlim[1], xlim[2], .5)
    }
  } else if(!is.null(xpar)) {
    if(!xpar %in% N) {
      stop("Incorrect values for argument 'xpar'", call. = FALSE)
    }
    if(is.null(xlim)) {
      stop("Argument 'xlim' is NULL", call. = FALSE)
    } else {
      if(!is.numeric(xlim)) {
        stop("Incorrect values for argument 'xlim'")
      }
      if(length(xlim) > 2) {
        xlim <- range(xlim)
      }
      xseq <- seq(xlim[1], xlim[2], .5)
    }
  }

  if(x$parms$dists == "empirical") {
    ifelse(is.null(score),
           stop("Score object is missing", call. = FALSE),
           x$parms$score <- score)
  }
  names.parms <-  names(x$parms)
  idx <- match(xpar, names.parms)
  x0 <- x$parms[[idx]]
  ifelse(ypar == "mdes", y0 <- x$mdes[1], y0 <- x$power)

  yout <- matrix(NA, nrow = length(xseq), ncol = 3)
  capture.output({
      for(i in 1:length(xseq)) {
        #x$parms["p"]<- list(NULL)
        x$parms[[idx]] <- xseq[i]
        if(tolower(ypar) == "mdes") {
          yout[i,] <- do.call(paste("mdes", class(x)[2], sep = "."), x$parms)$mdes
        } else if(tolower(ypar) == "power") {
          yout[i,1] <- do.call(paste("power", class(x)[2], sep = "."), x$parms)$power
        }
      }
  })

  if(is.null(ylim)) {
    ifelse(ypar == "mdes",
           ylim <- range(yout) + c(-.30, .30),
           ylim <- c(.01, .99))
  } else {
    if(!is.numeric(ylim)) {
      stop("Incorrect values for argument 'ylim'", call. = FALSE)
    }
    if(length(ylim) > 2) {
      ylim <- range(ylim)
    }
    if(ypar == "power") {
      if(min(ylim) < 0 | max(ylim) > 1) {
        stop("Incorrect values for argument 'ylim'", call. = FALSE)
      }
    }
  }

  ifelse(!is.null(ylab), ylab,
         ifelse(ypar == "mdes",
                ylab <- "Minimum Detectable Effect Size",
                ylab <- "Statistical Power"))

  plot.new()
  plot.window(xlim = range(xseq),
              ylim = ylim, ...)
  polygon(c(rev(xseq), xseq), c(rev(yout[,3]), yout[,2]), border = NA,
          col = adjustcolor(4, alpha.f = 0.2))
  lines(xseq, yout[,1], col = adjustcolor(4, alpha.f = 0.5), lty = 1, lwd = 2)
  if(ypar == "mdes") {
    lines(xseq, yout[,2], col = adjustcolor(4, alpha.f = 0.2), lty = 1, lwd = 1.5)
    lines(xseq, yout[,3], col = adjustcolor(4, alpha.f = 0.2), lty = 1, lwd = 1.5)
  }
  title(main = main, sub = sub,
        xlab = ifelse(!is.null(xlab), xlab, xpar),
        ylab = ylab)
  axis(1)
  axis(2)
  box()

  if(locate) {
    points(x0, y0, pch=21, bg = adjustcolor(2, alpha.f = 0.5), cex=1.5)
    abline(v = x0, lty = 5, col = adjustcolor(2, alpha.f = 0.5))
  }
  if(is.null(benchmark)) {
    abline(h = ifelse(ypar == "mdes", .20, .80), lty = 5, col = adjustcolor(2, alpha.f = 0.5))
  } else {
    abline(h = benchmark, lty = 5, col = adjustcolor(2, alpha.f = 0.5))
  }

}
