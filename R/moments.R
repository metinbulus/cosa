tnorm.moment <- function(mu = 0, sigma = 1, k1 = -Inf, k2 = Inf, order = 1, central = FALSE) {
  if(order < 0) stop("'order' should be 0 or greater", call. = FALSE)
  # do not use very large or small numbers for bounds, instead define infinite bounds explicitly
  if(central) {
    m1 <- integrate(function(x) {x * msm::dtnorm(x = x, mean = mu, sd = sigma, lower = k1, upper = k2)}, lower = k1, upper = k2, subdivisions = 1000L)$value
    mgf <- function(x) {(x-m1)^order * msm::dtnorm(x = x, mean = mu, sd = sigma, lower = k1, upper = k2)}
    return(integrate(mgf, lower = k1, upper = k2, subdivisions = 1000L)$value)
  } else {
    mgf <- function(x) {x^order * msm::dtnorm(x = x, mean = mu, sd = sigma, lower = k1, upper = k2)}
    return(integrate(mgf, lower = k1, upper = k2, subdivisions = 1000L)$value)
  }
}

# moments of uniform
unif.moment <- function(k1 = 0, k2 = 1, order = 1, central = FALSE) {
  if(order < 0) stop("'order' should be 0 or greater", call. = FALSE)
  if(central) {
    M1 <- (k1 + k2)/2
    k1 <- k1 - M1
    k2 <- k2 - M1
  }
  #return((k2^(order + 1) - k1^(order + 1)) / ((order + 1) * (k2 - k1)))
  mgf <- function(x) {x^order * dunif(x = x, min = k1, max = k2)}
  return(integrate(mgf, lower = k1, upper = k2)$value)
}

# empirical moments
emp.moment <- function(x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE) {
  if(!is.numeric(x) | !is.vector(x)) stop("'x' should be a numeric vector", call. = FALSE)
  if(na.rm) {x <- x[!is.na(x)]}
  if(central) {
    x <- x - mean(x)
    if(absolute) {x <- abs(x)}
  }
  if(absolute) {x <- abs(x)}
  return(sum(x^order) / length(x))
}

# function to compute RDDE for arbitrary order
# w/ or w/o interaction up to ~ 9 or 13 for normal dist. with cutoff = 0
.rdde <- function(score = NULL, order = 1, interaction = FALSE,
                  treat.lower = FALSE, cutoff = 0, center = TRUE,
                  mu = 0, sigma = 1, k1 = -Inf, k2 =  Inf,
                  dists = "normal", ndraw = 1000, nsim = 1000) {

  pars <- c(cutoff, mu, sigma, k1, k2, ndraw, nsim)
  if(any(!is.numeric(pars) || length(pars) > 7)) stop("Incorrect value for one of the numeric argument", call. = FALSE)
  if(!is.logical(interaction)) stop("Non-logical value for argument 'interaction'", call. = FALSE)
  if(order %% 1 != 0) order <- round(order)
  if(order > 5) stop("Too many terms in the model may cause singularity issues, try 'order' < 6", call. = FALSE)

  if(is.null(score)) {
    sim <- TRUE
    temp <- matrix(NA, nrow = nsim, ncol = 6)
    for(i in 1:nsim) {
      ifelse(dists == "normal",
             score <- msm::rtnorm(ndraw, mean = mu, sd = sigma, lower = k1, upper = k2),
             score <- runif(ndraw, k1, k2))
      ifelse(isFALSE(treat.lower),
             treatment <- ifelse(score > cutoff, 1, 0),
             treatment <- ifelse(score < cutoff, 1, 0))
      if(center == TRUE) score <- score - cutoff
      smat <- matrix(NA, nrow = length(score), ncol = order)
      for(j in 1:order) {
        smat[, j] <- score^j
      }
      ifelse(interaction,
             xmat <- cbind(1, treatment, smat, treatment * smat),
             xmat <- cbind(1, treatment, smat)
      )
      pi <- mean(treatment)
      rdde.temp <- try(solve(t(xmat) %*% xmat)[2,2] * pi * (1 - pi) * (ndraw - 1))
      if(inherits(rdde.temp, "try-error")) rdde.temp <- NA
      temp[i,1] <- pi
      temp[i,2] <- mean(score)
      temp[i,3] <- var(score)
      temp[i,4] <- min(score)
      temp[i,5] <- max(score)
      temp[i,6] <- rdde.temp
    }
    p <- mean(temp[,1])
    mu <- mean(temp[,2])
    sigma <- sqrt(sum(temp[,3]) / nsim)
    k1 <- mean(temp[,4])
    k2 <- mean(temp[,5])
    rdde <- mean(temp[,6], na.rm = TRUE)
    if(length(temp[,6][is.na(temp[,6])]) / ndraw > .05) warning("Computational singularity issues, interpret results with caution", call. = FALSE)
    message(cat("RDDE for", dists, "distribution \n based on simulation (approx. population moments)"))
  } else{
    if(!is.vector(score)) stop("Score variable should be a vector", call. = FALSE)
    sim <- FALSE
    score <- score[score < k2 & score > k1]
    ifelse(isFALSE(treat.lower),
           treatment <- ifelse(score > cutoff, 1, 0),
           treatment <- ifelse(score < cutoff, 1, 0))
    p <- mean(treatment)
    if(center == TRUE) score <- score - cutoff
    smat <- matrix(NA, nrow = length(score), ncol = order)
    for(i in 1:order) {
      smat[, i] <- score^i
    }
    ifelse(interaction,
           xmat <- cbind(1, treatment, smat, treatment * smat),
           xmat <- cbind(1, treatment, smat)
    )
    mu <- mean(score)
    sigma <- sd(score)
    k1 <- min(score)
    k2 <- max(score)
    rdde <- try(solve(t(xmat) %*% xmat)[2,2] * p * (1 - p) * length(score))
    if(inherits(rdde, "try-error")) rdde <- NA
    message(cat("RDDE for empirical score distribution \n based on sample moments"))
  }

  if(!is.null(score)) dists <- "empirical"
  parms <- list(dists = dists, mu = mu, sigma = sigma, k1 = k1, k2 = k2, sim = sim, ndraw = ndraw, nsim = nsim)
  rdde.out <- list(parms = parms, cutoff = cutoff, treat.lower = treat.lower, p = p,
                    order = order, interaction = interaction, center = center,
                    rdde = rdde)
  class(rdde.out) <- "rdde"

  return(invisible(rdde.out))
}


# function to compute RDDE
# mu: uncentered mean of the score var
# k1: uncentered lower bound for the score var
# k2: uncentered upper bound for the score var
inspect.score <- function(score = NULL, p = NULL, cutoff = NULL,
                           treat.lower = FALSE, order = 1, interaction = FALSE,
                           mu = 0, sigma = 1, k1 = -Inf, k2 =  Inf,
                           dists = "normal", sim = FALSE, ndraw = 1000, nsim = 1000) {

  pars <- c(p, cutoff, mu, sigma, k1, k2, ndraw, nsim)
  if(any(!is.numeric(pars) || length(pars) > 8)) stop("Incorrect value for one of the numeric argument", call. = FALSE)
  if(is.null(cutoff) & is.null(p)) stop("Specify one of the 'cutoff' or 'p' arguments", call. = FALSE)
  if(!is.null(cutoff) & !is.null(p)) warning("Ignoring 'cutoff' using 'p'", call. = FALSE)
  if(!is.logical(interaction)) stop("Non-logical value for argument 'interaction'", call. = FALSE)
  if(order < 1) stop("'order' < 1?", call. = FALSE)
  if(order %% 1 != 0) {
    order <- round(order)
    warning("'order' is rounded to the closest integer", call. = FALSE)
  }

  if(!is.null(p)) {
    ifelse(treat.lower, q <- p, q <- 1 - p)
    if(is.null(score)) {
      ifelse(dists == "normal",
             cutoff <- msm::qtnorm(q, mean = mu, sd = sigma, lower = k1, upper = k2),
             cutoff <- qunif(q, min = k1, max = k2))
    } else {
      cutoff <- quantile(score, q)
    }
  }

  if(!is.null(cutoff) & is.null(p)) {
    ifelse(dists == "normal",
           q <- msm::ptnorm(cutoff, mean = mu, sd = sigma, lower = k1, upper = k2),
           q <- punif(cutoff, min = k1, max = k2))
    ifelse(treat.lower, p <- q, p <- 1 - q)
  }
  if(round(p,2) < .05 | round(p,2) > .95) stop("'cutoff' in the vicinity of bounds", call. = FALSE)

  if(order < 3) {
    if(is.null(score) & isTRUE(sim)) {
      temp <- matrix(NA, nrow = nsim, ncol = 15)
      for(i in 1:nsim) {
        ifelse(dists == "normal",
               score <- msm::rtnorm(ndraw, mean = mu, sd = sigma, lower = k1, upper = k2),
               score <- runif(ndraw, k1, k2))
        ifelse(isFALSE(treat.lower),
               treatment <- ifelse(score > cutoff, 1, 0),
               treatment <- ifelse(score < cutoff, 1, 0))
        score <- score - cutoff
        temp[i,1] <- mean(treatment)
        temp[i,2] <- cor(treatment, score)
        temp[i,3] <- cor(treatment, score^2)
        temp[i,4] <- cor(score, score^2)
        temp[i,5] <- cor(treatment, treatment*score)
        temp[i,6] <- cor(treatment, treatment*score^2)
        temp[i,7] <- cor(score, treatment*score)
        temp[i,8] <- cor(score, treatment*score^2)
        temp[i,9] <- cor(score^2, treatment*score)
        temp[i,10] <- cor(score^2, treatment*score^2)
        temp[i,11] <- cor(treatment*score, treatment*score^2)
        temp[i,12] <- mean(score)
        temp[i,13] <- var(score)
        temp[i,14] <- min(score)
        temp[i,15] <- max(score)
      }
      p <- mean(temp[,1])
      rhots <- mean(temp[,2])
      rhots2 <- mean(temp[,3])
      rhoss2 <- mean(temp[,4])
      rhotts <- mean(temp[,5])
      rhotts2 <- mean(temp[,6])
      rhosts <- mean(temp[,7])
      rhosts2 <- mean(temp[,8])
      rhos2ts <- mean(temp[,9])
      rhos2ts2 <- mean(temp[,10])
      rhotsts2 <- mean(temp[,11])
      mu <- mean(temp[,12])
      sigma <- sqrt(sum(temp[,13]) / nsim)
      k1 <- mean(temp[,14])
      k2 <- mean(temp[,15])
      message(cat("RDDE for", dists, "distribution \n based on simulation (approx. population moments)"))
    } else if(is.null(score) & isFALSE(sim)) {
      if(dists == "normal") {
        m1 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = k2-cutoff, central = FALSE, order = 1)
        m2 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = k2-cutoff, central = FALSE, order = 2)
        m3 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = k2-cutoff, central = FALSE, order = 3)
        m4 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = k2-cutoff, central = FALSE, order = 4)
        if(treat.lower) {
          m1s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = 0, central = FALSE, order = 1)
          m2s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = 0, central = FALSE, order = 2)
          m3s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = 0, central = FALSE, order = 3)
          m4s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = k1-cutoff, k2 = 0, central = FALSE, order = 4)
        } else {
          m1s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = 0, k2 = k2-cutoff, central = FALSE, order = 1)
          m2s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = 0, k2 = k2-cutoff, central = FALSE, order = 2)
          m3s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = 0, k2 = k2-cutoff, central = FALSE, order = 3)
          m4s0 <- tnorm.moment(mu = mu-cutoff, sigma = sigma, k1 = 0, k2 = k2-cutoff, central = FALSE, order = 4)
        }
        # mu <- m1 # biased
        # sigma <- sqrt(tnorm.moment(mu = mu, sigma = sigma, k1 = k1, k2 = k2, central = TRUE, order = 2)) # biased
        mu <- mu-cutoff
        sigma <- sigma
        k1 <- k1-cutoff
        k2 <- k2-cutoff
      } else if(dists == "uniform") {
        m1 <- unif.moment(k1 = k1-cutoff, k2 = k2-cutoff, order = 1, central = FALSE)
        m2 <- unif.moment(k1 = k1-cutoff, k2 = k2-cutoff, order = 2, central = FALSE)
        m3 <- unif.moment(k1 = k1-cutoff, k2 = k2-cutoff, order = 3, central = FALSE)
        m4 <- unif.moment(k1 = k1-cutoff, k2 = k2-cutoff, order = 4, central = FALSE)
        if(treat.lower) {
          m1s0 <- unif.moment(k1 = k1-cutoff, k2 = 0, order = 1, central = FALSE)
          m2s0 <- unif.moment(k1 = k1-cutoff, k2 = 0, order = 2, central = FALSE)
          m3s0 <- unif.moment(k1 = k1-cutoff, k2 = 0, order = 3, central = FALSE)
          m4s0 <- unif.moment(k1 = k1-cutoff, k2 = 0, order = 4, central = FALSE)
        } else {
          m1s0 <- unif.moment(k1 = 0, k2 = k2-cutoff, order = 1, central = FALSE)
          m2s0 <- unif.moment(k1 = 0, k2 = k2-cutoff, order = 2, central = FALSE)
          m3s0 <- unif.moment(k1 = 0, k2 = k2-cutoff, order = 3, central = FALSE)
          m4s0 <- unif.moment(k1 = 0, k2 = k2-cutoff, order = 4, central = FALSE)
        }
        mu <- m1
        sigma <- sqrt(unif.moment(k1 = k1, k2 = k2, order = 2, central = TRUE))
        k1 <- k1-cutoff
        k2 <- k2-cutoff
      } else if(dists == "empirical"){
        stop("Empirical score variable is missing", call. = FALSE)
      } else {
        stop("Distribution type not supported", call. = FALSE)
      }
      rhots <- p*(m1s0 - m1) /  sqrt(p*(1-p)*(m2 - m1^2))
      rhots2 <- p*(m2s0 - m2) /  sqrt(p*(1-p)*(m4 - (m2)^2))
      rhoss2 <- (m3 - m1*m2) / sqrt((m2 - m1^2) * (m4 - m2^2))
      rhotts <- sqrt(1-p) * m1s0 / sqrt(m2s0 - p * m1s0^2)
      rhotts2 <- sqrt(1-p) * m2s0 / sqrt(m4s0 - p * m2s0^2)
      rhosts <- sqrt(p) * (m2s0 - m1 * m1s0) / sqrt((m2 - m1^2) * (m2s0 - p * m1s0^2))
      rhosts2 <- sqrt(p) * (m3s0 - m1 * m2s0) / sqrt((m2 - m1^2) * (m4s0 - p * m2s0^2))
      rhos2ts <- sqrt(p) * (m3s0 - m2 * m1s0) / sqrt((m4 - m2^2) * (m2s0 - p * m1s0^2))
      rhos2ts2 <- sqrt(p) * (m4s0 - m2 * m2s0) / sqrt((m4 - m2^2) * (m4s0 - p * m2s0^2))
      rhotsts2 <- (m3s0 - p * m1s0 * m2s0) / sqrt((m2s0 - p * m1s0^2) * (m4s0 - p * m2s0^2))
      message(cat("RDDE for", dists, "distribution \n based on population moments"))
    } else if(!is.null(score)){
      if(!is.vector(score)) stop("Score variable should be a vector", call. = FALSE)
      score <- score[score < k2 & score > k1]
      ifelse(isFALSE(treat.lower),
             treatment <- ifelse(score > cutoff, 1, 0),
             treatment <- ifelse(score < cutoff, 1, 0))
      p <- mean(treatment)
      score <- score - cutoff
      rhots <- cor(treatment, score)
      rhots2 <- cor(treatment, score^2)
      rhoss2 <- cor(score, score^2)
      rhotts <- cor(treatment, treatment*score)
      rhotts2 <- cor(treatment, treatment*score^2)
      rhosts <- cor(score, treatment*score)
      rhosts2 <- cor(score, treatment*score^2)
      rhos2ts <- cor(score^2, treatment*score)
      rhos2ts2 <- cor(score^2, treatment*score^2)
      rhotsts2 <- cor(treatment*score, treatment*score^2)
      mu <- mean(score)
      sigma <- sd(score)
      k1 <- min(score)
      k2 <- max(score)
      dists <- "empirical"
      message(cat("RDDE for empirical score distribution \n based on sample moments"))
    }
    corMat <- matrix(c(1, rhots, rhots2, rhotts, rhotts2,
                       rhots, 1, rhoss2, rhosts, rhosts2,
                       rhots2, rhoss2, 1, rhos2ts, rhos2ts2,
                       rhotts, rhosts, rhos2ts, 1, rhotsts2,
                       rhotts2, rhosts2, rhos2ts2, rhotsts2, 1),
                     nrow = 5, ncol = 5, byrow = TRUE)
    if(interaction) {
      rdde <- switch(order,
                     "1" = (1 - rhosts^2) / (1 - rhots^2 - rhotts^2 - rhosts^2 + 2*rhots*rhotts*rhosts),
                     "2" = solve(corMat)[1,1])
    } else {
      rdde <- switch(order,
                     "1" = 1 / (1 - rhots^2),
                     "2" = (1 - rhoss2^2) / (1 - rhots^2 - rhots2^2 - rhoss2^2 + 2*rhots*rhots2*rhoss2))
    }
    parms <- list(dists = dists, mu = mu, sigma = sigma, k1 = k1, k2 = k2, sim = sim, ndraw = ndraw, nsim = nsim)
    score.out <- list(parms = parms, cutoff = cutoff, treat.lower = treat.lower, p = p,
                      order = order, interaction = interaction,
                      rdde = rdde)
  } else {
    score.out <- .rdde(score = score, treat.lower = treat.lower, dists = dists,
                       order = order, interaction = interaction, cutoff = cutoff, center = TRUE,
                       mu = mu, sigma = sigma, k1 = k1, k2 = k2, ndraw = ndraw, nsim = nsim)
    rdde <- score.out$rdde
  }

  class(score.out) <- "score"
  cat("\n---------------------------------------",
      "\nPolynomial order =", order, "\nInteraction w/ treatment =",  interaction,
      "\nTreat if score < cutoff =", treat.lower,
      "\nCutoff =", round(cutoff,3), "| p =", round(p,3), "\nRDDE =", round(rdde,3), "\n\n")
  return(invisible(score.out))
}
