mdes.bcrd4f3 <- function(score = NULL, dists = "normal", k1 = -6, k2 = 6,
                         order = 1, interaction = FALSE, treat.lower = TRUE, cutoff = 0, p = NULL,
                         power = .80, alpha = .05, two.tailed = TRUE,
                         df = n4 * (n3 - 2) - g3 - order * (1 + interaction),
                         rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
                         rate.tp = 1, rate.cc = 0, n1, n2, n3, n4) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(df < 1) stop("Insufficient degrees of freedom", call. = FALSE)
  if(!is.null(score) & order == 0) warning("Ignoring information from the 'score' object \n", call. = FALSE)
  if(order == 0) {
    d <- 1
    if(is.null(p)) stop("'p' cannot be NULL in random assignment designs", call. = FALSE)
    idx.score <- intersect(c("dists", "k1", "k2", "interaction", "treat.lower", "cutoff"),  names(user.parms))
    if(length(idx.score) > 0)  cat("\nCAUTION: Ignoring argument(s):",
                                   sQuote(names(user.parms[idx.score])), "\n")
    ifelse(treat.lower, cutoff <- p, cutoff <- 1 - p)
    interaction <- FALSE
    dists <- "uniform"
    k1 <- 0
    k2 <- 1
  } else if(order %in% 1:8) {
    if(is.null(score)) {
      score <- inspect.score(order = order, interaction = interaction,
                             treat.lower = treat.lower, cutoff = cutoff,
                             p = p, k1 = k1, k2 = k2, dists = dists)
    } else {
      if("p" %in% names(user.parms)) warning("Using 'p' from the 'score' object, ignoring 'p' in the function call", call. = FALSE)
      if(!inherits(score, "score")) {
        score <- inspect.score(score = score, order = order, interaction = interaction,
                               treat.lower = treat.lower, cutoff = cutoff,
                               p = p, k1 = k1, k2 = k2, dists = dists)
      } else {
        idx.score <- intersect(c("dists", "k1", "k2", "order", "interaction", "treat.lower", "p", "cutoff"),  names(user.parms))
        if(length(idx.score) > 0)  cat("\nCAUTION: 'score' object overwrites argument(s):",
                                       sQuote(names(user.parms[idx.score])), "\n")
      }
    }
    d <- score$rdde
    p <- score$p
    cutoff <- score$cutoff
    treat.lower <- score$treat.lower
    order <- score$order
    interaction <- score$interaction
    dists <- score$parms$dists
    k1 <- score$parms$k1
    k2 <- score$parms$k2
  } else if(order > 8) {
    stop("'order' > 8 is not allowed", call. = FALSE)
  }

  sse <- (1/(rate.tp - rate.cc)) * sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
               (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(dists = dists, k1 = k1, k2 = k2,
                                order = order, interaction = interaction,
                                treat.lower = treat.lower, p = p, cutoff = cutoff,
                                power = power, alpha = alpha, two.tailed = two.tailed,
                                rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                                g3 = g3, rate.tp = rate.tp, rate.cc = rate.cc,
                                n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                   df = df,
                   sse = sse,
                   mdes = mdes)
  class(mdes.out) <- c("mdes", "bcrd4f3")
  .summary.mdes(mdes.out)
  return(invisible(mdes.out))
}

power.bcrd4f3 <- function(score = NULL, dists = "normal", k1 = -6, k2 = 6,
                          order = 1, interaction = FALSE, treat.lower = TRUE, cutoff = 0, p = NULL,
                          es = .25, alpha = .05, two.tailed = TRUE,
                          df = n4 * (n3 - 2) - g3 - order * (1 + interaction),
                          rho2, rho3, r21 = 0, r22 = 0, r23 = 0, g3 = 0,
                          rate.tp = 1, rate.cc = 0, n1, n2, n3, n4) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms)

  if(df < 1) stop("Insufficient degrees of freedom", call. = FALSE)
  if(!is.null(score) & order == 0) warning("Ignoring information from the 'score' object \n", call. = FALSE)
  if(order == 0) {
    d <- 1
    if(is.null(p)) stop("'p' cannot be NULL in random assignment designs", call. = FALSE)
    idx.score <- intersect(c("dists", "k1", "k2", "interaction", "treat.lower", "cutoff"),  names(user.parms))
    if(length(idx.score) > 0)  cat("\nCAUTION: Ignoring argument(s):",
                                   sQuote(names(user.parms[idx.score])), "\n")
    ifelse(treat.lower, cutoff <- p, cutoff <- 1 - p)
    interaction <- FALSE
    dists <- "uniform"
    k1 <- 0
    k2 <- 1
  } else if(order %in% 1:8) {
    if(is.null(score)) {
      score <- inspect.score(order = order, interaction = interaction,
                             treat.lower = treat.lower, cutoff = cutoff,
                             p = p, k1 = k1, k2 = k2, dists = dists)
    } else {
      if("p" %in% names(user.parms)) warning("Using 'p' from the 'score' object, ignoring 'p' in the function call", call. = FALSE)
      if(!inherits(score, "score")) {
        score <- inspect.score(score = score, order = order, interaction = interaction,
                               treat.lower = treat.lower, cutoff = cutoff,
                               p = p, k1 = k1, k2 = k2, dists = dists)
      } else {
        idx.score <- intersect(c("dists", "k1", "k2", "order", "interaction", "treat.lower", "p", "cutoff"),  names(user.parms))
        if(length(idx.score) > 0)  cat("\nCAUTION: 'score' object overwrites argument(s):",
                                       sQuote(names(user.parms[idx.score])), "\n")
      }
    }
    d <- score$rdde
    p <- score$p
    cutoff <- score$cutoff
    treat.lower <- score$treat.lower
    order <- score$order
    interaction <- score$interaction
    dists <- score$parms$dists
    k1 <- score$parms$k1
    k2 <- score$parms$k2
  } else if(order > 8) {
    stop("'order' > 8 is not allowed", call. = FALSE)
  }

  sse <- (1/(rate.tp - rate.cc)) * sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                     rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                     (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1)))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(dists = dists, k1 = k1, k2 = k2,
                                  order = order, interaction = interaction,
                                  treat.lower = treat.lower, p = p, cutoff = cutoff,
                                  es = es, alpha = alpha, two.tailed = two.tailed,
                                  rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                                  g3 = g3, rate.tp = rate.tp, rate.cc = rate.cc,
                                  n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                     df = df,
                     sse = sse,
                     power = power)
  class(power.out) <- c("power", "bcrd4f3")
  .summary.power(power.out)
  return(invisible(power.out))
}

cosa.bcrd4f3 <- function(score = NULL, dists = "normal", k1 = -6, k2 = 6, rhots = NULL,
                         order = 1, interaction = FALSE,
                         treat.lower = TRUE, cutoff = 0, p = NULL,
                         cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
                         n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL,
                         n0 = c(10, 3, 100 + g3 + order * (1 + interaction), 5), p0 = .499,
                         constrain = "power", round = TRUE, max.power = FALSE,
                         local.solver = c("LBFGS", "SLSQP"),
                         power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                         rho2, rho3, g3 = 0, r21 = 0, r22 = 0, r23 = 0) {

  user.parms <- as.list(match.call())
  .error.handler(user.parms, fun = "cosa")

  if(!is.null(rhots)) {
    if(rhots == 0) {
      if(order != 0) {
        order <- 0
        warning("'order' argument is ignored, forcing 'order = 0' because 'rhots = 0'", call. = FALSE)
      }
    } else {
      stop("'rhots' argument will be removed in the future, arbitrary correlations are not allowed,
              use inspect.score() function instead", call. = FALSE)
    }
  }

  if(!is.null(score) & order == 0) warning("Ignoring information from the 'score' object \n", call. = FALSE)
  if(order == 0) {
    d <- 1
    idx.score <- intersect(c("dists", "k1", "k2", "interaction", "treat.lower", "cutoff"),  names(user.parms))
    if(length(idx.score) > 0)  cat("\nCAUTION: Ignoring argument(s):",
                                   sQuote(names(user.parms[idx.score])), "\n")
    cutoff <- NA
    interaction <- FALSE
    dists <- "uniform"
    k1 <- 0
    k2 <- 1
  } else if(order %in% 1:8) {
    if(is.null(score)) {
      score <- inspect.score(order = order, interaction = interaction,
                             treat.lower = treat.lower, cutoff = cutoff,
                             p = p, k1 = k1, k2 = k2, dists = dists)
    } else {
      if("p" %in% names(user.parms)) warning("Using 'p' from the 'score' object, ignoring 'p' in the function call", call. = FALSE)
      if(!inherits(score, "score")) {
        score <- inspect.score(score = score, order = order, interaction = interaction,
                               treat.lower = treat.lower, cutoff = cutoff,
                               p = p, k1 = k1, k2 = k2, dists = dists)
      } else {
        idx.score <- intersect(c("dists", "k1", "k2", "order", "interaction", "treat.lower", "p", "cutoff"),  names(user.parms))
        if(length(idx.score) > 0)  cat("\nCAUTION: 'score' object overwrites argument(s):",
                                       sQuote(names(user.parms[idx.score])), "\n")
      }
    }
    d <- score$rdde
    p <- score$p
    cutoff <- score$cutoff
    treat.lower <- score$treat.lower
    order <- score$order
    interaction <- score$interaction
    dists <- score$parms$dists
    k1 <- score$parms$k1
    k2 <- score$parms$k2
  } else if(order > 8) {
    stop("'order' > 8 is not allowed", call. = FALSE)
  }

  fun <- "cosa.bcrd4f3"
  lb <- c(1, 1, g3 + order * (1 + interaction) + 3, 1)
  if(!is.null(n3)) {
    if(n3[1] < lb[3]) warning("Lower bound for 'n3' may violate minimum degrees of freedom requirement", call. = FALSE)
  }

  .df <- quote(n4 * (n3 - 2) - g3 - order * (1 + interaction))
  .sse <- quote(sqrt(d * (rho3 * (1 - r23) / (p * (1 - p) * n4 * n3) +
                            rho2 * (1 - r22) / (p * (1 - p) * n4 * n3 * n2) +
                            (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n4 * n3 * n2 * n1))))
  .cost <- quote(n4 * cn4 +
                   n4 * n3 * (cn3[2] + p * (cn3[1] - cn3[2])) +
                   n4 * n3 * n2 * (cn2[2] + p * (cn2[1] - cn2[2])) +
                   n4 * n3 * n2 * n1 * (cn1[2] + p * (cn1[1] - cn1[2])))

  # NOTE: numerical derivatives for equality constraint (on power rate) are not precise
  # possibly due to numerical instability resulting from small step size
  .var.jacob <- expression(
    c(
      -d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3 * n4 * n1^2),

      -d * rho2 * (1 - r22) / (p * (1 - p) * n2^2 * n3 * n4) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2^2 * n3 * n4 * n1),

      -d * rho3 * (1 - r23) / (p * (1 - p) * n3^2 * n4) -
        d * rho2 * (1 - r22) / (p * (1 - p) * n2 * n3^2 * n4) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3^2 * n4 * n1),

      -d * rho3 * (1 - r23) / (p * (1 - p) * n3 * n4^2) -
        d * rho2 * (1 - r22) / (p * (1 - p) * n2 * n3 * n4^2) -
        d * (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * n2 * n3 * n4^2 * n1),

      -(1 - 2 * p) * d * rho3 * (1 - r23) / ((1 - p)^2 * p^2 * n3 * n4) -
        (1 - 2 * p) * d * rho2 * (1 - r22) / ((1 - p)^2 * p^2 * n2 * n3 * n4) -
        (1 - 2 * p) * d * (1 - rho3 - rho2) * (1 - r21) / ((1 - p)^2 * p^2 * n2 * n3 * n4 * n1)
    )
  )

  .cost.jacob <- expression(
    c(
      n4 * n3 * n2 * (p * cn1[1] + (1 - p) * cn1[2]),

      n4 * n3 * (p * cn2[1] + (1 - p) * cn2[2]) +
        n4 * n3 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      n4 * (p * cn3[1] + (1 - p) * cn3[2]) +
        n4 * n2 * (p * cn2[1] + (1 - p) * cn2[2]) +
        n4 * n2 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      cn4 +
        n3 * (p * cn3[1] + (1 - p) * cn3[2]) +
        n3 * n2 * (p * cn2[1] + (1 - p) * cn2[2]) +
        n3 * n2 * n1 * (p * cn1[1] + (1 - p) * cn1[2]),

      n4 * n3 * (cn3[1] - cn3[2]) +
        n4 * n3 * n2 * (cn2[1] - cn2[2]) +
        n4 * n3 * n2 * n1 * (cn1[1] - cn1[2])
    )
  )

  if(all(cn1 == 0) & all(cn2 == 0) & all(cn3 == 0) & is.null(p)) p <- .50

  cosa <- .cosa(order = order, interaction = interaction,
                cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cost = cost,
                constrain = constrain, round = round,
                max.power = max.power, local.solver = local.solver,
                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                g3 = g3, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2, n3 = n3, n4 = n4)
  cosa.out <- list(parms = list(dists = dists, k1 = k1, k2 = k2,
                                order = order, interaction = interaction,
                                treat.lower = treat.lower, cutoff = cutoff,
                                cn1 = cn1, cn2 = cn2, cn3 = cn3, cn4 = cn4, cost = cost,
                                constrain = constrain, round = round,
                                max.power = max.power, local.solver = local.solver,
                                power = power, es = es, alpha = alpha, two.tailed = two.tailed,
                                rho2 = rho2, rho3 = rho3, r21 = r21, r22 = r22, r23 = r23,
                                g3 = g3, p0 = p0, p = p, n0 = n0, n1 = n1, n2 = n2, n3 = n3, n4 = n4),
                   cosa = cosa)
  class(cosa.out) <- c("cosa", "bcrd4f3")
  .summary.cosa(cosa.out)
  return(invisible(cosa.out))
}
