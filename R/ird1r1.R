mdes.ird1r1 <- function(score = NULL, dists = "normal", k1 = -6, k2 = 6,
                        order = 1, interaction = FALSE, treat.lower = TRUE, cutoff = 0, p = NULL,
                        power = .80, alpha = .05, two.tailed = TRUE,
                        df = n1 - g1 - order * (1 + interaction) - 2,
                        r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, n1) {

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

  sse <- (1/(rate.tp - rate.cc)) * sqrt(d * (1 - r21) / (p * (1 - p) * n1))

  mdes <- .mdes(power, alpha, sse, df, two.tailed)
  colnames(mdes) <- c("mdes", paste0(100 * (1 - round(alpha, 2)), "%lcl"),
                      paste0(100 * (1 - round(alpha, 2)), "%ucl"))
  mdes.out <- list(parms = list(dists = dists, k1 = k1, k2 = k2,
                                order = order, interaction = interaction,
                                treat.lower = treat.lower, p = p, cutoff = cutoff,
                                power = power, alpha = alpha, two.tailed = two.tailed,
                                r21 = r21, g1 = g1, rate.tp = rate.tp, rate.cc = rate.cc, n1 = n1),
                   df = df,
                   sse = sse,
                   mdes = mdes)

  class(mdes.out) <- c("mdes", "ird1r1")
  .summary.mdes(mdes.out)
  return(invisible(mdes.out))
}
mdes.ird <- mdes.ird1r1

power.ird1r1 <- function(score = NULL, dists = "normal", k1 = -6, k2 = 6,
                         order = 1, interaction = FALSE, treat.lower = TRUE, cutoff = 0, p = NULL,
                         es = .25, alpha = .05, two.tailed = TRUE,
                         df = n1 - g1 - order * (1 + interaction) - 2,
                         r21 = 0, g1 = 0, rate.tp = 1, rate.cc = 0, n1) {

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

  sse <- (1/(rate.tp - rate.cc)) * sqrt(d * (1 - r21) / (p * (1 - p) * n1))

  power <- .power(es, alpha, sse, df, two.tailed)
  power.out <-  list(parms = list(dists = dists, k1 = k1, k2 = k2,
                                  order = order, interaction = interaction,
                                  treat.lower = treat.lower, p = p, cutoff = cutoff,
                                  es = es, alpha = alpha, two.tailed = two.tailed,
                                  r21 = r21, g1 = g1, rate.tp = rate.tp, rate.cc = rate.cc, p = p, n1 = n1),
                     df = df,
                     sse = sse,
                     power = power)
  class(power.out) <- c("power", "ird1r1")
  .summary.power(power.out)
  return(invisible(power.out))
}
power.ird <- power.ird1r1
