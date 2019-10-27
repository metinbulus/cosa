# minimum detectable effect size
.mdes <- function(power, alpha, sse, df, two.tailed){
  t1 <- ifelse(two.tailed == TRUE, abs(qt(alpha / 2, df)), abs(qt(alpha, df)))
  t2 <- abs(qt(power, df))
  m <- ifelse(power >= 0.5, t1 + t2, t1 - t2)
  mdes <- m * sse
  lcl <- mdes * (1 - t1 / m)
  ucl <- mdes * (1 + t1 / m)
  mlu <- cbind(mdes, lcl, ucl)
  return(mlu)
}

# statistical power
.power <- function(es, alpha, sse, df, two.tailed){
  lambda <- es/sse
  power <- ifelse(two.tailed == FALSE,
                  1 - pt(qt(alpha, df, lower.tail = FALSE), df, lambda),
                  1 - pt(qt(alpha / 2, df, lower.tail = FALSE), df, lambda) +
                    pt(-qt(alpha / 2, df, lower.tail = FALSE), df, lambda))
  return(power)
}

# constrained optimal sample allocation
.cosa <- function(order = 2, interaction = FALSE,
                  cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0, cost = NULL,
                  n1 = NULL, n2 = NULL, n3 = NULL, n4 = NULL, p = NULL, n0, p0,
                  constrain, local.solver = c("LBFGS", "SLSQP", "MMA", "COBYLA"),
                  power = .80, es = .25, alpha = .05, two.tailed = TRUE,
                  rho2, rho3, rho4, omega2, omega3, omega4,
                  g1 = 0, g2 = 0, g3 = 0, g4 = 0,
                  r21 = 0, r22 = 0, r23 = 0, r24 = 0,
                  r2t2 = 0, r2t3 = 0, r2t4 = 0, round = TRUE,
                  jacobian = TRUE, max.power = FALSE,
                  localtol = 1e-8, fabstol = 1e-8, maxeval = 500) {

  .df <- get(".df", parent.frame())
  .sse <- get(".sse", parent.frame())
  .cost <- get(".cost", parent.frame())
  .cost.jacob <- get(".cost.jacob", parent.frame())
  .var.jacob <- get(".var.jacob", parent.frame())
  d <- get("d", parent.frame())
  lb <- get("lb", parent.frame())
  fun <- get("fun", parent.frame())

  fun_parsed <- scan(text = fun, what = "character", sep=".", quiet = TRUE)
  rlevel <- as.numeric(substr(fun, nchar(fun), nchar(fun)))
  nlevels <- as.numeric(substr(fun, nchar(fun)-2, nchar(fun)-2))

  # equality constraints on cost
  .eq.cost <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(eval(.cost) - cost)
  }

  # equality constraints on power
  .eq.power <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(abs(.power(es = es,  alpha = alpha,
                      sse = eval(.sse), df = eval(.df),
                      two.tailed = two.tailed) - power))
  }

  # equality constraints on ES
  .eq.es <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(abs(.mdes(power = power, alpha = alpha,
                     sse = eval(.sse), df = eval(.df),
                     two.tailed = two.tailed)[1] - es))
  }

  #  MDES and power
  .mlu.pwr <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(c(.mdes(power = power, alpha = alpha,
                   sse = eval(.sse), df = eval(.df),
                   two.tailed = two.tailed),
             .power(es = es, alpha = alpha,
                    sse = eval(.sse), df = eval(.df),
                    two.tailed = two.tailed)))
  }

  # minimize cost given power or mdes
  .min.cost <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(eval(.cost))
  }

  # minimize treatment variance given cost
  .min.var <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(eval(.sse)^2)
  }

  # maximize power rate
  .max.pwr <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(abs(.power(es = es,  alpha = alpha,
                      sse = eval(.sse), df = eval(.df),
                      two.tailed = two.tailed) - 1))
  }

  # starting values
  .start <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(.power(es = es,  alpha = alpha,
                  sse = eval(.sse), df = eval(.df),
                  two.tailed = two.tailed) - .95)
  }

  # jacobian for cost
  .cost.jac <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(
      eval(.cost.jacob)
    )
  }

  # jacobian for variance
  .var.jac <- function(ss){
    n1 <- ss[1]
    if(nlevels >= 2){
      n2 <- ss[2]
    }
    if(nlevels >= 3){
      n3 <- ss[3]
    }
    if(nlevels == 4){
      n4 <- ss[4]
    }
    p <- ss[nlevels + 1]
    return(
      eval(.var.jacob)
    )
  }

  cost.list <- list(cn1, cn2, cn3, cn4)
  cost.names <- c("cn1", "cn2", "cn3", "cn4")
  cost.lengths <- unlist(lapply(cost.list, length))
  if(rlevel < nlevels & any(cost.lengths[(rlevel + 1):nlevels] == 2)) {
    stop("Unequal cost applies to levels at or below randomization (or discontinuity) level", call.=FALSE)
  }

  if(any(cost.lengths[1:rlevel] > 2) || !is.numeric(unlist(cost.list))) {
    stop("Incorrect value for marginal cost arguments", call.=FALSE)
  }

  if(!is.null(cost)) {
    if(length(cost) > 1 || !is.numeric(cost) || cost < 0) {
      stop("Incorrect value for argument 'cost'", call.=FALSE)
    }
  }

  if(max.power == TRUE) {
    if(tolower(constrain) %in% c("es", "power")) {
      stop("Power cannot be maximized with primary constraint on effect size or power",
           call. = FALSE)
    }
    .min.var <- .max.pwr
    jacobian <- FALSE
  }
  fn.min <- switch(tolower(constrain),
                   "power" = .min.cost,
                   "es" = .min.cost,
                   "cost" = .min.var,
                   stop("Incorrect constraint",
                        call. = FALSE))
  fn.constr <- switch(tolower(constrain),
                      "power" = .eq.power,
                      "es" = .eq.es,
                      "cost" = .eq.cost)
  if(jacobian) {
    fn.min.jacob <- switch(tolower(constrain),
                           "power" = .cost.jac,
                           "es" = .cost.jac,
                           "cost" = .var.jac)
  } else {
    fn.min.jacob <- NULL
  }
  fn.constr.jacob <- switch(tolower(constrain),
                            "power" = NULL,
                            "es" = NULL,
                            "cost" = .cost.jac)

  if(rlevel >= 1 & length(cn1) == 1) {
    cn1 <- c(cn1, cn1)
  }
  if(rlevel >= 2 & length(cn2) == 1) {
    cn2 <- c(cn2, cn2)
  }
  if(rlevel >= 3 & length(cn3) == 1) {
    cn3 <- c(cn3, cn3)
  }
  if(rlevel == 4 & length(cn4) == 1) {
    cn4 <- c(cn4, cn4)
  }

  if(p0 < .01 || p0 > .99 || !is.numeric(p0) || length(p0) != 1) {
    stop("Incorrect value for argment 'p0'",  call. = FALSE)
  }
  if(any(n0 < 0) || !is.numeric(n0) || length(n0) != nlevels) {
    stop("Incorrect value for argument 'n0'",  call. = FALSE)
  }

  p0  <- ifelse(!is.null(p), mean(p), p0)
  plb <- ifelse(!is.null(p), min(p), .10)
  pub <- ifelse(!is.null(p), max(p), .90)
  n10  <- ifelse(!is.null(n1), mean(n1), n0[1])
  n1lb <- ifelse(!is.null(n1), min(n1), lb[1])
  n1ub <- ifelse(!is.null(n1), max(n1), 1e+6)
  if(nlevels == 1){
    ss0  <- c(n10,p0)
    sslb <- c(n1lb,plb)
    ssub <- c(n1ub,pub)
  }else if(nlevels == 2){
    n20  <- ifelse(!is.null(n2), mean(n2), n0[2])
    n2lb <- ifelse(!is.null(n2), min(n2), lb[2])
    n2ub <- ifelse(!is.null(n2), max(n2), 1e+6)
    ss0  <- c(n10,n20,p0)
    sslb <- c(n1lb,n2lb,plb)
    ssub <- c(n1ub,n2ub,pub)
  }else if(nlevels == 3){
    n20  <- ifelse(!is.null(n2), mean(n2), n0[2])
    n2lb <- ifelse(!is.null(n2), min(n2), lb[2])
    n2ub <- ifelse(!is.null(n2), max(n2), 1e+6)
    n30  <- ifelse(!is.null(n3), mean(n3), n0[3])
    n3lb <- ifelse(!is.null(n3), min(n3), lb[3])
    n3ub <- ifelse(!is.null(n3), max(n3), 1e+6)
    ss0  <- c(n10,n20,n30,p0)
    sslb <- c(n1lb,n2lb,n3lb,plb)
    ssub <- c(n1ub,n2ub,n3ub,pub)
  }else if(nlevels == 4){
    n20  <- ifelse(!is.null(n2), mean(n2), n0[2])
    n2lb <- ifelse(!is.null(n2), min(n2), lb[2])
    n2ub <- ifelse(!is.null(n2), max(n2), 1e+6)
    n30  <- ifelse(!is.null(n3), mean(n3), n0[3])
    n3lb <- ifelse(!is.null(n3), min(n3), lb[3])
    n3ub <- ifelse(!is.null(n3), max(n3), 1e+6)
    n40  <- ifelse(!is.null(n4), mean(n4), n0[4])
    n4lb <- ifelse(!is.null(n4), min(n4), lb[4])
    n4ub <- ifelse(!is.null(n4), max(n4), 1e+6)
    ss0  <- c(n10,n20,n30,n40,p0)
    sslb <- c(n1lb,n2lb,n3lb,n4lb,plb)
    ssub <- c(n1ub,n2ub,n3ub,n4ub,pub)
  }

  if(any(ss0 < sslb) || any(ss0 > ssub)) {
    stop("Starting values are out of bounds",  call. = FALSE)
  }

  local.solver <- toupper(local.solver)
  if(constrain %in% c("es","power") || max.power == TRUE) {
    solver.power.es <- c("SLSQP", "LBFGS", "MMA", "COBYLA")
    local.solver <- local.solver[order(match(local.solver, solver.power.es))]
  }

  i <- 1; conv <- FALSE
  while(i <= length(local.solver) & conv == FALSE){

    if(!local.solver[i] %in% c("LBFGS", "SLSQP",  "MMA", "COBYLA")) {
      stop("Incorrect value for argument 'local.solver'",  call. = FALSE)
    }

    if(local.solver[i] %in% c("MMA", "COBYLA")) {
      warning("Possibility of a local solution",  call. = FALSE)
    }

    if(constrain == "cost") {
      par0 <- try(nloptr::auglag(x0 = ss0, fn = .min.cost, gr = .cost.jac, heq = .start,
                          localsolver = local.solver[i], localtol = localtol,
                          lower = sslb, upper = ssub,
                          control = list(ftol_abs = fabstol, maxeval = maxeval))$par)
      if(!inherits(par0, "try-error")) {
        ss00 <- par0
      } else {
        ss00 <- ss0
      }
    } else {
      ss00 <- ss0
    }

    nlopt.ss <- nloptr::auglag(x0 = ss00, fn = fn.min, heq = fn.constr,
                       gr = fn.min.jacob, heqjac = fn.constr.jacob,
                       localsolver = local.solver[i], localtol = localtol,
                       lower = sslb, upper = ssub,
                       control = list(ftol_abs = fabstol, maxeval = maxeval))

    if(round) {
      if(length(intersect(ss0, round(nlopt.ss$par, 1))) < nlevels - 1){
        if(nlopt.ss$par[nlevels] %% 1 != 0) {
          sslb.round <- sslb
          ssub.round <- ssub
          ss00.round <- ss00
          ss00.round[nlevels] <- sslb.round[nlevels] <- ssub.round[nlevels] <- round(nlopt.ss$par[nlevels])
          nlopt.ss <- nloptr::auglag(x0 = ss00.round, fn = fn.min, heq = fn.constr,
                             gr = fn.min.jacob, heqjac = fn.constr.jacob,
                             localsolver = local.solver[i], localtol = localtol,
                             lower = sslb.round, upper = ssub.round,
                             control = list(ftol_abs = fabstol, maxeval = maxeval))
        }
      }
    }

    if(nlopt.ss$convergence < 0 | all(nlopt.ss$par == ss0) | any(nlopt.ss$par <= 0)) {
      conv <- FALSE
      cat("Solution is not feasible with ", local.solver[i], "\n")
    } else {
      conv <- TRUE
      cat("Solution converged with", local.solver[i], "\n")
      ss1 <- nlopt.ss$par
    }
    i <- i+1

  }

  if(nlopt.ss$convergence < 0 | all(nlopt.ss$par == ss0) | any(nlopt.ss$par <= 0)) {
    stop("Solution is not feasible. Change default settings", call.=FALSE)
  }

  col.names <- c(c(ifelse(!is.null(n1), ifelse(length(n1) >= 2, "<n1<", "[n1]"),"n1"),
                   ifelse(!is.null(n2), ifelse(length(n2) >= 2, "<n2<", "[n2]"),"n2"),
                   ifelse(!is.null(n3), ifelse(length(n3) >= 2, "<n3<", "[n3]"),"n3"),
                   ifelse(!is.null(n4), ifelse(length(n4) >= 2, "<n4<", "[n4]"),"n4"))[1:nlevels],
                 ifelse(!is.null(p), ifelse(length(p) == 2, "<p<", "[p]"),"p"),
                 ifelse(constrain == "cost","[cost]","cost"),
                 ifelse(constrain == "es","[mdes]","mdes"),
                 paste0(100 * round((1 - alpha), 2), "%lcl"),
                 paste0(100 * round((1 - alpha), 2), "%ucl"),
                 ifelse(constrain == "power","[power]","power"))

  if(round) {
    ss1[nlevels] <- round(ss1[nlevels])
    if(nlevels==2){
      ss1[1] <- round(prod(ss1[1:2]))/ss1[2]
    }else if(nlevels==3){
      ss1[2] <- round(prod(ss1[2:3]))/ss1[3]
      ss1[1] <- round(prod(ss1[1:3]))/prod(ss1[2:3])
    }else if(nlevels==4){
      ss1[3] <- round(prod(ss1[3:4]))/ss1[4]
      ss1[2] <- round(prod(ss1[2:4]))/prod(ss1[3:4])
      ss1[1] <- round(prod(ss1[1:4]))/prod(ss1[2:4])
    }
    if(nlevels==rlevel){
      ss1[nlevels+1] <- round(ss1[nlevels+1]*ss1[rlevel]) / ss1[rlevel]
    }else{
      ss1[nlevels+1] <- round(ss1[nlevels+1]*round(prod(ss1[rlevel:nlevels]))) / round(prod(ss1[rlevel:nlevels]))
    }

    est.mlu.power <- .mlu.pwr(ss1)
    cosa <- cbind(t(ss1), .min.cost(ss1), est.mlu.power[1], est.mlu.power[2], est.mlu.power[3], est.mlu.power[4])
    colnames(cosa) <-  col.names
  } else {
    est.mlu.power <- .mlu.pwr(ss1)
    cosa <- cbind(t(ss1), .min.cost(ss1), est.mlu.power[1], est.mlu.power[2], est.mlu.power[3], est.mlu.power[4])
    colnames(cosa) <-  col.names
  }

  return(invisible(cosa))
}

# object conversion
.cosa2mdes <- function(x, score = NULL){
  if(inherits(x, "cosa")){
    design <- class(x)[2]
    nlevels <- as.numeric(substr(design, nchar(design) - 2, nchar(design) - 2))
    fun <- paste("mdes", class(x)[2], sep = ".")
    parms <- x$parms[intersect(names(x$parms), names(formals(fun)))]
    if(x$parms$dists == "empirical") {
      ifelse(is.null(score),
             stop("Score object is missing", call. = FALSE),
             parms$score <- score)
    }
    parms$p <- x$cosa[nlevels + 1]
    parms$n1 <- x$cosa[1]
    if(nlevels >= 2) {
      parms$n2 <- x$cosa[2]
    }
    if(nlevels >= 3) {
      parms$n3 <- x$cosa[3]
    }
    if(nlevels == 4) {
      parms$n4 <- x$cosa[4]
    }
    return(invisible(do.call(fun, parms)))
  }else{
    stop("x should be an object returned from BCOSSA functions", call.=FALSE)
  }
}

.cosa2power <- function(x, score = NULL){
  if(inherits(x, "cosa")){
    design <- class(x)[2]
    nlevels <- as.numeric(substr(design, nchar(design) - 2, nchar(design) - 2))
    fun <- paste("power", class(x)[2], sep = ".")
    parms <- x$parms[intersect(names(x$parms), names(formals(fun)))]
    if(x$parms$dists == "empirical") {
      ifelse(is.null(score),
             stop("Score object is missing", call. = FALSE),
             parms$score <- score)
    }
    parms$p <- x$cosa[nlevels + 1]
    parms$n1 <- x$cosa[1]
    if(nlevels >= 2) {
      parms$n2 <- x$cosa[2]
    }
    if(nlevels >= 3) {
      parms$n3 <- x$cosa[3]
    }
    if(nlevels == 4) {
      parms$n4 <- x$cosa[4]
    }
    return(invisible(do.call(fun, parms)))
  }else{
    stop("x should be an object returned from BCOSSA functions", call.=FALSE)
  }
}

.power2mdes <- function(x, score = NULL){
  if(inherits(x, "power")){
    fun <- paste("mdes", class(x)[2], sep = ".")
    parms <- x$parms[intersect(names(x$parms), names(formals(fun)))]
    if(x$parms$dists == "empirical") {
      ifelse(is.null(score),
             stop("Score object is missing", call. = FALSE),
             parms$score <- score)
    }
    parms$power <- x$power
    return(invisible(do.call(fun, parms)))
  }else{
    stop("x should be an object returned from statistical power functions", call.=FALSE)
  }
}

.mdes2power <- function(x, score = NULL){
  if(inherits(x, "mdes")){
    fun <- paste("power", class(x)[2], sep = ".")
    parms <- x$parms[intersect(names(x$parms), names(formals(fun)))]
    if(x$parms$dists == "empirical") {
      ifelse(is.null(score),
             stop("Score object is missing", call. = FALSE),
             parms$score <- score)
    }
    parms$es<- x$mdes[1]
    return(invisible(do.call(fun, parms)))
  }else{
    stop("x should be an object returned from MDES functions", call.=FALSE)
  }
}
