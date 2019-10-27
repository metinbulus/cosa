.error.handler <- function(x) {

  names.x <- names(x)
  if(any(!names.x %in% c("", "cost", "cn1", "cn2", "cn3", "cn4", "n0", "p0",
                         "constrain", "round", "max.power", "local.solver",
                         "score", "dists",  "k1", "k2", "rhots",
                         "order", "interaction", "treat.lower", "cutoff",
                         "df", "n1", "n2", "n3", "n4", "g1", "g2", "g3", "g4",
                         "r21","r22","r23","r24", "r2t2", "r2t3", "r2t4",
                         "rho2", "rho3", "rho4", "omega2", "omega3", "omega4",
                         "rate.tp", "rate.cc", "p", "alpha", "power", "mdes",
                         "es", "two.tailed"))) {
    stop("Unused arguments", call. = FALSE)
  }

  # exclude NULL arguments and redefine the check list
  idx.notnull <- match(names(lapply(x, is.null)[!lapply(x, is.null) == TRUE]),
                       names.x)
  parms.notnull <- x[idx.notnull]
  names.x <- names(parms.notnull)
  x <- lapply(parms.notnull, eval)


  # validity check for sample sizes
  idx.n <- intersect(c("n1","n2","n3","n4", "df"),  names.x)
  length.list.n <- length(unlist(x[idx.n]))
  length.unlist.n <- length(x[idx.n])
  if(length.list.n == length.unlist.n){
    if(any(x[idx.n] <= 0) ||
       any(lapply(x[idx.n], function(x)!is.numeric(x)) == TRUE) ||
       any(lapply(x[idx.n], length) > 1)) {
      stop("Incorrect / insufficient sample size or degrees of freedom", call.=FALSE)
    }
  }

  # validity check for number of covariates
  idx.g <- intersect(c("g1", "g2", "g3", "g4"),  names.x)
  if(length(idx.g) > 0) {
    if(!is.numeric(x[[idx.g]]) ||
       length(x[idx.g]) > 1 ||
       any(x[idx.g] < 0)) {
      stop("Incorrect number of covariates", call.=FALSE)
    }
  }

  # validity check for variance parameters, proportions, and probabilities
  idx.var <- intersect(c("r21","r22","r23","r24", "r2t2", "r2t3", "r2t4",
                         "r2m1", "r2m2", "r2m3","rhom2", "rhom3", "omegam2", "omegam3",
                         "rho2", "rho3", "rho4", "omega2", "omega3", "omega4",
                         "rate.tp", "rate.cc", "alpha", "power"),  names.x)
  if(any(lapply(x[idx.var], function(x)!is.numeric(x)) == TRUE) ||
     any(lapply(x[idx.var], length) > 1) ||
     any(x[idx.var] < 0) ||
     any(x[idx.var] > 1)) {
    stop("Incorrect value for [0, 1] or (0, 1) bounded arguments", call.=FALSE)
  }

  # validity check for R-squared value and number of covariate consistency
  idx.r2 <- intersect(c("r21", "r22", "r23", "r24", "r2t2",
                        "r2t3", "r2t4"), names.x)
  if(length(idx.g) != 0 & length(idx.r2) != 0) {
    if (any(x[idx.r2] > 0) & x[idx.g] == 0) {
      x.r2 <- x[idx.r2]
      x.g <- x[idx.g]
      err.r2 <- names(x.r2[x.r2 > 0])
      err.g <- names(x.g[x.g == 0])
      if (any(substr(err.r2, nchar(err.r2), nchar(err.r2))== substr(err.g, 2, 2))){
        warning("R-squared value for a level may not be greater than zero",
                call. = FALSE)
      }
    } else if (any(x[idx.r2] == 0) & x[idx.g] > 0) {
      x.r2 <- x[idx.r2]
      x.g <- x[idx.g]
      err.r2 <- names(x.r2[x.r2 == 0])
      err.g <- names(x.g[x.g > 0])
      if (any(substr(err.r2, nchar(err.r2), nchar(err.r2)) == substr(err.g, 2, 2))) {
        warning("R-squared value for a level may not be zero",
                call. = FALSE)
      }
    }
  }

  if("es" %in% names.x) {
    if(is.na(x$es) ||
       any(lapply(x$es, function(x)!is.numeric(x)) == TRUE) ||
       any(lapply(x$es, length) > 1) ||
       any(x$es < 0)) {
      stop("Incorrect value for effect size", call.=FALSE)
    }
    if(any(x$es > 6)) {
      warning("Extreme value for effect size (es > 6?)", call.=FALSE)
    }
  }

  if("two.tailed" %in% names.x){
    if(!is.logical(x$two.tailed) || length(x$two.tailed) > 1 ){
      stop("Non-logical value for argument 'two.tailed'", call.=FALSE)
    }
  }

  if("interaction" %in% names.x){
    if(!is.logical(x$interaction) || length(x$interaction) > 1 ){
      stop("Non-logical value for argument 'interaction'", call.=FALSE)
    }
  }

  if("treat.lower" %in% names.x){
    if(!is.logical(x$treat.lower) || length(x$treat.lower) > 1 ){
      stop("Non-logical value for argument 'treat.lower'", call.=FALSE)
    }
  }

  if("order" %in% names.x){
    if(x$order %% 1 != 0 || x$order < 0 || x$order > 8) {
      stop("'order' argument can take values:
         0 (for random assignment designs)
         1 to 8 (for regression discontinuity designs)", call. = FALSE)
    }
  }

  if("constrain" %in% names.x & !"cost" %in% names.x){
    if(x$constrain == "cost") {
      stop("Primary constraint is placed on total cost but 'cost' argument is NULL", call. = FALSE)
    }
  }

  if("p" %in% names.x){
    if(any(x$p < .01) || any(x$p > .99) || !is.numeric(x$p) || length(x$p) > 2){
      stop("Incorrect value for [.01, .99] bounded argument 'p'", call. = FALSE)
      }
  }

  if(any(c("rate.tp", "rate.cc") %in% names.x)){
    ifelse(!"rate.cc" %in% names.x, rate.cc <- 0, rate.cc <- x$rate.cc)
    ifelse(!"rate.tp" %in% names.x, rate.tp <- 1, rate.tp <- x$rate.tp)
    if(rate.tp != 1 | rate.cc != 0) {
      message(cat("\nLocal average treatment effect (LATE)",
                  "\nTreatment group participant rate =", rate.tp,
                  "\nControl group cross-over rate =", rate.cc, "\n"))
    }
    if(rate.cc >= rate.tp) {
      stop("'rate.cc' >= 'rate.tp' ?!", call. = FALSE)
    }
    if(rate.cc >= .98) {
      stop("'rate.cc' =~ 1 ?!", call. = FALSE)
    }
  }

  # if(any(c("k1", "k2", "dists") %in% names.x)) {
  #   warning("'k1', 'k2', 'dists' arguments will be removed in the next release,
  #           use inspect.score() function instead", call. = FALSE)
  # }

} #.error.handler()

