# summarize mdes output
.summary.mdes <- function(object, ...) {

  order <- as.character(get("order", parent.frame()))
  score.info <- switch(order,
                     "0" = "Random assignment design",
                     "Regression discontinuity design")

  cat("\nMinimum detectable effect size: \n--------------------------------------- \n",
      round(object$mdes[1], 3), " ", 100 * (1 - round(object$parms$alpha, 2)),
      "% CI [", round(object$mdes[2], 3), ",", round(object$mdes[3], 3), "]",
      "\n---------------------------------------\nDegrees of freedom: ", object$df,
      "\nStandardized standard error: ", round(object$sse, 3), "\nType I error rate: ", object$parms$alpha,
      "\nType II error rate: ", round(1 - object$parms$power, 3), "\nTwo-tailed test: ", object$parms$two.tailed,
      "\n--------------------------------------- \n", score.info,
      "\n", sep = "")
}

# summarize power output
.summary.power <- function(object, ...) {

   order <- as.character(get("order", parent.frame()))
   score.info <- switch(order,
                        "0" = "Random assignment design",
                        "Regression discontinuity design")

   mlu <- .mdes(power = object$power, alpha = object$parms$alpha,
               sse = object$sse, df = object$df, two.tailed = object$parms$two.tailed)

   cat("\nStatistical power: \n--------------------------------------- \n ",
      round(object$power, 3),
      "\n--------------------------------------- \nDegrees of freedom: ", object$df,
      "\nStandardized standard error: ", round(object$sse, 3), "\nType I error rate: ", object$parms$alpha,
      "\nType II error rate: ", round(1 - object$power, 3), "\nTwo-tailed test: ", object$parms$two.tailed,
      "\n--------------------------------------- \n", score.info,
      "\n", sep = "")
}

.summary.cosa <- function(object, ...) {

  order <- as.character(get("order", parent.frame()))
  score.info <- switch(order,
                       "0" = "Random assignment design",
                       "Regression discontinuity design")

  design <- class(object)[2]
  rlevel <- as.numeric(substr(design, nchar(design), nchar(design)))
  nlevels <- as.numeric(substr(design, nchar(design) - 2, nchar(design) - 2))
  if(object$parms$round) {
    cat("\nRounded solution: \n--------------------------------------------------- \n")
  } else {
    cat("\nExact solution: \n--------------------------------------------------- \n")
  }
  print(as.data.frame(round(object$cosa, 3)), row.names=FALSE)
  if(length(object$parms$cn1) == 1) {
    object$parms$cn1 <- rep(object$parms$cn1, 2)
  }
  cat("--------------------------------------------------- \nPer unit marginal costs: \n",
      "Level 1 treatment:", object$parms$cn1[1] ,
      "\n Level 1 control:", object$parms$cn1[2], "\n")
  if(nlevels >= 2 & rlevel >= 2){
    if(length(object$parms$cn2) == 1) {
      object$parms$cn2 <- rep(object$parms$cn2, 2)
    }
    cat(" Level 2 treatment:", object$parms$cn2[1] ,
        "\n Level 2 control:", object$parms$cn2[2], "\n")
  }else if(nlevels >= 2 & rlevel < 2){
    cat(" Level 2:", object$parms$cn2[1], "\n")
  }
  if(nlevels >= 3 & rlevel >= 3){
    if(length(object$parms$cn3) == 1) {
      object$parms$cn3 <- rep(object$parms$cn3, 2)
    }
    cat(" Level 3 treatment:", object$parms$cn3[1] ,
        "\n Level 3 control:", object$parms$cn3[2], "\n")
  }else if(nlevels >= 3 & rlevel < 3){
    cat(" Level 3:", object$parms$cn3[1], "\n")
  }
  if(nlevels >= 4 & rlevel >= 4){
    if(length(object$parms$cn4) == 1) {
      object$parms$cn4 <- rep(object$parms$cn4, 2)
    }
    cat(" Level 4 treatment:", object$parms$cn4[1] ,
        "\n Level 4 control:", object$parms$cn4[2], "\n")
  }else if(nlevels >= 4 & rlevel < 4){
    cat(" Level 4:", object$parms$cn4[1], "\n")
  }
  cat("--------------------------------------------------- \nMDES = ", round(object$cosa[nlevels + 3], 3),
      " (with power = ", round(object$parms$power, 3)*100, ") \npower = ", round(object$cosa[nlevels + 6], 3),
      " (for ES = ", round(object$parms$es, 3), ") \n--------------------------------------------------- \n",
      "[]: point constrained (fixed) \n<<: bound constrained \n", score.info, "\n", sep = "")
}
