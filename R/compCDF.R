compCDF <- function(obj, scores = c("wmw", "logrank1","logrank2"),
                    B = 199, type = c("naive", "permutation"),
                    units = NULL, upperl, lowerl, display = TRUE){

  if(!inherits(obj, "drcte") | inherits(obj, "drcteList")) {
     stop("Method works only with 'drcte' objects in case of simultaneous fitting (i.e. it does not work when the 'separate = T' option has been used)")
  }
  ncol <- length(obj$data[1,])
  group <- as.character(rep(obj$data[,ncol - 1], obj$data[,ncol - 3]))
  nlev <- length(levels(factor(group)))
  if(nlev < 2 ) {
     stop("Nothing to compare")
  }

  # Manage unit variable
  unitName <- deparse(substitute(units))
  if(unitName != "NULL"){
    test <- obj$origData[[unitName]]
    if(is.null(test)) units <- units else units <- test
  }

  if(obj$fit$method == "NPMLE"){
    scores <- match.arg(scores)
    compCDFnp(obj = obj, scores = scores, B = B, units = units,
              display = display)
  } else if(obj$fit$method == "KDE"){
    compCDFkde(obj = obj, B = B, units = units, display = display)
  } else {
    type <- match.arg(type)
    compCDFpar(obj = obj, B = B, units = units, type = type,
               display = display)
  }
}

compCDFpar <- function(obj, B, units, type = type,
                       display){
  # Fit a null model
  args <- getCall(obj)
  call <- args[names(args) != "curveid" & names(args) != "pmodels" &
        names(args) != "upperl" & names(args) != "lowerl"]

  obj2 <- eval(call, parent.frame())

  # Edited 6/3/2023. For certain models, there is a form of
  # control to avoid that 'd' goes in an unreasonable range
  fctName <- deparse(substitute(obj$fct$name))
  if(grepl("L.3(",  fctName, fixed=TRUE) |
         grepl("LN.3(", fctName, fixed=TRUE) |
         grepl("W1.3(", fctName, fixed=TRUE) |
         grepl("W2.3(", fctName, fixed=TRUE) |
         grepl("G.3(",  fctName, fixed=TRUE) |
         grepl("loglogistic",  fctName, fixed=TRUE) ){
    d <- coef(obj2)[2]
    if(d > 1 | d < 0) obj2 <- update(obj2, upperl = c(NA, 1, NA))
  }

  # Naive lrt (wrong in two ways!)
  LRT <- as.numeric(2 * (logLik(obj) - logLik(obj2)))
  n.df <- length(coef(obj)) - length(coef(obj2))
  pval <- pchisq(LRT, n.df, lower.tail = F)

  if(type == "naive"){
    if(display){
      cat("\n")
      cat("\n")
      TEST <- "Likelihood ratio test"
      null.phrase <- "NULL: time-to-event curves are equal"
      alt.phrase <- "ALTERNATIVE: time-to-event curves are not equal"
      cat(TEST)
      cat("\n")
      cat(null.phrase)
      cat("\n")
      cat("\n")
      cat(paste("Observed LR value: ", round(LRT, 4)))
      cat("\n")
      cat(paste("Degrees of freedom: ", n.df))
      cat("\n")
      cat(paste("P-value: ", ifelse(pval < 0.00001, format(pval, scientific = T), pval)))
      cat("\n")
    }
    pout <- list()
    pout$method <- "LRT"
    pout$scores <- NULL
    pout$val0 <- LRT
    pout$vali <- NULL
    pout$pval <- pval
    pout$pvalb <- NULL
    pout$U <- NULL
    pout$N <- NULL

  } else {
    # permutation approach
    LRTb <- rep(NA, B)
    if(is.null(units)){

      # Individual based permutation
      ncol <- length(obj$data[1,])
      group <- as.character(rep(obj$data[,ncol - 1], obj$data[,ncol - 3]))
      L <- rep(obj$data[, 1], obj$data[,ncol - 3])
      R <- rep(obj$data[, 2], obj$data[,ncol - 3])
      df <- data.frame(L = L, R = R, group = group, count = rep(1, length(L)))
      message("Permuting individuals")
      for(b in 1:B){
        # If sampling individuals, the fit may not succeed!
        txtMes <- paste(round(b/B*100, 0), "%\r", sep ="")
        message(txtMes, appendLF = F)
        newList <- sample(group, replace = FALSE)
        df$group <- newList
        objNew <- try(update(obj, formula = count ~ L + R,
                         curveid = group, data = df), silent = T)
        if(any(class(objNew) == "try-error")){
          # Update 18/08/2022. In some cases, a parametric model
          # cannot be fit to the permuted sample
          message(" \r", appendLF = F)
          stop("\r Permutation based inference was not successful. See the documentation for more detail")
        }
        # obj2New <- update(objNew, curveid = NULL)
        # obj2New <- update(obj2, formula = count ~ L + R, curveid = NULL, data = df)
        args <- getCall(objNew)
        call <- args[names(args) != "curveid" & names(args) != "pmodels" &
        names(args) != "upperl" & names(args) != "lowerl"]
        obj2New <- eval(call)
        LRTb[b] <- as.numeric(2 * (logLik(objNew) - logLik(obj2New)))
        }
    } else {
      # Group-based permutation
      cluster <- units
      df <- data.frame(id = 1:length(cluster),
                   obj$data,
                   cluster = cluster)
      message("Permuting groups")
        for(b in 1:B){
          txtMes <- paste(round(b/B*100, 0), "%\r", sep ="")
          message(txtMes, appendLF = F)
          pr <- resample.cens(df, cluster = df$cluster, replace = c(F, F))
          pr$group <- obj$data[[5]]
          pr <- pr[order(pr$cluster, pr$timeBef), ]
          newList <- as.character(pr$group)
          objNew <- try(update(obj,
                   curveid = group, data = pr), silent = T)
          if(any(class(objNew) == "try-error")){
            # Update 18/08/2022. In some cases, a parametric model
            # cannot be fit to the permuted sample
            message("\r", appendLF = F)
            stop("\r Permutation based inference was not successful. See the documentation for more detail")
            }
          # obj2New <- update(objNew, curveid = NULL)
          args <- getCall(objNew)
          call <- args[names(args) != "curveid" & names(args) != "pmodels" &
                       names(args) != "upperl" & names(args) != "lowerl"]
          obj2New <- eval(call)
          LRTb[b] <- as.numeric(2 * (logLik(objNew) - logLik(obj2New)))
        }
    }
    pvalb <- (sum(LRTb > LRT) + 1)/(B + 1)
    if(display){
      cat("\r")
      # cat("\n")
      TEST <- "Likelihood ratio test (permutation based)"
      null.phrase <- "NULL: time-to-event curves are equal"
      alt.phrase <- "ALTERNATIVE: time-to-event curves are not equal"
      cat(TEST)
      cat("\n")
      cat(null.phrase)
      cat("\n")
      cat("\n")
      cat(paste("Observed LR value: ", round(LRT, 4)))
      cat("\n")
      cat(paste("Degrees of freedom: ", n.df))
      cat("\n")
      cat(paste("Naive P-value: ", ifelse(pval < 0.00001, format(pval, scientific = T), pval)))
      cat("\n")
      cat(paste("Permutation P-value (B = ", B,"): ", round(pvalb,6), sep = ""))
      cat("\n")
    }
    pout <- list()
    pout$method <- "LRT"
    pout$scores <- NULL
    pout$val0 <- LRT
    pout$vali <- LRTb
    pout$pval <- pval
    pout$pvalb <- pvalb
    pout$U <- NULL
    pout$N <- NULL
  }
  return(invisible(pout))
}

# compCDFkde.old <- function(obj,
#                        alternative= c("two.sided", "less", "greater"),
#                        B = 50){
#
#   obj2 <- update(obj, curveid = NULL) # Fitting cumulato
#
#   # Recuperare oggetti rilevanti
#   ug <- obj$dataList$names$rNames
#   recObj <- obj2$ICfit$naiveStart
#   recObj2 <- obj2$ICfit$naiveEnd
#   t <- (recObj$time + recObj2$time)/2
#   y <- unique(c(recObj$time, recObj2$time))
#   wpooled <- recObj$pdf
#   hpooled <- as.numeric(obj2$coefficients)
#   hvals <- as.numeric(obj$coefficients)
#   ntrt <- length(hvals)
#
#   ni <- tapply(obj$data[[3]], obj$data[[5]], sum)
#
#   calcDstat <- function(obj, obj2){
#     # Calculate the Cramer-von Mises distance
#
#     # bandwidths
#     hpooled <- as.numeric(obj2$coefficients)
#     hvals <- as.numeric(obj$coefficients)
#
#     t <- (obj2$ICfit$naiveStart$time + obj2$ICfit$naiveEnd$time)/2
#     wpooled <- obj2$ICfit$naiveStart$pdf
#     lim2 <- t[length(t)] # max
#     lim1 <- t[1] # min
#     limInf <- lim1 + hpooled * stats::qnorm(0.99)
#     limSup <- lim2 + hpooled * stats::qnorm(0.01)
#     Fhi <- obj$curve[[1]]
#     Fh <- obj2$curve[[1]][[1]]
#     fh <- Vectorize(function(x) 1 / hpooled * sum(wpooled * stats::dnorm((x-t)/hpooled)))
#     ntrt <- length(hvals)
#     Dval <- rep(NA, ntrt)
#
#     for(i in 1:ntrt){
#       Dval[i] <- stats::integrate(function(x)(Fhi[[i]](x)-Fh(x))^2*fh(x), limInf, limSup)$value
#       Dval[i] <- Dval[i] * ni[i]
#     }
#     obsD <- sum(Dval)/ntrt
#     return(list(Dvals = Dval, obsD=obsD))
#   }
#
#   Dlist <- calcDstat(obj, obj2)
#   D <- Dlist$obsD
#
#   Db <- rep(NA, B)
#
#   for(b in 1:B){
#     # Get a resample
#     # cat("Resampling: \t")
#     cat("\r", b)
#     newCounts <- c()
#     newId <- c()
#     for(i in 1:ntrt){
#       .tmp1 <- sample(t, size = ni[i], prob = wpooled, replace=TRUE)
#       .tmp2 <- hpooled * stats::rnorm(ni[i])
#       # xbi <- sort(sample(t, size = ni[i], prob = wpooled, replace=TRUE) + hpooled*stats::rnorm(ni[i]))
#       xbi <- sort(.tmp1 + .tmp2)
#       cbi <- table(cut(xbi, breaks = y))
#       wbi <- cbi/ni[i]
#       tti <- t
#       # wbi2 <- binnednp:::calcw_cpp(xbi, y)
#       newCounts <- c(newCounts, cbi)
#     }
#     df <- data.frame(start = rep(recObj$time, ntrt),
#                      end = rep(recObj2$time, ntrt),
#                      count = newCounts,
#                      Id = rep(1:ntrt, each = length(recObj$time)))
#     obj <- drmte(count ~ start + end, fct = KDE(),
#                  curveid = Id, data = df) # Fitting separato
#     obj2 <- update(obj, curveid = NULL) # Fitting cumulato
#     Db[b] <- calcDstat(obj, obj2)$obsD
#   }
#   pval <- mean(Db > D)
#
#   # Try a permutation approach
#
#   ## Describes the results
#   cat("\n")
#   TEST<-"Bootstrap test based on a Cramer-von-Mises type distance (Barreiro-Ures et al., 2019)"
#   null.phrase <- "NULL: time-to-event curves are equal"
#   alt.phrase <- "ALTERNATIVE: time-to-event curves are not equal"
#
#   tabRes <- data.frame("level" = as.character(ug), n = ni, D = Dlist$Dvals)
#   cat(TEST)
#   cat("\n")
#   cat(null.phrase)
#   cat("\n")
#   cat("\n")
#   print(tabRes)
#   cat("\n")
#   cat(paste("Observed D value = ", round(Dlist$obsD, 4)))
#   cat("\n")
#   cat(paste("P value = ", round(pval, 6)))
#   cat("\n")
#
#
#   # pout$data.name <- paste("{",L.name,",",R.name,"}"," by ",group.name,sep="")
#     pout <- list()
#     pout$method <- "Kramer - Von Mises distance"
#     pout$scores <- NULL
#     pout$val0 <- Dlist$obsD
#     pout$vali <- Db
#     pout$pval <- NULL
#     pout$pvalb <- pval
#     pout$U <- NULL
#     pout$N <- NULL
#   return(pout)
# }

compCDFkde <- function(obj,
                       alternative = c("two.sided", "less", "greater"),
                       B = 50, units, display){

  obj2 <- update(obj, curveid = NULL) # Fitting cumulato

  # Recuperare oggetti rilevanti
  ug <- obj$dataList$names$rNames
  recObj <- obj2$ICfit$naiveStart
  recObj2 <- obj2$ICfit$naiveEnd
  t <- (recObj$time + recObj2$time)/2
  y <- unique(c(recObj$time, recObj2$time))
  wpooled <- recObj$pdf
  hpooled <- as.numeric(obj2$coefficients)
  hvals <- as.numeric(obj$coefficients)
  ntrt <- length(hvals)
  ni <- tapply(obj$data[[3]], obj$data[[5]], sum)

  calcDstat <- function(obj, obj2){
    # Calculate the Cramer-von Mises distance

    # bandwidths
    hpooled <- as.numeric(obj2$coefficients)
    hvals <- as.numeric(obj$coefficients)

    t <- (obj2$ICfit$naiveStart$time + obj2$ICfit$naiveEnd$time)/2
    wpooled <- obj2$ICfit$naiveStart$pdf
    lim2 <- t[length(t)] # max
    lim1 <- t[1] # min
    limInf <- lim1 + hpooled * stats::qnorm(0.99)
    limSup <- lim2 + hpooled * stats::qnorm(0.01)
    Fhi <- obj$curve[[1]]
    Fh <- obj2$curve[[1]][[1]]
    fh <- Vectorize(function(x) 1 / hpooled * sum(wpooled * stats::dnorm((x-t)/hpooled)))
    ntrt <- length(hvals)
    Dval <- rep(NA, ntrt)

    for(i in 1:ntrt){
      Dval[i] <- stats::integrate(function(x)(Fhi[[i]](x)-Fh(x))^2*fh(x), limInf, limSup)$value
      Dval[i] <- Dval[i] * ni[i]
    }
    obsD <- sum(Dval)/ntrt
    return(list(Dvals = Dval, obsD=obsD))
  }

  Dlist <- calcDstat(obj, obj2)
  D <- Dlist$obsD

  Db <- rep(NA, B)
  # Try a permutation approach

  if(is.null(units)){
    message("Permuting individuals")
    for(b in 1:B){
    # Get a permutation sample (individual based)
    # cat("Resampling: \t")
    txtMes <- paste(round(b/B*100, 0), "%\r", sep ="")
    message(txtMes, appendLF = F)
    ncol <- length(obj$data[1,])
    group <- as.character(rep(obj$data[,ncol - 1], obj$data[,ncol - 3]))
    L <- rep(obj$data[, 1], obj$data[,ncol - 3])
    R <- rep(obj$data[, 2], obj$data[,ncol - 3])
    newGroup <- sample(group, replace = FALSE)
    df <- data.frame(L = L, R = R, group = newGroup, count = rep(1, length(L)))
    df <- dplyr::arrange(group_te(df), group)
    objNew <- try(drmte(count ~ L + R, fct = KDE(),
                 curveid = group, data = df), silent = T) # Fitting separato
    if(any(class(objNew) == "try-error")){
      # Update 18/08/2022. In some cases, a KDE model
      # cannot be fit to the permuted sample
      message("\r", appendLF = F)
      stop("\r Permutation based inference was not successful. See the documentation for more detail")
      }
    obj2New <- update(objNew, curveid = NULL) # Fitting cumulato
    Db[b] <- calcDstat(objNew, obj2New)$obsD
  }
  # pval <- mean(Db > D)
  } else {
    # Group-based permutation
    cluster <- units
    df <- data.frame(id = 1:length(cluster),
                   obj$data,
                   cluster = cluster)
    message("Permuting groups")
    for(b in 1:B){
      txtMes <- paste(round(b/B*100, 0), "%\r", sep ="")
      message(txtMes, appendLF = F)
      pr <- resample.cens(df, cluster = df$cluster, replace = c(F, F))
      pr$group <- obj$data[[5]]
      pr <- pr[order(pr$cluster, pr$timeBef), ]
      newList <- as.character(pr$group)
      objNew <- try(update(obj,
                   curveid = group, data = pr), silent = T)
      objNew <- try(update(obj,
                   curveid = group, data = pr), silent = T)
      if(any(class(objNew) == "try-error")){
        # Update 18/08/2022. In some cases, a parametric model
        # cannot be fit to the permuted sample
        message("\r", appendLF = F)
        stop("\r Permutation based inference was not successful. See the documentation for more detail")
        }
      obj2New <- update(objNew, curveid = NULL) # Fitting cumulato
      Db[b] <- calcDstat(objNew, obj2New)$obsD
  }
  # pval <- mean(Db > D)
  }

  pval <- (sum(Db > D) + 1)/(B + 1)
  ## Describes the results
  if(display){
    cat("\r")
    TEST<-"Permutation test based on a Cramer-von-Mises type distance (Barreiro-Ures et al., 2019)"
    null.phrase <- "NULL HYPOTHESIS: time-to-event curves are equal"
    alt.phrase <- "ALTERNATIVE: time-to-event curves are not equal"
    tabRes <- data.frame("level" = as.character(ug), n = ni, D = Dlist$Dvals)
    row.names(tabRes) <- 1:length(tabRes[,1])
    cat(TEST)
    cat("\n")
    cat(null.phrase)
    cat("\n")
    cat("\n")
    print(tabRes)
    cat("\n")
    cat(paste("Observed D value = ", round(Dlist$obsD, 4)))
    cat("\n")
    cat(paste("P value = ", round(pval, 6)))
    cat("\n")
  }
  # pout$data.name <- paste("{",L.name,",",R.name,"}"," by ",group.name,sep="")
    pout <- list()
    pout$method <- "Kramer - Von Mises distance"
    pout$scores <- NULL
    pout$val0 <- Dlist$obsD
    pout$vali <- Db
    pout$pval <- NULL
    pout$pvalb <- pval
    pout$U <- NULL
    pout$N <- NULL
  return(invisible(pout))
}
compCDFnp <- function(obj,
    scores,
    rho = NULL,
    alternative = c("two.sided", "less", "greater"),
    B, units, display) {

  # L.name <- This is missing. To be added
  cluster <- units
  R.name <- obj$dataList$names$dName
  group.name <- obj$dataList$names$cNames

  ## find NPMLE based on all the data (ignoring group membership)
  objNull <- update(obj, curveid = NULL)
  icFIT <- objNull$fit$icfitObjFull[[1]]
  # icFIT

  ## calculate scores from fitted model
  ## code taken from ictest
  A <- icFIT$A
  k <- dim(A)[[2]]
  n <- dim(A)[[1]]
  phat <- c(0, icFIT$pf)
  Shat <- 1 - cumsum(phat)

  if (scores=="logrank1"){
    Stilde <- exp( - c(0,cumsum(phat[2:(k+1)]/Shat[1:k])) )
    ckstar <- (1/phat[2:(k+1)])*( Shat[1:k]* log(Stilde[1:k]) - Shat[2:(k+1)]* log(Stilde[2:(k+1)]) )
  } else if (scores=="logrank2"){
      ckstar <- (1/phat[2:(k+1)])*c(
      Shat[1:(k-1)]* log( Shat[1:(k-1)]) - Shat[2:(k)]* log( Shat[2:(k)]),
      Shat[k]* log( Shat[k]))
  } else if (scores=="wmw"){
    ckstar <- Shat[1:k] + Shat[2:(k+1)] - 1
  }

  p <- phat[2:(k+1)]
  tempfunc<-function(Arow){
    sum(Arow * p * ckstar)/sum(Arow * p)
  }

  cc <- apply(A, 1, tempfunc)
  # print(cc)

  ## if group is numeric but only two unique levels then treat as two sample case
  group <- as.character(rep(obj$data[,5], obj$data[,3]))

  ## for 2- or k-sample: calculate efficient score statistics, U,
  ## and sample size per group, N
  # U is the group score sum (see below)
  # tapply(cc, group, sum)

  ug <- unique(group)
  ng <- length(ug)
  U <- rep(NA, ng)
  names(U) <- ug
  N <- U
  for (j in 1:ng){
    U[j] <- sum(cc[group==ug[j]])
    N[j] <- length(cc[group==ug[j]])
  }

  if (ng < 2){
    stop("Only one group found. Nothing to compare")
  # } else if(ug == 20){
  #   # 2-samples test
  #   # X<-cc[group==ug[1]]
  #   # Y<-cc[group==ug[2]]
  #   # # pout <- do.call("permTS", list(x=X, y=Y, alternative=alternative,
  #   #                                exact=exact,method=method, control=mcontrol))
  } else {
    # K-samples test
    # pout <- perm::permKS(x = cc, g = group,
    #                      exact = T) #, method = method) #, control = mcontrol))
    # pout
    # print(pout)
    # Calculate the resampling stat (i.e. the deviance for groups)
    cc <- cc - mean(cc) # centering
    calcTestStat<-function(x, g, Ng=ng, Ug=ug){
      mean.scores <- N <- rep(NA, Ng)
      for (j in 1:Ng){
        mean.scores[j] <- mean(x[g==Ug[j]])
        N[j] <- length(x[g==Ug[j]])
      }
      sum( N*(mean.scores^2) )
    }

    t0 <- calcTestStat(cc, group) # observed stat. quadrato scores medi * sample sizes
    # t0
    # anova(lm(cc ~ group))
    # sum(N * tapply(cc, group, mean)^2)


    # Build permutation distribution
    #    set.seed(1234321)
    ti <- rep(NA, B)

    for (i in 1:B) {
      if(is.null(cluster)){
        newList <- sample(group, replace = FALSE)
        ti[i] <- calcTestStat(cc, newList)
      } else {
        df <- data.frame(id = 1:length(cluster),
                       obj$data,
                       cluster = cluster)
        pr <- resample.cens(df, cluster = df$cluster, replace = c(F, F))
        pr$group <- obj$data[[5]]
        pr <- pr[order(pr$cluster, pr$timeBef), ]
        newList <- as.character(pr$group)
        ti[i] <- calcTestStat(cc, newList)
      }
    }
    nval <- sum(ti > t0)
    pval <- (nval + 1)/(B + 1)
  }
  # print(t0)
  ## TEST describes results
  if(display){
    TEST <- "Exact"
    if (scores=="logrank1" || scores=="logrank2") TEST <- paste(TEST,"Logrank")
    else if (scores=="wmw")  TEST <- paste(TEST,"Wilcoxon")
    TEST <- paste(TEST,"test (permutation form)")
    if (scores == "logrank1") TEST <- paste(TEST, " (Sun's scores)")
    if (scores == "logrank2") TEST <- paste(TEST, " (Finkelstein's scores)")
    null.phrase <- "NULL: time-to-event curves are equal"
    alt.phrase <- "ALTERNATIVE: time-to-event curves are not equal"
    tabRes <- data.frame("level" = as.character(ug), n = N, Scores = U)
    row.names(tabRes) <- 1:length(tabRes[,1])
    cat(TEST)
    cat("\n")
    cat(null.phrase)
    cat("\n")
    cat("\n")
    print(tabRes)
    cat("\n")
    cat(paste("Observed T value: ", round(t0, 4)))
    cat("\n")
    cat(paste("Permutation P-value (B = ", B,"): ", round(pval, 6), sep = ""))
    cat("\n")
  }
    # pout$data.name <- paste("{",L.name,",",R.name,"}"," by ",group.name,sep="")
    pout <- list()
    pout$method <- scores
    pout$scores <- cc
    pout$val0 <- t0
    pout$vali <- ti
    pout$pval <- NULL
    pout$pvalb <- pval
    pout$U <- U
    pout$N <- N
  return(invisible(pout))
  }

