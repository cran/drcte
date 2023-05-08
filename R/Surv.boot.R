# This module contains functions to retreive bootstrap standard errors
# for predictions from NPMLE fits
confint.plot <- function(L, R, conf.level=.95,
                          B = 200,
                          timeEpsilon = 10^-8,
                          seed = 1234,
                          messages = FALSE,...){
    if (B < 10) stop("B must be at least 10")
    if (!is.null(seed)) set.seed(seed)
    if (messages){
        message("Bootstrap Confidence intervals can be very time consuming.")
        utils::flush.console()
    }
    fit0 <- getNPMLE(survival::Surv(L, R, type = "interval2") ~ 1)

    ## just get the bootstrap value at a little bit before and a little bit after
    ## each interval time
    ## otherwise we get the points were the likelihood is unique
    ## use timeEpsilon to define the little bit
    eps <- timeEpsilon
    times <- unique(c(0, as.vector(fit0$intmap)))
    times <- times[is.finite(times)]
    num <- length(fit0$intmap[1, ])
    times <- c(0, times-eps, times+eps)
    times <- sort(times[times >= 0])
    n <-length(L)
    nt <- length(times)
# B <- 200
    LOWER <- UPPER <- matrix(NA, B, nt)
    df <- data.frame(L = L, R = R)
    # print(B)
    for (i in 1:B){
        cat("\r", i)
        dfi <- resample.cens(df, replace = list(TRUE))
        fiti <- getNPMLE(survival::Surv(dfi$L, dfi$R, type = "interval2") ~ 1)
        LOWER[i,] <- predictCDF(fiti, times, method="right")$S
        UPPER[i,] <- predictCDF(fiti, times, method="left")$S
    }

    percci<-function(Ti, conf.level=.95){
        ### get percentile bootstrap confidence intervals
        ### see Efron and Tibshirani, p. 160 bottom
        alpha <- (1 - conf.level)/2
        B <- length(Ti)
        k <- floor((B+1)*alpha)
        if (k==0){
             warning("increase number of bootstrap samples")
             ci <- c(-Inf,Inf)
        } else {
            oTi<-Ti[order(Ti)]
            ci<-oTi[c(k,B+1-k)]
        }
        ci
    }
    # pr <- 1:100
    # percci(pr)[1]
    # quantile(pr, probs = c(0.025, 0.975), type = 9)
    calclower <- function(x, CL=conf.level){ percci(x, conf.level=CL)[1] }
    calcupper<-function(x, CL=conf.level){ percci(x,conf.level=CL)[2] }
    lower <- apply(LOWER, 2, calclower)
    upper <- apply(UPPER, 2, calcupper)
    #
    maxlower <- binom.test(n, n, conf.level=conf.level)$conf.int[1]
    lower[lower>maxlower]<-maxlower
    lower[lower<0]<-0

    minupper<-binom.test(0,n,conf.level=conf.level)$conf.int[2]
    upper[upper<minupper]<-minupper
    upper[upper>1]<-1

    # boot.se = apply(LOWER2, 2, sd)
    # boot.se2 = apply(UPPER2, 2, sd)

    list(time=times, lower = 1 - upper, upper = 1 - lower, confMethod = "modboot",
         conf.level = conf.level)
}

confint.predict <- function(L, R, pred.times = NULL,
                            conf.level=.95,
                            B = 200,
                            # seed = 1234,
                            groups = NULL, ...){
  # Obtains standard errors for predictions with a NPMLE
  # fit
  if (B < 10) stop("B must be at least 10")
  # if (!is.null(seed)) set.seed(seed)
  fit0 <- getNPMLE(survival::Surv(L, R, type = "interval2") ~ 1)
  # times <- sort(pred.times[pred.times >= 0])
  times <- pred.times[pred.times >= 0]
  n <-length(L)
  nt <- length(times)
  LOWER <- UPPER <- INTERPOL <- matrix(NA, B, nt)
  df <- data.frame(L = L, R = R)

  for (i in 1:B){

    if(!is.null(groups)){
      message("\r Cluster Resampling:", i, appendLF = F)
      dfi <- resample.cens(df, cluster = groups, replace = list(TRUE, TRUE))
    } else {
      message("\r Resampling:", i, appendLF = F)
      dfi <- resample.cens(df, replace = list(TRUE))
    }
      fiti <- getNPMLE(survival::Surv(dfi$L, dfi$R, type = "interval2") ~ 1)
      INTERPOL[i,] <- predictCDF(fiti, times)$cdf
  }


  # percci<-function(Ti, conf.level=.95){
  #       ### get percentile bootstrap confidence intervals
  #       ### see Efron and Tibshirani, p. 160 bottom
  #       alpha <- (1 - conf.level)/2
  #       B <- length(Ti)
  #       k <- floor((B+1)*alpha)
  #       if (k==0){
  #            warning("increase number of bootstrap samples")
  #            ci <- c(-Inf,Inf)
  #       } else {
  #           oTi<-Ti[order(Ti)]
  #           ci<-oTi[c(k,B+1-k)]
  #       }
  #       ci
  #   }
    # # # pr <- 1:100
    # # # percci(pr)[1]
    # # # quantile(pr, probs = c(0.025, 0.975), type = 9)
    # # calclower <- function(x, CL=conf.level){ percci(x, conf.level=CL)[1] }
    # # calcupper<-function(x,CL=conf.level){ percci(x,conf.level=CL)[2] }
    # # lower <- apply(LOWER, 2, calclower)
    # # upper <- apply(UPPER, 2, calcupper)
    # # #
    # # maxlower <- binom.test(n, n, conf.level=conf.level)$conf.int[1]
    # # lower[lower>maxlower] <- maxlower
    # # lower[lower<0]<-0
    #
    # minupper<-binom.test(0,n,conf.level=conf.level)$conf.int[2]
    # upper[upper<minupper]<-minupper
    # upper[upper>1]<-1
  message("\n")
  boot.mean <- apply(INTERPOL, 2, mean)
  boot.se <- apply(INTERPOL, 2, sd)
  boot.upper <- apply(INTERPOL, 2, quantile, probs = 0.975, type = 9)
  boot.lower <- apply(INTERPOL, 2, quantile, probs = 0.025, type = 9)
  list(time = times, mean = boot.mean, se = boot.se,
      lower = boot.lower, upper = boot.upper)
}

confint.summary <- function(L, R, pred.times = NULL,
                            conf.level=.95,
                            B = 200,
                            # seed = 1234,
                            cluster = NULL, ...){
  groups <- cluster
  if (B < 10) stop("B must be at least 10")
  # if (!is.null(seed)) set.seed(seed)
  fit0 <- getNPMLE(survival::Surv(L, R, type = "interval2") ~ 1)
  times <- sort(pred.times[pred.times >= 0])
  n <-length(L)
  nt <- length(times)
  LOWER <- UPPER <- INTERPOL <- matrix(NA, B, nt)
  df <- data.frame(L = L, R = R)

  for (i in 1:B){

    if(!is.null(groups)){
      cat("\r Cluster Resampling:", i)
      dfi <- resample.cens(df, cluster = groups, replace = list(TRUE, TRUE))
    } else {
      cat("\r Resampling:", i)
      dfi <- resample.cens(df, replace = list(TRUE))
    }
      fiti <- getNPMLE(survival::Surv(dfi$L, dfi$R, type = "interval2") ~ 1)
      INTERPOL[i,] <- predictCDF(fiti, times)$cdf
  }

  # percci<-function(Ti, conf.level=.95){
  #       ### get percentile bootstrap confidence intervals
  #       ### see Efron and Tibshirani, p. 160 bottom
  #       alpha <- (1 - conf.level)/2
  #       B <- length(Ti)
  #       k <- floor((B+1)*alpha)
  #       if (k==0){
  #            warning("increase number of bootstrap samples")
  #            ci <- c(-Inf,Inf)
  #       } else {
  #           oTi<-Ti[order(Ti)]
  #           ci<-oTi[c(k,B+1-k)]
  #       }
  #       ci
  #   }
    # # # pr <- 1:100
    # # # percci(pr)[1]
    # # # quantile(pr, probs = c(0.025, 0.975), type = 9)
    # # calclower <- function(x, CL=conf.level){ percci(x, conf.level=CL)[1] }
    # # calcupper<-function(x,CL=conf.level){ percci(x,conf.level=CL)[2] }
    # # lower <- apply(LOWER, 2, calclower)
    # # upper <- apply(UPPER, 2, calcupper)
    # # #
    # # maxlower <- binom.test(n, n, conf.level=conf.level)$conf.int[1]
    # # lower[lower>maxlower] <- maxlower
    # # lower[lower<0]<-0
    #
    # minupper<-binom.test(0,n,conf.level=conf.level)$conf.int[2]
    # upper[upper<minupper]<-minupper
    # upper[upper>1]<-1

  boot.mean <- apply(INTERPOL, 2, mean)
  boot.se <- apply(INTERPOL, 2, sd)
  boot.upper <- apply(INTERPOL, 2, quantile, probs = 0.975, type = 9)
  boot.lower <- apply(INTERPOL, 2, quantile, probs = 0.025, type = 9)
  list(time = times, mean = boot.mean, se = boot.se,
      lower = boot.lower, upper = boot.upper)
}

