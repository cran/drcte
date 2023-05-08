# Function to get quantiles for drcte (time-to-event) objects
# with robust SEs. It is based on the ED.drcte function, that is
# only exposed for compatibility with drc.
# Last edited: 01/06/2021
# 1. for parametric models: the EDfct function is used
# 2. for NMPLE fit quantiles for grouped data are calculated (Farooq, 2005)
# 3. for KDE fit quantiles are calculated by using a bisection method
##########################################################################
quantile.drcte <- function(x, probs, restricted = FALSE, rate = FALSE,
                           interval = FALSE, level = ifelse(!(interval == "none"), 0.95, NULL),
                           robust = FALSE, B = 999, units = NULL,
                           display = TRUE, ...){
  object <- x

  # test whether units are in the original 'data.frame' or they
  # are given as an external vector
  if(!is.null(object$origData)){
    object$origData <- data.frame(object$origData)
  }

  if(!missing(units)){
    data <- object$origData
    if(!is.null(data)){
      tmp <- try(dplyr::select(data, {{ units }}), silent = T)
      if(!is(tmp, "try-error")){
        units <-  tmp[,1]
      }
    }
  }

  if(restricted == F){
    type = "absolute"
    respLev <- probs
  } else {
    type = "relative"
    respLev <- probs * 100
  }

  # Menage interval argument (for backward compatibility)
  if(is.logical(interval)){
    if(interval == F){
      interval <- "none"
    } else {
        if(object$fit$method == "Parametric") {
          interval <- "delta"
          }  else if(object$fit$method == "NPMLE") {
            interval <- "boot"
          } else {
            interval <- "none"
          }

    }}

  # Menage robust argument
  if(object$fit$method == "Parametric" & robust == T & is.null(units)){
    vcov <- sandwich(object)
  } else if(object$fit$method == "Parametric" & !is.null(units)){
    vcov <- vcovCL(object, cluster = units)
  }


  ED(object, respLev, interval = interval,
     clevel = NULL, level = level, type = type,
     # bound = bound,
     # od = od,
     # seed = seed,
     vcov. = vcov,
     display = display, units = units,
     rate = rate, B = B, ...)
}

# Other service functions #################################
quantileNPMLE <- Vectorize(function(start, end, count, probs,
                                 type = c("absolute", "relative"),
                                 rate = F){
# Function to calculate quantiles for grouped data. We need the start
# of interval, the end of interval and the count in each interval
  if (probs < 0 | probs > 1) stop("probs should be between 0 and 1.")
  L <- rep(start, count)
  R <- rep(end, count)
  obj <- getNPMLE(survival::Surv(L, R, type = "interval2") ~ 1)
  # print(obj)
  newStart <- obj$intmap[1,]
  newEnd <- obj$intmap[2,]
  n <- sum(count)
  newCount <- obj$pf*n
  cdf <- cumsum(obj$pf)
  cdfR <- cdf[is.finite(newEnd)]
  startR <- newStart[is.finite(newEnd)]
  endR <- newEnd[is.finite(newEnd)]

  type <- match.arg(type)
  if(type == "absolute") g <- probs else g <- probs * max(cdfR)
  if( g > max(cdfR)) {
    ret <- NA
    } else {
    pos <- head( which(cdf >= g), 1)
    t1 <- startR[pos] #ifelse(pos == 1, 0, startR[pos])
    t2 <- endR[pos] #ifelse(pos == 1, startR[pos], endR[pos])
    cdf1 <- ifelse(pos == 1, 0, cdf[pos - 1])
    cdf2 <- cdf[pos] # ifelse(pos == 1, cdf[pos - 1], cdf[pos])
    df <- data.frame(ti = c(t1, t2), q = c(cdf1, cdf2))
    # print(df)
    ret <- predict(lm(ti ~ q, data = df), newdata = data.frame(q = g))
    ret <- as.numeric(ret)
  }
  names(ret) <- paste(probs*100, "%", sep="")
  if(rate == F) {
    return(ret)
    } else {
      return( ifelse(is.na(ret), 0, 1/ret) )
      }
}, "probs")

quantileKDE <- Vectorize(function(fcn, lower = 0, upper = 1000, probs,
                                  type = c("absolute", "relative"), rate = F){
  # Function to calculate quantiles from a Kernel Distribution
  # Estimator (fcn). It uses the bisection method
  if (probs < 0 | probs > 1) stop("probs should be between 0 and 1.")
  type <- match.arg(type)
  Pmax <- fcn(10^9)
  if(type == "absolute") g <- probs else g <- probs * Pmax

  eqn <- function(x, g) fcn(x) - g
  if(g >= Pmax) {
    ret <- NA
  } else {
    ret <- uniroot(eqn, g = g, lower = lower, upper = upper)$root
    ret <- as.numeric(ret)
  }
  names(ret) <- paste(probs*100, "%", sep="")
  if(rate == F) {
    return(ret)
  } else {
    return( ifelse(is.na(ret), 0, 1/ret) )
  }
}, "probs")

quantileNPMLE.boot <- function(start, end, count, probs,
                                       type = c("absolute", "relative"),
                     rate = F, B = 200, level = 0.95,
                     cluster = NULL){

  if (any(probs < 0) | any(probs > 1)) stop("probs should be between 0 and 1.")

  # first estimate
  ED0 <- quantileNPMLE(start, end, count, probs,
                          type = type, rate = rate)

  # resample the data and populate a list
  L <- rep(start, count)
  R <- rep(end, count)
  df <- data.frame(L = L, R = R)
  if(!is.null(cluster)) gr <- rep(cluster, count) else gr <- NULL

  EDi <- list()

  for (i in 1:B){
    if(!is.null(cluster)){
      message("\r Cluster Resampling:", i, appendLF = F)
      # cat("\r Cluster Resampling:", i)
      tmp <- resample.cens(df, cluster = gr, replace = list(TRUE, TRUE))
    } else {
      message("\r Resampling:", i, appendLF = F)
      # cat("\r Resampling:", i)
      tmp <- resample.cens(df, replace = list(TRUE))
    }
    fiti <- getNPMLE(survival::Surv(tmp$L, tmp$R, type = "interval2") ~ 1)
    newStart <- fiti$intmap[1,]
    newEnd <- fiti$intmap[2,]
    n <- sum(count)
    newCount <- fiti$pf * n
    EDi[[i]] <- quantileNPMLE(newStart, newEnd, newCount, probs,
                         type = type,
                         rate = rate)
  }

  # fit NPMLE di base per aver gli intervalli
  cat("\n")
  temp <- do.call(rbind, EDi)
  level <- (1 - level)/2
  bootn <- apply(temp, 2, function(x) length(x[!is.na(x)]))
  bootq <- apply(temp, 2, mean, na.rm = TRUE)
  bootq2 <- apply(temp, 2, median, na.rm = TRUE)
  bootSE <- apply(temp, 2, sd, na.rm = TRUE)
  bootLow <- apply(temp, 2, quantile, na.rm = TRUE, probs = level)
  bootUp <- apply(temp, 2, quantile, na.rm = TRUE, probs = 1 - level)
  df <- list(n = bootn, Mean = bootq, Median = bootq2, SE = bootSE, Lower = bootLow,
                   Upper = bootUp)
  df <- as.data.frame(df)
  return(df)
  }

quantileKDE.boot <- function(start, end, count, probs,
                                       type = c("absolute", "relative"),
                     rate = F, B = 1000, seed = 1234, level = 0.95){

  if (any(probs < 0) | any(probs > 1)) stop("probs should be between 0 and 1.")
  # B <- 1; seed <- 1234; type = "absolute"; rate = F; probs = c(0.1, 0.2)
  # set.seed(seed)
  newdata <- simulateTE(start, end, count, B = B)
  fcns <- lapply(newdata, function(x) Kest.boot(x[,1], x[,2], x[,3])$Fh)
  values <- lapply(fcns, quantileKDE, lower = 0, upper = 1000, probs = probs,
         type = type, rate = rate)
  temp <- do.call(rbind, values)
  level <- (1 - level)/2
  bootn <- apply(temp, 2, function(x) length(x[!is.na(x)]))
  bootq <- apply(temp, 2, mean, na.rm = TRUE)
  bootq2 <- apply(temp, 2, median, na.rm = TRUE)
  bootSE <- apply(temp, 2, sd, na.rm = TRUE)
  bootLow <- apply(temp, 2, quantile, na.rm = TRUE, probs = level)
  bootUp <- apply(temp, 2, quantile, na.rm = TRUE, probs = 1 - level)
  df <- list(n = bootn, Mean = bootq, Median = bootq2, SE = bootSE, Lower = bootLow,
             Upper = bootUp)
  df <- as.data.frame(df)
  return(df)
  }


# quantileSG <- function(time, counts, probs, nSeeds,
#                        se.fit = F, rate = F, type = 1){
#   # Quantiles for grouped data (Farooq, 2005). It requires the end of
#   # interval (first interval begins at 0) the counts in each interval
#   # and the total number of individuals. Superseded !!!!!!!!!!!!!!!!!
#   ret <- c(); g <- probs
#   for(i in 1:length(g)){
#     # For each percentile
#     dec <- (nSeeds + type - 1) * g[i]/100
#     if( dec > sum(counts)) { ret[i] <- Inf; next() }
#     pos <- head( which(cumsum(counts) >= dec), 1)
#     t1 <- ifelse(pos == 1, 0, time[pos - 1])
#     t2 <- time[pos]
#     N1 <- ifelse(pos == 1, 0, cumsum(counts)[pos - 1])
#     N2 <- cumsum(counts)[pos]
#     ret[i] <- t1 +  (dec - N1) * abs( (t2 - t1)/(N2 - N1) )
#   }
#   names(ret) <- paste(probs, "%", sep="")
#   if(rate == F) {  return(ret)} else {
#     if(se.fit == F){return(1/ret)} else {
#       if(se.fit == T){
#         estimate <- 1/ret
#         se <- NULL
#           # se <- bootSGr.old(time, counts, probs, nSeeds, type)
#           return(list( estimate = estimate, se = se)) }
#       } } }



# bootSGr.old <- function(time, counts, probs, nSeeds, type=1){
# temp <- data.frame()
# for(i in 1:1000){
#   bySeed <- c( rep(time, counts), rep(10000, nSeeds - sum(counts)) )
#   # set.seed(1234)
#   resVec <- sample(bySeed, replace = T)
#   resCount <- as.numeric( table(cut(resVec, breaks=c(0, sort(unique(time)))) ) )
#   values <- quantileSG.old(time, resCount, probs, nSeeds, type)
#   #pMax <- sum(resCount)/nSeeds
#   report <- c(1/values)
#   temp <- rbind(temp, report)
# }
# names(temp) <- c(paste(probs, "%", sep = "") )
# #apply(temp, 2, mean)
# #apply(temp, 2, median)
# bootSE <- apply(temp, 2, sd)
# bootSE
# }

bootSGpMax <- function(time, counts, nSeeds){
  temp <- c()
  for(i in 1:1000){
    bySeed <- c( rep(time, counts), rep(10000, nSeeds - sum(counts)) )
    resVec <- sample(bySeed, replace = T)
    resCount <- as.numeric( table(cut(resVec, breaks=c(0, sort(unique(time)) )) ) )
    pMax <- sum(resCount)/nSeeds
    temp <- c(temp, pMax)
  }
  bootSE <- sd(temp)
  bootSE
}

pMaxFin <- function(time, counts, nSeeds, se.fit = F){
  pMax <- sum(counts)/nSeeds
  if(se.fit == T) {
    pMaxSE <- bootSGpMax(time, counts, nSeeds)
    returnList <- list(pMaxFin = pMax, se = pMaxSE)
    }else{ returnList <- list (pMaxFin = pMax) }
  return(returnList)
  }

