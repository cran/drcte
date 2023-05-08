# Robust standard errors for model parameters
# fully iterated jackknife - ungrouped/grouped version
# Only grouped, so far
# v.1.1 - 29/12/2021

jack.drcte <- function(obj, units = NULL) {
  if(is.null(units)) stop("Only the grouped jackknife is implemented")
  cluster <- factor(units)
  data <- obj$origData
  estim <- coef(obj)
  esN <- summary(obj)[[3]][,2]
  #numGroups <- length(levels(cluster))
  numParms <- length(estim)
  estim2 <- data.frame()
  cont <- 0
  for(i in 1:length(levels(cluster))){
    sel <- levels(cluster)[i]
    datS <- data[cluster!=sel,] #select data
    objNew <- try(update(obj, data=datS, start=estim), silent=T) #refit model
    if(any(class(objNew)=="try-error")) next
    estim2 <- rbind(estim2, coef(objNew))
    message(paste("\r Deleting group", i, "and refitting", sep=" "))
    cont <- cont + 1
  }
  names(estim2) <- names(estim)
  numGroups <- cont
  b <- sqrt( (numGroups - numParms)/numGroups * apply(((as.matrix(estim2) - matrix(rep(estim,numGroups),numGroups,numParms, byrow=T))^2), 2, sum) )
  cbind("Estimate"=estim, "SE"=esN, "Robust SE" = b)

  class(obj) <- "drc"
  sumObj <- summary(obj)
  sumObj$coefficients <- data.frame("Estimate" = sumObj$coefficients[,1],
                                    "Jackknife Mean" = apply(estim2, 2, mean),
                                    "Jackknife SE" = b)
  sumObj$varMat <- NA
  sumObj$robust <- paste("Delete-one fully-iterated jacknife")
  sumObj$resVar <- NULL
  sumObj$resamples <- estim2
  class(sumObj) <- c("summary.drcte", "summary.drc")
  return(sumObj)
}

jack.drcte2 <- function(obj, units = NULL) {
  # numFac <- length(attr(terms(as.formula(obj$call$formula)), "term.labels"))
  # counts <- obj$data[,numFac + 1]
  # L <- rep(obj$data[,1], counts)
  # R <- rep(obj$data[,2], counts)
  # numTreat <- length(levels(factor(obj$data[, numFac + 2])))
  # if(numTreat == 1) {
  #   df <- data.frame(L, R)
  # } else if(numTreat > 1) {
  #   treats <- rep(obj$data[,numFac + 2], counts)
  #   df <- data.frame(L, R, treats)
  # }
  # print(df)
  dfr <- obj$origData
  cnt <- obj$dataList$names$orName
  cnt <- dfr[,colnames(dfr) == cnt]
  dfr <- ungroup_te(dfr, cnt)

  if(!is.null(units)){
    cluster <- factor(units)
    B <- length(levels(cluster))
    res <- matrix(NA, nrow = B, ncol = length(coef(obj)))

    estim <- coef(obj)
    esN <- summary(obj)[[3]][,2]
    numParms <- length(estim)
    estim2 <- data.frame()
    cont <- 0
    for(i in 1:length(levels(cluster))){
      estim2 <- data.frame()
      sel <- levels(cluster)[i]
      datS <- dfr[cluster!=sel,] #select data
      objNew <- try(update(obj, data = datS, start = estim), silent=T) #refit model
      if(any(class(objNew) == "try-error" )) next
      estim2 <- rbind(estim2, coef(objNew))
      message(paste("Deleting group ", i, " and refitting", sep=""))
      cont <- cont + 1
      }
    names(estim2) <- names(estim)
    numGroups <- cont
    b <- sqrt( (numGroups - numParms)/numGroups * apply(((as.matrix(estim2) - matrix(rep(estim,numGroups),numGroups,numParms, byrow=T))^2), 2, sum) )
    cbind("Estimate"=estim, "SE"=esN, "Robust SE" = b)
  } else {
    estim <- coef(obj)
    esN <- summary(obj)[[3]][,2]
    numParms <- length(estim)
    estim2 <- data.frame()
    cont <- 0
    for(i in 1:length(dfr[1,])){
      datS <- dfr[-i,] # remove one datum
      print(datS)
      objNew <- try(update(obj, data = datS, start = estim), silent=T) #refit model
      print(objNew)
      stop()
      if(any(class(objNew) == "try-error" )) next
      estim2 <- rbind(estim2, coef(objNew))
      message(paste("Deleting subject ", i, " and refitting", sep=""))
      cont <- cont + 1
    }
    print(estim2)
    names(estim2) <- names(estim)
    numGroups <- cont
    b <- sqrt( (numGroups - numParms)/numGroups * apply(((as.matrix(estim2) - matrix(rep(estim,numGroups),numGroups,numParms, byrow=T))^2), 2, sum) )
    retRes <- cbind("Estimate"=estim, "SE"=esN, "Robust SE" = b)
  }
  retRes
}

