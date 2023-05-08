"predict.drcte" <- function(object, newdata, se.fit = FALSE,
                            interval = FALSE, # = c("none", "confidence", "boot"),
                            level = 0.95,
                            na.action = na.pass, # od = FALSE, vcov. = vcov,
                            npmle.type = c("interpolation", "left", "right", "midpoint"),
                            robust = FALSE, units = NULL,
                            B = 200, ...){
  ## Checking arguments
  # object <- mod
  # newdata <- 180
  # cluster <- NULL
  # test whether units are in the original 'data.frame' or they
  # are given as an external vector
  if(!is.null(object$origData)){
    object$origData <- data.frame(object$origData)
  }

  if(robust) se.fit <- TRUE
    if(!missing(units)){
      data <- object$origData
      if(!is.null(data)){
        tmp <- try(dplyr::select(data, {{ units }}), silent = T)
        if(!is(tmp, "try-error")){
          units <-  tmp[,1]
          }
      }
    }

  cluster <- units
  # interval <- match.arg(interval) # modified as TRUE/FALSE
  respType <- object[["type"]]
  dataList <- object[["dataList"]]
  dataSet <- object[["data"]]
  ncol <- length(dataSet[1,])
  cName <- dataList[["names"]][["cNames"]]
  cLevs <- dataList[["names"]][["rNames"]]

  if (missing(newdata))
  {
    # Case 1: no newdata
    dataSet <- dataSet[is.finite(dataSet[,2]) == T, ]
    doseVec <- dataSet[,-c(1,(ncol - 3):ncol)]
    groupLevels <- dataSet[,(ncol - 2)] # as.character(dataList[["plotid"]])
    varCurveId <- cLevs[groupLevels]
  } else {
    # Case 2: newdata is given.
    if(is.vector(newdata)) newdata <- data.frame(newdata)
    if(is.vector(dataList[["dose"]]) == T){
      # Only one predictor in the model
      dName <- dataList[["names"]][["dName"]]
      if (any(names(newdata) %in% dName))
      {
        doseVec <- newdata[, dName]
      } else {
        doseVec <- newdata[, 1]
        #         warning("Dose variable not in 'newdata'")
      }
    } else {
      # More than one predictor in the model (eg: HT models)
      doseDim <- ncol(dataList[["dose"]]) - 1
      if(length(newdata[1,]) <= doseDim) stop("The number of dependent variables in newdata is not equal to the number of dependent variables in the model")
      doseVec <- newdata
    }

    # With newdata, if more than one time-to-event curve is defined,
    # predictions are made for all curves
    nameLevels <- levels(factor(dataSet[,(ncol - 1)]))
    if(is.vector(doseVec)){
      groupLevels <- rep(nameLevels, each = length(doseVec))
      doseVec <- rep(doseVec, length(nameLevels))
    } else {
      groupLevels <- rep(nameLevels, each = length(doseVec[,1]))
      doseVec <- doseVec[rep(seq_len(nrow(doseVec)), length(nameLevels)),]
    }
    varCurveId <- groupLevels
  }

  # counts the number of predictions to be derived
  noNewData <- length(groupLevels)

  ## Retrieving matrix of parameter estimates
  parmMat <- object[["parmMat"]]
  if(!is.null(parmMat)) pm <- t(parmMat[, groupLevels, drop = FALSE])

  ## Retrieving variance-covariance matrix for parametric fits
  if(object$fit$method != "KDE" & object$fit$method != "NPMLE"){
    sumObj <- summary(object) #, od = od) #summary(object, od = od)

    if(is.null(cluster) & !robust){
      vcovMat <- vcov(object)
    } else if(is.null(cluster) & robust){
      vcovMat <- vcovCL(object)
    } else {
      vcovMat <- vcovCL(object, cluster = factor(cluster))
    }

  }

  ## Defining index matrix for parameter estimates
  indexMat <- object[["indexMat"]]

  ## Calculating predicted values
  retMat <- matrix(0, noNewData, 4)
  colnames(retMat) <- c("Prediction", "SE", "Lower", "Upper")
  objFct <- object[["fct"]]

  if(object$fit$method != "KDE" & object$fit$method != "NPMLE"){
    retMat[, 1] <- objFct$"fct"(doseVec, pm)
  } else if (object$fit$method == "KDE"){
    # KDE
    listaFun <- function(x, y) object$curve[[1]][[x]](y)
    retMat[, 1] <- as.numeric(mapply(listaFun, groupLevels, doseVec))
  } else {
    # NPMLE
    NPMLE.method = match.arg(npmle.type)
    listaFun <- function(x, y) object$curve[[1]][[x]](y, NPMLE.method)
    retMat[, 1] <- as.numeric(mapply(listaFun, groupLevels, doseVec))
    }

  ## Calculating the quantile to be used in the confidence intervals
  # if (!identical(interval, "none")) {
  if (interval == TRUE) {
    tquan <- qnorm(1 - (1 - level)/2)
  }

  ## Calculating standard errors and/or confidence intervals
  if(object$fit$method != "KDE" & object$fit$method != "NPMLE"){

    # SEs for Parametric time-to-event models
    ## Checking if derivatives are available
    deriv1 <- objFct$"deriv1"

    if (!is.null(deriv1) & (se.fit || interval == TRUE)) {
      sumObjRV <- 0
      piMat <- indexMat[, groupLevels, drop = FALSE]
      for (rowIndex in 1:noNewData){
        parmInd <- piMat[, rowIndex]
        varCov <- vcovMat[parmInd, parmInd]
        if(is.vector(doseVec)){
          dfEval <- deriv1(doseVec[rowIndex], pm[rowIndex, , drop = FALSE])
        } else{
          dfEval <- deriv1(doseVec[rowIndex, ], pm[rowIndex, , drop = FALSE])
        }

        varVal <- as.vector(dfEval) %*% varCov %*% as.vector(dfEval)
        retMat[rowIndex, 2] <- sqrt(varVal)

        if (!se.fit)
        {
          retMat[rowIndex, 3:4] <- retMat[rowIndex, 1] + (tquan * sqrt(varVal + sumObjRV)) * c(-1, 1)
        }
        # if (identical(interval, "confidence"))
        if (interval == TRUE)
        {
          retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal)
          retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal)
        }
        }
    }

    } else if(object$fit$method == "NPMLE"){
      # NPMLE fit
      if(se.fit || interval != FALSE){
        # return(retMat[, 1])
        # } else if(se.fit == T){
        if(is.null(cluster)) df <- object$data else df <- data.frame(object$data, group = cluster)
        splitData <- by(df, object$data[,5], function(x) x)

        confCal <- function(x, times.pred, B, groups = NULL) {
          L <- rep(x[,1], x[,3])
          R <- rep(x[,2], x[,3])
          if(!is.null(groups)) gr <- rep(x[,length(x[1,])], x[,3]) else gr <- NULL
          cis <- confint.predict(L, R, times.pred, B = B, groups = gr)
          # function confint.predict is in module surv.boot.R
          cis
        }
        # splitDoseVec <- by(doseVec, dataSet[,5], function(x) x)
        splitDoseVec <- by(doseVec, groupLevels, function(x) x)

        cis <- list()
        for(i in 1:length(splitData)){
          # i <- 1
          tmp <- confCal(splitData[[i]], splitDoseVec[[i]], B = B, groups = cluster)
          cis[[i]] <- as.data.frame(tmp)
          # cis <- append(cis, as.data.frame(tmp))
        }

        tmp <- data.frame(sapply(cis, function(x) x$se))
        tmp2 <- data.frame(sapply(cis, function(x) x$lower))
        tmp3 <- data.frame(sapply(cis, function(x) x$upper))

        if(length(tmp[1,]) == 1) {
          retMat[, 2] <- as.numeric(unlist(tmp))
          indCol <- 3
          if (interval == TRUE){
            retMat[, indCol] <- as.numeric(unlist(tmp2))
            retMat[, indCol + 1] <- as.numeric(unlist(tmp3))
            colnames(retMat)[indCol] <- "Lower"
            colnames(retMat)[indCol + 1] <- "Upper"
            indCol <- indCol + 2
          }
          retMat <- data.frame(retMat)
          retMat[, indCol] <- object$dataList$names$rNames
        } else {
          colnames(tmp) <- object$dataList$names$rNames
          tmp <- stack(tmp)
          tmp2 <- stack(tmp2)
          tmp3 <- stack(tmp3)

          retMat[, 2] <- tmp[,1]
          retMat <- data.frame(retMat)
          indCol <- 3

          if (interval == TRUE){
            retMat[, indCol] <- tmp2[, 1]
            retMat[, indCol + 1] <- tmp3[, 1]
            colnames(retMat)[indCol] <- "Lower"
            colnames(retMat)[indCol + 1] <- "Upper"
            indCol <- indCol + 2
          }
          retMat[, indCol] <- tmp[,2]
        }
        retMat <- data.frame(retMat)
        retMat[, indCol + 1] <- doseVec
        colnames(retMat)[indCol] <- cName
        colnames(retMat)[indCol + 1] <- "Time"
      }

  } else if(object$fit$method == "KDE"){
    retMat[,2:4] <- NA
  }

  ## Keeping relevant indices
  keepInd <- 1
  if (se.fit) {keepInd <- c(keepInd, 2)}
  if (interval == T) {
      keepInd <- c(keepInd, 3, 4)
      retMat[, 3] <- pmax(retMat[, 3], 0)}
  retMat <- retMat[, keepInd]
  if(!is.vector(retMat)) retMat <- data.frame(retMat)

  if(!missing(newdata)) {
      if(is.vector(retMat)) {
        retMat <- cbind(newdata, "Prediction" = retMat, row.names = NULL)
      } else {
        retMat <- cbind(newdata, retMat, row.names = NULL)
      } }
  if(length(levels(factor(groupLevels))) > 1 & !is.vector(retMat)){
    retMat <- data.frame(varCurveId, retMat)
    colnames(retMat)[1] <- cName
  }

  return(retMat)
}

"predict.llogistic" <- function(object, newdata, coefs, vcov. = NULL,
                                ...){
  if (missing(newdata)) stop ("No newdata have been given")
  if (missing(coefs)) stop ("No coefficients have been given")
  class(object) <- "list"
  predict(object = object, newdata = newdata, coefs = coefs, vcov. = vcov.)
}


"predict.list" <- function(object, newdata, coefs, vcov. = NULL,
                           ...){

  if (is.null(object$fct)) stop ("No 'fct' component exist in object")
  if (missing(newdata)) stop ("No newdata have been given")
  if (missing(coefs)) stop ("No coefficients have been given")
  if(length(coefs) != length(object$names)) stop("The number of parameters is not correct")
  parmMat <- matrix(coefs, 1, length(coefs))

  if(!is.null(vcov.)) vcovMat <- vcov.

  ## Calculating predicted values

  if(is.vector(newdata)) {
    noNewData <- length(newdata)
    retMat <- matrix(0, noNewData, 2)
    colnames(retMat) <- c("Prediction", "SE")
    retMat[, 1] <- object$"fct"(newdata, parmMat)
  } else {
    noNewData <- length(newdata[,1])
    retMat <- matrix(0, noNewData, 2)
    colnames(retMat) <- c("Prediction", "SE")
    objFct <- object[["fct"]]
    if(length(newdata[1,]) == 1) newdata <- newdata[,1]
    retMat[, 1] <- object$"fct"(newdata, parmMat)
    }

  ## Calculating standard errors and/or confidence intervals
  deriv1 <- object$"deriv1"
  if (!is.null(deriv1) & !is.null(vcov.)){
    if(is.vector(newdata)){
        dfEval <- deriv1(newdata, parmMat)
    } else {
        dfEval <- deriv1(newdata, parmMat)
    }
  varVal <- dfEval %*% vcovMat %*% t(dfEval)
  retMat[, 2] <- sqrt(diag(varVal))
  }

  ## Keeping relevant indice
  if(is.null(vcov.)){
    retMat <- data.frame(newdata, prediction = retMat[,1])
  } else {
    retMat <- data.frame(newdata, retMat)
  }
  return(retMat)
}


