"plotData" <- function(x, xlim, confidence.level = 0.95, gridsize = 100,
                     type = c("average", "all", "bars", "none", "obs", "confidence"),
                     npmle.type = c("interpolation", "midpoint", "right", "left", "none")
                     ){
  # This function is used to explore 'drcte' objects
  # and extract the data to be used for ggplots
  # Date of editing: 19/01/2022
  # x <- mod
  object <- x
  type <- match.arg(type)
  npmle.type <- match.arg(npmle.type)
  if(class(object)[1] == "drc"){
    stop("This method only works for 'drcte' objects")
    # method <- "drc"
  } else {
    method <- object$fit$method
  }
  dataList <- object[["dataList"]]
  dlNames <- dataList[["names"]]
  doseName <- dlNames[["dName"]]
  respName <- dlNames[["orName"]]
  curveidName <- dlNames[["cNames"]]

  ## Determining logarithmic scales
  logX <- FALSE

  ## Constructing the observed data
  if(method == "NPMLE"){
    # To display the NPMLE (pkg. interval)
    dose <- object$ICfit$npmle$time
    resp <- object$ICfit$npmle$cdf

    # Edited on 4/7/22: in order to avoid that curveid levels are scrambled
    # curveid <- object$ICfit$npmle[,1]
    curveid <- object$ICfit$npmle[,1]
    curveid <- factor(curveid, levels = unique(as.character(curveid)))
    plotid <- as.numeric(curveid)

    # Several fitting methods
    x1 <- dose
    y1 <- resp
    y2 <- as.numeric(unlist( tapply(y1, plotid, function(i) c(i[-1], i[length(i)])) ))
    y3 <- as.numeric(unlist( tapply(y1, plotid, function(i) c(0, i[-1], i[length(i)])) ))
    x3 <- as.numeric(unlist( tapply(x1, plotid, function(i) c(0, i + 0.001)) ))

    xmid <- as.numeric(unlist( tapply(x1, plotid, function(i) c(0, (i[-1] + i[-length(i)])/2, i[length(i)])) ))
    ymid <- as.numeric(unlist( tapply(y2, plotid, function(i) c(0, i)) ))
    curveidmid <- tapply(y2, plotid, function(i) c(0, i))
    tms <- as.numeric(unlist(lapply(curveidmid, length)) )
    curveidmid <- as.numeric(rep(names(curveidmid), tms))
    curveidmid <- unique(curveid)[curveidmid]

    plotidmid <- tapply(y2, plotid, function(i) c(0, i))
    tms <- as.numeric(unlist(lapply(plotidmid, length)) )
    plotidmid <- rep(names(plotidmid), tms)
    # print(plotidmid)

    if(npmle.type == "interpolation"){
        plotPoints <- data.frame(x1, y1, curveid, plotid)
    } else if(npmle.type == "midpoint"){
        plotPoints <- data.frame(xmid, ymid, curveidmid, plotidmid)
    } else if(npmle.type == "right"){
        plotPoints <- data.frame(x1, y1, curveid, plotid)
    } else if(npmle.type == "left"){
        plotPoints <- data.frame(x3, y3, curveidmid, plotidmid) #, curveid, plotid)
    }

  } else {
    # To display naive end-point estimator (upd. 21/12/21)
    dose <- dataList[["dose"]]
    resp <- dataList[["origResp"]]
    if(method == "drc"){
      curveid <- dataList[["curveid"]]
    } else {
      curveid <- dlNames$rNames[dataList$curveid]
    }

    plotid <- dataList[["plotid"]]
    plotPoints <- data.frame(dose, resp, curveid, plotid)
  }

  # Capire se voglio plottare tutti i dati...

  # print(plotPoints)
  # print(c(doseName, "CDF", curveidName))
  if(is.null(curveidName)) curveidName <- "curveid"
  names(plotPoints)[1:(length(doseName) + 2)] <- c(doseName, "CDF", curveidName)

  if (!is.null(plotid))
  {  # used for event-time data
      assayNoOld <- as.vector(plotid)
  } else {
      assayNoOld <- as.vector(curveid)
  }
  uniAss <- unique(assayNoOld)
  numAss <- length(uniAss)

  # Get the prediction function
  plotFct <- (object$"curve")[[1]]
  if(method != "NPMLE" && method != "KDE"){
    logDose <- (object$"curve")[[2]]
    naPlot <- ifelse(is.null(object$"curve"$"naPlot"), FALSE, TRUE)
  }

  # Handling multiple covariates
  doseDim <- 1
  if(!is.vector(dose)){
    doseDim <- length(dose[1,])
    doseOld <- dose
    addVars <- dose[,-1]
    dose <- dose[,1]
  }

  ## Determining range of dose values
  if (missing(xlim))
  {
    if(method != "NPMLE" && method != "KDE")  xLimits <- c(min(dose), max(dose)) else xLimits <- c(0, max(dose))
  } else {
    xLimits <- xlim  # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
  }



  ## Handling small dose values
  conLevel <- round(min(dose[is.finite(dose)])) - 1
  # print(xLimits[1])
  # print(xLimits[2])


  if ((xLimits[1] < conLevel) && (logX || (!is.null(logDose))))
  {
      xLimits[1] <- conLevel
      smallDoses <- (dose < conLevel)
      dose[smallDoses] <- conLevel
      if (is.null(conName))
      {
          if (is.null(logDose)) {conName <- expression(0)} else {conName <- expression(-infinity)}
      }
  } else {
      conName <- NULL
  }
  if (xLimits[1] >= xLimits[2]) {stop("Argument 'conLevel' is set too high")}

  ## Constructing dose values for plotting
   if (doseDim == 1)
     {
     dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)
     } else {
     dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)

     ## addVars <- as.numeric(apply(as.data.frame(addVars), 2, function(col) unique(as.character(col))))
     # Corrected on 16/3/2022
     addVars <- apply(as.data.frame(addVars), 2, function(col) as.numeric(unique(as.character(col))), simplify = F)
     # print(names(addVars))
     numAdd <- doseDim - 1
     allVars <- list(time = dosePts)
     for(i in 1:numAdd) allVars[[i + 1]] <- addVars[[i]]
     names(allVars)[2:(2 + numAdd - 1)] <- names(addVars)
     dosePts <- expand.grid(allVars)
     }
  if(method == "Parametric"){
    plotMat <- plotFct(dosePts)
    } else if(method == "KDE"){
    plotMat <- lapply(plotFct, function(x) x(dosePts))
    } else {
    plotMat <- NULL
  }

  plotMat <- as.data.frame(plotMat)

   if(method == "NPMLE"){
     retData <- NULL
   } else {
     doseName <- dlNames[["dName"]]
     if(length(dlNames[["rNames"]]) > 1){
       respName <- as.character(dlNames[["rNames"]])
     } else {
       respName <- "CDF"
     }
     retData <- data.frame(dosePts, as.data.frame(plotMat))
     colnames(retData) <- c(doseName, respName)
   }

   # Melt plot data
  if(method != "NPMLE" && numAss > 1){
    # print(dlNames$cNames)
    retData <- tidyr::pivot_longer(retData, names_to = dlNames$cNames,
                          values_to = "CDF",
                          cols = c((doseDim + 1):length(retData[1,])))
    retData <- as.data.frame(retData)
    # print(retData)
  }

  returnList <- list(plotPoints = plotPoints, plotFits = retData)
  return(returnList)
  # retData <- data.frame(dosePts, as.data.frame(plotMat))
  # colnames(retData) <- c(doseName, dlNames[["cNames"]])
  # returnList <- list(plotPoints = plotPoints, plotFits = retData)
  # return(invisible(returnList))

}




