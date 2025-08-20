# plot.drcte <- function(x, ..., add = FALSE, level = NULL,
#    shading = TRUE, type = "all",
#    npmle.type = c("interpolation", "midpoint", "right", "left", "none"),
#    npmle.points = FALSE, kde.points = TRUE,
#    broken = FALSE, bp, bcontrol = NULL, conName = NULL, axes = TRUE, gridsize = 100,
#    log = "", xtsty, xttrim = TRUE, xt = NULL, xtlab = NULL, xlab, xlim,
#    yt = NULL, ytlab = NULL, ylab, ylim,
#    cex, cex.axis = 1, col = FALSE, lty, pch,
#    legend, legendText, legendPos, cex.legend = 1,
#    normal = FALSE, normRef = 1, confidence.level = 0.95)
plot.drcte <- function(x, ...,
     npmle.type = c("interpolation", "midpoint", "right", "left", "none"),
     npmle.points = FALSE, kde.points = TRUE, shading = TRUE, bp) {
  obj <- x

  if(obj$fit$method == "NPMLE" | obj$fit$method == "KDE") {
    # Non-parametric fits ###########

    # Save options to restore
    oldOpt <- options()
    oldPar <- par(no.readonly = TRUE)
    on.exit(options(oldOpt))
    on.exit(par(oldPar), add = T)

    object <- obj
    type = "all"
    npmle.type <- match.arg(npmle.type)

    ## Determining logarithmic scales
    listArg <- as.list(match.call(expand.dots = FALSE)[-1])
    if(is.null(listArg$...$type)) type <- "all"
    if(is.null(listArg$...$conName)) conName <- NULL

    if(is.null(listArg$...$log)) log <- "" else log <- listArg$...$log
    if ((log == "") || (log == "y"))
    {
        logX <- FALSE
    } else {
        logX <- TRUE
    }

    ## Determining the tick mark style for the dose axis
    if(is.null(listArg$...$xtsty))
    {
        if (logX)
        {
            xtsty <- "base10"
        } else {
            xtsty <- "standard"
        }
    }

    # Selection of data to be plotted #########
    dataList <- object[["dataList"]]

    if(obj$fit$method == "NPMLE"){
      # To display the NPMLE (pkg. interval)
      dose <- obj$ICfit$npmle$time
      resp <- obj$ICfit$npmle$cdf
      curveid <- obj$ICfit$npmle[,1]
      plotid <- obj$ICfit$npmle[,1]
    } else {
      # To display naive end-point estimator (upd. 21/12/21)
      dose <- dataList[["dose"]]
      resp <- dataList[["origResp"]]
      curveid <- dataList$names$rNames[dataList$curveid]
      plotid <- dataList$names$rNames[dataList$plotid]
    }

    assayNoOld <- as.vector(curveid)
    uniAss <- unique(assayNoOld)
    numAss <- length(uniAss)
    if(is.null(listArg$...$level)) level <- NULL
    doPlot <- is.null(level) || any(uniAss %in% level)
    # if (!doPlot) {stop("Nothing to plot")}
    if (!doPlot & !is.numeric(level)) {stop("Nothing to plot")} # Change: 31/10/21

    plotFct <- (object$"curve")[[1]]
    logDose <- NULL #(object$"curve")[[2]]
    # naPlot <- ifelse(is.null(object$"curve"$"naPlot"), FALSE, TRUE)

    ## Assigning axis names
    dlNames <- dataList[["names"]]
    doseName <- "" #dlNames[["dName"]]
    respName <- "" #dlNames[["orName"]]
    # axis names are the names of the dose variable and response variable in the original data set
    if(is.null(listArg$...$xlab))

    # if (missing(xlab)) {if (doseName == "") {xlab <- "Time"} else {xlab <- doseName}}
    # if (missing(ylab)) {if (respName == "") {ylab <- "CDF"} else {ylab <- respName}}
    if(is.null(listArg$...$xlab))
    if(is.null(listArg$...$xlab)) {if (doseName == "") {xlab <- "Time"} else {xlab <- doseName}}
    if(is.null(listArg$...$ylab)) {if (respName == "") {ylab <- "CDF"} else {ylab <- respName}}


    ## Determining range of dose values
    # print(listArg)
    if(is.null(listArg$...$xlim))
    # if (missing(xlim))
    {
        xLimits <- c(0, max(dose)) #c(min(dose), max(dose))
    } else {
        xLimits <- eval(listArg$...$xlim)
    }

    ## Handling small dose values
    if(is.null(listArg$...$bp))
    # if (missing(bp))
    {

        ## Constructing appropriate break on dose axis
        if (!is.null(logDose))  # natural logarithm
        {
            conLevel <- round(min(dose[is.finite(dose)])) - 1
        } else {
            log10cl <- round(log10(min(dose[dose > 0]))) - 1
            conLevel <- 10^(log10cl)
        }
    } else {
        conLevel <- bp
    }
    if ((xLimits[1] < conLevel) && (logX || (!is.null(logDose))))
    {
        xLimits[1] <- conLevel
        smallDoses <- (dose < conLevel)
        dose[smallDoses] <- conLevel
        if (is.null(conName))
        {
            if (is.null(logDose)) {conName <- expression(0)} else {conName <- expression(-infinity)}
        }
#        conNameYes <- TRUE
    } else {
#        conNameYes <- FALSE
        conName <- NULL
    }
    if (xLimits[1] >= xLimits[2]) {stop("Argument 'conLevel' is set too high")}

    ## Constructing dose values for plotting ########
#    if (doseDim == 1)
    if(is.null(listArg$...$gridsize)) gridsize <- 100
     dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)
#    } else {}  # No handling of multi-dimensional dose values


    ## Finding minimum and maximum on response scale
    if(obj$fit$method == "KDE"){
      if (is.null(logDose)){
        plotMat <- lapply(plotFct, function(x) x(dosePts))
      } else {
        plotMat <- lapply(plotFct, function(x) x(logDose^(dosePts)))
      }
    } else {
      plotMat <- NULL
    }

    maxR <- max(resp)
    # Edited on 4/5/2023
    # options(warn = -1)  # suppressing warning in case maximum of NULL is taken
    suppressWarnings({
    maxPM <- unlist(lapply(plotMat, max, na.rm = TRUE))
    if (max(maxPM) > maxR) {maxPM <- maxPM[which.max(maxPM)]} else {maxPM <- maxR}
    })
    # options(warn=0)

    if(is.null(listArg$...$ylim)){
      if(is.null(listArg$...$xlim)){
            yLimits <- c(0, 1) # c(min(resp), maxPM)
        } else {
            logVec <- ((dose >= xLimits[1]) & (dose <= xLimits[2]))
            yLimits <- range(resp[logVec])
        }
    } else {
        yLimits <- eval(listArg$...$ylim)
    }

    ## Setting a few graphical parameters
    par(las = 1)

    ## Cutting away original x values outside the limits
    eps1 <- 1e-8
    logVec <- !( (dose < xLimits[1] - eps1) | (dose > xLimits[2] + eps1) )
    dose <- dose[logVec]
    resp <- resp[logVec]
    assayNoOld <- assayNoOld[logVec]

    # Helper functions
    barFct <- function(plotPoints){invisible(NULL)}
    ciFct <- function(level, ...){invisible(NULL)}
    pointFct <- function(plotPoints, cexVal, colVal, pchVal, ...){
            points(plotPoints, cex = cexVal, col = colVal, pch = pchVal, ...)
        }
    # }

    ## Determining levels to be plotted
    if (is.null(level))
    {
        level <- uniAss
    } else {
        # level <- intersect(level, uniAss) Corrected 31/10/21
        level2 <- intersect(level, uniAss)
        if(length(level2) == 0) level2 <- uniAss[level]
        level <- level2
        if(any(is.na(level))) stop("One or more levels are not found") # Added: 31/10/21
    }
    lenlev <- length(level)

    ## Determining presence of legend
    if(is.null(listArg$...$legend))
    # if (missing(legend))
    {
      # Setting the legend parameter to true when there are more
      # than 1 level, by default
      if (lenlev == 1) {legend <- FALSE} else {legend <- TRUE}
    } else {
      if (listArg$...$legend == FALSE) {legend <- FALSE} else {legend <- TRUE}
    }

    ## Setting graphical parameters
    colourVec <- rep(1, lenlev)
    if(is.null(listArg$...$col)) col <- FALSE
    if (is.logical(col) && col)
    {
        colourVec <- 1:lenlev
    }
    if (!is.logical(col) && (length(col) == lenlev) )
    {
        colourVec <- col
    }
    if (!is.logical(col) && (!(length(col) == lenlev)) )
    {
        colourVec <- rep(col, lenlev)
    }

    # Helper function
    parFct <- function(gpar, lenlev, defVal = NULL) {
      # if (!missing(gpar)){
      if(!is.null(gpar)){
        if (length(gpar) == 1)
        {
            return(rep(gpar, lenlev))
        } else {
            return(gpar)
        }
      } else {
        if (is.null(defVal)) {return(1:lenlev)} else {rep(defVal, lenlev)}
            }
    }
    if(is.null(listArg$...$cex)) cex <- NULL
    if(is.null(listArg$...$lty)) lty <- NULL
    if(is.null(listArg$...$pch)) pch <- NULL
    if(is.null(listArg$...$axes)) axes <- TRUE
    cexVec <- parFct(cex, lenlev, 1)
    ltyVec <- parFct(lty, lenlev)
    pchVec <- parFct(pch, lenlev)

    ## Plotting non parametric data ######################
    levelInd <- 1:lenlev
    ## Plotting data for the first curve id
    # print(class(yLimits))
    plot(0, type = "n", xlab = xlab, ylab = ylab, log = log,
         xlim = xLimits, ylim = yLimits,
         axes = axes,
         frame.plot = TRUE)


    if(shading & obj$fit$method == "NPMLE"){
      # Plotting grey areas, if requested
      for (i in levelInd)
      {
        indVec <- level[i] == assayNoOld
        plotPoints <- cbind(dose[indVec], resp[indVec])
        # print(plotPoints)
        plotPoints <- rbind(c(0, 0), plotPoints)
        ti <- plotPoints[, 1]
        pi <- plotPoints[, 2]
        coord <- lapply(1:(length(ti) - 1), function(i){
                res <- list(xx = c(ti[i],ti[1 + i],ti[1 + i],ti[i]),
                yy = c(pi[i],pi[i],pi[1 + i],pi[1 + i]))
                return(res)
                })

        # greyCol <- 1 - i/levelInd[i]
        greyCol <- 0.9
        p <- lapply(coord, function(x) {graphics::polygon(x$xx, x$yy, col = gray(greyCol), border = NA)})
        # points(plotPoints, type = plotType, col = colourVec[i], pch = pchVec[i],
        #            cex = cexVec[i], lty = ltyVec[i], ...)

        ## Adding error bars (inutile?)
        barFct(plotPoints)

        ## Add confidence regions
        # ciFct(level=i, col = alpha(colourVec[i],0.25))

        ## Adding axes
        # drc:::addAxes(axes, cex.axis, conName, xt, xtlab, xtsty, xttrim, logX, yt, ytlab, conLevel) #, logDose)

        ## Adding axis break
        # ivMid <- drc:::brokenAxis(bcontrol, broken, conLevel, dosePts, gridsize, log, logX) #, logDose)

        ## Plotting in the case "add = TRUE" and for all remaining curve ids
        }
    }

    ## Setting the plot type
    if(obj$fit$method == "NPMLE"){
      plotType <- "b"
    } else {
      plotType <- "p"
    }

    for (i in levelInd) {
      indVec <- level[i] == assayNoOld
      x1 <- dose[indVec]
      y1 <- resp[indVec]
      y2 <- c(y1[-1], y1[length(y1)])
      xmid <- c(0, (x1[-1] + x1[-length(x1)])/2, x1[length(x1)])
      ymid <- c(0, y2)

      if(obj$fit$method == "NPMLE"){
        if(npmle.type == "interpolation"){
          plotPoints <- cbind(x1, y1)
          if(npmle.points == T) plotType <- "b" else plotType <- "l"
        } else if(npmle.type == "midpoint"){
          plotPoints <- cbind(xmid, ymid)
          plotType <- "s"
        } else if(npmle.type == "right"){
          plotPoints <- cbind(x1, y1)
          plotType <- "s"
        } else if(npmle.type == "left"){
          plotPoints <- cbind(x1, y2)
          plotType <- "s"
        }
        # print(plotPoints)
      } else {
        plotPoints <- cbind(x1, y1)
        # plotType <- "p"
        if(kde.points == T) plotType <- "p" else plotType <- "n"
      }
      # print(colourVec)
      pointFct(plotPoints, cexVec[i], colourVec[i], pchVec[i], type = plotType,
                         lty = ltyVec[i], ...)
      }

    ## Plotting fitted curves #### I need this for KDEs
    plotMat <- as.data.frame(plotMat)
    noPlot <- rep(FALSE, lenlev)
    if(obj$fit$method == "KDE")
    {
        for (i in levelInd)
        {
            indVal <- uniAss %in% level[i]

            # Da verificare se necessario
            # if ( (!naPlot) && (any(is.na(plotMat[, indVal]))) )
            # {
            #     noPlot[i] <- TRUE
            #     next
            # }
            # lines(dosePts[ivMid], plotMat[ivMid, indVal], lty = ltyVec[i],
            #       col = colourVec[i], ...)
            lines(dosePts, plotMat[, indVal], lty = ltyVec[i],
                  col = colourVec[i], ...)
            # print(plotMat); stop()
        }
    }


    ## Adding legend
    # level <- dlNames[["rNames"]] Removed on 31/10/21: error
    makeLegend <- function(colourVec, legend, legendCex, legendPos, legendText,
                         lenlev, level, ltyVec, noPlot, pchVec, type,
                         xLimits, yLimits) {
      if (!legend) {return(invisible(NULL))}

        legendLevels <- as.character(level)
        if (!is.null(legendText))
        {
            lenLT <- length(legendText)

            if (lenLT == lenlev) {legendLevels <- legendText}

            if (lenLT == 1) {legendLevels <- rep(legendText, lenlev)}
        }
        levInd <- 1:lenlev

        ## Removing line types when lines are not drawn
        # Penso sia inutile, qui
        # ltyVec[noPlot] <- 0
        # if (identical(type, "all"))
        # {
        #     ltyVec[levInd] <- 0
        # }
        # print("OK")
        # ## Removing plot symbol when no points are drawn
        # Corrected: 21/12/21
        if ( obj$fit$method == "NPMLE" & (npmle.points == F | npmle.type != "interpolation"))
        {
            pchVec[levInd] <- NA
        } else if (obj$fit$method == "KDE" & kde.points == F ) pchVec[levInd] <- NA

        ## Defining position of legend

        if(!is.null(listArg$...$legendPos)) legendPos <- eval(listArg$...$legendPos) else legendPos <- NULL

        if (!is.null(legendPos))
        {
            if ( (is.numeric(legendPos)) && (length(legendPos) == 2) )
            xlPos <- legendPos[1]
            ylPos <- legendPos[2]
        } else {
            xlPos <- xLimits[2]
            ylPos <- yLimits[2]
        }

        ## Adding legend
        legend(xlPos, ylPos, legendLevels, lty = ltyVec[levInd], pch = pchVec[levInd],
        col = colourVec[levInd], bty = "n", xjust = 1, yjust = 1, cex = legendCex)
    }
    if(is.null(listArg$...$cex.legend)) cex.legend <- 1
    if(is.null(listArg$...$legendText)) legendText <- NULL
    if(is.null(listArg$...$legendPos)) legendPos <- NULL
    # print(legend)
    # stop()
    makeLegend(colourVec, legend, cex.legend, legendPos, legendText,
               lenlev, level, ltyVec, noPlot, pchVec = pchVec,
               type, xLimits, yLimits)

    ## Resetting graphical parameter
    par(las = 0)


    if(obj$fit$method == "NPMLE"){
      retData <- NULL
    } else {
      doseName <- dlNames[["dName"]]
      if(length(dlNames[["rNames"]]) > 1){
        respName <- as.character(dlNames[["rNames"]])
      } else {
        respName <- "Prop"
      }
      retData <- data.frame(dosePts, as.data.frame(plotMat))
      colnames(retData) <- c(doseName, respName)
    }

    }
   else if(obj$fit$method == "KDEfun2") {
      # Da fare?
    } else {

  # case parametric curve ############
  # Handling multiple covariates
  dose <- obj$dataList[["dose"]]
  if(!is.vector(dose)) doseDim <- length(dose[1,]) else doseDim <- 1
  if(doseDim > 1) stop("Plotting is not possible with additional covariates (apart from time)")

  obj$dataList$plotid <- obj$dataList$names$rNames[obj$dataList$plotid]
  class(obj) <- "drc"
  # Corrected on 24/2/25 and 19/8/2025
  # retData <- plot(obj, log = log, add = add, level = level, type = type, broken = broken,
  #             bp = bp, bcontrol = bcontrol, conName = conName, axes = axes, gridsize = gridsize,
  #             xtsty = xtsty, xttrim = xttrim, xt = xt, xtlab = xtlab, xlab = xlab, xlim = xlim,
  #             yt = yt, ytlab = ytlab, ylab = ylab, ylim = ylim, cex = cex,
  #             cex.axis = cex.axis, col = col, lty = lty, pch = pch, legend = legend,
  #             legendText = legendText, legendPos = legendPos, cex.legend = cex.legend,
  #             normal = normal, normRef = normRef, confidence.level = confidence.level)
  listArg <- as.list(match.call(expand.dots = FALSE)[-1])
  # print(listArg)
  if(is.null(listArg$...$log)) {
    retData <- plot(obj, log = "", ...)
  } else {
    retData <- plot(obj, ...)
  }
  dlNames <- obj$dataList$names
  doseName <- dlNames[["dName"]]
  if(length(dlNames[["rNames"]]) > 1){
     respName <- as.character(dlNames[["rNames"]])
   } else {
     respName <- "Prop"
   }
   colnames(retData) <- c(doseName, respName)

  }
  return(invisible(retData))
}


