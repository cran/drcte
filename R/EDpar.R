EDpar <- function(object, respLev, interval = c("none", "delta", "fls", "tfls"), clevel = NULL,
level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"),
type = c("relative", "absolute"), lref, uref, bound = TRUE, od = FALSE, vcov. = vcov, # robust = false,
display = TRUE, pool = TRUE, logBase = NULL, multcomp = FALSE, ...)
{

    # Save options to restore
    oldOpt <- options()
    # oldPar <- par(no.readonly = TRUE)
    on.exit(options(oldOpt))
    # on.exit(par(oldPar), add = T)

    interval <- match.arg(interval)
    reference <- match.arg(reference)
    type <- match.arg(type)

    ## Checking 'respLev' vector ... should be numbers between 0 and 100
    # print(respLev)
    if ( (type == "relative") && (bound) )
    {
        if (any(respLev <= 0 | respLev >= 100))
        {
            stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
        }
    }

    ## Retrieving relevant quantities
    EDlist <- object$fct$"edfct"
    if (is.null(EDlist)) {stop("ED values cannot be calculated")}
    indexMat <- object$"indexMat"
    parmMat <- object$"parmMat"
    strParm0 <- sort(colnames(object$"parmMat"))

    curveNames <- colnames(object$"parmMat")
    # options(warn = -1) # edited on 4/5/2023
    # switching off warnings caused by coercion in the if statement
    suppressWarnings({
    if (any(is.na(as.numeric(curveNames))))
    {
        curveOrder <- order(curveNames)
    } else { # if names are numbers then skip re-ordering
        curveOrder <- 1:length(curveNames)
    }
    })
    # options(warn = 0)  # normalizing behaviour of warnings

    strParm0 <- curveNames[curveOrder]
    indexMat <- indexMat[, curveOrder, drop = FALSE]
    parmMat <- parmMat[, curveOrder, drop = FALSE]

    strParm <- strParm0

    #Modified on 15/2/19#############################
    if (is.function(vcov.))
            vcMat <- vcov.(object, ...)
        else vcMat <- vcov.

    ## Defining vectors and matrices
    ncolIM <- ncol(indexMat)
    indexVec <- 1:ncolIM
#    lenEB <- ncolIM
    lenPV <- length(respLev)  # used twice below
    noRows <- ncolIM * lenPV
    dimNames <- rep("", noRows)  # lenEB*lenPV, 2)
    EDmat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)
    oriMat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)

    ## Skipping curve id if only one curve is present
#    lenIV <- lenEB  # ncol(indexMat)
    if (identical(length(unique(strParm)), 1))
    {
#        strParm[1:lenIV] <- rep("", lenIV)
        strParm[indexVec] <- rep("", ncolIM)
    } else {
        strParm <- paste(strParm, ":", sep = "")
    }

    ## Calculating estimates and estimated standard errors
    rowIndex <- 1
    lenIV <- length(indexVec)
    dEDmat <- matrix(0, lenPV * lenIV, nrow(vcMat))
    for (i in indexVec)
    {
        parmChosen <- parmMat[, i]
#        print(parmChosen)
        parmInd <- indexMat[, i]
        varCov <- vcMat[parmInd, parmInd]

        if ((is.null(clevel)) || (strParm0[i] %in% clevel))
        {
        for (j in 1:lenPV)
        {
            EDeval <- EDlist(parmChosen, respLev[j], reference = reference, type = type, ...)
            EDval <- EDeval[[1]]
            dEDval <- EDeval[[2]]
#            print(c(i,j,parmInd))
#            print(dEDval)
            dEDmat[(i-1)*lenPV + j, parmInd] <- dEDval

            oriMat[rowIndex, 1] <- EDval
            oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)

            if (!is.null(logBase))
            {
                EDval <- logBase^(EDval)
                dEDval <- EDval * log(logBase) * dEDval
            }
            EDmat[rowIndex, 1] <- EDval
            EDmat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
            # dimNames[rowIndex] <- paste(strParm[i], respLev[j], sep = "")
            if(type == "absolute") rowProb <- respLev[j]*100 else rowProb <- respLev[j]
            dimNames[rowIndex] <- paste(strParm[i], rowProb, sep = "")

            rowIndex <- rowIndex + 1
        }
        } else {
            rowsToRemove <- rowIndex:(rowIndex + lenPV - 1)
            EDmat <- EDmat[-rowsToRemove, , drop = FALSE]
            dimNames <- dimNames[-rowsToRemove]
        }
    }

    ## Defining column names
    colNames <- c("Estimate", "SE")

    ## Calculating the confidence intervals
    if (interval == "delta")
    {
        intMat <- confint.basic(EDmat, level, object$"type", df.residual(object), FALSE)
        intLabel <- "Delta method"
    }

    if (interval == "tfls")
    {
#        colNames <- c( colNames, "Lower", "Upper")
#        lsVal <- log(oriMat[, 1])
#        lsdVal <- oriMat[, 2]/oriMat[, 1]
#        ciMat <- matrix(0, lenEB*lenPV, 2)
#        tquan <- qFct(1 - (1 - level)/2)
#        ciMat[, 1] <- exp(lsVal - tquan * lsdVal)
#        ciMat[, 2] <- exp(lsVal + tquan * lsdVal)

        intMat <- exp(confint.basic(matrix(c(log(oriMat[, 1]), oriMat[, 2] / oriMat[, 1]), ncol = 2),
        level, object$"type", df.residual(object), FALSE))
        intLabel <- "To and from log scale"
    }

    if (interval == "fls")
    {
#        ciMat <- matrix(0, lenEB*lenPV, 2)
#        tquan <- qFct(1 - (1 - level)/2)
        if (is.null(logBase))
        {
            logBase <- exp(1)
            EDmat[, 1] <- exp(EDmat[, 1])  # back-transforming log ED values
        }
#        ciMat[, 1] <- logBase^(oriMat[, 1] - tquan * oriMat[, 2])
#        ciMat[, 2] <- logBase^(oriMat[, 1] + tquan * oriMat[, 2])

        intMat <- logBase^(confint.basic(oriMat, level, object$"type", df.residual(object), FALSE))
        intLabel <- "Back-transformed from log scale"

        ## Dropping estimated standard errors (not relevant after back transformation)
        EDmat <- EDmat[, -2, drop = FALSE]
        colNames <- colNames[-2]
#        colNames <- c(colNames[-2], "Lower", "Upper")  # standard errors not relevant
    }

    if (identical(interval, "none"))
    {
        intLabel <- NULL
    } else {
        EDmat <- as.matrix(cbind(EDmat, intMat))
        colNames <- c(colNames, "Lower", "Upper")
    }
    dimnames(EDmat) <- list(dimNames, colNames)

    # Edited on 19/08/2022
    #rownames(EDmat) <- paste("e", rownames(EDmat), sep = ":")
    rownames(EDmat) <- paste(rownames(EDmat), "%", sep = "")

    if(type == "relative") rownames(EDmat) <- paste(rownames(EDmat), "_R", sep = "")
    resPrint(EDmat, "Estimated effective doses", interval, intLabel, display = display)

    ## require(multcomp, quietly = TRUE)
#    invisible(list(EDdisplay = EDmat, EDmultcomp = list(EDmat[, 1], dEDmat %*% vcMat %*% t(dEDmat))))

    if(multcomp)
    {
        EDmat1 <- EDmat[, 1]
        namesVec <- names(EDmat1)  # paste("e", names(EDmat1), sep = ":")
#        names(EDmat1) <- namesVec

        EDmat1VC <- (dEDmat %*% vcMat %*% t(dEDmat))[1:nrow(EDmat), 1:nrow(EDmat), drop = FALSE]
        colnames(EDmat1VC) <- namesVec
        rownames(EDmat1VC) <- namesVec


        invisible(list(#EDdisplay = EDmat,
#                       EDmultcomp = parm(EDmat[, 1], (dEDmat %*% vcMat %*% t(dEDmat))[1:nrow(EDmat), 1:nrow(EDmat), drop = FALSE])))
                        EDmultcomp = parm(EDmat1, EDmat1VC)))
    } else {
      invisible(EDmat)
    }
}

"resPrint" <- function(resMat, headerText, interval, intervalLabel, display)
{
# Note: arguments "interval", "intervalLabel" no longer used
    if (display)
    {
        cat("\n")
        cat(paste(headerText, "\n", sep = ""))
        # if (!identical(interval, "none"))
        # {
        #     intervalText <- paste("(", intervalLabel, "-based confidence interval(s))\n", sep = "")
        #     cat(intervalText)
        # }
        cat("\n")
        printCoefmat(resMat, cs.ind = 1:ncol(resMat), tst.ind = NULL, has.Pvalue = FALSE)
    }
#    invisible(resMat)
}

## Defining basic function for providing confidence intervals
"confint.basic" <- function(estMat, level, intType, dfres, formatting = TRUE)
{
  alphah <- (1 - level)/2
  #    if (type == "u") {two <- qnorm(1 - alphah)}
  #    if (type == "t") {two <- qt(1 - alphah, df.residual(object))}
  tailPercentile <- switch(intType,
                           binomial = qnorm(1 - alphah),
                           continuous = qt(1 - alphah, dfres),
                           event = qnorm(1 - alphah),
                           Poisson = qnorm(1 - alphah),
                           negbin1 = qnorm(1 - alphah),
                           negbin2 = qnorm(1 - alphah))

  estVec <- estMat[, 1]
  halfLength <- tailPercentile * estMat[, 2]
  confMat <- matrix(c(estVec -  halfLength, estVec +  halfLength), ncol = 2)

  ## Formatting matrix
  if (formatting)
  {
    colnames(confMat) <- c(paste(format(100 * alphah), "%", sep = " "), paste(format(100*(1 - alphah)), "%", sep = " "))
    rownames(confMat) <- rownames(estMat)
  }

  return(confMat)
}
