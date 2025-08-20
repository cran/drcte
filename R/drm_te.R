"drmte" <- function(formula, curveid, pmodels, data = NULL, subset, fct,
start, na.action = na.omit, logDose = NULL, type = "event",
control = drmteControl(), lowerl = NULL, upperl = NULL, separate = FALSE,
pshifts = NULL, varcov = NULL){

  ## Get all arguments from call ############################
  bcVal <- NULL
  bcAdd <- 0
  robust = "mean"

  # options(na.action = deparse(substitute(na.action))) - Edited 4/5/23
  op1 <- options(na.action = deparse(substitute(na.action)))
  on.exit(options(op1), add=TRUE)

  useD <- control$"useD"
  constrained <- control$"constr"
  maxIt <- control$"maxIt"
  optMethod <- control$"method"
  relTol <- control$"relTol"
  warnVal <- control$"warnVal"
  rmNA <- control$"rmNA"  # in drmEM...
  errorMessage <- control$"errorm"  # in drmOpt
  noMessage <- control$"noMessage"  # reporting finding control measurements?
  dscaleThres <- control$"dscaleThres"
  rscaleThres <- control$"rscaleThres"
  conCheck <- control$"conCheck"
  # KDEmethod <- control$"KDEmethod"

  ## Setting warnings policy
  # options(warn = warnVal) edited 4/5/2023
  op2 <- options(warn = warnVal)
  on.exit(options(op2), add=TRUE)

  ## Handling 'start' argument
  if (missing(start)) {selfStart <- TRUE} else {selfStart <- FALSE}


  ## Handling 'fct' argument #######################
  ## fct is a function, but it returns a list, when it is passed with parentheses
  fctName <- deparse(substitute(fct))
  if(substr(fctName, 1, 4) == "KDE(") fctName <- "KDE()"
  if(substr(fctName, 1, 6) == "NPMLE(") fctName <- "NPMLE()"

  if ( (!is.list(fct)) && (!is.function(fct)) ) {stop("No function or list given in argument 'fct'")}
  if (is.function(fct))
  {
    fct <- fct2list(fct, 2)
  }

  if (is.null(names(fct))) {fct$"fct" <- fct[[1]]; fct$"ssfct" <- fct[[2]]; fct$"names" <- fct[[3]]}

  if (!is.function(fct$"fct"))
  {
      stop("First entry in list to 'fct' NOT a function")
  } else {
      drcFct <- fct$"fct"
  }

  if (is.null(fct$"ssfct")) {noSSfct <- TRUE} else {noSSfct <- FALSE}
  if ((!is.function(fct$"ssfct")) && selfStart)
  {
      stop("Neither self starter function nor starting values provided")
  } else {
      ssfct <- fct$"ssfct"
  }
  if (is.null(fct$"names") || (!is.character(fct$"names")))
  {
    if(fctName == "NPMLE()"){
      parNames <- NULL
      numNames <- 1
    } else stop("Parameter names (as vector a strings) are NOT supplied")
  } else {
      parNames <- fct$"names"
      numNames <- length(parNames)
  }

  ## Checking whether or not first derivates are supplied
  isDF <- is.function(fct$"deriv1")
  if ( (useD) && (isDF) )
  {
      dfct1 <- fct$"deriv1"  # deriv1  # [[4]]
  } else {
      dfct1 <- NULL
  }

  ## Checking whether or not second derivates are supplied
  if ( (useD) && (is.function(fct$"deriv2")) )
  {
      dfct2 <- fct$"deriv2"
  } else {
      dfct2 <- NULL
  }

  ## Storing call details
  callDetail <- match.call()

  ## Handling the 'formula', 'curveid' and 'data' arguments ##########
  anName <- deparse(substitute(curveid))  # storing name for later use MOVED UP
  if (length(anName) > 1) {anName <- anName[1]}  # to circumvent the behaviour of 'substitute' in do.call("multdrc", ...)
  if (nchar(anName) < 1) {anName <- "1"}  # in case only one curve is analysed

  mf <- match.call(expand.dots = FALSE)
  nmf <- names(mf)
  mnmf <- match(c("formula", "curveid", "data", "subset", "na.action", "weights"), nmf, 0)

  mf[[1]] <- as.name("model.frame")

  mf <- eval(mf[c(1, mnmf)], parent.frame())  #, globalenv())
  mt <- attr(mf, "terms")
  varNames <- names(mf)[c(2, 1)]  # Rigido x1 + y
  varNames0 <- names(mf) # tutte le variabili


  # only used once, but mf is overwritten later on
  dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
  xDim <- ncol(as.matrix(dose))
  resp <- model.response(mf, "numeric")

  if (is.null(resp))
  {
      # Non so a cosa serva..... da togliere, credo
      if (xDim > 1) {doseForResp <- dose[, 1]} else {doseForResp <- dose}
      resp <- ppoints(doseForResp, 0.5)[order(doseForResp)]  # just one option
      varNames[1] <- varNames[2]
      varNames[2] <- "proportion"
  }

  origDose <- dose
  origResp <- resp  # in case of transformation of the response
  lenData <- length(resp)
  numObs <- length(resp)

  ## Retrieving weights (useless...)
  wVec <- model.weights(mf)
  if (is.null(wVec))
  {
      wVec <- rep(1, numObs)
  }


  ## Finding indices for missing values
  missingIndices <- attr(mf, "na.action")
  if (is.null(missingIndices)) {removeMI <- function(x){x}} else {removeMI <- function(x){x[-missingIndices,]}}

  ## Handling "curveid" argument
  assayNo <- model.extract(mf, "curveid")
  if (is.null(assayNo))  # in case not supplied
  {
      assayNo <- rep(1, numObs)
  }
  uniqueNames <- unique(assayNo)
  colOrder <- order(uniqueNames)
  uniqueNames <- as.character(uniqueNames)

  ## Re-enumerating the levels in 'assayNo' and 'pmodels'
  assayNoOld <- assayNo

  ## Separate fitting is only possible for parametric curves
  if(fctName == "NPMLE()" || fctName == "KDE()") separate <- FALSE

  ## Detecting control measurements

  ## Defining helper function
  colConvert <- function(vec)
  {
      len <- length(vec)
      assLev <- unique(vec)

      retVec <- rep(0,len)
      j <- 1
      for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

      return(retVec)
  }

  assayNo <- colConvert(assayNoOld)
  assayNames <- as.character(unique(assayNoOld))
  numAss <- length(assayNames)

  if (xDim > 1) {tempDoseVec <- dose[, 1]} else {tempDoseVec <- dose}

  uniqueDose <- lapply(tapply(tempDoseVec, assayNoOld, unique), length)
  udNames <- names(uniqueDose[uniqueDose == 1])

  if ( (conCheck) && (length(udNames) > 0) ){

      cm <- udNames
      if (!noMessage)
      {
          cat(paste("Control measurements detected for level: ", udNames, "\n", sep = ""))

          if (separate)
          {
              stop("Having a common control when fitting separate models does not make sense!\n")
          }
      }
      conInd <- assayNoOld %in% udNames
      assayNo[conInd] <- (assayNo[!conInd])[1]
      ciOrigIndex <- unique(assayNo)
      ciOrigLength <- numAss

  ## Updating names, number of curves and the enumeration (starting from 1)

  assayNames <- as.character(unique(assayNoOld[!conInd]))
  numAss <- length(assayNames)
  assayNo <- colConvert(assayNo)
  cm <- NULL

  } else {

    cm <- NULL
    ciOrigIndex <- unique(assayNo)
    ciOrigLength <- numAss
  }

  ## Pooling data from different curves
  if ((separate) && (numAss < 2))
  {
#        warning("Nothing to pool", call. = FALSE)
        warning("Only one level: separate = TRUE has no effect", call. = FALSE)
        separate <- FALSE
    }
    if ((separate) && (!missing(pmodels)))
    {
        warning("Separate fitting switched off", call. = FALSE)
        separate <- FALSE
    }

    if (separate) {
      # Separate fitting #####
       # return(idrm(dose, resp, assayNo, wVec, fct, type))
       # return(idrm(dose, resp, assayNoOld, wVec, fct, type, control))
       # First of all, check whether separate = TRUE. In this case,
       # call the internal function drmte_sep
      if(grepl("L.3(",  fctName, fixed=TRUE) |
         grepl("LN.3(", fctName, fixed=TRUE) |
         grepl("W1.3(", fctName, fixed=TRUE) |
         grepl("W2.3(", fctName, fixed=TRUE) |
         grepl("G.3(",  fctName, fixed=TRUE) |
         grepl("lognormal",  fctName, fixed=TRUE) |
         grepl("exponential",  fctName, fixed=TRUE) |
         grepl("loglogistic",  fctName, fixed=TRUE) ){

        ## edited on 18/08/22: retransform into a factor assayNoOld
        # to avoid problems when the dataset is subsetted

         returnList <- by(data, factor(assayNoOld),
                         function(g)  drmte_sep1(formula = formula,
                                               data = g, subset = subset,
                                               fct = fct, start = start, na.action = na.action,
                                               control = control,
                                               lowerl = lowerl, upperl = upperl))
      } else {
        ## edited on 18/08/22: retransform into a factor assayNoOld
        # to avoid problems when the dataset is subsetted
        returnList <- by(data, factor(assayNoOld),
                         function(g) drmte_sep2(formula = formula,
                                                data = g, subset = subset,
                                                fct = fct, start = start, na.action = na.action,
                                                control = control,
                                                lowerl = lowerl, upperl = upperl)
                           )

      }
      sepList <- returnList
      # return(returnList)
      if( all(unlist(lapply(returnList, inherits, what = "character")), T) ){
        stop("Selected model could not be fit to any of the curveid levels!")
      }
      returnList <- sepFit2obj(returnList)
      returnList$"dataList"$dose <- as.numeric(returnList$"dataList"$dose)
      returnList$"dataList"$names$dName <- varNames0[3:(3+xDim-2)]
      returnList$"dataList"$names$orName <- varNames0[1]
      returnList$"dataList"$names$cNames <- anName
      # class(sepList) <- "list"
      returnList$separateFit <- sepList
      class(returnList) <- c("drcteList", "drcte", "drc")
      # class(returnList) <- c("drcte", "drc")
      returnList$call <- NULL
      return(returnList)
  }


    ## Handling "pmodels" argument ###########
    pmodelsList <- list()
    if (missing(pmodels))
    {
        if (length(unique(assayNo)) == 1)
        {
            for (i in 1:numNames)
            {
                pmodelsList[[i]] <- matrix(1, numObs, 1)
            }
        } else {
            modelMat <- model.matrix(~ factor(assayNo) - 1, level = unique(assayNo))  # no intercept term
            colnames(modelMat) <- assayNames
            for (i in 1:numNames)
            {
                pmodelsList[[i]] <- modelMat
            }
        }
    } else {
        ## Handling a list or data.frame argument of "pmodels"
        if (is.null(data))
        {
          # pmodels <- eval(substitute(pmodels), envir = .GlobalEnv)
          # Edited on 8/5/2023
          pmodels <- eval(substitute(pmodels))
          # print(pmodels)
        } else {
          pmodels <- eval(substitute(pmodels), envir = data, enclos = parent.frame())
          # print(pmodels)
        }

        if (is.data.frame(pmodels))
        {
            lenCol <- ncol(pmodels)
            pmodelsMat <- matrix(0, numObs, lenCol)

            for (i in 1:lenCol)
            {
                if (length(unique(pmodels[,i])) == 1)
                {
                    pmodelsList[[i]] <- matrix(1, numObs, 1)
                    pmodelsMat[,i] <- rep(1, numObs)
                } else {
                    mf <- eval(model.frame(~factor(pmodels[,i]) - 1), parent.frame())  # converting to factors
                    mt <- attr(mf, "terms")

                    mf2 <- model.matrix(mt, mf)
                    ncmf2 <- ncol(mf2)

                    mf3 <- removeMI(mf2)
                    pmodelsList[[i]] <- mf3
                    pmodelsMat[, i] <- mf3 %*% c(1:ncmf2)
                }
            }
        } else {

            if (is.list(pmodels))
            {
                lenCol <- length(pmodels)
                pmodelsMat <- matrix(0, length(resp), lenCol)

                for (i in 1:lenCol)
                {
                    if (paste(as.character(pmodels[[i]]), collapse = "") == "~1")
                    {
                        pmodelsList[[i]] <- matrix(1, numObs, 1)
                        pmodelsMat[,i] <- rep(1, numObs)
                    } else {
                        mf <- eval(model.frame(pmodels[[i]], data=data), parent.frame())
                        mt <- attr(mf, "terms")

                        mf2 <- model.matrix(mt, mf)
                        ncmf2 <- ncol(mf2)

                        mf3 <- removeMI(mf2)
                        pmodelsList[[i]] <- mf3

                        pmodelsMat[,i] <- mf3%*%c(1:ncmf2)
                    }
                }
            }
        }
    }


    ## Re-setting na.action
    # options(na.action = "na.omit") edited 4/5/2023
    op3 <- options(na.action = "na.omit")  # the default
    on.exit(options(op3), add=TRUE)

    ## Transforming dose value if they are provided as log dose
    if ( !is.null(logDose) && is.numeric(logDose) )
    {
       origDose <- dose
       dose <- logDose^dose
    }

    ## Finding parameters for the control measurements which will not be estimated
    pmodelsList2 <- list()
    for (i in 1:numNames)
    {
        colNames <- colnames(pmodelsList[[i]])

        if ( (!is.null(cm)) && (!is.null(colNames)) )
        {
            accm <- as.character(cm)
            pos <- grep(accm, colNames)
            if (length(pos) == 0)
            {
                candCol <- pmodelsList[[i]][, 1]
                if ( !(length(assayNoOld[candCol==1])==0) && (all(assayNoOld[candCol==1] == accm)) )
                {
                    pos <- 1  # the control measurements correspond to the "Intercept" term
                }
            }
        } else {pos <- numeric(0)}


        ## Defining 'pmodelsList2' from 'pmodelsList'
        if ((length(pos) > 0) && !(upperPos == i) )
        {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]][, -pos])  # column is removed
        } else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])  # column is kept
        }
    }

    for (i in 1:numNames)
    {
        if (ncol(pmodelsList[[i]]) > numAss)
        {
            pmodelsList2[[i]] <- model.matrix(~factor(assayNo) - 1)
            colnames(pmodelsList2[[i]]) <- assayNames
        } else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])  # columns are kept
        }
    }



    ## Constructing vectors 'ncclVec' and 'parmPos' used below
    ncclVec <- rep(0, numNames)
    for (i in 1:numNames)
    {
        ncclVec[i] <- ncol(pmodelsList2[[i]])  # ncol(as.matrix(pmodelsList2[[i]]))
    }
    parmPos <- c(0, cumsum(ncclVec)[-numNames])

    ## Constructing parameter names
    pnList <- drmParNames(numNames, parNames, pmodelsList2)
    parmVec <- pnList[[1]]
    parmVecA <- pnList[[2]]
    parmVecB <- pnList[[3]]

    ## Defining with indices for the individual parameters in the model
    parmIndex <- list()
    for (i in 1:numNames)
    {
        parmIndex[[i]] <- parmPos[i] + 1:ncclVec[i]
    }

    ## Scaling of dose and response values
    scaleFct <- fct$"scaleFct"
    origDoseVal <- origDose
    if (!is.null(scaleFct))  # && (is.null(lowerl)) && (is.null(upperl)) )
    # currently the scaling interferes with constraining optimization
    {
        # Defining scaling for dose and response values

      doseScaling <- 10^(floor(log10(median(dose))))
        if ( (is.na(doseScaling)) || (doseScaling < dscaleThres) )
        {
            doseScaling <- 1
        }

        respScaling <- 10^(floor(log10(median(resp))))
        if ( (is.na(respScaling)) || (respScaling < rscaleThres) || (!identical(type, "continuous")) || (!is.null(bcVal)) )
        {
            respScaling <- 1
        }

        ## Starting values need to be calculated after BC transformation!!!

        # Retrieving scaling vector
        longScaleVec <- rep(scaleFct(doseScaling, respScaling),
                            as.vector(unlist(lapply(parmIndex, length))))

    } else {
        doseScaling <- 1
        respScaling <- 1
        longScaleVec <- 1
    }

    ## Constructing vector of initial parameter values
    startVecList <- list()

    # if(!noSSfct)
    # {
    #     startMat <- matrix(0, numAss, numNames)
    #     lenASS <- length(formals(ssfct))
    #     if (lenASS > 1)
    #     # in case doseScaling and respScaling arguments are available
    #     # scaling is done inside ssfct()
    #     {
    #         doseresp <- data.frame(x = dose, y = origResp)
    #         ssFct <- function(dframe){ssfct(dframe, doseScaling, respScaling)}
    #
    #     } else {
    #     # scaling is explicitly applied to the dose and response values
    #         doseresp <- data.frame(x = dose / doseScaling, y = origResp / respScaling)
    #         ssFct <- ssfct
    #     }
    #     isfi <- is.finite(dose)  # removing infinite dose values

    ## Combining curves (NPMLE fit) ###############
        if (identical(type, "event") | identical(type, "event2"))
        {
          # If necessary, several curves are combined by using different methods
          # based on the NPcdf function.

          df <- data.frame(x = dose, idVar = assayNo, y = origResp) # dataset non scalato
          nColdf <- length(df[1,])

          # Lista di lists
          # print(df[,nColdf])
          # stop()
          NPcdfList <- list()

          NPcdfList <- plyr::dlply(df, 3:(nColdf - 1), function(x) NPcdf(x[,1], x[,2], x[,nColdf]))
          # naiveStart <- plyr::ddply(df, 3:(nColdf - 1), function(x) NPcdf(x[,1], x[,2], x[,nColdf])$Type4)
          # naiveEnd <- plyr::ddply(df, 3:(nColdf - 1), function(x) NPcdf(x[,1], x[,2], x[,nColdf])$Type3)
          # npmle <- plyr::ddply(df, 3:(nColdf - 1), function(x) NPcdf(x[,1], x[,2], x[,nColdf])$Type1plot)
          # icFit <- plyr::ddply(df, 3:(nColdf - 1), function(x) NPcdf(x[,1], x[,2], x[,nColdf])$Type1)

          naiveStart <- plyr::ldply(NPcdfList, function(x) x$Type4)
          naiveEnd <- plyr::ldply(NPcdfList, function(x) x$Type3)
          npmle <- plyr::ldply(NPcdfList, function(x) x$Type1plot)
          icFit <- plyr::ldply(NPcdfList, function(x) x$Type1)

          # New curveid, combining the environmental variables, if any
          ssSel2 <- naiveStart
          ssSel <- naiveEnd
          assayNoNew <- ssSel$idVar

          if(xDim <= 2){
              # assayNoNew <- ssSel$idVar
              doseresp <- data.frame(time = ssSel2$time, cdf = ssSel2$cdf)

          } else {
            # assayNoComb <- data.frame(ssSel[,1:(xDim - 1)])
            # assayNoComb <- do.call(paste, c(assayNoComb, sep=":"))
            doseresp <- data.frame(time = ssSel$time, ot = ssSel[,1:(xDim - 2)],
                                   cdf = ssSel$cdf)
          }
        } else {
          isFinite <- is.finite(doseresp[, 2])
        }


    ## Calculating initial estimates for the parameters #########
    ## using the self starter
    if(!noSSfct)
    {
      startMat <- matrix(0, numAss, numNames)

      lenASS <- length(formals(ssfct))

        if (lenASS > 1)

        # in case doseScaling and respScaling arguments are available
        # scaling is done inside ssfct()
        {
            # doseresp <- data.frame(x = dose, y = origResp)
            ssFct <- function(dframe){ssfct(dframe, doseScaling, respScaling)}

        } else {
        # scaling is explicitly applied to the dose and response values
            doseresp[,1] <- doseresp[,1]/doseScaling
            ssFct <- ssfct

        }
        isfi <- is.finite(dose)  # removing infinite dose values

        ## Finding starting values for each curve
        for (i in 1:numAss)
        {
            indexT1 <- (assayNoNew == i)

            if (any(indexT1))
            {

              logVec <- indexT1
              # print(doseresp[logVec, ])
              startMat[i, ] <- ssFct(doseresp[logVec, ])

            } else {
                 startMat[i, ] <- rep(NA, numNames)
            }

            ## Identifying a dose response curve only consisting of control measurements
            if (sum(!is.na(startMat[i, ])) == 1)
            {
                upperPos <- (1:numNames)[!is.na(startMat[i, ])]
#                print(upperPos)
            }
        }

        ## Transforming matrix of starting values into a vector
        nrsm <- nrow(startMat)
        for (i in 1:numNames)
        {
            sv <- rep(0, max(nrsm, ncclVec[i]))
            indVec <- 1:ncclVec[i]
            sv[1:nrsm] <- startMat[, i]
            sv <- sv[!is.na(sv)]

            isZero <- (sv == 0)
            sv[isZero] <- mean(sv)

            startVecList[[i]] <- sv[indVec]
#            print(startVecList[[i]])
        }
        startVec <- unlist(startVecList)
    } else {

        startVec <- start  # no checking if no self starter function is provided!!!
    }

    ## Checking the number of start values provided
    if (!selfStart && !noSSfct)
    {
        lenReq <- length(startVec)  # generated from self starter
        #print(lenReq); print(start)
        if (length(start) == lenReq)
        {
            startVec <- start / longScaleVec
        } else {
            stop(paste("Wrong number of initial parameter values. ", lenReq, " values should be supplied", sep = ""))
        }
    }

    ## Converting parameters
    if (selfStart)
    {
        startVec <- drmConvertParm(startVec, startMat, assayNo, pmodelsList2)
    }

    # Scaling starting values (currently not done in drmEMls)
    startVecSc <- startVec

    ## Defining function which converts parameter vector to parameter matrix
    parmMatrix <- matrix(0, numObs, numNames)
    parm2mat <- function(parm)
    {
#        parmMatrix <- matrix(0, lenData, numNames)
        for (i in 1:numNames)
        {
           parmMatrix[, i] <- pmodelsList2[[i]] %*% parm[parmIndex[[i]]]
        }
        return(parmMatrix)
    }

    ## Defining non-linear function ---------------------------------

    ## Defining model function
    # Dose: contains all independent variables (original unscaled)
    # parm2mat: function to convert parameter vector into matrix
    # drcFct: mean function
    # cm: NULL ?
    multCurves <- modelFunction(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, fct$"retFct",
                                doseScaling, respScaling, isFinite = rep(TRUE, lenData), pshifts)
    ## Defining first derivative (if available) ... used once in drmEMls()
    if (!is.null(dfct1))
    {
        dmatfct <- function(dose, parm)
        {
            dfct1(dose, parm2mat(parm))
        }
    } else {
        dmatfct <-NULL
    }

    ## Box-Cox transformation is applied
    if (!is.null(bcVal))  # (boxcox)
    {
#        varPower <- FALSE  # not both boxcox and varPower at the same time

        ## Defining Box-Cox transformation function
        bcfct <- function(x, lambda, bctol, add = bcAdd)
        {
            if (abs(lambda) > bctol)
            {
                return(((x + add)^lambda - 1)/lambda)
            } else {
                return(log(x + add))
            }
        }

        ## Setting the tolerance for Box-Cox transformation being the logarithm transformation
        ##  (same as in boxcox.default in MASS package)
        bcTol <- 0.02

#        resp <- bcfct(resp, lambda, bcTol)
        resp <- bcfct(resp, bcVal, bcTol)

        multCurves2 <- function(dose, parm)
        {
            bcfct(multCurves(dose, parm), bcVal, bcTol)
        }
    } else {multCurves2 <- multCurves}

    ## Defining estimation method -- perhaps working for continuous data
    # robustFct <- drmRobust(robust, callDetail, numObs, length(startVec))

    # if (type == "continuous")
    # {
    #     ## Ordinary least squares estimation
    #     estMethod <- drmEMls(dose, resp, multCurves2, startVecSc, robustFct, wVec, rmNA, dmf = dmatfct,
    #     doseScaling = doseScaling, respScaling = respScaling, varcov = varcov)

    # }
    # if (identical(type, "binomial"))
    # {
    #     estMethod <- drmEMbinomial(dose, resp, multCurves2, startVecSc, robustFct, wVec, rmNA,
    #     doseScaling = doseScaling)
    # }
    # if (identical(type, "Poisson"))
    # {
    #     estMethod <- drmEMPoisson(dose, resp, multCurves2, startVecSc, weightsVec = wVec,
    #                               doseScaling = doseScaling)
    # }
    # if (identical(type, "negbin1") || identical(type, "negbin2"))
    # {
    #     estMethod <- drmEMnegbin(dose, resp, multCurves2, startVecSc, weightsVec = wVec,
    #                              doseScaling = doseScaling,
    #                              dist.type = ifelse(type == "negbin1", 1, 2))
    # }

    if (identical(type, "event")){
      estMethod <- drmEMeventtime(dose, resp, multCurves2,
                                  doseScaling = doseScaling)
    } else if (identical(type, "event2")) {
      # Yet to be programmed, for uncensored observations
      stop("Uncensored observations are not supported, yet")
      # estMethod <- drmEMeventtimeUnc(dose, resp, multCurves2,
      #                             doseScaling = doseScaling)
    }

    opfct <- estMethod$opfct

    ## Defining lower and upper limits of parameters
#    if (constrained)
#    {
    if (!is.null(lowerl))
    {
        if (!is.numeric(lowerl) || !((length(lowerl) == sum(ncclVec)) || (length(lowerl) == numNames)))
        {
            stop("Not correct 'lowerl' argument")
        } else {
            if (length(lowerl) == numNames)
            {
                lowerLimits <- rep(lowerl, ncclVec)
            } else {
                lowerLimits <- lowerl
            }
        }
        constrained <- TRUE

    } else {  ## In case lower limits are not specified
        lowerLimits <- rep(-Inf, length(startVec))
    }

    if (!is.null(upperl))
    {
        if (!is.numeric(upperl) || !((length(upperl) == sum(ncclVec)) || (length(upperl) == numNames)))
        {
            stop("Not correct 'upperl' argument")
        } else {
            if (length(upperl) == numNames)
            {
                upperLimits <- rep(upperl, ncclVec)
            } else {
                upperLimits <- upperl
            }
        }
        constrained <- TRUE

    } else {  ## In case upper limits are not specified
        upperLimits <- rep(Inf, length(startVec))
    }

    lowerLimits <- lowerLimits  / longScaleVec
    upperLimits <- upperLimits  / longScaleVec

    ## Optimisation

    ## Setting derivatives
    opdfctTemp <- estMethod$"opdfct1"
    appFct <- function(x, y){tapply(x, y, sum)}

    if (!is.null(opdfctTemp))
    {
        opdfct1 <- function(parm)
        {
#            print(as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo)))
            as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo))
        }
    } else {
        opdfct1 <- NULL
    }


    ## Updating starting values
    startVecSc <- as.vector(startVecSc)  # removing names

    ## Optimising the objective function previously defined ################
    # opfct: loglik function
    # opdfct1: NULL
    # startVecSc: starting values
    # parmVec: parameter names with curveid value
    # drmOpt: function that performs the optimisation

    if(fctName == "NPMLE()"){
      ## this is not optimising, more returning the NPMLE fit
      nlsFit <- list()
      retObj <- NPcdfList
      retList1 <- lapply(retObj, function(x) x$SurvObjSum)
      retList2 <- lapply(retObj, function(x) x$SurvObj)
      # Corrected on 4/7/22: bug that scrambled the naming of curve levels
      # names(retList1) <- levels(factor(assayNoOld))
      # names(retList2) <- levels(factor(assayNoOld))
      names(retList1) <- unique(assayNoOld)
      names(retList2) <- unique(assayNoOld)

      plotFctList <- lapply(retObj, function(x) x$Fh) #x$Type1plot)
      names(plotFctList) <- levels(factor(assayNoOld))

      nlsFit$convergence <- TRUE
      nlsFit$value <- NULL
      nlsFit$par <- NULL
      nlsFit$counts <- NULL
      nlsFit$message <- NULL
      nlsFit$hessian <- NULL
      nlsFit$ovalue <- NULL
      nlsFit$method <- "NPMLE"
      nlsFit$icfitObj <- retList1
      nlsFit$icfitObjFull <- retList2



      } else if(fctName == "KDE()"){
      ## KDE fit
      nlsFit <- list()
      count <- origResp #data[,varNames0[1]]
      timeBef <- origDose[,1] #data[,varNames0[2]]
      timeAf <- origDose[,2] # data[,varNames0[3]]

      fitData <- data.frame(count, timeBef, timeAf, groups = factor(assayNo))

      if(fct$bw == "boot"){
        retObj <- by(fitData, fitData$groups, function(x) Kest.boot(x[,2], x[,3], x[,1]))
      } else {
        retObj <- by(fitData, fitData$groups, function(x) Kest(x[,2], x[,3], x[,1]))
      }
      # print(retObj$)
      pars <- unlist(lapply(retObj, function(x) x$h))
      plotFctList <- lapply(retObj, function(x) x$Fh)
      # Corrected on 4/7/22: bug that scrambled the naming of curve levels
      # names(plotFctList) <- levels(factor(assayNoOld))
      names(plotFctList) <- unique(assayNoOld)
      # print(plotFctList)
      nlsFit$convergence <- TRUE
      nlsFit$par <- pars
      nlsFit$value <- NULL
      nlsFit$counts <- NULL
      nlsFit$message <- NULL
      nlsFit$hessian <- NULL
      nlsFit$ovalue <- NULL
      nlsFit$method <- "KDE"
      nlsFit$KDEmethod <- fct$bw

      } else {
      ## parametric fit
      nlsFit <- drmOpt(opfct, opdfct1, startVecSc, optMethod, constrained, warnVal,
      upperLimits, lowerLimits, errorMessage, maxIt, relTol, parmVec = parmVec,
      traceVal = control$"trace",
      matchCall = callDetail, silentVal = control$"otrace")
    }

    if (!nlsFit$convergence) {return(nlsFit)}

    ## Manipulating after optimisation ###################
    if (identical(type, "event") | identical(type, "event2"))
    {
    # # 9/4/19 - Correction by Andrea
    # assayNo0 <- assayNo[isFinite]
    # # dose00 <- dose[, -1] #Original timeBef
    # dose0 <- dose[, -1] #Original timeAf
    # # print(doseresp)
    # #Ok se dose0 is vector
    # #Also Create an id for experimental units (Petri dishes or other)
    # if(is.vector(dose0)==T){
    #     dose <- dose0[is.finite(dose0)==T]
    #     dose1 <- dose0[is.finite(dose0)==T]
    #     dose <- as.vector(unlist(tapply(dose1, assayNo0, function(x){unique(sort(x))})))
    #     Dish <- c(); cont <- 1
    #     for(i in 1:length(dose0)) {Dish[i] <- cont; if(is.finite(dose0[i]) == F ) cont <- cont+1}
    #   } else {
    #     dose <- dose0[is.finite(dose0[,1])==T,]
    #     dose1 <- dose0[is.finite(dose0[,1])==T,]
    #     dose <- dose[order(dose[,1]),]
    #     Dish <- c(); cont <- 1
    #     for(i in 1:length(dose0[,1])) {Dish[i] <- cont; if(is.finite(dose0[i,1]) == F ) cont <- cont+1}
    #     }
    #
    #     ## Rescaling per curve id
    #     idList <- split(data.frame(dose0, resp, assayNoOld), assayNo) #For assay
    #     idList2 <- split(data.frame(dose0, resp, assayNoOld), Dish) #For Dish
    #     # print(idList)
    #     respFct <- function(idListElt)
    #     {
    #         doseVec <- idListElt[, 1]
    #         respIdx <- length(idListElt[1, ]) - 1
    #         numAssay <- idListElt[, 3]
    #         dose2 <- unique(sort(doseVec))
    #         orderDose <- order(doseVec)
    #         resp1 <- tapply(idListElt[orderDose, respIdx], doseVec[orderDose], sum)  # obtaining one count per time interval
    #         numAssayRed <- idListElt[duplicated(idListElt[orderDose, 1]) == F, 3]
    #         resp2 <- cumsum(resp1) / sum(resp1)
    #         cbind(dose2, resp2, numAssayRed)[is.finite(dose2), , drop = FALSE]
    #     }
    #
    #     #These functions here do not work properly with replicates.
    #     drList <- lapply(idList, respFct) #dose/resp per assay
    #     drList2 <- lapply(idList2, respFct) #dose/resp per dish
    #     # print(drList)
    #
    #     lapList <- lapply(drList, function(x){x[, 1]})
    #     curveidList <- as.vector(unlist(lapply(drList, function(x){x[, 3]}) ))
    #     doseList <- as.vector(unlist(lapply(drList, function(x){x[, 1]}) ))
    #     resp <- as.vector(unlist(lapply(drList, function(x){x[, 2]}))) #This are the means???
    #     resp2 <- as.vector(unlist(lapply(drList2, function(x){x[, 2]}))) #This are the data???
    #
    #     splitFactor <- factor(assayNo, exclude = NULL)
    #     listCI <- split(splitFactor, splitFactor)
    #     lenVec <- as.vector(unlist(lapply(lapList, length)))
    #     plotid <- as.factor(as.vector(unlist(mapply(function(x,y){x[1:y]}, listCI, lenVec))))
    #     if(is.vector(dose0)==T){
    #       plotid2 <- as.factor(as.vector(unlist(listCI))[is.finite(dose0)])
    #     }else{
    #       plotid2 <- as.factor(as.vector(unlist(listCI))[is.finite(dose0[,1])])
    #     }
    #     levels(plotid) <- unique(assayNoOld)
      plotid <- factor(npmle$idAssay)
    } else {
        plotid <- NULL
    }

    ## Adjusting for pre-fit scaling
    if (!is.null(scaleFct))
    {
        # Scaling the sums of squares value back
        nlsFit$value <- nlsFit$value * (respScaling^2)

        # Scaling estimates and Hessian back
#        print(longScaleVec)
        nlsFit$par <- nlsFit$par * longScaleVec
        nlsFit$hessian <- nlsFit$hessian * (1/outer(longScaleVec/respScaling, longScaleVec/respScaling))
    }

    # Testing against the ANOVA (F-test)
    nlsSS <- nlsFit$value
    nlsDF <- numObs - length(startVec)

    ## Constructing a plot function #########################
    if(fctName != "NPMLE()"){
    ## Picking parameter estimates for each curve.
    # Does only work for factors not changing within a curve!
    if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames==cm)]} else {iVec <- 1:numAss}
    pickCurve <- rep(0, length(iVec))
    for (i in iVec)
    {
       pickCurve[i] <- (1:numObs)[assayNo == i][1]
    }
    parmMat <- matrix(NA, numAss, numNames)
    fixedParm <- (estMethod$"parmfct")(nlsFit)
    parmMat[iVec, ] <- (parm2mat(fixedParm))[pickCurve, ]
    indexMat2 <- parm2mat(1:length(fixedParm))
    indexMat2 <- indexMat2[!duplicated(indexMat2), ]

    if (!is.null(cm))
    {
#        conPos <- upperPos
#        print(conPos)
        parmMat[-iVec, upperPos] <- (parm2mat(fixedParm))[assayNoOld == cm, , drop = FALSE][1, upperPos]
        # 1: simply picking the first row
    }
    rownames(parmMat) <- assayNames


    pmFct <- function(fixedParm)
    {
        if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames == cm)]} else {iVec <- 1:numAss}

        if (!is.null(cm))
        {
#            conPos <- conList$"pos"
            parmMat[-iVec, upperPos] <- (parm2mat(fixedParm))[assayNoOld == cm, , drop = FALSE][1, upperPos]
            # 1: simply picking the first row
        }
        rownames(parmMat) <- assayNames

        return(parmMat)
    }
    parmMat <- pmFct(fixedParm)  # (estMethod$"parmfct")(nlsFit) )

    ## Defining the plot function
    # print(obj)
    obj <- parmMat # Added on 1/12/2022. To check!
    pfFct <- function(obj)
    {
        plotFct <- function(dose)
        {
            if (is.vector(dose))
            {
                lenPts <- length(dose)
            } else {
                lenPts <- nrow(dose)
            }

            curvePts <- matrix(NA, lenPts, ciOrigLength)  # numAss)
            for (i in 1:numAss)
            {
                if (i %in% iVec)
                {
#                    parmChosen <- parmMat[i, ]
                    parmChosen <- parmMat[i, complete.cases(parmMat[i, ])]  # removing NAs
#                    print(parmChosen)

                    parmMat2 <- matrix(parmChosen, lenPts, numNames, byrow = TRUE)
#                    print(parmMat2)
                    curvePts[, ciOrigIndex[i]] <- drcFct(dose, parmMat2)
                } else { curvePts[, i] <- rep(NA, lenPts)}
            }
            return(curvePts)
        }

        return(plotFct)
    }
    plotFct <- pfFct(obj)
    } else {
      plotFct <- NULL
      pfFct <- NULL
      pmFct <- NULL
      indexMat2 <- NULL
    }

    ## Computation of fitted values and residuals #####
    ## 9/4/2019 - Andrea Onofri. The original routine did
    ## not appear to work with type = "event". Therefore I
    ## parted the two routines
    if (identical(type, "event") | identical(type, "event2"))
    {
        dose <- dose[,-1]
        if(fctName == "KDE()"){
          tmp_df <- data.frame(origDose, assayNo)
          predVec <- apply(tmp_df, 1, function(x) plotFctList[[ x[length(x)] ]](as.numeric(x[2])))
          # predVec <- plotFuncList[assayNo](dose)
          resVec <- rep(NA, length(origResp))
          diagMat <- matrix(c(predVec, resVec), length(origResp), 2)
          colnames(diagMat) <- c("Predicted values", "Residuals")
        } else if(fctName == "NPMLE()"){
          diagMat <- NULL
        } else {

        multCurves2 <- modelFunction(dose, parm2mat, drcFct, cm, assayNoOld, upperPos,
                                           fct$"retFct",
                                           doseScaling, respScaling,
                                           isFinite = rep(T, length(origResp)))

        predVec <- multCurves2(dose, fixedParm)

        resVec <- rep(NA, length(origResp))
        diagMat <- matrix(c(predVec, resVec), length(origResp), 2)
        colnames(diagMat) <- c("Predicted values", "Residuals")

        }
    }
    else{
    #Everything else, but time-to-event models
    # if (identical(type, "ssd"))
    # {
    #     multCurves2 <- modelFunction(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, fct$"retFct", doseScaling, respScaling, isFinite)
    # }

    predVec <- multCurves2(dose, fixedParm)
    resVec <- resp - predVec
    resVec[is.nan(predVec)] <- 0

    diagMat <- matrix(c(predVec, resVec), length(dose), 2)
    colnames(diagMat) <- c("Predicted values", "Residuals")
    }

    ## Adjusting for robust estimation: MAD based on residuals, centered at 0, is used as scale estimate
    # if (robust%in%c("median", "trimmed", "tukey", "winsor"))
    # {
    #     nlsFit$value <- (mad(resVec, 0)^2)*nlsDF
    # }
#    if (robust=="winsor")
#    {
#        K <- 1 + length(startVec)*var(psi.huber(resVec/s, deriv=1))
#    }
    # if (robust%in%c("lms", "lts"))  # p. 202 i Rousseeuw and Leroy: Robust Regression and Outlier Detection
    # {
    #     scaleEst <- 1.4826*(1+5/(numObs-length(nlsFit$par)))*sqrt(median(resVec^2))
    #     w <- (resVec/scaleEst < 2.5)
    #     nlsFit$value <- sum(w*resVec^2)/(sum(w)-length(nlsFit$par))
    # }


    ## Adding meaningful names for robust methods
    # robust <- switch(robust, median="median", trimmed="metric trimming", tukey="Tukey's biweight",
    #                          winsor="metric Winsorizing", lms="least median of squares",
    #                          lts="least trimmed squares")


    ## Collecting summary output
    sumVec <- c(NA, NA, NA, nlsSS, nlsDF, numObs)  # , alternative)
    sumList <- list(lenData = numObs,
    alternative = NULL,  # alternative,
    df.residual = numObs - length(startVec))


    ## The data set
    if (!is.null(logDose))
    {
        dose <- origDose
    }

    dataSet <- data.frame(origDose, origResp, assayNo, assayNoOld, wVec)


    # 17/11/20 CORRECTED HERE
    lengX <- ifelse(anName != "1", length(varNames0) - 1, length(varNames0))
    if(anName == "1") anName <- "curveid"
    if (identical(type, "event") | identical(type, "event2"))
    {
        #names(dataSet) <- c(varNames0[c(2, 3, 1)], anName, paste("orig.", anName, sep = ""), "weights")
        names(dataSet) <- c(varNames0[c(2:lengX, 1)], anName, paste("orig.", anName, sep = ""),
                            "weights")
      } else {
        # names(dataSet) <- c(varNames0[c(2, 1)], anName, paste("orig.", anName, sep = ""), "weights")
        names(dataSet) <- c(varNames0[c(2:lengX, 1)], anName, paste("orig.", anName, sep = ""),
                            "weights")
    }


    ## Matrix of first derivatives evaluated at the parameter estimates
    if (isDF)
    {
#        print((parmMat[assayNo, , drop = FALSE])[isFinite, , drop = FALSE])
        isFinite <- rep(T, length(origResp))
        deriv1Mat <- fct$"deriv1"(dose, (parmMat[assayNo, , drop = FALSE])[isFinite, , drop = FALSE])

    } else {
        deriv1Mat <- NULL
    }
#    deriv1Mat <- NULL

    ## Box-Cox information
    if (!is.null(bcVal))
    {
        bcVec <- list(lambda = bcVal, ci = c(NA, NA), bcAdd = bcAdd)
    } else {
        bcVec <- NULL
    }

    ## Parameter estimates
    coefVec <- nlsFit$par
    if(!is.null(coefVec)){
      names(coefVec) <- parmVec
      indexMat <- apply(t(parmMat), 2, function(x){match(x, coefVec)})
    } else {
      indexMat <- NULL
      parmMat <- NULL
    }

    ## Constructing data list ... where is it used? And it gives problems (AO)!!!!!!
    wName <- callDetail[["weights"]]
    if (is.null(wName))
    {
        wName <- "weights"
    } else {
        wName <- deparse(wName)
    }
#    dataList <- list(dose = as.vector(origDose), origResp = as.vector(origResp), weights = wVec,
    # dataList <- list(dose = origDose, origResp = as.vector(origResp), weights = wVec,
    # curveid = assayNoOld, resp = as.vector(resp),
    # names = list(dName = varNames[1], orName = varNames[2], wName = wName, cNames = anName, rName = ""))

    if (identical(type, "event") | identical(type, "event2"))
    {
      # For compatibility with 'drm()' returns the data for plotting
      # i.e., the end-point estimator

      if(xDim <= 2){
        dataListDose <- naiveEnd$time
        adVarNames <- NULL
      } else {
        dataListDose <- data.frame(naiveEnd$time, naiveEnd[,1:(xDim - 2)])
        # dataListDose <- data.frame(naiveEnd$time, naiveEnd[,xDim - 2]) # Corrected 14/02/2022
        # colnames(dataListDose) <- c("dose", varNames0[4:length(varNames0)]) # Corrected 15/03/2022
        colnames(dataListDose) <- c("dose", varNames0[4:(4 + xDim - 2 - 1)])
        # adVarNames <- c(varNames0[4:length(varNames0)])
        adVarNames <- c(varNames0[4:(4 + xDim - 2 - 1)])
      }

      dataList <- list(dose = dataListDose, origResp = naiveEnd$cdf,
                         weights = NA, curveid = naiveEnd$idVar,
                         plotid = naiveEnd$idVar, resp = naiveEnd$cdf,
                  names = list(dName = varNames0[3:(3+xDim-2)], orName = varNames0[1],
                               wName = wName, cNames = anName,
                               rNames = unique(assayNoOld)),
                  adVarNames = adVarNames
                  ) #rName = levels(factor(assayNoOld))


      if(length(varNames0) > 3) {
        for(i in 4:length(varNames0)){
          names(naiveStart)[i - 3] <- varNames0[i]
          names(naiveEnd)[i - 3] <- varNames0[i]
          names(npmle)[i - 3] <- varNames0[i]
          names(icFit)[i - 3] <- varNames0[i]
        }}
      # This is new and returns several CDF estimators
      # The first commands take the numerical 'curveid' and transform it into
      # the corresponding namings for id levels
      colVal <- xDim - 2 + 1
      naiveStart[,colVal] <- dataList$names$rNames[naiveStart[,colVal]]; names(naiveStart)[colVal] <- dataList$names$cNames
      naiveEnd[,colVal] <- dataList$names$rNames[naiveEnd[,colVal]]; names(naiveEnd)[colVal] <- dataList$names$cNames
      npmle[,colVal] <- dataList$names$rNames[npmle[,colVal]]; names(npmle)[colVal] <- dataList$names$cNames
      icFit[,colVal] <- dataList$names$rNames[icFit[,colVal]]; names(icFit)[colVal] <- dataList$names$cNames
      dataList2 <- list(naiveStart=naiveStart, naiveEnd=naiveEnd, npmle = npmle, icfit = icFit)
    }


    ## What about naming the vector of weights?
    # print(dataSet); print("OK"); stop()

    if(fctName == "KDE()" | fctName == "NPMLE()") plotFct <- plotFctList
    # if(fctName == "KDE()" | fctName == "NPMLE()") callDetail$fct <- fctName
    ## Returning the fit ######

    if(is.null(parmMat)) {
      parmMatOut <- NULL
      df.res <- NULL
      resultMat <- do.call(rbind, lapply(retList1, function(x) x[,-c(1,2)]))
      coefVec <- resultMat[ ,3]
      names(coefVec) <- row.names(resultMat)

    } else {
      parmMatOut <- t(parmMat)
      df.res <- numObs - length(startVec)
    }

    if(fctName == "KDE()") names(coefVec) <- sub("Intercept", "bandwidth", names(coefVec))
    if(fctName == "KDE()") parmVec <- sub("Intercept", "bandwidth", parmVec)
    if(fctName == "KDE()") parmVecB <- sub("Intercept", "bandwidth", parmVecB)

    isFinite <- is.finite(origDose[,2]) == T
    returnList <- list(NULL, nlsFit, list(plotFct, logDose), sumVec,
                        startVecSc * longScaleVec,
                       list(parmVec, parmVecA, parmVecB), diagMat, callDetail,
                       dataSet, parmMatOut, fct, robust, estMethod, df.res,
                       sumList, NULL, pmFct, pfFct, type, indexMat, logDose, cm, deriv1Mat[isFinite,],
                       anName, data, wVec, dataList, coefVec, bcVec,
                       indexMat2, dataList2)


    names(returnList) <- c("varParm", "fit", "curve", "summary", "start", "parNames",
                           "predres", "call", "data",
    "parmMat", "fct", "robust", "estMethod", "df.residual",
    "sumList", "scaleFct", "pmFct", "pfFct", "type", "indexMat", "logDose", "cm", "deriv1",
    "curveVarNam", "origData", "weights",
    "dataList", "coefficients", "boxcox", "indexMat2", "ICfit")
    class(returnList) <- c("drcte", "drc")
    return(returnList)
}

