"drmEMeventtime" <- function(dose, resp, multCurves,
                                 doseScaling = 1){
    ## Defining the objective function
    opfct <- function(c)  # dose, resp and weights are fixed. resp are the counts
    {
      Fstart <- multCurves(dose[, -2] / doseScaling, c)
      dose2 <- dose[, -1]
      Fend <- multCurves(dose2 / doseScaling, c)

      ifelse(is.matrix(dose2)==T,
             Fend[is.finite(dose2[, 1])==F] <- 1,
             Fend[!is.finite(dose2)] <- 1)

      # Edited on 24/4/24: to support uncensored observations
      # needs to be tested!
      Func <- multCurves(dose2 / doseScaling - 1e-8, c)
      temp <- ifelse(dose[, 1] == dose[, 2], (Fend - Func)/1e-08, Fend - Fstart)
      temp[temp==0] <- 10e-8

      retVal <- -sum(resp * log(temp))
      return( retVal )
      # minus in front of sum() as maximization is done as minimization
    }

    ## Defining self starter function
    ssfct <- NULL

    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]
#        c( sum(log(choose(total, success))) - object$"fit"$"ofvalue", object$"sumList"$"df.residual" )

        c(
        -object$"fit"$value,  # oops a constant is missing!
        object$"sumList"$"df.residual"
        )
    }


    ## Defining functions returning the residual variance, the variance-covariance and the fixed effects estimates
    rvfct <- NULL

    vcovfct <- function(object)
    {
        solve(object$fit$hessian)
    }

    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct,
    parmfct = parmfct))
}

"drmEMeventtime_old" <-
function(dose, resp, multCurves, doseScaling = 1)
{
    ## Defining the objective function
    opfct <- function(c)  # dose, resp and weights are fixed
    {

      Fstart <- multCurves(dose[, -2] / doseScaling, c)
      dose2 <- dose[, -1]
      Fend <- multCurves(dose2 / doseScaling, c)

      ifelse(is.matrix(dose2)==T,
             Fend[is.finite(dose2[, 1])==F] <- 1,
             Fend[!is.finite(dose2)] <- 1)

      # Edited on 24/4/24: to support uncensored observations
      # needs to be tested!
      temp <- Fend - Fstart
      temp[temp==0] <- 10e-8
      # print(temp)
      # FendEx <- multCurves(dose2 / doseScaling - 10e-8, c)
      # temp[temp==0] <- (Fend - FendEx)/10e-8

      return( -sum(resp * log(temp)) )
      # minus in front of sum() as maximization is done as minimization
    }


    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]
#        c( sum(log(choose(total, success))) - object$"fit"$"ofvalue", object$"sumList"$"df.residual" )

        c(
        -object$"fit"$value,  # oops a constant is missing!
        object$"sumList"$"df.residual"
        )
    }


    ## Defining functions returning the residual variance, the variance-covariance and the fixed effects estimates
    rvfct <- NULL

    vcovfct <- function(object)
    {
        solve(object$fit$hessian)
    }

    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct,
    parmfct = parmfct))
}


"drmLOFeventtime" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}

"drmOpt" <- function(opfct, opdfct1, startVec, optMethod, constrained,
                     warnVal, upperLimits, lowerLimits, errorMessage,
                     maxIt, relTol, opdfct2 = NULL, parmVec, traceVal,
                     silentVal = TRUE, matchCall){
    ## Controlling the warnings
    op1 <- options(warn = warnVal)
    on.exit(options(op1), add=TRUE)

    ## Calculating hessian
    if (is.null(opdfct2)) {hes <- TRUE} else {hes <- FALSE}

    ## Setting scaling parameters for optim()
    psVec <- abs(startVec)
    psVec[psVec < 1e-4] <- 1

    ## Derivatives are used
    {if (!is.null(opdfct1))
    {
        if (constrained)
        {
            nlsObj <- try(optim(startVec, opfct, opdfct1, hessian = hes, method = "L-BFGS-B",
            lower = lowerLimits, upper = upperLimits,
            control = list(maxit = maxIt, reltol = relTol, parscale = psVec)), silent = silentVal)
        } else {
            nlsObj <- try(optim(startVec, opfct, opdfct1, hessian = hes, method = optMethod,
            control = list(maxit = maxIt, reltol = relTol, parscale = psVec)), silent = silentVal)
        }
        options(warn = 0)

        if (!inherits(nlsObj, "try-error"))
        {
            nlsFit <- nlsObj
            nlsFit$convergence <- TRUE
        } else {
#            stop("Convergence failed")
            warning("Convergence failed. The model was not fitted!", call. = FALSE)

#            callDetail <- match.call()
#            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}
            return(list(call = matchCall, parNames = parmVec, startVal = startVec, convergence = FALSE))
        }
        if (!hes) {nlsFit$hessian <- opdfct2(nlsFit$par)}

    ## Derivatives are not used
    } else {

        if (constrained)
        {
            nlsObj <- try(optim(startVec, opfct, hessian = TRUE, method = "L-BFGS-B",
            lower = lowerLimits, upper = upperLimits,
            control = list(maxit = maxIt, parscale = psVec, reltol = relTol, trace = traceVal)), silent = silentVal)
            # parscale is needed for the example in methionine.Rd
        } else {
#            psVec <- abs(startVec)
#            psVec[psVec<1e-4] <- 1

            nlsObj <- try(optim(startVec, opfct, hessian = TRUE, method = optMethod,
            control = list(maxit = maxIt, reltol = relTol, parscale = psVec, trace = traceVal)), silent = silentVal)

#            nlsObj0 <- try(optim(startVec, opfct, method=optMethod,
#            control=list(maxit=maxIt, reltol=relTol, parscale=psVec)), silent=TRUE)
#            nlsObj <- try(optim(nlsObj0$par, opfct, hessian=TRUE, method=optMethod,
#            control=list(maxit=maxIt, reltol=relTol)), silent=TRUE)
        }
        options(warn = 0)

        if (!inherits(nlsObj, "try-error"))
        {
            nlsFit <- nlsObj
            nlsFit$convergence <- TRUE
        } else {  # to avoid an error if used in a loop
            if (errorMessage)
            {
                stop("Convergence failed")
            } else {
                warning("Convergence failed. The model was not fitted!", call. = FALSE)
            }

#            callDetail <- match.call()
#            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(LL.4())}
            return(list(call = matchCall, parNames = parmVec, startVal = startVec, convergence = FALSE))
        }
    }}

#    nlsFit$ofvalue <- nlsFit$value
    nlsFit$ovalue <- nlsFit$value  # used in the var-cov matrix ... check
#    nlsFit$value <- opfct(nlsFit$par, scaling = FALSE)  # used in the residual variance ... check
    nlsFit$value <- opfct(nlsFit$par)
    nlsFit$method <- "Parametric"

    ## Returning the fit
    return(nlsFit)
}
