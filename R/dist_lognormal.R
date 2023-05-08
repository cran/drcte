"lognormal" <- function(
fixed = c(NA, NA, NA), names = c("b", "d", "e"))
{

    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the model function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        #parmMat[, 2]/(1 + exp(- parmMat[, 1]*(log(dose + 0.000001) - log(parmMat[, 3]))))
        parmMat[, 2] * pnorm(parmMat[, 1]*(log(dose + 0.000001) - log(parmMat[, 3])))
    }

    ## Defining the self starter function
    ssfct <- function(data){
          x <- data[, 1]
          y <- data[, 2]
          y <- y[x > 0]
          x <- x[x > 0]


          d <- max(y) * 1.01

          ## Linear regression on pseudo y values
          pseudoY <- log((d - y)/(y + 0.000001))
          coefs <- coef( lm(pseudoY ~ log(x)))
          b <- - coefs[2]
          k <- coefs[1];
          e <- exp(k/b)
          value <- c(b, ifelse(d>=1, 0.999, d), e)

          return(value[notFixed])
    }


    ## Defining names
    names <- names[notFixed]

    ##Defining the first derivatives (in the parameters)
    deriv1 <- function(dose, parm)
    {
        # ~d * pnorm(b * (log(x + 0.000001) - log(e)))
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        fctb <- deriv(~d * pnorm(b * (log(x + 0.000001) - log(e))), "b",
                     function.arg = c("b", "d", "e", "x"))
        fctd <- deriv(~d * pnorm(b * (log(x + 0.000001) - log(e))), "d",
                     function.arg = c("b", "d", "e", "x"))
        fcte <- deriv(~d * pnorm(b * (log(x + 0.000001) - log(e))), "e",
                     function.arg = c("b", "d", "e", "x"))
        derb <- as.numeric( attr(fctb(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        derd <- as.numeric( attr(fctd(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        dere <- as.numeric( attr(fcte(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )

        cbind(derb, derd, dere)[, notFixed]
    }

    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        retVec <- deriv(~d * pnorm(b * (log(x + 0.000001) - log(e))), "x",
                     function.arg = c("b", "d", "e", "x"))
        retVec
    }


    ## Defining the ED function
    edfct <- function(parm, respl = 50, reference, type, ...)
    {
       respl <- respl
       parmVec[notFixed] <- parm
       if (type == "absolute")
       {
           tempVal <- log( (parmVec[2] - respl)/respl )
           dVal <- 1
       } else {
           # respl <- respl * 100
           tempVal <- log( (1 - respl)/respl )
           dVal <- 0
       }
        EDp <- exp( - tempVal/parmVec[1] + log(parmVec[3]))
        EDder1 <- - EDp * 1/(parmVec[1]^2) * tempVal
        EDder2 <- dVal * (EDp * (- 1/parmVec[1] * (100/respl/((100 * parmVec[2] - respl)/respl))))
        EDder3 <- EDp * 1/parmVec[3]

        # D(expression(exp(- 1/b * log((100 - p)/p) + log(e))), "d")
        # D(expression(exp(- 1/b * log((100*d - p)/p) + log(e))), "d")

        EDder <- c(EDder1, EDder2, EDder3)
        return(list(EDp, EDder[notFixed]))
    }

    ## Defining the inverse function
    invfct <- function(y, parm)
    {
        parmVec[notFixed] <- parm
        exp(log((parmVec[2] - y)/(y))/parmVec[1] + log(parmVec[3]))
    }


    ## Returning the function with self starter and names
    returnList <-
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, inversion = invfct,
    name = "lognormal",
    text = "Log-normal distribution of event times",
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "drcMean"
    invisible(returnList)
}

"lognormalSurv" <- function(
fixed = c(NA, NA, NA), names = c("b", "d", "e"))
{

    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the model function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        parmMat[, 2]/(1 + exp(- (log(dose + 0.000001) - parmMat[, 3])/exp(parmMat[, 1])))
    }

    ## Defining the self starter function
    ssfct <- function(data){
          x <- data[, 1]
          y <- data[, 2]
          y <- y[x > 0]
          x <- x[x > 0]


          d <- max(y) * 1.01

          ## Linear regression on pseudo y values
          pseudoY <- log((d - y)/(y + 0.000001))
          # print(pseudoY)
          coefs <- coef( lm(pseudoY ~ log(x)))
          b <- - 1/coefs[2]
          k <- coefs[1];
          e <- k * b
          value <- c(log(b), ifelse(d>=1, 0.999, d), e)

          return(value[notFixed])
    }


    ## Defining names
    names <- names[notFixed]

    ##Defining the first derivatives (in the parameters)
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        meanFct <- ~d/(1 + exp(-(log(x) - e)/exp(b)))
        fctb <- deriv(meanFct, "b",
                     function.arg = c("b", "d", "e", "x"))
        fctd <- deriv(meanFct, "d",
                     function.arg = c("b", "d", "e", "x"))
        fcte <- deriv(meanFct, "e",
                     function.arg = c("b", "d", "e", "x"))
        # print(fctb(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose)); stop()
        derb <- as.numeric( attr(fctb(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        derd <- as.numeric( attr(fctd(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        dere <- as.numeric( attr(fcte(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )

        cbind(derb, derd, dere)[, notFixed]
    }

    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        meanFct <- ~ d/(1 + exp(-(log(x) - e)/exp(b)))
        retVec <- deriv(meanFct, "x",
                     function.arg = c("b", "d", "e", "x"))
        retVec
    }


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
       parmVec[notFixed] <- parm
       if (type == "absolute")
       {
           tempVal <- log( (100 * parmVec[2] - respl)/respl )
           dVal <- 1
       } else {
           tempVal <- log( (100 - respl)/respl )
           dVal <- 0
       }
        EDp <- exp(tempVal * exp(parmVec[1]) + parmVec[3])
        EDder1 <- - EDp * exp(parmVec[1]) * tempVal
        EDder2 <- dVal * (EDp * (exp(parmVec[1]) * (100/respl/((100 * parmVec[2] - respl)/respl))))
        EDder3 <- EDp

        # D(expression(exp(exp(b) * log((100 - p)/p) + e)), "e")
        # D(expression(exp(exp(b) * log((100*d - p)/p) + e)), "e")

        EDder <- c(EDder1, EDder2, EDder3)
        return(list(EDp, EDder[notFixed]))
    }


    ## Defining the inverse function
    invfct <- function(y, parm)
    {
        parmVec[notFixed] <- parm

        exp(log(((parmVec[3] - parmVec[2])/(y - parmVec[2])) - 1)/parmVec[1] + log(parmVec[3]))
    }

    ## Returning the function with self starter and names
    returnList <-
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, inversion = invfct,
    # name = fctName,
    text = "Log-logistic distribution of germination times",
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "drcMean"
    invisible(returnList)
}
