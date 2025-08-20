"exponential" <- function(
fixed = c(NA, NA, NA), names = c("b", "d", "shift")){
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
        retVal <- parmMat[, 2] * (1 - exp(- parmMat[, 1] * (dose - parmMat[, 3])))
        # print(retVal)
        ifelse(retVal < 0, 0, retVal)
        # d * (1 - exp(- b * (x - shift)))
    }

    ## Defining the self starter function
    ssfct <- function(data){
          x <- data[, 1]
          y <- data[, 2]
          y <- y[x > 0]
          x <- x[x > 0]
          y <- y[!is.na(x)]
          x <- x[!is.na(x)]

          d <- max(y) * 1.01

          ## Linear regression on pseudo y values
          pseudoY <- log((d - y)/d)
          coefs <- coef( lm(pseudoY ~ x))
          b <- - coefs[2]
          shift <- coefs[1]/b
          # e <- exp(k/b)
          value <- c(b, ifelse(d>=1, 0.999, d), shift)
          # print(c(b,d,shift))
          return(value[notFixed])
    }


    ## Defining names
    names <- names[notFixed]

    ##Defining the first derivatives (in the parameters)
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        fctb <- deriv(~ d * (1 - exp(- b * (x - shift))), "b",
                     function.arg = c("b", "d", "shift", "x"))
        fctd <- deriv(~ d * (1 - exp(- b * (x - shift))), "d",
                     function.arg = c("b", "d", "shift", "x"))
        fcte <- deriv(~d * (1 - exp(- b * (x - shift))), "shift",
                     function.arg = c("b", "d", "shift", "x"))
        derb <- as.numeric( attr(fctb(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        derd <- as.numeric( attr(fctd(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        dere <- as.numeric( attr(fcte(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        ret <- cbind(derb, derd, dere)[, notFixed]
        # print(as.matrix(ret))
        # cbind(derb, derd, dere)[, notFixed]
        as.matrix(ret)
    }

    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        retVec <- deriv(~d * (1 - exp(- b * (x - shift))), "x",
                     function.arg = c("b", "d", "shift", "x"))
        retVec
    }


    ## Defining the ED function
    edfct <- function(parm, respl = 50, reference, type, ...)
    {
       respl <- respl
       parmVec[notFixed] <- parm
       if (type == "absolute")
       {
           tempVal <- log( (parmVec[2] - respl)/parmVec[2] )
           dVal <- 1
       } else {
           # respl <- respl * 100
           tempVal <- log( 1 - respl )
           dVal <- 0
       }
        EDp <-  - tempVal/parmVec[1]  + parmVec[3]
        EDder1 <- tempVal/(parmVec[1]^2)
        EDder2 <- dVal * -((1/parmVec[2] - (parmVec[2] - respl)/parmVec[2]^2)/((parmVec[2] - respl)/parmVec[2])/parmVec[1])
        EDder3 <- 1

        # D(expression(-log((d - respl)/d)/b + shift), "d")
        # D(expression(-log(1 - respl)/b + shift), "d")
        # D(expression(-log((d - respl)/d)/b + shift), "b")
        # D(expression(-log(1 - respl)/b + shift), "b")
        # D(expression(-log((d - respl)/d)/b + shift), "shift")
        # D(expression(-log(1 - respl)/b + shift), "shift")


        EDder <- c(EDder1, EDder2, EDder3)
        return(list(EDp, EDder[notFixed]))
    }

    ## Defining the inverse function
    invfct <- function(y, parm)
    {
        parmVec[notFixed] <- parm
        -log((parmVec[2] - y)/parmVec[2])/parmVec[1] + parmVec[3]
    }


    ## Returning the function with self starter and names
    returnList <-
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1,
         deriv2 = deriv2, derivx = derivx,
         edfct = edfct, inversion = invfct,
    name = "exponential",
    text = "exponential distribution of event times",
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "drcMean"
    invisible(returnList)
}

"exponentialSurv" <- function(
fixed = c(NA, NA, NA), names = c("b", "d", "shift")){

    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the model function
    fct <- function(dose, parm){
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        rate <- exp(-parmMat[, 1])
        retVal <- parmMat[, 2] * (1 - exp(- rate * (dose - parmMat[, 3])))
        ifelse(retVal < 0, 0, retVal)
        # d * (1 - exp(- b * (x - shift)))
    }

    ## Defining the self starter function
    ssfct <- function(data){
          x <- data[, 1]
          y <- data[, 2]
          y <- y[x > 0]
          x <- x[x > 0]

          d <- max(y) * 1.01

          ## Linear regression on pseudo y values
          pseudoY <- log((d - y)/d)
          coefs <- coef( lm(pseudoY ~ x))
          rate <- - coefs[2]
          shift <- coefs[1]/rate
          # e <- exp(k/b)
          b <- -log(rate)
          value <- c(b, ifelse(d>=1, 0.999, d), shift)
          return(value[notFixed])
    }


    ## Defining names
    names <- names[notFixed]

    ##Defining the first derivatives (in the parameters)
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        fctb <- deriv(~ d * (1 - exp(-exp(-b) * (x - shift))), "b",
                     function.arg = c("b", "d", "shift", "x"))
        fctd <- deriv(~ d * (1 - exp(-exp(-b) * (x - shift))), "d",
                     function.arg = c("b", "d", "shift", "x"))
        fcte <- deriv(~d * (1 - exp(-exp(-b) * (x - shift))), "shift",
                     function.arg = c("b", "d", "shift", "x"))
        derb <- as.numeric( attr(fctb(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        derd <- as.numeric( attr(fctd(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        dere <- as.numeric( attr(fcte(parmMat[, 1], parmMat[, 2], parmMat[, 3], dose), "gradient") )
        ret <- cbind(derb, derd, dere)[, notFixed]
        # print(as.matrix(ret))
        # cbind(derb, derd, dere)[, notFixed]
        as.matrix(ret)
    }

    ## Defining the second derivative (in the parameters)
    deriv2 <- NULL

    ## Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        retVec <- deriv(~d * (1 - exp(-exp(-b) * (x - shift))), "x",
                     function.arg = c("b", "d", "shift", "x"))
        retVec
    }

    ## Defining the ED function
    edfct <- function(parm, respl = 50, reference, type, ...)
    {
       respl <- respl
       parmVec[notFixed] <- parm
       if (type == "absolute")
       {
           tempVal <- log( (parmVec[2] - respl)/parmVec[2] )
           dVal <- 1
       } else {
           # respl <- respl * 100
           tempVal <- log( 1 - respl )
           dVal <- 0
       }
        EDp <-  - tempVal/exp(-parmVec[1])  + parmVec[3]
        EDder1 <- -(tempVal * exp(-parmVec[1])/exp(-parmVec[1])^2)
        EDder2 <- dVal * -((1/parmVec[2] - (parmVec[2] - respl)/parmVec[2]^2)/((parmVec[2] - respl)/parmVec[2])/exp(-parmVec[1]))
        EDder3 <- 1

        # D(expression(-log((d - respl)/d)/exp(-b) + shift), "d")
        # D(expression(-log(1 - respl)/exp(-b) + shift), "d")
        # D(expression(-log((d - respl)/d)/exp(-b) + shift), "b")
        # D(expression(-log(1 - respl)/exp(-b) + shift), "b")
        # D(expression(-log((d - respl)/d)/exp(-b) + shift), "shift")
        # D(expression(-log(1 - respl)/exp(-b) + shift), "shift")


        EDder <- c(EDder1, EDder2, EDder3)
        return(list(EDp, EDder[notFixed]))
    }

    ## Defining the inverse function
    invfct <- function(y, parm)
    {
        parmVec[notFixed] <- parm
        -log((parmVec[2] - y)/parmVec[2])/exp(-parmVec[1]) + parmVec[3]
    }

    ## Defining the density function
    densy <- function(y, parm)
    {
        # parmVec[notFixed] <- parm
        # -log((parmVec[2] - y)/parmVec[2])/exp(-parmVec[1]) + parmVec[3]
    }


    ## Returning the function with self starter and names
    returnList <-
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1,
         deriv2 = deriv2, derivx = derivx,
         edfct = edfct, inversion = invfct,
    name = "exponentialSurv",
    text = "exponential distribution of event times",
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "drcMean"
    invisible(returnList)
}
