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
        # parmMat[, 2]/(1 + exp(- parmMat[, 1]*(log(dose + 0.000001) - log(parmMat[, 3]))))
        parmMat[, 2] * pnorm(parmMat[, 1] * (log(dose + 0.000001) - log(parmMat[, 3])))

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
    edfct <- function(parm, respl, reference, type, ...)
    {
      edfct.abs <- function(p1, p2, p3, respl){
        tempVal <- qnorm(respl/p2)
        exp( tempVal * 1/p1 + log(p3) )
      }
      edfct.rel <- function(p1, p2, p3, respl){
         tempVal <- qnorm(respl)
         exp( tempVal * 1/p1 + log(p3) )
      }
      respl <- respl
      parmVec[notFixed] <- parm
      if (type == "absolute")
       {
         EDp <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.1 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.abs(parmVec[1] +  10e-7, parmVec[2], parmVec[3], respl)
         EDder1 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.abs(parmVec[1], parmVec[2] +  10e-7, parmVec[3], respl)
         EDder2 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3] +  10e-7, respl)
         EDder3 <- (d1.2 - d1.1)/10e-7
       } else {
         EDp <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.1 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.rel(parmVec[1] +  10e-7, parmVec[2], parmVec[3], respl)
         EDder1 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.rel(parmVec[1], parmVec[2] +  10e-7, parmVec[3], respl)
         EDder2 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3] +  10e-7, respl)
         EDder3 <- (d1.2 - d1.1)/10e-7
       }

        EDder <- c(EDder1, EDder2, EDder3)
        return(list(EDp, EDder[notFixed]))
    }

    ## Defining the inverse function
    invfct <- function(y, parm)
    {
        parmVec[notFixed] <- parm
        # exp(log((parmVec[2] - y)/(y))/parmVec[1] + log(parmVec[3]))
        exp(1/parmVec[1]*qnorm(y/parmVec[2]) + log(parmVec[3]))

    }


    ## Returning the function with self starter and names
    returnList <-
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, inversion = invfct,
    name = "lognormal.2",
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
    if (!(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the model function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        if(!is.na(fixed[2])){
          d <- 1
        } else {
          logitd <- parmMat[,2] # ifelse(parmMat[,2] > 20, 20, parmMat[,2])
          d <- exp(logitd)/(1 + exp(logitd))
          # d <- ifelse(is.nan(d), 1, d)
        }

        d * pnorm((log(dose + 0.000001) - parmMat[, 3])/exp(parmMat[, 1]))
        # d * pnorm(parmMat[, 1]*(log(dose + 0.000001) - log(parmMat[, 3])))
    }

    ## Defining the self starter function
    ssfct <- function(data){
          x <- data[, 1]
          y <- data[, 2]
          y <- y[x > 0]
          x <- x[x > 0]
          # print(data)
          d <- max(y) * 1.01

          ## Linear regression on pseudo y values
          pseudoY <- log((d - y)/(y + 0.000001))
          coefs <- coef( lm(pseudoY ~ log(x)))
          b <- - 1/coefs[2]
          k <- coefs[1]
          e <- k * b
          d <- ifelse(d > 1, 0.99, d)
          d <- log(d/(1 - d)) # exp(d)/(1 + exp(d))
          # print(d); print(exp(d)/(1 + exp(d)))
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
        meanFct <- ~ exp(d)/(1 + exp(d)) * pnorm((log(x + 0.000001) - e)/b)
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
        meanFct <- ~ exp(d)/(1 + exp(d)) * pnorm((log(x + 0.000001) - e)/b)
        retVec <- deriv(meanFct, "x",
                     function.arg = c("b", "d", "e", "x"))
        retVec
    }


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
      edfct.abs <- function(p1, p2, p3, respl){
        if(!is.na(fixed[2])) d <- 1 else d <- exp(p2)/(1 + exp(p2))
         tempVal <- qnorm(respl/d)
         exp( tempVal * exp(p1) + p3 )
      }
      edfct.rel <- function(p1, p2, p3, respl){
         # d <- exp(p2)/(1 + exp(p2))
         tempVal <- qnorm(respl)
         exp( tempVal * exp(p1) + p3 )
      }
      respl <- respl
      parmVec[notFixed] <- parm
      if (type == "absolute")
       {
         EDp <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.1 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.abs(parmVec[1] +  10e-7, parmVec[2], parmVec[3], respl)
         EDder1 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.abs(parmVec[1], parmVec[2] +  10e-7, parmVec[3], respl)
         EDder2 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.abs(parmVec[1], parmVec[2], parmVec[3] +  10e-7, respl)
         EDder3 <- (d1.2 - d1.1)/10e-7
       } else {
         EDp <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.1 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.rel(parmVec[1] +  10e-7, parmVec[2], parmVec[3], respl)
         EDder1 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.rel(parmVec[1], parmVec[2] +  10e-7, parmVec[3], respl)
         EDder2 <- (d1.2 - d1.1)/10e-7
         d1.1 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3], respl)
         d1.2 <- edfct.rel(parmVec[1], parmVec[2], parmVec[3] +  10e-7, respl)
         EDder3 <- (d1.2 - d1.1)/10e-7
       }

        EDder <- c(EDder1, EDder2, EDder3)
        print(EDder)
        return(list(EDp, EDder[notFixed]))
    }


    ## Defining the inverse function
    invfct <- function(y, parm)
    {
        parmVec[notFixed] <- parm
        exp(exp(parmVec[1])*qnorm(y/parmVec[2]) + parmVec[3])
        # exp(log(((parmVec[3] - parmVec[2])/(y - parmVec[2])) - 1)/parmVec[1] + log(parmVec[3]))
    }
    linkFct <- function(){
      link1 <- "1/exp(b)"
      link2 <- "exp(d)/(1 + exp(d))"
      link3 <- "exp(e)"
      link <- c(link1, link2, link3)
      names(link) <- c("b", "d", "e")
      return(link[notFixed])
    }
    ## Returning the function with self starter and names
    returnList <-
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx,
    edfct = edfct, inversion = invfct, linkFct = linkFct,
    name = "lognormal",
    text = "Log-normal distribution of event times",
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "drcMean"
    invisible(returnList)
}
