"estfun.drcte" <- function (x, ...)
#"estfun.drc" <- function (x, cvar = NULL, ...)
{
    ## Extending the matrix of derivatives to have one column per parameter (not only per model parameter)
    ##  only relevant in case several curves (with different parameters) were fitted
    xderiv1 <- x$deriv1
    indexMat0 <- x$indexMat2
    if (is.matrix(indexMat0))
    {
        indMat <- t(indexMat0)
    } else {
        indMat <- t(matrix(indexMat0, nrow = 1))
    }
    colnames(indMat) <- colnames(x$indexMat)
    curveID <- x$dataList[["curveid"]]

    ## Defining a helper function for making a wide version of the matrix
    # of parameter estimates
    ## One column per parameter, with 0s
    xderiv2Fct <- function(xderiv1, indMat, curveID)
    {
        xderiv2 <- xderiv1[, rep(1:ncol(xderiv1), apply(indMat, 1, function(x){length(unique(x))}))]
        cnInd <- colnames(indMat)
        for (i in 1:ncol(indMat))
        {
            xderiv2[curveID != cnInd[i], indMat[, i]] <- 0
        }
        xderiv2
    }

    if (identical(x$type, "event"))
    {  # not handling missing values

        doses <- NULL
        if(is.vector(x$dataList$dose) == F) doses <- 3:(length(x$dataList[["dose"]][1,])+1)
        # pr <- x[["fct"]]$"deriv1"(x$data[, c(1, doses)], t(x[["parmMat"]][, as.character(x$data[, 5])]))
        # pr2 <- x[["fct"]]$"deriv1"(x$data[, c(2, doses)], t(x[["parmMat"]][, as.character(x$data[, 5])]))
        # Verificare posizioni colonne. It depends whether we have additional predictors
        # With respect to timeBef and timeAf

        # Correction on 18/2/19 ################################
        # ids <- ifelse(is.null(doses), 5, 5 + length(x$dataList[["dose"]][1,]) - 2)
        ids <- ifelse(is.null(doses), 5, 5 + length(x$dataList[["dose"]][1,]) - 1)
        countCol <- ifelse(is.null(doses), 3, 3 + length(x$dataList[["dose"]][1,]) - 1)
        Ft1 <- predict(x, data.frame(x$data[, c(1, doses)], x$data[, ids]))
        Ft2 <- predict(x, data.frame(x$data[, c(2, doses)], x$data[, ids]))
        # 31/3/2022: needed this correction after modifying the output of predict
        Ft1 <- Ft1[,length(Ft1[1,])]
        Ft2 <- Ft2[,length(Ft2[1,])]
        Ft2[!is.finite(x$data[, 2])] <- 1
        diffF <- Ft2 - Ft1 #likelihood for individual seed
        diffF[diffF==0] <- 10e-6


        #Correction on 18/2/19 (went back) ################################
        DFt1.0 <- x[["fct"]]$"deriv1"(x$data[, c(1, doses)], t(x[["parmMat"]][, as.character(x$data[, ids])]))
        DFt2.0 <- x[["fct"]]$"deriv1"(x$data[, c(2, doses)], t(x[["parmMat"]][, as.character(x$data[, ids])]))

        # DFt1.0 <- x[["fct"]]$"deriv1"(x$data[, c(1, doses)], t(x[["parmMat"]][, as.character(x$data[, ids + 1])]))
        # DFt2.0 <- x[["fct"]]$"deriv1"(x$data[, c(2, doses)], t(x[["parmMat"]][, as.character(x$data[, ids + 1])]))

        DFt2.0[!is.finite(x$data[, 2]), ] <- 0
        diffDF.0 <- DFt2.0 - DFt1.0 # Derivative of likelihood funct

        #Modified on 12/2/19 - Went back on 18/2/19##########################################
        diffDF <- xderiv2Fct(diffDF.0, indMat, x$data[, ids])
        #diffDF <- xderiv2Fct(diffDF.0, indMat, x$data[, ids + 1])
        #print(diffDF); stop()
        rval <- (x$data[, countCol] / diffF) * diffDF #Derivative composite function D(log(f(x)))
        #print(diffDF)
        if (FALSE)
        {
        dataList <- x[["data"]]
        xderiv1 <- x[["fct"]]$"deriv1"(dataList[, 1], t(x[["parmMat"]][, as.character(dataList[, (ids - 1)])]))
        #print(xderiv1)
        xderiv2 <- xderiv2Fct(xderiv1, indMat, data[, (ids-1)])
#        print(xderiv2)
        resp <- c(0, dataList[, countCol])
        fittedVal <- c(0, predict(dataList[, 2]))
#        print(resp)
#        print(fittedVal)
#        fittedVal[length(fittedVal)] <- 1
#        fittedVal2 <- c(fittedVal[-1], 1)  # assuming data ordered according to time
#        rval0 <- diff(resp) / diff(fittedVal)
#        print(rval0)
##        lagDeriv1 <- apply(x$deriv1, 2, function(x){diff(x)})
#        lagDeriv1 <- apply(xderiv2, 2, function(x){diff(x)})
#        rval <- lagDeriv1 * rval0
#        rval <- xderiv2 * rval0
#        rval
        rval <- xderiv2 * (diff(resp) / diff(fittedVal))
        }
    }
#    ## Summing up according to specified clusters (the variable "cvar")
#    if (!is.null(cvar))
#    {
#        rval <- rowsum(rval, cvar)
#    }

#    print(rval)
    colnames(rval) <- names(coef(x))
#    colnames(rval) <- x$fct$names
    rval
}

"bread.drcte" <- function (x, ...)
{
#    if (identical(x$type, "binomial"))
#    {
#        breadMat <- vcov(x) * unlist(x$sumList[1])
#    }
    if (identical(x$type, "continuous"))
    {
#        breadMat <- summary(x)$cov.unscaled * unlist(x$sumList[1])
        breadMat <- vcov(x) / (summary(x)$rse[1]^2) * unlist(x$sumList[1])
    } else { ## Note: not checked for event time data!
        breadMat <- vcov(x) * unlist(x$sumList[1]) #Wrong to me
        #breadMat <- vcov(x) * sum(x$data[,3])
    }
    return(breadMat)
}

