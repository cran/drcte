"drmConvertParm" <- 
function(startVec, startMat, factor1, colList)
{
    startMat2 <- startMat
    if (length(unique(factor1)) == 1) {return(startVec)}
    mmat <- model.matrix(~factor(factor1) - 1)
    pm <- list()
    for (i in 1:length(colList))
    {
        clElt <- colList[[i]]
        ncclElt <- dim(clElt)[2]   
            
        indVec <- !is.na(startMat2[, i, drop = FALSE])
        indVal <- min(c(sum(indVec), dim(clElt)[2]))
        indVec2 <- (1:ncclElt)[indVec]
        if (length(indVec2) > ncclElt) {indVec2 <- 1:ncclElt}
        pm[[i]] <- (ginv(t(clElt)%*%clElt)%*%t(clElt))[1:indVal, ,drop = FALSE]%*%mmat[,indVec]%*%startMat2[indVec, i, drop = FALSE]
    }  
    tempVec <- unlist(pm)
    tempVec <- tempVec[!is.na(tempVec)]

    ## Checking whether the intercept column has been removed
    indVec3 <- ( abs(tempVec)<1e-10 )
    if (any(indVec3)) 
    {
        tempVec <- startVec
    }

    return(tempVec)
}

"drmParNames" <- 
function(numNames, parNames, collapseList2, repStr1 = "factor(pmodels[, i])", repStr2 = "factor(assayNo)")
{
    ## Retrieving names for parameters
    parmVecList <- list()
    for (i in 1:numNames)
    {
        colNames1 <- colnames(collapseList2[[i]])
        if (is.null(colNames1)) 
        {
            parmVecList[[i]] <- paste(parNames[i], "(Intercept)", sep = ":")
        } else { 
            parmVecList[[i]] <- paste(parNames[i], colNames1, sep = ":")
        }        
        parmVecList[[i]] <- (parmVecList[[i]])[1:ncol(collapseList2[[i]])]  # min(maxParm[i], length(colNames1))]
    }
    parmVec <- unlist(parmVecList)
        
    parmVec2 <- parmVec
#    print(parmVec2)
    for (i in 1:length(parmVec))
    {
        pos <- regexpr(repStr1, parmVec[i], fixed = TRUE)
        if (pos > 0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i], 1, pos-1), substring(parmVec[i], pos + 20), sep = "")
        }
            
        pos <- regexpr(repStr2, parmVec[i], fixed = TRUE)
        if (pos > 0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i], 1, pos-1), substring(parmVec[i], pos + 15), sep = "")
        }
    }
    return(drmPNsplit(parmVec2, ":"))
}

"vec2mat" <- function(fct, no)
{
    parName <- names(formals(fct)[no])
    if (is.na(parName)) {stop("Argument number does not exist")}
    parName1 <- paste(parName, "[", sep = "")
    parName2 <- paste(parName, "[,", sep = "")

#    bodyStr <- as.character(body(fct))
#    if (bodyStr[1] == "{") {bodyStr <- bodyStr[-1]}  #  else {bodyStr <- bodyStr[1]}
    bodyStr <- deparse(body(fct))
    if (bodyStr[1] == "{") {bodyStr <- paste(head(tail(bodyStr, -1), -1), collapse = "")}  #  else {bodyStr <- bodyStr[1]}
    
    lenbs <- length(bodyStr)
    bsList <- list()
#    options(warn = -1)
    for (i in 1:lenbs)
    {
        tempText <- gsub(parName1, parName2, bodyStr[i], fixed = TRUE)
        bsList[[i]] <- paste(tempText, ";", sep = "")
    }
#    options(warn = 0)
    bodyStr2 <- paste("{", paste(as.vector(unlist(bsList)), collapse = ""), "}")
    bodyStr3 <- parse(text = bodyStr2)

    newfct <- fct
    body(newfct) <- bodyStr3  # bodyStr2

    return(list(newfct, bodyStr2, parName))
}

"nParm" <- function(bodyStr)
{
    gregObj <- gregexpr("(\\[,){1}[[:digit:]]+\\]{1}", bodyStr)[[1]]
    posVec <- gregObj
    lenVec <- attr(gregObj, "match.length")

    lenpv <- length(posVec)
    numVec <- rep(0, lenpv)
    for (i in 1:lenpv)
    {
        numVec[i] <- as.numeric(substr(bodyStr, posVec[i] + 2, posVec[i] + lenVec[i] - 2))
    }
    return(length(unique(numVec)))
}

"fParm" <- function(fct, no, fixed)
{
    v2m <- vec2mat(fct, no)
    if (all(is.na(fixed))) {return(v2m[[1]])}    
    
    bodyStr <- v2m[[2]]

    fStr <- paste("{ f <- c(", paste(fixed[!is.na(fixed)], collapse = ", "), ");")
    bodyStr <- paste(fStr, substr(bodyStr, 2, nchar(bodyStr)))
    parName <- v2m[[3]]

    gregObj <- gregexpr("(\\[,){1}[[:digit:]]+\\]{1}", bodyStr)[[1]]
    posVec <- gregObj
    lenVec <- attr(gregObj, "match.length")

    lenpv <- length(posVec)
    numVec <- rep(0, lenpv)
    
    realPos <- rep(NA, length(fixed))
    realPos[is.na(fixed)] <- 1:sum(is.na(fixed)) 
    
    fixPos <- rep(NA, length(fixed))
    fixPos[!is.na(fixed)] <- 1:sum(!is.na(fixed))    
    
#    options(warn = -1)    
    for (i in 1:lenpv)
    {
        numVec[i] <- as.numeric(substr(bodyStr, posVec[i] + 2, posVec[i] + lenVec[i] - 2))
        
        if (is.na(realPos[numVec[i]]))
        {
            inStr0 <- paste("f[", as.character(fixPos[numVec[i]]), "]", sep = "")
           
            lenBl <- nchar(parName) + lenVec[numVec[i]] - nchar(inStr0)
            inStr1 <- paste(rep(" ", nchar(parName) + lenVec[numVec[i]] - nchar(inStr0)), collapse = "")
            inStr2 <- paste(inStr0, inStr1, sep = "")
  
            substr(bodyStr, posVec[i] - nchar(parName), posVec[i] + lenVec[i] - 1) <- inStr2



#           inStr0 <- as.character(fixed[numVec[i]])
#           
#           lenBl <- nchar(parName) + lenVec[numVec[i]] - nchar(inStr0)
#           if (lenBl > 0) 
#           {
#               inStr1 <- paste(rep(" ", lenBl), collapse = "")
#               inStr2 <- paste(inStr0, inStr1, sep = "")
#           } else {
#               inStr2 <- inStr0
#           }         
#           substr(bodyStr, posVec[i] - nchar(parName), posVec[i] + lenVec[i] - 1) <- inStr2

        } else {
        
            inStr3 <- as.character(realPos[numVec[i]])
            
            ## In case the number changes from 10 to 9, 100 to 99 and so on
            if (nchar(inStr3) < nchar(as.character(numVec[i])) )
            {
                numSpaces <- nchar(as.character(numVec[i])) - nchar(inStr3)
                tempStr <- paste(rep(" ", numSpaces), collapse = "")
            
                inStr3 <- paste(tempStr, inStr3, sep = "")
            }
            substr(bodyStr, posVec[i] + 2, posVec[i] + lenVec[i] - 3) <- inStr3           
        }
    }

    bodyStr2 <- parse(text = bodyStr)
    newfct <- fct
    body(newfct) <- bodyStr2

    return(newfct)
}


fct2list <- function(fct, no)
{
    v2m <- vec2mat(fct, no) 
    list(v2m[[1]], NULL, letters[1:nParm(v2m[[2]])])
}

modelFunction <- function(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, retFct, 
                          doseScaling, respScaling, isFinite, pshifts = NULL)
{
    if (!is.null(retFct))
    {
        drcFct <- retFct(doseScaling, respScaling)
    }
    drcFct1 <- function(dose, parm)
    {
        parmVal <- parm2mat(parm)
#        print(c(dim(pshifts), dim(parmVal)))
        if ((!is.null(pshifts)) & all(dim(pshifts) == dim(parmVal))) 
        {
            parmVal <- parmVal + pshifts
        }     
#        drcFct(dose, (parm2mat(parm))[isFinite, , drop = FALSE])
        drcFct(dose, parmVal[isFinite, , drop = FALSE])
    }

    if (is.null(cm))
    {
        multCurves <- function(dose, parm)
        {
           drcFct1(dose, parm)
        }
    } else {  # not adapting to scaling (not using drcFct1)!!!
        iv <- isFinite & (assayNoOld == cm)
        niv <- !iv
        fctEval <- rep(0, length(dose))

        multCurves <- function(dose, parm)
        {
            parmVal <- (parm2mat(parm))[isFinite, , drop = FALSE]
#            print(c(dim(pweights), dim(parmVal)))
            if ((!is.null(pshifts)) & all(dim(pshifts) == dim(parmVal))) 
            {
                parmVal <- parmVal + pshifts
            }
            fctEval[iv] <- parmVal[iv, upperPos, drop = FALSE]
            fctEval[niv] <- drcFct(dose[niv], parmVal[niv, , drop = FALSE])

            fctEval
        }
    }
    
    multCurves
}

"drmPNsplit" <- 
function(parmVec, sep) 
{
    lenPV <- length(parmVec)
    parmVecA <- rep(0, lenPV)
    parmVecB <- rep(0, lenPV)

    splitList <- strsplit(parmVec, sep, fixed = TRUE)
    for (i in 1:lenPV)
    {
        parmVecA[i] <- splitList[[i]][1]
        
        lenSL <- length(splitList[[i]])
        parmVecB[i] <- paste(splitList[[i]][2:lenSL], collapse = "")  # 'paste' is needed in case several ":" occur
    }
    return(list(parmVec, parmVecA, parmVecB))
}
