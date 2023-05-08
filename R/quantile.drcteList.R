##########################################################################
quantile.drcteList <- function(x, probs, restricted = FALSE,
                           interval = c("none", "delta", "boot"),
                           clevel = NULL, level = ifelse(!(interval == "none"), 0.95, NULL),
                           bound = TRUE, od = FALSE,
                           vcov. = vcov, robust = FALSE, units = NULL,
                           display = FALSE, rate = F, B = 999, ...){

  # Edited 18/08/2022: Added robust and cluster robust SEs
  unName <- deparse(substitute(units))

  quantileList <- function(el, probs, rate, display, robust, unName){

    if(unName != "NULL"){
    EDlist <- try(quantile(el, probs = probs,
                           rate = rate, display = F,
                           robust = robust,
                           units = el$origData[[unName]]),
                  silent = T)
    } else {

    EDlist <- try(quantile(el, probs = probs,
                         rate = rate, display = F,
                         robust = robust),
                silent = T)
    }
    if(any(class(EDlist) == "try-error") == T) {
    estim <- ifelse(rate == T, 0, NA)
    estim <- rep(estim, length(probs))
    SE <- rep(NA, length(probs))
    retDF <- data.frame(estim,
                        rLab = SE)
    colnames(retDF) <- c("Estimate", "SE")
    row.names(retDF) <- paste("1:", probs*100, "%", sep = "")
    EDlist <- retDF
    }
    EDlist
  }

  ret <- lapply(x$separateFit, quantileList, probs = probs,
              rate = rate, display = F, robust=robust, unName = unName)
  # Edited on 7/3/2023. In order to avoid problems on row.names
  # when probs is of length 1
  if(length(probs) == 1){
      retNames <- do.call(rbind, lapply(ret, row.names))
      retFin <- do.call(rbind, ret)
      row.names(retFin) <- paste(row.names(retFin), retNames, sep = "")
    } else {
      retFin <- do.call(rbind, ret)
    }
  row.names(retFin) <- sub(".1:", ":", row.names(retFin), fixed = T)
  row.names(retFin) <- sub(".:", ":", row.names(retFin), fixed = T)
  return(retFin)
}
