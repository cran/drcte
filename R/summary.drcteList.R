summary.drcteList <- function(object,
                          robust = FALSE,
                          units = NULL, ...)
{

  # Edited 18/08/2022: Added robust and cluster robust SEs
  unName <- deparse(substitute(units))

  if(unName != "NULL"){
    tmp <- lapply(object$objVal, function(el) summary(el, robust = robust,
                                                    units = el$origData[[unName]])[[2]])
  } else {
    tmp <- lapply(object$objVal, function(el) summary(el, robust = robust)$coefficients)
  }

  tab <- do.call(rbind, tmp)

  row.names(tab) <- names(object$coefficients)
  sumObj <- list()
  sumObj$coefficients <- tab
  sumObj$varMat <- NA
  if(unName != "NULL") {
    sumObj$robust <- "cluster-robust sandwich"
  } else if(robust == T){
    sumObj$robust <- "sandwich"
  } else {
    sumObj$robust <- "no"
  }

  sumObj$resVar <- NULL
  sumObj$fctName <- "Parametric"
  sumObj$text <- "Separate fitting of several time-to-event curves"
  class(sumObj) <- c("summary.drcte", "summary.drc")
  return(sumObj)
}
