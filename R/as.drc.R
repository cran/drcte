as.drc <- function(obj){
  # This function transform an lm object into a drc object
  # All slots are ok (19/3/2022)
  nlsFit <- list()
  nlsFit$par <- coef(obj)
  nlsFit$value <- deviance(obj) # as.numeric(logLik(obj))
  nlsFit$counts <- NULL
  nlsFit$convergence <- TRUE
  nlsFit$message <- NULL

  tmp <- vcov(obj)
  tmp[tmp==0] <- 1E-6
  nlsFit$hessian <- matrix(solve(tmp)) * 2 * summary(obj)$sigma^2
  nlsFit$method <- "Parametric"

  # varParm
  varParm <- NULL

  ## curve function
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

      curvePts <- predict(obj, newdata = data.frame(dose))
      curvePts <- matrix(curvePts, lenPts, 1)
      return(curvePts)
    }
    return(plotFct)
  }
  plotFct <- pfFct(obj)

  # summary
  sumVec <- c(NA, NA, NA, deviance(obj), obj$df.residuals, length(obj$residuals))

  # start
  startVec <- NULL

  # parmVec
  parmVec <- list(names(obj$coefficients),
                  names(obj$coefficients),
                  names(obj$coefficients))

  # predres
  diagMat <- cbind(fitted(obj), residuals(obj))
  callDetail <- obj$call
  origData <- eval(obj$call$data)
  dataSet <- data.frame(obj$model, curveid = rep(1, length(residuals(obj))),
                        orig.curveid = rep(1, length(residuals(obj))),
                        weights = rep(1, length(residuals(obj))))
  parmMatOut <- matrix(coef(obj), length(coef(obj)), 1)
  fct <- NULL
  robust <- NULL
  estMethod <- NULL
  df.residual <- obj$df.residual

  sumList <- list(lenData = length(obj$residuals),
                  alternative = NULL,
                  df.residual = obj$df.residuals)
  scaleFct <- NULL
  fixedParm <- coef(obj)
  pmFct <- function(fixedParm = coef(obj))
  {
    parmMat <- matrix(fixedParm, 1, length(fixedParm))
    return(parmMat)
  }
  type <- "continuous"
  indexMat <- matrix(c(1:length(coef(obj))), length(coef(obj)), 1)
  logDose <- NULL
  cm <- NULL
  deriv1 <- model.matrix(obj)
  anName <- NULL
  wVec <- rep(1, length(residuals(obj)))

  dataList <- list()
  if(length(obj$model[1,]) == 1) {
    dataList$dose <- rep(1, length(obj$model[,1]))
  } else {
    dataList$dose <- obj$model[,-1]
  }
  dataList$origResp <- obj$model[,1]
  dataList$weights <- NA
  dataList$curveid <- rep(1, length(residuals(obj)))
  dataList$plotid <- rep(1, length(residuals(obj)))
  dataList$resp <- obj$model[,1]
  dataList$names <- list()
  dataList$names$dName <- names(obj$model)[-1]
  dataList$names$orName <- names(obj$model)[1]
  dataList$names$wName <- NULL
  dataList$names$cNames <- NULL
  dataList$names$rNames <- 1
  dataList$names$adVarNames <- NULL
  coefVec <- coef(obj)
  bcVec <- NULL
  newObj <- list(NULL, nlsFit, list(plotFct, NULL), sumVec,
                     startVec, parmVec, diagMat, callDetail,
                     dataSet, parmMatOut, fct, robust, estMethod, df.residual,
                     sumList, NULL, pmFct, pfFct, type, indexMat, logDose, cm, deriv1, #[isFinite,],
                     anName, origData, wVec,
                     dataList, coefVec, bcVec)


  names(newObj) <- c("varParm", "fit", "curve", "summary",
                         "start", "parNames", "predres", "call",
                         "data", "parmMat", "fct", "robust", "estMethod", "df.residual",
                         "sumList", "scaleFct", "pmFct", "pfFct", "type", "indexMat", "logDose", "cm", "deriv1",
                         "curveVarNam", "origData", "weights",
                         "dataList", "coefficients", "boxcox")
  class(newObj) <- c("drc")
  invisible(newObj)
}
