"vcov.drcte" <-
function(object, ..., corr = FALSE)
{
  ## Defining function for calculating variance-covariance matrix
  vcMat <- solve(object$fit$hessian)
  if(length(object$dataList$names$rNames) > 1) vcMat[vcMat < 1E-8] <- 0
  vcMat
}