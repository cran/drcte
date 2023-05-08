"coef.drcte" <-
function(object, ...)
{
    if (!is.null(object$"coefficients"))
    {
        return(object$"coefficients")
    } else {
        retVec <- c()
        # names(retVec) <- object$parNames[[1]]
        return(retVec)
    }
}
