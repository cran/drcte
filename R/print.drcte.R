"print.drcte" <- function(x, ..., digits = max(3, getOption("digits") - 3))
{
    object <- x

    classList <- class(object)
    cat(paste("\n", "A time-to-event model fitted with drcte", "\n", sep=""))
    if(object$fit$method == "NPMLE") fitMess = "NPMLE estimator for interval censored data (Fay and Shaw, 2010)"
    else if(object$fit$method == "KDE") fitMess = "Kernel density estimator for grouped data (Barreiro-Ures et al., 2019)"
    else fitMess <- "Parametric ML estimator"
    cat(paste("\n", "Method: ", fitMess, "\n", sep=""))


    ## Borrowing from print.lm
    ## Corrected on 6/3/2023 to produce a different message with drcteList objects
    if(inherits(object$"call", "call"))
      {cat("\nCall:\n", deparse(object$"call"), "\n\n", sep = "")
    }else {
      cat("\nCall:\n", "This is a list of multiple calls. See the 'separateFit' slot for more detail.", "\n\n", sep = "")
    }
    if (length(coef(object))>0)
    {
        cat("Coefficients:\n")
        print.default(format(coef(object), digits = digits), print.gap = 2, quote = FALSE)
    } else {
        cat("No coefficients\n")
    }
    cat("\n")

    invisible(object)
}
