#' Predict ann Object oc class "qrsvm"
#'
#' @param model An object of class "qrsvm"
#' @param newdata The predictors of the predictable data in an n X m Matrix
#' @return A numeric vector of predicted values
#' @import kernlab
#' @export
predict.qrsvm <- function(model, newdata) {

  library(kernlab)
  xold <- model[[2]]
  alpha <- model[[1]]
  kern <- model[[3]]
  b <- model[[6]]

  if (ncol(newdata) != ncol(xold)) {
    cat("Newdata has different number of columns than xtrain please check consistency!",
        fill = TRUE)
  }

  pred <- kernelMult(kern, newdata, xold, alpha)
  pred <- pred + b
  return(pred)
}
