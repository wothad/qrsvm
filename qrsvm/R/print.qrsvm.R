
print.qrsvm <- function(x) {
  cat(paste0("Quantile Regression SVM. Cost is ", x$cost),
      fill = TRUE)
  cat(paste0("Estimated Quantile is ", x$tau), fill = TRUE)
}
