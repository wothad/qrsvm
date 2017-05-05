
summary.qrsvm <- function(x) {
  for (i in c(4, 5, 6, 8, 9, 10, 11, 12)) {
    cat(paste(names(x[[i]]), " is ", x[[i]]), fill = TRUE)
  }
}
