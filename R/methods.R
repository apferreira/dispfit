#' @export
summary.dispfit <- function(x, ...) {
  stopifnot(inherits(x, "dispfit"))
  res <- list()
  res$parameters <- x$distribution.parameters
  res$selection <- x$distribution.selection
  class(res) <- "summary.dispfit"
  res
}

#' @export
print.summary.dispfit <- function(x, ...)
{
  message("\nDistribution selection")
  print(x$selection)
  message("\nDistribution parameters")
  tmp <- x$parameters
  tmp[tmp == Inf] <- "Infinite value"
  #tmp[is.na(tmp)] <- "in progress"
  print(tmp)
}

#' @export
print.dispfit <- function(x, ...)
{
  message("\nDistribution selection")
  print(x$distribution.selection)
  message("\nDistribution parameters")
  tmp <- x$distribution.parameters
  tmp[tmp == Inf] <- "Infinite value"
  #tmp[is.na(tmp)] <- "in progress"
  print(tmp)
}

#' @export
AIC.dispfit <- function(x, ...)
{
  print(x$distribution.selection[1:7])
}
