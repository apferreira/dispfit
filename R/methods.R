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
  print("Distribution selection")
  print(x$selection)
  print("Distribution parameters")
  tmp <- x$parameters
  tmp[tmp == Inf] <- "Infinite valuesss"
  tmp[is.na(tmp)] <- "in progress"
  print(tmp)
}

#' @export
print.dispfit <- function(x, ...)
{
  print("Distribution selection")
  print(x$distribution.selection)
  print("Distribution parameters")
  tmp <- x$distribution.parameters
  tmp[tmp == Inf] <- "Infinite valuesss"
  tmp[is.na(tmp)] <- "in progress"
  print(tmp)
}

#' @export
AIC.dispfit <- function(x, ...)
{
  print(x$distribution.selection[1:7])
}
