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
  print(x$parameters)
}

#' @export
print.dispfit <- function(x, ...)
{
  print("Distribution selection")
  print(x$distribution.selection)
  print("Distribution parameters")
  print(x$distribution.parameters)
}

#' @export
AIC.dispfit <- function(x, ...)
{
  print(x$distribution.selection[1:7])
}
