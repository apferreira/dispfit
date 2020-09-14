#' @export
summary.dispfit <- function(object, ...) {
  # stopifnot(inherits(x, "dispfit"))
  res <- list()
  res$parameters <- object$distribution.parameters
  res$selection <- object$distribution.selection
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
