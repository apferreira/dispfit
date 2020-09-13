summary.dispfit <- function(x, ...) {
  stopifnot(inherits(x, "dispfit"))
  print("Distribution parameters")
  print(x$distribution.parameters)
  print("Distribution selection")
  print(x$distribution.selection)
}
