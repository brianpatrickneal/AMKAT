
# Apply phimr filter (permutation-based high-dimensional multiple response);
# return list of selected column indices for x
phimr <- function(y, x) {
  .checkNonEmpty("y", y)
  .checkNonEmpty("x", x)
  if (!is.matrix(y) | !is.numeric(y)) y <- .convertToNumericMatrix(y)
  if (!is.matrix(x) | !is.numeric(x)) x <- .convertToNumericMatrix(x)
  .checkYX(y, x)
  return(1 + .Call(`_AMKAT_applyAmkatFilter`, y, x)) # C++ index offset
}
