# Generate Empirical Centralized Kernel Matrix
generateKernelMatrix <- function(x, kernel_function = "gau") {
  .checkNonEmpty("x", x)
  if (!is.matrix(x) | !is.numeric(x)) x <- .convertToNumericMatrix(x)
  .checkX(x)
  .checkKernelFunction(kernel_function)
  .Call(`_AMKAT_generateKernelMatrix`, x, kernel_function)
}