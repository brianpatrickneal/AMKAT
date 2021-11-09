# Supplemental functions for AMKAT


# Apply AMKAT's filter method; return list of selected column indices for x
applyAmkatFilter <- function(y, x) {
  .checkNonEmpty("y", y)
  .checkNonEmpty("x", x)
  if (!is.matrix(y) | !is.numeric(y)) y <- .convertToNumericMatrix(y)
  if (!is.matrix(x) | !is.numeric(x)) x <- .convertToNumericMatrix(x)
  .checkYX(y, x)
  return(1 + .Call(`_AMKAT_applyAmkatFilter`, y, x)) # C++ index offset
}

# Generate Empirical Centralized Kernel Matrix
generateKernelMatrix <- function(x, kernel_function = "gau") {
  .checkNonEmpty("x", x)
  if (!is.matrix(x) | !is.numeric(x)) x <- .convertToNumericMatrix(x)
  .checkX(x)
  .checkKernelFunction(kernel_function)
  .Call(`_AMKAT_generateKernelMatrix`, x, kernel_function)
}

# Internal functions -----------------------------------------------------------
# These functions are included for performing AMKAT faster within loops
# by foregoing argument checks and conversion; since they can cause R to crash
# if mishandled, they are left as internal functions and not exported. They
# can still be accessed using the AMKAT::: qualification prefix.
.applyAmkatFilter <- function(y, x) {
  return(1 + .Call(`_AMKAT_applyAmkatFilter`, y, x)) # C++ index offset
}
.generateKernelMatrix <- function(x, kernel_function) {
  .Call(`_AMKAT_generateKernelMatrix`, x, kernel_function)
}
.generatePermStats <- function(y, y_variances, x, candidate_kernels,
                               num_permutations) {
  .Call(`_AMKAT_generatePermStats`, y, y_variances, x,
        candidate_kernels, num_permutations)
}

.generatePermStatsNoFilter <- function(y, y_variances, x, candidate_kernels,
                                       num_permutations) {
  .Call(`_AMKAT_generatePermStatsNoFilter`, y, y_variances, x,
        candidate_kernels, num_permutations)
}

.generateTestStat <- function(y, y_variances, x, candidate_kernels) {
  .Call(`_AMKAT_generateTestStat`, y, y_variances, x, candidate_kernels)
}

.generateTestStatMultiple <- function(y, y_variances, x, candidate_kernels,
                                      num_test_statistics) {
  .Call(`_AMKAT_generateTestStatMultiple`, y, y_variances, x,
        candidate_kernels, num_test_statistics)
}

.generateTestStatNoFilter <- function(y, y_variances, x, candidate_kernels) {
  .Call(`_AMKAT_generateTestStatNoFilter`, y, y_variances, x, candidate_kernels)
}

.generateTestStatsAllResults <- function(y, y_variances, x, candidate_kernels,
                                         num_test_statistics) {
  .Call(`_AMKAT_generateTestStatsAllResults`, y, y_variances, x,
        candidate_kernels, num_test_statistics)
}