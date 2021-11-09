# Functions to check inputs for type errors and property constraints
#
# AMKAT package for R
# Copyright (C) 2021, Brian Neal
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# attempts to coerce input to numeric matrix
.convertToNumericMatrix <- function(x) {
  if (!is.numeric(x)) {
    return(as.matrix(sapply(x, as.numeric)))
  } else {
    return(as.matrix(x))
  }
}

# checks that the argument is nonempty
.checkNonEmpty <- function(arg_name, arg_value) {
  if (length(arg_value) == 0) stop(paste0("'", arg_name, "' has zero length"))
}

# checks that an argument is TRUE or FALSE
.checkTrueOrFalse <- function(arg_name, arg_value) {
  if (length(arg_value) != 1 | sum(arg_value %in% c(TRUE, FALSE)) != 1) {
    stop(paste0("'", arg_name, "' must evaluate to \"TRUE\" or \"FALSE\""))
  }
}

# checks that an argument is a strictly-positive integer
.checkPositiveInteger <- function(arg_name, arg_value) {
  if (length(arg_value) != 1) {
    stop(paste0("'", arg_name, "' must be a finite, strictly-positive integer"))
  }
  if (is.na(arg_value) | !is.numeric(arg_value)) {
    stop(paste0("'", arg_name, "' must be a finite, strictly-positive integer"))
  }
  if (!is.finite(arg_value)) {
    stop(paste0("'", arg_name, "' must be a finite, strictly-positive integer"))
  }
  if (arg_value %% 1 != 0 | arg_value < 1) {
    stop(paste0("'", arg_name, "' must be a finite, strictly-positive integer"))
  }
}

# checks that y and x have matching row dimension and no NA/Inf values while
# enforcing minimum feature dimension and sample size
# y and x are numeric matrices
.checkYX <- function(y, x) {

  min_sample_size <- 16

  if (nrow(y) < min_sample_size) {
    stop(paste0("'y' must have more than ", min_sample_size, " rows to ensure ",
                "numerical stability"))
  }
  if (nrow(x) != nrow(y)) stop("'y' and 'x' must have the same number of rows")
  if (sum(is.na(y)) != 0) stop("'y' contains NA/NaN values")
  if (sum(is.na(x)) != 0) stop("'x' contains NA/NaN values")
  if (sum(is.finite(y)) != length(y)) stop("'y' contains Inf/-Inf values")
  if (sum(is.finite(x)) != length(x)) stop("'x' contains Inf/-Inf values")
}

# checks that x meets minimum dimensions and has no missing/infinite values
.checkX <- function(x) {

  min_sample_size <- 16

  if (nrow(x) < min_sample_size) {
    stop(paste0("'x' must have more than ", min_sample_size, " rows to ensure ",
                "numerical stability"))
  }
  if (sum(is.na(x)) != 0) stop("'x' contains NA/NaN values")
  if (sum(is.finite(x)) != length(x)) stop("'x' contains Inf/-Inf values")
}

# checks that covariates are either null or have positive length
.checkCovariateArgument <- function(covariates) {
  if (!is.null(covariates) & length(covariates) == 0) {
    stop(paste0("'covariates' must either be NULL or have positive length"))
  }
}

# checks for character vector with valid kernel ID string entries
.checkCandidateKernels <- function(candidate_kernels) {
  .checkNonEmpty('candidate_kernels', candidate_kernels)
  if (!is.character(candidate_kernels) | !is.vector(candidate_kernels) |
      sum(candidate_kernels %in% listAmkatKernelFunctions()) <
      length(candidate_kernels)) {
    stop(paste0(
      "'candidate_kernels' must be a character vector containing ",
      "one or more of the following values: \"",
      paste(listAmkatKernelFunctions(), collapse = "\", \""), "\""))
  }
}

# checks for valid kernel ID string
.checkKernelFunction <- function(kernel_function) {
  if (!is.character(kernel_function) | length(kernel_function) != 1 |
      sum(kernel_function %in% listAmkatKernelFunctions()) != 1) {
    stop(paste0(
      "'kernel_function' must be one of the following values: \"",
      paste(listAmkatKernelFunctions(), collapse = "\", \""), "\""))
  }
}

# checks for character vector with valid entries for 'p_value_adjustment'
.checkPValueAdjustment <- function(p_value_adjustment) {
  if (length(p_value_adjustment) != 1 |
      sum(p_value_adjustment %in% c('pseudocount', 'floor', 'none')) != 1) {
    stop(paste0(
      "value of 'p_value_adjustment' must be ",
      "either \"pseudocount\", \"floor\" or \"none\""))
  }
}

# checks covariates for dimension and missing values
# covariates is a nonempty numeric matrix
# n is a positive integer
.checkCovariateContent <- function(covariates, n) {
  if (nrow(covariates) != n) {
    stop(paste0(
      "'covariates' must have the same row dimension as 'y' and 'x'"))
  }
  if (sum(is.na(covariates)) != 0) {
    stop("'covariates' contains NA/NaN values")
  }
  if (sum(is.finite(covariates)) != length(covariates)) {
    stop("'covariates' contains Inf/-Inf values")
  }
  if (ncol(covariates) > n - 2) {
    stop(paste0("cannot fit null model when ncol(covariates) > nrow(y) - 2"))
  }
}