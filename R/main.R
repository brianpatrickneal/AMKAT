# AMKAT: Adaptive Multivariate Kernel-Based Association Test
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


# amkat ------------------------------------------------------------------------
amkat <-
  function(y, x, covariates = NULL, filter_x = TRUE,
           candidate_kernels = c("lin", "quad", "gau", "exp"),
           num_permutations = 1000, p_value_adjustment = "pseudocount",
           num_test_statistics = 1, output_test_statistics = TRUE,
           output_selected_kernels = TRUE, output_selected_x_columns = TRUE,
           output_null_residuals = TRUE, output_p_value_only = FALSE) {

    .checkNonEmpty("y", y)
    .checkNonEmpty("x", x)
    if (!is.matrix(y) | !is.numeric(y)) y <- .convertToNumericMatrix(y)
    if (!is.matrix(x) | !is.numeric(x)) x <- .convertToNumericMatrix(x)
    .checkAmkatInputs(
      y, x, covariates, filter_x, candidate_kernels, num_permutations,
      p_value_adjustment, num_test_statistics, output_test_statistics,
      output_selected_kernels, output_selected_x_columns,
      output_null_residuals, output_p_value_only)

    null_fit <- .fitAmkatNullModel(y, x, covariates)

    if (ncol(x) == 1) filter_x <- FALSE
    if (output_p_value_only) {
      output <-
        .generateAmkatPvalue(null_fit, x, candidate_kernels, num_permutations,
                             filter_x, num_test_statistics, p_value_adjustment)
    } else {
      test_results <- .generateAmkatResults(
        null_fit, x, candidate_kernels, num_permutations, filter_x,
        num_test_statistics, output_selected_kernels, output_selected_x_columns,
        p_value_adjustment)
      output <- .formatAmkatOutput(
        nrow(y), ncol(y), ncol(x), null_fit, test_results,
        output_null_residuals, filter_x, output_selected_x_columns,
        candidate_kernels, output_selected_kernels, num_test_statistics,
        output_test_statistics, num_permutations)
    }
    return(output)
  }

# listAmkatKernelFunctions -----------------------------------------------------
# Acceptable values for kernel_function parameter of generateKernelMatrix()
# located in AMKAT/src/GenerateKernelMatrix.cpp

listAmkatKernelFunctions <- function() {
  return(c("Linear kernel" = "lin",
           "Quadratic kernel" = "quad",
           "Gaussian kernel" = "gau",
           "Exponential kernel" = "exp",
           "Identical-by-State kernel" = "IBS"))
}

# Internal helpers for amkat ---------------------------------------------------
# Helper function to check inputs
.checkAmkatInputs <-
  function(y, x, covariates, filter_x, candidate_kernels, num_permutations,
           p_value_adjustment, num_test_statistics, output_test_statistics,
           output_selected_kernels, output_selected_x_columns,
           output_null_residuals, output_p_value_only) {
    .checkYX(y, x)
    .checkCovariateArgument(covariates)
    .checkTrueOrFalse("filter_x", filter_x)
    .checkCandidateKernels(candidate_kernels)
    .checkPositiveInteger("num_permutations", num_permutations)
    .checkPValueAdjustment(p_value_adjustment)
    .checkPositiveInteger("num_test_statistics", num_test_statistics)
    .checkTrueOrFalse("output_test_statistics", output_test_statistics)
    .checkTrueOrFalse("output_selected_kernels", output_selected_kernels)
    .checkTrueOrFalse("output_selected_x_columns", output_selected_x_columns)
    .checkTrueOrFalse("output_null_residuals", output_null_residuals)
    .checkTrueOrFalse("output_p_value_only", output_p_value_only)
  }

# Helper function to fit null model
.fitAmkatNullModel <- function(y, x, covariates) {
  n <- nrow(y)
  if (is.null(covariates)) {
    covdim <- 0
    hat_matrix <- matrix(data = 1 / n, nrow = n, ncol = n)
  } else {
    if (!is.matrix(covariates) | !is.numeric(covariates)) {
      covariates <- .convertToNumericMatrix(covariates)
    }
    .checkCovariateContent(covariates, n)
    covdim <- ncol(covariates)
    w_aug <- cbind(rep(1, times = n), covariates)
    hat_matrix <- w_aug %*% tcrossprod(solve(crossprod(w_aug, w_aug)), w_aug)
  }
  null_residuals <- (diag(n) - hat_matrix) %*% y
  null_standard_errors <-
    diag(crossprod(y, null_residuals)) / (n - covdim - 1)
  return(list("num_covariates" = covdim,
              "residuals" = null_residuals,
              "standard_errors" = null_standard_errors))
}

# Helper function to generate P-value
.generateAmkatPvalue <-
  function(null_fit, x, candidate_kernels, num_permutations,
           filter_x, num_test_statistics, p_value_adjustment) {
    if (filter_x) {
      test_statistic <- mean(
        .Call(`_AMKAT_generateTestStatMultiple`,
              null_fit$residuals, null_fit$standard_errors, x,
              candidate_kernels, num_test_statistics))
      permutation_statistics <-
        .Call(`_AMKAT_generatePermStats`,
              null_fit$residuals, null_fit$standard_errors, x,
              candidate_kernels, num_permutations)
    } else {
      test_statistic <-
        .Call(`_AMKAT_generateTestStatNoFilter`,
              null_fit$residuals, null_fit$standard_errors, x,
              candidate_kernels)$test_statistic
      permutation_statistics <-
        .Call(`_AMKAT_generatePermStatsNoFilter`,
              null_fit$residuals, null_fit$standard_errors, x,
              candidate_kernels, num_permutations)
    }
    p_value <-
      mean(test_statistic <= permutation_statistics)
    if (p_value_adjustment == 'pseudocount') {
      p_value <- p_value + 1 / num_permutations
    } else if (p_value_adjustment == 'floor') {
      p_value <- max(p_value, 1 / num_permutations)
    }
    return(p_value)
  }

# Helper function to generate full test results
.generateAmkatResults <- function(
  null_fit, x, candidate_kernels, num_permutations, filter_x,
  num_test_statistics, output_selected_kernels, output_selected_x_columns,
  p_value_adjustment) {

  if (filter_x) {
    if (num_test_statistics == 1) {
      test_results <-
        .Call(`_AMKAT_generateTestStat`,
              null_fit$residuals, null_fit$standard_errors, x,
              candidate_kernels)
      test_results$using_mean_observed_stat <- FALSE
    } else {
      if (output_selected_kernels | output_selected_x_columns) {
        test_results <-
          .Call(`_AMKAT_generateTestStatsAllResults`,
                null_fit$residuals, null_fit$standard_errors, x,
                candidate_kernels, num_test_statistics)
      } else {
        test_results <- list(
          "test_statistics" =
            .Call(`_AMKAT_generateTestStatMultiple`,
                  null_fit$residuals, null_fit$standard_errors, x,
                  candidate_kernels, num_test_statistics))
      }
      test_results$test_statistic <-
        mean(test_results$test_statistics)
      test_results$using_mean_observed_stat <- TRUE
    }
    test_results$permutation_statistics <-
      .Call(`_AMKAT_generatePermStats`,
            null_fit$residuals, null_fit$standard_errors, x,
            candidate_kernels, num_permutations)
  } else {
    test_results <-
      .Call(`_AMKAT_generateTestStatNoFilter`,
            null_fit$residuals, null_fit$standard_errors, x,
            candidate_kernels)
    test_results$using_mean_observed_stat <- FALSE
    test_results$permutation_statistics <-
      .Call(`_AMKAT_generatePermStatsNoFilter`,
            null_fit$residuals, null_fit$standard_errors, x,
            candidate_kernels, num_permutations)
  }
  p_value <-
    mean(test_results$test_statistic <= test_results$permutation_statistics)
  if (p_value_adjustment == 'pseudocount') {
    test_results$p_value <- min(1, p_value + 1 / num_permutations)
    test_results$pv_adjust_desc <-
      paste0('Pseudocount value of (1 / ', num_permutations, ') added')
  } else if (p_value_adjustment == 'floor') {
    test_results$p_value <- max(p_value, 1 / num_permutations)
    test_results$pv_adjust_desc <-
      paste0('Floor value of (1 / ', num_permutations, ') applied')
  } else {
    test_results$p_value <- p_value
    test_results$pv_adjust_desc <- 'No adjustment'
  }
  return(test_results)
}

# Helper function to format list output
.formatAmkatOutput <- function(
  n, y_dim, p, null_fit, test_results, output_null_residuals, filter_x,
  output_selected_x_columns, candidate_kernels, output_selected_kernels,
  num_test_statistics, output_test_statistics, num_permutations) {

  out <- list(sample_size = n, y_dimension = y_dim,
              x_dimension = p, number_of_covariates = null_fit$num_covariates)
  if (output_null_residuals) {
    out$null_residuals <- null_fit$residuals
    out$null_standard_errors <- null_fit$standard_errors
  }
  out$filter_x <- filter_x
  if (filter_x & output_selected_x_columns) {
    out$selected_x_columns <- test_results$selected_x_columns
  }
  out$candidate_kernels <- candidate_kernels
  if (output_selected_kernels) {
    out$selected_kernels <- test_results$selected_kernels
  }
  if (test_results$using_mean_observed_stat) {
    out$test_statistic_type <-
      paste0('Mean of ', num_test_statistics, ' values, each generated by ',
             'repeating the feature selection and kernel selection process. ',
             'The variation in the test statistic introduced by feature ',
             'selection (for a particular set of data) decreases as more ',
             'values are used.')
    if (output_test_statistics) {
      out$generated_test_statistics <- test_results$test_statistics
    }
  }
  if (output_test_statistics) {
    out$test_statistic_value <- test_results$test_statistic
  }
  out$number_of_permutations <- num_permutations
  if (output_test_statistics) {
    out$permutation_statistics <- test_results$permutation_statistics
  }
  out$p_value_adjustment <- test_results$pv_adjust_desc
  out$p_value <- test_results$p_value
  return(out)
}