\name{amkat}
\alias{amkat}
\title{Adaptive Multivariate Kernel-based Association Test}
\description{
Tests nonparametrically for joint association between the columns of \code{x} and those of \code{y} while optionally adjusting for covariates, accommodating high-dimensional data for \code{x} without distributional assumptions using a kernel-based test statistic. Includes supervised methods for feature selection and kernel selection.
}
\usage{
amkat(y, x, covariates = NULL, filter_x = TRUE,
      candidate_kernels = c("lin", "quad", "gau", "exp"),
      num_permutations = 1000,
      p_value_adjustment = "pseudocount",
      num_test_statistics = 1,
      output_test_statistics = TRUE,
      output_selected_kernels = TRUE,
      output_selected_x_columns = TRUE,
      output_null_residuals = TRUE,
      output_p_value_only = FALSE)
}
\arguments{
  \item{y}{a numeric matrix containing data on the dependent variables, with  observations indexed by row.}

  \item{x}{a numeric matrix with the same number of rows as \code{y} containing data on the independent variables.}

  \item{covariates}{an optional numeric matrix with the same number of rows as \code{y} containing data on the covariates. The number of columns cannot exceed \code{nrow(y) - 2}.}

  \item{filter_x}{logical, indicating whether to apply AMKAT's permutation-based filter method for feature selection. Requires \code{x} to have at least two columns.}

  \item{candidate_kernels}{an optional character vector specifying the kernel functions to be considered during kernel selection. For a list of valid character strings and the kernel function corresponding to each, use \code{listAmkatKernelFunctions()}.}

  \item{num_permutations}{an optional strictly-positive integer specifying the number of test statistics from the permutation null distribution that are to be used in approximating the \emph{P}-value for the test.}

  \item{p_value_adjustment}{an optional character string specifying the method of adjustment to apply to the permutation-based \emph{P}-value for the test. Acceptable values are \code{"pseudocount"}, \code{"floor"} and \code{"none"}.}

  \item{num_test_statistics}{an optional strictly-positive integer indicating the number of observed test statistic values to generate. If greater than 1, the mean value of the statistics is used as the observed value for testing; this stabilizes variation in the test statistic introduced by AMKAT's permutation-based filter method for feature selection. Has no effect if \code{filter_x = FALSE}.}

  \item{output_test_statistics}{logical, indicating whether output should  include the values of all observed test statistics and permutation test statistics generated during testing. Has no effect if \code{output_p_value_only = TRUE}.}

  \item{output_selected_kernels}{logical, indicating whether output should include the kernel function selected for each column of \code{y}. Has no effect if \code{output_p_value_only = TRUE}.}

  \item{output_selected_x_columns}{logical, indicating whether output should include a record of which columns of \code{x} were selected for use in testing by AMKAT's filter method. Has no effect if \code{filter_x = FALSE} or if  \code{output_p_value_only = TRUE}.}

  \item{output_null_residuals}{logical, indicating whether output should include the residuals and standard errors from the fitted null model (after covariate adjustment, if applicable). Has no effect if \code{output_p_value_only = TRUE}.}

  \item{output_p_value_only}{logical; if \code{TRUE}, the function returns only the \emph{P}-value for the test rather than a list of results.}
}
\details{
A minimum requirement of 16 observations is enforced to avoid \code{NaN} values when estimating the asymptotic variance of the test statistic.

The kernel-based AMKAT test statistic incorporates individual kernel selection for each variable in \code{y} using a maximum statistic method. The kernel functions considered during kernel selection can be specified via the argument \code{candidate_kernels}. The function \code{listAmkatKernelFunctions()} retrieves the strings used to identify the available kernel functions, which include the Linear, Quadratic, Gaussian, Exponential and Identical-By-State (IBS) kernels. By default, all of these are considered except for the IBS kernel. \code{?listAmkatKernelFunctions} provides details on the individual kernel functions.

Prior to performing kernel selection and computing the test statistic, feature selection is performed using a filter method in which permuted data is compared to original data using tests of Spearman's Rho between columns of \code{x} and {y}.

An option is included to repeat the feature selection and kernel selection process in order to generate multiple test statistic values and use their mean value for testing, in order to reduce variation in the test statistic introduced by the feature selection method.

The \emph{P}-value for the test is computed by drawing a sample of test statistics from the permutation null distribution and comparing them to the value of the test statistic obtained from the original data. Permutation statistics are generated with feature selection and kernel selection reapplied to each permuted copy of the data. By default, the calculation of the \emph{P}-value includes a positive adjustment of \code{1/num_permutations} to avoid \emph{P}-values of 0, which are never possible for a permutation test using all possible permutations of the data (due to the identity permutation). Alternatively, \code{p_value_adjustment = "floor"} may be used to apply a floor of \code{1/num_permutations} to the \emph{P}-value in place of the adjustment, while \code{p_value_adjustment = "none"} will forego the adjustment and allow for \emph{P}-values of 0.

Covariate adjustment is performed prior to testing by using ordinary least squares to fit a null model in which the covariate effects are modeled as linear effects. The residuals and standard errors from this model are used in place of the raw values and estimated variances for \code{y} during testing.
}

\value{
  When invoked with \code{output_p_value_only = TRUE}, \code{amkat} returns a \code{double} containing the \emph{P}-value for the test; otherwise, it returns a \code{list} whose components vary depending on the arguments supplied to \code{amkat}. Components listed below are always present in the list unless otherwise specified:

  \item{sample_size}{the row dimension of \code{x} and of \code{y}.}

  \item{y_dimension}{the column dimension of \code{y}.}

  \item{x_dimension}{the column dimension of \code{x}.}

  \item{number_of_covariates}{the column dimension of \code{covariates}, or \code{0} if no covariates were used.}

  \item{null_residuals}{a numeric matrix of the same dimensions as \code{y} containing the residuals from the null model (including covariates, if applicable). Not included when \code{output_null_residuals = FALSE}.}

  \item{null_standard_errors}{a numeric vector of length \code{ncol(y)} containing the standard error for each column of \code{null_residuals}. Not included when \code{output_null_residuals = FALSE}.}

  \item{filter_x}{a logical value indicating whether AMKAT's permutation-based filter method for feature selection was applied during testing.}

  \item{selected_x_columns}{if \code{num_test_statistics = 1}, a numeric vector containing the indices of the columns of \code{x} selected by AMKAT's filter method. Otherwise, a numeric matrix with  \code{num_test_statistics} rows and \code{ncol(x)} columns, where the \emph{(i,j)}th entry is \code{1} if the \emph{j}th column of \code{x} was selected by AMKAT's filter when generating the \emph{i}th test statistic and \code{0} otherwise. Not included if \code{filter_x = FALSE} or \code{output_selected_x_columns = FALSE}.}

  \item{candidate_kernels}{a character vector indicating the candidate kernels that were used during AMKAT's kernel selection process.}

  \item{selected_kernels}{if \code{num_test_statistics = 1}, a character vector of length \code{ncol(y)} containing the kernel function selected for each \code{y} variable by AMKAT's kernel selection method. If \code{num_test_statistics > 1}, a character matrix with  \code{num_test_statistics} rows and \code{ncol(y)} columns, where the \emph{(i, j)}th entry is \code{1} if the \emph{j}th column of \code{x} was selected by AMKAT's filter when generating the \emph{i}th test statistic. Not included if \code{filter_x = FALSE} or \code{output_selected_x_columns = FALSE}.}

  \item{test_statistic_type}{a description of the test statistic used. Only included when using the mean of multiple observed test statistic values, i.e., when \code{filter_x = TRUE} and \code{num_test_statistics > 1}.}

  \item{generated_test_statistics}{a numeric vector containing the values of the observed test statistics. Only included when \code{output_test_statistics = TRUE} and \code{num_test_statistics > 1}.}

  \item{test_statistic_value}{the value of the observed test statistic, or the mean value of the observed test statistics when \code{num_test_statistics > 1}. Only included when \code{output_test_statistics = TRUE}.}

  \item{number_of_permutations}{the number of permutation test statistics used in testing.}

  \item{permutation_statistics}{a numeric vector containing the values of the permutation test statistics. Only included when \code{output_test_statistics = TRUE}.}

  \item{p_value_adjustment}{a character string describing the method of adjustment used for the \emph{P}-value.}

  \item{p_value}{the \emph{P}-value for the test.}
}

\references{Neal, Brian and He, Tao. \dQuote{An adaptive multivariate kernel-based test for association with multiple quantitative traits in high-dimensional data.} \emph{Genetic Epidemiology} (not yet submitted).}

\examples{
y <- matrix(rnorm(4 * 25), nrow = 25, ncol = 4)
x <- matrix(rnorm(200 * 25), nrow = 25, ncol = 200)
test_results <- amkat(y, x)

x <- matrix(rbinom(n = 200 * 25, size = 2, prob = 0.25),
            nrow = 25, ncol = 200)
w <- matrix(rnorm(2 * 25), nrow = 25, ncol = 2)
test_results <-
amkat(y, x, covariates = w,
      candidate_kernels = listAmkatKernelFunctions(),
      num_permutations = 500,
      num_test_statistics = 50,
      p_value_adjustment = "floor")
}
\author{Brian Neal}