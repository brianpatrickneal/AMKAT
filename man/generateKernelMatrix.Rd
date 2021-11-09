\name{generateKernelMatrix}
\alias{generateKernelMatrix}
\title{Centralized Kernel Similarity Matrix}
\description{
Given a \eqn{n \times p} matrix of sample observations indexed by row and a choice of kernel function, returns an empirical estimate of the \eqn{n \times n} matrix whose \eqn{(i,j)}th entry is the centralized version of the kernel function evaluated at the pair of rows \eqn{i} and \eqn{j} of \code{x}.
}

\usage{
generateKernelMatrix(x, kernel_function = "gau")
}

\arguments{
  \item{x}{a numeric matrix containing the sample data with observations indexed by row.}
  \item{kernel_function}{a character string identifying the kernel function to use. For a list of valid strings, use \code{listAmkatKernelFunctions()}. }
}

\details{
If the kernel function is expressed as \eqn{f(x_1, x_2)} for two \eqn{p}-dimensional vectors \eqn{x_1} and \eqn{x_2}, then \code{generateKernelMatrix} first computes the non-centralized empirical kernel matrix \eqn{K} whose \eqn{(i,j)}th entry is \eqn{f(x_i, x_j)}, where \eqn{x_i} and \eqn{x_i} are the \eqn{i}th and \eqn{j}th rows of \code{x}, respectively. Letting \eqn{K*} denote a copy of \eqn{K} with zero entries along the main diagonal and \eqn{J} denote a \eqn{n \times n} matrix of ones, the empirical centralized kernel matrix is then computed as
\deqn{K - (JK* + K*J - JK*J / n) / (n - 1)}.

Details on the available kernel functions can be found in the help topic for \code{listAmkatKernelFunctions}.
}

\value{
A square numeric matrix with the same row dimension as \code{x}, whose \eqn{(i,j)}th entry is the empirically centralized value for kernel-based pairwise similarity between rows \eqn{i} and \eqn{j} of \code{x}.
}

\references{Neal, Brian and He, Tao. \dQuote{An adaptive multivariate kernel-based test for association with multiple quantitative traits in high-dimensional data.} \emph{Genetic Epidemiology} (not yet submitted).}

\examples{
x <- matrix(rnorm(50 * 25), nrow = 25, ncol = 50)
ker_mat_lin <- generateKernelMatrix(x, "lin")
ker_mat_quad <- generateKernelMatrix(x, "quad")
ker_mat_gau <- generateKernelMatrix(x, "gau")
ker_mat_exp <- generateKernelMatrix(x, "exp")

x <- matrix(rbinom(n = 50 * 25, size = 2, prob = 0.25),
            nrow = 25, ncol = 50)
ker_mat_ibs <- generateKernelMatrix(x, "IBS")
}

\author{Brian Neal}