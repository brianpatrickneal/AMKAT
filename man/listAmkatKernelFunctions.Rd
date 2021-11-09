\name{listAmkatKernelFunctions}
\alias{listAmkatKernelFunctions}
\title{List Available Kernel Functions In the AMKAT Package}
\description{
Returns a vector containing the character strings that identify each of the
    kernel functions included in the \code{AMKAT} package.
}

\usage{
listAmkatKernelFunctions()
}

\details{
The main function \code{amkat} includes an argument \code{candidate_kernels} that specifies the candidate kernel functions to consider during kernel selection. Its value must be a subset of the character vector returned by \code{listAmkatKernelFunctions}.

The function \code{generateKernelMatrix} in the \code{AMKAT} package includes an argument \code{kernel_function} that specifies the kernel function used to generate the empirical kernel matrix. The value of \code{kernel_function} must be a character string matching one of the entries in the value returned by \code{listAmkatKernelFunctions}.

Five kernel functions are available as of the current version: linear, quadratic, Gaussian, exponential, and Identical-by-State. Details on the individual kernel functions can be found below. For kernel functions that include a tuning parameter, the value of the parameter is set at \eqn{p}, the dimension of each of the kernel function's arguments (i.e., the column dimension of the argument \code{x} in the function \code{amkat}). Please note that the IBS kernel is only applicable to very specific types of data.
}

\value{
A character vector containing the following components, listed along with a description of the kernel function corresponding to each one:

  \item{"lin"}{Linear kernel, defined as \deqn{f(x_1, x_2) = x_1^T x_2 / p} for two \eqn{p}-dimensional column vectors \eqn{x_1} and \eqn{x_2}.}

  \item{"quad"}{Quadratic kernel, defined as \deqn{f(x_1, x_2) = (x_1^T x_2 / p + 1)^2} for two \eqn{p}-dimensional column vectors \eqn{x_1} and \eqn{x_2}.}

  \item{"gau"}{Gaussian kernel, defined as \deqn{f(x_1, x_2) = exp(-|| x_1 - x_2 ||^2 / p)} for two \eqn{p}-dimensional column vectors \eqn{x_1} and \eqn{x_2}, where \eqn{|| x_1 - x_2 ||} is the Euclidean distance.}

  \item{"exp"}{Exponential kernel, defined as \deqn{f(x_1, x_2) = exp(-(|| x_1 ||^2 + 3|| x_1 - x_2 ||^2 + || x_2 ||^2) / p)} for two \eqn{p} -dimensional column vectors \eqn{x_1} and \eqn{x_2}, where \eqn{|| . ||} is the Euclidean norm.}

  \item{"IBS"}{Identical-by-State kernel, defined as \deqn{(2p)^(-1) \sum{k=1 to k=p} (2 - |x_{1,k} - x_{2,k}|)} for two \eqn{p}-dimensional column vectors \eqn{x_1} and \eqn{x_2} written as \eqn{x_1 = (x_{1,1}, x_{1,2}, \ldots, x_{1,p})^T} and \eqn{x_2 = (x_{2,1}, x_{2,2}, \ldots, x_{2,p})^T}. Only applicable when the entries of \eqn{x_1} and \eqn{x_2} take on values in the set \eqn{{0, 1, 2}}, e.g., for additively-encoded single-nucleotide polymorphism (SNP) genotype data.}
}

\references{Neal, Brian and He, Tao. \dQuote{An adaptive multivariate kernel-based test for association with multiple quantitative traits in high-dimensional data.} \emph{Genetic Epidemiology} (not yet submitted).}


\author{Brian Neal}