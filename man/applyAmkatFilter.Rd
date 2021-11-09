\name{applyAmkatFilter}
\alias{applyAmkatFilter}
\title{Permutation-Based Filter Method For Feature Selection}
\description{
Selects columns of \code{x} related to one or more columns of \code{y} using a filter method in which permuted data is compared to original data using nonparametric tests of association between columns from \code{x} and columns from \code{y}.
}
\usage{
applyAmkatFilter(y, x)
}
\arguments{
  \item{y}{a numeric matrix containing data on the dependent variables, with at least two observations indexed by row.}

  \item{x}{a numeric matrix with the same number of rows as \code{y}, containing data on the independent variables.}
}

\details{
A two-sided test of the hypothesis that Spearman's Rho equals 0 is performed for each pair of columns, one from \code{x} and one from \code{y}. A row-permuted copy of \code{x} is created, and the tests are repeated using the permuted copy of \code{x}. A column of \code{x} passes the filter (i.e., is deemed to be related to \code{y}) if the minimum of the \emph{P}-values for that column from the original data is less than the minimum of the \emph{P}-values from the permuted data.

If no columns of \code{x} pass the filter, then the column with the lowest
minimum \emph{P}-value is selected by default.
}

\value{
An integer vector of variable length, containing the indices of the columns of \code{x} which passed the filter.
}

\references{Neal, Brian and He, Tao. \dQuote{An adaptive multivariate kernel-based test for association with multiple quantitative traits in high-dimensional data.} \emph{Genetic Epidemiology} (not yet submitted).}

\examples{
y <- matrix(rnorm(4 * 25), nrow = 25, ncol = 4)
x <- matrix(rnorm(50 * 25), nrow = 25, ncol = 50)
selected_x_columns <- applyAmkatFilter(y, x)
}

\author{Brian Neal}