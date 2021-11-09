
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AMKAT

<!-- badges: start -->
<!-- badges: end -->

AMKAT (short for Adaptive Multivariate Kernel-based Association Test) is
a flexible, nonparametric test for joint association between a
multivariate quantitative response and a high-dimensional set of
features. It includes the following notable features:

-   A feature selection method for filtering noise variables out of the
    feature set
-   A kernel-based test statistic that uses a different kernel function
    for each response variable, with a method for selecting each kernel
    function individually
-   A method of adjusting for covariate effects on the response
    variables
-   All significant computations performed using compiled C++ routines

## Installation

You can install the development version of AMKAT using the following
commands in `R`:

    install.packages("devtools")
    devtools::install_github("brianpatrickneal/AMKAT")

## Example

Testing is straightforward using the `amkat` function. Below, we
generate toy data for a set of feature variables (`x`) and for a set of
response variables (`y`), as well as for a set of covariates (`w`), and
test for joint association between the features and response:

``` r
library(AMKAT)
set.seed(42)
n <- 20 # sample size
y <- matrix(rnorm(3 * n), nrow = n, ncol = 3) # 3 response variables
x <- matrix(rnorm(50 * n), nrow = n, ncol = 50) # 50 feature variables
w <- matrix(rnorm(2 * n), nrow = n, ncol = 2) # 2 covariates
test_results <- amkat(y, x, covariates = w)
str(test_results)
#> List of 15
#>  $ sample_size           : int 20
#>  $ y_dimension           : int 3
#>  $ x_dimension           : int 50
#>  $ number_of_covariates  : int 2
#>  $ null_residuals        : num [1:20, 1:3] 1.0686 -0.9554 -0.0843 0.6126 1.785 ...
#>  $ null_standard_errors  : num [1:3] 1.38 1.37 1.13
#>  $ filter_x              : logi TRUE
#>  $ selected_x_columns    : num [1:20, 1] 2 4 5 7 8 11 14 16 17 18 ...
#>  $ candidate_kernels     : chr [1:4] "lin" "quad" "gau" "exp"
#>  $ selected_kernels      : chr [1:3] "lin" "gau" "exp"
#>  $ test_statistic_value  : num 8.14
#>  $ number_of_permutations: num 1000
#>  $ permutation_statistics: num [1:1000, 1] 6.03 4.35 4.62 5.66 4.88 ...
#>  $ p_value_adjustment    : chr "Pseudocount value of (1 / 1000) added"
#>  $ p_value               : num 0.099
```

## Other Information

More details on the main function `amkat` can be found in its help file.
Type `?AMKATpackage` for an index of the other contents in the package’s
namespace.
