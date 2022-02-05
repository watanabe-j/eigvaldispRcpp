## Documentation of the package
#' eigvaldispRcpp: Rcpp Extension for 'eigvaldisp'
#'
#' This package extends the main package \code{eigvaldisp} with
#' fast C++ functions via the \code{Rcpp} interface.
#'
#' The DESCRIPTION file:
#' \packageDESCRIPTION{eigvaldispRcpp}
#' \packageIndices{eigvaldispRcpp}
#'
#' @section Author/Maintainer:
#' Junya Watanabe <jw2098@cam.ac.uk>
#'
#' @references
#' Watanabe, J. (2021). eigvaldisp: R package for statistics of
#'   eigenvalue dispersion. https://github.com/watanabe-j/eigvaldisp.
#'
#' Watanabe, J. (2022). Statistics of eigenvalue dispersion indices:
#'  quantifying the magnitude of phenotypic integration. *Evolution*,
#'  **76**, 4--28. doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
#'
#' @seealso
#' \code{\link[eigvaldisp:eigvaldisp-package]{eigvaldisp}}: The main package
#'
#' \code{\link[eigvaldisp]{Var.VRR}}, \code{\link[eigvaldisp]{AVar.VRR_pfc}}:
#'   The C++ functions are internally used in these outer functions
#'
#' @examples
#' ## Generate a (moderately) large population correlation matrix
#' Rho <- GenCov(p = 50L, shape = "linearly_decreasing", evectors = "Givens")
#'
#' ## Arbitrary sample sizes
#' Ns <- c(10, 20, 50)
#'
#' ## Calculate the variance of Vrel(Rho)
#' ## When the argument cppfun is provided, Var.VRR() automatically calls
#' ## eigvaldisp:::AVar.VRR_pfc():
#' eigvaldisp::Var.VRR(Rho, Ns - 1, cppfun = "Cov_r2C") # Base Rcpp
#' eigvaldisp::Var.VRR(Rho, Ns - 1, cppfun = "Cov_r2A") # RcppArmadillo + OpenMP
#' eigvaldisp::Var.VRR(Rho, Ns - 1, cppfun = "Cov_r2E") # RcppEigen + OpenMP
#' eigvaldisp::Var.VRR(Rho, Ns - 1, cppfun = "Cov_r2P") # RcppParallel
#' ## These use four alternative functions provided and yield identical results
#' ## The last three parallelize the computation when possible
#' ## Speed may vary, depending on environment
#'
#' @docType package
#' @name eigvaldispRcpp-package
#' @aliases eigvaldispRcpp
#'
#' @import eigvaldisp
#'
# ## Setting Rcpp-related import
#' @useDynLib eigvaldispRcpp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom RcppParallel RcppParallelLibs
NULL

## Documentation of C++ functions
##
#' Sum of covariances of squared correlations
#'
#' Calculate (twice) sum of covariances of
#' all non-redundant pairs of squared correlation coefficients.
#' Used to evaluate \eqn{Var[Vrel(R)]} with Pan & Frank's approximation.
#'
#' These are internal C++ functions to be used in within
#' \code{eigvaldisp:::AVar.VRR_pfc()} (which in turn is typically called by
#' the front-end \code{eigvaldisp::Var.VRR()}). The choice among these C++
#' functions is made by the argument \code{cppfun} in the outside R functions.
#'
#' All these functions implement essentially the same algorithm.  The default
#' is \code{Cov_r2C}, which is basic but should be fast enough for correlation
#' matrices with \eqn{p} up to 100 or so. The others aim to speed-up
#' the evaluation for larger matrices via parallelization.
#'
#' \code{Cov_r2A} and \code{Cov_r2E} are parallelized with OpenMP, whereas
#' \code{Cov_r2P} is with IntelTBB via \code{RcppParallel}
#' (when the environment allows). They tend to show very similar performance,
#' though substantial environment-specific variation has been noted.
#' \code{Cov_r2P} seems to show similar performance across various environments,
#' but tends to be slightly slower than the others when OpenMP is properly
#' supported. It also does not allow the user to specify the number of threads,
#' hence typically ends up using all processors available. On the other hand,
#' the number of threads can be specified in \code{Cov_r2A} and \code{Cov_r2E}
#' thanks to OpenMP. Unfortunately, however, OpenMP does not seem to have
#' a native support in mac environments ([R Installation and Administration
#' manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html)).
#'
#' For \code{Cov_r2A} and \code{Cov_r2E}, the number of parallel threads is
#' controlled by the argument \code{nthreads} in the front-end R functions.
#' The default value \code{0L} automatically sets itself into one-half of
#' the number of processors detected, i.e., \code{omp_get_num_procs()/2}.
#' Setting this value beyond the number of physical cores typically does not
#' lead to better performance (or can actually compromise it).
#'
#' Technically, the argument \code{E} is assumed to have the form of
#' \code{matrix(sapply(n, eigvaldisp:::Exv.r1, x = R[upper.tri(R)]),
#' ncol = length(n))}.
#' Internally, \code{n}, \code{R} and \code{E} are converted into
#' \code{Rcpp::NumericVector} in \code{Cov_r2C()} and \code{Cov_r2P()},
#' \code{arma::vec} in \code{Cov_r2A()},
#' and \code{Eigen::ArrayXd} in \code{Cov_r2E()}.
#'
#' @name Cov_r2
#' @rdname Cov_r2
#'
#' @param n
#'   Degrees of freedom (not sample sizes); numeric of length 1 or more.
#' @param R
#'   Population correlation matrix; assumed to be validly constructed;
#'   numeric of length \code{p * p}
#' @param E
#'   Matrix of expectation of correlation coefficients corresponding to
#'   the upper triangular of \code{R}; numeric of length
#'   \code{p * (p - 1) / 2 * length(n)}.
#' @param dummy_nthreads
#'   Dummy argument to be ignored in \code{Cov.r2C} and \code{Cov.r2P}.
#'   For coding compatibility with other functions that use \code{nthreads}.
#' @param nthreads
#'   Integer to specify the number of threads in \code{Cov.r2A}
#'   and \code{Cov.r2E}. The default is one-half of the number of processors
#'   detected. See Details.
#'
#' @return
#'   A numeric vector of \eqn{\sum 2 Cov(r_{ij}^2, r_{kl}^2)},
#'   corresponding to \code{n}.
#'
#' @seealso
#'   \code{\link[eigvaldisp]{Exv.VXX}}, \code{\link[eigvaldisp]{AVar.VRR_xx}}
#'
NULL
