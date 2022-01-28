// C++ code to be used within eigvaldisp:::AVAR.VRR_pfc(),
// with the package RcppArmadillo.
// Parallelized with OpenMP; this seems thread-safe, unlike Rcpp::NumericVector.

#define ARMA_USE_BLAS
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#define _USE_MATH_DEFINES
#include <cmath>
#ifdef _OPENMP
    #include <omp.h>
    // [[Rcpp::plugins(openmp)]]
    // #define ARMA_USE_OPENMP
#endif
using namespace arma;

//*** Cov_r2A ***
//' Sum of covariances of squared correlations via \code{RcppArmadillo}.
//'
//' @rdname Cov_r2
//'
//' @importFrom Rcpp evalCpp
//'
// [[Rcpp::export]]
arma::rowvec Cov_r2A(const arma::rowvec n, const arma::vec R,
                     const arma::vec E, int nthreads = 0) {
    const int ln = n.size();
    const int p = sqrt(R.size());
    const int nrE = p * (p - 1) / 2;
    const vec E2 = M_SQRT2 * E;
    const rowvec ni = 1 / n;
    rowvec out(ln);
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#pragma omp declare reduction(+: rowvec: omp_out += omp_in) \
        initializer(omp_priv = zeros<rowvec>(omp_orig.size()))
#pragma omp parallel reduction(+: out)
{
#endif
    double rij, rkl, rik, ril, rjk, rjl, C;
    double Eij, Ekl, Cijkl;
    for(int i = 0; i < p - 1; i++) {
#ifdef _OPENMP
#pragma omp for
#endif
        for(int j = i + 1; j < p; j++) {
            for(int k = i; k < p - 1; k++) {
                for(int l = ((i == k) ? (j + 1) : (k + 1)); l < p; l++) {
                    rij = R[i + j * p];
                    rkl = R[k + l * p];
                    rik = R[i + k * p];
                    rjl = R[j + l * p];
                    ril = R[i + l * p];
                    rjk = R[j + k * p];
                    C = rij * rkl * (rik * rik + ril * ril +
                                     rjk * rjk + rjl * rjl) / 2.0 +
                        rik * rjl + ril * rjk - (rij * rik * ril
                        + rij * rjk * rjl + rik * rjk * rkl + ril * rjl * rkl);
                    for(int m = 0; m < ln; m++) {
                        Eij = E2[i + j * p - j * (2 * p - j + 1) / 2 + m * nrE];
                        Ekl = E2[k + l * p - l * (2 * p - l + 1) / 2 + m * nrE];
                        Cijkl = C * ni[m];
                        out[m] += (Eij * Ekl + Cijkl) * Cijkl;
                    }
                }
            }
        }
    }
#ifdef _OPENMP
}
#endif
    out *= 4.0;
    return out;
}
