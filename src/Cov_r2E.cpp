// C++ code to be used within eigvaldisp:::AVAR.VRR_pfc(),
// with the package RcppEigen.
// Parallelized with OpenMP; this seems thread-safe, unlike Rcpp::NumericVector.
// The right place for "#pragma omp for" (in terms of speed) seems to be
// above the for loop for j.
// Usually, omp_get_max_threads() gets the environment variable OMP_NUM_THREADS,
// which cannot be overwritten once an R session is initiated
// (Sys.setenv() does not seem to work for this).
// However, when omp_set_num_threads() is called, this can override
// the environment variable (and even num_threads() in the "#pragma omp" line)
// When using clang++, conflict with GOMP_CPU_AFFINITY seems to exist
// (correct result, but essentially single-threaded execution)

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#define _USE_MATH_DEFINES
#include <cmath>
#ifdef _OPENMP
    #include <omp.h>
    // [[Rcpp::plugins(openmp)]]
#endif
using Eigen::ArrayXd;

//*** Cov_r2E ***
//' Sum of covariances of squared correlations via \code{RcppEigen}.
//'
//' @rdname Cov_r2
//'
//' @importFrom Rcpp evalCpp
//'
// [[Rcpp::export]]
Eigen::ArrayXd Cov_r2E(const Eigen::ArrayXd n, const Eigen::ArrayXd R,
                       const Eigen::ArrayXd E, int nthreads = 0) {
    const int ln = n.size();
    const int p = sqrt(R.size());
    const int nrE = p * (p - 1) / 2;
    const ArrayXd E2 = M_SQRT2 * E;
    const ArrayXd ni = n.cwiseInverse();
    ArrayXd out = ArrayXd::Zero(ln);
#ifdef _OPENMP
    // Rcpp::Rcout << "Max no. of threads: " << omp_get_num_procs() << "\n";
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#pragma omp declare reduction(+: ArrayXd: omp_out += omp_in) \
        initializer(omp_priv = ArrayXd::Zero(omp_orig.size()))
    // Rcpp::Rcout << "No. of threads: " << nthreads << "\n";
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
