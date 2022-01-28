// C++ code to be used within eigvaldisp:::AVar.VRR_pfc(),
// with the package RcppParallel.
// Parallelized with IntelTBB; this seems thread-safe.
// The environment variable GOMP_CPU_AFFINITY should be unset to enable
// multi-threading (this should be done before launching R)

#include <Rcpp.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#define _USE_MATH_DEFINES
#include <cmath>
using namespace Rcpp;
using namespace RcppParallel;

struct Compute_Cov : public Worker {
    const size_t p, ln, nrE;
    const RVector<double> R;
    const RVector<double> E2;
    const RVector<double> ni;
    std::vector<double> out;
    Compute_Cov(const size_t p_, const size_t ln_, const size_t nrE_,
        const NumericVector R_, const NumericVector E2_,
        const NumericVector ni_): p(p_), ln(ln_), nrE(nrE_),
        R(R_), E2(E2_), ni(ni_), out() { out.resize(ln, 0.0); }
    Compute_Cov(const Compute_Cov& Cov, Split):
            p(Cov.p), ln(Cov.ln), nrE(Cov.nrE),
            R(Cov.R), E2(Cov.E2), ni(Cov.ni), out() { out.resize(ln, 0.0); }
    void operator()(std::size_t begin, std::size_t end) {
        double rij, rkl, rik, ril, rjk, rjl, C;
        double Eij, Ekl, Cijkl;
        for(std::size_t i = begin; i < end; i++) { // Note i < end, not end - 1
            for(std::size_t j = i + 1; j < p; j++) {
                for(std::size_t k = i; k < p - 1; k++) {
                    for(std::size_t l = ((i == k) ? (j + 1) : (k + 1)); l < p; l++) {
                        rij = R[i + j * p];
                        rkl = R[k + l * p];
                        rik = R[i + k * p];
                        rjl = R[j + l * p];
                        ril = R[i + l * p];
                        rjk = R[j + k * p];
                        C = rij * rkl * (rik * rik + ril * ril +
                                         rjk * rjk + rjl * rjl) / 2 +
                            rik * rjl + ril * rjk - (rij * rik * ril
                            + rij * rjk * rjl + rik * rjk * rkl + ril * rjl * rkl);
                        for(std::size_t m = 0; m < ln; m++) {
                            Eij = E2[i + j * p - j * (2 * p - j + 1) / 2 + m * nrE];
                            Ekl = E2[k + l * p - l * (2 * p - l + 1) / 2 + m * nrE];
                            Cijkl = C * ni[m];
                            out[m] += (Eij * Ekl + Cijkl) * Cijkl;
                        }
                    }
                }
            }
        }
    }
    void join(const Compute_Cov& rhs) {
        for(size_t ii = 0; ii < ln; ii++){
            out[ii] += rhs.out[ii];
        }
    }
};

//*** Cov_r2P ***
//' Sum of covariances of squared correlations via \code{RcppParallel}.
//'
//' @rdname Cov_r2
//'
// [[Rcpp::export]]
NumericVector Cov_r2P(NumericVector n, NumericVector R,
                      NumericVector E, int dummy_nthreads) {
    const size_t ln  = static_cast<size_t>(n.size());
    const size_t p   = static_cast<size_t>(sqrt(R.size()));
    const size_t nrE = p * (p - 1) / 2;
    const NumericVector E2 = M_SQRT2 * E;
    const NumericVector ni = 1 / n;
    Compute_Cov Computing_Cov(p, ln, nrE, R, E2, ni);
    parallelReduce(0, p, Computing_Cov);
    for(size_t i = 0; i < ln; i++){
        Computing_Cov.out[i] *= 4.0;
    }
    return wrap(Computing_Cov.out);
}
