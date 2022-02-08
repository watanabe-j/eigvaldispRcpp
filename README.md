# eigvaldispRcpp
Rcpp extension of eigvaldisp

This package extends the R package
[`eigvaldisp`](https://github.com/watanabe-j/eigvaldisp) by
fast C++ functions using `Rcpp`-related packages.
This is to speed-up evaluation of the approximate variance of
the relative eigenvalue variance of correlation matrices with
`eigvaldisp::Var.VRR()`, which takes prohibitively long time in R, even
with vectorization and parallelization (see `eigvaldisp:::AVar.VRR_pfd()`).

This extension is provided as a separate package to minimize the dependency
of `eigvaldisp`. Practically, one would not need this extension unless
interested in applying `eigvaldisp::Var.VRR()` to large correlation
matrices (p > 100 or so).


## Installation
```
# install.packages("devtools")
devtools::install_github("watanabe-j/eigvaldispRcpp")
```

This package has the following dependencies:
```
Depends: eigvaldisp
Imports: Rcpp, RcppParallel
LinkingTo: Rcpp, RcppArmadillo, RcppEigen, RcppParallel
Suggests: testthat (>= 3.0.0)
SystemRequirements: GNU make
```

In case you are unfamiliar with `Rcpp`, it requires a working C++ compiler.
Most Linux and macOS users would already have one. Windows users will require
[Rtools](https://cran.r-project.org/bin/windows/Rtools/),
which comes with GNU `make`.
See http://www.rcpp.org/ or https://cran.r-project.org/package=Rcpp for details.


## See also
Watanabe, J. (2021). eigvaldisp: R package for statistics of eigenvalue dispersion. https://github.com/watanabe-j/eigvaldisp.

Watanabe, J. (2022). Statistics of eigenvalue dispersion indices: quantifying the magnitude of phenotypic integration. *Evolution*, **76**, 4&ndash;28. doi:[10.1111/evo.14382](https://doi.org/10.1111/evo.14382).
