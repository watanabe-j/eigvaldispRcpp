test_that("behavior of AVar.VRR_pfx under arbitrary conditions", {
    p <- sample(3:20, 1)
    N <- sample(3:20, 1)
    n <- N - 1
    Lambda <- sort(stats::rchisq(p, 5), decreasing = TRUE)
    Lambda <- Lambda / sum(Lambda) * p
    A <- GenCov(evalues = Lambda, evectors = "Givens")

    pfd  <- eigvaldisp:::AVar.VRR_pfd(A, n)
    pfcC <- eigvaldisp:::AVar.VRR_pfc(A, n, cppfun = "Cov_r2C")
    pfcA <- eigvaldisp:::AVar.VRR_pfc(A, n, cppfun = "Cov_r2A")
    pfcE <- eigvaldisp:::AVar.VRR_pfc(A, n, cppfun = "Cov_r2E")
    pfcP <- eigvaldisp:::AVar.VRR_pfc(A, n, cppfun = "Cov_r2P")
    expect_equal(pfd, pfcC)
    expect_equal(pfd, pfcA)
    expect_equal(pfd, pfcE)
    expect_equal(pfd, pfcP)
})
