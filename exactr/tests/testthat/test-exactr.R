test_that("dbinom works", {
  ## Based on dbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R .
  n0 <- 50; n1 <- 16; n2 <- 20; n3 <- 8
  for (n in rbinom(n1, size = 2*n0, p = .4)) {
    for (p in c(0,1,rbeta(n2, 2, 4))) {
      for (k in rbinom(n3, size = n, prob = runif(1))) {
        if (k!=n && p!=0) {
          expect_equal(stats::pf((k+1)/(n-k)*(1-p)/p, df1=2*(n-k), df2=2*(k+1)),
                       sum(dbinom(0:k, size = n, prob = p)),
                       label=paste0("stats::pf((k+1)/(n-k)*(1-p)/p, df1=2*(n-k), df2=2*(k+1)) for n=", n, ", p=", p, ", k=", k))
        }
      }
    }
  }
  ## (more tests coming)
})
