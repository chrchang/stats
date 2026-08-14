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

  x0 <- -2 * 10^-c(22,10,7,5)
  suppressWarnings(expect_warning(fx0 <- dbinom(x0, size = 3, prob = 0.1)))
  expect_all_equal(fx0, 0)

  expect_no_error(lapply(sample(10000, size=1000), function(M) {
    ## Range reduced for now due to 2^52 limit.
    ## n <- (M/100)*10^(2:20)
    n <- (M/100)*10^(2:11)
    if (anyNA(P <- dbinom(1,n,0.5))) {
      stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
    }
  }))
})
