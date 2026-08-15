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
    ## Domain reduced for now due to 2^52 limit.
    # n <- (M/100)*10^(2:20)
    n <- (M/100)*10^(2:11)
    if (anyNA(P <- dbinom(1,n,0.5))) {
      # stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
      stop("NA for M=", M, "; 10ex=",paste((2:11)[is.na(P)], collapse=", "))
    }
  }))

  expect_all_equal(dbinom(2^c(0:1023, 1023.999), size=Inf, prob = .1), 0)

  ## Uncomment if/when x domain expanded.
  # x. <- 1.20e308
  # N <- 1.72e308
  # prb <- print(seq(13, 127, by = 6))/128
  # db <- dbinom(x., N, prob = prb, log = TRUE)
  # expect_equal(-2^1012 * c(3978.52477729, 3004.42235841, 2321.14764068, 1804.10617471, 1395.99909831,
  #                          1065.91552786, 795.509574314, 573.263664821, 391.782927199, 246.423578896,
  #                          134.598198071, 55.4909088367, 10.1000515546, 1.6625043907, 36.710476479,
  #                          127.476528317, 297.82558994, 600.93919587, 1195.32578926, 3368.52998705),
  #              db,
  #              tolerance=1e-11)

  ## Uncomment if/when denormal prob handled correctly.
  # x <- c(12:20, 100*c(1,3,10,20,50))
  # db <- dbinom(x, size=x+1, prob = 2^-1024.1, log = TRUE)
  # expect_equal(c(-8515.66, -9225.44, -9935.22, -10645, -11354.8, -12064.6, -12774.4,
  #              -13484.2, -14194, -70980.6, -212950, -709845, -1419700, -3549250),
  #            db,
  #            tolerance=1e-5)
})

test_that("pbinom works", {
  ## Based on pbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R .
  n0 <- 50; n1 <- 16; n2 <- 20; n3 <- 8
  for (n in rbinom(n1, size = 2*n0, p = .4)) {
    for (p in c(0,1,rbeta(n2, 2, 4))) {
      for (k in rbinom(n3, size = n, prob = runif(1))) {
        if (k!=n && p!=0) {
          expect_equal(       pbinom(0:k, size=n, prob=p),
                       cumsum(dbinom(0:k, size=n, prob=p)))
        }
      }
    }
  }

  set.seed(123)
  n <- 20
  Rbinom <- sort(unique(rbinom(n, size = 55, prob = pi/16)))
  Pbinom <- pbinom(Rbinom, size = 55, prob = pi/16)
  expect_equal(log1p(-Pbinom), pbinom(Rbinom, size=55, prob=pi/16, lower=F, log=T))

  pr <- 1e-23
  expect_equal(pr^12, pbinom(11, 12, prob=pr, lower=FALSE), tolerance=1e-12)

  x0 <- -2 * 10^-c(22,10,7,5)
  expect_all_equal(pbinom(x0, size=3, prob=0.1), 0)

  expect_no_error(lapply(sample(10000, size=1000), function(M) {
    ## Domain reduced for now due to 2^52 limit.
    # n <- (M/100)*10^(2:20)
    n <- (M/100)*10^(2:11)
    if (anyNA(P <- pbinom(1,n,0.5))) {
      # stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
      stop("NA for M=", M, "; 10ex=",paste((2:11)[is.na(P)], collapse=", "))
    }
  }))
})

test_that("qbinom works", {
  ## Based on qbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R .
  set.seed(123)
  n <- 20
  Rbinom <- sort(unique(rbinom(n, size = 55, prob = pi/16)))
  Pbinom <- pbinom(Rbinom, size = 55, prob = pi/16)
  ep <- 1e-7
  f1 <- 1 - 1e-7
  expect_equal(Rbinom, qbinom(Pbinom*f1, size=55, prob=pi/16))
  p1 <- 1 + ep
  expect_equal(Rbinom, qbinom(p1 - Pbinom, size=55, prob=pi/16, lower=F))
  expect_equal(Rbinom, qbinom(log(Pbinom)-ep, size=55, prob=pi/16, log=TRUE))
  expect_equal(Rbinom, qbinom(log1p(-Pbinom)+ep, size=55, prob=pi/16, lower=F, log=T))

  M <- .Machine$integer.max
  k <- 0:32
  for (n in c((M+1)/2, M, 2*M, 10*M)) {
    for (pr in c(1e-8, 1e-9, 1e-10)) {
      nDup <- !duplicated(pb <- pbinom(k, n, pr))
      qb <- qbinom(pb[nDup], n, pr)
      pn1 <- pb[nDup] < 1  # todo: check if MKL-on-RHEL guard is still needed
      expect_equal(k[nDup][pn1], qb[pn1])
    }
  }

  p.s <- c(.01, .001, .1, .25)
  pr <- (2:20)*1e-7
  sizes <- 1000*(5000 + c(0,6,16)) + 279
  k.s <- 0:15
  q.xct <- rep(k.s, each=length(pr))
  for (sz in sizes) {
    for (p in p.s) {
      qb <- qbinom(p=p, size=sz, prob=pr)
      pb <- stats::qpois(p=p, lambda=sz*pr)
      expect_equal(qb, pb)
    }
    pp.x <- outer(pr, k.s, function(pr, q) pbinom(q, size=sz, prob=pr))
    qq.x <- apply(pp.x, 2, function(p)     qbinom(p, size=sz, prob=pr))
    expect_equal(as.vector(qq.x), q.xct)
  }

  sz <- 6040:6045
  prb <- 0.995
  qb6 <- qbinom(p=0.05, size=sz, prob=prb)
  pqb6 <- pbinom(qb6, size=sz, prob=prb)
  pqb6_1 <- pbinom(qb6-1, size=sz, prob=prb)
  expect_equal(qb6, c(6001:6004, 6004:6005))
  expect_true(all((1 > pqb6) & (pqb6 >= 0.05)))
  expect_true(all((0.05 > pqb6_1) & (pqb6_1 >= 0.035)))
})

test_that("binom.test works", {
  ## From R 4.6.1 tests/reg-test-2.{R,Rout.save} .
  expect_snapshot({
    binom.test(c(800,10))
  })
})

test_that("dhyper works", {
  ## Based on dhyper() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R and
  ## tests/print-tests.{R,Rout.save} .
  expect_no_error(lapply(sample(10000, size=1000), function(M) {
    ## Domain reduced for now due to 2^52 limit.
    # n <- (M/100)*10^(2:20)
    n <- (M/100)*10^(2:11)
    if (anyNA(P <- dhyper(n+1,n+5,n+5,n))) {
      # stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
      stop("NA for M=", M, "; 10ex=",paste((2:11)[is.na(P)], collapse=", "))
    }
  }))
  # todo: more tests
})

test_that("phyper works", {
  ## Based on phyper() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R and
  ## tests/print-tests.{R,Rout.save} .
  .suppHyper <- function(m,n,k) max(0, k-n) : min(k, m)
  hyp.mn <- rbind(m = c(10, 15, 999),
                  n = c( 7,  0,   0))
  for (j in 1:ncol(hyp.mn)) {
    mn <- hyp.mn[,j]
    m <- mn[["m"]]
    n <- mn[["n"]]
    for (k in 2:m) {
      x <- .suppHyper(m,n,k)
      x <- c(x[1]-1L, x)
      expect_equal(phyper(x, m, n, k), cumsum(dhyper(x, m, n, k)))
      expect_equal(phyper(x, m, n, k, log.p=TRUE),
        log(cumsum(dhyper(x, m, n, k))))
    }
  }
  # todo: more tests
})
