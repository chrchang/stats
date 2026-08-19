test_that("dbinom works", {
  ## Based on dbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R .
  n0 <- 50; n1 <- 16; n2 <- 20; n3 <- 8
  for (n in stats::rbinom(n1, size = 2*n0, p = .4)) {
    for (p in c(0,1,rbeta(n2, 2, 4))) {
      for (k in stats::rbinom(n3, size = n, prob = runif(1))) {
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

  ## This test is quite slow, since exactr::dbinom uses Rmpfr for the larger
  ## exponents.  May want to skip most of the Rmpfr-requiring cases in some
  ## circumstances.
  expect_no_error(lapply(sample(10000, size=1000), function(M) {
    n <- (M/100)*10^(2:20)
    if (anyNA(P <- dbinom(1,n,0.5))) {
      stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
    }
  }))

  expect_all_equal(dbinom(2^c(0:1023, 1023.999), size=Inf, prob = .1), 0)

  x. <- 1.20e308
  N <- 1.72e308
  prb <- print(seq(13, 127, by = 6))/128
  db <- dbinom(x., N, prob = prb, log = TRUE)
  expect_equal(-2^1012 * c(3978.52477729, 3004.42235841, 2321.14764068, 1804.10617471, 1395.99909831,
                           1065.91552786, 795.509574314, 573.263664821, 391.782927199, 246.423578896,
                           134.598198071, 55.4909088367, 10.1000515546, 1.6625043907, 36.710476479,
                           127.476528317, 297.82558994, 600.93919587, 1195.32578926, 3368.52998705),
               db,
               tolerance=1e-11)

  x <- c(12:20, 100*c(1,3,10,20,50))
  db <- dbinom(x, size=x+1, prob = 2^-1024.1, log = TRUE)
  expect_equal(c(-8515.66, -9225.44, -9935.22, -10645, -11354.8, -12064.6, -12774.4,
                 -13484.2, -14194, -70980.6, -212950, -709845, -1419700, -3549250),
               db,
               tolerance=1e-5)
})

test_that("pbinom works", {
  ## Based on pbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R .
  n0 <- 50; n1 <- 16; n2 <- 20; n3 <- 8
  for (n in stats::rbinom(n1, size = 2*n0, p = .4)) {
    for (p in c(0,1,rbeta(n2, 2, 4))) {
      for (k in stats::rbinom(n3, size = n, prob = runif(1))) {
        if (k!=n && p!=0) {
          expect_equal(       pbinom(0:k, size=n, prob=p),
                       cumsum(dbinom(0:k, size=n, prob=p)))
        }
      }
    }
  }

  set.seed(123)
  n <- 20
  Rbinom <- sort(unique(stats::rbinom(n, size = 55, prob = pi/16)))
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
  Rbinom <- sort(unique(stats::rbinom(n, size = 55, prob = pi/16)))
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
  ## Based on dhyper() tests in R 4.6.1 tests/d-p-q-r-tst-2.R and
  ## tests/print-tests.R .
  expect_no_error(lapply(sample(10000, size=1000), function(M) {
    ## Domain reduced for now due to 2^52 limit.
    # n <- (M/100)*10^(2:20)
    n <- (M/100)*10^(2:11)
    if (anyNA(P <- dhyper(n+1,n+5,n+5,n))) {
      # stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
      stop("NA for M=", M, "; 10ex=",paste((2:11)[is.na(P)], collapse=", "))
    }
  }))

  N1 <- 10
  N2 <- 7
  n <- 8
  x <- 0:N1
  Mhyp2 <- dhyper(x, N1, N2, n)
  expect_equal(Mhyp2, c(0, 0.0004113534, 0.01295763, 0.103661, 0.3023447,
                        0.3628137, 0.1814068, 0.03455368, 0.00185109, 0, 0),
               tolerance=6e-7)
})

test_that("phyper works", {
  ## Based on phyper() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R and
  ## tests/print-tests.R .
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

  set.seed(123)
  n <- 20
  Rhyper <- sort(unique(stats::rhyper(n, m=40, n=30, k=20)))
  Phyper <- phyper(Rhyper, m=40, n=30, k=20)
  expect_equal(log1p(-Phyper), phyper(Rhyper, m=40, n=30, k=20, lower=F, log=T))

  expect_all_equal(phyper(c(0:3, 1e67), 0, 0, 0), 1)

  N1 <- 10
  N2 <- 7
  n <- 8
  x <- 0:N1
  Mhyp1 <- phyper(x, N1, N2, n)
  expect_equal(Mhyp1, c(0, 0.0004113534, 0.01336898, 0.117030, 0.4193747,
                        0.7821884, 0.9635952, 0.99814891, 1.00000000, 1, 1),
               tolerance=6e-7)
})

test_that("qhyper works", {
  ## Based on qhyper() tests in R 4.6.1 tests/d-p-q-r-tests.R .
  set.seed(123)
  n <- 20
  Rhyper <- sort(unique(stats::rhyper(n, m=40, n=30, k=20)))
  Phyper <- phyper(Rhyper, m=40, n=30, k=20)
  ep <- 1e-7
  f1 <- 1 - 1e-7
  expect_equal(Rhyper, qhyper(Phyper*f1, m=40, n=30, k=20))
  p1 <- 1 + ep
  expect_equal(Rhyper, qhyper(p1-Phyper, m=40, n=30, k=20, lower=F))
  expect_equal(Rhyper, qhyper(log(Phyper)-ep, m=40, n=30, k=20, log=TRUE))
  expect_equal(Rhyper, qhyper(log1p(-Phyper)+ep, m=40, n=30, k=20, lower=F, log=T))
})

test_that("fisher.test works", {
  ## Based on fisher.test() tests in R 4.6.1 tests/reg-tests-1{a,b,d,e}.R .
  x <- matrix(c(2, 2, 4, 8, 6, 0, 1, 1, 7, 8, 1, 3, 1, 3, 7, 4, 2, 2, 2,
                1, 1, 0, 0, 0, 0, 0, 1, 1, 2, 0, 1, 1, 0, 2, 1, 0, 0, 0),
              nc = 2)
  fisher.test(x)

  fisher.test(cbind(0, c(0,0,0,1)))

  ## R PR#4688: once gave p.value=Inf, now gives FEXACT error 501
  ## (After new algorithm is fully implemented, this case should no longer
  ## yield an error at all.  But we'll probably want to skip it under some
  ## circumstances, it looks like it'll be slow.)
  reli <- cbind(Si = c(2121, 100, 27, 0),
                av = c(4700, 216, 67, 0),
                Nc = c(6234,2461,502,14))
  expect_error(fisher.test(reli, workspace=2000000))

  a <- diag(1:3)
  p <- fisher.test(a, simulate.p.value=TRUE)$p.value
  expect_gt(p, 0.001)

  dd <- data.frame(group=1, score=c(rep(0,14), rep(1,29), rep(2, 16)))[rep(1:59, 2),]
  dd[,"group"] <- c(rep("DOG", 59), rep("kitty", 59))
  Pv <- with(dd, fisher.test(score, group)$p.value)
  expect_gte(Pv, 0)
  expect_lte(Pv, 1)

  set.seed(7)
  t44 <- table(sample(LETTERS[1:4], size=50, replace=TRUE),
               sample(letters[1:4], size=50, replace=TRUE))
  ft44 <- do.call(fisher.test, list(t44))
  expect_equal(t44, eval(str2lang(ft44$data.name)))

  ## (This case might also become tractable soon.)
  d <- matrix(c(1,0,5,2,1,90
               ,2,1,0,2,3,89
               ,0,0,0,1,0,14
               ,0,0,0,0,0, 5
               ,0,0,0,0,0, 2
               ,0,0,0,0,0, 2
                ), nrow=6, byrow=TRUE)
  expect_error(fisher.test(d))
})
