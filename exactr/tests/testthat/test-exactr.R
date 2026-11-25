test_that("dbinom works", {
  ## Based on dbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R and
  ## tests/reg-tests-1{a,e}.R .
  set.seed(123)
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
  expect_equal(-2^1012 * c(3978.52477729, 3004.42235841, 2321.14764068,
                           1804.10617471, 1395.99909831, 1065.91552786,
                           795.509574314, 573.263664821, 391.782927199,
                           246.423578896, 134.598198071, 55.4909088367,
                           10.1000515546, 1.6625043907, 36.710476479,
                           127.476528317, 297.82558994, 600.93919587,
                           1195.32578926, 3368.52998705),
               db,
               tolerance=1e-11)

  dbArg <- cbind(x = c(1.465e+308, 1.4715e+308, 1.4762e+308, 1.4822e+308,
                       1.4869e+308, 1.4949e+308, 1.5034e+308, 1.5137e+308,
                       1.523e+308, 1.5305e+308, 1.5416e+308, 1.5486e+308,
                       1.5653e+308, 1.5853e+308, 1.639e+308),
                 size=c(1.6574e+308, 1.6514e+308, 1.7035e+308, 1.679e+308,
                        1.6531e+308, 1.6285e+308, 1.6993e+308, 1.6661e+308,
                        1.6801e+308, 1.6873e+308, 1.7506e+308, 1.7052e+308,
                        1.6752e+308, 1.6885e+308, 1.731e+308),
                 prob=c(27, 27, 26, 27, 28, 29, 28, 29, 29, 30, 29, 30, 32, 33,
                        35)/128)
  db <- apply(dbArg, 1, \(v3) do.call(dbinom, c(v3, list(log=TRUE))))
  expect_equal(db,
               c(-1.73031e+308, -1.76399e+308, -1.73534e+308, -1.74653e+308,
                 -1.7615e+308, -1.79181e+308, -1.7259e+308,  -1.77688e+308,
                 -1.77981e+308, -1.74055e+308, -1.70236e+308, -1.76548e+308,
                 -1.79598e+308, -1.79126e+308, -1.79514e+308),
                 tolerance=1e-5)

  x <- c(12:20, 100*c(1,3,10,20,50))
  db <- dbinom(x, size=x+1, prob = 2^-1024.1, log = TRUE)
  expect_equal(c(-8515.66, -9225.44, -9935.22, -10645, -11354.8, -12064.6,
                 -12774.4, -13484.2, -14194, -70980.6, -212950, -709845,
                 -1419700, -3549250),
               db,
               tolerance=1e-5)

  Pr <- stats::dmultinom(c(0,0,3), prob=c(1,1,14))
  expect_equal(Pr, dbinom(3, 3, p=14/16))

  ## get value of `expr` and keep warning as attribute (if there is one)
  getVaW <- function(expr, obj=FALSE) {
    W <- NULL
    withCallingHandlers(val <- expr,
                        warning = function(w) {
                          W <<- if(obj) w else conditionMessage(w)
                          invokeRestart("muffleWarning")
                        })
    structure(val %||% quote(._NULL_()), warning = W)
  }
  onWindows <- .Platform$OS.type == "windows"
  englishMsgs <- {
    if (nzchar(lang <- Sys.getenv("LANGUAGE"))) {
      lang == "en"
    } else {
      if (!onWindows) {
        lc.msgs <- sub("\\..*", "", print(Sys.getlocale("LC_MESSAGES")))
        lc.msgs == "C" || substr(lc.msgs, 1, 2) == en
      } else {
        lc.type <- sub("\\..*", "", sub("_.*", "", print(Sys.getlocale("LC_CTYPE"))))
        lc.type == "English" || lc.type == "C"
      }
    }
  }
  gd <- getVaW(dbinom(1234560:1234570, 9876543.2, .5))
  if (englishMsgs) {
    expect_identical(gd, structure(rep(NaN, 11), warning="NaNs produced"))
  }
})

test_that("pbinom works", {
  ## Based on pbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R and
  ## tests/reg-tests-1{a,b}.R .
  set.seed(123)
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

  n <- 20
  Rbinom <- sort(unique(stats::rbinom(n, size = 55, prob = pi/16)))
  Pbinom <- pbinom(Rbinom, size = 55, prob = pi/16)
  expect_equal(log1p(-Pbinom), pbinom(Rbinom, size=55, prob=pi/16, lower=F, log=T))

  pr <- 1e-23
  expect_equal(pr^12, pbinom(11, 12, prob=pr, lower=FALSE), tolerance=1e-12)

  x0 <- -2 * 10^-c(22,10,7,5)
  expect_all_equal(pbinom(x0, size=3, prob=0.1), 0)

  expect_no_error(lapply(sample(10000, size=1000), function(M) {
    n <- (M/100)*10^(2:20)
    if (anyNA(P <- pbinom(1,n,0.5))) {
      stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
    }
  }))

  x <- c(-1,0,1,2)
  expect_identical(pbinom(x, size=0, p=0.5), c(0,1,1,1))
  expect_identical(pbinom(x, size=0, p=0.5, lower.tail=FALSE), c(1,0,0,0))
  expect_identical(pbinom(c(-1, 20, 25), 20, 0.3, lower.tail=FALSE), c(1, 0, 0))
  expect_identical(pbinom(c(-1, 20, 25), 20, 0.3, lower.tail=FALSE, log.p=TRUE), c(0, -Inf, -Inf))

  expect_false(is.nan(pbinom(10, 1e6, 0.01, log.p=TRUE)))
})

test_that("qbinom works", {
  ## Based on qbinom() tests in R 4.6.1 tests/d-p-q-r-{tests,tst-2}.R and
  ## tests/reg-tests-1a.R .
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

  expect_equal(qbinom(0.95, 10, 1), 10)
  expect_equal(qbinom(0, 10, 1), 0)
  expect_equal(qbinom(0.95, 10, 0), 0)
  expect_equal(qbinom(0, 10, 0), 0)
  expect_equal(qbinom(0.95, 0, 0.5), 0)

  z <- qbinom(-Inf, 1, 0.5, log.p=TRUE)
  expect_true(is.finite(z))
})

test_that("binom.test works", {
  ## From R 4.6.1 tests/reg-test-2.{R,Rout.save} .
  expect_snapshot({
    binom.test(c(800,10))
  })

  ## One-sided p-values, including midp, and the x=0 and x=n boundaries.
  ## (alternative="greater", midp=TRUE previously returned
  ## P[X >= x-1] - 0.5 * P[X = x-1] instead of P[X > x] + 0.5 * P[X = x].)
  n <- 20
  p <- 0.3
  d <- stats::dbinom(0:n, n, p)
  for (x in 0:n) {
    less_want <- sum(d[1:(x+1)])
    greater_want <- sum(d[(x+1):(n+1)])
    half_px <- 0.5 * d[x+1]
    expect_equal(binom.test(x, n, p, alternative="less")$p.value, less_want)
    expect_equal(binom.test(x, n, p, alternative="greater")$p.value, greater_want)
    expect_equal(binom.test(x, n, p, alternative="less", midp=TRUE)$p.value,
                 less_want - half_px)
    expect_equal(binom.test(x, n, p, alternative="greater", midp=TRUE)$p.value,
                 greater_want - half_px)
    expect_equal(binom.test(x, n, p, alternative="greater", midp=TRUE, log.p=TRUE)$log.p.value,
                 log(greater_want - half_px))
  }
  expect_identical(binom.test(0, n, p, alternative="greater")$p.value, 1)
  expect_identical(binom.test(n, n, p, alternative="less")$p.value, 1)
  expect_equal(binom.test(7, 20, 0.3, alternative="greater", midp=TRUE)$p.value,
               stats::pbinom(7, 20, 0.3, lower.tail=FALSE) + 0.5 * stats::dbinom(7, 20, 0.3))

  ## method label
  expect_identical(binom.test(7, 20, 0.3)$method, "Exact binomial test")
  expect_identical(binom.test(7, 20, 0.3, log.p=TRUE)$method, "Exact binomial test")
  expect_identical(binom.test(7, 20, 0.3, midp=TRUE)$method, "Exact binomial test (midp)")
  expect_identical(binom.test(7, 20, 0.3, midp=TRUE, log.p=TRUE)$method, "Exact binomial test (midp)")
})

test_that("ppois works", {
  q <- c(-1, 0, 1, 3, 10, 50, Inf)
  for (lambda in c(0, 1e-5, 0.5, 4, 30, 1234.5)) {
    expect_equal(ppois(q, lambda), stats::ppois(q, lambda), tolerance=1e-13)
    expect_equal(ppois(q, lambda, lower.tail=FALSE),
                 stats::ppois(q, lambda, lower.tail=FALSE), tolerance=1e-13)
    expect_equal(ppois(q, lambda, log.p=TRUE),
                 stats::ppois(q, lambda, log.p=TRUE), tolerance=1e-13)
  }
  expect_equal(cumsum(stats::dpois(0:100, 37.5)), ppois(0:100, 37.5))
  expect_true(is.nan(suppressWarnings(ppois(1, -1))))
  ## Deep tail where non-log result underflows.
  expect_equal(ppois(10, 2000, log.p=TRUE),
               stats::ppois(10, 2000, log.p=TRUE), tolerance=1e-13)

  ## approx=FALSE (the default) promises <1 ULP error relative to the exact
  ## (infinite-precision) value, so usually only the two doubles bracketing it
  ## are acceptable.  (As of this writing, all of these cases are correctly
  ## rounded, but we don't want to fail if a future speed optimization gives up
  ## unpromised accuracy.)  Several of these cases are within 0.001 ULP of a
  ## rounding boundary.
  ##
  ## Each exact value is represented as hi + lo, where hi is the 320-bit MPFR
  ## result correctly rounded to double, and lo is the remainder rounded to
  ## double.  (In the subnormal case, lo is only accurate to 2^{-1074}, which
  ## is still enough to identify the two acceptable results.)  Inputs and
  ## reference values are written in hex since R's decimal parser isn't always
  ## correctly rounded (e.g. it reads 244837.89991574993 as
  ## 0x1.de32f33070758p+17 instead of 0x1.de32f33070759p+17); and it parses
  ## subnormal hex literals like 0x1p-1074 as 0, so that lo is written as
  ## 2^-1074.  The cases
  ## exercise direct summation, Temme's expansion, both continued fractions,
  ## and a subnormal result.
  ##
  ## Spacing between consecutive doubles in the binade containing the exact
  ## value hi + lo.
  ulp_of_exact <- function(hi, lo) {
    x <- abs(hi)
    if (x < .Machine$double.xmin) {
      return(2^-1074)
    }
    e <- floor(log2(x))
    if (2^e > x) {
      e <- e - 1
    } else if (2^(e + 1) <= x) {
      e <- e + 1
    }
    if ((x == 2^e) && (sign(lo) == -sign(hi)) && (e > -1022)) {
      ## hi is a power of 2 and the exact value is just below it in magnitude.
      e <- e - 1
    }
    2^(e - 52)
  }
  ## Entries are (q, lambda, lower.tail, log.p, hi, lo).
  ulp_cases <- list(
    list(141, 0x1.26674275d3e35p+7, TRUE, FALSE, 0x1.4acd640bbc5dbp-2, 0x1.ff0dd00803c4fp-56),
    list(12663, 0x1.8c776f1914c96p+13, FALSE, FALSE, 0x1.29e65f0a030f1p-1, -0x1.fcf6333abbe22p-55),
    list(1076091, 0x1.0816fd59c0965p+20, TRUE, FALSE, 0x1.1117594e845e3p-25, -0x1.fcd7607725959p-79),
    list(1846969520, 1846972037, FALSE, FALSE, 0x1.0bf3eb68c07ccp-1, 0x1.f5cd59ee281fcp-55),
    list(67534192, 0x1.01a0a7ec0556bp+26, TRUE, FALSE, 0x1.be51b0cf35846p-2, 0x1.ebbcdbe5e4246p-56),
    list(168122, 0x1.3f44f4c5caecbp+17, FALSE, FALSE, 0x1.4582b85633a33p-100, -0x1.ecaeefc0b9f13p-154),
    list(236634, 0x1.de32f33070759p+17, FALSE, TRUE, -0x1.11d642801bf89p-206, -0x1.ebeee96035b19p-260),
    list(1197, 3000.75, TRUE, FALSE, 0x1.196432abd9e51p-1021, 2^-1074),
    list(960000, 1000000.5, TRUE, TRUE, -0x1.97bf23c291c95p+9, -0x1.db65f41e8705p-45),
    list(1040000, 1000000.5, FALSE, TRUE, -0x1.8d132acf660e3p+9, 0x1.30a5633a6172ap-45))
  for (cc in ulp_cases) {
    hi <- cc[[5]]
    lo <- cc[[6]]
    got <- ppois(cc[[1]], cc[[2]], lower.tail=cc[[3]], log.p=cc[[4]])
    ## (got - hi) is computed exactly whenever got is within a few ULPs of hi.
    expect_lt(abs((got - hi) - lo), ulp_of_exact(hi, lo),
              label=paste0("ppois(", cc[[1]], ", ", sprintf("%a", cc[[2]]), ") error"))
    got_approx <- ppois(cc[[1]], cc[[2]], lower.tail=cc[[3]], log.p=cc[[4]], approx=TRUE)
    expect_equal(got_approx, hi, tolerance=2e-15)
  }
})

test_that("poisson.test works", {
  ## From R 4.6.1 poisson.test() examples.
  expect_snapshot({
    poisson.test(137, 24.19893)
    poisson.test(c(11, 6+8+7), c(800, 1083+1050+878))
  })

  ## Agreement with stats::poisson.test() on ordinary cases, excluding the
  ## integer-(r * T) tie cases where rounding error matters.
  set.seed(1)
  for (i in 1:200) {
    m <- 10^runif(1, -1, 3.5)
    x <- max(0, round(m + rnorm(1) * 3 * sqrt(m)))
    T <- 10^runif(1, -1, 1)
    for (alt in c("two.sided", "less", "greater")) {
      expect_equal(poisson.test(x, T, m/T, alternative=alt)[c("statistic", "parameter", "p.value", "conf.int", "estimate", "null.value", "alternative", "method")],
                   stats::poisson.test(x, T, m/T, alternative=alt)[c("statistic", "parameter", "p.value", "conf.int", "estimate", "null.value", "alternative", "method")],
                   tolerance=1e-12)
    }
  }
  x <- c(7, 13)
  T <- c(1.5, 2.25)
  for (alt in c("two.sided", "less", "greater")) {
    expect_equal(poisson.test(x, T, 0.75, alternative=alt)[c("statistic", "parameter", "p.value", "conf.int", "estimate", "null.value", "alternative", "method")],
                 stats::poisson.test(x, T, 0.75, alternative=alt)[c("statistic", "parameter", "p.value", "conf.int", "estimate", "null.value", "alternative", "method")],
                 tolerance=1e-12)
  }

  ## Ties are exact when (r * T) is an exact integer.
  expect_identical(poisson.test(4, r=5)$p.value, 1)
  expect_identical(poisson.test(5, r=5)$p.value, 1)
  d <- stats::dpois(0:200, 6)
  expect_equal(poisson.test(9, r=6)$p.value, sum(d[d <= stats::dpois(9, 6) * (1 + 1e-12)]))

  ## midp
  m <- 3.3
  d <- stats::dpois(0:200, m)
  expect_equal(poisson.test(7, r=m, midp=TRUE)$p.value,
               sum(d[d < stats::dpois(7, m)]) + 0.5 * stats::dpois(7, m))
  expect_equal(poisson.test(7, r=m, alternative="less", midp=TRUE)$p.value,
               stats::ppois(6, m) + 0.5 * stats::dpois(7, m))
  expect_equal(poisson.test(7, r=m, alternative="greater", midp=TRUE)$p.value,
               stats::ppois(7, m, lower.tail=FALSE) + 0.5 * stats::dpois(7, m))
  ## lambda=2: pmf(1) = pmf(2) is an exact tie between the two modes, and no
  ## other pmf value ties pmf(0).
  d <- stats::dpois(0:100, 2)
  expect_equal(poisson.test(1, r=2, midp=TRUE)$p.value, 1 - stats::dpois(1, 2))
  expect_equal(poisson.test(0, r=2, midp=TRUE)$p.value,
               sum(d[d < 0.999 * stats::dpois(0, 2)]) + 0.5 * stats::dpois(0, 2))
  expect_identical(poisson.test(0, T=0, midp=TRUE)$p.value, 0.5)

  ## log.p
  res <- poisson.test(10, r=2000, log.p=TRUE)
  ## Checked against a 320-bit MPFR implementation.
  expect_equal(res$log.p.value, -1938.6127473836495, tolerance=1e-14)

  ## Large cases, checked against a 320-bit MPFR implementation.
  expect_equal(poisson.test(1000123456, r=1e9)$p.value, 9.4638293091475092e-05, tolerance=1e-14)
  expect_equal(poisson.test(999876544, r=1e9)$p.value, 9.4582633468018055e-05, tolerance=1e-14)
  expect_equal(poisson.test(1000123456, r=1e9, log.p=TRUE)$log.p.value, -9.2654483742611529, tolerance=1e-14)

  expect_error(poisson.test(1, r=2^51))
  expect_error(poisson.test(c(1, 2, 3)))
  expect_error(poisson.test(-1))
})

test_that("dhyper works", {
  ## Based on dhyper() tests in R 4.6.1 tests/d-p-q-r-tst-2.R and
  ## tests/print-tests.R .
  expect_no_error(lapply(sample(10000, size=1000), function(M) {
    n <- (M/100)*10^(2:20)
    if (anyNA(P <- dhyper(n+1,n+5,n+5,n))) {
      stop("NA for M=", M, "; 10ex=",paste((2:20)[is.na(P)], collapse=", "))
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

  ## R PR#4688: once gave p.value=Inf, then FEXACT error 501.  Now handled
  ## by the 3xk code; the p-value underflows to zero, which is detected
  ## quickly.  (log.p = TRUE would be very slow here.)
  reli <- cbind(Si = c(2121, 100, 27, 0),
                av = c(4700, 216, 67, 0),
                Nc = c(6234,2461,502,14))
  expect_identical(fisher.test(reli, workspace=2000000)$p.value, 0)

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

  ## stats::fisher.test() fails on this one (FEXACT hash key overflow, even
  ## with a large workspace); now handled by the network code.  Regression
  ## value from this implementation; a 10^6-replicate simulation gave
  ## 0.63275 +- 0.00048.
  d <- matrix(c(1,0,5,2,1,90
               ,2,1,0,2,3,89
               ,0,0,0,1,0,14
               ,0,0,0,0,0, 5
               ,0,0,0,0,0, 2
               ,0,0,0,0,0, 2
                ), nrow=6, byrow=TRUE)
  expect_equal(fisher.test(d)$p.value, 0.63221605315178575, tolerance = 1e-12)
})

test_that("fisher.test works on 2xk and 3x3 tables", {
  ## Reference values were computed by exact enumeration with rational
  ## arithmetic.
  t24 <- matrix(c(3, 10, 4, 7, 9, 2, 4, 1), nrow = 2, byrow = TRUE)
  expect_equal(fisher.test(t24)$p.value, 0.007673597307624796, tolerance = 1e-12)
  expect_equal(fisher.test(t24, midp = TRUE)$p.value, 0.00715611580507905, tolerance = 1e-12)
  expect_equal(fisher.test(t(t24), midp = TRUE)$p.value, 0.00715611580507905, tolerance = 1e-12)
  t52 <- matrix(c(5, 1, 2, 8, 3, 3, 6, 0, 1, 4), ncol = 2, byrow = TRUE)
  expect_equal(fisher.test(t52)$p.value, 0.004501551251436071, tolerance = 1e-12)
  t33 <- matrix(c(3, 1, 4, 1, 5, 9, 2, 6, 5), nrow = 3, byrow = TRUE)
  expect_equal(fisher.test(t33)$p.value, 0.282360228768354, tolerance = 1e-12)
  expect_equal(fisher.test(t33, midp = TRUE)$p.value, 0.2760652637970296, tolerance = 1e-12)
  ft33 <- fisher.test(matrix(c(8, 0, 1, 0, 7, 2, 1, 1, 9), nrow = 3), log.p = TRUE)
  expect_equal(ft33$log.p.value, log(6.011127956704347e-07), tolerance = 1e-12)
  expect_equal(fisher.test(matrix(2, 3, 3), midp = TRUE)$p.value, 0.9787502413552833, tolerance = 1e-12)

  ## Agreement with stats::fisher.test() (which has a relative tolerance of
  ## 1e-7 when deciding which tables are at least as extreme).
  set.seed(11)
  for (i in 1:5) {
    m <- matrix(rpois(8, 30), nrow = 2)
    expect_equal(fisher.test(m)$p.value, stats::fisher.test(m, workspace = 2e7)$p.value, tolerance = 1e-6)
    m <- matrix(rpois(9, 20), nrow = 3)
    expect_equal(fisher.test(m)$p.value, stats::fisher.test(m, workspace = 2e7)$p.value, tolerance = 1e-6)
  }

  ## Large table; regression value from this implementation.
  big33 <- matrix(c(400, 520, 580, 510, 470, 520, 600, 480, 420), nrow = 3, byrow = TRUE)
  expect_equal(fisher.test(big33, log.p = TRUE)$log.p.value, -31.139221226132207, tolerance = 1e-10)

  ## p-values which underflow (these use the early-exit path when log.p is
  ## FALSE).
  for (m in list(matrix(c(300, 0, 300, 0, 0, 300, 0, 300), nrow = 2, byrow = TRUE),
                 diag(300, 3))) {
    expect_identical(fisher.test(m)$p.value, 0)
    expect_lt(fisher.test(m, log.p = TRUE)$log.p.value, -745.2)
  }

  ## Tables with at least 4 rows and columns support midp now (exact value).
  expect_equal(fisher.test(matrix(1, 4, 4), midp = TRUE)$p.value, exp(-0.0026339782217517601), tolerance = 1e-13)
  expect_identical(fisher.test(matrix(1, 4, 4))$p.value, 1)
  ## Still forwarded to stats::fisher.test() (more than 8 rows and columns).
  expect_error(fisher.test(matrix(1, 9, 9), midp = TRUE))
  expect_error(fisher.test(matrix(1, 2, 9), midp = FALSE, simulate.p.value = TRUE, B = 10), NA)
})

test_that("fisher.test works on larger 2xk and 3xk tables", {
  ## Reference values were computed by (pruned) exact enumeration with
  ## 256-bit MPFR arithmetic.  Each case is (column-major entries, number of
  ## rows, log-p-value, log-midp-value).
  cases <- list(
    list(c(59, 64, 79, 94, 38, 45, 90, 86, 82, 82, 60, 72, 70, 79), 2,
         -0.071938761276373928, -0.071938815194004191),
    list(c(69, 85, 66, 75, 90, 86, 28, 37, 68, 70, 90, 84, 77, 75), 2,
         -0.25190433769013532, -0.25190440603994091),
    list(c(1, 2, 1, 4, 3, 1, 2, 5, 0, 1, 0, 1, 0, 2, 3, 4, 3, 2), 2,
         -0.35329748516790693, -0.41242615795654952),
    list(c(1, 1, 3, 1, 1, 1, 1, 0, 0, 2, 3, 1, 3, 7, 5, 2, 3, 0), 2,
         -1.4843038728469717, -1.4918907420522329),
    list(c(1, 1, 0, 1, 1, 0, 3, 0, 1, 0, 1, 1, 2, 2, 2, 2, 3, 1), 3,
         -0.11365018305619758, -0.17781090707359118),
    list(c(4, 1, 1, 0, 0, 1, 0, 2, 1, 1, 0, 2, 1, 0, 3, 1, 3, 1), 3,
         -1.4070116119685145, -1.5536083327901693),
    list(c(0, 1, 0, 1, 2, 0, 4, 1, 2, 0, 6, 1, 7, 4, 1), 3,
         -2.8640803039982159, -2.9425671320127691),
    list(c(3, 0, 0, 2, 1, 2, 3, 0, 0, 1, 1, 1, 1, 1, 2), 3,
         -0.54059281778090373, -0.64751570884718135),
    list(c(1, 0, 2, 0, 0, 1, 1, 2, 2, 4, 1, 5), 3,
         -0.18301893558768439, -0.21874741757827967))
  for (cc in cases) {
    m <- matrix(cc[[1]], nrow = cc[[2]])
    for (mm in list(m, t(m))) {
      expect_equal(fisher.test(mm)$p.value, exp(cc[[3]]), tolerance = 1e-13)
      expect_equal(fisher.test(mm, log.p = TRUE)$log.p.value, cc[[3]], tolerance = 1e-13)
      expect_equal(fisher.test(mm, midp = TRUE)$p.value, exp(cc[[4]]), tolerance = 1e-13)
    }
  }

  ## Agreement with stats::fisher.test().
  set.seed(13)
  for (i in 1:3) {
    m <- matrix(rpois(12, 8), nrow = 3)
    expect_equal(fisher.test(m)$p.value, stats::fisher.test(m, workspace = 2e7)$p.value, tolerance = 1e-6)
    m <- matrix(rpois(16, 6), nrow = 2)
    expect_equal(fisher.test(m)$p.value, stats::fisher.test(m, workspace = 2e7)$p.value, tolerance = 1e-6)
  }

  ## Too large for stats::fisher.test() with its default workspace.
  big <- matrix(c(500, 400, 300, 200, 100, 260, 210, 160, 110, 60), nrow = 2)
  expect_gt(fisher.test(big)$p.value, 0)
})

test_that("fisher.test works on tables with 4+ rows and columns", {
  ## Reference values were computed by exact rational enumeration.  Each case
  ## is (row-major entries, number of rows, log-p-value, log-midp-value).
  cases <- list(
    list(c(3, 0, 1, 0, 0, 2, 1, 1, 1, 0, 0, 3, 0, 2, 1, 0), 4,
         -3.1065721504237533, -3.4441981430055453),
    list(c(2, 1, 0, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 1, 0, 2, 1, 0, 1, 1), 5,
         -0.4361435607576663, -0.67329362509793245))
  for (cc in cases) {
    m <- matrix(cc[[1]], nrow = cc[[2]], byrow = TRUE)
    for (mm in list(m, t(m))) {
      expect_equal(fisher.test(mm, log.p = TRUE)$log.p.value, cc[[3]], tolerance = 1e-13)
      expect_equal(fisher.test(mm, midp = TRUE)$p.value, exp(cc[[4]]), tolerance = 1e-13)
    }
  }
  ## Agreement with stats::fisher.test() (which has a relative tolerance of
  ## ~1e-7 for ties).
  set.seed(17)
  for (i in 1:3) {
    m <- matrix(rpois(16, 6), nrow = 4)
    expect_equal(fisher.test(m)$p.value, stats::fisher.test(m, workspace = 2e7)$p.value, tolerance = 1e-6)
    m <- matrix(rpois(20, 4), nrow = 5)
    expect_equal(fisher.test(m)$p.value, stats::fisher.test(m, workspace = 2e7)$p.value, tolerance = 1e-6)
  }
  ## Empty rows/columns are dropped; data.name is preserved.
  m0 <- rbind(matrix(c(3, 0, 1, 0, 0, 2, 1, 1, 1, 0, 0, 3, 0, 2, 1, 0), nrow = 4, byrow = TRUE), 0)
  expect_equal(fisher.test(m0, log.p = TRUE)$log.p.value, -3.1065721504237533, tolerance = 1e-13)
})

test_that("fisher.test workspace limits the 3xk network algorithm", {
  ## The network's limit is max(4 * workspace, 256 MiB); Inf removes it.
  ## Results must agree.
  for (m in list(matrix(c(102, 65, 39, 19, 59, 42, 34, 30, 32, 26, 25, 27), nrow = 3, byrow = TRUE),
                 matrix(c(23, 14, 8, 3, 0, 9, 7, 6, 6, 6, 3, 3, 4, 5, 7), nrow = 3, byrow = TRUE))) {
    p0 <- fisher.test(m, log.p = TRUE)$log.p.value
    expect_equal(fisher.test(m, workspace = 1000, log.p = TRUE)$log.p.value, p0, tolerance = 1e-12)
    expect_equal(fisher.test(t(m), workspace = Inf, log.p = TRUE)$log.p.value, p0, tolerance = 1e-12)
    expect_equal(fisher.test(m, workspace = Inf, midp = TRUE)$p.value, fisher.test(m, midp = TRUE)$p.value, tolerance = 1e-12)
  }
  expect_error(fisher.test(matrix(1:9, 3), workspace = -1), "workspace")
  expect_error(fisher.test(matrix(1:9, 3), workspace = NA), "workspace")
})

test_that("fisher.test respects time.limit", {
  ## Near-modal 2x4 and 3x3 tables with huge totals take many minutes.
  m <- 2500000
  d <- m %/% 250
  x24 <- matrix(c(m, m, m, m, m - d, m + d, m + d / 2, m - d / 2), nrow = 2,
                byrow = TRUE)
  expect_error(fisher.test(x24, time.limit = 0.05), "time.limit")
  expect_error(fisher.test(t(x24), time.limit = 0.05), "time.limit")
  m <- 3000000
  d <- m %/% 1000
  x33 <- matrix(c(m, m, m, m, m, m, m, m, m - d), nrow = 3, byrow = TRUE)
  expect_error(fisher.test(x33, time.limit = 0.05), "time.limit")
  ## Fast tables are unaffected by the limit.
  x35 <- matrix(c(10, 12, 8, 15, 9, 7, 11, 13, 6, 10, 9, 8, 12, 7, 11),
                nrow = 3)
  expect_identical(fisher.test(x35, time.limit = Inf)$p.value,
                   fisher.test(x35)$p.value)
  for (bad in list(0, -1, NA, "60", c(1, 2))) {
    expect_error(fisher.test(x35, time.limit = bad), "time.limit")
  }
})
