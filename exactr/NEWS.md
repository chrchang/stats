# exactr 0.1.2 (2026-11-18)
* New poisson.test().
* binom.test(alternative="greater", midp=TRUE) now returns
  P[X > x] + 0.5 * P[X = x]; it previously returned
  P[X >= x-1] - 0.5 * P[X = x-1].
* binom.test(alternative="greater") now returns 1 when x = 0 (previously 0).
* binom.test(midp=TRUE) now reports method "Exact binomial test (midp)".

# exactr 0.1.0 (2026-10-28)
Initial CRAN submission.
