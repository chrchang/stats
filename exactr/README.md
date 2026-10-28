exactr
======

This package contains accurate and efficient implementations of the binomial
and Fisher's (2x2..2x6 and 3x3 so far) exact tests.

It also has unexported implementations of the related discrete distribution
functions (dbinom, dhyper, pbinom, qbinom, phyper, qhyper), which are slower
but more accurate than their R standard library counterparts.

### Comparison with R 4.6.1 stats

The MPFR-comparison (requires the Rmpfr package) and benchmark scripts under
utils/ in the GitHub repository (https://github.com/chrchang/stats ) compare
these functions against their R 4.6.1 stats counterparts.  Errors are RMS
relative errors vs. a 256-bit MPFR implementation; n in [2^35, 2^36) MPFR
tiers are omitted since Rmpfr is too slow there.

stats::binom.test()'s two-sided p-value is accurate, but it takes O(n) time
and memory; it is skipped for n >= 2^27:

    $ utils/binomtest_accuracy.R --z-score 1 -e 5,20,25,30
    n in [2^5, 2^6): errRMS=1.15e-16  statsErrRMS=1.79e-15
    n in [2^20, 2^21): errRMS=7.36e-16  statsErrRMS=6.62e-16
    n in [2^25, 2^26): errRMS=4.85e-16  statsErrRMS=8e-16
    n in [2^30, 2^31): errRMS=7.28e-16
    $ utils/binomtest_benchmark.R --z-score 1  # ~650x speedup at n~=1m
    n=(2^5)-1: base=2.92e-05  stats=2.7e-05 sec/iter
    n=(2^20)-1: base=2.86e-05  stats=0.0189 sec/iter
    n=(2^35)-1: base=3.36e-05 sec/iter

For 2x2 tables, stats::fisher.test()'s p-value is accurate, but it takes O(n)
time and memory (and requires entries < 2^31); it is skipped for n >= 2^27.
Its odds-ratio estimate and confidence interval are only accurate to ~5-7
significant digits, since they come from uniroot() with default tolerance:

    $ utils/fisher_exact_22_accuracy.R --z-score -1 -e 5,20,25,30
    n in [2^5, 2^6): errRMS=7.56e-17  statsErrRMS=3.15e-16
    n in [2^20, 2^21): errRMS=1.48e-15  statsErrRMS=3.28e-15
    n in [2^25, 2^26): errRMS=1.55e-14  statsErrRMS=5.33e-15
    n in [2^30, 2^31): errRMS=6.47e-14
    $ utils/fisher_exact_22_benchmark.R --z-score -1  # conf.int=FALSE; ~1200x speedup at n~=1m
    n=(2^5)-1: base=0.000102  stats=9.13e-05 sec/iter
    n=(2^20)-1: base=7.14e-05  stats=0.0851 sec/iter
    n=(2^35)-1: base=0.00189 sec/iter

    $ utils/odds_ratio_accuracy.R --z-score 1 -e 5,10,15,20
    n in [2^5, 2^6): estErrRMS=4.58e-16  ciLowErrRMS=8.87e-16  ciHighErrRMS=6.81e-16  statsEstErrRMS=1.82e-05  statsCiLowErrRMS=2.06e-05  statsCiHighErrRMS=8.15e-05
    n in [2^10, 2^11): estErrRMS=1.43e-16  ciLowErrRMS=2.01e-15  ciHighErrRMS=1.94e-15  statsEstErrRMS=1.21e-06  statsCiLowErrRMS=1.21e-05  statsCiHighErrRMS=9.96e-06
    n in [2^15, 2^16): estErrRMS=1.2e-16  ciLowErrRMS=8.32e-15  ciHighErrRMS=8.37e-15  statsEstErrRMS=5.28e-06  statsCiLowErrRMS=1.16e-05  statsCiHighErrRMS=9.57e-06
    n in [2^20, 2^21): estErrRMS=7e-17  ciLowErrRMS=4.57e-14  ciHighErrRMS=4.57e-14  statsEstErrRMS=1.59e-07  statsCiLowErrRMS=2.41e-05  statsCiHighErrRMS=1.72e-05
    $ utils/odds_ratio_benchmark.R --z-score 1  # conf.int=TRUE
    n=(2^5)-1: base=0.000103  stats=0.000171 sec/iter
    n=(2^10)-1: base=6.49e-05  stats=0.000405 sec/iter
    n=(2^15)-1: base=6.95e-05  stats=0.00683 sec/iter
    n=(2^35)-1: base=0.00362 sec/iter
