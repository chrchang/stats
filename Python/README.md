exact_tests
===========

This package implements the binomial, Hardy-Weinberg equilibrium, and Fisher's
exact tests (2x2 and 2x3 so far), along with corresponding discrete
distributions (binomial, hypergeometric).  scipy- and R-style interfaces are
provided.  cdf-like functions offer a fast 'approx' mode.

As of this writing, these functions (even in 'approx' mode) are more accurate,
and often simultaneously more efficient, than their scipy counterparts: see the
MPFR-comparison and benchmark scripts under utils/ .  E.g.

    $ utils/pbinom_accuracy.py
    n in [2^5, 2^6): errRMS=0  approxErrRMS=7.31e-17  scipyErrRMS=0
    n in [2^20, 2^21): errRMS=0  approxErrRMS=1.25e-15  scipyErrRMS=2.18e-15
    n in [2^35, 2^36): errRMS=0  approxErrRMS=6.77e-15  scipyErrRMS=1.15e-14
    $ utils/pbinom_accuracy.py -p 0.1  # scipy accuracy degrades more quickly with most values of p != 0.5, or k further from np
    n in [2^5, 2^6): errRMS=0  approxErrRMS=1.06e-16  scipyErrRMS=4.22e-16
    n in [2^20, 2^21): errRMS=0  approxErrRMS=1.62e-15  scipyErrRMS=1.23e-11
    n in [2^35, 2^36): errRMS=0  approxErrRMS=5.48e-15  scipyErrRMS=3.47e-07
    $ utils/pbinom_accuracy.py --z-score -2
    n in [2^5, 2^6): errRMS=6.51e-17  approxErrRMS=2.48e-16  scipyErrRMS=1.96e-16
    n in [2^20, 2^21): errRMS=0  approxErrRMS=5.08e-16  scipyErrRMS=1.01e-13
    n in [2^35, 2^36): errRMS=0  approxErrRMS=6.87e-16  scipyErrRMS=1.66e-11
    $ utils/pbinom_benchmark.py  # yes, default mode becomes slower than scipy for large n, that's why approx= exists
    n=(2^5)-1: base=6.92e-07  approx=2.92e-07  scipy=2.99e-05 sec/iter
    n=(2^20)-1: base=4.6e-05  approx=3.52e-06  scipy=2.99e-05 sec/iter
    n=(2^35)-1: base=0.00113  approx=7.96e-05  scipy=9.5e-05 sec/iter

    $ utils/binomtest_accuracy.py --z-score 1
    n in [2^5, 2^6): errRMS=1.44e-16  scipyErrRMS=1.96e-16
    n in [2^20, 2^21): errRMS=5.4e-16  scipyErrRMS=6e-14
    n in [2^35, 2^36): errRMS=5.83e-16  scipyErrRMS=7.53e-12
    $ utils/binomtest_benchmark.py --z-score 1  # ~100-200x speedup, better accuracy
    n=(2^5)-1: base=8.42e-07  scipy=0.000189 sec/iter
    n=(2^20)-1: base=4.03e-06  scipy=0.000464 sec/iter
    n=(2^35)-1: base=5.67e-06  scipy=0.000774 sec/iter

    $ utils/phyper_accuracy.py --z-score -0.5  # scipy accuracy degrades quickly
    n in [2^5, 2^6): errRMS=0  approxErrRMS=1.72e-16  scipyErrRMS=2.84e-16
    n in [2^20, 2^21): errRMS=0  approxErrRMS=3.01e-15  scipyErrRMS=1.73e-10
    n in [2^35, 2^36): errRMS=0  approxErrRMS=3.93e-13  scipyErrRMS=1.05e-05
    $ utils/phyper_benchmark.py --z-score -0.5  # scipy speed is ok, except...
    n=(2^5)-1: base=7.54e-07  approx=2.08e-07  scipy=3.51e-05 sec/iter
    n=(2^20)-1: base=6.98e-05  approx=5.37e-06  scipy=4.35e-05 sec/iter
    n=(2^35)-1: base=0.00697  approx=0.000401  scipy=0.00195 sec/iter
    $ utils/phyper_benchmark.py --z-score 0.0001  # ...it blows up for large n when z approaches 0.
    n=(2^5)-1: base=7.12e-07  approx=2.12e-07  scipy=3.51e-05 sec/iter
    n=(2^20)-1: base=6.88e-05  approx=5.32e-06  scipy=0.000498 sec/iter
    n=(2^35)-1: base=0.0147  approx=0.000873  scipy=2.69 sec/iter
    $ utils/phyper_benchmark.py
    n=(2^5)-1: base=7.33e-07  approx=2.12e-07  scipy=3.49e-05 sec/iter
    n=(2^20)-1: base=6.87e-05  approx=5.32e-06  scipy=0.000492 sec/iter
    n=(2^35)-1: base=0.0147  approx=0.000869  scipy=14.8 sec/iter

    $ utils/fisher_exact_22_accuracy.py --z-score -1
    n in [2^5, 2^6): errRMS=9.71e-17  scipyErrRMS=2.39e-16
    n in [2^20, 2^21): errRMS=3.17e-15  scipyErrRMS=1.85e-10
    n in [2^35, 2^36): errRMS=6.63e-13
    $ utils/fisher_exact_22_benchmark.py --z-score -1  # ~75-250x speedup, better accuracy
    n=(2^5)-1: base=1.07e-06  scipy=0.000258 sec/iter
    n=(2^20)-1: base=8.34e-06  scipy=0.000667 sec/iter
    n=(2^35)-1: base=0.000742 sec/iter

    $ utils/odds_ratio_concordance.py --z-score 1  # we have high agreement with scipy
    n in [2^5, 2^6): estErrRMS=1.73e-15  ciLowErrRMS=2.32e-16  ciHighErrRMS=4.42e-16
    n in [2^10, 2^11): estErrRMS=1.15e-14  ciLowErrRMS=3.46e-16  ciHighErrRMS=2.83e-15
    n in [2^15, 2^16): estErrRMS=1.52e-12  ciLowErrRMS=1.6e-12  ciHighErrRMS=1.26e-14
    $ utils/odds_ratio_benchmark.py --z-score 1  # 200x and larger speedups
    n=(2^5)-1: est=6.38e-07  ci=1.07e-06  scipy_est=0.000445  scipy_ci=0.00257 sec/iter
    n=(2^10)-1: est=6.79e-07  ci=1.79e-06  scipy_est=0.00051  scipy_ci=0.0281 sec/iter
    n=(2^15)-1: est=2.98e-06  ci=6.4e-06  scipy_est=0.000702  scipy_ci=2.27 sec/iter
    n=(2^35)-1: est=0.000859  ci=0.00176 sec/iter

    $ utils/HWE_accuracy.py --maf 0.05 --z-score 1
    n in [2^10, 2^11): errRMS=1.95e-16  snphweErrRMS=2.04e-16
    n in [2^15, 2^16): errRMS=2.9e-16  snphweErrRMS=6.48e-16
    n in [2^20, 2^21): errRMS=6.25e-16  snphweErrRMS=1.6e-15
    $ utils/HWE_benchmark.py --maf 0.05 --z-score 1  # snphwe.snphwe() is accurate, but slow for biobank-scale cases (we're >200x as fast for n~=1m)
    n=(2^10)-1: base=2.83e-07  snphwe=8.08e-07 sec/iter
    n=(2^15)-1: base=4.5e-07  snphwe=1.24e-05 sec/iter
    n=(2^20)-1: base=1.48e-06  snphwe=0.000388 sec/iter

The central building block is a high-precision log-factorial function utilizing
the QD library (https://github.com/BL-highprecision/QD ).

### Build instructions
PyPI:
```
python -m pip install 'pip>=20.3'
python -m pip install exact_tests
```

GitHub:
```
python -m pip install 'pip>=20.3'
python -m pip install -e 'git+https://github.com/chrchang/stats.git#egg=exact_tests&subdirectory=Python'
```

Or install from a cloned copy:
```
# clone repo
git clone https://github.com/chrchang/stats
# go to python folder
cd stats/Python
# install the package
python -m pip install -e .
```

You can test the package with `pytest`.

##### Example usage:
```
import exact_tests

exact_tests.binomtest(3, n=15, p=0.1, alternative="greater")
exact_tests.binomtest(2, 29, 0.1)
exact_tests.binomtest(2, 29, "0.1")
exact_tests.binomtest(2, 29, "0.1", midp=True)
exact_tests.binomtest(100, 10000, "0.000001", logp=True)
exact_tests.fisher_exact([[4, 0], [0, 4]], alternative="greater", midp=True)
exact_tests.fisher_exact([[10000, 20000], [30000, 40000], [50000, 60000]], logp=True)
```
