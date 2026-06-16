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
    n=(2^5)-1: base=2.17e-06  approx=4.86e-07  scipy=3.73e-05 sec/iter
    n=(2^20)-1: base=4.8e-05  approx=3.64e-06  scipy=3.33e-05 sec/iter
    n=(2^35)-1: base=0.00114  approx=7.95e-05  scipy=9.98e-05 sec/iter

    $ utils/binomtest_accuracy.py --z-score 1
    n in [2^5, 2^6): errRMS=1.44e-16  scipyErrRMS=1.96e-16
    n in [2^20, 2^21): errRMS=5.4e-16  scipyErrRMS=6e-14
    n in [2^35, 2^36): errRMS=5.83e-16  scipyErrRMS=7.53e-12
    $ utils/binomtest_benchmark.py --z-score 1  # ~50-100x speedup, better accuracy
    n=(2^5)-1: base=3.72e-06  scipy=0.00021 sec/iter
    n=(2^20)-1: base=5.17e-06  scipy=0.000478 sec/iter
    n=(2^35)-1: base=7.36e-06  scipy=0.000775 sec/iter

    $ utils/phyper_accuracy.py --z-score -0.5  # scipy accuracy degrades quickly
    n in [2^5, 2^6): errRMS=0  approxErrRMS=1.72e-16  scipyErrRMS=2.84e-16
    n in [2^20, 2^21): errRMS=0  approxErrRMS=3.01e-15  scipyErrRMS=1.73e-10
    n in [2^35, 2^36): errRMS=0  approxErrRMS=3.93e-13  scipyErrRMS=1.05e-05
    $ utils/phyper_benchmark.py --z-score -0.5  # scipy speed is ok, except...
    n=(2^5)-1: base=2.74e-06  approx=4.31e-07  scipy=4.45e-05 sec/iter
    n=(2^20)-1: base=6.76e-05  approx=5.51e-06  scipy=5.36e-05 sec/iter
    n=(2^35)-1: base=0.0069  approx=0.000402  scipy=0.00195 sec/iter
    $ utils/phyper_benchmark.py --z-score 0.0001  # ...it blows up for large n when z approaches 0.
    n=(2^5)-1: base=2.53e-06  approx=4.44e-07  scipy=4.21e-05 sec/iter
    n=(2^20)-1: base=7.16e-05  approx=5.67e-06  scipy=0.000502 sec/iter
    n=(2^35)-1: base=0.0146  approx=0.000928  scipy=2.69 sec/iter
    $ utils/phyper_benchmark.py
    n=(2^5)-1: base=6.37e-05  approx=4.17e-07  scipy=4.57e-05 sec/iter
    n=(2^20)-1: base=7.15e-05  approx=5.78e-06  scipy=0.000509 sec/iter
    n=(2^35)-1: base=0.0146  approx=0.000854  scipy=14.7 sec/iter

    $ utils/fisher_exact_22_accuracy.py --z-score -1
    n in [2^5, 2^6): errRMS=9.71e-17  scipyErrRMS=2.39e-16
    n in [2^20, 2^21): errRMS=3.17e-15  scipyErrRMS=1.85e-10
    n in [2^35, 2^36): errRMS=6.63e-13
    $ utils/fisher_exact_22_benchmark.py --z-score -1  # ~60-75x speedup, better accuracy
    n=(2^5)-1: base=4.4e-06  scipy=0.000284 sec/iter
    n=(2^20)-1: base=9.21e-06  scipy=0.00069 sec/iter
    n=(2^35)-1: base=0.000743 sec/iter

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
