exact_tests
===========

This package implements the binomial, Hardy-Weinberg equilibrium, and Fisher's
exact tests (2x2 and 2x3 so far), along with corresponding discrete
distributions (binomial, hypergeometric).  scipy- and R-style interfaces are
provided.

As of this writing, these functions are more accurate, and often simultaneously
more efficient, than their scipy counterparts: see the MPFR-comparison and
benchmark scripts under utils/ .  E.g.

    $ utils/pbinom_accuracy.py
    n in [2^5, 2^6): errRMS=0  approxErrRMS=7.30607e-17  scipyErrRMS=0
    n in [2^20, 2^21): errRMS=0  approxErrRMS=1.25261e-15  scipyErrRMS=2.18032e-15
    n in [2^35, 2^36): errRMS=0  approxErrRMS=6.76522e-15  scipyErrRMS=1.15131e-14
    $ utils/pbinom_accuracy.py --z-score -2
    n in [2^5, 2^6): errRMS=6.50688e-17  approxErrRMS=2.47829e-16  scipyErrRMS=1.96176e-16
    n in [2^20, 2^21): errRMS=0  approxErrRMS=5.0815e-16  scipyErrRMS=1.01243e-13
    n in [2^35, 2^36): errRMS=0  approxErrRMS=6.87091e-16  scipyErrRMS=1.66324e-11
    $ utils/pbinom_accuracy.py -p 0.1
    n in [2^5, 2^6): errRMS=0  approxErrRMS=1.0628e-16  scipyErrRMS=4.22428e-16
    n in [2^20, 2^21): errRMS=0  approxErrRMS=1.62409e-15  scipyErrRMS=1.23353e-11
    n in [2^35, 2^36): errRMS=0  approxErrRMS=5.48266e-15  scipyErrRMS=3.46685e-07
    $ utils/pbinom_benchmark.py
    n=(2^5)-1: base=6.91666e-06  approx=5.5534e-07  scipy=3.8514e-05 sec/iter
    n=(2^20)-1: base=0.000119333  approx=3.62502e-06  scipy=3.55417e-05 sec/iter
    n=(2^35)-1: base=0.00115099  approx=7.9375e-05  scipy=9.75833e-05 sec/iter
    $ utils/binomtest_accuracy.py -p 0.7 --z-score 0.5
    n in [2^5, 2^6): errRMS=1.19694e-16  scipyErrRMS=1.53967e-16
    n in [2^20, 2^21): errRMS=6.4373e-16  scipyErrRMS=1.26149e-11
    n in [2^35, 2^36): errRMS=1.48472e-15  scipyErrRMS=3.29018e-06
    $ utils/binomtest_benchmark.py -p 0.7 --z-score 0.5
    n=(2^5)-1: base=4.76401e-06  scipy=0.000209583 sec/iter
    n=(2^20)-1: base=6.90267e-06  scipy=0.000506403 sec/iter
    n=(2^35)-1: base=1.8389e-05  scipy=0.000786333 sec/iter
    $ utils/phyper_accuracy.py
    n in [2^5, 2^6): errRMS=0  approxErrRMS=1.27832e-16  scipyErrRMS=3.26681e-16
    n in [2^20, 2^21): errRMS=0  approxErrRMS=7.5947e-16  scipyErrRMS=2.76054e-10
    n in [2^35, 2^36): errRMS=0  approxErrRMS=1.64e-14  scipyErrRMS=1.27419e-05
    $ utils/phyper_benchmark.py
    n=(2^5)-1: base=6.63901e-06  approx=4.72336e-07  scipy=4.50973e-05 sec/iter
    n=(2^20)-1: base=6.7625e-05  approx=5.44434e-06  scipy=0.000506611 sec/iter
    n=(2^35)-1: base=0.0146301  approx=0.000858139  scipy=14.8541 sec/iter

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
