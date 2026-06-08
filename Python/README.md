exact_tests
===========

This package provides binomial, Hardy-Weinberg equilibrium, and Fisher's exact
test (2x2 and 2x3 so far) functions, along with binomial and hypergeometric
distribution functions with scipy- and R-style interfaces.  The functions are
accurate and efficient:

- High-precision and interval arithmetic are used to ensure likelihood
  near-ties are correctly resolved.  Log-(mid)p-values are available when you
  don't want to just round tiny p-values down to zero.

- It's relatively expensive to compute a contingency table's likelihood from
  scratch, but it's cheap when the likelihood of an adjacent contingency table
  is known.  And it's even cheaper to skip a likelihood calculation when the
  value must be too small to matter.  These functions take advantage of the
  latter possibilities.


### Build instructions
PyPI:
```
pip install 'pip>=20.3'
pip install exact_tests
```

GitHub:
```
pip install 'pip>=20.3'
pip install -e 'git+https://github.com/chrchang/stats.git#egg=exact_tests&subdirectory=Python'
```

Or install from a cloned copy:
```
# clone repo
git clone https://github.com/chrchang/stats
# go to python folder
cd stats/Python
# install the package
pip install -e .
```

You can test the package with `pytest`.

##### Example usage:
```
import exact_tests

exact_tests.binom_test(3, n=15, p=0.1, alternative="greater")
exact_tests.binom_test(2, 29, 0.1)
exact_tests.binom_test(2, 29, "0.1")
exact_tests.binom_test(2, 29, "0.1", midp=True)
exact_tests.binom_test(100, 10000, "0.000001", logp=True)
exact_tests.fisher_test([[4, 0], [0, 4]], alternative="greater", midp=True)
exact_tests.fisher_test([[10000, 20000], [30000, 40000], [50000, 60000]], logp=True)
```
