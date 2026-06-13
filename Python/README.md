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
    n in [2^5, 2^6): RMS=0.0  approxRMS=7.306071560616266e-17  scipyRMS=0.0
    n in [2^20, 2^21): RMS=0.0  approxRMS=1.2526064965951066e-15  scipyRMS=2.1803173153679657e-15
    n in [2^35, 2^36): RMS=0.0  approxRMS=6.765221412952766e-15  scipyRMS=1.151310810012349e-14
    $ utils/pbinom_accuracy.py --z-score -2
    n in [2^5, 2^6): RMS=6.50688277690206e-17  approxRMS=2.4782896885670045e-16  scipyRMS=1.961761908748987e-16
    n in [2^20, 2^21): RMS=0.0  approxRMS=5.081502940565382e-16  scipyRMS=1.0124311233166742e-13
    n in [2^35, 2^36): RMS=0.0  approxRMS=6.870911531925403e-16  scipyRMS=1.6632421178076063e-11
    $ utils/pbinom_accuracy.py -p 0.1
    n in [2^5, 2^6): RMS=0.0  approxRMS=1.0628000483584193e-16  scipyRMS=4.2242773368500245e-16
    n in [2^20, 2^21): RMS=0.0  approxRMS=1.6240941533512943e-15  scipyRMS=1.2335258925049512e-11
    n in [2^35, 2^36): RMS=0.0  approxRMS=5.482663890316387e-15  scipyRMS=3.466853038093633e-07
    $ utils/pbinom_benchmark.py
    n=(2^5)-1: base=9.99998883344233e-07  approx=2.582994056865573e-07  scipy=3.0166699434630573e-05 sec/iter
    n=(2^20)-1: base=4.6941702021285894e-05  approx=3.274998744018376e-06  scipy=2.860420208889991e-05 sec/iter
    n=(2^35)-1: base=0.0011365917016519234  approx=7.932920125313103e-05  scipy=9.544169879518449e-05 sec/iter
    $ utils/phyper_accuracy.py
    n in [2^5, 2^6): RMS=0.0  approxRMS=1.2783240627107975e-16  scipyRMS=3.2668147214494056e-16
    n in [2^20, 2^21): RMS=0.0  approxRMS=7.594696303752788e-16  scipyRMS=2.7605400523075697e-10
    n in [2^35, 2^36): RMS=0.0  approxRMS=1.639999166953556e-14  scipyRMS=1.2741882006508117e-05
    $ utils/phyper_benchmark.py
    n=(2^5)-1: base=2.847324746350447e-06  approx=4.306687818219264e-07  scipy=4.3361001492788397e-05 sec/iter
    n=(2^20)-1: base=7.218066214894255e-05  approx=5.430668049181501e-06  scipy=0.0004969029978383332 sec/iter
    n=(2^35)-1: base=0.014636972337029874  approx=0.0008547363298324248  scipy=14.766715555665238 sec/iter

The primary building block is a high-precision log-factorial function utilizing
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
