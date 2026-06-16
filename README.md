stats
=====

Implementations of the binomial, Hardy-Weinberg equilibrium,
and Fisher's exact tests (2x2 and 2x3 so far), along with corresponding
discrete distributions (binomial, hypergeometric).

As of this writing, these functions are more accurate, and often simultaneously
more efficient, than their scipy counterparts; see the MPFR-comparison and
benchmark scripts under Python/utils/ .  (R package to be provided soon.)

The central building block is a high-precision log-factorial function utilizing
the QD library (https://github.com/BL-highprecision/QD ).


Test program compilation:

make binom_demo
make fisher_demo
make hwe_demo

Test program usage examples:

    $ ./binom_demo 2900 5000 0.52
    Two-sided p-value: 1.86026e-17
    $ ./binom_demo 2900000000000 5000000000000 0.52
    Two-sided p-value: 1.04869e-15748396144
    $ ./fisher_demo 20 30 40 50 60 70
    P-value: 0.780903
    $ ./hwe_demo 16 3 81
    P-value: 0.0899446

Licensing:
- include/plink2_highprec's adaptation of a subset of the QD library is under
  the include/LICENSE.QD BSD-3-Clause-LBNL license.
- Other code is covered by LGPL-3.0-only (see the COPYING.LESSER file).
