stats
=====

Implementations of the binomial, Hardy-Weinberg equilibrium,
and Fisher's exact tests (2x2 and 2x3 so far), along with corresponding
discrete distributions (binomial, hypergeometric).

As of this writing, these functions are more accurate, and often simultaneously
more efficient, than their scipy counterparts; see the MPFR-comparison and
benchmark scripts under Python/utils/ .  (R package to be provided soon.)

The primary building block is a high-precision log-factorial function utilizing
the QD library (https://github.com/BL-highprecision/QD ).


Test program compilation:

make binom_demo
make fisher_demo
make hwe_demo

Test program usage examples:

  ./binom_demo binom_demo.txt

  ./binom_demo 2900 5000 0.52

  ./fisher_demo fisher_demo.txt

  ./fisher_demo 20 30 40 50 60 70

  ./hwe_demo hwe_demo.txt

  ./hwe_demo 16 3 81

Licensing:
- include/plink2_highprec's adaptation of a subset of the QD library is under
  the include/LICENSE.QD BSD-3-Clause-LBNL license.
- Other code is covered by LGPL-3.0-only (see the COPYING.LESSER file).
