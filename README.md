stats
=====

Functions for calculating binomial, Hardy-Weinberg equilibrium, and Fisher's
2x2/2x3 exact test log-p and log-midp values.  They are accurate and efficient:

- High-precision and interval arithmetic are used to ensure likelihood
  near-ties are correctly resolved.  Log-(mid)p-values are available when you
  don't want to just round tiny p-values down to zero.

- It's relatively expensive to compute a contingency table's likelihood from
  scratch, but it's cheap when the likelihood of an adjacent contingency table
  is known.  And it's even cheaper to skip a likelihood calculation when the
  value must be too small to matter.  These functions take advantage of the
  latter possibilities.


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
- include/plink2_highprec's adaptation of a small subset of the QD library
  (https://github.com/BL-highprecision/QD ) is under the include/LICENSE.QD
  BSD-3-Clause-LBNL license.
- binom_demo.cc, fisher_demo.cc, hwe_demo.cc, and all of include/ except the
  aforementioned bit of plink2_highprec are covered by LGPLv3 (see the
  COPYING.LESSER file).
- mini-gmp/ is dual-licensed under LGPLv3 and GPLv2.
