stats
=====

Standalone functions for calculating binomial, Hardy-Weinberg equilibrium, and
Fisher 2x2/2x3 exact test statistics.

Test program compilation:

gcc -O2 binom.c binom_test.c -o binom_test
make fisher_demo
make hwe_demo

Test program usage examples:

  ./binom_test binom_test.txt

  ./binom_test 2900 5000 0.52

  ./fisher_demo fisher_demo.txt

  ./fisher_demo 20 30 40 50 60 70

  ./hwe_demo hwe_demo.txt

  ./hwe_demo 16 3 81

Licensing:
- include/plink2_highprec's adaptation of a small subset of the QD library
  (https://github.com/BL-highprecision/QD ) is under the include/LICENSE.QD
  BSD-3-Clause-LBNL license.
- binom.c, fisher_demo.cc, hwe_demo.cc, and all of include/ except the
  aforementioned bit of plink2_highprec are covered by LGPLv3 (see the
  COPYING.LESSER file).
- mini-gmp/ is dual-licensed under LGPLv3 and GPLv2.
