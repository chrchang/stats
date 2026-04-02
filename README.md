stats
=====

Standalone functions for calculating binomial, Hardy-Weinberg equilibrium, and
Fisher 2x2/2x3 exact test statistics.

Test program compilation:

gcc -O2 binom.c binom_test.c -o binom_test

make hwe_test

gcc -O2 fisher.c fisher_test.c -o fisher_test

Test program usage examples:

  ./binom_test binom_test.txt

  ./binom_test 2900 5000 0.52

  ./fisher_test fisher_test.txt

  ./fisher_test 20 30 40 50 60 70

  ./hwe_test hwe_test.txt

  ./hwe_test 16 3 81

binom.c, fisher.c, hwe_test.cc, and all files under include/ except
plink2_highprec.* are covered by LGPLv3 (see the COPYING.LESSER file).
plink2_highprec.*'s adaptation of a small subset of the QD library
(https://github.com/BL-highprecision/QD ) is under the include/LICENSE.QD
BSD-3-Clause-LBNL license.  mini-gmp/* is dual-licensed under LGPLv3 and GPLv2.
