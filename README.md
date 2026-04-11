stats
=====

Standalone functions for calculating binomial, Hardy-Weinberg equilibrium, and
Fisher 2x2/2x3 exact test statistics.

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
- binom.c, binom_demo.cc, fisher_demo.cc, hwe_demo.cc, and all of include/
  except the aforementioned bit of plink2_highprec are covered by LGPLv3 (see
  the COPYING.LESSER file).
- mini-gmp/ is dual-licensed under LGPLv3 and GPLv2.
