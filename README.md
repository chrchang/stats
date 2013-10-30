stats
=====

Standalone functions for calculating binomial exact test, Hardy-Weinberg
equilibrium, and Fisher 2x2/2x3 exact test statistics.

Test program compilation:

gcc -O2 binom.c binom_test.c -o binom_test

gcc -O2 snphwe2.c snphwe_test.c -o snphwe_test

gcc -O2 fisher.c fisher_test.c -o fisher_test

Test program usage examples:

  binom_test binom_test.txt

  binom_test 2900 5000 0.52

  snphwe_test snphwe_test.txt

  snphwe_test 16 3 81

  fisher_test fisher_test.txt

  fisher_test 20 30 40 50 60 70

binom.c and fisher.c are covered by GPLv3 (see the LICENSE file).  For
snphwe2.c usage rules, refer to the Abecasis Lab's SNP-HWE page at
http://sph.umich.edu/csg/abecasis/Exact .