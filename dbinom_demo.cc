#include "include/binom.h"
#include "include/plink2_base.h"

#include <math.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  using namespace plink2;
  PglErr reterr = kPglRetSuccess;
  {
    if (argc != 4) {
      fputs(
"Binomial log-probability                      https://github.com/chrchang/stats\n"
"(C) 2013-2026 Christopher Chang     GNU Lesser General Public License version 3\n\n"
"  dbinom_demo <succ ct> <total obs ct> <expected succ rate>\n", stdout);
      reterr = kPglRetSkipped;
    } else {
      char* endptr;
      const double k = strtod(argv[1], &endptr);
      if (unlikely((endptr == argv[1]) || (k < 0) || (k != nearbyint(k)))) {
        fprintf(stderr, "Error: Invalid success count '%s'.\n", argv[1]);
        goto main_ret_INVALID_CMDLINE;
      }
      const double n = strtod(argv[2], &endptr);
      if (unlikely((endptr == argv[2]) || (n < 0) || (n != nearbyint(n)))) {
        fprintf(stderr, "Error: Invalid observation count '%s'.\n", argv[2]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(k > n)) {
        fputs("Error: # successes > # observations.\n", stderr);
        goto main_ret_INVALID_CMDLINE;
      }
      const double p = strtod(argv[3], &endptr);
      if ((endptr == argv[3]) || (p <= 0) || (p >= 1)) {
        fprintf(stderr, "Error: Invalid or unsupported rate '%s'.\n", argv[3]);
        goto main_ret_INVALID_CMDLINE;
      }
      const double ln_pval = BinomMassExtrange(k, n, p, 1);
      printf("Log-probability: %.17g\n", ln_pval);
    }
  }
  while (0) {
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
  return S_CAST(int32_t, reterr);
}
