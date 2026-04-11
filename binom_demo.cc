// #include "include/binom.h"
#include "include/plink2_base.h"
#include "include/plink2_ln.h"

#include <errno.h>
#include <string.h>

#if defined(__cplusplus)
extern "C" {
#endif

double binom_2sided(uint32_t succ, uint32_t obs, double rate, uint32_t midp);
double binom_1sided(int32_t succ, int32_t obs, double rate);

#if defined(__cplusplus)
}
#endif

int main(int argc, char** argv) {
  using namespace plink2;
  FILE* test_file = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    uint32_t midp = 0;
    if (!strcmp(argv[argc - 1], "midp")) {
      midp = 1;
      argc--;
    }
    if ((argc == 4) || (argc == 5)) {
      uint32_t succ;
      uint32_t obs;
      double rate;
      if (unlikely(ScanUintDefcap(argv[1], &succ))) {
        fprintf(stderr, "Error: Invalid success count '%s'.\n", argv[1]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScanPosintDefcap(argv[2], &obs))) {
        fprintf(stderr, "Error: Invalid observation count '%s'.\n", argv[2]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely((sscanf(argv[3], "%lg", &rate) != 1) || (rate <= 0) || (rate >= 1))) {
        fprintf(stderr, "Error: Invalid expected success rate '%s'.\n", argv[3]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(succ > obs)) {
        fputs("Error: # successes > # observations.\n", stderr);
        goto main_ret_INVALID_CMDLINE;
      }
      double p_value;
      if (argc == 4) {
        p_value = binom_2sided(succ, obs, rate, midp);
      } else {
        if (argv[4][1] || ((argv[4][0] != '+') && (argv[4][0] != '-'))) {
          printf("Error: Invalid alternative hypothesis ('+' = more successes, '-' = fewer)\n");
          return 3;
        }
        if (midp) {
          printf("Error: 1-sided demo does not currently support mid-p adjustment.\n");
          return 3;
        }
        if (argv[4][0] == '+') {
          p_value = binom_1sided(obs - succ, obs, 1 - rate);
        } else {
          p_value = binom_1sided(succ, obs, rate);
        }
      }
      printf("P-value: %g\n", p_value);
    } else if (argc != 2) {
      fputs(
"Binomial test                                 https://github.com/chrchang/stats\n"
"(C) 2013-2026 Christopher Chang     GNU Lesser General Public License version 3\n\n"
"  binom_demo <success ct> <total obs ct> <expected succ rate> <+ | - | midp>\n"
"  binom_demo <filename>\n\n"
"If a filename is provided, each line of the file is expected to contain an ID\n"
"in the first column, and then 3 values (in succ-obs-rate order).\n", stdout);
      reterr = kPglRetSkipped;
    } else {
      test_file = fopen(argv[1], "r");
      if (!test_file) {
        fprintf(stderr, kErrprintfFopen, argv[1], strerror(errno));
        goto main_ret_OPEN_FAIL;
      }
      while (!feof(test_file)) {
        char name[100];
        int32_t succ;
        int32_t obs;
        double rate;
        fscanf(test_file, "%s %d %d %lg\n", name, &succ, &obs, &rate);
        const double p_value = binom_2sided(succ, obs, rate, midp);
        printf("P-value for %s: %g\n", name, p_value);
      }
    }
  }
  while (0) {
  main_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
  if (test_file) {
    fclose(test_file);
  }
  return S_CAST(int32_t, reterr);
}
