#include "include/binom.h"
#include "include/plink2_base.h"
#include "include/plink2_ln.h"

#include <errno.h>
#include <string.h>

#if defined(__cplusplus)
extern "C" {
#endif

double binom_1sided(int32_t succ, int32_t obs, double rate);

#if defined(__cplusplus)
}
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

BoolErr ParseProbFrac(const char* fracstr, int64_t* numerp, int64_t* denomp) {
  // Returns error on parse failure, value isn't in (0, 1), or the value can't
  // be represented with in-range denom.
  // todo: allow e.g. "1/10"
  double dxx;
  if (sscanf(fracstr, "%lg", &dxx) != 1) {
    return 1;
  }
  // ensure we fail here on nan
  if ((dxx <= 0.0) || (!(dxx < 1.0))) {
    return 1;
  }
  int pow;
  double value = frexp(dxx, &pow);
  int64_t numer = S_CAST(int64_t, scalbn(value, 53));
  int rshift = ctzu64(numer);
  pow = 53 - rshift - pow;
  if (pow > 62) {
    return 1;
  }
  assert(pow >= 0);
  *numerp = numer >> rshift;
  *denomp = 1LL << S_CAST(uint32_t, pow);
  return 0;
}

#ifdef __cplusplus
}  // namespace plink2
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
      if (unlikely(ScanUintDefcap(argv[1], &succ))) {
        fprintf(stderr, "Error: Invalid success count '%s'.\n", argv[1]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScanPosintDefcap(argv[2], &obs))) {
        fprintf(stderr, "Error: Invalid observation count '%s'.\n", argv[2]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(succ > obs)) {
        fputs("Error: # successes > # observations.\n", stderr);
        goto main_ret_INVALID_CMDLINE;
      }
      if (argc == 4) {
        int64_t rate_numer;
        int64_t rate_denom;
        if (unlikely(ParseProbFrac(argv[3], &rate_numer, &rate_denom))) {
          fprintf(stderr, "Error: Invalid rate '%s'.\n", argv[3]);
          goto main_ret_MALFORMED_INPUT;
        }
        double logp;
        if (unlikely(BinomLnP(succ, obs, rate_numer, rate_denom - rate_numer, midp, &logp))) {
          goto main_ret_NOMEM;
        }
        fputs("Two-sided ", stdout);
        if (midp) {
          fputs("mid", stdout);
        }
        fputs("p-value: ", stdout);
        char buf[80];
        char* write_iter = lntoa_g(logp, buf);
        memcpy_k2(write_iter, "\n");
        fputs(buf, stdout);
      } else {
        double rate;
        if (unlikely((sscanf(argv[3], "%lg", &rate) != 1) || (rate <= 0) || (rate >= 1))) {
          fprintf(stderr, "Error: Invalid expected success rate '%s'.\n", argv[3]);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(argv[4][1] || ((argv[4][0] != '+') && (argv[4][0] != '-')))) {
          fputs("Error: Invalid alternative hypothesis ('+' = more successes, '-' = fewer)\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(midp)) {
          fputs("Error: 1-sided demo does not currently support mid-p adjustment.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        double p_value;
        if (argv[4][0] == '+') {
          p_value = binom_1sided(obs - succ, obs, 1 - rate);
        } else {
          p_value = binom_1sided(succ, obs, rate);
        }
        printf("P-value: %g\n", p_value);
      }
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
      char buf[kMaxMediumLine];
      buf[kMaxMediumLine - 1] = ' ';
      uintptr_t line_idx = 0;
      while (fgets(buf, kMaxMediumLine, test_file)) {
        ++line_idx;
        if (!buf[kMaxMediumLine - 1]) {
          fprintf(stderr, "Error: Line %" PRIuPTR " of %s is too long.\n", line_idx, argv[1]);
          goto main_ret_MALFORMED_INPUT;
        }
        char* bufptr = buf;
        while ((*bufptr == ' ') || (*bufptr == '\t')) {
          ++bufptr;
        }
        if (*bufptr < ' ') {
          continue;
        }
        char name[kMaxMediumLine];
        char ratestr[kMaxMediumLine];
        int32_t succ;
        int32_t obs;
        if (unlikely(sscanf(bufptr, "%s %d %d %s", name, &succ, &obs, ratestr) < 4)) {
          fprintf(stderr, "Error: Line %" PRIuPTR " of %s has fewer fields than expected.\n", line_idx, argv[1]);
          goto main_ret_MALFORMED_INPUT;
        }
        int64_t rate_numer;
        int64_t rate_denom;
        if (unlikely(ParseProbFrac(ratestr, &rate_numer, &rate_denom))) {
          fprintf(stderr, "Error: Invalid rate '%s' on line %" PRIuPTR " of %s.\n", ratestr, line_idx, argv[1]);
          goto main_ret_MALFORMED_INPUT;
        }
        double logp;
        if (unlikely(BinomLnP(succ, obs, rate_numer, rate_denom - rate_numer, midp, &logp))) {
          goto main_ret_NOMEM;
        }
        fputs("Two-sided ", stdout);
        if (midp) {
          fputs("mid", stdout);
        }
        fputs("p-value for ", stdout);
        fputs(name, stdout);
        char buf[80];
        char* write_iter = memcpya_k2(buf, ": ");
        write_iter = lntoa_g(logp, write_iter);
        memcpy_k2(write_iter, "\n");
        fputs(buf, stdout);
      }
    }
  }
  while (0) {
  main_ret_NOMEM:
    fputs("Error: Out of memory.\n", stderr);
    reterr = kPglRetNomem;
    break;
  main_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  if (test_file) {
    fclose(test_file);
  }
  return S_CAST(int32_t, reterr);
}
