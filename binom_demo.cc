#include "include/binom.h"
#include "include/plink2_base.h"
#include "include/plink2_highprec.h"
#include "include/plink2_ln.h"

#include <errno.h>
#include <string.h>

#ifdef __cplusplus
namespace plink2 {
#endif

// This may move to plink2_base.
static inline BoolErr ScanmovU64Finish(const char** str_iterp, uint64_t* valp) {
  const char* str_iter = *str_iterp;
  // limit is 20 digits, we've already read one
  const char* str_limit = &(str_iter[20]);
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = ctou64(*str_iter++) - 48;
    if (cur_digit >= 10) {
      break;
    }
    const uint64_t cur_digit2 = ctou64(*str_iter++) - 48;
    if (str_iter == str_limit) {
      if (unlikely((cur_digit2 < 10) || ((val >= (UINT64_MAX / 10)) && ((val > (UINT64_MAX / 10)) || (cur_digit > (UINT64_MAX % 10)))))) {
        return 1;
      }
      val = val * 10 + cur_digit;
      break;
    }
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
  }
  *valp = val;
  *str_iterp = &(str_iter[-1]);
  return 0;
}

BoolErr ScanmovU64(const char** str_iterp, uint64_t* valp) {
  const char* str_iter = *str_iterp;
  *valp = ctou32(*str_iter++) - 48;
  if (*valp >= 10) {
    if (unlikely(*valp != 0xfffffffbU)) {
      return 1;
    }
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  while (!(*valp)) {
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely((*valp) >= 10)) {
      return 1;
    }
  }
  *str_iterp = str_iter;
  return ScanmovU64Finish(str_iterp, valp);
}

BoolErr ParseProbStr(const char* fracstr, qd_real* p_qdr_ptr) {
  // Returns error on parse failure, value isn't in (0, 1), or the value can't
  // be represented with in-range denom.
  const char* slash_ptr = strchr(fracstr, '/');
  if (slash_ptr) {
    const char* numer_iter = fracstr;
    const char* denom_iter = &(slash_ptr[1]);
    uint64_t numer;
    uint64_t denom;
    if (ScanmovU64(&numer_iter, &numer) || ScanmovU64(&denom_iter, &denom) || (numer_iter != slash_ptr) || (denom_iter[0] != '\0') || (numer == 0) || (numer >= denom)) {
      return 1;
    }
    *p_qdr_ptr = qdr_accurate_div(qdr_makeu64(numer), qdr_makeu64(denom));
    return 0;
  }
  double dxx;
  if (sscanf(fracstr, "%lg", &dxx) != 1) {
    return 1;
  }
  // ensure we fail here on nan
  if ((dxx <= 0.0) || (!(dxx < 1.0))) {
    return 1;
  }
  *p_qdr_ptr = qdr_make1(dxx);
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
      double ln_pval;
      if (argc == 4) {
        qd_real p_qdr;
        if (unlikely(ParseProbStr(argv[3], &p_qdr))) {
          fprintf(stderr, "Error: Invalid or unsupported rate '%s'.\n", argv[3]);
          goto main_ret_INVALID_CMDLINE;
        }
        ln_pval = BinomTwoSidedP(succ, obs, p_qdr, midp, 1);
        fputs("Two-sided ", stdout);
      } else {
        qd_real p_qdr;
        if (unlikely(ParseProbStr(argv[3], &p_qdr))) {
          fprintf(stderr, "Error: Invalid or unsupported rate '%s'.\n", argv[3]);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(argv[4][1] || ((argv[4][0] != '+') && (argv[4][0] != '-')))) {
          fputs("Error: Invalid alternative hypothesis ('+' = more successes, '-' = fewer)\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        const dd_real q_ddr = ddr_negate(ddr_make_qd(qdr_addd(p_qdr, -1.0)));
        ln_pval = BinomOneSidedP(succ, obs, ddr_make_qd(p_qdr), q_ddr, (argv[4][0] == '+'), midp, 1);
        fputs("One-sided ", stdout);
      }
      if (midp) {
        fputs("mid", stdout);
      }
      fputs("p-value: ", stdout);
      char buf[80];
      char* write_iter = lntoa_g(ln_pval, buf);
      memcpy_k2(write_iter, "\n");
      fputs(buf, stdout);
    } else if (argc != 2) {
      fputs(
"Binomial test                                 https://github.com/chrchang/stats\n"
"(C) 2013-2026 Christopher Chang     GNU Lesser General Public License version 3\n\n"
"  binom_demo <succ ct> <total obs ct> <expected succ rate> ['+' | '-'] ['midp']\n"
"  binom_demo <filename> ['midp']\n\n"
"Rates can be entered as fractions (e.g. '3/10', no spaces allowed).  When a\n"
"rate is entered as a decimal, it is first parsed as a float64 (with the usual\n"
"potential for rounding error), and then converted to a fraction with a\n"
"power-of-2 denominator.\n\n"
"If the optional 4th parameter is '+', a 1-sided test is performed where the\n"
"alternative hypothesis is that #succ is greater than expected; similarly, '-'\n"
"invokes the 1-sided test with the #succ-less-than-expected alternative.\n"
"With neither, the 2-sided test is performed.\n\n"
"If a filename is provided, each line of the file is expected to contain an ID\n"
"in the first column, and then 3 values (in succ-obs-rate order; one-sided test)\n"
"not currently supported here).\n", stdout);
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
        qd_real p_qdr;
        if (unlikely(ParseProbStr(ratestr, &p_qdr))) {
          fprintf(stderr, "Error: Invalid or unsupported rate '%s' on line %" PRIuPTR " of %s.\n", ratestr, line_idx, argv[1]);
          goto main_ret_MALFORMED_INPUT;
        }
        const double ln_pval = BinomTwoSidedP(succ, obs, p_qdr, midp, 1);
        fputs("Two-sided ", stdout);
        if (midp) {
          fputs("mid", stdout);
        }
        fputs("p-value for ", stdout);
        fputs(name, stdout);
        char buf[80];
        char* write_iter = memcpya_k2(buf, ": ");
        write_iter = lntoa_g(ln_pval, write_iter);
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
