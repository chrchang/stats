#include "include/fisher.h"
#include "include/hypergeom.h"
#include "include/plink2_base.h"
#include "include/plink2_ln.h"

#include <errno.h>
#include <string.h>

#ifdef __cplusplus
namespace plink2 {
#endif

// Assumes cap < 2^64 / 100.
// This may move to plink2_base.
static inline BoolErr ScanmovU64CappedFinish(const char** str_iterp, uint64_t cap, uint64_t* valp) {
  const char* str_iter = *str_iterp;
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = ctou64(*str_iter++) - 48;
    if (cur_digit >= 10) {
      break;
    }
    const uint64_t cur_digit2 = ctou64(*str_iter++) - 48;
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      if (unlikely(val > cap)) {
        return 1;
      }
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
    if (unlikely(val > cap)) {
      return 1;
    }
  }
  *valp = val;
  *str_iterp = &(str_iter[-1]);
  return 0;
}

BoolErr ScanmovU64Capped(const char** str_iterp, uint64_t cap, uint64_t* valp) {
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
  return ScanmovU64CappedFinish(str_iterp, cap, valp);
}

static inline BoolErr ScantokU64Capped(const char* str_iter, uint64_t cap, uint64_t* valp) {
  return ScanmovU64Capped(&str_iter, cap, valp) || (!IsSpaceOrEoln(*str_iter));
}

#ifdef __cplusplus
}  // namespace plink2
#endif

int32_t main(int argc, char** argv) {
  using namespace plink2;
  FILE* test_file = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    uint32_t midp = 0;
    if (!strcmp(argv[argc - 1], "midp")) {
      midp = 1;
      argc--;
    }
    const uint64_t kFisher22Max = (1LLU << 52) - 1;
    if ((argc >= 5) && (argc <= 7)) {
      uint64_t m11;
      uint64_t m12;
      uint64_t m21;
      uint64_t m22;
      if (unlikely(ScantokU64Capped(argv[1], kFisher22Max, &m11))) {
        fprintf(stderr, "Error: Invalid m11 value '%s'.\n", argv[1]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScantokU64Capped(argv[2], kFisher22Max, &m12))) {
        fprintf(stderr, "Error: Invalid m12 value '%s'.\n", argv[2]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScantokU64Capped(argv[3], kFisher22Max, &m21))) {
        fprintf(stderr, "Error: Invalid m21 value '%s'.\n", argv[3]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScantokU64Capped(argv[4], kFisher22Max, &m22))) {
        fprintf(stderr, "Error: Invalid m22 value '%s'.\n", argv[4]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (argc == 7) {
        uint32_t m31;
        uint32_t m32;
        if (unlikely(ScanUintDefcap(argv[5], &m31))) {
          fprintf(stderr, "Error: Invalid m11 value '%s'.\n", argv[5]);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(ScanUintDefcap(argv[6], &m32))) {
          fprintf(stderr, "Error: Invalid m12 value '%s'.\n", argv[6]);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(m11 + m12 + m21 + m22 + m31 + m32 > 0x7fffffff)) {
          fputs("Error: Problem instance too large.\n", stderr);
          reterr = kPglRetNotYetSupported;
          goto main_ret_1;
        }
        double ln_pval = Fisher23LnP(m11, m21, m31, m12, m22, m32, midp);
        char buf[80];
        char* write_iter = strcpya(buf, "P-value: ");
        write_iter = lntoa_g(ln_pval, write_iter);
        memcpy_k2(write_iter, "\n");
        fputs(buf, stdout);
      } else {
        if (unlikely(m11 + m12 + m21 + m22 > kFisher22Max)) {
          fputs("Error: Problem instance too large.\n", stderr);
          reterr = kPglRetNotYetSupported;
          goto main_ret_1;
        }
        if (argc == 5) {
          double ln_pval = Fisher22TwoSidedP(m11, m12, m21, m22, midp, 1);
          char buf[80];
          char* write_iter = strcpya(buf, "P-value: ");
          write_iter = lntoa_g(ln_pval, write_iter);
          memcpy_k2(write_iter, "\n");
          fputs(buf, stdout);
        } else {
          if ((argv[5][1] != '\0') || ((argv[5][0] != '+') && (argv[5][0] != '-'))) {
            goto main_std_help;
          }
          const double ln_pval = PhyperApprox(m11, m12, m21, m22, (argv[5][0] == '+')? 1 : 0, midp, 1);
          char buf[80];
          char* write_iter = strcpya(buf, "P-value: ");
          write_iter = lntoa_g(ln_pval, write_iter);
          memcpy_k2(write_iter, "\n");
          fputs(buf, stdout);
        }
      }
    } else if (argc != 2) {
    main_std_help:
      fputs(
"Fisher 2x2 and 2x3 exact test                 https://github.com/chrchang/stats\n"
"(C) 2013-2026 Christopher Chang     GNU Lesser General Public License version 3\n\n"
"  fisher_demo <m11> <m12> <m21> <m22> ['+' | '-'] ['midp']\n"
"  fisher_demo <m11> <m12> <m21> <m22> <m31> <m32> ['midp']\n"
"  fisher_demo <filename> ['midp']\n\n"
"For the 2x2 case, if the optional 5th parameter is '+', a 1-sided test is\n"
"performed where the alternative hypothesis is that m11 is greater than\n"
"expected; similarly, '-' invokes the 1-sided test with the m11-less-than-exp.\n"
"alternative.  With neither, the 2-sided test is performed.\n\n"
"If a filename is provided, each line of the file is expected to contain an ID\n"
"in the first column, and then either 4 or 6 values (in m11-m12-m21-m22-m31-m32\n"
"order).\n\n"
"If 'midp' is the last parameter, Lancaster's mid-p correction is applied.\n", stdout);
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
          bufptr++;
        }
        if (*bufptr < ' ') {
          continue;
        }
        char idstr[kMaxMediumLine];
        uint64_t m11;
        uint64_t m12;
        uint64_t m21;
        uint64_t m22;
        uint32_t m31;
        uint32_t m32;
        double ln_pval;
        if (sscanf(bufptr, "%s %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %u %u", idstr, &m11, &m12, &m21, &m22, &m31, &m32) < 7) {
          if (sscanf(bufptr, "%s %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64, idstr, &m11, &m12, &m21, &m22) < 5) {
            // skip improperly formatted line
            continue;
          }
          if (unlikely((m11 > kFisher22Max) || (m12 > kFisher22Max) || (m21 > kFisher22Max) || (m22 > kFisher22Max) || (m11 + m12 + m21 + m22 > kFisher22Max))) {
            fputs("Error: Problem instance too large.\n", stderr);
            reterr = kPglRetNotYetSupported;
            goto main_ret_1;
          }
          ln_pval = Fisher22TwoSidedP(m11, m12, m21, m22, midp, 1);
        } else {
          if (unlikely((m11 > 0x7fffffff) || (m12 > 0x7fffffff) || (m21 > 0x7fffffff) || (m22 > 0x7fffffff) || (m31 > 0x7fffffff) || (m32 > 0x7fffffff) || (m11 + m12 + m21 + m22 + m31 + m32 > 0x7fffffff))) {
            fputs("Error: Problem instance too large.\n", stderr);
            reterr = kPglRetNotYetSupported;
            goto main_ret_1;
          }
          ln_pval = Fisher23LnP(m11, m21, m31, m12, m22, m32, midp);
        }
        char buf2[80];
        fputs("P-value for ", stdout);
        fputs(idstr, stdout);
        char* write_iter = strcpya(buf2, ": ");
        write_iter = lntoa_g(ln_pval, write_iter);
        memcpy_k2(write_iter, "\n");
        fputs(buf2, stdout);
      }
      if (!feof(test_file)) {
        fputs("Error: File read failure.\n", stderr);
        goto main_ret_READ_FAIL;
      }
    }
  }
  while (0) {
  main_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  main_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  main_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 main_ret_1:
  if (test_file) {
    fclose(test_file);
  }
  return S_CAST(int32_t, reterr);
}
