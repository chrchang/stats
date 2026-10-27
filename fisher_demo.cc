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
      // currently doesn't accept "-0", etc.
      return 1;
    }
    // accept leading '+'
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  *str_iterp = str_iter;
  return ScanmovU64CappedFinish(str_iterp, cap, valp);
}

static inline BoolErr ScantokU64Capped(const char* str_iter, uint64_t cap, uint64_t* valp) {
  return ScanmovU64Capped(&str_iter, cap, valp) || (!IsSpaceOrEoln(*str_iter));
}

static const uint32_t kFisher2xkDemoMaxK = 256;

// Parses up to max_ct whitespace-separated nonnegative integers (each at most
// cap), stopping at the end of the line or at the first token that isn't
// one.  Returns the number of values parsed.
static uint32_t ScanU64Tokens(const char* str_iter, uint32_t max_ct, uint64_t cap, uint64_t* vals) {
  uint32_t ct = 0;
  for (; ct != max_ct; ++ct) {
    while ((*str_iter == ' ') || (*str_iter == '\t')) {
      ++str_iter;
    }
    if (ScanmovU64Capped(&str_iter, cap, &(vals[ct])) || (!IsSpaceOrEoln(*str_iter))) {
      break;
    }
  }
  return ct;
}

// vals[] is in m11-m12-m21-m22-... order (k rows of 2), and must sum to less
// than 2^31.
static double Fisher2xkDemoLnP(const uint64_t* vals, uint32_t k, uint32_t midp) {
  int32_t obs[2 * kFisher2xkDemoMaxK];
  for (uint32_t row_idx = 0; row_idx != k; ++row_idx) {
    obs[row_idx] = S_CAST(int32_t, vals[2 * row_idx]);
    obs[k + row_idx] = S_CAST(int32_t, vals[2 * row_idx + 1]);
  }
  return Fisher2xkLnP(obs, k, midp, -DBL_MAX);
}

// vals[] is a 3x3 table in row-major order, and must sum to less than 2^31.
static double Fisher33DemoLnP(const uint64_t* vals, uint32_t midp) {
  int32_t obs[9];
  for (uint32_t cell_idx = 0; cell_idx != 9; ++cell_idx) {
    obs[cell_idx] = S_CAST(int32_t, vals[cell_idx]);
  }
  return Fisher33LnP(obs, midp, -DBL_MAX);
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
    } else if (argc == 10) {
      uint64_t vals[9];
      uint64_t tot = 0;
      for (uint32_t val_idx = 0; val_idx != 9; ++val_idx) {
        uint32_t cur_val;
        if (unlikely(ScanUintIcap(argv[val_idx + 1], &cur_val))) {
          fprintf(stderr, "Error: Invalid m%u%u value '%s'.\n", val_idx / 3 + 1, (val_idx % 3) + 1, argv[val_idx + 1]);
          goto main_ret_INVALID_CMDLINE;
        }
        vals[val_idx] = cur_val;
        tot += cur_val;
      }
      if (unlikely(tot > 0x7fffffff)) {
        fputs("Error: Problem instance too large.\n", stderr);
        reterr = kPglRetNotYetSupported;
        goto main_ret_1;
      }
      const double ln_pval = Fisher33DemoLnP(vals, midp);
      if (unlikely(ln_pval != ln_pval)) {
        goto main_ret_NOMEM;
      }
      char buf[80];
      char* write_iter = strcpya(buf, "P-value: ");
      write_iter = lntoa_g(ln_pval, write_iter);
      memcpy_k2(write_iter, "\n");
      fputs(buf, stdout);
    } else if ((argc >= 9) && (argc % 2) && (argc <= S_CAST(int32_t, 2 * kFisher2xkDemoMaxK + 1))) {
      const uint32_t k = (argc - 1) / 2;
      uint64_t vals[2 * kFisher2xkDemoMaxK];
      uint64_t tot = 0;
      for (uint32_t val_idx = 0; val_idx != 2 * k; ++val_idx) {
        uint32_t cur_val;
        if (unlikely(ScanUintIcap(argv[val_idx + 1], &cur_val))) {
          fprintf(stderr, "Error: Invalid m%u%u value '%s'.\n", val_idx / 2 + 1, (val_idx % 2) + 1, argv[val_idx + 1]);
          goto main_ret_INVALID_CMDLINE;
        }
        vals[val_idx] = cur_val;
        tot += cur_val;
      }
      if (unlikely(tot > 0x7fffffff)) {
        fputs("Error: Problem instance too large.\n", stderr);
        reterr = kPglRetNotYetSupported;
        goto main_ret_1;
      }
      const double ln_pval = Fisher2xkDemoLnP(vals, k, midp);
      if (unlikely(ln_pval != ln_pval)) {
        goto main_ret_NOMEM;
      }
      char buf[80];
      char* write_iter = strcpya(buf, "P-value: ");
      write_iter = lntoa_g(ln_pval, write_iter);
      memcpy_k2(write_iter, "\n");
      fputs(buf, stdout);
    } else if (argc != 2) {
    main_std_help:
      fputs(
"Fisher 2x2, 2xk, and 3x3 exact test           https://github.com/chrchang/stats\n"
"(C) 2013-2026 Christopher Chang            GNU General Public License version 3\n\n"
"  fisher_demo <m11> <m12> <m21> <m22> ['+' | '-'] ['midp']\n"
"  fisher_demo <m11> <m12> <m21> <m22> <m31> <m32> [<m41> <m42> ...] ['midp']\n"
"  fisher_demo <m11> <m12> <m13> <m21> <m22> <m23> <m31> <m32> <m33> ['midp']\n"
"  fisher_demo <filename> ['midp']\n\n"
"For the 2x2 case, if the optional 5th parameter is '+', a 1-sided test is\n"
"performed where the alternative hypothesis is that m11 is greater than\n"
"expected; similarly, '-' invokes the 1-sided test with the m11-less-than-exp.\n"
"alternative.  With neither, the 2-sided test is performed.\n\n"
"With 6 or an even number of 8 or more values, the table has one row per pair of\n"
"values (e.g. m31 and m32 are the third row); up to 256 rows are supported.\n"
"With exactly 9 values, the table is 3x3.\n\n"
"If a filename is provided, each line of the file is expected to contain an ID\n"
"in the first column, and then 4 or more values (in m11-m12-m21-m22-m31-m32-...\n"
"order), or 9 values for a 3x3 table.  On lines with any other odd number of\n"
"values, only the first 4 or 6 are used.\n\n"
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
        double ln_pval;
        {
          // Check for 3x3 and 4+-row tables first.
          const char* vals_start = bufptr;
          while (!IsSpaceOrEoln(*vals_start)) {
            ++vals_start;
          }
          uint64_t vals[2 * kFisher2xkDemoMaxK + 1];
          const uint32_t val_ct = ScanU64Tokens(vals_start, 2 * kFisher2xkDemoMaxK + 1, kFisher22Max, vals);
          if (val_ct == 9) {
            uint64_t tot = 0;
            for (uint32_t val_idx = 0; val_idx != 9; ++val_idx) {
              tot += vals[val_idx];
            }
            if (unlikely(tot > 0x7fffffff)) {
              fputs("Error: Problem instance too large.\n", stderr);
              reterr = kPglRetNotYetSupported;
              goto main_ret_1;
            }
            const uintptr_t id_slen = vals_start - bufptr;
            memcpy(idstr, bufptr, id_slen);
            idstr[id_slen] = '\0';
            ln_pval = Fisher33DemoLnP(vals, midp);
            if (unlikely(ln_pval != ln_pval)) {
              goto main_ret_NOMEM;
            }
            goto main_print_pval;
          }
          if ((val_ct >= 8) && (!(val_ct % 2))) {
            if (unlikely(val_ct > 2 * kFisher2xkDemoMaxK)) {
              fprintf(stderr, "Error: Line %" PRIuPTR " of %s has too many values.\n", line_idx, argv[1]);
              goto main_ret_MALFORMED_INPUT;
            }
            uint64_t tot = 0;
            for (uint32_t val_idx = 0; val_idx != val_ct; ++val_idx) {
              tot += vals[val_idx];
            }
            if (unlikely(tot > 0x7fffffff)) {
              fputs("Error: Problem instance too large.\n", stderr);
              reterr = kPglRetNotYetSupported;
              goto main_ret_1;
            }
            const uintptr_t id_slen = vals_start - bufptr;
            memcpy(idstr, bufptr, id_slen);
            idstr[id_slen] = '\0';
            ln_pval = Fisher2xkDemoLnP(vals, val_ct / 2, midp);
            if (unlikely(ln_pval != ln_pval)) {
              goto main_ret_NOMEM;
            }
            goto main_print_pval;
          }
        }
        {
          uint64_t m11;
          uint64_t m12;
          uint64_t m21;
          uint64_t m22;
          uint32_t m31;
          uint32_t m32;
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
        }
      main_print_pval:
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
  main_ret_NOMEM:
    fputs("Error: Out of memory.\n", stderr);
    reterr = kPglRetNomem;
    break;
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
