#include "include/fisher.h"
#include "include/plink2_base.h"
#include "include/plink2_ln.h"

#include <errno.h>
#include <string.h>

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(__cplusplus)
}
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
    if ((argc >= 5) && (argc <= 7)) {
      uint32_t m11;
      uint32_t m12;
      uint32_t m21;
      uint32_t m22;
      if (unlikely(ScanUintDefcap(argv[1], &m11))) {
        fprintf(stderr, "Error: Invalid m11 value '%s'.\n", argv[1]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScanUintDefcap(argv[2], &m12))) {
        fprintf(stderr, "Error: Invalid m12 value '%s'.\n", argv[2]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScanUintDefcap(argv[3], &m21))) {
        fprintf(stderr, "Error: Invalid m21 value '%s'.\n", argv[3]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScanUintDefcap(argv[4], &m22))) {
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
        if (unlikely(S_CAST(uint64_t, m11) + m12 + m21 + m22 + m31 + m32 > 0x7fffffff)) {
          fputs("Error: Problem instance too large.\n", stderr);
          reterr = kPglRetNotYetSupported;
          goto main_ret_1;
        }
        double logp;
        if (unlikely(Fisher23LnP(m11, m21, m31, m12, m22, m32, midp, &logp))) {
          goto main_ret_NOMEM;
        }
        char buf[80];
        char* write_iter = strcpya(buf, "P-value: ");
        write_iter = lntoa_g(logp, write_iter);
        memcpy_k2(write_iter, "\n");
        fputs(buf, stdout);
      } else {
        if (unlikely(S_CAST(uint64_t, m11) + m12 + m21 + m22 > 0x7fffffff)) {
          fputs("Error: Problem instance too large.\n", stderr);
          reterr = kPglRetNotYetSupported;
          goto main_ret_1;
        }
        if (argc == 5) {
          double logp;
          if (unlikely(Fisher22LnP(m11, m12, m21, m22, midp, &logp))) {
            goto main_ret_NOMEM;
          }
          char buf[80];
          char* write_iter = strcpya(buf, "P-value: ");
          write_iter = lntoa_g(logp, write_iter);
          memcpy_k2(write_iter, "\n");
          fputs(buf, stdout);
        } else {
          if ((argv[5][1] != '\0') || ((argv[5][0] != '+') && (argv[5][0] != '-'))) {
            goto main_std_help;
          }
          const double logp = Fisher22OneSidedLnP(m11, m12, m21, m22, (argv[5][0] == '+')? 1 : 0, midp);
          char buf[80];
          char* write_iter = strcpya(buf, "P-value: ");
          write_iter = lntoa_g(logp, write_iter);
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
"For the 2x2 case, if the optional 5th parameter is '+', a 1-sided test is used\n"
"where the alternative hypothesis is that m11 is greater than expected;\n"
"similarly, '-' invokes the 1-sided test with the m11-is-less-than-expected\n"
"alternative.  With neither, a 2-sided test is employed.\n\n"
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
        uint32_t m11;
        uint32_t m12;
        uint32_t m21;
        uint32_t m22;
        uint32_t m31;
        uint32_t m32;
        double logp;
        if (sscanf(bufptr, "%s %u %u %u %u %u %u", idstr, &m11, &m12, &m21, &m22, &m31, &m32) < 7) {
          if (sscanf(bufptr, "%s %u %u %u %u", idstr, &m11, &m12, &m21, &m22) < 5) {
            // skip improperly formatted line
            continue;
          }
          if (unlikely(Fisher22LnP(m11, m12, m21, m22, midp, &logp))) {
            goto main_ret_NOMEM;
          }
        } else {
          if (unlikely(Fisher23LnP(m11, m21, m31, m12, m22, m32, midp, &logp))) {
            goto main_ret_NOMEM;
          }
        }
        char buf[80];
        fputs("P-value for ", stdout);
        fputs(idstr, stdout);
        char* write_iter = strcpya(buf, ": ");
        write_iter = lntoa_g(logp, write_iter);
        memcpy_k2(write_iter, "\n");
        fputs(buf, stdout);
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
