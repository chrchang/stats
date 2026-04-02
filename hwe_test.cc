#include "include/plink2_base.h"
#include "include/plink2_hwe.h"
#include "include/plink2_ln.h"

#include <errno.h>
#include <string.h>

int32_t main(int argc, char** argv) {
  using namespace plink2;
  FILE* test_file = nullptr;
  PglErr reterr = kPglRetSuccess;
  uint32_t midp = 0;
  {
    uint32_t midp = 0;
    if (!strcmp(argv[argc - 1], "midp")) {
      midp = 1;
      argc--;
    }
    if ((argc == 4) || (argc == 5)) {
      uint32_t hets;
      uint32_t homs1;
      uint32_t homs2;
      if (unlikely(ScanUintDefcap(argv[1], &hets))) {
        fprintf(stderr, "Error: Invalid het count '%s'.\n", argv[1]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScanUintDefcap(argv[2], &homs1))) {
        fprintf(stderr, "Error: Invalid hom1 count '%s'.\n", argv[2]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely(ScanUintDefcap(argv[3], &homs2))) {
        fprintf(stderr, "Error: Invalid hom2 count '%s'.\n", argv[3]);
        goto main_ret_INVALID_CMDLINE;
      }
      if (argc == 4) {
        double logp;
        if (unlikely(HweLnP(hets, homs1, homs2, midp, &logp))) {
          goto main_ret_NOMEM;
        }
        char buf[80];
        char* write_iter = strcpya(buf, "P-value: ");
        write_iter = lntoa_g(logp, write_iter);
        memcpy_k2(write_iter, "\n");
        fputs(buf, stdout);
      } else {
        double ln_thresh;
        if (unlikely((!ScantokLn(argv[4], &ln_thresh)) || (ln_thresh > 0.0))) {
          fprintf(stderr, "Error: Invalid p-value threshold '%s'.\n", argv[4]);
          goto main_ret_INVALID_CMDLINE;
        }
        uint32_t passed;
        if (unlikely(HweThreshLn(hets, homs1, homs2, midp, exp(ln_thresh), ln_thresh, &passed))) {
          goto main_ret_NOMEM;
        }
        if (passed) {
          printf("Test failed (p-value below threshold).\n");
        } else {
          printf("Test passed (p-value above threshold).\n");
        }
      }
    } else if (argc != 2) {
      fputs(
"HweLnP() demo                                 https://github.com/chrchang/stats\n\n"
"  hwe_test <het count> <hom count 1> <hom count 2> [threshold] ['midp']\n"
"  hwe_test <filename> ['midp']\n\n"
"If a filename is provided, the file is expected to contain marker names in the\n"
"first column, heterozygote counts in the second, and homozygote counts in the\n"
"third and fourth.\n", stdout);
      reterr = kPglRetSkipped;
    } else {
      test_file = fopen(argv[1], "r");
      if (!test_file) {
        fprintf(stderr, kErrprintfFopen, argv[1], strerror(errno));
        goto main_ret_OPEN_FAIL;
      }
      uintptr_t line_idx = 0;
      while (!feof(test_file)) {
        char name[kMaxMediumLine];
        int32_t hets;
        int32_t homs1;
        int32_t homs2;
        ++line_idx;
        if (unlikely(fscanf(test_file, "%s %d %d %d\n", name, &hets, &homs1, &homs2) != 4)) {
          fprintf(stderr, "Error: Failed to read or parse line %" PRIuPTR " of %s .\n", line_idx, argv[1]);
          goto main_ret_READ_FAIL;
        }
        double logp;
        if (unlikely(HweLnP(hets, homs1, homs2, midp, &logp))) {
          goto main_ret_NOMEM;
        }
        fputs("P-value for ", stdout);
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
  main_ret_READ_FAIL:
    reterr = kPglRetReadFail;
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
