#include <stdint.h>
#include <inttypes.h>

#include <stdio.h>
#include <stdlib.h>

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2);
int32_t SNPHWE_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);

int main(int argc, char** argv) {
  FILE* test_file;
  int32_t hets;
  int32_t homs1;
  int32_t homs2;
  double p_value;
  char name[100];
  if ((argc == 4) || (argc == 5)) {
    hets = atoi(argv[1]);
    homs1 = atoi(argv[2]);
    homs2 = atoi(argv[3]);
    if ((hets < 0) || (homs1 < 0) || (homs2 < 0)) {
      printf("Error: Negative parameter.\n");
      return 3;
    }
    if (argc == 4) {
      p_value = SNPHWE2(hets, homs1, homs2);
      printf("P-value: %lg\n", p_value);
    } else {
      if (sscanf(argv[4], "%lg", &p_value) != 1) {
	printf("Error: Invalid p-value '%s'.\n", argv[4]);
	return 3;
      }
      if (SNPHWE_t(hets, homs1, homs2, p_value)) {
	printf("Test failed (p-value below threshold).\n");
      } else {
	printf("Test passed (p-value above threshold).\n");
      }
    }
    return 0;
  } else if (argc != 2) {
    printf(
"SNPHWE2 demo    http://www.sph.umich.edu/csg/abecasis/Exact/\n\n"
"  snphwe_test [het count] [hom count 1] [hom count 2] {threshold}\n"
"  snphwe_test [filename]\n\n"
"If a filename is provided, the file is expected to contain marker names in the\n"
"first column, heterozygote counts in the second, and homozygote counts in the\n"
"third and fourth.\n");
    return 1;
  }
  test_file = fopen(argv[1], "r");
  if (!test_file) {
    printf("Error: Unable to open file.\n");
    return 2;
  }
  while (!feof(test_file)) {
    fscanf(test_file, "%s %d %d %d\n", name, &hets, &homs1, &homs2);
    p_value = SNPHWE2(hets, homs1, homs2);
    printf("P-value for %s: %lg\n", name, p_value);
  }
  fclose(test_file);
  return 0;
}
