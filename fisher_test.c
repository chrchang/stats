#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

// Fisher 2x2 and 2x3 exact test command line utility
// Copyright (C) 2013  Christopher Chang  chrchang@alumni.caltech.edu

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22);

double fisher22_1sided(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t m11_is_greater_alt);

void fisher22_precomp_thresh(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiep);

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23);

#define MAXLINELEN 131072

int main(int argc, char** argv) {
  FILE* test_file;
  char buf[MAXLINELEN];
  char idstr[MAXLINELEN];
  char* bufptr;
  uint32_t m11;
  uint32_t m12;
  uint32_t m21;
  uint32_t m22;
  uint32_t m31;
  uint32_t m32;
  uint32_t tie;
  if (argc == 5) {
    m11 = atoi(argv[1]);
    m12 = atoi(argv[2]);
    m21 = atoi(argv[3]);
    m22 = atoi(argv[4]);
    printf("p-value: %g\n", fisher22(m11, m12, m21, m22));
    fisher22_precomp_thresh(m11, m12, m21, m22, &m31, &m32, &tie);
    if (!m32) {
      printf("(This is maximal.)\n");
    } else {
      printf("%u <= m11 < %u results in a higher p-value.\n", m31, m32);
    }
    if (tie != m11) {
      printf("m11 == %u results in the same p-value.\n", tie);
    }
    return 0;
  } else if (argc == 6) {
    if ((argv[5][1] != '\0') || ((argv[5][0] != '+') && (argv[5][0] != '-'))) {
      goto main_std_help;
    }
    m11 = atoi(argv[1]);
    m12 = atoi(argv[2]);
    m21 = atoi(argv[3]);
    m22 = atoi(argv[4]);
    printf("p-value: %g\n", fisher22_1sided(m11, m12, m21, m22, (argv[5][0] == '+')? 1 : 0));
    return 0;
  } else if (argc == 7) {
    printf("p-value: %g\n", fisher23(atoi(argv[1]), atoi(argv[3]), atoi(argv[5]), atoi(argv[2]), atoi(argv[4]), atoi(argv[6])));
    return 0;
  } else if (argc != 2) {
  main_std_help:
    printf(
"Fisher 2x2 and 2x3 exact test    https://github.com/chrchang/stats\n"
"(C) 2013 Christopher Chang, GNU General Public License version 3\n\n"
"Usage: fisher_test [m11] [m12] [m21] [m22] <+ | ->\n"
"       fisher_test [m11] [m12] [m21] [m22] [m31] [m32]\n"
"       fisher_test [filename]\n\n"
"For the 2x2 case, if the optional last parameter is '+', a 1-sided test is used\n"
"where the alternative hypothesis is that m11 is greater than expected;\n"
"similarly, '-' invokes the 1-sided test with the m11-is-less-than-expected\n"
"alternative.  With neither, a 2-sided test is employed.\n\n"
"If a filename is provided, each line of the file is expected to contain an ID\n"
"in the first column, and then either 4 or 6 values (in m11-m12-m21-m22-m31-m32\n"
"order).\n"
	   );
    return 1;
  }
  test_file = fopen(argv[1], "r");
  if (!test_file) {
    printf("Error: Unable to open file.\n");
    return 2;
  }
  buf[MAXLINELEN - 1] = ' ';
  while (fgets(buf, MAXLINELEN, test_file)) {
    if (!buf[MAXLINELEN - 1]) {
      printf("Error: Excessively long line in input file.\n");
      fclose(test_file);
      return 3;
    }
    bufptr = buf;
    while ((*bufptr == ' ') || (*bufptr == '\t')) {
      bufptr++;
    }
    if (*bufptr < ' ') {
      continue;
    }
    if (sscanf(bufptr, "%s %u %u %u %u %u %u", idstr, &m11, &m12, &m21, &m22, &m31, &m32) < 7) {
      if (sscanf(bufptr, "%s %u %u %u %u", idstr, &m11, &m12, &m21, &m22) < 5) {
        // skip improperly formatted line
        continue;
      }
      printf("p-value for %s: %g\n", idstr, fisher22(m11, m12, m21, m22));
    } else {
      printf("p-value for %s: %g\n", idstr, fisher23(m11, m21, m31, m12, m22, m32));
    }
  }
  if (!feof(test_file)) {
    printf("Error: File read failure.\n");
    fclose(test_file);
    return 4;
  }
  fclose(test_file);
  return 0;
}
