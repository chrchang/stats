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
  if (argc == 5) {
    printf("p-value: %g\n", fisher22(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])));
    return 0;
  } else if (argc == 7) {
    printf("p-value: %g\n", fisher23(atoi(argv[1]), atoi(argv[3]), atoi(argv[5]), atoi(argv[2]), atoi(argv[4]), atoi(argv[6])));
    return 0;
  } else if (argc != 2) {
    printf(
"Fisher 2x2 and 2x3 exact test    https://www.cog-genomics.org/software/stats\n"
"(C) 2013 Christopher Chang, GNU General Public License version 3\n\n"
"Usage: fisher [m11] [m12] [m21] [m22]\n"
"       fisher [m11] [m12] [m21] [m22] [m31] [m32]\n"
"       fisher [filename]\n\n"
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
