#ifndef __PLINK2_LN_H__
#define __PLINK2_LN_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Just the Scan{adv,tok}Ln() and lntoa_g() functions, and what they depend on,
// from plink2_string.

#include "plink2_base.h"

#ifdef _WIN32
#  define EOLN_STR "\r\n"
#else
#  define EOLN_STR "\n"
#endif

// generic maximum line byte length, currently also used as a default I/O
// buffer size.  .ped/.vcf/etc. lines can of course be longer.
CONSTI32(kMaxMediumLine, 131072);

#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef __cplusplus
#  define CXXCONST_CP const char*
#  define CXXCONST_VOIDP const void*
#  define TO_CONSTCPCONSTP(char_pp) (char_pp)

HEADER_INLINE const char* DowncastToXC(const void* pp) {
  return S_CAST(const char*, pp);
}
#else
#  define CXXCONST_CP char*
#  define CXXCONST_VOIDP void*
#  define TO_CONSTCPCONSTP(char_pp) ((const char* const*)(char_pp))

HEADER_INLINE char* DowncastToXC(const void* pp) {
  return R_CAST(char*, pp);
}
#endif

HEADER_INLINE int32_t IsSpaceOrEoln(unsigned char ucc) {
  return (ucc <= 32);
}

CXXCONST_CP ScanadvLn(const char* str_iter, double* ln_ptr);

HEADER_INLINE CXXCONST_CP ScantokLn(const char* str_iter, double* ln_ptr) {
  CXXCONST_CP parsed_end = ScanadvLn(str_iter, ln_ptr);
  if ((!parsed_end) || (!IsSpaceOrEoln(*parsed_end))) {
    return nullptr;
  }
  return parsed_end;
}

#ifdef __cplusplus
HEADER_INLINE char* ScanadvLn(char* str_iter, double* ln_ptr) {
  return const_cast<char*>(ScanadvLn(const_cast<const char*>(str_iter), ln_ptr));
}

HEADER_INLINE char* ScantokLn(char* str_iter, double* ln_ptr) {
  return const_cast<char*>(ScantokLn(const_cast<const char*>(str_iter), ln_ptr));
}
#endif

char* lntoa_g(double ln_val, char* start);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_LN_H__
