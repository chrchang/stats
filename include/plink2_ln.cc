// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_ln.h"

#include <math.h>
#include <string.h>

#include "plink2_float.h"

#ifdef __cplusplus
namespace plink2 {
#endif

CXXCONST_CP ScanadvLn(const char* str_iter, double* ln_ptr) {
  // revised ScanadvDouble() which currently requires number to be nonnegative
  // returns -DBL_MAX on 0
  uint32_t cur_char_code = ctou32(*str_iter);
  const uint32_t is_negative = (cur_char_code == 45);
  if (is_negative || (cur_char_code == 43)) {
    cur_char_code = ctou32(*(++str_iter));
  }
  uint32_t cur_digit = cur_char_code - 48;
  intptr_t e10 = 0;
  const char* dot_ptr;
  int64_t digits;
#ifdef __LP64__
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvLn_parse_decimal;
        }
        goto ScanadvLn_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    // (could keep ~19 instead, but if we're systematically losing the last two
    // bits of precision anyway...)
    const char* last_sig_fig_ptr = str_iter;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
    } while (cur_digit < 10);
    e10 = S_CAST(intptr_t, str_iter - last_sig_fig_ptr) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
    }
    goto ScanadvLn_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = str_iter;
  cur_digit = ctou32(*(++str_iter)) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits = cur_digit;
 ScanadvLn_parse_decimal:
  while (1) {
    cur_digit = ctou32(*(++str_iter)) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
      break;
    }
    digits = digits * 10 + cur_digit;
    if (digits >= 10000000000000000LL) {
      e10 = -S_CAST(intptr_t, str_iter - dot_ptr);
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
      break;
    }
  }
 ScanadvLn_parse_exponent:
  if (is_negative && (digits != 0)) {
    return nullptr;
  }
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = ctou32(*(++str_iter));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = ctou32(*(++str_iter));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 214748364) {
        // may as well guard against exponent overflow
        if (!exp_is_negative) {
          return nullptr;
        }
        *ln_ptr = -DBL_MAX;
        do {
          cur_digit = ctou32(*(++str_iter)) - 48;
        } while (cur_digit < 10);
        return S_CAST(CXXCONST_CP, str_iter);
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ctou32(*(++str_iter)) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#else  // not __LP64__
  int32_t digits_short;
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits_short = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvLn_parse_decimal;
        }
        digits = digits_short;
        goto ScanadvLn_parse_exponent;
      }
      digits_short = digits_short * 10 + cur_digit;
    } while (digits_short < 100000000);
    digits = digits_short;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvLn_parse_decimal_long;
        }
        goto ScanadvLn_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    const char* last_sig_fig_ptr = str_iter;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
    } while (cur_digit < 10);
    e10 = S_CAST(intptr_t, str_iter - last_sig_fig_ptr) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
    }
    goto ScanadvLn_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = str_iter;
  cur_digit = ctou32(*(++str_iter)) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits_short = cur_digit;
 ScanadvLn_parse_decimal:
  while (1) {
    cur_digit = ctou32(*(++str_iter)) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
      digits = digits_short;
      break;
    }
    digits_short = digits_short * 10 + cur_digit;
    if (digits_short >= 100000000) {
      digits = digits_short;
    ScanadvLn_parse_decimal_long:
      while (1) {
        cur_digit = ctou32(*(++str_iter)) - 48;
        if (cur_digit >= 10) {
          e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
          goto ScanadvLn_parse_exponent;
        }
        digits = digits * 10 + cur_digit;
        if (digits >= 10000000000000000LL) {
          e10 = -S_CAST(intptr_t, str_iter - dot_ptr);
          do {
            cur_digit = ctou32(*(++str_iter)) - 48;
          } while (cur_digit < 10);
          goto ScanadvLn_parse_exponent;
        }
      }
    }
  }
 ScanadvLn_parse_exponent:
  if (is_negative && (digits != 0)) {
    return nullptr;
  }
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = ctou32(*(++str_iter));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = ctou32(*(++str_iter));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 107374182) {
        // may as well guard against exponent overflow
        if (!exp_is_negative) {
          return nullptr;
        }
        *ln_ptr = -DBL_MAX;
        do {
          cur_digit = ctou32(*(++str_iter)) - 48;
        } while (cur_digit < 10);
        return S_CAST(CXXCONST_CP, str_iter);
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ctou32(*(++str_iter)) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#endif
  if (digits == 0) {
    *ln_ptr = -DBL_MAX;
    return S_CAST(CXXCONST_CP, str_iter);
  }
  double ln_val = log(S_CAST(double, digits));
  if (e10) {
    // I don't expect log() to be bit-identical between FMA and non-FMA, so
    // this doesn't make floating-point variation meaningfully worse...
    ln_val += e10 * kLn10;
  }
  *ln_ptr = ln_val;
  return S_CAST(CXXCONST_CP, str_iter);
}

static inline char* uitoa_trunc6(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  memcpy_k2(start, &(kDigitPair[quotient]));
  uii -= 10000 * quotient;
  if (uii) {
    quotient = uii / 100;
    start += 2;
    memcpy_k2(start, &(kDigitPair[quotient]));
    uii -= 100 * quotient;
    if (uii) {
      start += 2;
      memcpy_k2(start, &(kDigitPair[uii]));
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* rtoa_p5(uint32_t remainder, char* start) {
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  uint32_t quotient = remainder / 1000;
  memcpy_k2(start, &(kDigitPair[quotient]));
  remainder -= 1000 * quotient;
  if (remainder) {
    quotient = remainder / 10;
    start += 2;
    memcpy_k2(start, &(kDigitPair[quotient]));
    remainder -= 10 * quotient;
    if (remainder) {
      start[2] = '0' + remainder;
      return &(start[3]);
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* qrtoa_1p5(uint32_t quotient, uint32_t remainder, char* start) {
  *start++ = '0' + quotient;
  return rtoa_p5(remainder, start);
}

static const STD_ARRAY_DECL(double, 2, kBankerRound8) = STD_ARRAY_INIT_START() {0.499999995, 0.500000005} STD_ARRAY_INIT_END();

static inline uint32_t BankerRoundD(double dxx, STD_ARRAY_KREF(double, 2) banker_round) {
  uint32_t result = S_CAST(int32_t, dxx);
  return result + S_CAST(int32_t, (dxx - u31tod(result)) + banker_round[result & 1]);
}

// These are separate functions so the compiler can optimize the integer
// divisions.
static inline void BankerRoundD1(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10;
  *remainderp = remainder - (*quotientp) * 10;
}

static inline void BankerRoundD2(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100;
  *remainderp = remainder - (*quotientp) * 100;
}

static inline void BankerRoundD3(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 1000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 1000;
  *remainderp = remainder - (*quotientp) * 1000;
}

static inline void BankerRoundD4(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10000;
  *remainderp = remainder - (*quotientp) * 10000;
}

static inline void BankerRoundD5(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100000;
  *remainderp = remainder - (*quotientp) * 100000;
}

char* dtoa_so6(double dxx, char* start) {
  // 6 sig fig number, 0.999995 <= dxx < 999999.5
  // 'so' = "significand only"
  // Just hardcoding all six cases, in the absence of a better approach...
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.999949999999) {
    if (dxx < 9.9999949999999) {
      BankerRoundD5(dxx, kBankerRound8, &quotient, &remainder);
      return qrtoa_1p5(quotient, remainder, start);
    }
    BankerRoundD4(dxx, kBankerRound8, &quotient, &remainder);
    start = memcpya_k2(start, &(kDigitPair[quotient]));
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 100;
    memcpy_k2(start, &(kDigitPair[quotient]));
    remainder -= 100 * quotient;
    if (remainder) {
      start += 2;
    dtoa_so6_pretail:
      memcpy_k2(start, &(kDigitPair[remainder]));
    }
  dtoa_so6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  }
  if (dxx < 9999.9949999999) {
    if (dxx < 999.99949999999) {
      BankerRoundD3(dxx, kBankerRound8, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya_k2(start, &(kDigitPair[quotient]));
      if (!remainder) {
        return start;
      }
      *start++ = '.';
      quotient = remainder / 10;
      memcpy_k2(start, &(kDigitPair[quotient]));
      remainder -= quotient * 10;
      if (!remainder) {
        goto dtoa_so6_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    BankerRoundD2(dxx, kBankerRound8, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya_k2(start, &(kDigitPair[quotient]));
    quotient = uii - (100 * quotient);
    start = memcpya_k2(start, &(kDigitPair[quotient]));
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so6_pretail;
  }
  if (dxx >= 99999.949999999) {
    return u32toa_z6(BankerRoundD(dxx, kBankerRound8), start);
  }
  BankerRoundD1(dxx, kBankerRound8, &uii, &remainder);
  quotient = uii / 10000;
  *start = '0' + quotient;
  uii -= 10000 * quotient;
  quotient = uii / 100;
  start = memcpya_k2(&(start[1]), &(kDigitPair[quotient]));
  uii = uii - 100 * quotient;
  start = memcpya_k2(start, &(kDigitPair[uii]));
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  *start = '0' + remainder;
  return &(start[1]);
}

// todo: benchmark the exponential-notation part of this vs. dtoa_g(); maybe
// dtoa_g() should actually call this (or at least the exponential-notation
// part, put into its own function) for the most extreme values?
char* lntoa_g(double ln_val, char* start) {
  // log(999999.49999999)
  if (ln_val < 13.81551005796414) {
    // log(9.9999949999999e-5)
    if (ln_val > -9.210340871976317) {
      // No exponential notation.

      // log(0.99999949999999)
      if (ln_val > -5.000001349509205e-7) {
        // may as well fast-path x=1; since the most common use-case for this
        // function is p-value printing, x=1 should happen a lot more than x>1.
        // log(1.0000050000001)
        if (ln_val < 4.999987599993995e-6) {
          *start++ = '1';
          return start;
        }
        return dtoa_so6(exp(ln_val), start);
      }
      double dxx = exp(ln_val);
      // 6 sig fig decimal, no less than ~0.0001
      start = memcpya_k2(start, "0.");
      if (dxx < 9.9999949999999e-3) {
        dxx *= 100;
        start = memcpya_k2(start, "00");
      }
      if (dxx < 9.9999949999999e-2) {
        dxx *= 10;
        *start++ = '0';
      }
      return uitoa_trunc6(BankerRoundD(dxx * 1000000, kBankerRound8), start);
    }
    // if exponent is in danger of overflowing int32, just print '0'
    if (ln_val < 0x7ffffffb * (-kLn10)) {
      *start++ = '0';
      return start;
    }
  } else {
    // if exponent is in danger of overflowing int32, just print 'inf'
    if (ln_val > 0x7ffffffb * kLn10) {
      memcpy(start, "inf", 4);
      return &(start[3]);
    }
  }
  int32_t xp10 = S_CAST(int32_t, prefer_fma(ln_val, kRecipLn10, 5.000001349509205e-7 * kRecipLn10));
  double mantissa = exp(prefer_fma(xp10, -kLn10, ln_val));
  // mantissa will usually be in [.9999995, 9.999995], but |ln_val| can be
  // larger than 2^32, and floating point errors in either direction are
  // definitely possible (<20 bits of precision).
  if (mantissa < 0.99999949999999) {
    mantissa *= 10;
    xp10 -= 1;
  } else if (mantissa > 9.9999949999999) {
    mantissa *= 0.1;
    xp10 += 1;
  }
  uint32_t quotient;
  uint32_t remainder;
  BankerRoundD5(mantissa, kBankerRound8, &quotient, &remainder);
  start = qrtoa_1p5(quotient, remainder, start);
  if (xp10 < 0) {
    start = memcpya_k2(start, "e-");
    if (xp10 > -10) {
      *start++ = '0';
    }
    return u32toa(-xp10, start);
  }
  start = memcpya_k2(start, "e+");
  if (xp10 < 10) {
    *start++ = '0';
  }
  return u32toa(xp10, start);
}

#ifdef __cplusplus
}  // namespace plink2
#endif
