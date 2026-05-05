# cython: language_level=3
from libc.stdint cimport int64_t, uint32_t, int32_t
from libc.math cimport exp
import fractions

__version__ = "0.1.0"

cdef extern from "../include/binom.h" namespace "plink2":
    ctypedef uint32_t BoolErr

    BoolErr BinomLnP(int32_t obs_succ, int32_t obs_tot, int64_t succ_odds_ratio_numer, int64_t succ_odds_ratio_denom, int32_t midp, double* resultp) nogil

    double BinomOneSidedLnP(int32_t obs_succ, int32_t obs_tot, double succ_odds_ratio, uint32_t succ_is_greater_alt, int32_t midp)


cdef extern from "../include/fisher.h" namespace "plink2":
    BoolErr Fisher22LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, double* resultp) nogil

    double Fisher22OneSidedLnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp) nogil

    BoolErr Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp, double* resultp) nogil


# n must be in [1, 2^31 - 1], and k must be in [0, n].
#
# For the two-sided test, to enable precise classification of likelihood ties
# and near-ties, p is expected to be a fractions.Fraction in (2^{-63}, 1) with
# positive denominator < 2^63, or valid input to the fractions.Fraction()
# constructor (e.g. float, string, decimal.Decimal).  Note that float values
# smaller than 1/512 may convert to a fraction with too large of a denominator;
# a workaround is to pass in
#   fractions.Fraction(p).limit_denominator(2**63 - 1)
# when that's a problem.  (This conversion isn't done automatically because it
# would water down what this function promises about likelihood ties and
# near-ties.)
#
# For the one-sided tests, p is expected to be a float in (2^{-63}, 1), or
# valid input to the float() constructor.  (Likelihood-ties don't matter in
# this case.)
def binom(int32_t k, int32_t n, object p=0.5, str alternative = "two-sided", bint midp = 0, bint logp = 0):
    if k < 0 or k > n:
        raise RuntimeError("k must be nonnegative and <= n.")
    cdef double ln_result
    if alternative == "two-sided":
        if not isinstance(p, fractions.Fraction):
            p = fractions.Fraction(p)
        numer, denom = p.as_integer_ratio()
        if numer <= 0 or numer >= denom:
            raise RuntimeError("p must be in (2^{-63}, 1).")
        if denom >= 2**63:
            raise RuntimeError("for the two-sided test, p must be precisely representable as a fraction with denominator < 2^63; one workaround is to pass in fractions.Fraction(p).limit_denominator(2**63 - 1).")
        if BinomLnP(k, n, <int64_t>(numer), <int64_t>(denom - numer), midp, &ln_result):
            raise MemoryError()
        if logp:
            return ln_result
        else:
            return exp(ln_result)
    cdef bint k_is_greater_alt = (alternative == "greater")
    if alternative != "less" and not k_is_greater_alt:
        raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
    cdef double pfloat = float(p)
    # strictly speaking, BinomOneSidedLnP() still works with p as small as
    # 10^{-298}, but let's minimize the difference between the one-sided and
    # two-sided interface for now.
    if pfloat <= 1.0842021724855044e-19 or pfloat >= 1:
        raise RuntimeError("p must be in (2^{-63}, 1).")
    ln_result = BinomOneSidedLnP(k, n, pfloat / (1 - pfloat), k_is_greater_alt, midp)
    if logp:
        return ln_result
    else:
        return exp(ln_result)


# table must be a 2x2 or larger matrix, represented as a list of equal-length
# lists.  Values must be nonnegative integers which add up to <2^31.
#
# alternative must be one of the following:
#   "two-sided": default, must be this if table is larger than 2x2.
#   "less": alt hypothesis is that table[0][0] is smaller than expected.
#   "greater": alt hypothesis is that table[0][0] is larger than expected.
def fisher(list table, str alternative="two-sided", bint midp=0, bint logp=0):
    cdef uint32_t nrow = len(table)
    if nrow < 2:
        raise RuntimeError("table has less than 2 rows.")
    cdef uint32_t ncol = len(table[0])
    cdef uint32_t row_idx
    for row_idx in range(1, nrow):
        if len(table[row_idx]) != ncol:
            raise RuntimeError("table rows have unequal lengths.")
    if ncol < 2:
        raise RuntimeError("table has less than 2 columns.")
    cdef int32_t m11 = table[0][0]
    cdef int32_t m12 = table[0][1]
    cdef int32_t m21 = table[1][0]
    cdef int32_t m22 = table[1][1]
    cdef int64_t total = <int64_t>(m11) + <int64_t>(m12) + <int64_t>(m21) + <int64_t>(m22)
    if m11 < 0 or m12 < 0 or m21 < 0 or m22 < 0 or total > 0x7fffffff:
        raise RuntimeError("table entries must be nonnegative and sum to <2^31")
    cdef bint m11_is_greater_alt = 0
    cdef int32_t m13
    cdef int32_t m23
    cdef double ln_result
    if alternative != "two-sided":
        if nrow > 2 or ncol > 2:
            raise RuntimeError("alternative must be 'two-sided' for tables larger than 2x2.")
        m11_is_greater_alt = (alternative == "greater")
        if alternative != "less" and not m11_is_greater_alt:
            raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
        ln_result = Fisher22OneSidedLnP(m11, m12, m21, m22, m11_is_greater_alt, midp)
    else:
        if nrow == 2 and ncol == 2:
            if Fisher22LnP(m11, m12, m21, m22, midp, &ln_result):
                raise MemoryError()
        elif (nrow == 2 and ncol == 3) or (nrow == 3 and ncol == 2):
            if nrow == 2:
                m13 = table[0][2]
                m23 = table[1][2]
            else:
                m12, m21 = m21, m12
                m13 = table[2][0]
                m23 = table[2][1]
            total += <int64_t>(m13) + <int64_t>(m23)
            if m13 < 0 or m23 < 0 or total > 0x7fffffff:
                raise RuntimeError("table entries must be nonnegative and sum to <2^31")
            with nogil:
                if Fisher23LnP(m11, m12, m13, m21, m22, m23, midp, &ln_result):
                    raise MemoryError()
        else:
             raise RuntimeError("tables larger than 2x3 not yet supported")
    if logp:
        return ln_result
    return exp(ln_result)
