# cython: language_level=3
from libc.stdint cimport int64_t, uint32_t, int32_t
from libc.math cimport NAN
import fractions

__version__ = "0.3.3"

cdef extern from "../include/plink2_highprec.h" namespace "plink2":
    cdef struct dd_real_struct:
        double x[2]

    dd_real_struct ddr_maked(const double a) nogil

    dd_real_struct ddr_ieee_sub(const dd_real_struct a, const dd_real_struct b) nogil

    dd_real_struct ddr_accurate_div(const dd_real_struct a, const dd_real_struct b) nogil

    int32_t ddr_is_zero(const dd_real_struct a) nogil

    int32_t ddr_is_one(const dd_real_struct a) nogil

    int32_t ddr_ltd(const dd_real_struct a, double b) nogil

    int32_t ddr_leqd(const dd_real_struct a, double b) nogil


cdef extern from "../include/binom.h" namespace "plink2":
    ctypedef uint32_t BoolErr

    double BinomMass(int64_t k, int64_t n, dd_real_struct p_ddr, uint32_t logp) nogil

    BoolErr BinomLnP(int32_t obs_succ, int32_t obs_tot, int64_t succ_odds_ratio_numer, int64_t succ_odds_ratio_denom, int32_t midp, double* resultp) nogil

    double BinomOneSidedLnP(int64_t obs_succ, int64_t obs_tot, double succ_odds_ratio, uint32_t succ_is_greater_alt, int32_t midp) nogil

    double Pbinom(int64_t obs_k, int64_t n, dd_real_struct p_ddr, uint32_t logp) nogil


cdef extern from "../include/fisher.h" namespace "plink2":
    BoolErr Fisher22LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, double* resultp) nogil

    double Fisher22OneSidedLnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp) nogil

    BoolErr Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp, double* resultp) nogil


cdef extern from "../include/plink2_float.h" namespace "plink2":
    cdef enum:
        kLn2

    double flush_if_denormal(double xx) nogil

    double exp_flush(double xx) nogil


cdef extern from "../include/plink2_hwe.h" namespace "plink2":
    BoolErr HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp, double* resultp) nogil


# For dbinom() and pbinom(), we want to be able to deliver <1 ULP relative
# error (at least for n < ~10^10, p >= 2^{-53}) for the calculation the caller
# actually wants; and similarly, qbinom() should be based on a cmf
# approximation with at least 53-bit accuracy.
#
# Of course, a calculation with e.g. p=1/3 is a reasonable thing to want, yet
# 1/3 cannot be precisely represented by a float64.  We won't achieve our goal
# if we force the caller to misrepresent p by up to 0.5 ULPs up front.
#
# So these functions allow p to be a fractions.Fraction, and in that case we
# convert to a dd_real ("double-double") with ~106-bit precision instead of the
# usual 53.  Then we perform the underlying calculation with >53 bits of
# precision, unless n is too large for that to be practical.
cdef dd_real_struct DdrMake(object p):
    cdef dd_real_struct p_ddr
    if isinstance(p, float):
        p_ddr.x[0] = p
        p_ddr.x[1] = 0.0
        return p_ddr
    if not isinstance(p, fractions.Fraction):
        p = fractions.Fraction(p)
    p_ddr.x[0] = float(p)
    p_ddr.x[1] = float(p - fractions.Fraction(p_ddr.x[0]))
    return p_ddr


# n must be in [0, 2^31) for the two-sided test, and [0, 2^52) for one-sided
# tests.  k must be in [0, n].  p must be in [0, 1].
#
# For the two-sided test, fractions.Fraction(p) should have denominator < 2^63.
# This may require you to clarify your intent when passing in p < 1/512:
# - If an exact decimal value (e.g. 0.0001) is intended, you can pass it in as
#   a string.
# - If another sort of exact fraction is intended, pass in
#     fractions.Fraction(numer, denom)
# - Otherwise, you can use the generic approximation
#     fractions.Fraction(p).limit_denominator(2**63 - 1)
# The benefit of this interface is that it enables *perfect* handling of
# likelihood near-ties; we aren't limited to even qbinom()'s ~53 bits.
#
# Since there is negligible practical difference between ~45-bit and 53-bit
# p-value precision, the implementation aims for the former rather than the
# latter.  (This may seem inconsistent with the finicky handling of likelihood
# near-ties, but there's an important distinction: in the rare scenarios where
# a likelihood near-tie does show up, handling it incorrectly results in
# *large* relative error.)
def binom(int64_t k, int64_t n, object p=0.5, str alternative="two-sided", bint midp=0, bint logp=0):
    if k < 0 or k > n:
        raise RuntimeError("k must be nonnegative and <= n.")
    cdef double ln_result
    if alternative == "two-sided":
        if n > 0x7fffffff:
            raise RuntimeError("For alternative='two-sided', n must be less than 2^31.")
        numer, denom = p.as_integer_ratio()
        if denom == 1:
            # Degenerate cases.
            if numer != 0 and numer != 1:
                raise RuntimeError("For alternative='two-sided', p must be 0 or in (2^{-63}, 1].")
            if (numer == 0 and k == 0) or (numer == 1 and k == n):
                if midp:
                    ln_result = -kLn2
                else:
                    ln_result = 0.0
            else:
                if not logp:
                    return 0.0
                ln_result = NAN
        else:
            if numer * (2**63) <= denom or numer >= denom:
                raise RuntimeError("For alternative='two-sided', p must be 0 or in (2^{-63}, 1].")
            if denom >= 2**63:
                raise RuntimeError("For alternative='two-sided', fractions.Fraction(p) must have denominator < 2^63.  If an exact decimal value (e.g. 0.0001) is intended, it can be passed in as a string, otherwise a generic workaround is to pass in fractions.Fraction(p).limit_denominator(2**63 - 1).")
            if BinomLnP(k, n, <int64_t>(numer), <int64_t>(denom - numer), midp, &ln_result):
                raise MemoryError()
        if logp:
            return ln_result
        return exp_flush(ln_result)
    cdef bint k_is_greater_alt = (alternative == "greater")
    if alternative != "less" and not k_is_greater_alt:
        raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
    if n >= (1LL << 52):
        raise RuntimeError("n must be less than 2^52.")
    cdef double pfloat = float(p)
    if pfloat == 0.0:
        # Degenerate cases.
        if k == 0 and midp:
            ln_result = -kLn2
        elif k == 0 or not k_is_greater_alt:
            ln_result = 0.0
        else:
            if not logp:
                return 0.0
            ln_result = NAN
    elif pfloat == 1.0:
        if k == n and midp:
            ln_result = -kLn2
        elif k == n or k_is_greater_alt:
            ln_result = 0.0
        else:
            if not logp:
                return 0.0
            ln_result = NAN
    else:
        if pfloat < 0.5**960 or not pfloat <= 1.0:
            # Don't think it's worth supporting p in (0, 2^{-960}); we'd need
            # to either slow down the common case to avoid underflow/overflow,
            # or bloat the library with a separate function just to handle that
            # case.
            raise RuntimeError("p must be 0, or in [2^{-960}, 1].")
        ln_result = BinomOneSidedLnP(k, n, pfloat / (1 - pfloat), k_is_greater_alt, midp)
    if logp:
        return ln_result
    return exp_flush(ln_result)


# Returns likelihood of exactly k successes.  Relative error should be <0.6 ULP
# unless n is huge.
def dbinom(int64_t k, int64_t n, object p=0.5, bint logp=0):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    cdef dd_real_struct p_ddr = DdrMake(p)
    if ddr_ltd(p_ddr, 0.0) or not ddr_leqd(p_ddr, 1.0):
        raise RuntimeError("p must be in [0, 1].")
    if k < 0 or k > n:
        if logp:
            return NAN
        return 0.0
    if ddr_is_zero(p_ddr) or ddr_is_one(p_ddr):
        if (ddr_is_zero(p_ddr) and k == 0) or (k == n and not ddr_is_zero(p_ddr)):
            if logp:
                return 0.0
            return 1.0
        if logp:
            return NAN
        return 0.0
    return BinomMass(k, n, p_ddr, logp)

# Returns cumulative mass function, e.g. pbinom(n, n) is 1.
#
# If approx=True, this is essentially equivalent to a binom() call with
# alternative="less", which uses a faster algorithm that doesn't try to get the
# last few mantissa bits right.
# Otherwise, relative error should be <0.6 ULP unless n is huge.
def pbinom(int64_t k, int64_t n, object p=0.5, bint logp=0, bint approx=0):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    cdef dd_real_struct p_ddr = DdrMake(p)
    if not ddr_is_zero(p_ddr) and (ddr_ltd(p_ddr, 0.5**960) or not ddr_leqd(p_ddr, 1.0)):
        raise RuntimeError("p must be 0, or in [2^{-960}, 1].")
    cdef double ln_result
    if k < 0:
        # Degenerate cases.
        if not logp:
            return 0.0
        ln_result = NAN
    elif k >= n or ddr_is_zero(p_ddr):
        ln_result = 0.0
    elif ddr_is_one(p_ddr):
        if not logp:
            return 0.0
        ln_result = NAN
    else:
        if not approx:
            return flush_if_denormal(Pbinom(k, n, p_ddr, logp))
        ln_result = BinomOneSidedLnP(k, n, ddr_accurate_div(p_ddr, ddr_ieee_sub(ddr_maked(1.0), p_ddr)).x[0], 0, 0)
    if logp:
        return ln_result
    return exp_flush(ln_result)

# Todo: qbinom function.


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
    return exp_flush(ln_result)


# "HWE" is short for Hardy-Weinberg Equilibrium.
#
# hom1, hets, and hom2 must be nonnegative, and add up to <2^31.
#
# alternative="less" and alternative="greater" refer to the heterozygote count.
# (alternative="greater" has more practical value in identifying
# variant-calling errors without throwing out variants affected by the Wahlund
# effect.)  These one-sided tests are not implemented yet.
#
# Variants with k>2 alleles can be evaluated with k one-vs.-rest tests.
def HWE(int32_t hom1, int32_t hets, int32_t hom2, str alternative="two-sided", bint midp=0, bint logp=0):
    cdef int64_t total = <int64_t>(hom1) + <int64_t>(hets) + <int64_t>(hom2)
    if hom1 < 0 or hets < 0 or hom2 < 0 or total > 0x7fffffff:
        raise RuntimeError("hom1, hets, and hom2 must be nonnegative and sum to <2^31")
    cdef bint hets_is_greater_alt = 0
    cdef double ln_result
    if alternative != "two-sided":
        hets_is_greater_alt = (alternative == "greater")
        if alternative != "less" and not hets_is_greater_alt:
            raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
        raise RuntimeError("one-sided tests not implemented yet")
    else:
        # note different parameter order
        if HweLnP(hets, hom1, hom2, midp, &ln_result):
            raise MemoryError()
    if logp:
        return ln_result
    return exp_flush(ln_result)
