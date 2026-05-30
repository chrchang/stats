# cython: language_level=3
from libc.stdint cimport int64_t, uint32_t, int32_t
from libc.math cimport NAN
import fractions
import random, scipy
from libc.math cimport sqrt

__version__ = "0.4.5"

cdef extern from "../include/plink2_highprec.h" namespace "plink2":
    cdef struct qd_real_struct:
        double x[2]

    cdef struct dd_real_struct:
        double x[2]

    dd_real_struct ddr_maked(const double a) nogil

    dd_real_struct ddr_make_qd(const qd_real_struct a) nogil

    dd_real_struct ddr_negate(const dd_real_struct a) nogil

    dd_real_struct ddr_subd(const dd_real_struct a, double b) nogil

    int32_t ddr_is_zero(const dd_real_struct a) nogil

    int32_t ddr_is(const dd_real_struct a, double b) nogil

    int32_t ddr_ltd(const dd_real_struct a, double b) nogil

    int32_t ddr_leqd(const dd_real_struct a, double b) nogil

    qd_real_struct qdr_addd(const qd_real_struct a, double b) nogil

    int32_t qdr_is_zero(const qd_real_struct a) nogil

    int32_t qdr_is(const qd_real_struct a, double b) nogil

    int32_t qdr_ltd(const qd_real_struct a, double b) nogil

    int32_t qdr_leqd(const qd_real_struct a, double b) nogil


cdef extern from "../include/binom.h" namespace "plink2":
    double BinomMass(int64_t k, int64_t n, dd_real_struct p_ddr, dd_real_struct q_ddr, uint32_t logp) nogil

    double PbinomApprox(int64_t obs_k, int64_t n, dd_real_struct p_ddr, dd_real_struct q_ddr, uint32_t complement, int32_t midp, uint32_t logp) nogil

    double Pbinom(int64_t obs_k, int64_t n, dd_real_struct p_ddr, dd_real_struct q_ddr, uint32_t complement, uint32_t logp) nogil

    int64_t QbinomHalfUlp(dd_real_struct targetp_or_lnp_ddr, int64_t n, dd_real_struct distp_ddr, uint32_t log_targetp) nogil

    double BinomTwoSidedP(int32_t obs_succ, int32_t obs_tot, qd_real_struct p_qdr, int32_t midp, uint32_t logp) nogil


cdef extern from "../include/fisher.h" namespace "plink2":
    double HypergeomMass(int64_t m11, int64_t m12, int64_t m21, int64_t m22, uint32_t logp)

    double Fisher22TwoSidedP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, uint32_t logp) nogil

    double PhyperApprox(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp, uint32_t logp) nogil

    double Phyper(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t logp) nogil

    int64_t QhyperHalfUlp(dd_real_struct p_or_lnp_ddr, int64_t ac, int64_t bd, int64_t ab, uint32_t logp) nogil

    double Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp) nogil


cdef extern from "../include/plink2_float.h" namespace "plink2":
    cdef enum:
        kLn2

    double flush_if_denormal(double xx) nogil

    double exp_flush(double xx) nogil


cdef extern from "../include/plink2_hwe.h" namespace "plink2":
    double HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp) nogil


# For dbinom() and pbinom(), we want to be able to deliver <1 ULP relative
# error for the calculation the caller actually wants; and similarly, qbinom()
# should be based on a cmf approximation with at least 53-bit accuracy.
#
# Of course, a calculation with e.g. p=1/3 is a reasonable thing to want, yet
# 1/3 cannot be precisely represented by a float64.  We're poorly positioned to
# achieve our goal if we force the caller to misrepresent p by up to 0.5 ULPs
# up front.
#
# So these functions allow p to be a fractions.Fraction, and in that case we
# convert to a dd_real ("double-double") with ~106-bit precision instead of the
# usual 53.  Then we perform the underlying calculation with >53 bits of
# precision, unless n is too large for that to be practical.
#
# We also return 1-p as a dd_real, since this is occasionally important for p
# near 1.
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

cdef dd_real_struct DdrMake2(object p, dd_real_struct* q_ddr_ptr):
    cdef dd_real_struct p_ddr
    if isinstance(p, float):
        p_ddr.x[0] = p
        p_ddr.x[1] = 0.0
        q_ddr_ptr[0].x[0] = 1.0 - p
        q_ddr_ptr[0].x[1] = 1.0 - (1.0 - p) - p
        return p_ddr
    if not isinstance(p, fractions.Fraction):
        p = fractions.Fraction(p)
    p_ddr.x[0] = float(p)
    p_ddr.x[1] = float(p - fractions.Fraction(p_ddr.x[0]))
    q = fractions.Fraction(fractions.Fraction(1, 1) - p)
    q_ddr_ptr[0].x[0] = float(q)
    q_ddr_ptr[0].x[1] = float(q - fractions.Fraction(q_ddr_ptr[0].x[0]))
    return p_ddr

# For two-sided binom(), we convert p to a qd_real with ~212-bit precision, so
# we can use a very small epsilon for near-tie resolution.
cdef qd_real_struct QdrMake(object p):
    cdef qd_real_struct p_qdr
    if isinstance(p, float):
        p_qdr.x[0] = p
        p_qdr.x[1] = 0.0
        p_qdr.x[2] = 0.0
        p_qdr.x[3] = 0.0
        return p_qdr
    if not isinstance(p, fractions.Fraction):
        p = fractions.Fraction(p)
    p_qdr.x[0] = float(p)
    rem1 = p - fractions.Fraction(p_qdr.x[0])
    p_qdr.x[1] = float(rem1)
    rem2 = rem1 - fractions.Fraction(p_qdr.x[1])
    p_qdr.x[2] = float(rem2)
    p_qdr.x[3] = float(rem2 - fractions.Fraction(p_qdr.x[2]))
    return p_qdr

cdef double zeroval(bint logp):
    if logp:
        return NAN
    return 0.0

cdef double oneval(bint logp):
    return 1.0 - logp

cdef double half_or_oneval(bint midp, bint logp):
    if midp:
        if logp:
            return -kLn2
        return 0.5
    return 1.0 - logp

cdef double pbinom_p01(int64_t k, int64_t n, double p, bint complement, bint midp, bint logp):
    if complement:
        k = n + midp - k - 1
        p = 1 - p
    if k < 0:
        return zeroval(logp)
    if k > n:
        return oneval(logp)
    if p == 0:
        if k == 0:
            return half_or_oneval(midp, logp)
        return oneval(logp)
    # p == 1
    if k == n:
        return half_or_oneval(midp, logp)
    return zeroval(logp)


# n must be in [0, 2^52), k must be in [0, n], p must be in [0, 1].
#
# If p is a fractions.Fraction(), it's expanded to a "quad-double" with ~212
# bit accuracy.  This enables very accurate handling of near-ties.
def binom(int64_t k, int64_t n, object p=0.5, str alternative="two-sided", bint midp=0, bint logp=0):
    if k < 0 or k > n:
        raise RuntimeError("k must be nonnegative and <= n.")
    if n >= (1LL << 52):
        raise RuntimeError("n must be less than 2^52.")
    cdef bint complement = (alternative == "greater")
    if alternative != "two-sided" and alternative != "less" and not complement:
        raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
    if complement and (not midp):
        k -= 1
    cdef qd_real_struct p_qdr = QdrMake(p)
    if qdr_is_zero(p_qdr) or qdr_is(p_qdr, 1):
        # Degenerate cases.
        if alternative == "two-sided":
            if (p_qdr.x[0] == 0 and k == 0) or (p_qdr.x[0] == 1 and k == n):
                return half_or_oneval(midp, logp)
            return zeroval(logp)
        return pbinom_p01(k, n, p_qdr.x[0], complement, midp, logp)
    if (qdr_ltd(p_qdr, 0.5**960) or not qdr_leqd(p_qdr, 1.0)):
        # TODO: these functions should allow p in (0, 2^{-960}).  Deferred
        # since, as of this writing, there are much higher-priority problems to
        # solve; but this is straightforward to get right, just need to be
        # careful about underflow.
        raise RuntimeError("p must be 0, or in [2^{-960}, 1].")
    if alternative == "two-sided":
        return flush_if_denormal(BinomTwoSidedP(k, n, p_qdr, midp, logp))
    cdef dd_real_struct q_ddr = ddr_negate(ddr_make_qd(qdr_addd(p_qdr, -1.0)))
    return flush_if_denormal(PbinomApprox(k, n, ddr_make_qd(p_qdr), q_ddr, complement, midp, logp))


# Returns likelihood of exactly k successes.  Relative error should be <0.6 ULP
# unless n is huge.
def dbinom(int64_t k, int64_t n, object p=0.5, bint logp=0):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    cdef dd_real_struct q_ddr
    cdef dd_real_struct p_ddr = DdrMake2(p, &q_ddr)
    if ddr_ltd(p_ddr, 0.0) or not ddr_leqd(p_ddr, 1.0):
        raise RuntimeError("p must be in [0, 1].")
    if k < 0 or k > n:
        if logp:
            return NAN
        return 0.0
    if ddr_is_zero(p_ddr) or ddr_is(p_ddr, 1):
        if (ddr_is_zero(p_ddr) and k == 0) or (k == n and not ddr_is_zero(p_ddr)):
            return oneval(logp)
        return zeroval(logp)
    return flush_if_denormal(BinomMass(k, n, p_ddr, q_ddr, logp))


# Returns cumulative mass function, e.g. pbinom(n, n) is 1.
#
# If approx=True, this is essentially equivalent to a binom() call with
# alternative="less", which uses a faster algorithm that doesn't try to get the
# last few mantissa bits right.
# Otherwise, relative error should be <0.6 ULP unless n is huge.
def pbinom(int64_t k, int64_t n, object p=0.5, bint complement=0, bint logp=0, bint approx=0):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    cdef dd_real_struct q_ddr
    cdef dd_real_struct p_ddr = DdrMake2(p, &q_ddr)
    if not ddr_is_zero(p_ddr) and (ddr_ltd(p_ddr, 0.5**960) or not ddr_leqd(p_ddr, 1.0)):
        raise RuntimeError("p must be 0, or in [2^{-960}, 1].")
    if ddr_is_zero(p_ddr) or ddr_is(p_ddr, 1):
        return pbinom_p01(k, n, p_ddr.x[0], complement, 0, logp)
    if approx:
        return flush_if_denormal(PbinomApprox(k, n, p_ddr, q_ddr, complement, 0, logp))
    return flush_if_denormal(Pbinom(k, n, p_ddr, q_ddr, complement, logp))


# Returns smallest nonnegative k for which cdf(k) >= targetP if logTarget is
# is False, and cdf(k) >= exp(targetP) if logTarget is True.
#
# Implementation is *not* built on top of pbinom() in a way that e.g.
# guarantees qbinom(pbinom(k, n, succP), n, succP) == k or
# qbinom(pbinom(k, n, succP) * (1 + 0.5**52), n, succP) > k in non-degenerate
# cases.  However, it is designed to make these outcomes very likely:
# - Qbinom() is designed for <0.6 ULP relative error (except when n is huge),
#   and achieves <0.5 ULP the vast majority of the time.
# - The internal Qbinom() call is made with 0.5 ULP subtracted off of q.
def qbinom(object targetP, int64_t n, object succP=0.5, bint logTarget=0):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    cdef dd_real_struct distp_ddr = DdrMake(succP)
    if not ddr_is_zero(distp_ddr) and (ddr_ltd(distp_ddr, 0.5**960) or not ddr_leqd(distp_ddr, 1.0)):
        raise RuntimeError("succP must be 0, or in [2^{-960}, 1].")
    cdef dd_real_struct targetp_or_lnp_ddr = DdrMake(targetP)
    if logTarget:
        if not ddr_leqd(targetp_or_lnp_ddr, 0.0):
            raise RuntimeError("targetP must be <= 0 when logTarget is True.")
    else:
        if ddr_ltd(targetp_or_lnp_ddr, 0.0) or not ddr_leqd(targetp_or_lnp_ddr, 1.0):
            raise RuntimeError("targetP must be in [0, 1] when logTarget is False.")
    return QbinomHalfUlp(targetp_or_lnp_ddr, n, distp_ddr, logTarget)


# table must be a 2x2 or larger matrix, represented as a list of equal-length
# lists.  For two-sided tests, values must be nonnegative integers which add up
# to <2^31.  For one-sided tests, they must add up to <2^52.
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
    cdef int64_t m11 = table[0][0]
    cdef int64_t m12 = table[0][1]
    cdef int64_t m21 = table[1][0]
    cdef int64_t m22 = table[1][1]
    cdef int64_t total = m11 + m12 + m21 + m22
    if m11 < 0 or m12 < 0 or m21 < 0 or m22 < 0:
        raise RuntimeError("table entries must be nonnegative")
    if m11 >= (1LL << 52) or m12 >= (1LL << 52) or m21 >= (1LL << 52) or m22 >= (1LL << 52):
        raise RuntimeError("table entries must be <2^52")
    cdef bint m11_is_greater_alt = 0
    cdef double ln_result
    if alternative != "two-sided":
        if nrow > 2 or ncol > 2:
            raise RuntimeError("alternative must be 'two-sided' for tables larger than 2x2.")
        if total >= (1LL << 52):
            raise RuntimeError("table entries must sum to <2^52")
        m11_is_greater_alt = (alternative == "greater")
        if alternative != "less" and not m11_is_greater_alt:
            raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
        return flush_if_denormal(PhyperApprox(m11, m12, m21, m22, m11_is_greater_alt, midp, logp))
    if nrow == 2 and ncol == 2:
        if total > 0x7fffffff:
            raise RuntimeError("table entries must sum to <2^31")
        ln_result = Fisher22TwoSidedP(m11, m12, m21, m22, midp, logp)
        return flush_if_denormal(ln_result)
    cdef int32_t m13
    cdef int32_t m23
    if (nrow == 2 and ncol == 3) or (nrow == 3 and ncol == 2):
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
            ln_result = Fisher23LnP(m11, m12, m13, m21, m22, m23, midp)
    else:
         raise RuntimeError("tables larger than 2x3 not yet supported")
    if logp:
        return ln_result
    return exp_flush(ln_result)


# This mirrors R dhyper()'s parameters.
def dhyper(int64_t a, int64_t ac, int64_t bd, int64_t ab, bint logp=0):
    if a < 0 or ac < 0 or bd < 0 or ab < 0 or a >= (1LL << 52) or ac >= (1LL << 52) or bd >= (1LL << 52) or ab >= (1LL << 52):
        raise RuntimeError("Parameters must be in [0, 2^52).")
    cdef int64_t b = ab - a
    cdef int64_t c = ac - a
    cdef int64_t d = bd - b
    if b < 0 or c < 0 or d < 0:
        raise RuntimeError("(ab-a), (ac-a), and (bd-(ab-a)) must be nonnegative.")
    if ac + bd >= (1LL << 52):
        raise RuntimeError("ac+bd must be <2^52.")
    return flush_if_denormal(HypergeomMass(a, b, c, d, logp))


def phyper(int64_t a, int64_t ac, int64_t bd, int64_t ab, bint complement=0, bint logp=0, bint approx=0):
    if a < 0 or ac < 0 or bd < 0 or ab < 0 or a >= (1LL << 52) or ac >= (1LL << 52) or bd >= (1LL << 52) or ab >= (1LL << 52):
        # Unlike pbinom(), we don't bother with returning NAN/0/1 in some of
        # these cases.
        raise RuntimeError("Parameters must be in [0, 2^52).")
    cdef int64_t b = ab - a
    cdef int64_t c = ac - a
    cdef int64_t d = bd - b
    if b < 0 or c < 0 or d < 0:
        raise RuntimeError("(ab-a), (ac-a), and (bd-(ab-a)) must be nonnegative.")
    if ac + bd >= (1LL << 52):
        raise RuntimeError("ac+bd must be <2^52.")
    if complement:
        a, b = b, a
        c, d = d, c
        if a == 0 or d == 0:
            if logp:
                return NAN
            return 0
        a -= 1
        b += 1
        c += 1
        d -= 1
    if approx:
        return flush_if_denormal(PhyperApprox(a, b, c, d, 0, 0, logp))
    return flush_if_denormal(Phyper(a, b, c, d, logp))


# Returns smallest 'a' in the distribution support for which cdf(a) >= p.
#
# Implementation is *not* built on top of phyper() in a way that e.g.
# guarantees qhyper(phyper(a, ac, bd, ab), ac, bd, ab) == a or
# qhyper(phyper(a, ac, bd, ab) * (1 + 0.5**52), ac, bd, ab) > a in
# non-degenerate cases.  However, it is designed to make these outcomes very
# likely:
# - Phyper() is designed for <0.6 ULP relative error (except when n is huge),
#   and achieves <0.5 ULP the vast majority of the time.
# - The internal Qhyper() call is made with 0.5 ULP subtracted off of q.
def qhyper(object p, int64_t ac, int64_t bd, int64_t ab, bint logp=0):
    if ac < 0 or bd < 0 or ab < 0 or ac >= (1LL << 52) or bd >= (1LL << 52) or ab >= (1LL << 52):
        raise RuntimeError("ac, bd, and ab must be in [0, 2^52).")
    if ac + bd >= (1LL << 52):
        raise RuntimeError("ac+bd must be <2^52.")
    if ab > ac + bd:
        raise RuntimeError("ab must be <= ac+bd.")
    cdef dd_real_struct p_ddr = DdrMake(p)
    if logp:
        if not ddr_leqd(p_ddr, 0.0):
            raise RuntimeError("p must be <= 0 when logp is True.")
    else:
        if ddr_ltd(p_ddr, 0.0) or not ddr_leqd(p_ddr, 1.0):
            raise RuntimeError("p must be in [0, 1] when logp is False.")
    return QhyperHalfUlp(p_ddr, ac, bd, ab, logp)


# TODO: scipy-style entry points.


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
        ln_result = HweLnP(hets, hom1, hom2, midp)
    if logp:
        return ln_result
    return exp_flush(ln_result)


def pbinom_approx_accuracy_test(float p, uint32_t min_pow2, uint32_t max_pow2, uint32_t trials_per_pow2, bint logp=0):
    cdef uint32_t pow2
    cdef uint32_t trial_idx
    cdef int64_t min_n
    cdef int64_t max_n
    cdef int64_t k
    cdef int64_t n
    cdef double got
    cdef double want
    cdef double relerr
    cdef double relerr_ssq
    cdef double relerr_scipy_ssq
    cdef double approx_rms
    cdef double scipy_rms
    for pow2 in range(min_pow2, max_pow2 + 1):
        min_n = 1LL << pow2
        max_n = min_n * 2 - 1
        relerr_ssq = 0.0
        relerr_scipy_ssq = 0.0
        for trial_idx in range(0, trials_per_pow2):
            n = random.randint(min_n, max_n)
            k = <int64_t>(float(n) * p)
            want = pbinom(k, n, p, logp)
            got = pbinom(k, n, p, logp, approx=True)
            relerr = (got - want) / want
            relerr_ssq += relerr * relerr
            if logp:
                got = scipy.stats.binom.logcdf(k, n, p)
            else:
                got = scipy.stats.binom.cdf(k, n, p)
            relerr = (got - want) / want
            relerr_scipy_ssq += relerr * relerr
        approx_rms = sqrt(relerr_ssq / float(trials_per_pow2))
        scipy_rms = sqrt(relerr_scipy_ssq / float(trials_per_pow2))
        print("[2^" + str(pow2) + ", 2^" + str(pow2+1) + " - 1): approxRMS=" + str(approx_rms) + "  scipyRMS=" + str(scipy_rms))
