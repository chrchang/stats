# cython: language_level=3, boundscheck=False, wraparound=False
from libc.stdint cimport int64_t, uintptr_t, uint32_t, int32_t
from libc.math cimport INFINITY
import fractions
import numpy as np
cimport numpy as cnp

__version__ = "0.9.1"

cdef extern from "../include/plink2_highprec.h" namespace "plink2":
    cdef struct td_real_struct:
        double x[3]

    cdef struct dd_real_struct:
        double x[2]


    dd_real_struct ddr_negate(const dd_real_struct a) nogil

    dd_real_struct ddr_add2d(double a, double b) nogil

    dd_real_struct ddr_addd(const dd_real_struct a, double b) nogil

    int32_t ddr_ltd(const dd_real_struct a, double b) nogil

    int32_t ddr_leqd(const dd_real_struct a, double b) nogil

    dd_real_struct ddr_maked(const double a) nogil

    dd_real_struct ddr_exp(const dd_real_struct a) nogil

    td_real_struct tdr_make1(const double a) nogil

    int32_t tdr_ltd(const td_real_struct a, double b) nogil

    int32_t tdr_leqd(const td_real_struct a, double b) nogil


cdef extern from "../include/binom_detail.h" namespace "plink2":
    void BinomMassMultiKPrecomp(int64_t n, td_real_struct p_tdr, uint32_t* p_is_half_ptr, td_real_struct* lfact_n_tdr_ptr, td_real_struct* lnp_tdr_ptr, td_real_struct* lnq_tdr_ptr) nogil

    double BinomMassJustK(int64_t k, int64_t n, uint32_t p_is_half, const td_real_struct lfact_n_tdr, const td_real_struct lnp_tdr, const td_real_struct lnq_tdr, uint32_t logp) nogil

    void BinomMassMultiPPrecomp(int64_t k, int64_t n, td_real_struct* lfact_n_tdr_ptr, td_real_struct* neg_lfact_k_tdr_ptr, td_real_struct* neg_lfact_nmk_tdr_ptr)

    double BinomMassJustP(td_real_struct p_tdr, int64_t k, int64_t n, const td_real_struct lfact_n_tdr, const td_real_struct neg_lfact_k_tdr, const td_real_struct neg_lfact_nmk_tdr, uint32_t logp)


cdef extern from "../include/binom.h" namespace "plink2":
    double BinomMass(int64_t k, int64_t n, td_real_struct p_tdr, uint32_t logp) nogil

    double PbinomApprox(int64_t obs_k, int64_t n, td_real_struct p_tdr, uint32_t complement, int32_t midp, uint32_t logp) nogil

    double Pbinom(int64_t obs_k, int64_t n, td_real_struct p_tdr, uint32_t complement, uint32_t logp) nogil

    int64_t QbinomHalfUlp(dd_real_struct targetp_or_lnp_ddr, int64_t n, td_real_struct distp_tdr, uint32_t log_targetp) nogil

    double BinomTwoSidedP(int32_t obs_succ, int32_t obs_tot, td_real_struct p_tdr, int32_t midp, uint32_t logp) nogil


cdef extern from "../include/hypergeom_detail.h" namespace "plink2":
    void HypergeomMassMultiKPrecomp(int64_t mxx, int64_t m1x, int64_t mx1, td_real_struct* lfact_m1x_tdr_ptr, td_real_struct* lfact_m2x_tdr_ptr, td_real_struct* lfact_mx1_tdr_ptr, td_real_struct* lfact_mx2_tdr_ptr, td_real_struct* lfact_mxx_tdr_ptr) nogil

    double HypergeomMassJustK(int64_t m11, int64_t mxx, int64_t m1x, int64_t mx1, const td_real_struct lfact_m1x_tdr, const td_real_struct lfact_m2x_tdr, const td_real_struct lfact_mx1_tdr, const td_real_struct lfact_mx2_tdr, const td_real_struct lfact_mxx_tdr, uint32_t logp) nogil


cdef extern from "../include/hypergeom.h" namespace "plink2":
    double HypergeomMass(int64_t m11, int64_t m12, int64_t m21, int64_t m22, uint32_t logp)

    double PhyperApprox(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp, uint32_t logp) nogil

    double Phyper(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t logp) nogil

    int64_t QhyperHalfUlp(dd_real_struct p_or_lnp_ddr, int64_t ac, int64_t bd, int64_t ab, uint32_t logp) nogil


cdef extern from "../include/fisher.h" namespace "plink2":
    double Fisher22TwoSidedP(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, int32_t midp, uint32_t logp) nogil

    double Fisher22OddsRatio(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22) nogil

    void Fisher22OddsRatioCI(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double low_p, double high_p, double* low_resultp, double* high_resultp) nogil

    double Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp) nogil


cdef extern from "../include/plink2_float.h" namespace "plink2":
    cdef enum:
        kLn2

    double flush_if_denormal(double xx) nogil

    double exp_flush(double xx) nogil


cdef extern from "../include/plink2_hwe.h" namespace "plink2":
    double HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp) nogil


# The pmf and cmf functions here default to delivering <1 ULP relative error;
# and similarly, the quantile functions are designed to correspond to cmf
# approximations with at least 53-bit accuracy everywhere.
#
# However, there can be additional error from imprecise representation of p,
# e.g. the user actually wants to perform the calculation for p=1/3, but the
# closest float64 is off by 1/3 ULP.  Because the underlying C-library
# functions already need high-precision arithmetic to handle log-factorials
# well, it doesn't take much of an additional stretch to support a
# high-precision representation of p.  So they actually accept p as a td_real
# ("triple-double") with ~159-bit precision instead of the usual 53; and this
# TdrMake() function converts a fractions.Fraction-compatible Python object to
# a td_real.
#
# For now, this capability is only exposed in the Python API by the binomtest()
# function.  It has extra value there because the usual alternative (for the
# default two-sided test) is to have a wide tie-detection epsilon, which makes
# it impossible to handle some large cases accurately.
cdef td_real_struct TdrMake(object p):
    if isinstance(p, (float, np.float64, np.float32)):
        return tdr_make1(p)
    if not isinstance(p, fractions.Fraction):
        p = fractions.Fraction(p)
    cdef td_real_struct p_tdr
    p_tdr.x[0] = float(p)
    rem1 = p - fractions.Fraction(p_tdr.x[0])
    p_tdr.x[1] = float(rem1)
    p_tdr.x[2] = float(rem1 - fractions.Fraction(p_tdr.x[1]))
    return p_tdr

cdef double zeroval(bint logp) nogil:
    if logp:
        return -INFINITY
    return 0.0

cdef double oneval(bint logp) nogil:
    return 1.0 - logp


cdef double dbinom_internal(int64_t k, int64_t n, double p, bint logp) nogil:
    if k < 0 or k > n:
        if logp:
            return -INFINITY
        else:
            return 0.0
    if p == 0.0 or p == 1.0:
        if (p == 0.0 and k == 0) or (k == n and p == 1.0):
            return oneval(logp)
        return zeroval(logp)
    return flush_if_denormal(BinomMass(k, n, tdr_make1(p), logp))

cdef dbinom_vectorize_k(object k_obj, int64_t n, double p, bint logp):
    ka = np.asarray(k_obj, dtype=np.int64)
    cdef int64_t [::1] kar = ka.ravel()
    cdef uintptr_t kar_size = kar.size
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(kar_size, dtype=np.float64)
    cdef double out_of_support_result
    if logp:
        out_of_support_result = -INFINITY
    else:
        out_of_support_result = 0.0
    cdef uintptr_t idx
    cdef int64_t ki
    cdef double result
    if p == 0.0 or p == 1.0:
        with nogil:
            for idx in range(kar_size):
                ki = kar[idx]
                if ki < 0 or ki > n:
                    result = out_of_support_result
                elif (p == 0.0 and ki == 0) or (ki == n and p == 1.0):
                    result = oneval(logp)
                else:
                    result = zeroval(logp)
                results[idx] = result
        return np.reshape(results, ka.shape)
    cdef uint32_t p_is_half
    cdef td_real_struct lfact_n_tdr
    cdef td_real_struct lnp_tdr
    cdef td_real_struct lnq_tdr
    with nogil:
        BinomMassMultiKPrecomp(n, tdr_make1(p), &p_is_half, &lfact_n_tdr, &lnp_tdr, &lnq_tdr)
        for idx in range(kar_size):
            ki = kar[idx]
            if ki < 0 or ki > n:
                result = out_of_support_result
            else:
                result = flush_if_denormal(BinomMassJustK(ki, n, p_is_half, lfact_n_tdr, lnp_tdr, lnq_tdr, logp))
            results[idx] = result
    return np.reshape(results, ka.shape)

cdef dbinom_v_internal(object k_obj, int64_t n, double p, bint logp):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    cdef int64_t ki
    try:
        ki = k_obj
    except TypeError:
        return dbinom_vectorize_k(k_obj, n, p, logp)
    return dbinom_internal(ki, n, p, logp)

cdef dbinom_vectorize_all(object k_obj, object n_obj, object p_obj, bint logp):
    ka = np.asarray(k_obj, dtype=np.int64)
    na = np.asarray(n_obj, dtype=np.int64)
    pa = np.asarray(p_obj, dtype=np.float64)
    it = np.nditer([ka, na, pa], flags=['c_index'])
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.float64)
    # important to declare these types, otherwise loops are a lot slower
    cdef int64_t ki
    cdef int64_t ni
    cdef double pd
    for ki, ni, pd in it:
        if ni < 0 or ni >= (1LL << 52):
            raise RuntimeError("n must be in [0, 2^52).")
        results[it.index] = dbinom_internal(ki, ni, pd, logp)
    return np.reshape(results, np.broadcast_shapes(ka.shape, na.shape, pa.shape))

cdef dbinom_vectorize_kp(object k_obj, int64_t n, object p_obj, bint logp):
    cdef int64_t ki
    try:
        ki = k_obj
    except TypeError:
        return dbinom_vectorize_all(k_obj, n, p_obj, logp)
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    # This case isn't that uncommon, and is straightforward to optimize.
    pa = np.asarray(p_obj, dtype=np.float64)
    cdef double [::1] par = pa.ravel()
    cdef uintptr_t par_size = par.size
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(par_size, dtype=np.float64)
    cdef double out_of_support_result
    if logp:
        out_of_support_result = -INFINITY
    else:
        out_of_support_result = 0.0
    cdef double pd
    if ki < 0 or ki > n:
        for idx in range(par_size):
            pd = par[idx]
            if pd < 0.0 or not pd <= 1.0:
                raise RuntimeError("p must be in [0, 1].")
            results[idx] = out_of_support_result
        return np.reshape(results, pa.shape)
    cdef td_real_struct lfact_n_tdr
    cdef td_real_struct neg_lfact_k_tdr
    cdef td_real_struct neg_lfact_nmk_tdr
    BinomMassMultiPPrecomp(ki, n, &lfact_n_tdr, &neg_lfact_k_tdr, &neg_lfact_nmk_tdr)
    cdef double result
    for idx in range(par_size):
        pd = par[idx]
        if pd < 0.0 or not pd <= 1.0:
            raise RuntimeError("p must be in [0, 1].")
        if pd == 0.0 or pd == 1.0:
            if (pd == 0.0 and ki == 0) or (ki == n and pd == 1.0):
                result = oneval(logp)
            else:
                result = zeroval(logp)
        else:
            result = flush_if_denormal(BinomMassJustP(tdr_make1(pd), ki, n, lfact_n_tdr, neg_lfact_k_tdr, neg_lfact_nmk_tdr, logp))
        results[idx] = result
    return np.reshape(results, pa.shape)

cdef dbinom_vv_internal(object k_obj, object n_obj, object p_obj, bint logp):
    cdef int64_t ni
    try:
        ni = n_obj
    except TypeError:
        return dbinom_vectorize_all(k_obj, n_obj, p_obj, logp)
    cdef double pd
    try:
        pd = p_obj
    except TypeError:
        return dbinom_vectorize_kp(k_obj, ni, p_obj, logp)
    return dbinom_v_internal(k_obj, ni, pd, logp)

# Returns likelihood of exactly k successes.  Relative error should be <1 ULP.
# This R-style entry point can only broadcast k, not n or p.
def dbinom(object k, int64_t n, double p=0.5, bint logp=0):
    return dbinom_v_internal(k, n, p, logp)


cdef double pbinom_internal(int64_t k, int64_t n, td_real_struct p_tdr, bint complement, bint logp, bint approx) nogil:
    cdef double result
    if approx:
        result = PbinomApprox(k, n, p_tdr, complement, 0, logp)
    else:
        result = Pbinom(k, n, p_tdr, complement, logp)
    return flush_if_denormal(result)

cdef pbinom_vectorize_k(object k_obj, int64_t n, double p, bint complement, bint logp, bint approx):
    ka = np.asarray(k_obj, dtype=np.int64)
    cdef int64_t [::1] kar = ka.ravel()
    cdef uintptr_t kar_size = kar.size
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(kar_size, dtype=np.float64)
    cdef td_real_struct p_tdr = tdr_make1(p)
    cdef uintptr_t idx
    cdef int64_t ki
    with nogil:
        for idx in range(kar_size):
            ki = kar[idx]
            results[idx] = pbinom_internal(ki, n, p_tdr, complement, logp, approx)
    return np.reshape(results, ka.shape)

cdef pbinom_v_internal(object k_obj, int64_t n, double p, bint complement, bint logp, bint approx):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    if p < 0.0 or not p <= 1.0:
        raise RuntimeError("p must be in [0, 1].")
    cdef int64_t ki
    try:
        ki = k_obj
    except TypeError:
        return pbinom_vectorize_k(k_obj, n, p, complement, logp, approx)
    return pbinom_internal(ki, n, tdr_make1(p), complement, logp, approx)

cdef pbinom_vectorize_all(object k_obj, object n_obj, object p_obj, bint complement, bint logp, bint approx):
    ka = np.asarray(k_obj, dtype=np.int64)
    na = np.asarray(n_obj, dtype=np.int64)
    pa = np.asarray(p_obj, dtype=np.float64)
    it = np.nditer([ka, na, pa], flags=['c_index'])
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.float64)
    cdef int64_t ki
    cdef int64_t ni
    cdef double pd
    for ki, ni, pd in it:
        if ni < 0 or ni >= (1LL << 52):
            raise RuntimeError("n must be in [0, 2^52).")
        if pd < 0.0 or not pd <= 1.0:
            raise RuntimeError("p must be in [0, 1].")
        results[it.index] = pbinom_internal(ki, ni, tdr_make1(pd), complement, logp, approx)
    return np.reshape(results, np.broadcast_shapes(ka.shape, na.shape, pa.shape))

cdef pbinom_vv_internal(object k_obj, object n_obj, object p_obj, bint complement, bint logp, bint approx):
    cdef int64_t ni
    cdef double pd
    try:
        ni = n_obj
        pd = p_obj
    except TypeError:
        return pbinom_vectorize_all(k_obj, n_obj, p_obj, complement, logp, approx)
    return pbinom_v_internal(k_obj, ni, pd, complement, logp, approx)

# Returns cumulative mass function, e.g. pbinom(n, n) is 1.
#
# If approx=True, this is essentially equivalent to a binom() call with
# alternative="less", which uses a faster algorithm that doesn't try to get the
# last few mantissa bits right.
# Otherwise, relative error should be <0.6 ULP unless n is huge.
def pbinom(object k, int64_t n, double p=0.5, bint complement=0, bint logp=0, bint approx=0):
    return pbinom_v_internal(k, n, p, complement, logp, approx)


cdef int64_t qbinom_internal(double q, int64_t n, double succP, bint invert, bint logTarget) except? -1:
    cdef dd_real_struct q_or_lnq_ddr
    if invert:
        if logTarget:
            q_or_lnq_ddr = ddr_addd(ddr_negate(ddr_exp(ddr_maked(q))), 1.0)
            logTarget = False
        else:
            q_or_lnq_ddr = ddr_add2d(1.0, -q)
    else:
        q_or_lnq_ddr = ddr_maked(q)
    if logTarget:
        if not ddr_leqd(q_or_lnq_ddr, 0.0):
            raise RuntimeError("targetP must be <= 0 when logTarget is True.")
    else:
        if ddr_ltd(q_or_lnq_ddr, 0.0) or not ddr_leqd(q_or_lnq_ddr, 1.0):
            raise RuntimeError("targetP must be in [0, 1] when logTarget is False.")
    return QbinomHalfUlp(q_or_lnq_ddr, n, tdr_make1(succP), logTarget)

cdef qbinom_vectorize_targetp(object targetP_obj, int64_t n, double succP, bint invert, bint logTarget):
    qa = np.asarray(targetP_obj, dtype=np.float64)
    cdef double [::1] qar = qa.ravel()
    cdef uintptr_t qar_size = qar.size
    cdef cnp.ndarray[cnp.int64_t,mode="c",ndim=1] results = np.empty(qar_size, dtype=np.int64)
    cdef uintptr_t idx
    cdef double qd
    for idx in range(qar_size):
        qd = qar[idx]
        results[idx] = qbinom_internal(qd, n, succP, invert, logTarget)
    return np.reshape(results, qa.shape)

cdef qbinom_v_internal(object targetP_obj, int64_t n, double succP, bint invert, bint logTarget):
    if n < 0 or n >= (1LL << 52):
        raise RuntimeError("n must be in [0, 2^52).")
    if succP < 0.0 or not succP <= 1.0:
        raise RuntimeError("succP must be in [0, 1].")
    cdef double qd
    try:
        qd = targetP_obj
    except TypeError:
        return qbinom_vectorize_targetp(targetP_obj, n, succP, invert, logTarget)
    return qbinom_internal(qd, n, succP, invert, logTarget)

cdef qbinom_vectorize_all(object targetP_obj, object n_obj, object succP_obj, bint invert, bint logTarget):
    qa = np.asarray(targetP_obj, dtype=np.float64)
    na = np.asarray(n_obj, dtype=np.int64)
    pa = np.asarray(succP_obj, dtype=np.float64)
    it = np.nditer([qa, na, pa], flags=['c_index'])
    cdef cnp.ndarray[cnp.int64_t,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.int64)
    cdef double qd
    cdef int64_t ni
    cdef double pd
    for qd, ni, pd in it:
        if ni < 0 or ni >= (1LL << 52):
            raise RuntimeError("n must be in [0, 2^52).")
        if pd < 0.0 or not pd <= 1.0:
            raise RuntimeError("succP must be in [0, 1].")
        results[it.index] = qbinom_internal(qd, ni, pd, invert, logTarget)
    return np.reshape(results, np.broadcast_shapes(qa.shape, na.shape, pa.shape))

cdef qbinom_vv_internal(object targetP_obj, object n_obj, object succP_obj, bint invert, bint logTarget):
    cdef int64_t ni
    cdef double pd
    try:
        ni = n_obj
        pd = succP_obj
    except TypeError:
        return qbinom_vectorize_all(targetP_obj, n_obj, succP_obj, invert, logTarget)
    return qbinom_v_internal(targetP_obj, ni, pd, invert, logTarget)

# Returns smallest nonnegative k for which cdf(k) >= targetP if logTarget is
# is False, and cdf(k) >= exp(targetP) if logTarget is True.
#
# Implementation is *not* built on top of pbinom() in a way that e.g.
# guarantees qbinom(pbinom(k, n, succP), n, succP) == k or
# qbinom(pbinom(k, n, succP) * (1 + 0.5**52), n, succP) > k in non-degenerate
# cases.  However, it is designed to make these outcomes very likely:
# - Qbinom() is designed for <0.6 ULP relative error and achieves <0.5 ULP the
#   vast majority of the time.
# - The internal Qbinom() call is made with 0.5 ULP subtracted off of q.
def qbinom(object targetP, int64_t n, double succP=0.5, bint logTarget=0):
    return qbinom_v_internal(targetP, n, succP, False, logTarget)


# scipy-style interfaces.  Straightforward to fill in the missing methods (e.g.
# .stats()) if it matters.
class _BinomDist:
    @staticmethod
    def cdf(object k, object n, object p=0.5, bint approx=False):
        return pbinom_vv_internal(k, n, p, complement=False, logp=False, approx=approx)

    @staticmethod
    def isf(object q, object n, object p=0.5):
        return qbinom_vv_internal(q, n, p, invert=True, logTarget=False)

    @staticmethod
    def logcdf(object k, object n, object p=0.5, bint approx=False):
        return pbinom_vv_internal(k, n, p, complement=False, logp=True, approx=approx)

    @staticmethod
    def logpmf(object k, object n, object p=0.5):
        return dbinom_vv_internal(k, n, p, logp=True)

    @staticmethod
    def logsf(object k, object n, object p=0.5, bint approx=False):
        return pbinom_vv_internal(k, n, p, complement=True, logp=True, approx=approx)

    @staticmethod
    def median(object n, object p):
        # silly to have a p=0.5 default here
        return qbinom_vv_internal(0.5, n, p, invert=False, logTarget=False)

    @staticmethod
    def pmf(object k, object n, object p=0.5):
        return dbinom_vv_internal(k, n, p, logp=False)

    @staticmethod
    def ppf(object q, object n, object p=0.5, bint logQ=False):
        return qbinom_vv_internal(q, n, p, invert=False, logTarget=logQ)

    @staticmethod
    def sf(object k, object n, object p=0.5, bint approx=False):
        return pbinom_vv_internal(k, n, p, complement=True, logp=True, approx=approx)

binom = _BinomDist()

cdef class Binomial:
    cdef object n
    cdef object p
    cdef bint _approx
    # further optimizations possible: can special-case scalar n and p, cache
    # lfact_n_tdr / lnp_tdr / lnq_tdr, call n.ravel() and p.ravel() during
    # initialization and split them into cython memoryview + shape-tuple, etc.

    def __cinit__(self, n, p, **kwargs):
        # Straightforward (though a bit tedious) to improve 'tol' support by
        # adding parallel functions that opportunistically stop using
        # dd_real/td_real.
        # Validation seems cheap enough that 'validation_policy' probably isn't
        # worth implementing.
        # Could suppory 'cache_policy' by imitating scipy caching behavior for
        # the slower functions.
        if 'tol' in kwargs:
            # possible todo: prove an upper bound for PbinomApprox()'s relative
            # error and make that the threshold; 2^{-32} is just a guess based
            # on utils/pbinom_accuracy.py results I've seen
            self._approx = (kwargs['tol'] >= (0.5 ** 32))
        else:
            self._approx = False
        self.n = np.asarray(n, dtype=np.int64)
        self.p = np.asarray(p, dtype=np.float64)

    def ccdf(self, k):
        return pbinom_vv_internal(k, self.n, self.p, complement=True, logp=False, approx=self._approx)

    def cdf(self, k):
        return pbinom_vv_internal(k, self.n, self.p, complement=False, logp=False, approx=self._approx)

    def iccdf(self, q):
        return qbinom_vv_internal(q, self.n, self.p, invert=True, logTarget=False)

    def icdf(self, q):
        return qbinom_vv_internal(q, self.n, self.p, invert=False, logTarget=False)

    def ilogccdf(self, q):
        return qbinom_vv_internal(q, self.n, self.p, invert=True, logTarget=True)

    def ilogcdf(self, q):
        return qbinom_vv_internal(q, self.n, self.p, invert=False, logTarget=True)

    def logccdf(self, k):
        return pbinom_vv_internal(k, self.n, self.p, complement=True, logp=True, approx=self._approx)

    def logcdf(self, k):
        return pbinom_vv_internal(k, self.n, self.p, complement=False, logp=True, approx=self._approx)

    def logpmf(self, k):
        return dbinom_vv_internal(k, self.n, self.p, logp=True)

    def median(self):
        return qbinom_vv_internal(0.5, self.n, self.p, invert=False, logTarget=False)

    def pmf(self, k):
        return dbinom_vv_internal(k, self.n, self.p, logp=False)


# doesn't check whether p is in [0, 1]
cdef double binomtest_internal(int64_t k, int64_t n, td_real_struct p_tdr, bint twosided, bint complement, bint midp, bint logp) except? 2.0:
    if k < 0 or k > n:
        raise RuntimeError("k must be nonnegative and <= n.")
    if n >= (1LL << 52):
        raise RuntimeError("n must be less than 2^52.")
    if complement and (not midp):
        k -= 1
    cdef double result
    if twosided:
        result = BinomTwoSidedP(k, n, p_tdr, midp, logp)
    else:
        result = PbinomApprox(k, n, p_tdr, complement, midp, logp)
    return flush_if_denormal(result)

cdef binomtest_vectorize_all(object k_obj, object n_obj, object pa, bint twosided, bint complement, bint midp, bint logp):
    ka = np.asarray(k_obj, dtype=np.int64)
    na = np.asarray(n_obj, dtype=np.int64)
    it = np.nditer([ka, na, pa], flags=['c_index', 'refs_ok'])
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.float64)
    cdef int64_t ki
    cdef int64_t ni
    cdef object poa
    cdef td_real_struct p_tdr
    for ki, ni, poa in it:
        p_tdr = TdrMake(poa.item(0))
        if tdr_ltd(p_tdr, 0) or not tdr_leqd(p_tdr, 1):
            raise RuntimeError("p must be in [0, 1].")
        results[it.index] = binomtest_internal(ki, ni, p_tdr, twosided, complement, midp, logp)
    return np.reshape(results, np.broadcast_shapes(ka.shape, na.shape, pa.shape))

cdef binomtest_vectorize_kn(object k_obj, object n_obj, td_real_struct p_tdr, bint twosided, bint complement, bint midp, bint logp):
    ka = np.asarray(k_obj, dtype=np.int64)
    na = np.asarray(n_obj, dtype=np.int64)
    it = np.nditer([ka, na], flags=['c_index'])
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.float64)
    cdef int64_t ki
    cdef int64_t ni
    for ki, ni in it:
        results[it.index] = binomtest_internal(ki, ni, p_tdr, twosided, complement, midp, logp)
    return np.reshape(results, np.broadcast_shapes(ka.shape, na.shape))

# n must be in [0, 2^52), k must be in [0, n], p must be in [0, 1].
#
# If p is a fractions.Fraction(), it's expanded to a "triple-double" with ~159
# bit accuracy.  This option lets us justify use of a very small tie-detection
# epsilon: we aren't effectively forced to treat p=0.30000000000000004 as if
# the user may have intended p=3/10, because the user has a proper way to
# specify the latter when they care about the difference.  Which, in turn, is
# actually a hard prerequisite for handling some large cases accurately: see
# e.g. the n=21000000 test cases where it has been taken for granted for
# decades that R sometimes has relative error >1e-4.
cpdef binomtest(object k, object n, object p=0.5, str alternative="two-sided", bint midp=0, bint logp=0):
    cdef bint twosided = (alternative == "two-sided")
    cdef bint complement = (alternative == "greater")
    if alternative != "less" and not (twosided or complement):
        raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
    cdef td_real_struct p_tdr
    if isinstance(p, (float, np.float32, np.float64)):
        p_tdr = tdr_make1(p)
    else:
        pa = np.asarray(p)
        if pa.size != 1:
            return binomtest_vectorize_all(k, n, pa, twosided, complement, midp, logp)
        p_tdr = TdrMake(pa.item(0))
    if tdr_ltd(p_tdr, 0) or not tdr_leqd(p_tdr, 1):
        raise RuntimeError("p must be in [0, 1].")
    cdef int64_t ki
    cdef int64_t ni
    try:
        ki = k
        ni = n
    except TypeError:
        return binomtest_vectorize_kn(k, n, p_tdr, twosided, complement, midp, logp)
    return binomtest_internal(ki, ni, p_tdr, twosided, complement, midp, logp)


# scipy and R don't use the same parameters for the hypergeometric
# distribution:
#   R (m, n, k) = (a+c, b+d, a+b)
#   scipy (M, n, N) = (a+b+c+d, a+b, a+c)
#                   = (m+n, k, m) from R
# After trying a few possibilities, I think the least-bad internal
# parameterization is scipy's, since that naturally supports R-entry-point
# vectorization while the reverse isn't true.
cdef double hypergeom_pmf_internal(int64_t k, int64_t M, int64_t n, int64_t N, bint logp) except? 2.0:
    if k < 0 or k >= (1LL << 52):
        raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
    cdef int64_t b = n - k
    cdef int64_t c = N - k
    cdef int64_t d = M - n - N + k
    if b < 0 or c < 0 or d < 0:
        raise RuntimeError("Observation is outside of distribution support.")
    return flush_if_denormal(HypergeomMass(k, b, c, d, logp))

cdef hypergeom_pmf_vectorize_k(object k_obj, int64_t M, int64_t n, int64_t N, bint logp):
    ka = np.asarray(k_obj, dtype=np.int64)
    cdef int64_t [::1] kar = ka.ravel()
    cdef uintptr_t kar_size = kar.size
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(kar_size, dtype=np.float64)
    cdef td_real_struct lfact_m1x_tdr
    cdef td_real_struct lfact_m2x_tdr
    cdef td_real_struct lfact_mx1_tdr
    cdef td_real_struct lfact_mx2_tdr
    cdef td_real_struct lfact_mxx_tdr
    cdef uintptr_t idx
    cdef int64_t ki
    with nogil:
        HypergeomMassMultiKPrecomp(M, n, N, &lfact_m1x_tdr, &lfact_m2x_tdr, &lfact_mx1_tdr, &lfact_mx2_tdr, &lfact_mxx_tdr)
        for idx in range(kar_size):
            ki = kar[idx]
            results[idx] = HypergeomMassJustK(ki, M, n, N, lfact_m1x_tdr, lfact_m2x_tdr, lfact_mx1_tdr, lfact_mx2_tdr, lfact_mxx_tdr, logp)
    return np.reshape(results, ka.shape)

cdef hypergeom_pmf_v_internal(object k_obj, int64_t M, int64_t n, int64_t N, bint logp):
    if M < 0 or n < 0 or N < 0 or M >= (1LL << 52) or n > M or N > M:
        raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
    cdef int64_t ki
    try:
        ki = k_obj
    except TypeError:
        return hypergeom_pmf_vectorize_k(k_obj, M, n, N, logp)
    return hypergeom_pmf_internal(ki, M, n, N, logp)

cdef hypergeom_pmf_vectorize_all(object k_obj, object M_obj, object n_obj, object N_obj, bint logp):
    ka = np.asarray(k_obj, dtype=np.int64)
    Ma = np.asarray(M_obj, dtype=np.int64)
    na = np.asarray(n_obj, dtype=np.int64)
    Na = np.asarray(N_obj, dtype=np.int64)
    it = np.nditer([ka, Ma, na, Na], flags=['c_index'])
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.float64)
    cdef int64_t ki
    cdef int64_t Mi
    cdef int64_t ni
    cdef int64_t Ni
    for ki, Mi, ni, Ni in it:
        if Mi < 0 or ni < 0 or Ni < 0 or Mi >= (1LL << 52) or ni > Mi or Ni > Mi:
            raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
        results[it.index] = hypergeom_pmf_internal(ki, Mi, ni, Ni, logp)
    return np.reshape(results, np.broadcast_shapes(ka.shape, Ma.shape, na.shape, Na.shape))

cdef hypergeom_pmf_vv_internal(object k_obj, object M_obj, object n_obj, object N_obj, bint logp):
    cdef int64_t Mi
    cdef int64_t ni
    cdef int64_t Ni
    try:
        Mi = M_obj
        ni = n_obj
        Ni = N_obj
    except TypeError:
        return hypergeom_pmf_vectorize_all(k_obj, M_obj, n_obj, N_obj, logp)
    return hypergeom_pmf_v_internal(k_obj, Mi, ni, Ni, logp)

def dhyper(object x, int64_t m, int64_t n, int64_t k, bint logp=0):
    return hypergeom_pmf_v_internal(x, m+n, k, m, logp)


cdef double hypergeom_cdf_internal(int64_t k, int64_t M, int64_t n, int64_t N, bint lowertail, bint logp, bint approx) except? 2.0:
    if k < 0 or k >= (1LL << 52):
        # Unlike pbinom(), we don't bother with returning -INFINITY/0/1 in some
        # of these cases.
        raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
    cdef int64_t b = n - k
    cdef int64_t c = N - k
    cdef int64_t d = M - n - N + k
    if b < 0 or c < 0 or d < 0:
        raise RuntimeError("Observation is outside of distribution support.")
    if not lowertail:
        k, b = b, k
        c, d = d, c
        if k == 0 or d == 0:
            if logp:
                return -INFINITY
            return 0
        k -= 1
        b += 1
        c += 1
        d -= 1
    cdef double result
    if approx:
        result = PhyperApprox(k, b, c, d, 0, 0, logp)
    else:
        result = Phyper(k, b, c, d, logp)
    return flush_if_denormal(result)

cdef hypergeom_cdf_vectorize_k(object k_obj, int64_t M, int64_t n, int64_t N, bint lowertail, bint logp, bint approx):
    ka = np.asarray(k_obj, dtype=np.int64)
    cdef int64_t [::1] kar = ka.ravel()
    cdef uintptr_t kar_size = kar.size
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(kar_size, dtype=np.float64)
    cdef uintptr_t idx
    cdef int64_t ki
    for idx in range(kar_size):
        ki = kar[idx]
        results[idx] = hypergeom_cdf_internal(ki, M, n, N, lowertail, logp, approx)
    return np.reshape(results, ka.shape)

cdef hypergeom_cdf_v_internal(object k_obj, int64_t M, int64_t n, int64_t N, bint lowertail, bint logp, bint approx):
    if M < 0 or n < 0 or N < 0 or M >= (1LL << 52) or n > M or N > M:
        raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
    cdef int64_t ki
    try:
        ki = k_obj
    except TypeError:
        return hypergeom_cdf_vectorize_k(k_obj, M, n, N, lowertail, logp, approx)
    return hypergeom_cdf_internal(ki, M, n, N, lowertail, logp, approx)

cdef hypergeom_cdf_vectorize_all(object k_obj, object M_obj, object n_obj, object N_obj, bint lowertail, bint logp, bint approx):
    ka = np.asarray(k_obj, dtype=np.int64)
    Ma = np.asarray(M_obj, dtype=np.int64)
    na = np.asarray(n_obj, dtype=np.int64)
    Na = np.asarray(N_obj, dtype=np.int64)
    it = np.nditer([ka, Ma, na, Na], flags=['c_index'])
    cdef cnp.ndarray[double,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.float64)
    cdef int64_t ki
    cdef int64_t Mi
    cdef int64_t ni
    cdef int64_t Ni
    for ki, Mi, ni, Ni in it:
        if Mi < 0 or ni < 0 or Ni < 0 or Mi >= (1LL << 52) or ni > Mi or Ni > Mi:
            raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
        results[it.index] = hypergeom_cdf_internal(ki, Mi, ni, Ni, lowertail, logp, approx)
    return np.reshape(results, np.broadcast_shapes(ka.shape, Ma.shape, na.shape, Na.shape))

cdef hypergeom_cdf_vv_internal(object k_obj, object M_obj, object n_obj, object N_obj, bint lowertail, bint logp, bint approx):
    cdef int64_t Mi
    cdef int64_t ni
    cdef int64_t Ni
    try:
        Mi = M_obj
        ni = n_obj
        Ni = N_obj
    except TypeError:
        return hypergeom_cdf_vectorize_all(k_obj, M_obj, n_obj, N_obj, lowertail, logp, approx)
    return hypergeom_cdf_v_internal(k_obj, Mi, ni, Ni, lowertail, logp, approx)

def phyper(object x, int64_t m, int64_t n, int64_t k, bint lowertail=1, bint logp=0, bint approx=0):
    return hypergeom_cdf_v_internal(x, m+n, k, m, lowertail, logp, approx)


cdef int64_t hypergeom_ppf_internal(double p, int64_t M, int64_t n, int64_t N, bint invert, bint logp) except? -1:
    cdef dd_real_struct p_ddr
    if invert:
        if logp:
            p_ddr = ddr_addd(ddr_negate(ddr_exp(ddr_maked(p))), 1.0)
            logp = False
        else:
            p_ddr = ddr_add2d(1.0, -p)
    else:
        p_ddr = ddr_maked(p)
    if logp:
        if not ddr_leqd(p_ddr, 0.0):
            raise RuntimeError("p must be <= 0 when logp is True.")
    else:
        if ddr_ltd(p_ddr, 0.0) or not ddr_leqd(p_ddr, 1.0):
            raise RuntimeError("p must be in [0, 1] when logp is False.")
    return QhyperHalfUlp(p_ddr, N, M-N, n, logp)

cdef hypergeom_ppf_vectorize_p(object p_obj, int64_t M, int64_t n, int64_t N, bint invert, bint logp):
    pa = np.asarray(p_obj, dtype=np.float64)
    cdef double [::1] par = pa.ravel()
    cdef uintptr_t par_size = par.size
    cdef cnp.ndarray[cnp.int64_t,mode="c",ndim=1] results = np.empty(par_size, dtype=np.int64)
    cdef uintptr_t idx
    cdef double pd
    for idx in range(par_size):
        pd = par[idx]
        results[idx] = hypergeom_ppf_internal(pd, M, n, N, invert, logp)
    return np.reshape(results, pa.shape)

cdef hypergeom_ppf_v_internal(object p_obj, int64_t M, int64_t n, int64_t N, bint invert, bint logp):
    if M < 0 or n < 0 or N < 0 or M >= (1LL << 52) or n > M or N > M:
        raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
    cdef double pd
    try:
        pd = p_obj
    except TypeError:
        return hypergeom_ppf_vectorize_p(p_obj, M, n, N, invert, logp)
    return hypergeom_ppf_internal(pd, M, n, N, invert, logp)

cdef hypergeom_ppf_vectorize_all(object p_obj, object M_obj, object n_obj, object N_obj, bint invert, bint logp):
    pa = np.asarray(p_obj, dtype=np.float64)
    Ma = np.asarray(M_obj, dtype=np.int64)
    na = np.asarray(n_obj, dtype=np.int64)
    Na = np.asarray(N_obj, dtype=np.int64)
    it = np.nditer([pa, Ma, na, Na], flags=['c_index'])
    cdef cnp.ndarray[cnp.int64_t,mode="c",ndim=1] results = np.empty(it.itersize, dtype=np.int64)
    cdef double pd
    cdef int64_t Mi
    cdef int64_t ni
    cdef int64_t Ni
    for pd, Mi, ni, Ni in it:
        if Mi < 0 or ni < 0 or Ni < 0 or Mi >= (1LL << 52) or ni > Mi or Ni > Mi:
            raise RuntimeError("Parameters, row/column sums, and population size must be in [0, 2^52).")
        results[it.index] = hypergeom_ppf_internal(pd, Mi, ni, Ni, invert, logp)
    return np.reshape(results, np.broadcast_shapes(pa.shape, Ma.shape, na.shape, Na.shape))

cdef hypergeom_ppf_vv_internal(object p_obj, object M_obj, object n_obj, object N_obj, bint invert, bint logp):
    cdef int64_t Mi
    cdef int64_t ni
    cdef int64_t Ni
    try:
        Mi = M_obj
        ni = n_obj
        Ni = N_obj
    except TypeError:
        return hypergeom_ppf_vectorize_all(p_obj, M_obj, n_obj, N_obj, invert, logp)
    return hypergeom_ppf_v_internal(p_obj, Mi, ni, Ni, invert, logp)

# Returns smallest x in the distribution support for which cdf(x) >= p.
#
# Implementation is *not* built on top of phyper() in a way that e.g.
# guarantees qhyper(phyper(x, m, n, k), m, n, k) == x or
# qhyper(phyper(x, m, n, k) * (1 + 0.5**52), m, n, k) > a in non-degenerate
# cases.  However, it is designed to make these outcomes very likely:
# - Phyper() is designed for <0.6 ULP relative error (except when n is huge),
#   and achieves <0.5 ULP the vast majority of the time.
# - The internal Qhyper() call is made with 0.5 ULP subtracted off of q.
def qhyper(object p, int64_t m, int64_t n, int64_t k, bint logp=0):
    return hypergeom_ppf_v_internal(p, m+n, k, m, False, logp)


# scipy-style interface.  Straightforward to fill in the missing methods (e.g.
# .stats()) if it matters.
class _HypergeomDist:
    @staticmethod
    def cdf(object k, object M, object n, object N, bint approx=False):
        return hypergeom_cdf_vv_internal(k, M, n, N, lowertail=True, logp=False, approx=approx)

    @staticmethod
    def isf(object q, object M, object n, object N):
        return hypergeom_ppf_vv_internal(q, M, n, N, invert=True, logp=False)

    @staticmethod
    def logcdf(object k, object M, object n, object N, bint approx=False):
        return hypergeom_cdf_vv_internal(k, M, n, N, lowertail=True, logp=True, approx=approx)

    @staticmethod
    def logpmf(object k, object M, object n, object N):
        return hypergeom_pmf_vv_internal(k, M, n, N, logp=True)

    @staticmethod
    def logsf(object k, object M, object n, object N, bint approx=False):
        return hypergeom_cdf_vv_internal(k, M, n, N, lowertail=False, logp=True, approx=approx)

    @staticmethod
    def median(object M, object n, object N):
        return hypergeom_ppf_vv_internal(0.5, M, n, N, invert=False, logp=False)

    @staticmethod
    def pmf(object k, object M, object n, object N):
        return hypergeom_pmf_vv_internal(k, M, n, N, logp=False)

    @staticmethod
    def ppf(object q, object M, object n, object N):
        return hypergeom_ppf_vv_internal(q, M, n, N, invert=False, logp=False)

    @staticmethod
    def sf(object k, object M, object n, object N, bint approx=False):
        return hypergeom_cdf_vv_internal(k, M, n, N, lowertail=False, logp=False, approx=approx)

hypergeom = _HypergeomDist()


# table must be a 2x2 or larger matrix, represented as a list of equal-length
# lists.  For 2x2 tests, values must be nonnegative integers which add up
# to <2^52.  For tests on larger tables, they must add up to <2^31.
#
# alternative must be one of the following:
#   "two-sided": default, must be this if table is larger than 2x2.
#   "less": alt hypothesis is that table[0][0] is smaller than expected.
#   "greater": alt hypothesis is that table[0][0] is larger than expected.
def fisher_exact(list table, str alternative="two-sided", bint midp=0, bint logp=0):
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
            raise RuntimeError("2x2 table entries must sum to <2^52")
        m11_is_greater_alt = (alternative == "greater")
        if alternative != "less" and not m11_is_greater_alt:
            raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
        return flush_if_denormal(PhyperApprox(m11, m12, m21, m22, m11_is_greater_alt, midp, logp))
    if nrow == 2 and ncol == 2:
        if total >= (1LL << 52):
            raise RuntimeError("2x2 table entries must sum to <2^52")
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
            raise RuntimeError(f"{nrow}x{ncol} table entries must be nonnegative and sum to <2^31")
        with nogil:
            ln_result = Fisher23LnP(m11, m12, m13, m21, m22, m23, midp)
    else:
         raise RuntimeError("tables larger than 2x3 not yet supported")
    if logp:
        return ln_result
    return exp_flush(ln_result)


# Point estimate provided by R fisher.test.
def cond_odds_ratio(int64_t m11, int64_t m12, int64_t m21, int64_t m22):
    if m11 < 0 or m12 < 0 or m21 < 0 or m22 < 0:
        raise RuntimeError("table entries must be nonnegative")
    if m11 >= (1LL << 52) or m12 >= (1LL << 52) or m21 >= (1LL << 52) or m22 >= (1LL << 52):
        raise RuntimeError("table entries must be <2^52")
    cdef int64_t total = m11 + m12 + m21 + m22
    if total >= (1LL << 52):
        raise RuntimeError("table entries must sum to <2^52")
    return Fisher22OddsRatio(m11, m12, m21, m22)


# CI provided by R fisher.test.
def cond_odds_ratio_ci(int64_t m11, int64_t m12, int64_t m21, int64_t m22, double low=0.025, double high=0.975):
    if m11 < 0 or m12 < 0 or m21 < 0 or m22 < 0:
        raise RuntimeError("table entries must be nonnegative")
    if m11 >= (1LL << 52) or m12 >= (1LL << 52) or m21 >= (1LL << 52) or m22 >= (1LL << 52):
        raise RuntimeError("table entries must be <2^52")
    cdef int64_t total = m11 + m12 + m21 + m22
    if total >= (1LL << 52):
        raise RuntimeError("table entries must sum to <2^52")
    if low > high:
        raise RuntimeError("low can't be greater than high")
    if ((low < 0.5**24) and (low != 0.0)) or ((low > 1 - 0.5**24) and (low != 1.0)) or ((high < 0.5**24) and (high != 0.0)) or ((high > 1 - 0.5**24) and (high != 1.0)):
        raise RuntimeError("low and high must be 0, 1, or in [2^{-24}, 1 - 2^{-24}].")
    cdef double low_result
    cdef double high_result
    Fisher22OddsRatioCI(m11, m12, m21, m22, low, high, &low_result, &high_result)
    return (low_result, high_result)


cdef double HWE_exact_2sided_internal(int32_t hom1, int32_t hets, int32_t hom2, bint midp, bint logp) except? 2.0:
    cdef int64_t total = <int64_t>(hom1) + <int64_t>(hets) + <int64_t>(hom2)
    if hom1 < 0 or hets < 0 or hom2 < 0 or total > 0x7fffffff:
        raise RuntimeError("hom1, hets and hom2 must be nonnegative and sum to <2^31")
    cdef bint hets_is_greater_alt = 0
    cdef double ln_result = HweLnP(hets, hom1, hom2, midp)
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
def HWE_exact(int32_t hom1, int32_t hets, int32_t hom2, str alternative="two-sided", bint midp=0, bint logp=0):
    if alternative == "two-sided":
        return HWE_exact_2sided_internal(hom1, hets, hom2, midp, logp)
    cdef bint hets_is_greater_alt = (alternative == "greater")
    if alternative != "less" and not hets_is_greater_alt:
        raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
    raise RuntimeError("one-sided tests not implemented yet")

def snphwe(int32_t hets, int32_t hom1, int32_t hom2, bint midp=0, bint logp=0):
    return HWE_exact_2sided_internal(hom1, hets, hom2, midp, logp)
