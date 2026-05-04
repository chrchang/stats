# cython: language_level=3

from libc.stdint cimport uint32_t, int32_t
from libc.math cimport exp

__version__ = "0.1.0"

cdef extern from "../include/fisher.h" namespace "plink2":
    ctypedef uint32_t BoolErr

    BoolErr Fisher22LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, double* resultp)

    double Fisher22OneSidedLnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp)

    BoolErr Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp, double* resultp)


# table must be a 2x2 or larger matrix, represented as a list of equal-length
# lists.  Values must be nonnegative integers which add up to <2^31.
#
# alternative must be one of the following:
#   "two-sided": default, must be this if table is larger than 2x2.
#   "less": alt hypothesis is that table[0][0] is smaller than expected.
#   "greater": alt hypothesis is that table[0][0] is larger than expected.
cpdef double fisher(table, str alternative = "two-sided", bint midp = 0, bint logp = 0):
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
    cdef bint m11_is_greater_alt = 0
    cdef int32_t m13_or_31
    cdef int32_t m23_or_32
    cdef double ln_result
    if alternative != "two-sided":
        if nrow > 2 or ncol > 2:
            raise RuntimeError("alternative must be 'two-sided' for tables larger than 2x2.")
        m11_is_greater_alt = (alternative == "greater")
        if alternative != "less" and not m11_is_greater_alt:
            raise RuntimeError("alternative is not in {'two-sided', 'less', 'greater'}.")
        ln_result = Fisher22OneSidedLnP(m11, m12, m21, m22, m11_is_greater_alt, midp)
    else:
        if nrow == 2:
            if ncol == 2:
                if Fisher22LnP(m11, m12, m21, m22, midp, &ln_result):
                    raise MemoryError()
            else:
                if ncol > 3:
                    raise RuntimeError("tables larger than 2x3 not yet supported")
                m13_or_31 = table[0][2]
                m23_or_32 = table[1][2]
                if Fisher23LnP(m11, m12, m13_or_31, m21, m22, m23_or_32, midp, &ln_result):
                    raise MemoryError()
        else:
            if nrow > 3 or ncol > 2:
                raise RuntimeError("tables larger than 2x3 not yet supported")
            m13_or_31 = table[2][0]
            m23_or_32 = table[2][1]
            if Fisher23LnP(m11, m21, m13_or_31, m12, m22, m23_or_32, midp, &ln_result):
                raise MemoryError()
    if logp:
        return ln_result
    return exp(ln_result)
