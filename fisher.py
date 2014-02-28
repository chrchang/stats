#!/usr/bin/env python

# Fisher 2x2 exact test p-value calculator, Python port
# Copyright (C) 2014  Christopher Chang  chrchang@alumni.caltech.edu

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses>.

EXACT_TEST_BIAS = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125
# 2^{-40} for now, since 2^{-44} was too small on real data
FISHER_EPSILON = 0.0000000000009094947017729282379150390625

# Python 2.x does not have native support for comparison to infinity
INFINITY = 1e5000

def fisher22(m11, m12, m21, m22, midp=False):
    # Basic 2x2 Fisher exact test p-value calculation.  m11, m12, m21, and m22
    # must be nonnegative integers smaller than 2^31 (higher-precision
    # arithmetic would be needed to go beyond that).
    # See https://github.com/chrchang/stats for C/C++ extensions of this
    # algorithm to 2x3, Hardy-Weinberg, and binomial tests.
    tprob = (1 - FISHER_EPSILON) * EXACT_TEST_BIAS
    cur_prob = tprob
    cprob = 0.0
    tie_ct = 1
    # Ensure we are left of the distribution center, m11 <= m22, and
    # m12 <= m21.
    if m12 > m21:
        uii = m12
        m12 = m21
        m21 = uii
    if m11 > m22:
        uii = m11
        m11 = m22
        m22 = uii
    if m11 * m22 > m12 * m21:
        uii = m11
        m11 = m12
        m12 = uii
        uii = m21
        m21 = m22
        m22 = uii
    cur11 = m11
    cur12 = m12
    cur21 = m21
    cur22 = m22
    while cur12 > 0.5:
        # must force floating-point division here
        cur11 += 1.0
        cur22 += 1.0
        cur_prob *= (cur12 * cur21) / (cur11 * cur22)
        cur12 -= 1.0
        cur21 -= 1.0
        if cur_prob == INFINITY:
            return 0.0
        if cur_prob < EXACT_TEST_BIAS:
            if cur_prob > (1 - 2 * FISHER_EPSILON) * EXACT_TEST_BIAS:
                tie_ct += 1
            tprob += cur_prob
            break
        cprob += cur_prob
    if cprob == 0.0 and not midp:
        return 1.0
    while cur12 > 0.5:
        cur11 += 1.0
        cur22 += 1.0
        cur_prob *= (cur12 * cur21) / (cur11 * cur22)
        cur12 -= 1.0
        cur21 -= 1.0
        preaddp = tprob
        tprob += cur_prob
        if tprob <= preaddp:
	    break
    if m11 > 0.5:
        cur11 = m11
        cur12 = m12
        cur21 = m21
        cur22 = m22
        cur_prob = (1 - FISHER_EPSILON) * EXACT_TEST_BIAS
        while cur11 > 0.5:
            cur12 += 1.0
            cur21 += 1.0
            cur_prob *= (cur11 * cur22) / (cur12 * cur21)
            cur11 -= 1.0
            cur22 -= 1.0
            preaddp = tprob
            tprob += cur_prob
            if tprob <= preaddp:
                if midp:
                    return (preaddp - ((1 - FISHER_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (cprob + preaddp)
                else:
                    return preaddp / (cprob + preaddp)
    if midp:
        return (tprob - ((1 - FISHER_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (cprob + tprob)
    else:
        return tprob / (cprob + tprob)
