#include <stdint.h>
#include <inttypes.h>

// for INFINITY
#include <math.h>

#define EPSILON 0.0000000000001

// A bias of this sort is needed to enable calculation of p-values down to the
// minimum representable positive number.
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

double binom_2sided(uint32_t succ, uint32_t obs, double rate, uint32_t midp) {
  // straightforward to generalize this to any success probability
  double cur_succ_t2 = (int32_t)succ;
  double cur_fail_t2 = (int32_t)(obs - succ);
  double tailp = (1 - EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp;
  double lastp1 = tailp;
  int32_t tie_ct = 1;
  double rate_mult_incr = rate / (1 - rate);
  double rate_mult_decr = (1 - rate) / rate;
  double cur_succ_t1;
  double cur_fail_t1;
  double preaddp;
  if (!obs) {
    return midp? 0.5 : 1;
  }
  if (obs * rate < succ) {
    while (cur_succ_t2 > 0.5) {
      cur_fail_t2 += 1;
      lastp2 *= rate_mult_decr * cur_succ_t2 / cur_fail_t2;
      cur_succ_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * EPSILON) * EXACT_TEST_BIAS) {
	  tie_ct++;
	}
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (cur_succ_t2 > 0.5) {
      cur_fail_t2 += 1;
      lastp2 *= rate_mult_decr * cur_succ_t2 / cur_fail_t2;
      cur_succ_t2 -= 1;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    cur_succ_t1 = (int32_t)(succ + 1);
    cur_fail_t1 = (int32_t)(obs - succ);
    while (cur_fail_t1 > 0.5) {
      lastp1 *= rate_mult_incr * cur_fail_t1 / cur_succ_t1;
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      cur_succ_t1 += 1;
      cur_fail_t1 -= 1;
    }
  } else {
    while (cur_fail_t2 > 0.5) {
      cur_succ_t2++;
      lastp2 *= rate_mult_incr * cur_fail_t2 / cur_succ_t2;
      cur_fail_t2--;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * EPSILON) * EXACT_TEST_BIAS) {
	  tie_ct++;
	}
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (cur_fail_t2 > 0.5) {
      cur_succ_t2 += 1;
      lastp2 *= rate_mult_incr * cur_fail_t2 / cur_succ_t2;
      cur_fail_t2 -= 1;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    cur_succ_t1 = (int32_t)succ;
    cur_fail_t1 = (int32_t)(obs - succ);
    while (cur_succ_t1 > 0.5) {
      cur_fail_t1 += 1;
      lastp1 *= rate_mult_decr * cur_succ_t1 / cur_fail_t1;
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      cur_succ_t1 -= 1;
    }
  }
  if (!midp) {
    return tailp / (tailp + centerp);
  } else {
    return (tailp - ((1 - EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (tailp + centerp);
  }
}

double binom_1sided(int32_t succ, int32_t obs, double rate) {
  double cur_prob = EXACT_TEST_BIAS;
  double left_prob = cur_prob;
  double right_prob = 0;
  double rate_mult_incr = rate / (1 - rate);
  double rate_mult_decr = (1 - rate) / rate;
  double cur_succ = succ;
  double cur_fail = obs - succ;
  double preaddp;
  if (obs * rate < succ) {
    while (cur_succ > 0.5) {
      cur_fail += 1;
      cur_prob *= rate_mult_decr * cur_succ / cur_fail;
      cur_succ -= 1;
      preaddp = left_prob;
      left_prob += cur_prob;
      if (left_prob <= preaddp) {
	break;
      }
      if (left_prob >= 1.0) {
	return 1;
      }
    }
    cur_succ = succ;
    cur_fail = obs - succ;
    cur_prob = EXACT_TEST_BIAS;
    right_prob = left_prob;
    while (cur_fail > 0.5) {
      cur_succ += 1;
      cur_prob *= rate_mult_incr * cur_fail / cur_succ;
      cur_fail -= 1;
      preaddp = right_prob;
      right_prob += cur_prob;
      if (right_prob <= preaddp) {
	break;
      }
    }
    return left_prob / right_prob;
  } else {
    while (cur_fail > 0.5) {
      cur_succ += 1;
      cur_prob *= rate_mult_incr * cur_fail / cur_succ;
      cur_fail -= 1;
      preaddp = right_prob;
      right_prob += cur_prob;
      if (right_prob == INFINITY) {
	return 0;
      }
      if (right_prob <= preaddp) {
	break;
      }
    }
    cur_succ = succ;
    cur_fail = obs - succ;
    cur_prob = EXACT_TEST_BIAS;
    while (cur_succ > 0.5) {
      cur_fail += 1;
      cur_prob *= rate_mult_decr * cur_succ / cur_fail;
      cur_succ -= 1;
      preaddp = left_prob;
      left_prob += cur_prob;
      if (left_prob <= preaddp) {
	break;
      }
    }
    return left_prob / (left_prob + right_prob);
  }
}
