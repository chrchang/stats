#ifndef __FISHER_NETWORK_H__
#define __FISHER_NETWORK_H__

// Fisher's Exact Test library, copyright (C) 2013-2026 Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Internal interface to the network algorithm used by Fisher3xkLnPEx().

#include "fisher.h"

#ifdef __cplusplus
namespace plink2 {
#endif

enum {
  kFisherNetOk = 0,
  // The tracked workspace would exceed mem_limit_bytes (or an allocation
  // failed).
  kFisherNetMemLimit = 1,
  // More than 8 rows and columns, a collision between hashed node keys (only
  // used when the row totals don't fit in a 62-bit key), a detected hash
  // collision or unresolvable near-tie.
  kFisherNetUnsupported = 2,
  // check_fn returned nonzero.
  kFisherNetCanceled = 3
};

// Natural log of the two-sided Fisher's exact test [mid-]p-value for a
// row_ct x col_ct table (obs[] in row-major order; rows and columns may have
// zero totals).  Tables count as at least as extreme when their probability is
// <= the observed table's; ties are resolved exactly, and count half with
// midp.  Returns NaN with *status_out != kFisherNetOk on failure.
// mem_limit_bytes bounds the main workspace (see Fisher3xkLnPEx()); 0 means no
// limit, which is dangerous (see the warning there).  check_fn may be null.
// Thread-safe: each call has its own workspace (so concurrent calls can use up
// to mem_limit_bytes each).
double FisherNetworkLnP(const int32_t* obs, uint32_t row_ct, uint32_t col_ct, uint32_t midp, uint64_t mem_limit_bytes, FisherCheckFn check_fn, void* check_arg, uint32_t* status_out);

#ifdef __cplusplus
}
#endif

#endif  // __FISHER_NETWORK_H__
