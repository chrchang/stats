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

// Network algorithm for r x c Fisher's exact tests; see fisher_network.h.
//
// Columns are processed one at a time in ascending-total order.  A node is the
// sorted vector of remaining row totals; it holds the list of distinct "past
// values" (log-probability contributions of the columns fixed so far), merged
// on equality of factorial products (an additive hash of prime exponents,
// confirmed by value, so that a hash collision makes the call fall back
// instead of being miscounted).  Each node gets cheap bounds on its most and
// least probable completions (an exact transportation mode is computed lazily
// when the cheap bound is inconclusive), so whole prefixes/suffixes of a
// mother's sorted list are counted in closed form or dropped.  Daughters are
// enumerated in groups (closed-form Vandermonde totals behind an even-split
// upper bound, and an interval search on two-free-row groups, whose best
// completion value is concave), and the last two columns are enumerated
// outward from the mode with ratio recurrences.  Under memory pressure,
// pending merges are flushed early and second-to-last-stage nodes switch to a
// list-free mode.
// Decisions use double sums of log-factorials with a margin; every added weight
// is computed from double-double exponents, and tables inside the margin are
// resolved exactly (factorial-product hash confirmed by value, then
// double-double comparison); what that can't resolve makes the call return
// kFisherNetUnsupported.
//
// Like the rest of this library, this compiles as either C99 or C++.  All of a
// call's memory is owned by its Net workspace and released at the end of
// FisherNetworkLnP().  Failures (memory limit or allocation failure,
// cancellation, an unresolvable tie or detected collision) set Net.err, and
// callers return as soon as they see it.

#include "fisher_network.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#else
#  include <sys/mman.h>
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

// Large buffers are mapped directly from the OS and unmapped when freed.  The
// network's big buffers (hash tables, merge/sort scratch, past-value arenas)
// grow and shrink from stage to stage; through malloc, freed blocks stayed
// resident (macOS kept ~0.5 GiB of empty large-allocation regions), and peak
// resident memory reached ~1.5-2x the tracked workspace with 1 GiB limits.
// These return nullptr on failure.
static const size_t kBigAllocMin = S_CAST(size_t, 1) << 20;
static void* BigAllocRaw(size_t sz) {
  if (sz < kBigAllocMin) {
    return malloc(sz);
  }
#ifdef _WIN32
  return VirtualAlloc(nullptr, sz, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
#else
  void* p = mmap(nullptr, sz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
  return (p == MAP_FAILED)? nullptr : p;
#endif
}
static void BigFreeRaw(void* p, size_t sz) {
  if (!p) {
    return;
  }
  if (sz < kBigAllocMin) {
    free(p);
    return;
  }
#ifdef _WIN32
  VirtualFree(p, 0, MEM_RELEASE);
#else
  munmap(p, sz);
#endif
}
// always OS-mapped (any size)
static void* MapRaw(size_t sz) {
#ifdef _WIN32
  return VirtualAlloc(nullptr, sz, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
#else
  void* p = mmap(nullptr, sz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
  return (p == MAP_FAILED)? nullptr : p;
#endif
}
static void UnmapRaw(void* p, size_t sz) {
#ifdef _WIN32
  (void)sz;
  VirtualFree(p, 0, MEM_RELEASE);
#else
  munmap(p, sz);
#endif
}

// Tuning switches.  Normally compile-time constants (so disabled code paths
// fold away); with FISHER_DEBUG_OVERRIDES they are mutable, set from
// environment variables by ApplyDebugEnv() on every call, and debug builds
// should then not run the network concurrently on several threads.
#ifdef FISHER_DEBUG_OVERRIDES
#  define NET_TUNE static
#elif defined(__cplusplus)
#  define NET_TUNE static constexpr
#else
#  define NET_TUNE static const
#endif
NET_TUNE uint32_t g_bsort = 1;
NET_TUNE size_t g_bsort_min = 2048;
NET_TUNE uint32_t g_bucket = 1;
NET_TUNE int g_bucket_min = 64;
// Second-to-last-stage nodes whose runs average >= g_defer_avg entries are
// deferred (see Core()); g_defer: -1 = by row count, 0 = off, 1 = on.
NET_TUNE int g_defer = -1;
NET_TUNE double g_defer_avg = 16.;
NET_TUNE uint32_t g_gskip = 1;
NET_TUNE uint32_t g_interval = 1;
#ifdef FISHER_DEBUG_OVERRIDES
static size_t g_seg_t = 0;          // sort concatenated segments outright below this many entries per segment
static int g_seg_w = 0;             // weight-balanced segment merge (measured neutral-to-slower)
static uint32_t g_gcheck = 0;       // verify group-skip tails by brute force
static const char* g_order_spec = nullptr;
#endif
// Last-stage tails stop once the remaining terms are below this fraction of
// the tail sum; previous tails are reused for at most this many steps.
static const double kTailEps = 1e-17;
static const int kTailReuseMax = 16;

// Maximum number of rows (the smaller dimension, after transposing).  Per-node
// row totals and the per-column scratch arrays below (x[], dr[], hi[], rem[],
// caps[], md[], S[], sumlf[]) are sized for this.  It's not a correctness
// limit, but costs grow quickly with the number of rows: subproblem counts
// like R^(r-1), column configurations per subproblem like c^(r-1) / (r-1)!,
// and GroupMin() checks j * 2^(j-1) vertices for j free rows.  See
// kFisherRxcMaxMinDim in fisher.h, which must not exceed this.
enum { kFisherNetMaxRows = 8 };
#ifdef __cplusplus
static_assert(kFisherRxcMaxMinDim <= S_CAST(uint32_t, kFisherNetMaxRows), "kFisherRxcMaxMinDim exceeds kFisherNetMaxRows");
#endif

// Double-double values (plink2_highprec.h's dd_real): decisions use the high
// parts; every weight added to the result is computed from double-double
// exponents (ddr_exp_to_double()), so rounding of the large log-factorial
// sums doesn't reach the p-value.  Sums use the non-renormalizing
// ddr_sloppier_add()/ddr_sloppier_sub(); the kTieEpsDD margin accounts for
// that.

// compensated sum
typedef struct AccStruct {
  double s, comp;
} Acc;
static inline void AccAdd(Acc* a, double v) {
  neumaier_add(v, &a->s, &a->comp);
}
static inline double AccGet(const Acc* a) {
  return a->s + a->comp;
}

typedef struct EntryStruct {
  double v;
  uint64_t h;
  double f;
  double vl;  // low part of v
#ifdef __cplusplus
  bool operator<(const struct EntryStruct& rhs) const {
    return v < rhs.v;
  }
#endif
} Entry;

// A node's past-value list: a span in its stage's arena (see Arena).
typedef struct EntryListStruct {
  Entry* p;
  size_t n;
} EntryList;

// deferred leaf: threshold (strictly-less side), weight, completion hash,
// completion value
typedef struct LeafRecStruct {
  double q;
  double w;
  uint64_t h;
  double tv, tvl;
#ifdef __cplusplus
  bool operator<(const struct LeafRecStruct& rhs) const {
    return q < rhs.q;
  }
#endif
} LeafRec;

// A slice [lo, hi) of mother m's list, shifted by the column value colv (+
// colvl) and the column's factorial-product hash ch.
typedef struct RunStruct {
  int m;
  uint32_t lo, hi;
  double colv;
  uint64_t ch;
  double colvl;
} Run;

// MergeRuns() hash slot: slots hold (h, index) so the 92% of inserts that
// merge need no load of out[ix].h.
typedef struct SlotStruct {
  uint64_t h;
  int32_t ix;
} Slot;

// flush order: largest run sets first (unique keys, so the order is fully
// determined)
typedef struct FlushOrdStruct {
  size_t cap, d;
#ifdef __cplusplus
  bool operator<(const struct FlushOrdStruct& rhs) const {
    return (cap > rhs.cap) || ((cap == rhs.cap) && (d > rhs.d));
  }
#endif
} FlushOrd;

// ModeCertify() residual-graph edge
typedef struct CertEdgeStruct {
  int u, v, q, up;
} CertEdge;

#ifndef __cplusplus
static int CmpEntryV(const void* aa, const void* bb) {
  const double x = S_CAST(const Entry*, aa)->v;
  const double y = S_CAST(const Entry*, bb)->v;
  return (x > y) - (x < y);
}
static int CmpLeafQ(const void* aa, const void* bb) {
  const double x = S_CAST(const LeafRec*, aa)->q;
  const double y = S_CAST(const LeafRec*, bb)->q;
  return (x > y) - (x < y);
}
static int CmpFlushOrd(const void* aa, const void* bb) {
  const FlushOrd* a = S_CAST(const FlushOrd*, aa);
  const FlushOrd* b = S_CAST(const FlushOrd*, bb);
  if (a->cap != b->cap) {
    return (a->cap > b->cap)? -1 : 1;
  }
  return (a->d > b->d)? -1 : ((a->d < b->d)? 1 : 0);
}
static int CmpInt(const void* aa, const void* bb) {
  const int x = *S_CAST(const int*, aa);
  const int y = *S_CAST(const int*, bb);
  return (x > y) - (x < y);
}
#endif

typedef struct RunVecStruct {
  Run* p;
  size_t n, cap;
} RunVec;

typedef struct NodeStruct {
  int r[kFisherNetMaxRows];
  EntryList ents;
  uint32_t have_bounds, have_lp;
  // max completion (lazy), lower bound on min completion, log total, cheap
  // upper bound on lp, low part of logt
  double lp, splb, logt, ub, logtl;
  // deferred last stage: per-entry weights / prefix sums (scaled by pmaxv)
  double* wv;
  double* cumv;
  size_t wv_n, cumv_n;
  double pmaxv, pmaxvl;
} Node;

typedef struct NetStruct Net;
static BoolErr NetFail(Net* net, uint32_t code);
static void* NetRegrow(Net* net, void* old, size_t used_bytes, size_t old_bytes, size_t new_bytes);

// Growable array of T with std::vector-like capacity behavior.  Reserve()
// grows geometrically and keeps the contents; Assign() is like assign(n, x)
// or resize(n) on an empty vector: it reallocates exactly n when the capacity
// is too small, and leaves the contents unspecified.  (Not every element type
// uses every function.)
#define NET_VEC(Name, T) \
typedef struct Name##Struct { \
  T* p; \
  size_t n, cap; \
} Name; \
__maybe_unused static inline BoolErr Name##Reserve(Net* net, Name* v, size_t need) { \
  if (need <= v->cap) { \
    return 0; \
  } \
  const size_t nc = MAXV(2 * v->cap, need); \
  void* q = NetRegrow(net, v->p, v->n * sizeof(T), v->cap * sizeof(T), nc * sizeof(T)); \
  if (!q) { \
    return 1; \
  } \
  v->p = S_CAST(T*, q); \
  v->cap = nc; \
  return 0; \
} \
__maybe_unused static inline BoolErr Name##Push(Net* net, Name* v, T e) { \
  if ((v->n == v->cap) && Name##Reserve(net, v, v->n + 1)) { \
    return 1; \
  } \
  v->p[v->n++] = e; \
  return 0; \
} \
__maybe_unused static inline BoolErr Name##Assign(Net* net, Name* v, size_t n) { \
  if (n > v->cap) { \
    BigFreeRaw(v->p, v->cap * sizeof(T)); \
    v->p = nullptr; \
    v->n = 0; \
    v->cap = 0; \
    void* q = BigAllocRaw(n * sizeof(T)); \
    if (!q) { \
      return NetFail(net, kFisherNetMemLimit); \
    } \
    v->p = S_CAST(T*, q); \
    v->cap = n; \
  } \
  v->n = n; \
  return 0; \
} \
__maybe_unused static inline void Name##Free(Name* v) { \
  BigFreeRaw(v->p, v->cap * sizeof(T)); \
  v->p = nullptr; \
  v->n = 0; \
  v->cap = 0; \
}

NET_VEC(DoubleVec, double)
NET_VEC(IntVec, int)
NET_VEC(U32Vec, uint32_t)
NET_VEC(U64Vec, uint64_t)
NET_VEC(EntryVec, Entry)
NET_VEC(LeafVec, LeafRec)
NET_VEC(SlotVec, Slot)
NET_VEC(NodeVec, Node)
NET_VEC(RunVecVec, RunVec)
NET_VEC(FlushOrdVec, FlushOrd)
NET_VEC(CertEdgeVec, CertEdge)

// Per-stage arena for past-value lists.  Lists are packed into 8 MiB blocks
// (a list larger than a quarter block gets a block of its own), and a stage's
// whole arena is freed at once when its lists are no longer needed.
typedef struct ArenaStruct {
  Entry** blocks;
  size_t* sizes;
  size_t nblk, blk_cap;
  Entry* cur;
  size_t left, bytes;
} Arena;
static const size_t kArenaBlock = (S_CAST(size_t, 8) << 20) / sizeof(Entry);

// Pool for the per-daughter run vectors.  They grow by doubling and are freed
// after their daughter is merged; through malloc, the freed blocks stayed
// resident as fragmented free space (up to ~125 MiB with a 250 MiB limit, on
// top of the tracked workspace).  Here small blocks come in power-of-two size
// classes carved from OS-mapped slabs and are reused within their class,
// larger blocks are mapped individually, and the pool's whole footprint counts
// toward the memory limit.  The pool is emptied whenever no run vector is
// live.
typedef struct RunPoolStruct {
  char** slabs;
  size_t nslab, slab_cap;
  void* free_list[64];
  char* cur;
  size_t left, bytes, live;
} RunPool;
static const size_t kRunSlab = S_CAST(size_t, 2) << 20;
// blocks this large are mapped individually and unmapped when freed (freed
// large blocks parked in free lists were ~2/3 of the pool after flushes)
static const size_t kRunDirectMin = S_CAST(size_t, 16) << 10;

// Simple open-addressing map from 64-bit key to int (key 0 reserved).
typedef struct NetMapStruct {
  uint64_t* k;
  int* v;
  size_t cnt, mask, cap;
} NetMap;

// All state of one network computation.  FisherNetworkLnP() allocates one per
// call (zero-initialized), so concurrent calls are independent.
struct NetStruct {
  double* LF;   // lfact, high part (double-double with LFL)
  double* LFL;  // lfact, low part
  uint64_t* FH; // additive hash of x! (prime-exponent based); equal products <=> equal hashes (w.h.p.)
  size_t lf_ct;
  double kTieEps;    // decision margin; set per call from the magnitude of the log-factorials
  double kTieEpsDD;  // double-double resolution margin (see Core())
  FisherCheckFn check_fn;
  void* check_arg;
  uint64_t poll_work;
  uint32_t err;   // kFisherNetOk, or the status of the first failure
  uint32_t midp;
  int NR;         // number of rows (after transposing so rows <= columns)

  // ---- memory limit ----
  // Tracked: past-value lists, run descriptors, node records, deferred-mode
  // weight arrays, hash tables, merge/sort scratch buffers, and the
  // log-factorial tables (see MemScratch()).
  size_t mem_limit, mem_peak;
  size_t mem_mothers, mem_runs, mem_nodes, mem_new, flush_at;
  uint64_t push_ctr;

  Arena arena_cur, arena_nxt, arena_prev;   // lists of cur, nxt, and (deferred mode) prevm
  Arena* out_arena;                         // where MergeRuns() puts new lists
  RunPool pool;
  SlotVec tab;           // MergeRuns() hash table
  U32Vec seg, seg2;      // MergeRuns() segment bounds
  U32Vec bcnt, bsort_start;
  U64Vec cnt_tab;        // CountNewDistinct() table
  EntryVec mbuf;         // MergeRuns() scratch
  NodeVec cur, nxt, prevm;
  RunVecVec druns, prevruns;
  NetMap nodemap;        // this stage's daughter nodes
  DoubleVec cum;         // this mother's prefix sums
  U32Vec bidx;           // last-stage bucket index
  LeafVec lvbuf;         // deferred last stage: leaf records
  DoubleVec wsuf;
  U32Vec qb;
  FlushOrdVec ord;
  IntVec mode_i;         // MaxCompletion() scratch
  DoubleVec mode_d;
  IntVec cb;             // ConcentrateBound() scratch
  int* ibuf;             // Core()'s copies of the table and its margins
};

static BoolErr NetFail(Net* net, uint32_t code) {
  if (!net->err) {
    net->err = code;
  }
  return 1;
}

static void* NetRegrow(Net* net, void* old, size_t used_bytes, size_t old_bytes, size_t new_bytes) {
  void* q = BigAllocRaw(new_bytes);
  if (!q) {
    NetFail(net, kFisherNetMemLimit);
    return nullptr;
  }
  if (used_bytes) {
    memcpy(q, old, used_bytes);
  }
  BigFreeRaw(old, old_bytes);
  return q;
}

static inline dd_real DLF(const Net* net, long long x) {
  return ddr_make(net->LF[x], net->LFL[x]);
}

// ---- cancellation hook ----
static BoolErr PollSlow(Net* net) {
  net->poll_work = 0;
  if (net->check_fn && net->check_fn(net->check_arg)) {
    return NetFail(net, kFisherNetCanceled);
  }
  return 0;
}
// roughly every 1-10 ms
static inline BoolErr Poll(Net* net, uint64_t work) {
  net->poll_work += work;
  if (net->poll_work >= (S_CAST(uint64_t, 1) << 18)) {
    return PollSlow(net);
  }
  return 0;
}

// ---- arena ----
static Entry* ArenaNewBlock(Net* net, Arena* a, size_t n) {
  if (a->nblk == a->blk_cap) {
    const size_t nc = a->blk_cap? (2 * a->blk_cap) : 8;
    Entry** nb = S_CAST(Entry**, malloc(nc * sizeof(Entry*)));
    size_t* ns = S_CAST(size_t*, malloc(nc * sizeof(size_t)));
    if ((!nb) || (!ns)) {
      free(nb);
      free(ns);
      NetFail(net, kFisherNetMemLimit);
      return nullptr;
    }
    if (a->nblk) {
      memcpy(nb, a->blocks, a->nblk * sizeof(Entry*));
      memcpy(ns, a->sizes, a->nblk * sizeof(size_t));
    }
    free(a->blocks);
    free(a->sizes);
    a->blocks = nb;
    a->sizes = ns;
    a->blk_cap = nc;
  }
  Entry* b = S_CAST(Entry*, BigAllocRaw(n * sizeof(Entry)));
  if (!b) {
    NetFail(net, kFisherNetMemLimit);
    return nullptr;
  }
  a->blocks[a->nblk] = b;
  a->sizes[a->nblk] = n * sizeof(Entry);
  ++a->nblk;
  a->bytes += n * sizeof(Entry);
  return b;
}
// Returns nullptr for n == 0, or on failure (net->err set).
static Entry* ArenaAlloc(Net* net, Arena* a, size_t n) {
  if (!n) {
    return nullptr;
  }
  if (n > kArenaBlock / 4) {
    return ArenaNewBlock(net, a, n);
  }
  if (n > a->left) {
    Entry* b = ArenaNewBlock(net, a, kArenaBlock);
    if (!b) {
      return nullptr;
    }
    a->cur = b;
    a->left = kArenaBlock;
  }
  Entry* p = a->cur;
  a->cur += n;
  a->left -= n;
  return p;
}
// returns a list's space if it had a block of its own (smaller lists' space is
// reclaimed with the arena)
static void ArenaFree(Arena* a, Entry* p, size_t n) {
  if (n <= kArenaBlock / 4) {
    return;
  }
  for (size_t i = a->nblk; i--; ) {
    if (a->blocks[i] != p) {
      continue;
    }
    BigFreeRaw(p, a->sizes[i]);
    a->bytes -= a->sizes[i];
    memmove(&a->blocks[i], &a->blocks[i + 1], (a->nblk - i - 1) * sizeof(Entry*));
    memmove(&a->sizes[i], &a->sizes[i + 1], (a->nblk - i - 1) * sizeof(size_t));
    --a->nblk;
    return;
  }
}
static void ArenaClear(Arena* a) {
  for (size_t i = 0; i < a->nblk; ++i) {
    BigFreeRaw(a->blocks[i], a->sizes[i]);
  }
  a->nblk = 0;
  a->cur = nullptr;
  a->left = 0;
  a->bytes = 0;
}
static void ArenaDestroy(Arena* a) {
  ArenaClear(a);
  free(a->blocks);
  free(a->sizes);
  a->blocks = nullptr;
  a->sizes = nullptr;
  a->blk_cap = 0;
}
static void ArenaSwap(Arena* a, Arena* b) {
  const Arena t = *a;
  *a = *b;
  *b = t;
}

// ---- run pool and run vectors ----
static inline int RunSizeClass(size_t n) {
  int c = 0;
  while ((S_CAST(size_t, 1) << c) < n) {
    ++c;
  }
  return c;
}
static Run* RunPoolAlloc(Net* net, RunPool* pl, size_t n) {
  const int c = RunSizeClass(n);
  const size_t sz = (S_CAST(size_t, 1) << c) * sizeof(Run);
  if (sz >= kRunDirectMin) {
    void* p = MapRaw(sz);
    if (!p) {
      NetFail(net, kFisherNetMemLimit);
      return nullptr;
    }
    ++pl->live;
    pl->bytes += sz;
    return S_CAST(Run*, p);
  }
  if (pl->free_list[c]) {
    void* p = pl->free_list[c];
    pl->free_list[c] = *S_CAST(void**, p);
    ++pl->live;
    return S_CAST(Run*, p);
  }
  if (sz > pl->left) {
    if (pl->nslab == pl->slab_cap) {
      const size_t nc = pl->slab_cap? (2 * pl->slab_cap) : 1;
      char** ns = S_CAST(char**, malloc(nc * sizeof(char*)));
      if (!ns) {
        NetFail(net, kFisherNetMemLimit);
        return nullptr;
      }
      if (pl->nslab) {
        memcpy(ns, pl->slabs, pl->nslab * sizeof(char*));
      }
      free(pl->slabs);
      pl->slabs = ns;
      pl->slab_cap = nc;
    }
    char* sl = S_CAST(char*, BigAllocRaw(kRunSlab));
    if (!sl) {
      NetFail(net, kFisherNetMemLimit);
      return nullptr;
    }
    pl->slabs[pl->nslab++] = sl;
    pl->cur = sl;
    pl->left = kRunSlab;
    pl->bytes += kRunSlab;
  }
  void* p = pl->cur;
  pl->cur += sz;
  pl->left -= sz;
  ++pl->live;
  return S_CAST(Run*, p);
}
static void RunPoolFree(RunPool* pl, void* p, size_t n) {
  const int c = RunSizeClass(n);
  const size_t sz = (S_CAST(size_t, 1) << c) * sizeof(Run);
  --pl->live;
  if (sz >= kRunDirectMin) {
    UnmapRaw(p, sz);
    pl->bytes -= sz;
    return;
  }
  *S_CAST(void**, p) = pl->free_list[c];
  pl->free_list[c] = p;
}
static inline size_t RunPoolBytes(const RunPool* pl) {
  return pl->bytes + pl->slab_cap * sizeof(char*);
}
static void RunPoolReset(RunPool* pl) {
  for (size_t i = 0; i < pl->nslab; ++i) {
    BigFreeRaw(pl->slabs[i], kRunSlab);
  }
  free(pl->slabs);
  pl->slabs = nullptr;
  pl->nslab = 0;
  pl->slab_cap = 0;
  for (int c = 0; c < 64; ++c) {
    pl->free_list[c] = nullptr;
  }
  pl->cur = nullptr;
  pl->left = 0;
  pl->bytes = 0;
  pl->live = 0;
}
static BoolErr RunVecPush(Net* net, RunVec* v, const Run* run) {
  if (v->n == v->cap) {
    const size_t nc = v->cap? (2 * v->cap) : 1;
    Run* q = RunPoolAlloc(net, &net->pool, nc);
    if (!q) {
      return 1;
    }
    if (v->n) {
      memcpy(q, v->p, v->n * sizeof(Run));
    }
    if (v->p) {
      RunPoolFree(&net->pool, v->p, v->cap);
    }
    v->p = q;
    v->cap = nc;
  }
  v->p[v->n++] = *run;
  return 0;
}
static void RunVecRelease(Net* net, RunVec* v) {
  if (v->p) {
    RunPoolFree(&net->pool, v->p, v->cap);
  }
  v->p = nullptr;
  v->n = 0;
  v->cap = 0;
}
static void RunVecVecRelease(Net* net, RunVecVec* vv) {
  for (size_t i = 0; i < vv->n; ++i) {
    RunVecRelease(net, &vv->p[i]);
  }
  RunVecVecFree(vv);
}

// ---- nodes ----
static void NodeVecRelease(NodeVec* nv) {
  for (size_t i = 0; i < nv->n; ++i) {
    Node* nd = &nv->p[i];
    BigFreeRaw(nd->wv, nd->wv_n * sizeof(double));
    BigFreeRaw(nd->cumv, nd->cumv_n * sizeof(double));
  }
  NodeVecFree(nv);
}
static inline size_t NodeBytes(const Node* nd) {
  // (lists are in the arenas)
  return sizeof(Node) + (nd->wv_n + nd->cumv_n) * sizeof(double);
}

// ---- node map ----
static inline uint64_t Mix(uint64_t x) {
  x ^= x >> 31;
  x *= 0x7FB5D329728EA185ULL;
  x ^= x >> 27;
  x *= 0x81DADEF4BC2DD44DULL;
  x ^= x >> 33;
  return x;
}
static void MapFree(NetMap* m) {
  BigFreeRaw(m->k, m->cap * sizeof(uint64_t));
  BigFreeRaw(m->v, m->cap * sizeof(int));
  m->k = nullptr;
  m->v = nullptr;
  m->cnt = 0;
  m->mask = 0;
  m->cap = 0;
}
// replaces *m's arrays (without freeing them) with c empty slots
static BoolErr MapAllocSlots(Net* net, NetMap* m, size_t c) {
  uint64_t* k = S_CAST(uint64_t*, BigAllocRaw(c * sizeof(uint64_t)));
  int* v = S_CAST(int*, BigAllocRaw(c * sizeof(int)));
  if ((!k) || (!v)) {
    BigFreeRaw(k, c * sizeof(uint64_t));
    BigFreeRaw(v, c * sizeof(int));
    return NetFail(net, kFisherNetMemLimit);
  }
  memset(k, 0, c * sizeof(uint64_t));
  for (size_t i = 0; i < c; ++i) {
    v[i] = -1;
  }
  m->k = k;
  m->v = v;
  m->cap = c;
  m->mask = c - 1;
  m->cnt = 0;
  return 0;
}
static BoolErr MapInit(Net* net, NetMap* m, size_t capn) {
  size_t c = 16;
  while (c < 2 * capn) {
    c <<= 1;
  }
  MapFree(m);
  return MapAllocSlots(net, m, c);
}
static BoolErr MapGrow(Net* net, NetMap* m) {
  NetMap old = *m;
  if (MapAllocSlots(net, m, old.cap * 2)) {
    return 1;
  }
  for (size_t i = 0; i < old.cap; ++i) {
    if (old.k[i]) {
      size_t p = Mix(old.k[i]) & m->mask;
      while (m->k[p]) {
        p = (p + 1) & m->mask;
      }
      m->k[p] = old.k[i];
      m->v[p] = old.v[i];
      ++m->cnt;
    }
  }
  MapFree(&old);
  return 0;
}
// returns the slot's value (created as -1 if absent), or nullptr on failure
static inline int* MapAt(Net* net, NetMap* m, uint64_t key) {
  if ((2 * (m->cnt + 1) > m->cap) && MapGrow(net, m)) {
    return nullptr;
  }
  size_t p = Mix(key) & m->mask;
  while (m->k[p] && (m->k[p] != key)) {
    p = (p + 1) & m->mask;
  }
  if (!m->k[p]) {
    m->k[p] = key;
    ++m->cnt;
  }
  return &m->v[p];
}

// ---- memory accounting ----
static size_t MemScratch(const Net* net) {
  size_t b = 3 * net->lf_ct * 8;
  b += net->tab.cap * sizeof(Slot) + net->mbuf.cap * sizeof(Entry) + (net->seg.cap + net->seg2.cap + net->bcnt.cap) * 4;
  b += net->nodemap.cap * 8 + net->nodemap.cap * 4;
  b += net->cum.cap * 8;
  b += net->bidx.cap * 4;
  b += net->lvbuf.cap * sizeof(LeafRec) + net->wsuf.cap * 8 + net->qb.cap * 4;
  return b;
}
static inline size_t MemNow(const Net* net) {
  return net->mem_mothers + RunPoolBytes(&net->pool) + net->mem_nodes + net->mem_new + MemScratch(net);
}
static inline BoolErr MemCheck(Net* net, size_t extra) {
  const size_t m = MemNow(net) + extra;
  if (m > net->mem_peak) {
    net->mem_peak = m;
  }
  if (m > net->mem_limit) {
    return NetFail(net, kFisherNetMemLimit);
  }
  return 0;
}

static BoolErr InitTables(Net* net, int n) {
  const size_t ct = S_CAST(size_t, n) + 2;
  net->lf_ct = ct;
  net->LF = S_CAST(double*, BigAllocRaw(ct * sizeof(double)));
  net->LFL = S_CAST(double*, BigAllocRaw(ct * sizeof(double)));
  net->FH = S_CAST(uint64_t*, BigAllocRaw(ct * sizeof(uint64_t)));
  if ((!net->LF) || (!net->LFL) || (!net->FH)) {
    return NetFail(net, kFisherNetMemLimit);
  }
  double* LF = net->LF;
  double* LFL = net->LFL;
  uint64_t* FH = net->FH;
  memset(FH, 0, ct * sizeof(uint64_t));
  for (int i = 0; i <= n + 1; ++i) {
    const dd_real t = ddr_lfact(S_CAST(double, i));
    LF[i] = t.x[0];
    LFL[i] = t.x[1];
    if ((!(i & 4095)) && Poll(net, 4096)) {
      return 1;
    }
  }
  uint64_t* R = S_CAST(uint64_t*, malloc(ct * sizeof(uint64_t)));
  if (!R) {
    return NetFail(net, kFisherNetMemLimit);
  }
  uint64_t s = 0x9E3779B97F4A7C15ULL;
  for (int i = 0; i <= n + 1; ++i) {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    R[i] = s * 0xD6E8FEB86659FD93ULL + S_CAST(uint64_t, i);
  }
  // hash(x!) = hash((x-1)!) + sum over prime factorization of x of e * R[p]
  for (int x = 2; x <= n + 1; ++x) {
    uint64_t h = 0;
    int y = x;
    for (int p = 2; p * p <= y; ++p) {
      while (y % p == 0) {
        h += R[p];
        y /= p;
      }
    }
    if (y > 1) {
      h += R[y];
    }
    FH[x] = FH[x - 1] + h;
  }
  free(R);
  return 0;
}

// first index in [lo, hi) whose v exceeds lim (std::upper_bound)
static inline size_t UpperBoundV(const Entry* a, size_t lo, size_t hi, double lim) {
  size_t len = hi - lo;
  while (len) {
    const size_t half = len >> 1;
    if (a[lo + half].v <= lim) {
      lo += half + 1;
      len -= half + 1;
    } else {
      len = half;
    }
  }
  return lo;
}

// Compares prod a[0..na) with prod b[0..nb) exactly (factors < 2^32); returns
// -1, 0, or 1.  work must have room for na + nb + 2 values.
static size_t ProdLimbs(const uint32_t* f, size_t n, uint32_t* pr) {
  size_t len = 1;
  pr[0] = 1;
  for (size_t k = 0; k < n; ++k) {
    if (!f[k]) {
      pr[0] = 0;
      return 1;
    }
    uint64_t carry = 0;
    for (size_t d = 0; d < len; ++d) {
      const uint64_t t = S_CAST(uint64_t, pr[d]) * f[k] + carry;
      pr[d] = S_CAST(uint32_t, t);
      carry = t >> 32;
    }
    if (carry) {
      pr[len++] = S_CAST(uint32_t, carry);
    }
  }
  return len;
}
static int CmpProducts(const uint32_t* a, size_t na, const uint32_t* b, size_t nb, uint32_t* work) {
  uint32_t* pa = work;
  uint32_t* pb = &work[na + 1];
  const size_t la = ProdLimbs(a, na, pa);
  const size_t lb = ProdLimbs(b, nb, pb);
  if (la != lb) {
    return (la < lb)? -1 : 1;
  }
  for (size_t k = la; k-- > 0; ) {
    if (pa[k] != pb[k]) {
      return (pa[k] < pb[k])? -1 : 1;
    }
  }
  return 0;
}

static BoolErr U32VecCopy(Net* net, U32Vec* dst, const U32Vec* src) {
  if (U32VecReserve(net, dst, src->n)) {
    return 1;
  }
  if (src->n) {
    memcpy(dst->p, src->p, src->n * sizeof(uint32_t));
  }
  dst->n = src->n;
  return 0;
}
static BoolErr U32VecAppend(Net* net, U32Vec* dst, const U32Vec* src) {
  if (U32VecReserve(net, dst, dst->n + src->n)) {
    return 1;
  }
  if (src->n) {
    memcpy(&dst->p[dst->n], src->p, src->n * sizeof(uint32_t));
  }
  dst->n += src->n;
  return 0;
}

// Exact optimality check for the R x m transportation table x (row-major) in
// ModeGeneral(): x is a mode iff no cycle in the residual graph (edge row i ->
// col j: cost log(x_ij + 1); col j -> row i if x_ij > 0: cost -log(x_ij)) has
// negative cost, i.e. prod (x + 1) over its increased cells < prod x over its
// decreased cells.  pot (the R + m potentials left by the double-precision
// Bellman-Ford, or null) and cup[q] = log(x[q] + 1) only narrow the search:
// every edge's reduced cost cost(u, v) + pot[u] - pot[v] is above ~-1e-12, so
// all edges of an improving cycle (at most 2R of them; the reduced costs sum to
// the cycle's cost) have reduced cost below 1e-9.  If those edges form no
// directed cycle, x is optimal; otherwise Bellman-Ford runs on them with exact
// path costs (lists of factors, compared with CmpProducts()).  Returns 1 iff
// x is optimal (or on failure, with net->err set); otherwise an improving
// cycle has been pushed (one unit) and 0 is returned.
static int ModeCertify(Net* net, int* x, int R, int m, const double* pot, const double* cup) {
  const int V = R + m;
  int result = 1;
  CertEdgeVec te = {nullptr, 0, 0};
  IntVec iw = {nullptr, 0, 0};
  U32Vec* lists = nullptr;
  U32Vec lhs = {nullptr, 0, 0};
  U32Vec rhs = {nullptr, 0, 0};
  U32Vec cn = {nullptr, 0, 0};
  U32Vec cd = {nullptr, 0, 0};
  U32Vec work = {nullptr, 0, 0};
  IntVec ups = {nullptr, 0, 0};
  IntVec downs = {nullptr, 0, 0};
  {
    for (int i = 0; i < R; ++i) {
      for (int j = 0; j < m; ++j) {
        const int q = i * m + j;
        if ((!pot) || (cup[q] + pot[i] - pot[R + j] < 1e-9)) {
          CertEdge e;
          e.u = i;
          e.v = R + j;
          e.q = q;
          e.up = 1;
          if (CertEdgeVecPush(net, &te, e)) {
            goto ModeCertify_ret;
          }
        }
        if (x[q] && ((!pot) || (-log(S_CAST(double, x[q])) + pot[R + j] - pot[i] < 1e-9))) {
          CertEdge e;
          e.u = R + j;
          e.v = i;
          e.q = q;
          e.up = 0;
          if (CertEdgeVecPush(net, &te, e)) {
            goto ModeCertify_ret;
          }
        }
      }
    }
    const int nte = S_CAST(int, te.n);
    if (pot) {
      // any directed cycle among them?  (iterative DFS)
      if (IntVecAssign(net, &iw, S_CAST(size_t, 5 * V + 1 + nte))) {
        goto ModeCertify_ret;
      }
      int* head = iw.p;
      int* adj = &head[V + 1];
      int* pos = &adj[nte];
      int* color = &pos[V];
      int* sv = &color[V];
      int* sp = &sv[V];
      for (int v = 0; v <= V; ++v) {
        head[v] = 0;
      }
      for (int k = 0; k < nte; ++k) {
        ++head[te.p[k].u + 1];
      }
      for (int v = 0; v < V; ++v) {
        head[v + 1] += head[v];
      }
      for (int v = 0; v < V; ++v) {
        pos[v] = head[v];
        color[v] = 0;
      }
      for (int k = 0; k < nte; ++k) {
        adj[pos[te.p[k].u]++] = te.p[k].v;
      }
      int cyc = 0;
      for (int s0 = 0; (s0 < V) && (!cyc); ++s0) {
        if (color[s0]) {
          continue;
        }
        color[s0] = 1;
        sv[0] = s0;
        sp[0] = head[s0];
        int top = 1;
        while (top && (!cyc)) {
          const int tv = sv[top - 1];
          if (sp[top - 1] == head[tv + 1]) {
            color[tv] = 2;
            --top;
            continue;
          }
          const int w = adj[sp[top - 1]++];
          if (color[w] == 1) {
            cyc = 1;
          } else if (!color[w]) {
            color[w] = 1;
            sv[top] = w;
            sp[top] = head[w];
            ++top;
          }
        }
      }
      if (!cyc) {
        goto ModeCertify_ret;
      }
    }
    // Bellman-Ford from a virtual source; a path's cost is log(prod num / prod den)
    lists = S_CAST(U32Vec*, calloc(2 * S_CAST(size_t, V), sizeof(U32Vec)));
    if ((!lists) || IntVecAssign(net, &iw, S_CAST(size_t, V))) {
      NetFail(net, kFisherNetMemLimit);
      goto ModeCertify_ret;
    }
    U32Vec* num = lists;
    U32Vec* den = &lists[V];
    int* pe = iw.p;
    for (int v = 0; v < V; ++v) {
      pe[v] = -1;
    }
    int lastv = -1;
    for (int it = 0; it <= V; ++it) {
      lastv = -1;
      for (int k = 0; k < nte; ++k) {
        const CertEdge* e = &te.p[k];
        const uint32_t f = e->up? (S_CAST(uint32_t, x[e->q]) + 1) : S_CAST(uint32_t, x[e->q]);
        // relax iff cost(u) + cost(e) < cost(v), i.e. num_u den_v (f if up) < num_v den_u (f if down)
        lhs.n = 0;
        rhs.n = 0;
        if (U32VecAppend(net, &lhs, &num[e->u]) || U32VecAppend(net, &lhs, &den[e->v]) ||
            U32VecAppend(net, &rhs, &num[e->v]) || U32VecAppend(net, &rhs, &den[e->u]) ||
            U32VecPush(net, e->up? (&lhs) : (&rhs), f) ||
            U32VecReserve(net, &work, lhs.n + rhs.n + 2)) {
          goto ModeCertify_ret;
        }
        if (CmpProducts(lhs.p, lhs.n, rhs.p, rhs.n, work.p) >= 0) {
          continue;
        }
        if (U32VecCopy(net, &num[e->v], &num[e->u]) || U32VecCopy(net, &den[e->v], &den[e->u]) ||
            U32VecPush(net, e->up? (&num[e->v]) : (&den[e->v]), f)) {
          goto ModeCertify_ret;
        }
        pe[e->v] = k;
        lastv = e->v;
      }
      if (lastv < 0) {
        goto ModeCertify_ret;
      }
    }
    // still relaxing after V + 1 passes: the predecessor graph has an improving cycle
    int v = lastv;
    for (int z = 0; z < V; ++z) {
      if (pe[v] < 0) {
        goto ModeCertify_ret;
      }
      v = te.p[pe[v]].u;
    }
    const int start = v;
    do {
      if (pe[v] < 0) {
        goto ModeCertify_ret;
      }
      const CertEdge* e = &te.p[pe[v]];
      if (e->up) {
        if (U32VecPush(net, &cn, S_CAST(uint32_t, x[e->q]) + 1) || IntVecPush(net, &ups, e->q)) {
          goto ModeCertify_ret;
        }
      } else {
        if (U32VecPush(net, &cd, S_CAST(uint32_t, x[e->q])) || IntVecPush(net, &downs, e->q)) {
          goto ModeCertify_ret;
        }
      }
      v = e->u;
    } while (v != start);
    if (U32VecReserve(net, &work, cn.n + cd.n + 2)) {
      goto ModeCertify_ret;
    }
    if (CmpProducts(cn.p, cn.n, cd.p, cd.n, work.p) >= 0) {
      goto ModeCertify_ret;  // (not expected)
    }
    for (size_t k = 0; k < ups.n; ++k) {
      ++x[ups.p[k]];
    }
    for (size_t k = 0; k < downs.n; ++k) {
      --x[downs.p[k]];
    }
    result = 0;
  }
 ModeCertify_ret:
  if (lists) {
    for (int v = 0; v < 2 * V; ++v) {
      U32VecFree(&lists[v]);
    }
    free(lists);
  }
  CertEdgeVecFree(&te);
  IntVecFree(&iw);
  U32VecFree(&lhs);
  U32VecFree(&rhs);
  U32VecFree(&cn);
  U32VecFree(&cd);
  U32VecFree(&work);
  IntVecFree(&ups);
  IntVecFree(&downs);
  return result;
}

// ---- bounds for the remaining NR x m problem (rows r, columns c[0..m-1]) ----
// Cells' factorial-product hash and double-double value (-sum lfact).
static void CellsHashVal(const Net* net, const int* cells, int n, uint64_t* h_out, dd_real* v_out) {
  uint64_t h = 0;
  dd_real v = ddr_make(0., 0.);
  for (int q = 0; q < n; ++q) {
    h += net->FH[cells[q]];
    v = ddr_sloppier_sub(v, DLF(net, cells[q]));
  }
  *h_out = h;
  *v_out = v;
}
// Mode of a two-way split: maximize -sum_i [lf(x_i) + lf(r_i - x_i)] subject
// to sum x = c0, 0 <= x_i <= r_i (an n x 2 table with row totals r and
// first-column total c0, or a 2 x n table with column totals r and first-row
// total c0).  Separable concave: proportional start, greedy fill, then pairwise
// exchanges while one is (exactly) an improvement.
static double ModeTwoWay(Net* net, const int* r, int n, long long c0, uint64_t* h_out, dd_real* v_out) {
  const double* LF = net->LF;
  int xbuf[64];
  int* x = xbuf;
  if (n > 64) {
    if (IntVecAssign(net, &net->mode_i, S_CAST(size_t, n))) {
      return 0.;
    }
    x = net->mode_i.p;
  }
  long long N = 0;
  for (int i = 0; i < n; ++i) {
    N += r[i];
  }
  long long s = 0;
  for (int i = 0; i < n; ++i) {
    x[i] = N? S_CAST(int, ((S_CAST(long long, r[i]) * c0) / N)) : 0;
    s += x[i];
  }
  for (; s < c0; ++s) {
    // add a unit where the marginal gain (r_i - x_i)/(x_i + 1) is largest
    int bi = -1;
    long long bn = 0;
    long long bd = 1;
    for (int i = 0; i < n; ++i) {
      const long long nn = r[i] - x[i];
      const long long dd = x[i] + 1;
      if ((nn > 0) && (nn * bd > bn * dd)) {
        bi = i;
        bn = nn;
        bd = dd;
      }
    }
    ++x[bi];
  }
  while (1) {
    // move a unit from j to i while it helps: (r_i - x_i) x_j > (x_i + 1)(r_j - x_j + 1)
    int bi = -1;
    int bj = -1;
    double br = 1.;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if ((i == j) || (x[j] == 0) || (x[i] == r[i])) {
          continue;
        }
        // both sides < 2^62: exact
        const long long nn = S_CAST(long long, (r[i] - x[i])) * x[j];
        const long long dd = S_CAST(long long, (x[i] + 1)) * (r[j] - x[j] + 1);
        if (nn <= dd) {
          continue;
        }
        const double ratio = S_CAST(double, nn) / S_CAST(double, dd);  // only ranks the improving exchanges
        if ((bi < 0) || (ratio > br)) {
          bi = i;
          bj = j;
          br = ratio;
        }
      }
    }
    if (bi < 0) {
      break;
    }
    ++x[bi];
    --x[bj];
  }
  double v = 0.;
  for (int i = 0; i < n; ++i) {
    v -= LF[x[i]] + LF[r[i] - x[i]];
  }
  if (h_out) {
    uint64_t h = 0;
    dd_real vd = ddr_make(0., 0.);
    for (int i = 0; i < n; ++i) {
      h += net->FH[x[i]] + net->FH[r[i] - x[i]];
      vd = ddr_sloppier_sub(ddr_sloppier_sub(vd, DLF(net, x[i])), DLF(net, r[i] - x[i]));
    }
    *h_out = h;
    *v_out = vd;
  }
  return v;
}
// Mode of a separable concave transportation problem (rows r[0..NR-1],
// columns c[0..m-1]): maximize -sum lf(x_ij).  Optimal iff the residual graph
// has no improving cycle; find one by Bellman-Ford in double precision (edge
// row i -> col j: cost log(x_ij + 1); col j -> row i if x_ij > 0: cost
// -log(x_ij)), push one unit around it, repeat; then certify optimality
// exactly (ModeCertify(), which pushes any improving cycle the tolerances
// hid).  Start from proportional rounding + northwest corner.
static double ModeGeneral(Net* net, const int* r, const int* c, int m, uint64_t* h_out, dd_real* v_out) {
  const int R = net->NR;
  const int V = R + m;
  if (IntVecAssign(net, &net->mode_i, S_CAST(size_t, R * m + m + V)) ||
      DoubleVecAssign(net, &net->mode_d, S_CAST(size_t, V + 2 * R * m))) {
    return 0.;
  }
  int* x = net->mode_i.p;
  int* cc = &x[R * m];
  int* pred = &cc[m];
  double* dist = net->mode_d.p;
  double* cup = &dist[V];
  double* cdn = &cup[R * m];
  for (int q = 0; q < R * m; ++q) {
    x[q] = 0;
  }
  for (int j = 0; j < m; ++j) {
    cc[j] = c[j];
  }
  int rr[kFisherNetMaxRows];
  long long N = 0;
  for (int i = 0; i < R; ++i) {
    rr[i] = r[i];
    N += r[i];
  }
  if (N == 0) {
    if (h_out) {
      *h_out = 0;
      *v_out = ddr_make(0., 0.);
    }
    return 0.;
  }
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < R; ++i) {
      int v = S_CAST(int, ((S_CAST(long long, r[i]) * c[j]) / N));
      v = MINV(v, MINV(rr[i], cc[j]));
      x[i * m + j] = v;
      rr[i] -= v;
      cc[j] -= v;
    }
  }
  for (int i = 0, j = 0; (i < R) && (j < m); ) {
    const int v = MINV(rr[i], cc[j]);
    x[i * m + j] += v;
    rr[i] -= v;
    cc[j] -= v;
    if (!rr[i]) {
      ++i;
    } else {
      ++j;
    }
  }
  for (int pot_ok = 0; ; pot_ok = 0) {
    while (1) {
      for (int q = 0; q < R * m; ++q) {
        cup[q] = log(x[q] + 1.0);
        cdn[q] = x[q]? (-log(S_CAST(double, x[q]))) : 1e300;
      }
      for (int v = 0; v < V; ++v) {
        dist[v] = 0.;
        pred[v] = -1;
      }
      int last = -1;
      for (int it = 0; it < V; ++it) {
        last = -1;
        for (int i = 0; i < R; ++i) {
          for (int j = 0; j < m; ++j) {
            const int q = i * m + j;
            if (dist[i] + cup[q] < dist[R + j] - 1e-12) {
              dist[R + j] = dist[i] + cup[q];
              pred[R + j] = i;
              last = R + j;
            }
            if (x[q] && (dist[R + j] + cdn[q] < dist[i] - 1e-12)) {
              dist[i] = dist[R + j] + cdn[q];
              pred[i] = R + j;
              last = i;
            }
          }
        }
        if (last < 0) {
          break;
        }
      }
      if (last < 0) {
        // no improving cycle beyond the tolerance; dist[] = potentials
        pot_ok = 1;
        break;
      }
      int v = last;
      for (int z = 0; z < V; ++z) {
        v = pred[v];  // now on the cycle
      }
      const int start = v;
      // walk the cycle backwards via pred; edge pred[v] -> v
      double gain = 0.;
      int u;
      do {
        u = pred[v];
        const int up = (u < R);
        const int i = up? u : v;
        const int j = (up? v : u) - R;
        gain -= up? cup[i * m + j] : cdn[i * m + j];
        v = u;
      } while (v != start);
      if (gain <= 1e-12) {
        break;  // too small to trust in double precision; decided exactly below
      }
      v = start;
      do {
        u = pred[v];
        if (u < R) {
          ++x[u * m + (v - R)];
        } else {
          --x[v * m + (u - R)];
        }
        v = u;
      } while (v != start);
    }
    if (ModeCertify(net, x, R, m, pot_ok? dist : nullptr, cup)) {
      break;
    }
  }
  double val = 0.;
  for (int q = 0; q < R * m; ++q) {
    val -= net->LF[x[q]];
  }
  if (h_out) {
    CellsHashVal(net, x, R * m, h_out, v_out);
  }
  return val;
}
// Value of the most probable completion (the mode).  If h_out is non-null, the
// mode table's factorial-product hash and double-double value are written to
// *h_out and *v_out.  On failure, net->err is set.
static double MaxCompletion(Net* net, const int* r, const int* c, int m, uint64_t* h_out, dd_real* v_out) {
  const int NR = net->NR;
  const double* LF = net->LF;
  if (m <= 2) {
    if (m == 1) {
      double v = 0.;
      for (int i = 0; i < NR; ++i) {
        v -= LF[r[i]];
      }
      if (h_out) {
        CellsHashVal(net, r, NR, h_out, v_out);
      }
      return v;
    }
    return ModeTwoWay(net, r, NR, c[0], h_out, v_out);
  }
  if (NR == 2) {
    return ModeTwoWay(net, c, m, r[0], h_out, v_out);
  }
  if (NR != 3) {
    return ModeGeneral(net, r, c, m, h_out, v_out);
  }
  // mode of a separable concave transportation problem; 4- and 6-cycles
  // cover all simple cycles with 3 rows, so local optimality is global.  A
  // cycle improves iff the product of its decreased cells' counts exceeds the
  // product of its increased cells' counts plus one; that is decided exactly,
  // so the search stops only at a mode.
  int xbuf[3 * 64];
  int ccbuf[64];
  int* x = xbuf;
  int* cc = ccbuf;
  if (m > 64) {
    if (IntVecAssign(net, &net->mode_i, 4 * S_CAST(size_t, m))) {
      return 0.;
    }
    x = net->mode_i.p;
    cc = &x[3 * m];
  }
  for (int q = 0; q < 3 * m; ++q) {
    x[q] = 0;
  }
  for (int j = 0; j < m; ++j) {
    cc[j] = c[j];
  }
  const long long N = S_CAST(long long, r[0]) + r[1] + r[2];
  if (N == 0) {
    if (h_out) {
      *h_out = 0;
      *v_out = ddr_make(0., 0.);
    }
    return 0.;
  }
  int rr[3];
  rr[0] = r[0];
  rr[1] = r[1];
  rr[2] = r[2];
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < 3; ++i) {
      int v = S_CAST(int, ((S_CAST(long long, r[i]) * c[j]) / N));
      v = MINV(v, MINV(rr[i], cc[j]));
      x[i * m + j] = v;
      rr[i] -= v;
      cc[j] -= v;
    }
  }
  for (int i = 0, j = 0; (i < 3) && (j < m); ) {
    // northwest corner on residual
    const int v = MINV(rr[i], cc[j]);
    x[i * m + j] += v;
    rr[i] -= v;
    cc[j] -= v;
    if (!rr[i]) {
      ++i;
    } else {
      ++j;
    }
  }
#define XL(ii, jj) x[(ii) * m + (jj)]
  static const int kPerm3[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  while (1) {
    int bt = 0;
    int a0 = 0;
    int a1 = 0;
    int a2 = 0;
    int b0 = 0;
    int b1 = 0;
    int b2 = 0;
    // gain > 0 iff num > den; improving cycles are ranked by num/den in double precision
    double br = 1.;
    for (int i1 = 0; i1 < 3; ++i1) {
      for (int i2 = 0; i2 < 3; ++i2) {
        if (i1 == i2) {
          continue;
        }
        for (int j1 = 0; j1 < m; ++j1) {
          for (int j2 = 0; j2 < m; ++j2) {
            if ((j1 == j2) || (!XL(i1, j2)) || (!XL(i2, j1))) {
              continue;
            }
            // both sides < 2^62: exact
            const long long nu = S_CAST(long long, XL(i1, j2)) * XL(i2, j1);
            const long long de = (XL(i1, j1) + 1LL) * (XL(i2, j2) + 1LL);
            if (nu <= de) {
              continue;
            }
            const double ratio = S_CAST(double, nu) / S_CAST(double, de);
            if ((!bt) || (ratio > br)) {
              br = ratio;
              bt = 4;
              a0 = i1;
              a1 = i2;
              b0 = j1;
              b1 = j2;
            }
          }
        }
      }
    }
    // 6-cycles: +(p0,j0) -(p0,j1) +(p1,j1) -(p1,j2) +(p2,j2) -(p2,j0), over row orders
    for (int pq = 0; pq < 6; ++pq) {
      const int p0 = kPerm3[pq][0];
      const int p1 = kPerm3[pq][1];
      const int p2 = kPerm3[pq][2];
      for (int j0 = 0; j0 < m; ++j0) {
        for (int j1 = 0; j1 < m; ++j1) {
          if (j1 == j0) {
            continue;
          }
          for (int j2 = 0; j2 < m; ++j2) {
            if ((j2 == j0) || (j2 == j1) || (!XL(p0, j1)) || (!XL(p1, j2)) || (!XL(p2, j0))) {
              continue;
            }
            // products of three counts: relative rounding error < 2^-51
            const double nu = S_CAST(double, XL(p0, j1)) * XL(p1, j2) * XL(p2, j0);
            const double de = (XL(p0, j0) + 1.0) * (XL(p1, j1) + 1.0) * (XL(p2, j2) + 1.0);
            if (nu < de * (1. - 1e-12)) {
              continue;
            }
            if (nu <= de * (1. + 1e-12)) {
              // too close for doubles: decide exactly
              uint32_t na[3];
              uint32_t da[3];
              uint32_t work[8];
              na[0] = S_CAST(uint32_t, XL(p0, j1));
              na[1] = S_CAST(uint32_t, XL(p1, j2));
              na[2] = S_CAST(uint32_t, XL(p2, j0));
              da[0] = S_CAST(uint32_t, XL(p0, j0)) + 1;
              da[1] = S_CAST(uint32_t, XL(p1, j1)) + 1;
              da[2] = S_CAST(uint32_t, XL(p2, j2)) + 1;
              if (CmpProducts(na, 3, da, 3, work) <= 0) {
                continue;
              }
            }
            const double ratio = nu / de;
            if ((!bt) || (ratio > br)) {
              br = ratio;
              bt = 6;
              a0 = p0;
              a1 = p1;
              a2 = p2;
              b0 = j0;
              b1 = j1;
              b2 = j2;
            }
          }
        }
      }
    }
    if (!bt) {
      break;
    }
    if (bt == 4) {
      ++XL(a0, b0);
      ++XL(a1, b1);
      --XL(a0, b1);
      --XL(a1, b0);
    } else {
      ++XL(a0, b0);
      --XL(a0, b1);
      ++XL(a1, b1);
      --XL(a1, b2);
      ++XL(a2, b2);
      --XL(a2, b0);
    }
  }
#undef XL
  double v = 0.;
  for (int q = 0; q < 3 * m; ++q) {
    v -= LF[x[q]];
  }
  if (h_out) {
    CellsHashVal(net, x, 3 * m, h_out, v_out);
  }
  return v;
}
// sum over tots of min over {x: sum = t, 0 <= x <= caps} of -sum lfact(x)
static double ConcentrateBound(Net* net, const int* caps_in, int ncap, const int* tots, int ntot) {
  const double* LF = net->LF;
  int cbuf[64];
  int* caps = cbuf;
  if (ncap > 64) {
    if (IntVecAssign(net, &net->cb, S_CAST(size_t, ncap))) {
      return 0.;
    }
    caps = net->cb.p;
    memcpy(caps, caps_in, ncap * sizeof(int));
    STD_SORT(ncap, CmpInt, caps);
    for (int i = 0, j = ncap - 1; i < j; ++i, --j) {
      // ascending -> descending
      const int t = caps[i];
      caps[i] = caps[j];
      caps[j] = t;
    }
  } else {
    for (int i = 0; i < ncap; ++i) {
      caps[i] = caps_in[ncap - 1 - i];  // ascending input -> descending
    }
    for (int i = 1; i < ncap; ++i) {
      const int v = caps[i];
      int q = i;
      while ((q > 0) && (caps[q - 1] < v)) {
        caps[q] = caps[q - 1];
        --q;
      }
      caps[q] = v;
    }
  }
  double v = 0.;
  for (int j = 0; j < ntot; ++j) {
    int t = tots[j];
    for (int i = 0; (i < ncap) && (t > 0); ++i) {
      const int x = MINV(t, caps[i]);
      v -= LF[x];
      t -= x;
    }
  }
  return v;
}
// Most even split of t over caps (water-filling): max of -sum lfact(x), sum x =
// t, x <= caps (caps sorted ascending).
static double EvenSplitMaxSorted(const double* LF, const int* caps, int n, long long t) {
  double v = 0.;
  for (int i = 0; i < n; ++i) {
    const long long left = n - i;
    const long long share = t / left;
    if (caps[i] <= share) {
      v -= LF[caps[i]];
      t -= caps[i];
      continue;
    }
    const long long q = t / left;
    const long long rmd = t % left;
    v -= S_CAST(double, rmd) * LF[q + 1] + S_CAST(double, (left - rmd)) * LF[q];
    return v;
  }
  return v;
}
// Cheap upper bound on the max completion of rows r[NR] (sorted) x columns
// c[0..m-1] (sorted ascending).
static double MaxCompletionUB(const Net* net, const int* r, const int* c, int m) {
  const double* LF = net->LF;
  double cw = 0.;
  double rw = 0.;
  for (int l = 0; l < m; ++l) {
    cw += EvenSplitMaxSorted(LF, r, net->NR, c[l]);
  }
  for (int i = 0; i < net->NR; ++i) {
    rw += EvenSplitMaxSorted(LF, c, m, r[i]);
  }
  return MINV(cw, rw);
}
static BoolErr BoundsLazy(Net* net, Node* nd, const int* c, int m) {
  const int NR = net->NR;
  nd->ub = MaxCompletionUB(net, nd->r, c, m);
  const double b1 = ConcentrateBound(net, nd->r, NR, c, m);
  const double b2 = ConcentrateBound(net, c, m, nd->r, NR);
  if (net->err) {
    return 1;
  }
  nd->splb = MAXV(b1, b2);
  long long N = 0;
  for (int i = 0; i < NR; ++i) {
    N += nd->r[i];
  }
  dd_real t = DLF(net, N);
  for (int i = 0; i < NR; ++i) {
    t = ddr_sloppier_sub(t, DLF(net, nd->r[i]));
  }
  nd->logt = t.x[0];
  nd->logtl = t.x[1];
  for (int j = 0; j < m; ++j) {
    nd->logt -= net->LF[c[j]];
  }
  nd->have_bounds = 1;
  return 0;
}

// ---- Last stage: outward-from-mode enumeration over the 2 x NR daughter
// problem ----
// Daughter d (sum c, 0 <= d_i <= r_i) has table value tv(d) = -sum_i g_i(d_i),
// g_i(x) = lf(x) + lf(r_i - x).  Past value p counts iff p + tv(d) <= thr.
typedef struct LastCtxStruct {
  LeafVec* leafrec;   // deferred mode: record (threshold, weight) per leaf
  const uint32_t* bidx;   // bucket index
  double bvmin, binv;
  int bnb;
  const int* r;
  const Entry* srt;
  size_t n;
  const double* cum;
  double cum0;   // subtracted from cum[] (for a sub-range of the mother list)
  double pmax, pmin, cum_all, vobs, thr;
  double pmaxl, vobsl;   // low parts
  // exact ties: factorial-product hash of the observed table; weight of a tie
  // (0.5 for mid-p)
  uint64_t obs_h;
  double tie_w;
  long long S[kFisherNetMaxRows + 1];
  double sumlf[kFisherNetMaxRows + 1];
  double sumlfl[kFisherNetMaxRows + 1];
  Acc* acc;
  // previous jump targets per direction
  int xc[2];
  uint32_t have_xc[2];
  // previous tail per direction: m, tail start, tail sum (w-space), fixed,
  // weight at the fully counted leaf
  uint32_t pt_ok[2];
  int pt_reuse[2];
  long long pt_m[2];
  int pt_s[2];
  double pt_sum[2], pt_fixed[2], pt_fixedl[2], pt_wfull[2];
} LastCtx;

static inline double Gi(const Net* net, const LastCtx* L, int i, int x) {
  return net->LF[x] + net->LF[L->r[i] - x];
}
static inline dd_real GiD(const Net* net, const LastCtx* L, int i, int x) {
  return ddr_sloppier_add(DLF(net, x), DLF(net, L->r[i] - x));
}
static inline dd_real LogTotD(const Net* net, const LastCtx* L, int j, long long m) {
  return ddr_sloppier_sub(ddr_sloppier_sub(ddr_sloppier_sub(DLF(net, L->S[j]), DLF(net, m)), DLF(net, L->S[j] - m)), ddr_make(L->sumlf[j], L->sumlfl[j]));
}
// Resolve the entries srt[q0..) whose value is inside the tie window at a leaf
// with completion value tvD and completion hash leafh; adds their
// contributions (scale w relative to the list's cum scaling).
static void ResolveWindow(Net* net, LastCtx* L, size_t q0, double lim_win, dd_real tvD, uint64_t leafh, double w) {
  const dd_real pmaxD = ddr_make(L->pmax, L->pmaxl);
  const dd_real vobsD = ddr_make(L->vobs, L->vobsl);
  const double eps = net->kTieEpsDD;
  for (size_t q = q0; (q < L->n) && (L->srt[q].v <= lim_win); ++q) {
    const Entry* e = &L->srt[q];
    const double wt = e->f * ddr_exp_to_double(ddr_sloppier_sub(ddr_make(e->v, e->vl), pmaxD));
    const dd_real diff = ddr_sloppier_sub(ddr_sloppier_add(ddr_make(e->v, e->vl), tvD), vobsD);
    const double d = diff.x[0] + diff.x[1];
    // outside the margin, the sign of d is reliable whatever the hashes say
    // (equal hashes there can only be a collision)
    if (fabs(d) > eps) {
      if (d < 0.) {
        AccAdd(L->acc, w * wt);
      }
      continue;
    }
    // inside the margin: a tie iff the factorial products are equal
    if (e->h + leafh != L->obs_h) {
      NetFail(net, kFisherNetUnsupported);
      return;
    }
    AccAdd(L->acc, w * L->tie_w * wt);
  }
}
// max over d[0..j-1] (sum m) of -sum g; mode written to md
static double GroupMax(const Net* net, const LastCtx* L, int j, long long m, int* md) {
  long long cur = 0;
  for (int i = 0; i < j; ++i) {
    md[i] = S_CAST(int, (L->S[j]? ((S_CAST(long long, L->r[i]) * m) / L->S[j]) : 0));
    cur += md[i];
  }
  for (; cur < m; ++cur) {
    int best = -1;
    long long bn = 0;
    long long bd = 1;
    for (int i = 0; i < j; ++i) {
      const long long nn = L->r[i] - md[i];
      const long long dd = md[i] + 1;
      if (nn * bd > bn * dd) {
        best = i;
        bn = nn;
        bd = dd;
      }
    }
    ++md[best];
  }
  while (1) {
    int bi = -1;
    int bdc = -1;
    long long in_ = 0;
    long long id_ = 1;
    long long dn = 0;
    long long dd_ = 1;
    for (int i = 0; i < j; ++i) {
      long long nn = L->r[i] - md[i];
      long long dd = md[i] + 1;
      if (nn * id_ > in_ * dd) {
        bi = i;
        in_ = nn;
        id_ = dd;
      }
      nn = md[i];
      dd = L->r[i] - md[i] + 1;
      if (nn * dd_ > dn * dd) {
        bdc = i;
        dn = nn;
        dd_ = dd;
      }
    }
    if ((bi < 0) || (bdc < 0) || (bi == bdc) || (in_ * dn <= id_ * dd_)) {
      break;
    }
    ++md[bi];
    --md[bdc];
  }
  double v = 0.;
  for (int i = 0; i < j; ++i) {
    v -= Gi(net, L, i, md[i]);
  }
  return v;
}
// exact min over d[0..j-1] (sum m) of -sum g: at a vertex (all but one
// coordinate at a bound)
static double GroupMin(const Net* net, const LastCtx* L, int j, long long m) {
  double best = 1e300;
  for (int fr = 0; fr < j; ++fr) {
    for (int mask = 0; mask < (1 << j); ++mask) {
      if (mask & (1 << fr)) {
        continue;
      }
      long long s = 0;
      double v = 0.;
      for (int i = 0; i < j; ++i) {
        if (i != fr) {
          const int x = ((mask >> i) & 1)? L->r[i] : 0;
          s += x;
          v -= Gi(net, L, i, x);
        }
      }
      const long long xf = m - s;
      if ((xf < 0) || (xf > L->r[fr])) {
        continue;
      }
      v -= Gi(net, L, fr, S_CAST(int, xf));
      if (v < best) {
        best = v;
      }
    }
  }
  return best;
}
static void LastTail(const Net* net, LastCtx* L, int j, long long m, double fixed, double fixedl, int xs, int dir) {
  const int xlo = S_CAST(int, MAXV(0LL, m - L->S[j - 1]));
  const int xhi = S_CAST(int, MINV(S_CAST(long long, L->r[j - 1]), m));
  double local = 0.;
  double prev = 0.;
  const dd_real baseD = ddr_sloppier_sub(ddr_sloppier_add(ddr_make(L->pmax, L->pmaxl), ddr_make(fixed, fixedl)), ddr_make(L->vobs, L->vobsl));
  // term(x) = exp(base - g_{j-1}(x) + LogTot(j-1, m-x)); consecutive terms
  // differ by a rational factor:
  // x -> x+1: (r-x)(m-x) / ((x+1)(S-m+x+1)), with r = r_{j-1}, S = S_{j-1}
  const double rr = L->r[j - 1];
  const double SS = S_CAST(double, L->S[j - 1]);
  const double mm = S_CAST(double, m);
  double term = 0.;
  int since = 64;
  for (int x = xs; (x >= xlo) && (x <= xhi); x += dir) {
    if (since >= 64) {
      term = ddr_exp_to_double(ddr_sloppier_add(ddr_sloppier_sub(baseD, GiD(net, L, j - 1, x)), LogTotD(net, L, j - 1, m - x)));
      since = 0;
    } else {
      const double xp = x - dir;  // previous x
      term *= (dir > 0)? ((rr - xp) * (mm - xp) / ((xp + 1.) * (SS - mm + xp + 1.))) : (xp * (SS - mm + xp) / ((rr - xp + 1.) * (mm - xp + 1.)));
    }
    ++since;
    local += term;
    if ((term < prev) && (term < 1e-20 * local)) {
      break;
    }
    prev = term;
  }
  AccAdd(L->acc, local * L->cum_all);
}
// count of entries with v <= lim via the bucket index: bucket(v) = floor((v -
// vmin) * inv) is monotone in v, so entries in lower buckets are all <= lim and
// entries in higher buckets all > lim.
static inline size_t BucketCount(const Entry* a, const uint32_t* bidx, double vmin, double inv, int nb, double lim) {
  if (lim < vmin) {
    return 0;
  }
  const double bq = (lim - vmin) * inv;
  const int b = (bq >= S_CAST(double, nb))? nb : S_CAST(int, bq);
  const size_t s0 = bidx[b];
  const size_t s1 = bidx[b + 1];
  if (s1 - s0 <= 8) {
    size_t q = s0;
    while ((q < s1) && (a[q].v <= lim)) {
      ++q;
    }
    return q;
  }
  return UpperBoundV(a, s0, s1, lim);
}
// galloping upper_bound on a[lo..n): first index with a[i].v > lim
static inline size_t Advance(const Entry* a, size_t n, size_t lo, double lim) {
  if ((lo >= n) || (a[lo].v > lim)) {
    return lo;
  }
  size_t prev = lo;
  size_t step = 1;
  while ((lo + step < n) && (a[lo + step].v <= lim)) {
    prev = lo + step;
    step <<= 1;
  }
  const size_t hi = MINV(lo + step, n);
  return UpperBoundV(a, prev + 1, hi, lim);
}
static inline double L2TV(const Net* net, const LastCtx* L, double fixed, long long m, int x) {
  return fixed - Gi(net, L, 1, x) - Gi(net, L, 0, S_CAST(int, (m - x)));
}
static inline dd_real L2TVD(const Net* net, const LastCtx* L, double fixed, double fixedl, long long m, int x) {
  return ddr_sloppier_sub(ddr_sloppier_sub(ddr_make(fixed, fixedl), GiD(net, L, 1, x)), GiD(net, L, 0, S_CAST(int, (m - x))));
}
// exact leaf weight
static inline double L2WX(const Net* net, const LastCtx* L, dd_real baseD, long long m, int x) {
  return ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_sub(baseD, GiD(net, L, 1, x)), GiD(net, L, 0, S_CAST(int, (m - x)))));
}
// j == 2 level: leaves along x = d[1] outward from the mode x0. tv(x) is
// concave in x, so moving outward tv falls, the count threshold thr - tv rises
// (lo only advances), and the weight w = exp(pmax + tv - vobs) follows
// w(x+1)/w(x) = (r1-x)(m-x) / ((x+1)(r0-m+x+1)).  The tail after the first
// fully counted leaf has the same weights times cum_all.
static void Level2Fast(Net* net, LastCtx* L, long long m, double fixed, double fixedl, uint64_t fixedh, int x0) {
  const uint64_t* FH = net->FH;
  const double kTieEps = net->kTieEps;
  const int r0 = L->r[0];
  const int r1 = L->r[1];
  if (Poll(net, 64)) {
    return;
  }
  const dd_real baseD = ddr_sloppier_sub(ddr_sloppier_add(ddr_make(L->pmax, L->pmaxl), ddr_make(fixed, fixedl)), ddr_make(L->vobs, L->vobsl));
  const int xlo = S_CAST(int, MAXV(0LL, m - r0));
  const int xhi = S_CAST(int, MINV(S_CAST(long long, r1), m));
  size_t lo0 = 0;
  for (int dir = 1; dir >= -1; dir -= 2) {
    int x = (dir > 0)? x0 : (x0 - 1);
    if ((x < xlo) || (x > xhi)) {
      continue;
    }
    const int xend = (dir > 0)? xhi : xlo;
    double tv = L2TV(net, L, fixed, m, x);
    // tv falls monotonically outward; skip leaves where even the smallest past
    // value fails
    if (L->pmin + tv > L->thr) {
      if (L->pmin + L2TV(net, L, fixed, m, xend) > L->thr) {
        continue;  // nothing counts
      }
      int bad = x;  // P(bad) false, P(good) true
      int good = xend;
      const int di = (dir > 0)? 0 : 1;
      if (L->have_xc[di]) {
        int g = L->xc[di];
        if ((g - x) * dir < 0) {
          g = x;
        }
        if ((g - xend) * dir > 0) {
          g = xend;
        }
        if ((g != x) && (L->pmin + L2TV(net, L, fixed, m, g) <= L->thr)) {
          // crossing at or before g: gallop back toward x
          good = g;
          for (int step = 1; ; step <<= 1) {
            const int c = good - step * dir;
            if ((c - x) * dir <= 0) {
              break;
            }
            if (L->pmin + L2TV(net, L, fixed, m, c) <= L->thr) {
              good = c;
            } else {
              bad = c;
              break;
            }
          }
        } else {
          // crossing after g: gallop forward
          bad = g;
          for (int step = 1; ; step <<= 1) {
            const int c = bad + step * dir;
            if ((c - xend) * dir >= 0) {
              break;
            }
            if (L->pmin + L2TV(net, L, fixed, m, c) <= L->thr) {
              good = c;
              break;
            }
            bad = c;
          }
        }
      }
      while ((good - bad) * dir > 1) {
        const int mid = bad + (good - bad) / 2;
        if (L->pmin + L2TV(net, L, fixed, m, mid) <= L->thr) {
          good = mid;
        } else {
          bad = mid;
        }
      }
      L->xc[di] = good;
      L->have_xc[di] = 1;
      x = good;
      tv = L2TV(net, L, fixed, m, x);
    }
    double w = L2WX(net, L, baseD, m, x);
    int since = 0;
    size_t lo = lo0;
    int tail = 0;
    double local = 0.;
    // leaves until the first fully counted one
    for (;;) {
      {
        const double lim_def = L->vobs - kTieEps - tv;  // entries <= lim_def: strictly less probable
        const uint64_t leafh = fixedh + FH[x] + FH[r1 - x] + FH[m - x] + FH[r0 - m + x];
        if (L->leafrec) {
          LeafRec lr;
          lr.q = lim_def;
          lr.w = w;
          lr.h = leafh;
          lr.tv = fixed - Gi(net, L, 1, x) - Gi(net, L, 0, S_CAST(int, (m - x)));
          lr.tvl = L2TVD(net, L, fixed, fixedl, m, x).x[1];
          if (LeafVecPush(net, L->leafrec, lr)) {
            return;
          }
        } else {
          lo = L->bidx? BucketCount(L->srt, L->bidx, L->bvmin, L->binv, L->bnb, lim_def) : Advance(L->srt, L->n, lo, lim_def);
          if ((dir > 0) && (x == x0)) {
            lo0 = lo;
          }
          if (lo) {
            AccAdd(L->acc, w * (L->cum[lo] - L->cum0));
          }
          if ((lo < L->n) && (L->srt[lo].v <= L->thr - tv)) {
            ResolveWindow(net, L, lo, L->thr - tv, L2TVD(net, L, fixed, fixedl, m, x), leafh, w);
            if (net->err) {
              return;
            }
          }
        }
      }
      if (L->pmax + tv <= L->vobs - kTieEps) {
        tail = 1;
        break;
      }
      const int xn = x + dir;
      if ((xn < xlo) || (xn > xhi)) {
        break;
      }
      const double ratio = (dir > 0)? (S_CAST(double, (r1 - x)) * (m - x) / (S_CAST(double, (x + 1)) * (r0 - m + x + 1))) : (S_CAST(double, x) * (r0 - m + x) / (S_CAST(double, (r1 - x + 1)) * (m - x + 1)));
      x = xn;
      tv = L2TV(net, L, fixed, m, x);
      if (++since >= 64) {
        w = L2WX(net, L, baseD, m, x);
        since = 0;
      } else {
        w *= ratio;
      }
    }
    if (!tail) {
      L->pt_ok[(dir > 0)? 0 : 1] = 0;
      continue;
    }
    const int di = (dir > 0)? 0 : 1;
    const int st = x + dir;  // first tail x; x is the fully counted leaf
    const long long N2 = S_CAST(long long, r0) + r1;
    const double wfull = w;  // weight at the fully counted leaf (the direct tail loop overwrites w)
    int reused = 0;
    if (L->pt_ok[di] && (L->pt_reuse[di] < kTailReuseMax)) {
      const long long dm = m - L->pt_m[di];
      const int ds = st - L->pt_s[di];
      if ((di == 0) && (dm == 1) && (ds <= 0)) {
        // U_{m}(c) = [(N - m') U_{m'}(c) + (r1 - c + 1) h_{m'}(c - 1)] / (m' + 1), m' = m - 1, c = previous start
        const long long mp = L->pt_m[di];
        const int c = L->pt_s[di];
        double U = ddr_exp_to_double(ddr_sloppier_sub(ddr_make(fixed, fixedl), ddr_make(L->pt_fixed[di], L->pt_fixedl[di]))) * (S_CAST(double, (N2 - mp)) * L->pt_sum[di] + S_CAST(double, (r1 - c + 1)) * L->pt_wfull[di]) / S_CAST(double, (mp + 1));
        double ww = w;
        for (int xx = x; xx + 1 <= c - 1; ++xx) {
          ww *= S_CAST(double, (r1 - xx)) * (m - xx) / (S_CAST(double, (xx + 1)) * (r0 - m + xx + 1));
          U += ww;
        }
        local = U;
        reused = 1;
      } else if ((di == 1) && (dm == -1) && (ds >= 0)) {
        // D_m(c) = [(m + 1) D_{m+1}(c) + (r1 - c) h_m(c)] / (N - m), m + 1 = previous m, c = previous start
        const int c = L->pt_s[di];
        double ww = w;
        double add = 0.;
        for (int xx = x; xx - 1 >= c; --xx) {
          ww *= S_CAST(double, xx) * (r0 - m + xx) / (S_CAST(double, (r1 - xx + 1)) * (m - xx + 1));
          if (xx - 1 > c) {
            add += ww;
          }
        }
        local = (S_CAST(double, (m + 1)) * ddr_exp_to_double(ddr_sloppier_sub(ddr_make(fixed, fixedl), ddr_make(L->pt_fixed[di], L->pt_fixedl[di]))) * L->pt_sum[di] + S_CAST(double, (r1 - c)) * ww) / S_CAST(double, (N2 - m)) + add;
        reused = 1;
      }
    }
    if (reused) {
      ++L->pt_reuse[di];
    } else {
      L->pt_reuse[di] = 0;
      // tight tail loop: integer-free float recurrence on (num1, num2, den1, den2)
      // dir>0: ratio(x) = (r1-x)(m-x) / ((x+1)(r0-m+x+1)); dir<0: x (r0-m+x) / ((r1-x+1)(m-x+1))
      double a1;  // factors at the current x; each changes by -1 / +1 per step
      double a2;
      double b1;
      double b2;
      if (dir > 0) {
        a1 = r1 - x;
        a2 = S_CAST(double, (m - x));
        b1 = x + 1;
        b2 = S_CAST(double, (r0 - m + x + 1));
      } else {
        a1 = x;
        a2 = S_CAST(double, (r0 - m + x));
        b1 = r1 - x + 1;
        b2 = S_CAST(double, (m - x + 1));
      }
      const long long nsteps = (dir > 0)? (xhi - x) : (x - xlo);  // steps available after x
      local = 0.;  // the leaf at x was already counted in full; the tail starts at x + dir
      // ratio(t) = n_t / d_t with n_t = a1 a2, d_t = b1 b2; ratios decrease
      // outward (log-concave), so once n < d the remaining tail after weight w
      // is <= w n / (d - n).
      long long t = 0;
      for (; t + 4 <= nsteps; t += 4) {
        const double n1 = a1 * a2;
        const double d1 = b1 * b2;
        const double n2 = (a1 - 1.) * (a2 - 1.);
        const double d2 = (b1 + 1.) * (b2 + 1.);
        const double n3 = (a1 - 2.) * (a2 - 2.);
        const double d3 = (b1 + 2.) * (b2 + 2.);
        const double n4 = (a1 - 3.) * (a2 - 3.);
        const double d4 = (b1 + 3.) * (b2 + 3.);
        const double d34 = d3 * d4;
        const double inv = w / (d1 * d2 * d34);
        const double n12 = n1 * n2;
        const double n123 = n12 * n3;
        const double n1234 = n123 * n4;
        const double w1 = n1 * (d2 * d34) * inv;
        const double w2 = n12 * d34 * inv;
        const double w3 = n123 * d4 * inv;
        const double w4 = n1234 * inv;
        local += (w1 + w2) + (w3 + w4);
        w = w4;
        a1 -= 4.;
        a2 -= 4.;
        b1 += 4.;
        b2 += 4.;
        const double n5 = a1 * a2;
        const double d5 = b1 * b2;
        if ((n5 < d5) && (w * n5 < kTailEps * local * (d5 - n5))) {
          goto Level2Fast_tail_done;
        }
      }
      for (; t < nsteps; ++t) {
        const double nn = a1 * a2;
        const double dd = b1 * b2;
        if ((nn < dd) && (w * nn < kTailEps * local * (dd - nn))) {
          break;
        }
        const double ratio = nn / dd;  // keep the divide off the w dependency chain
        w *= ratio;
        local += w;
        a1 -= 1.;
        a2 -= 1.;
        b1 += 1.;
        b2 += 1.;
      }
    Level2Fast_tail_done:
      ;
    }
    L->pt_ok[di] = 1;
    L->pt_m[di] = m;
    L->pt_s[di] = st;
    L->pt_sum[di] = local;
    L->pt_fixed[di] = fixed;
    L->pt_fixedl[di] = fixedl;
    L->pt_wfull[di] = wfull;
    AccAdd(L->acc, local * L->cum_all);
  }
}
// returns 1 iff the rank is fully counted (and has been added); on failure,
// net->err is set
static int LastLevel(Net* net, LastCtx* L, int j, long long m, double fixed, double fixedl, uint64_t fixedh) {
  const double kTieEps = net->kTieEps;
  if (j == 1) {
    const dd_real tvD = ddr_sloppier_sub(ddr_make(fixed, fixedl), GiD(net, L, 0, S_CAST(int, m)));
    const double tv = fixed - Gi(net, L, 0, S_CAST(int, m));
    const double w = ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_add(ddr_make(L->pmax, L->pmaxl), tvD), ddr_make(L->vobs, L->vobsl)));
    const uint64_t leafh = fixedh + net->FH[m] + net->FH[L->r[0] - m];
    const double lim_def = L->vobs - kTieEps - tv;
    if (L->leafrec) {
      LeafRec lr;
      lr.q = lim_def;
      lr.w = w;
      lr.h = leafh;
      lr.tv = tv;
      lr.tvl = tvD.x[1];
      if (LeafVecPush(net, L->leafrec, lr)) {
        return 0;
      }
    } else {
      const size_t lo = UpperBoundV(L->srt, 0, L->n, lim_def);
      if (lo) {
        AccAdd(L->acc, w * (L->cum[lo] - L->cum0));
      }
      if ((lo < L->n) && (L->srt[lo].v <= L->thr - tv)) {
        ResolveWindow(net, L, lo, L->thr - tv, tvD, leafh, w);
      }
    }
    return L->pmax + tv <= L->vobs - kTieEps;
  }
  int md[kFisherNetMaxRows + 1];
  const double gmax = GroupMax(net, L, j, m, md);
  if (L->pmax + fixed + gmax <= L->vobs - kTieEps) {
    AccAdd(L->acc, ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_add(ddr_sloppier_add(ddr_make(L->pmax, L->pmaxl), ddr_make(fixed, fixedl)), LogTotD(net, L, j, m)), ddr_make(L->vobs, L->vobsl))) * L->cum_all);
    return 1;
  }
  if ((j != 2) && (L->pmin + fixed + GroupMin(net, L, j, m) > L->thr)) {
    return 0;  // nothing counts
  }
  const int x0 = md[j - 1];
  if (j == 2) {
    Level2Fast(net, L, m, fixed, fixedl, fixedh, x0);
    return 0;
  }
  const int xlo = S_CAST(int, MAXV(0LL, m - L->S[j - 1]));
  const int xhi = S_CAST(int, MINV(S_CAST(long long, L->r[j - 1]), m));
  for (int dir = 1; dir >= -1; dir -= 2) {
    for (int x = (dir > 0)? x0 : (x0 - 1); (x >= xlo) && (x <= xhi); x += dir) {
      const dd_real fD = ddr_sloppier_sub(ddr_make(fixed, fixedl), GiD(net, L, j - 1, x));
      const int full = LastLevel(net, L, j - 1, m - x, fD.x[0], fD.x[1], fixedh + net->FH[x] + net->FH[L->r[j - 1] - x]);
      if (net->err) {
        return 0;
      }
      if (full) {
        LastTail(net, L, j, m, fixed, fixedl, x + dir, dir);
        break;
      }
    }
  }
  return 0;
}

// ---- merging runs into a daughter's list ----
static inline size_t Slot0(uint64_t h, int shift) {
  return S_CAST(size_t, ((h * 0x9E3779B97F4A7C15ULL) >> shift));
}
// std::merge by v (ties take the first range's entry)
static void MergeEntries(const Entry* a, const Entry* a_end, const Entry* b, const Entry* b_end, Entry* out) {
  while ((a != a_end) && (b != b_end)) {
    if (b->v < a->v) {
      *out++ = *b++;
    } else {
      *out++ = *a++;
    }
  }
  if (a != a_end) {
    memcpy(out, a, S_CAST(size_t, a_end - a) * sizeof(Entry));
    out += a_end - a;
  }
  if (b != b_end) {
    memcpy(out, b, S_CAST(size_t, b_end - b) * sizeof(Entry));
  }
}
#ifdef FISHER_DEBUG_OVERRIDES
// merge segments [i, j) (bounds seg[]); result goes to dst if to_dst else
// stays in src (ping-pong)
static void MergeSegW(const uint32_t* seg, size_t i, size_t j, Entry* src, Entry* dst, uint32_t to_dst) {
  if (j == i + 1) {
    if (to_dst) {
      memcpy(&dst[seg[i]], &src[seg[i]], (seg[j] - seg[i]) * sizeof(Entry));
    }
    return;
  }
  // split where the cumulative size is closest to half, so a big segment is
  // not re-copied at every level
  const uint32_t half = seg[i] + (seg[j] - seg[i]) / 2;
  size_t k = i + 1;
  while ((k < j) && (seg[k] <= half)) {
    ++k;
  }
  if ((k > i + 1) && (half - seg[k - 1] < seg[k] - half)) {
    --k;
  }
  if (k <= i) {
    k = i + 1;
  }
  if (k >= j) {
    k = j - 1;
  }
  MergeSegW(seg, i, k, src, dst, !to_dst);
  MergeSegW(seg, k, j, src, dst, !to_dst);
  Entry* from = to_dst? src : dst;
  Entry* into = to_dst? dst : src;
  MergeEntries(&from[seg[i]], &from[seg[k]], &from[seg[k]], &from[seg[j]], &into[seg[i]]);
}
#endif
// bucket sort net->mbuf by v into out: uniform buckets over [vmin, vmax], one
// histogram pass, one scatter, small per-bucket sorts
static BoolErr BucketSortV(Net* net, Entry* out) {
  const Entry* in = net->mbuf.p;
  const size_t n = net->mbuf.n;
  if (Poll(net, n)) {
    return 1;
  }
  double vmin = in[0].v;
  double vmax = in[0].v;
  for (size_t q = 0; q < n; ++q) {
    vmin = MINV(vmin, in[q].v);
    vmax = MAXV(vmax, in[q].v);
  }
  if (!(vmax > vmin)) {
    memcpy(out, in, n * sizeof(Entry));
    return 0;
  }
  const size_t nb = n;
  const double inv = S_CAST(double, nb) / (vmax - vmin);
  if (U32VecAssign(net, &net->bcnt, nb + 2) || U32VecAssign(net, &net->bsort_start, nb + 2)) {
    return 1;
  }
  uint32_t* bcnt = net->bcnt.p;
  uint32_t* start = net->bsort_start.p;
  memset(bcnt, 0, (nb + 2) * sizeof(uint32_t));
  for (size_t q = 0; q < n; ++q) {
    const double bq = (in[q].v - vmin) * inv;
    ++bcnt[((bq >= S_CAST(double, nb))? nb : S_CAST(size_t, bq)) + 1];
  }
  for (size_t b = 0; b <= nb; ++b) {
    bcnt[b + 1] += bcnt[b];
  }
  memcpy(start, bcnt, (nb + 2) * sizeof(uint32_t));
  for (size_t q = 0; q < n; ++q) {
    const double bq = (in[q].v - vmin) * inv;
    out[start[(bq >= S_CAST(double, nb))? nb : S_CAST(size_t, bq)]++] = in[q];
  }
  for (size_t b = 0; b <= nb; ++b) {
    Entry* lo = &out[bcnt[b]];
    Entry* hi = &out[bcnt[b + 1]];
    if (hi - lo <= 1) {
      continue;
    }
    if (hi - lo <= 16) {
      for (Entry* p = lo + 1; p < hi; ++p) {
        const Entry e = *p;
        Entry* q = p;
        while ((q > lo) && ((q - 1)->v > e.v)) {
          *q = *(q - 1);
          --q;
        }
        *q = e;
      }
    } else {
      STD_SORT(hi - lo, CmpEntryV, lo);
    }
  }
  return 0;
}
// sorts net->mbuf (a concatenation of v-sorted segments with bounds net->seg)
// by v into out; mbuf is clobbered
static BoolErr SortSegments(Net* net, Entry* out) {
  EntryVec* buf = &net->mbuf;
  if (Poll(net, buf->n)) {
    return 1;
  }
  // drop empty segments; seg holds starts, plus buf->n at the end
  {
    uint32_t* seg = net->seg.p;
    size_t w = 0;
    for (size_t i = 0; i + 1 < net->seg.n; ++i) {
      if (seg[i + 1] > seg[i]) {
        seg[w++] = seg[i];
      }
    }
    seg[w++] = S_CAST(uint32_t, buf->n);
    net->seg.n = w;
    // coalesce neighbours that are already in order (last of i <= first of i+1)
    size_t w2 = 1;
    for (size_t i = 1; i + 1 < net->seg.n; ++i) {
      if (buf->p[seg[i] - 1].v > buf->p[seg[i]].v) {
        seg[w2++] = seg[i];
      }
    }
    seg[w2++] = S_CAST(uint32_t, buf->n);
    net->seg.n = w2;
  }
#ifdef FISHER_DEBUG_OVERRIDES
  if ((net->seg.n > 2) && (buf->n < g_seg_t * (net->seg.n - 1))) {
    STD_SORT(buf->n, CmpEntryV, buf->p);
    memcpy(out, buf->p, buf->n * sizeof(Entry));
    return 0;
  }
#endif
  if (net->seg.n <= 2) {
    if (buf->n) {
      memcpy(out, buf->p, buf->n * sizeof(Entry));
    }
    return 0;
  }
  if (g_bsort && (buf->n >= g_bsort_min) && (net->seg.n > 3)) {
    return BucketSortV(net, out);
  }
  Entry* src = buf->p;
  Entry* dst = out;  // ping-pong between buf and out
#ifdef FISHER_DEBUG_OVERRIDES
  if (g_seg_w) {
    MergeSegW(net->seg.p, 0, net->seg.n - 1, src, dst, 1);
    return 0;
  }
#endif
  while (net->seg.n > 2) {
    const uint32_t* seg = net->seg.p;
    net->seg2.n = 0;
    size_t i = 0;
    for (; i + 2 < net->seg.n; i += 2) {
      if (U32VecPush(net, &net->seg2, seg[i])) {
        return 1;
      }
      MergeEntries(&src[seg[i]], &src[seg[i + 1]], &src[seg[i + 1]], &src[seg[i + 2]], &dst[seg[i]]);
    }
    if (i + 1 < net->seg.n) {
      if (U32VecPush(net, &net->seg2, seg[i])) {
        return 1;
      }
      memcpy(&dst[seg[i]], &src[seg[i]], (seg[i + 1] - seg[i]) * sizeof(Entry));
    }
    if (U32VecPush(net, &net->seg2, S_CAST(uint32_t, buf->n))) {
      return 1;
    }
    const U32Vec tmp = net->seg;
    net->seg = net->seg2;
    net->seg2 = tmp;
    Entry* t = src;
    src = dst;
    dst = t;
  }
  if (src != out) {
    memcpy(out, src, buf->n * sizeof(Entry));
  }
  return 0;
}

static inline uint32_t CntIns(uint64_t* tab, size_t mask, int shift, uint64_t h) {
  const uint64_t hk = h? h : 1;
  size_t p = Slot0(hk, shift);
  while (tab[p] && (tab[p] != hk)) {
    p = (p + 1) & mask;
  }
  if (tab[p]) {
    return 0;
  }
  tab[p] = hk;
  return 1;
}
// number of distinct new hashes the runs would add to the (distinct) list
// `have`; stops once it exceeds `cap`.  On failure, net->err is set.
static size_t CountNewDistinct(Net* net, const NodeVec* mothers, const Run* runs, size_t nr, const EntryList* have, size_t cap) {
  const size_t want = 4 * (have->n + MINV(cap, S_CAST(size_t, 1) << 26) + 16);
  int lg = 4;
  while ((S_CAST(size_t, 1) << lg) < want) {
    ++lg;
  }
  const size_t mask = (S_CAST(size_t, 1) << lg) - 1;
  const int shift = 64 - lg;
  // 0 = empty (a factorial-product hash of exactly 0 is treated as a
  // collision-free miss)
  if (U64VecAssign(net, &net->cnt_tab, mask + 1)) {
    return 0;
  }
  uint64_t* tab = net->cnt_tab.p;
  memset(tab, 0, (mask + 1) * sizeof(uint64_t));
  for (size_t q = 0; q < have->n; ++q) {
    CntIns(tab, mask, shift, have->p[q].h);
  }
  size_t added = 0;
  for (size_t i = 0; (i < nr) && (added <= cap); ++i) {
    const Run* R = &runs[i];
    const Entry* me = mothers->p[R->m].ents.p;
    if (Poll(net, R->hi - R->lo)) {
      break;
    }
    for (uint32_t q = R->lo; q < R->hi; ++q) {
      if (CntIns(tab, mask, shift, me[q].h + R->ch) && (++added > cap)) {
        break;
      }
    }
  }
  U64VecFree(&net->cnt_tab);
  return added;
}
// Merge v-sorted runs (each a slice of a mother list shifted by colv) into a
// v-sorted list, merging entries with identical factorial-product hash h
// (equal h => v equal up to rounding).  Slot index = Fibonacci hash (one
// multiply, top bits).  out is a concatenation of v-increasing segments (one
// per run: a run's new entries arrive in v order, and merges never move an
// entry), which SortSegments() merges.
static BoolErr MergeRuns(Net* net, const NodeVec* mothers, const Run* runs, size_t nr, EntryList* out, uint32_t append) {
  if ((!append) || (!out->n)) {
    append = 0;
    out->p = nullptr;
    out->n = 0;
  }
  if (nr == 0) {
    return 0;
  }
  if ((nr == 1) && (!append)) {
    // single run: already distinct and sorted
    const Run* R = &runs[0];
    const Entry* me = mothers->p[R->m].ents.p;
    Entry* dst = ArenaAlloc(net, net->out_arena, R->hi - R->lo);
    if (net->err) {
      return 1;
    }
    for (uint32_t q = R->lo; q < R->hi; ++q) {
      const dd_real v = ddr_sloppier_add(ddr_make(me[q].v, me[q].vl), ddr_make(R->colv, R->colvl));
      Entry* d = &dst[q - R->lo];
      d->v = v.x[0];
      d->h = me[q].h + R->ch;
      d->f = me[q].f;
      d->vl = v.x[1];
    }
    out->p = dst;
    out->n = R->hi - R->lo;
    return 0;
  }
  // Build in a reusable scratch buffer; the final sort pass writes into the arena.
  EntryVec* buf = &net->mbuf;
  buf->n = 0;
  if (append) {
    if (EntryVecAssign(net, buf, out->n)) {
      return 1;
    }
    memcpy(buf->p, out->p, out->n * sizeof(Entry));
    ArenaFree(net->out_arena, out->p, out->n);
  }
  size_t mx = buf->n;
  for (size_t i = 0; i < nr; ++i) {
    mx = MAXV(mx, S_CAST(size_t, runs[i].hi - runs[i].lo));
  }
  int lg = 4;
  while ((S_CAST(size_t, 1) << lg) < 4 * mx) {
    ++lg;
  }
  size_t cap = S_CAST(size_t, 1) << lg;
  size_t mask = cap - 1;
  int shift = 64 - lg;
  if (SlotVecAssign(net, &net->tab, cap)) {
    return 1;
  }
  Slot* tab = net->tab.p;
  for (size_t z = 0; z < cap; ++z) {
    tab[z].h = 0;
    tab[z].ix = -1;
  }
  net->seg.n = 0;
  if (append) {
    if (U32VecPush(net, &net->seg, 0)) {
      return 1;
    }
    for (size_t z = 0; z < buf->n; ++z) {
      size_t pp = Slot0(buf->p[z].h, shift);
      while (tab[pp].ix >= 0) {
        pp = (pp + 1) & mask;
      }
      tab[pp].h = buf->p[z].h;
      tab[pp].ix = S_CAST(int32_t, z);
    }
  }
  const double eps = net->kTieEpsDD;
  uint32_t collision = 0;
  for (size_t i = 0; i < nr; ++i) {
    const Run* R = &runs[i];
    const Entry* me = mothers->p[R->m].ents.p;
    if (Poll(net, R->hi - R->lo) || U32VecPush(net, &net->seg, S_CAST(uint32_t, buf->n))) {
      return 1;
    }
    const double colv = R->colv;
    const double colvl = R->colvl;
    const uint64_t ch = R->ch;
    for (uint32_t q = R->lo; q < R->hi; ++q) {
      const uint64_t h = me[q].h + ch;
      size_t p = Slot0(h, shift);
      while ((tab[p].ix >= 0) && (tab[p].h != h)) {
        p = (p + 1) & mask;
      }
      if (tab[p].ix >= 0) {
        // Equal hashes mean equal factorial products unless the hash collided,
        // so the values are checked too: a collision is detected (and the call
        // falls back) unless the values also agree to double-double accuracy.
        // sv.x[0] - o->v has only relative rounding error (and is exact when
        // they're within a factor of 2), so a tiny result is trustworthy.  The
        // failure is deferred to the end of the merge, which keeps it out of
        // the hot loop.
        Entry* o = &buf->p[tab[p].ix];
        const dd_real sv = ddr_add2d(me[q].v, colv);
        collision |= (fabs((sv.x[0] - o->v) + (((sv.x[1] + me[q].vl) + colvl) - o->vl)) > eps);
        o->f += me[q].f;
        continue;
      }
      tab[p].h = h;
      tab[p].ix = S_CAST(int32_t, buf->n);
      {
        const dd_real v = ddr_sloppier_add(ddr_make(me[q].v, me[q].vl), ddr_make(colv, colvl));
        Entry e;
        e.v = v.x[0];
        e.h = h;
        e.f = me[q].f;
        e.vl = v.x[1];
        if (EntryVecPush(net, buf, e)) {
          return 1;
        }
      }
      if (2 * buf->n > cap) {
        // grow and rehash
        ++lg;
        cap <<= 1;
        mask = cap - 1;
        shift = 64 - lg;
        if (SlotVecAssign(net, &net->tab, cap)) {
          return 1;
        }
        tab = net->tab.p;
        for (size_t z = 0; z < cap; ++z) {
          tab[z].h = 0;
          tab[z].ix = -1;
        }
        for (size_t z = 0; z < buf->n; ++z) {
          size_t pp = Slot0(buf->p[z].h, shift);
          while (tab[pp].ix >= 0) {
            pp = (pp + 1) & mask;
          }
          tab[pp].h = buf->p[z].h;
          tab[pp].ix = S_CAST(int32_t, z);
        }
      }
    }
  }
  if (collision) {
    return NetFail(net, kFisherNetUnsupported);
  }
  if (U32VecPush(net, &net->seg, S_CAST(uint32_t, buf->n))) {
    return 1;
  }
  Entry* dst = ArenaAlloc(net, net->out_arena, buf->n);
  if (net->err || SortSegments(net, dst)) {
    return 1;
  }
  out->p = dst;
  out->n = buf->n;
  return 0;
}

// incremental merge: when run descriptors take too much of the budget, merge
// nodes whose runs are short (a 32-byte run descriptor covering ~1 entry is
// bigger than the 24-byte entry it merges into)
static BoolErr Flush(Net* net) {
  // largest run sets first; stop once runs are back under 1/8 of the budget
  FlushOrdVec* ord = &net->ord;
  ord->n = 0;
  for (size_t d = 0; d < net->druns.n; ++d) {
    if (net->druns.p[d].n) {
      FlushOrd od;
      od.cap = net->druns.p[d].cap;
      od.d = d;
      if (FlushOrdVecPush(net, ord, od)) {
        return 1;
      }
    }
  }
  if (ord->n > 1) {
    STD_SORT(ord->n, CmpFlushOrd, ord->p);
  }
  for (size_t oi = 0; oi < ord->n; ++oi) {
    if (net->mem_runs <= net->mem_limit / 8) {
      break;
    }
    const size_t d = ord->p[oi].d;
    RunVec* rv = &net->druns.p[d];
    size_t tot = 0;
    for (size_t q = 0; q < rv->n; ++q) {
      tot += rv->p[q].hi - rv->p[q].lo;
    }
    const size_t run_bytes = rv->cap * sizeof(Run);
    if (tot * sizeof(Entry) > run_bytes) {
      // worst case would grow memory: count the new distinct entries first
      // (early exit), merge only if it shrinks
      const size_t budget = run_bytes / sizeof(Entry);
      const size_t nd = CountNewDistinct(net, &net->cur, rv->p, rv->n, &net->nxt.p[d].ents, budget);
      if (net->err) {
        return 1;
      }
      if (nd > budget) {
        continue;
      }
    }
    if (MergeRuns(net, &net->cur, rv->p, rv->n, &net->nxt.p[d].ents, 1)) {
      return 1;
    }
    net->mem_new = net->arena_nxt.bytes;
    net->mem_runs -= rv->cap * sizeof(Run);
    RunVecRelease(net, rv);
  }
  // whatever could not be merged waits until runs grow by half again
  net->flush_at = net->mem_runs + net->mem_runs / 2;
  return MemCheck(net, 0);
}

// ---- daughter enumeration (non-last stages) ----
typedef struct StageCtxStruct {
  Net* net;
  Acc* acc;
  const int* cols;  // column totals in processing order
  int nr, s, rem_after, bits;
  uint32_t packed_keys;
  double vobs, vobsl, thr;
  // current mother
  int mi;
  const Entry* srt;
  size_t nsrt;
  const double* cum;
  double pmax, v0;
  dd_real pmaxD;
  uint32_t single;
  const int* r;
  int x[kFisherNetMaxRows];
  // grouped enumeration
  const int* crem;
  int mcols;
  dd_real KD;
  long long S[kFisherNetMaxRows + 1];
  dd_real slf[kFisherNetMaxRows + 1];
  double cumall;
  int* gub_cc;   // GroupUB() scratch: the later columns plus one more
} StageCtx;

// Node keys: the sorted remaining row totals, packed exactly when they fit in
// 63 bits, and otherwise hashed (a lookup then checks the stored totals, and a
// genuine collision makes the call return kFisherNetUnsupported).
static inline uint64_t NodeKey(const StageCtx* sc, const int* dr) {
  const int nr = sc->nr;
  if (sc->packed_keys) {
    uint64_t key = 1ULL << 63;
    for (int q = 0; q < nr; ++q) {
      key |= S_CAST(uint64_t, dr[q]) << (q * sc->bits);
    }
    return key;
  }
  uint64_t h = 0x9E3779B97F4A7C15ULL;
  for (int q = 0; q < nr; ++q) {
    h = Mix(h ^ (S_CAST(uint64_t, dr[q]) + 0x632BE59BD9B4E019ULL * S_CAST(uint64_t, q + 1)));
  }
  return h | 1;
}
// Finds or creates the daughter node with sorted row totals dr; returns its
// index, or -1 on failure.
static int GetDaughter(StageCtx* sc, const int* dr) {
  Net* net = sc->net;
  const int nr = sc->nr;
  int* slot = MapAt(net, &net->nodemap, NodeKey(sc, dr));
  if (!slot) {
    return -1;
  }
  if (*slot < 0) {
    const int idx = S_CAST(int, net->nxt.n);
    Node nd;
    memset(&nd, 0, sizeof(Node));
    RunVec rv;
    rv.p = nullptr;
    rv.n = 0;
    rv.cap = 0;
    if (NodeVecPush(net, &net->nxt, nd) || RunVecVecPush(net, &net->druns, rv)) {
      return -1;
    }
    *slot = idx;
    Node* dn = &net->nxt.p[idx];
    memcpy(dn->r, dr, sizeof(int) * nr);
    net->mem_nodes = net->nxt.cap * sizeof(Node) + net->druns.cap * sizeof(RunVec);
    if (BoundsLazy(net, dn, &sc->cols[sc->s + 1], sc->rem_after)) {
      return -1;
    }
    return idx;
  }
  if ((!sc->packed_keys) && memcmp(net->nxt.p[*slot].r, dr, sizeof(int) * nr)) {
    NetFail(net, kFisherNetUnsupported);
    return -1;
  }
  return *slot;
}
static inline void SortRows(int* dr, int nr) {
  for (int a = 1; a < nr; ++a) {
    const int v = dr[a];
    int b = a;
    while ((b > 0) && (dr[b - 1] > v)) {
      dr[b] = dr[b - 1];
      --b;
    }
    dr[b] = v;
  }
}
// Processes the daughter column in sc->x; returns 1 iff the whole mother list
// counts (or on failure, with net->err set).
static int Body(StageCtx* sc) {
  Net* net = sc->net;
  const int nr = sc->nr;
  const int* r = sc->r;
  const int* x = sc->x;
  const double kTieEps = net->kTieEps;
  dd_real colvD = ddr_make(0., 0.);
  uint64_t ch = 0;
  int dr[kFisherNetMaxRows];
  for (int q = 0; q < nr; ++q) {
    colvD = ddr_sloppier_sub(colvD, DLF(net, x[q]));
    ch += net->FH[x[q]];
    dr[q] = r[q] - x[q];
  }
  const double colv = colvD.x[0];
  if (Poll(net, 1)) {
    return 1;
  }
  SortRows(dr, nr);
  const int di = GetDaughter(sc, dr);
  if (di < 0) {
    return 1;
  }
  Node* dn = &net->nxt.p[di];
  const Entry* srt = sc->srt;
  const size_t nsrt = sc->nsrt;
  const double vobs = sc->vobs;
  const double thr = sc->thr;
  size_t pre;
  size_t keep_end;
  if (sc->pmax + colv + dn->ub <= vobs - kTieEps) {
    pre = nsrt;
    keep_end = nsrt;
  } else if (srt[0].v + colv + dn->splb > thr) {
    pre = 0;
    keep_end = 0;
  } else {
    if (!dn->have_lp) {
      dn->lp = MaxCompletion(net, dn->r, &sc->cols[sc->s + 1], sc->rem_after, nullptr, nullptr);
      if (net->err) {
        return 1;
      }
      dn->have_lp = 1;
    }
    if (sc->single) {
      pre = (sc->v0 <= vobs - kTieEps - colv - dn->lp);
      keep_end = (sc->v0 <= thr - colv - dn->splb);
    } else {
      pre = UpperBoundV(srt, 0, nsrt, vobs - kTieEps - colv - dn->lp);
      keep_end = UpperBoundV(srt, 0, nsrt, thr - colv - dn->splb);
    }
  }
  if (pre) {
    AccAdd(sc->acc, ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_add(ddr_sloppier_add(sc->pmaxD, colvD), ddr_make(dn->logt, dn->logtl)), ddr_make(vobs, sc->vobsl))) * sc->cum[pre]);
  }
  if (keep_end > pre) {
    RunVec* rv = &net->druns.p[di];
    const size_t c0 = rv->cap;
    Run run;
    run.m = sc->mi;
    run.lo = S_CAST(uint32_t, pre);
    run.hi = S_CAST(uint32_t, keep_end);
    run.colv = colv;
    run.ch = ch;
    run.colvl = colvD.x[1];
    if (RunVecPush(net, rv, &run)) {
      return 1;
    }
    net->mem_runs += (rv->cap - c0) * sizeof(Run);
    if (!(++net->push_ctr & 0x7ff)) {
      if ((net->mem_limit != SIZE_MAX) && (net->mem_runs > MAXV(net->mem_limit / 4, net->flush_at))) {
        if (Flush(net)) {
          return 1;
        }
      }
      if (MemCheck(net, 0)) {
        return 1;
      }
    }
  }
  return pre == nsrt;
}
// exact best completion value colv(x) + lp(r - x) of the daughter in sc->x
// (creates its node if needed); on failure, net->err is set
static double ExactF(StageCtx* sc) {
  Net* net = sc->net;
  const int nr = sc->nr;
  double colv = 0.;
  int dr[kFisherNetMaxRows];
  for (int q = 0; q < nr; ++q) {
    colv -= net->LF[sc->x[q]];
    dr[q] = sc->r[q] - sc->x[q];
  }
  SortRows(dr, nr);
  const int di = GetDaughter(sc, dr);
  if (di < 0) {
    return 0.;
  }
  Node* dn = &net->nxt.p[di];
  if (!dn->have_lp) {
    dn->lp = MaxCompletion(net, dn->r, &sc->cols[sc->s + 1], sc->rem_after, nullptr, nullptr);
    if (net->err) {
      return 0.;
    }
    dn->have_lp = 1;
  }
  return colv + dn->lp;
}
// grouped enumeration bound: rows j..nr-1 fixed (sc->x), rows 0..j-1 free with
// sum mrem
static double GroupUB(const StageCtx* sc, int j, int mrem, double fixedv) {
  const double* LF = sc->net->LF;
  const int nr = sc->nr;
  const int* r = sc->r;
  const int* x = sc->x;
  const int* crem = sc->crem;
  const int mcols = sc->mcols;
  // column-wise: this column (fixed cells + even split of mrem over free rows) + each later column
  int caps[kFisherNetMaxRows];
  for (int q = 0; q < nr; ++q) {
    caps[q] = (q < j)? r[q] : (r[q] - x[q]);
  }
  SortRows(caps, nr);
  double cw = fixedv + EvenSplitMaxSorted(LF, r, j, mrem);
  for (int l = 0; l < mcols; ++l) {
    cw += EvenSplitMaxSorted(LF, caps, nr, crem[l]);
  }
  // row-wise: fixed rows split r_i - x_i over later columns; free rows split r_i over {mrem} + later columns
  int* cc = sc->gub_cc;
  int nc = 0;
  int ins = 0;
  for (int l = 0; l < mcols; ++l) {
    if ((!ins) && (mrem <= crem[l])) {
      cc[nc++] = mrem;
      ins = 1;
    }
    cc[nc++] = crem[l];
  }
  if (!ins) {
    cc[nc++] = mrem;
  }
  double rw = fixedv;
  for (int q = j; q < nr; ++q) {
    rw += EvenSplitMaxSorted(LF, crem, mcols, r[q] - x[q]);
  }
  for (int q = 0; q < j; ++q) {
    rw += EvenSplitMaxSorted(LF, cc, nc, r[q]);
  }
  return MINV(cw, rw);
}
static inline int Visit(StageCtx* sc, int mrem, int t) {
  sc->x[1] = t;
  sc->x[0] = mrem - t;
  return Body(sc);
}
static inline double Fval(StageCtx* sc, int mrem, int t) {
  sc->x[1] = t;
  sc->x[0] = mrem - t;
  return ExactF(sc);
}
// grouped enumeration: rows j..nr-1 fixed, rows 0..j-1 free with sum mrem
static void Rec(StageCtx* sc, int j, int mrem, double fixedv, dd_real fixedt) {
  Net* net = sc->net;
  int* x = sc->x;
  const int* r = sc->r;
  const double vobs = sc->vobs;
  if (j == 1) {
    x[0] = mrem;
    Body(sc);
    return;
  }
  if (sc->pmax + GroupUB(sc, j, mrem, fixedv) <= vobs - net->kTieEps) {
    const dd_real lt = ddr_sloppier_sub(ddr_sloppier_sub(ddr_sloppier_sub(DLF(net, sc->S[j]), DLF(net, mrem)), DLF(net, sc->S[j] - mrem)), sc->slf[j]);
    AccAdd(sc->acc, ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_add(ddr_sloppier_add(ddr_sloppier_add(sc->pmaxD, sc->KD), fixedt), lt), ddr_make(vobs, sc->vobsl))) * sc->cumall);
    return;
  }
  const int lo = S_CAST(int, MAXV(0LL, mrem - sc->S[j - 1]));
  const int hi2 = MINV(r[j - 1], mrem);
  if ((j == 2) && g_interval) {
    // f(t) = colv + lp for x1 = t, x0 = mrem - t is concave, so the daughters
    // that are not fully counted form one interval of t.  Visit a block [a, b]
    // containing it; sum the rest by recurrence.
    const int r0 = r[0];
    const int r1 = r[1];
    int t0 = S_CAST(int, (S_CAST(long long, r1) * mrem / MAXV(1LL, S_CAST(long long, r0) + r1)));
    t0 = MAXV(lo, MINV(hi2, t0));
    int a = t0;
    int b = t0;
    int found = !Visit(sc, mrem, t0);
    if (net->err) {
      return;
    }
    if (!found) {
      // climb toward the maximum of f; stop at a non-full daughter or at the
      // top (then all are full)
      const double f0 = Fval(sc, mrem, t0);
      int dir = 0;
      double fn = 0.;
      if ((t0 + 1 <= hi2) && ((fn = Fval(sc, mrem, t0 + 1)) > f0 + 1e-12)) {
        dir = 1;
      } else if ((t0 - 1 >= lo) && ((fn = Fval(sc, mrem, t0 - 1)) > f0 + 1e-12)) {
        dir = -1;
      }
      if (net->err) {
        return;
      }
      int t = t0;
      while (dir) {
        t += dir;
        if (dir > 0) {
          b = t;
        } else {
          a = t;
        }
        const int full = Visit(sc, mrem, t);
        if (net->err) {
          return;
        }
        if (!full) {
          found = 1;
          break;
        }
        const int tn = t + dir;
        if ((tn < lo) || (tn > hi2)) {
          break;
        }
        const double ft = fn;
        fn = Fval(sc, mrem, tn);
        if (net->err) {
          return;
        }
        if (!(fn > ft + 1e-12)) {
          break;
        }
      }
    }
    if (found) {
      while (b + 1 <= hi2) {
        ++b;
        const int full = Visit(sc, mrem, b);
        if (net->err) {
          return;
        }
        if (full) {
          break;
        }
      }
      while (a - 1 >= lo) {
        --a;
        const int full = Visit(sc, mrem, a);
        if (net->err) {
          return;
        }
        if (full) {
          break;
        }
      }
    }
    // tails [lo, a-1] and [b+1, hi2]: all fully counted
    const dd_real base = ddr_sloppier_sub(ddr_sloppier_add(ddr_sloppier_add(sc->pmaxD, sc->KD), fixedt), ddr_make(vobs, sc->vobsl));
#ifdef FISHER_DEBUG_OVERRIDES
    if (g_gcheck) {
      // debug builds only: keeps exit() and stderr out of normal builds
      for (int t = lo; t < a; ++t) {
        if (!Visit(sc, mrem, t)) {
          fprintf(stderr, "GCHECK violation\n");
          exit(5);
        }
      }
      for (int t = b + 1; t <= hi2; ++t) {
        if (!Visit(sc, mrem, t)) {
          fprintf(stderr, "GCHECK violation\n");
          exit(5);
        }
      }
      return;
    }
#endif
    for (int dir = 1; dir >= -1; dir -= 2) {
      int t = (dir > 0)? (b + 1) : (a - 1);
      if ((t < lo) || (t > hi2)) {
        continue;
      }
      double w = ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_sub(ddr_sloppier_sub(ddr_sloppier_sub(base, DLF(net, t)), DLF(net, r1 - t)), DLF(net, mrem - t)), DLF(net, r0 - mrem + t)));
      double sum = 0.;
      int since = 0;
      for (;;) {
        sum += w;
        const int tn = t + dir;
        if ((tn < lo) || (tn > hi2)) {
          break;
        }
        const double ratio = (dir > 0)? (S_CAST(double, (r1 - t)) * (mrem - t) / (S_CAST(double, (t + 1)) * (r0 - mrem + t + 1))) : (S_CAST(double, t) * (r0 - mrem + t) / (S_CAST(double, (r1 - t + 1)) * (mrem - t + 1)));
        if ((ratio < 1.) && (w * ratio < 1e-17 * sum * (1. - ratio))) {
          break;
        }
        t = tn;
        if (++since >= 64) {
          w = ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_sub(ddr_sloppier_sub(ddr_sloppier_sub(base, DLF(net, t)), DLF(net, r1 - t)), DLF(net, mrem - t)), DLF(net, r0 - mrem + t)));
          since = 0;
        } else {
          w *= ratio;
        }
      }
      AccAdd(sc->acc, sum * sc->cumall);
    }
    return;
  }
  for (int t = lo; t <= hi2; ++t) {
    x[j - 1] = t;
    Rec(sc, j - 1, mrem - t, fixedv - net->LF[t], ddr_sloppier_sub(ddr_sloppier_sub(fixedt, DLF(net, t)), DLF(net, r[j - 1] - t)));
    if (net->err) {
      return;
    }
  }
}

// releases the containers that live for one stage
static void StageRelease(Net* net) {
  NodeVecRelease(&net->nxt);
  RunVecVecRelease(net, &net->druns);
  MapFree(&net->nodemap);
  DoubleVecFree(&net->cum);
  U32VecFree(&net->bidx);
}
static void SwapNodeVecs(NodeVec* a, NodeVec* b) {
  const NodeVec t = *a;
  *a = *b;
  *b = t;
}

// (entry, leaf) pair inside the tie window in the deferred last stage: exact
// resolution.  Returns the pair's weight factor (tie_w for a tie, 1 if the
// table counts, 0 if not); on failure, net->err is set.
static double ResolveDeferred(Net* net, const Entry* e, const Run* R, const LeafRec* lr, uint64_t obs_h, double vobs, double vobsl, double tie_w) {
  const dd_real diff = ddr_sloppier_sub(ddr_sloppier_add(ddr_sloppier_add(ddr_make(e->v, e->vl), ddr_make(R->colv, R->colvl)), ddr_make(lr->tv, lr->tvl)), ddr_make(vobs, vobsl));
  const double d = diff.x[0] + diff.x[1];
  // outside the margin, the sign of d is reliable whatever the hashes say
  if (fabs(d) > net->kTieEpsDD) {
    return (d < 0.)? 1. : 0.;
  }
  // inside the margin: a tie iff the factorial products are equal
  if (e->h + R->ch + lr->h != obs_h) {
    NetFail(net, kFisherNetUnsupported);
    return 0.;
  }
  return tie_w;
}
static inline size_t LeafBucket(double x, double qmin, double qinv, size_t qnb) {
  const double t = (x - qmin) * qinv;
  return (t >= S_CAST(double, qnb))? qnb : S_CAST(size_t, t);
}
// first leaf index whose threshold is >= x (std::lower_bound)
static inline size_t FirstLeafGe(const LeafRec* lv, size_t nl, double x) {
  size_t lo = 0;
  size_t len = nl;
  while (len) {
    const size_t half = len >> 1;
    if (lv[lo + half].q < x) {
      lo += half + 1;
      len -= half + 1;
    } else {
      len = half;
    }
  }
  return lo;
}
static inline int Lg2(size_t v) {
  int b = 1;
  while (v > 1) {
    v >>= 1;
    ++b;
  }
  return b;
}

#ifdef FISHER_DEBUG_OVERRIDES
static void ApplyDebugEnv(void) {
  if (getenv("NO_GSKIP")) {
    g_gskip = 0;
  }
  if (getenv("NO_BUCKET")) {
    g_bucket = 0;
  }
  if (getenv("DEFER")) {
    g_defer = atoi(getenv("DEFER"))? 1 : 0;
  }
  if (getenv("DEFER_AVG")) {
    g_defer_avg = atof(getenv("DEFER_AVG"));
  }
  if (getenv("BUCKET_MIN")) {
    g_bucket_min = atoi(getenv("BUCKET_MIN"));
  }
  if (getenv("SEG_T")) {
    g_seg_t = S_CAST(size_t, atol(getenv("SEG_T")));
  }
  if (getenv("SEGW")) {
    g_seg_w = 1;
  }
  if (getenv("NO_BSORT")) {
    g_bsort = 0;
  }
  if (getenv("BSORT_MIN")) {
    g_bsort_min = S_CAST(size_t, atol(getenv("BSORT_MIN")));
  }
  if (getenv("NO_INTERVAL")) {
    g_interval = 0;
  }
  if (getenv("GCHECK")) {
    g_gcheck = 1;
  }
  g_order_spec = getenv("NET_ORDER");
}
#endif

// Returns the log-p-value; on failure, net->err is set (and the return value
// is meaningless).
static double Core(Net* net, const int32_t* tin, int nr, int nc) {
  // orient so rows <= columns (nodes are keyed on row totals)
  const int transpose = (nr > nc);
  if (transpose) {
    const int tmp = nr;
    nr = nc;
    nc = tmp;
  }
  if (nr > kFisherNetMaxRows) {
    NetFail(net, kFisherNetUnsupported);
    return 0.;
  }
  if (nr <= 1) {
    return 0.;
  }
  net->NR = nr;
  const int k = nc;
  const size_t ib_ct = S_CAST(size_t, nr) * k + S_CAST(size_t, nr) + 3 * S_CAST(size_t, k) + 1;
  net->ibuf = S_CAST(int*, malloc(ib_ct * sizeof(int)));
  if (!net->ibuf) {
    NetFail(net, kFisherNetMemLimit);
    return 0.;
  }
  int* t = net->ibuf;
  int* row = &t[nr * k];
  int* col = &row[nr];
  int* cols = &col[k];
  int* gub_cc = &cols[k];
  if (transpose) {
    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < nr; ++j) {
        t[j * k + i] = tin[i * nr + j];
      }
    }
  } else {
    for (int q = 0; q < nr * k; ++q) {
      t[q] = tin[q];
    }
  }
  long long N = 0;
  for (int i = 0; i < nr; ++i) {
    row[i] = 0;
  }
  for (int j = 0; j < k; ++j) {
    col[j] = 0;
  }
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < k; ++j) {
      row[i] += t[i * k + j];
      col[j] += t[i * k + j];
      N += t[i * k + j];
    }
  }
  int bits = 1;
  while ((1LL << bits) <= N) {
    ++bits;
  }
  // log-factorial (high, low) and factorial-hash tables
  if (MemCheck(net, 3 * sizeof(double) * (S_CAST(size_t, N) + 2)) || InitTables(net, S_CAST(int, N))) {
    return 0.;
  }
  const double* LF = net->LF;
  dd_real vobsD = ddr_make(0., 0.);
  uint64_t obs_h = 0;
  for (int q = 0; q < nr * k; ++q) {
    vobsD = ddr_sloppier_sub(vobsD, DLF(net, t[q]));
    obs_h += net->FH[t[q]];
  }
  const double vobs = vobsD.x[0];
  const double vobsl = vobsD.x[1];
  dd_real cnD = ddr_negate(DLF(net, N));
  for (int i = 0; i < nr; ++i) {
    cnD = ddr_sloppier_add(cnD, DLF(net, row[i]));
  }
  for (int j = 0; j < k; ++j) {
    cnD = ddr_sloppier_add(cnD, DLF(net, col[j]));
  }
  // decisions compare double sums of up to ~nr*k log-factorials; their
  // rounding stays well inside this margin
  net->kTieEps = MAXV(1e-11, 16. * (nr * k + 8) * DBL_EPSILON * MAXV(1., LF[N]));
  const double kTieEps = net->kTieEps;
  // Inside that margin, a table's log-probability difference from the observed
  // table is recomputed in (unnormalized) double-double: a signed sum of at
  // most 2 * nr * k + 8 log-factorials (each at most LF[N]), plus partial sums
  // carried over from earlier stages.  ddr_lfact() is accurate to ~2^-98
  // relative, and each low-part addition rounds at ~2^-106 of the running
  // magnitude, so the error is far below 2^-96 * (2 * nr * k + 8) * LF[N];
  // this margin is 4x that (and the factor counting every term at full size
  // LF[N] is itself generous).  Differences inside it, between tables with
  // unequal factorial products (unequal hashes), make the call return
  // kFisherNetUnsupported.
  net->kTieEpsDD = 4. * (2 * nr * k + 8) * ldexp(1., -96) * MAXV(1., LF[N]);
  const double tie_w = net->midp? 0.5 : 1.;
  memcpy(cols, col, k * sizeof(int));
  STD_SORT(k, CmpInt, cols);
#ifdef FISHER_DEBUG_OVERRIDES
  if (g_order_spec) {
    const char* ord = g_order_spec;
    if (!strcmp(ord, "desc")) {
      for (int i = 0, j = k - 1; i < j; ++i, --j) {
        const int tmp = cols[i];
        cols[i] = cols[j];
        cols[j] = tmp;
      }
    } else if (strcmp(ord, "asc")) {
      // comma-separated permutation of indices into the ascending order
      memcpy(gub_cc, cols, k * sizeof(int));
      int cnt = 0;
      const char* p = ord;
      while (*p) {
        const int idx = atoi(p);
        if ((cnt >= k) || (idx < 0) || (idx >= k)) {
          NetFail(net, kFisherNetUnsupported);
          return 0.;
        }
        cols[cnt++] = gub_cc[idx];
        while (*p && (*p != ',')) {
          ++p;
        }
        if (*p) {
          ++p;
        }
      }
      if (cnt != k) {
        NetFail(net, kFisherNetUnsupported);
        return 0.;
      }
    }
  }
#endif
  const double thr = vobs + kTieEps;
  Acc acc;
  acc.s = 0.;
  acc.comp = 0.;
  {
    Node nd;
    memset(&nd, 0, sizeof(Node));
    memcpy(nd.r, row, nr * sizeof(int));
    SortRows(nd.r, nr);
    Entry* p = ArenaAlloc(net, &net->arena_cur, 1);
    if (!p) {
      return 0.;
    }
    p[0].v = 0.;
    p[0].h = 0;
    p[0].f = 1.;
    p[0].vl = 0.;
    nd.ents.p = p;
    nd.ents.n = 1;
    if (NodeVecPush(net, &net->cur, nd)) {
      return 0.;
    }
  }
  if ((k == 1) || (N == 0)) {
    return 0.;
  }
  // If the observed table is a mode, every table is at most as probable, so
  // the (non-mid) p-value is exactly 1; summing would leave rounding error.
  // Being within the decision margin of the mode's value isn't enough (a
  // slightly more probable table would be miscounted), so the observed table
  // must have the same factorial product as the computed mode table (equal
  // hashes, confirmed by value); otherwise the sum is computed as usual.
  if (!net->midp) {
    uint64_t mode_h;
    dd_real mode_v;
    const double mode_val = MaxCompletion(net, net->cur.p[0].r, cols, k, &mode_h, &mode_v);
    if (net->err) {
      return 0.;
    }
    if (vobs >= mode_val - kTieEps) {
      const dd_real d = ddr_sloppier_sub(mode_v, vobsD);
      if ((mode_h == obs_h) && (fabs(d.x[0] + d.x[1]) <= net->kTieEpsDD)) {
        return 0.;
      }
    }
  }
  uint32_t deferred = 0;
  const uint32_t defer_choice = (g_defer < 0)? (nr <= 3) : (g_defer > 0);
  StageCtx sc;
  memset(&sc, 0, sizeof(StageCtx));
  sc.net = net;
  sc.acc = &acc;
  sc.cols = cols;
  sc.nr = nr;
  sc.bits = bits;
  sc.packed_keys = (nr * bits <= 62);
  sc.vobs = vobs;
  sc.vobsl = vobsl;
  sc.thr = thr;
  sc.gub_cc = gub_cc;
  for (int s = 0; s + 1 < k; ++s) {
    const int c = cols[s];
    const int rem_after = k - 1 - s;
    const uint32_t defer_here = (rem_after == 2) && (defer_choice || (net->mem_limit != SIZE_MAX));
    net->mem_mothers = net->arena_cur.bytes + net->arena_prev.bytes;
    for (size_t mi = 0; mi < net->cur.n; ++mi) {
      net->mem_mothers += NodeBytes(&net->cur.p[mi]);
    }
    net->out_arena = &net->arena_nxt;
    net->mem_runs = 0;
    net->mem_nodes = 0;
    net->mem_new = 0;
    net->flush_at = 0;
    if (!net->pool.live) {
      RunPoolReset(&net->pool);
    }
    if (MemCheck(net, 0) || MapInit(net, &net->nodemap, 1024)) {
      return 0.;
    }
    sc.s = s;
    sc.rem_after = rem_after;
    for (size_t mi = 0; mi < net->cur.n; ++mi) {
      Node* mo = &net->cur.p[mi];
      if (deferred && (rem_after == 1) && net->prevruns.p[mi].n) {
        // ---- deferred last stage: node = union of runs (slices of stage-2
        // mother lists shifted by colv) ----
        const RunVec* runs = &net->prevruns.p[mi];
        double pmin_n = 1e300;
        double pmax_n = -1e300;
        dd_real pmaxnD = ddr_make(-1e300, 0.);
        for (size_t ri = 0; ri < runs->n; ++ri) {
          const Run* R = &runs->p[ri];
          const Node* M = &net->prevm.p[R->m];
          pmin_n = MINV(pmin_n, M->ents.p[R->lo].v + R->colv);
          const dd_real top = ddr_sloppier_add(ddr_make(M->ents.p[R->hi - 1].v, M->ents.p[R->hi - 1].vl), ddr_make(R->colv, R->colvl));
          if (top.x[0] > pmax_n) {
            pmax_n = top.x[0];
            pmaxnD = top;
          }
        }
        Acc callA;
        callA.s = 0.;
        callA.comp = 0.;
        for (size_t ri = 0; ri < runs->n; ++ri) {
          const Run* R = &runs->p[ri];
          const Node* M = &net->prevm.p[R->m];
          AccAdd(&callA, (M->cumv[R->hi] - M->cumv[R->lo]) * ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_add(ddr_make(M->pmaxv, M->pmaxvl), ddr_make(R->colv, R->colvl)), pmaxnD)));
        }
        const double call = AccGet(&callA);
        net->lvbuf.n = 0;
        LastCtx L;
        memset(&L, 0, sizeof(LastCtx));
        L.r = mo->r;
        L.srt = nullptr;
        L.n = 0;
        L.cum = nullptr;
        L.leafrec = &net->lvbuf;
        L.pmax = pmax_n;
        L.pmaxl = pmaxnD.x[1];
        L.pmin = pmin_n;
        L.cum_all = call;
        L.vobs = vobs;
        L.vobsl = vobsl;
        L.thr = thr;
        L.acc = &acc;
        L.obs_h = obs_h;
        L.tie_w = tie_w;
        L.S[0] = 0;
        L.sumlf[0] = 0.;
        L.sumlfl[0] = 0.;
        for (int i = 0; i < nr; ++i) {
          L.S[i + 1] = L.S[i] + mo->r[i];
          const dd_real sl = ddr_sloppier_add(ddr_make(L.sumlf[i], L.sumlfl[i]), DLF(net, mo->r[i]));
          L.sumlf[i + 1] = sl.x[0];
          L.sumlfl[i + 1] = sl.x[1];
        }
        if (c <= L.S[nr]) {
          LastLevel(net, &L, nr, c, 0., 0., 0);
          if (net->err) {
            return 0.;
          }
        }
        if (!net->lvbuf.n) {
          continue;
        }
        LeafRec* lvbuf = net->lvbuf.p;
        const size_t nl = net->lvbuf.n;
        STD_SORT(nl, CmpLeafQ, lvbuf);
        if (DoubleVecAssign(net, &net->wsuf, nl + 1)) {
          return 0.;
        }
        double* wsuf = net->wsuf.p;
        wsuf[nl] = 0.;
        for (size_t q = nl; q-- > 0; ) {
          wsuf[q] = wsuf[q + 1] + lvbuf[q].w;
        }
        const double win = 2. * kTieEps;
        // bucket index over the thresholds: qb[b] = first threshold index whose
        // bucket is >= b (bucket monotone in q)
        const double qmin = lvbuf[0].q;
        const double qmax = lvbuf[nl - 1].q;
        const size_t qnb = nl;
        const double qinv = (qmax > qmin)? (S_CAST(double, qnb) / (qmax - qmin)) : 0.;
        if (U32VecAssign(net, &net->qb, qnb + 2)) {
          return 0.;
        }
        uint32_t* qb = net->qb.p;
        memset(qb, 0, (qnb + 2) * sizeof(uint32_t));
        for (size_t q = 0; q < nl; ++q) {
          ++qb[LeafBucket(lvbuf[q].q, qmin, qinv, qnb) + 1];
        }
        for (size_t b = 0; b <= qnb; ++b) {
          qb[b + 1] += qb[b];
        }
        for (size_t ri = 0; ri < runs->n; ++ri) {
          const Run* R = &runs->p[ri];
          const Node* M = &net->prevm.p[R->m];
          const Entry* me = M->ents.p;
          const double colv = R->colv;
          const double xmin = me[R->lo].v + colv;
          const double xmax = me[R->hi - 1].v + colv;
          if (Poll(net, R->hi - R->lo)) {
            return 0.;
          }
          const size_t i0 = FirstLeafGe(lvbuf, nl, xmin - win);
          if (i0 == nl) {
            continue;  // every leaf threshold (with its tie window) is below this run
          }
          const size_t i1 = FirstLeafGe(lvbuf, nl, xmax);  // leaves i1.. count the whole run
          double sum = 0.;
          const size_t len = R->hi - R->lo;
          const size_t nthr = i1 - i0;
          if (S_CAST(double, len) * Lg2(nthr + 1) <= S_CAST(double, (nthr + 1)) * Lg2(len)) {
            // walk the run's entries, galloping in the thresholds
            for (uint32_t e = R->lo; e < R->hi; ++e) {
              const double xv = me[e].v + colv;
              if (xv > qmax + win) {
                break;
              }
              size_t idx;
              if (xv > qmax) {
                idx = nl;
              } else {
                // thresholds in lower buckets are all < xv
                idx = (xv <= qmin)? 0 : qb[LeafBucket(xv, qmin, qinv, qnb)];
                while (lvbuf[idx].q < xv) {
                  ++idx;
                }
              }
              sum += M->wv[e] * wsuf[idx];
              for (size_t l = idx; (l-- > 0) && (lvbuf[l].q >= xv - win); ) {
                const double f = ResolveDeferred(net, &me[e], R, &lvbuf[l], obs_h, vobs, vobsl, tie_w);
                if (f != 0.) {
                  sum += f * lvbuf[l].w * M->wv[e];
                }
              }
              if (net->err) {
                return 0.;
              }
            }
          } else {
            // walk the thresholds inside the run's range, galloping in the run
            size_t pos = R->lo;
            for (size_t l = i0; l < i1; ++l) {
              const double q = lvbuf[l].q;
              if ((pos < R->hi) && (me[pos].v + colv <= q)) {
                size_t step = 1;
                size_t prev = pos;
                while ((pos + step < R->hi) && (me[pos + step].v + colv <= q)) {
                  prev = pos + step;
                  step <<= 1;
                }
                size_t lo2 = prev + 1;
                size_t hi2 = MINV(pos + step, S_CAST(size_t, R->hi));
                while (lo2 < hi2) {
                  const size_t mid = (lo2 + hi2) / 2;
                  if (me[mid].v + colv <= q) {
                    lo2 = mid + 1;
                  } else {
                    hi2 = mid;
                  }
                }
                pos = lo2;
              }
              sum += lvbuf[l].w * (M->cumv[pos] - M->cumv[R->lo]);
              for (size_t e = pos; (e < R->hi) && (me[e].v + colv <= q + win); ++e) {
                const double f = ResolveDeferred(net, &me[e], R, &lvbuf[l], obs_h, vobs, vobsl, tie_w);
                if (f != 0.) {
                  sum += f * lvbuf[l].w * M->wv[e];
                }
              }
              if (net->err) {
                return 0.;
              }
            }
            sum += (M->cumv[R->hi] - M->cumv[R->lo]) * wsuf[i1];
            // leaves i1.. take the whole run, except entries in a leaf's window
            // (q < x <= q + win), which are resolved
            for (size_t l = i1; (l < nl) && (lvbuf[l].q < xmax); ++l) {
              for (uint32_t e = R->hi; (e-- > R->lo) && (me[e].v + colv > lvbuf[l].q); ) {
                const double f = ResolveDeferred(net, &me[e], R, &lvbuf[l], obs_h, vobs, vobsl, tie_w);
                sum -= (1. - f) * lvbuf[l].w * M->wv[e];
              }
              if (net->err) {
                return 0.;
              }
            }
          }
          AccAdd(&acc, sum * ddr_exp_to_double(ddr_sloppier_sub(ddr_sloppier_add(ddr_make(M->pmaxv, M->pmaxvl), ddr_make(R->colv, R->colvl)), pmaxnD)));
        }
        if (MemCheck(net, 0)) {
          return 0.;
        }
        continue;
      }
      if (!mo->ents.n) {
        continue;
      }
      const Entry* srt = mo->ents.p;
      const size_t nsrt = mo->ents.n;
      const double pmax = srt[nsrt - 1].v;
      const dd_real pmaxD = ddr_make(pmax, srt[nsrt - 1].vl);
      if (DoubleVecAssign(net, &net->cum, nsrt + 1)) {
        return 0.;
      }
      double* cum = net->cum.p;
      cum[0] = 0.;
      {
        Acc ca;
        ca.s = 0.;
        ca.comp = 0.;
        for (size_t q = 0; q < nsrt; ++q) {
          AccAdd(&ca, srt[q].f * ddr_exp_to_double(ddr_sloppier_sub(ddr_make(srt[q].v, srt[q].vl), pmaxD)));
          cum[q + 1] = AccGet(&ca);
        }
      }
      if (Poll(net, nsrt)) {
        return 0.;
      }
      if (defer_here) {
        mo->cumv = S_CAST(double*, BigAllocRaw((nsrt + 1) * sizeof(double)));
        if (!mo->cumv) {
          NetFail(net, kFisherNetMemLimit);
          return 0.;
        }
        mo->cumv_n = nsrt + 1;
        memcpy(mo->cumv, cum, (nsrt + 1) * sizeof(double));
        mo->pmaxv = pmax;
        mo->pmaxvl = pmaxD.x[1];
        mo->wv = S_CAST(double*, BigAllocRaw(nsrt * sizeof(double)));
        if (!mo->wv) {
          NetFail(net, kFisherNetMemLimit);
          return 0.;
        }
        mo->wv_n = nsrt;
        for (size_t q = 0; q < nsrt; ++q) {
          mo->wv[q] = srt[q].f * ddr_exp_to_double(ddr_sloppier_sub(ddr_make(srt[q].v, srt[q].vl), pmaxD));
        }
        net->mem_mothers += (mo->wv_n + mo->cumv_n) * sizeof(double);
        if (MemCheck(net, 0)) {
          return 0.;
        }
      }
      const int* r = mo->r;
      if (rem_after == 1) {
        LastCtx L;
        memset(&L, 0, sizeof(LastCtx));
        L.r = r;
        L.srt = srt;
        L.n = nsrt;
        L.cum = cum;
        L.pmax = pmax;
        L.pmaxl = pmaxD.x[1];
        L.pmin = srt[0].v;
        L.cum_all = cum[nsrt];
        L.vobs = vobs;
        L.vobsl = vobsl;
        L.thr = thr;
        L.acc = &acc;
        L.obs_h = obs_h;
        L.tie_w = tie_w;
        L.S[0] = 0;
        L.sumlf[0] = 0.;
        L.sumlfl[0] = 0.;
        for (int i = 0; i < nr; ++i) {
          L.S[i + 1] = L.S[i] + r[i];
          const dd_real sl = ddr_sloppier_add(ddr_make(L.sumlf[i], L.sumlfl[i]), DLF(net, r[i]));
          L.sumlf[i + 1] = sl.x[0];
          L.sumlfl[i + 1] = sl.x[1];
        }
        if (g_bucket && (nsrt >= S_CAST(size_t, g_bucket_min)) && (srt[nsrt - 1].v > srt[0].v)) {
          const int nb = S_CAST(int, nsrt);
          const double vmin = srt[0].v;
          const double inv = nb / (srt[nsrt - 1].v - vmin);
          if (U32VecAssign(net, &net->bidx, S_CAST(size_t, nb) + 2)) {
            return 0.;
          }
          uint32_t* bidx = net->bidx.p;
          memset(bidx, 0, (S_CAST(size_t, nb) + 2) * sizeof(uint32_t));
          for (size_t q = 0; q < nsrt; ++q) {
            const double bq = (srt[q].v - vmin) * inv;
            ++bidx[((bq >= S_CAST(double, nb))? nb : S_CAST(int, bq)) + 1];
          }
          for (int b = 0; b <= nb; ++b) {
            bidx[b + 1] += bidx[b];
          }
          L.bidx = bidx;
          L.bvmin = vmin;
          L.binv = inv;
          L.bnb = nb;
        }
        if (c <= L.S[nr]) {
          LastLevel(net, &L, nr, c, 0., 0., 0);
          if (net->err) {
            return 0.;
          }
        }
        if (MemCheck(net, 0)) {
          return 0.;
        }
        continue;
      }
      sc.mi = S_CAST(int, mi);
      sc.srt = srt;
      sc.nsrt = nsrt;
      sc.cum = cum;
      sc.pmax = pmax;
      sc.pmaxD = pmaxD;
      sc.single = (nsrt == 1);
      sc.v0 = srt[0].v;
      sc.r = r;
      // enumerate daughter columns x (sum c, 0 <= x_i <= r_i) with an odometer; x[nr-1] forced
      int suf[kFisherNetMaxRows + 1];
      suf[nr] = 0;
      for (int i = nr - 1; i >= 0; --i) {
        suf[i] = suf[i + 1] + r[i];
      }
      if (c > suf[0]) {
        continue;
      }
      if (g_gskip) {
        // ---- grouped enumeration: rows j..nr-1 fixed, rows 0..j-1 free with sum mrem ----
        sc.crem = &cols[s + 1];
        sc.mcols = rem_after;
        const long long Nd = S_CAST(long long, suf[0]) - c;
        dd_real KD = DLF(net, Nd);
        for (int l = 0; l < sc.mcols; ++l) {
          KD = ddr_sloppier_sub(KD, DLF(net, sc.crem[l]));
        }
        sc.KD = KD;
        sc.S[0] = 0;
        sc.slf[0] = ddr_make(0., 0.);
        for (int q = 0; q < nr; ++q) {
          sc.S[q + 1] = sc.S[q] + r[q];
          sc.slf[q + 1] = ddr_sloppier_add(sc.slf[q], DLF(net, r[q]));
        }
        sc.cumall = cum[nsrt];
        Rec(&sc, nr, c, 0., ddr_make(0., 0.));
        if (net->err) {
          return 0.;
        }
        continue;
      }
      int* x = sc.x;
      int hi[kFisherNetMaxRows];
      int rem[kFisherNetMaxRows + 1];
      rem[0] = c;
      int i = 0;
      // each level starts at its lowest feasible value
      hi[0] = MINV(r[0], rem[0]);
      x[0] = MAXV(0, rem[0] - suf[1]);
      while (1) {
        if (i < nr - 1) {
          rem[i + 1] = rem[i] - x[i];
          ++i;
          if (i < nr - 1) {
            hi[i] = MINV(r[i], rem[i]);
            x[i] = MAXV(0, rem[i] - suf[i + 1]);
          }
          continue;
        }
        x[nr - 1] = rem[nr - 1];
        Body(&sc);
        if (net->err) {
          return 0.;
        }
        i = nr - 2;
        while ((i >= 0) && (x[i] == hi[i])) {
          --i;
        }
        if (i < 0) {
          break;
        }
        ++x[i];
      }
    }
    if (defer_here) {
      // per node: defer (keep runs) when runs are long on average, i.e. merging
      // would collapse a lot; otherwise merge into an ordinary list as before
      long long nruns = 0;
      for (size_t d = 0; d < net->nxt.n; ++d) {
        RunVec* rv = &net->druns.p[d];
        long long tot = 0;
        for (size_t q = 0; q < rv->n; ++q) {
          tot += rv->p[q].hi - rv->p[q].lo;
        }
        if (!rv->n) {
          continue;
        }
        // merge only if the worst case (every push distinct) still fits;
        // otherwise keep the runs (deferred)
        const uint32_t fits = (MemNow(net) + S_CAST(size_t, tot) * sizeof(Entry) <= net->mem_limit);
        const uint32_t has_list = (net->nxt.p[d].ents.n != 0);
        if ((!has_list) && ((defer_choice && (S_CAST(double, tot) >= g_defer_avg * S_CAST(double, rv->n))) || (!fits))) {
          nruns += S_CAST(long long, rv->n);
        } else {
          if (MergeRuns(net, &net->cur, rv->p, rv->n, &net->nxt.p[d].ents, 1)) {
            return 0.;
          }
          net->mem_new = net->arena_nxt.bytes;
          net->mem_runs -= rv->cap * sizeof(Run);
          RunVecRelease(net, rv);
          if (MemCheck(net, 0)) {
            return 0.;
          }
        }
      }
      if (nruns) {
        // keep the mothers (and their lists) for the deferred last stage
        SwapNodeVecs(&net->prevm, &net->cur);
        const RunVecVec tmp = net->prevruns;
        net->prevruns = net->druns;
        net->druns = tmp;
        deferred = 1;
        ArenaSwap(&net->arena_prev, &net->arena_cur);
      }
      // (otherwise nothing was deferred, and the mothers can go)
    } else if (rem_after > 1) {
      for (size_t d = 0; d < net->nxt.n; ++d) {
        RunVec* rv = &net->druns.p[d];
        if (MergeRuns(net, &net->cur, rv->p, rv->n, &net->nxt.p[d].ents, 1)) {
          return 0.;
        }
        net->mem_new = net->arena_nxt.bytes;
        net->mem_runs -= rv->cap * sizeof(Run);
        RunVecRelease(net, rv);
        if (MemCheck(net, 0)) {
          return 0.;
        }
      }
    }
    SwapNodeVecs(&net->cur, &net->nxt);
    ArenaSwap(&net->arena_cur, &net->arena_nxt);
    ArenaClear(&net->arena_nxt);
    StageRelease(net);  // the old mothers (now in nxt), this stage's runs, map, and buffers
    if (rem_after == 1) {
      break;
    }
  }
  net->out_arena = nullptr;
  const double Ssum = AccGet(&acc);
  const dd_real lnpD = ddr_sloppier_add(ddr_sloppier_add(ddr_make(log(Ssum), 0.), vobsD), cnD);
  return MINV(0., lnpD.x[0] + lnpD.x[1]);
}

static void NetDestroy(Net* net) {
  BigFreeRaw(net->LF, net->lf_ct * sizeof(double));
  BigFreeRaw(net->LFL, net->lf_ct * sizeof(double));
  BigFreeRaw(net->FH, net->lf_ct * sizeof(uint64_t));
  ArenaDestroy(&net->arena_cur);
  ArenaDestroy(&net->arena_nxt);
  ArenaDestroy(&net->arena_prev);
  // run vectors go back to the pool before it is emptied
  RunVecVecRelease(net, &net->druns);
  RunVecVecRelease(net, &net->prevruns);
  RunPoolReset(&net->pool);
  SlotVecFree(&net->tab);
  U32VecFree(&net->seg);
  U32VecFree(&net->seg2);
  U32VecFree(&net->bcnt);
  U32VecFree(&net->bsort_start);
  U64VecFree(&net->cnt_tab);
  EntryVecFree(&net->mbuf);
  NodeVecRelease(&net->cur);
  NodeVecRelease(&net->nxt);
  NodeVecRelease(&net->prevm);
  MapFree(&net->nodemap);
  DoubleVecFree(&net->cum);
  U32VecFree(&net->bidx);
  LeafVecFree(&net->lvbuf);
  DoubleVecFree(&net->wsuf);
  U32VecFree(&net->qb);
  FlushOrdVecFree(&net->ord);
  IntVecFree(&net->mode_i);
  DoubleVecFree(&net->mode_d);
  IntVecFree(&net->cb);
  free(net->ibuf);
}

double FisherNetworkLnP(const int32_t* obs, uint32_t row_ct, uint32_t col_ct, uint32_t midp, uint64_t mem_limit_bytes, FisherCheckFn check_fn, void* check_arg, uint32_t* status_out) {
#ifdef FISHER_DEBUG_OVERRIDES
  ApplyDebugEnv();
#endif
  Net* net = S_CAST(Net*, calloc(1, sizeof(Net)));
  if (!net) {
    *status_out = kFisherNetMemLimit;
    return S_CAST(double, NAN);
  }
  net->midp = (midp != 0);
  net->check_fn = check_fn;
  net->check_arg = check_arg;
  net->mem_limit = mem_limit_bytes? S_CAST(size_t, MINV(mem_limit_bytes, S_CAST(uint64_t, SIZE_MAX))) : SIZE_MAX;
  net->kTieEps = 1e-11;
  net->kTieEpsDD = 1e-25;
  net->NR = 3;
  double lnp = Core(net, obs, S_CAST(int, row_ct), S_CAST(int, col_ct));
  const uint32_t status = net->err;
  if (status) {
    lnp = S_CAST(double, NAN);
  }
  NetDestroy(net);
  free(net);
  *status_out = status;
  return lnp;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
