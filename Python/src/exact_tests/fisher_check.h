// Cancellation support for fisher_exact()'s 2xk/3xk searches: a wall-clock
// time limit, plus Python signal handling (so e.g. Ctrl-C raises
// KeyboardInterrupt instead of being ignored until the search finishes).
//
// The search runs without the GIL, and calls FisherPyCheck() every few
// milliseconds (see FisherCheckFn in fisher.h).

#ifndef FISHER_CHECK_H
#define FISHER_CHECK_H

#include <Python.h>
#include <stdint.h>

#include <chrono>

typedef std::chrono::steady_clock FisherPyClock;

typedef struct FisherPyCheckStateStruct {
  FisherPyClock::time_point deadline;
  FisherPyClock::time_point next_signal_check;
  int has_deadline;
  int timed_out;
  // Set when a signal handler raised an exception.  The exception is left
  // set in the calling thread's error indicator, to be propagated once the
  // search has unwound.
  int signalled;
} FisherPyCheckState;

// time_limit is in seconds; non-finite or very large values mean no limit.
static inline void FisherPyCheckInit(double time_limit, FisherPyCheckState* statep) {
  const FisherPyClock::time_point now = FisherPyClock::now();
  statep->has_deadline = (time_limit < 1e9);
  if (statep->has_deadline) {
    statep->deadline = now + std::chrono::duration_cast<FisherPyClock::duration>(std::chrono::duration<double>(time_limit));
  }
  statep->next_signal_check = now + std::chrono::milliseconds(100);
  statep->timed_out = 0;
  statep->signalled = 0;
}

// Must be called without the GIL held.
static inline uint32_t FisherPyCheck(void* arg) {
  FisherPyCheckState* statep = static_cast<FisherPyCheckState*>(arg);
  const FisherPyClock::time_point now = FisherPyClock::now();
  if (statep->has_deadline && (now >= statep->deadline)) {
    statep->timed_out = 1;
    return 1;
  }
  // Acquiring the GIL can block while other threads run Python code, so
  // signal checks are rate-limited.
  if (now >= statep->next_signal_check) {
    statep->next_signal_check = now + std::chrono::milliseconds(100);
    const PyGILState_STATE gil_state = PyGILState_Ensure();
    const int err = PyErr_CheckSignals();
    PyGILState_Release(gil_state);
    if (err) {
      statep->signalled = 1;
      return 1;
    }
  }
  return 0;
}

#endif  // FISHER_CHECK_H
