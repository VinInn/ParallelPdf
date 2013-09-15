
#ifndef OPENMP_H
#define OPENMP_H

#include "Partitioner.h"
#include "omp.h"
#include <algorithm>

namespace OpenMP {

  inline double GetTime() {
    return omp_get_wtime();
  }

  inline int GetNumThreads() {
    return std::max(1,omp_get_num_threads());
  }
  
  inline int GetRankThread() {
    return omp_get_thread_num();
  }

  // determine the number of elements per each thread (statically allocated)
  // determine also the first element per each thread
  inline int GetThreadElements(int nEvents, int& iStart, int& iEnd) {
    return Partitioner::GetElements(GetNumThreads(),GetRankThread(),nEvents,iStart,iEnd);
  }

  inline int GetThreadElements(int nEvents) {
    int iStart(0), iEnd(0);
    return GetThreadElements(nEvents,iStart,iEnd);
  }

  inline int GetMaxNumThreads() {
    return omp_get_max_threads();
  }

  inline bool IsInParallel() {
    return (omp_in_parallel() || GetMaxNumThreads()==1);
  }

}


#endif
