
#ifndef TBB_H
#define TBB_H

//#define TBB_ALLOCATOR

#include <iostream>
#include <cstdlib>
#include <vector>

#include "TMath.h"

class Results;


#ifdef TBB_ALLOCATOR

#ifdef USE_TBB

// TBB Allocator
#include "tbb/scalable_allocator.h"
#define VectorSTD(t) std::vector<t,tbb::scalable_allocator<t> >

#endif

#else

// STD Allocator
#define VectorSTD(t) std::vector<t>

#endif


#ifdef USE_TBB

#include <tbb/task_scheduler_init.h>

#define TBB_PREVIEW_DETERMINISTIC_REDUCE 1
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"

#define TBBRange            tbb::blocked_range<size_t>

#define TBBSafeCall(p)      p

#else

#define TBBSafeCall(p)      std::cerr << "Not compiled with TBB!!!" << std::endl; exit(1);

#endif


namespace TBB {

  template<class T> TMath::ValueAndError_t reduce(T functor, const Results *&results, size_t grainsize) {
    TMath::ValueAndError_t res;
    results = 0;
    TBBSafeCall(
		int size = functor.m_pdf->GetResultsSize();
		tbb::parallel_deterministic_reduce(TBBRange(0,size,grainsize),functor);
		res = functor.m_value;
		results = functor.m_results;
		);
    return res;
  }

}

#endif
