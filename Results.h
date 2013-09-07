#ifndef RESULTS
#define RESULTS

#include "TMath.h"
#include "openmp.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <malloc.h>

#include "tbb.h"

class Results {

public:
  Results(){}
  ~Results() { ClearAll(); }
  
  // no thread-safe
  void ClearAll();
  inline Bool_t IsEmpty() const { return m_resultsCPU.empty() && m_parallelResultsCPU.empty(); }
  inline UInt_t GetSize() const {
    UInt_t size(0);
    if (m_parallelResultsCPU.empty())
      size = m_resultsCPU.size();
    else {
      for (Int_t i = 0; i<m_parallelResultsCPU.size(); i++)
	size += m_parallelSize[i];
    }
    return size;
  }
  
  // Sequential reduction using Knuth for ALL events (no matter which thread)
  inline TMath::ValueAndError_t NegReduction() const {
    TMath::ValueAndError_t res;
    Double_t const * pResults(0);
    UInt_t iThread(0);
    
    do {

      pResults =   m_parallelResultsCPU.empty() ? &m_resultsCPU.front() : m_parallelResultsCPU[iThread];
      
      Int_t size = m_parallelResultsCPU.empty() ? m_resultsCPU.size() : m_parallelSize[iThread];

      for (Int_t idx = 0; idx<size; idx++) {
	TMath::KnuthAccumulationSub(res.value,res.error,pResults[idx]);
      }
      
      if (m_parallelResultsCPU.empty())
	break;
      
      iThread++;
      
    } while (iThread<m_parallelResultsCPU.size());
    
    return res;
  }
  
  inline Double_t *GetData(UInt_t &nEvents) const {
    UInt_t iEnd;
    return GetData(nEvents,iEnd);
  }
  
  inline Double_t *GetData(UInt_t &nEvents, UInt_t &iEnd, UInt_t iStart = 0, UInt_t nPartialEvents = 0) const {
    // sequential container
    Double_t *ret(0);
    
    if (m_parallelResultsCPU.empty()) {
      nEvents = m_resultsCPU.size();
      ret = &m_resultsCPU[0];
    }
    else {
      // parallel container
      UInt_t rank = OpenMP::GetRankThread();
      if (rank<m_parallelResultsCPU.size()) {
	nEvents = m_parallelSize[rank];
	ret = (Double_t*)__builtin_assume_aligned(m_parallelResultsCPU[rank],32);
	  }
      else {
	// Error
	nEvents = 0;
      }
    }

    iEnd = (nPartialEvents>0) ? TMath::Min(iStart+nPartialEvents,nEvents) : nEvents;
    if (iStart>=iEnd) 
      // Error
      nEvents = 0;

    return ret;
  }

  // no thread-safe
  inline void ResizeParallelContainer(Bool_t forceSequential = false) {
    if (OpenMP::GetMaxNumThreads()>1 && !forceSequential) {
      m_resultsCPU.clear();
      m_parallelResultsCPU.resize(OpenMP::GetMaxNumThreads(),0);
      m_parallelSize.resize(OpenMP::GetMaxNumThreads(),0);
      m_parallelCapacity.resize(OpenMP::GetMaxNumThreads(),0);
    }
    else {
      m_parallelResultsCPU.clear();
      m_parallelSize.clear();
      m_parallelCapacity.clear();
    }
  }

  inline Double_t *AllocateData(UInt_t nEvents) {
    // sequential container
    if (m_parallelResultsCPU.empty()) {
      m_resultsCPU.resize(nEvents);
      return &m_resultsCPU[0];
    }

    // create a new container local to the thread
    if (nEvents > m_parallelCapacity[OpenMP::GetRankThread()]) {
      free(m_parallelResultsCPU[OpenMP::GetRankThread()]);
      m_parallelResultsCPU[OpenMP::GetRankThread()] = (Double_t*)__builtin_assume_aligned(memalign(32,nEvents*sizeof(double)),32);
      m_parallelCapacity[OpenMP::GetRankThread()] = nEvents;
    }
    
    m_parallelSize[OpenMP::GetRankThread()] = nEvents;
    
    
    return (Double_t*)__builtin_assume_aligned(m_parallelResultsCPU[OpenMP::GetRankThread()],32);
    
  }

  private:

  mutable VectorSTD(Double_t) m_resultsCPU; // Global container of results (for not parallel execution)
  mutable std::vector<Double_t*> m_parallelResultsCPU; // results for parallel execution
  mutable std::vector<unsigned int>  m_parallelSize;
  mutable std::vector<unsigned int>  m_parallelCapacity;
};

#endif

