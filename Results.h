
#ifndef RESULTS
#define RESULTS

#include "TMath.h"
#include "openmp.h"

#include <iostream>
#include <vector>

#include "tbb.h"

class Results {

 public:
  Results();
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
	size += m_parallelResultsCPU[i]->size();
    }
    return size;
  }

  // Sequential reduction using Knuth for ALL events (no matter which thread)
  inline TMath::ValueAndError_t NegReduction() const {
    TMath::ValueAndError_t res;
    VectorSTD(Double_t) *results(0);
    Double_t *pResults(0);
    UInt_t iThread(0);

    do {

      results = m_parallelResultsCPU.empty() ? &(m_resultsCPU) : m_parallelResultsCPU[iThread];
      pResults = &((*results)[0]);

      Int_t size = results->size();
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
	nEvents = m_parallelResultsCPU[rank]->size();
	ret = &(*m_parallelResultsCPU[rank])[0];
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
    }
    else
      m_parallelResultsCPU.clear();
  }

  inline Double_t *AllocateData(UInt_t nEvents) {
    // sequential container
    if (m_parallelResultsCPU.empty()) {
      m_resultsCPU.resize(nEvents);
      return &m_resultsCPU[0];
    }

    // create a new container local to the thread
    if (0==m_parallelResultsCPU[OpenMP::GetRankThread()])
      m_parallelResultsCPU[OpenMP::GetRankThread()] = new VectorSTD(Double_t);

    m_parallelResultsCPU[OpenMP::GetRankThread()]->resize(nEvents);

    return &((*m_parallelResultsCPU[OpenMP::GetRankThread()])[0]);

  }

 private:

  mutable VectorSTD(Double_t) m_resultsCPU; // Global container of results (for not parallel execution)
  mutable std::vector<VectorSTD(Double_t)*> m_parallelResultsCPU; // results for parallel execution

};

#endif
