
#ifndef ABS_PDF
#define ABS_PDF

#include "mpisvc.h"

#include <vector>
#include <iostream>
#include <cassert>

#include "Named.h"
#include "TMath.h"
#include "Data.h"
#include "Results.h"
#include "openmp.h"
#include "cilk.h"

#define RooAbsPdf AbsPdf

class AbsPdf : public Named {
 public:

  enum DoCalculationBy { kVirtual=1, kOpenMP, kTBB, kCilk_for, kCilk_spawn, kMIC };

  AbsPdf(const Char_t* name, const Char_t* title);
  virtual ~AbsPdf();

  virtual void RandomizeFloatParameters();
  virtual void GetParameters(List<Variable>& parameters) { }

  // no thread-safe
  void Init(const Data& data, Bool_t doLog, DoCalculationBy doCalculationBy = kOpenMP) {
    Init(data,&doLog,doCalculationBy);
  }

  // no thread-safe
  virtual void Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy = kOpenMP);

  virtual Double_t GetLogVal() { return TMath::Log(GetVal()); }
  virtual const Results &GetLogValSIMD(UInt_t iStart = 0, UInt_t nPartialEvents = 0);

  virtual Double_t GetVal() { return evaluate()/GetIntegral(); }
  virtual const Results &GetValSIMD(UInt_t iStart = 0, UInt_t nPartialEvents = 0);
  
  // no thread-safe
  virtual Double_t GetIntegral() {
    if (m_cacheIntegral) {
      m_integral = integral();
      m_cacheIntegral = kFALSE;
    }

    return m_integral;
  }

  // no thread-safe
  virtual void CacheIntegral() { m_cacheIntegral = kTRUE; GetIntegral(); }

  virtual Double_t ExtendedTerm(UInt_t observed) const { return .0; }
  virtual Bool_t IsExtended() const { return kFALSE; }
  virtual Double_t ExpectedEvents() const { return .0; }

  // no thread-safe
  virtual void ClearResults(Bool_t recursive = kFALSE);

  inline Int_t GetResultsSize() const { return m_resultsCPU.GetSize(); }

 protected:

  virtual Bool_t evaluateSIMD(const UInt_t& /*iPartialStart*/, const UInt_t& /*nPartialEvents*/, 
			      const Double_t /*invIntegral = 1.*/) {
    std::cerr << "No SIMD Algo for pdf " << GetName() << "!!!" << std::endl;
    return kFALSE;
  }
  virtual Double_t evaluate() const = 0;
  virtual Double_t integral() const = 0;

  inline Double_t *GetDataResultsCPUThread(const Double_t *&dataCPU, UInt_t &iPartialEnd, 
					   const UInt_t iPartialStart, const UInt_t nPartialEvents) { 

    Int_t iStartMPI, iStartOpenMP(0), iEnd(0);
    UInt_t nEvents = m_data->GetEntries(); // Full sample

    // First decomposition in MPI
    nEvents = MPISvc::GetProcElements(nEvents,iStartMPI,iEnd);
    Double_t *resultsCPU(0);

    if (m_doCalculationBy==kOpenMP) {
      // Second decomposition in OpenMP
      nEvents = OpenMP::GetThreadElements(nEvents,iStartOpenMP,iEnd);

      resultsCPU = m_resultsCPU.AllocateData(nEvents);

    } 
    else { // case TBB or Cilk
      UInt_t nResults;
      resultsCPU = m_resultsCPU.GetData(nResults);
      assert(nEvents==nResults);
    }

    iPartialEnd = (nPartialEvents>0) ? TMath::Min(iPartialStart+nPartialEvents,nEvents) : nEvents;

    dataCPU += (iStartMPI+iStartOpenMP);

    return resultsCPU;
  }

  inline Int_t GetParallelNumElements(UInt_t nEvents) const {
    // First decomposition in MPI
    nEvents = MPISvc::GetProcElements(nEvents);

    if (m_doCalculationBy==kOpenMP) {
      // Second decomposition in OpenMP
      nEvents = OpenMP::GetThreadElements(nEvents);
    }

    return nEvents;
  }

  Results m_resultsCPU;
  Results m_logResultsCPU;

  // By default the calculation is done with Virtual algorithm
  DoCalculationBy m_doCalculationBy;
  const Data* m_data;

  Double_t m_integral;
  Bool_t m_cacheIntegral;

 private:

  Int_t RunGetLogVal(UInt_t iStart, UInt_t nPartialEvents);

};

#endif
