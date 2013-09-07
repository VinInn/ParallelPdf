
#ifndef NLL_H
#define NLL_H

#include "AbsPdf.h"

#include "Named.h"
#include "TMath.h"

#include "Data.h"

class NLL : public Named {
 public:
  NLL(const Char_t* name, const Char_t* title, Data &data, AbsPdf &pdf, 
      Bool_t externalLoop = kFALSE,
      AbsPdf::DoCalculationBy doCalculationBy = AbsPdf::kOpenMP);

  virtual ~NLL();
  
  Double_t GetVal();
  inline AbsPdf *GetPdf() { return m_pdf; }

  inline void SetBlockEventsSize(UInt_t nBlockEvents) { m_nBlockEvents = nBlockEvents; }

  // Sequential reduction using Knuth for events of a single thread
  static void PartialNegReduction(TMath::IntLog &value,
				  const Double_t *pResults,
				  UInt_t iEnd, UInt_t iStart = 0); // __attribute__((optimize("O2")));

 protected:

  Data *m_data;
  AbsPdf *m_pdf;
  // externalLoop parallelization only in case OpenMP
  // TBB uses always external loop
  Bool_t m_externalLoop;
  AbsPdf::DoCalculationBy m_doCalculationBy;

 private:

  struct RunTBB {
    TMath::ValueAndError_t m_value;
    AbsPdf *m_pdf;
    const Results *m_results;
    RunTBB(AbsPdf *pdf) : m_pdf(pdf), m_results(0) { }
#ifdef USE_TBB
    RunTBB(RunTBB &other, tbb::split) : m_pdf(other.m_pdf), m_results(0) { }
    void operator()(const TBBRange& range) {

      m_results = &(m_pdf->GetLogValSIMD(range.begin(),range.size()));
      
      UInt_t nEvents(0);
      const Double_t *pResults = m_results->GetData(nEvents);
      // Reduction using Knuth
      for (Int_t idx = range.begin(); idx<range.end(); ++idx) {
	TMath::KnuthAccumulationSub(m_value.value,m_value.error,pResults[idx]);
      }

    }
#endif
    void join(RunTBB& other) { 
      m_value += other.m_value;
    }
  };

#ifdef USE_CILK
  // Hook function to run _Cilk_spawn
  // It looks like Cilk doesn't like reference as return values
  const Results *RunCILK(UInt_t iStart, UInt_t nPartialEvents,
			 cilk::reducer_opadd<TMath::ValueAndError_t> &result) { 
    const Results *results = &(m_pdf->GetLogValSIMD(iStart,nPartialEvents));
    UInt_t nEvents, iEnd;
    const Double_t *pResults = results->GetData(nEvents,iEnd,iStart,nPartialEvents);
    // Do the reduction
    for (Int_t idx = (Int_t)iStart; idx<(Int_t)iEnd; idx++) {
      result -= TMath::ValueAndError_t(pResults[idx]);
    }
    return results;
  }
#endif


  TMath::ValueAndError_t m_result;
  std::vector<TMath::ValueAndError_t> m_sums; // result and error of the partial parallel sums (logs of below)
  std::vector<TMath::IntLog> m_logs; // result the partial parallel sums (still as exp,mant)

  UInt_t m_nBlockEvents;
  TMath::ValueAndError_t m_zeroValueAndError;

  enum ParallelReduction { kSequential=0, 
			   kParallel };
  // TBB runs alwyas the parallel reduction
  ParallelReduction m_parallelReduction;

  const Results *RunEvaluationBlockSplitting();
  const Results *RunEvaluation();
  Bool_t RunParallelReduction(const Results *results, UInt_t iStart = 0);
  Bool_t RunParallelReductionOpenMP(const Results *results, UInt_t iStart);
 
};

#endif
