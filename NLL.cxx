
#include "NLL.h"
#include "openmp.h"
#include "Results.h"

#include <iostream>

NLL::NLL(const Char_t* name, const Char_t* title, Data &data, AbsPdf &pdf,
	 bool externalLoop,
	 AbsPdf::DoCalculationBy doCalculationBy) :
  Named(name,title), m_data(&data), m_pdf(&pdf), 
  m_externalLoop(externalLoop),
  m_doCalculationBy(doCalculationBy),m_nBlockEvents(0),
  m_parallelReduction(NLL::kParallel)
  //  m_parallelReduction(NLL::kSequential)
{
  if (m_doCalculationBy!=AbsPdf::kVirtual)
    m_data->DoVectors();

}

NLL::~NLL()
{

}

 // Sequential reduction using Knuth for events of a single thread
void NLL::PartialNegReduction(TMath::IntLog &value,
                                  const Double_t *pResults,
                                  UInt_t iEnd, UInt_t iStart)  {
    value = IntLogAccumulate(value, pResults+iStart, iEnd-iStart);

}  

Double_t NLL::GetVal()
{
  Bool_t doVirtualAlgo(kTRUE);
  m_result = m_zeroValueAndError;
  
  m_pdf->CacheIntegral();
  m_pdf->Init(*m_data,kTRUE,m_doCalculationBy);

  if (m_doCalculationBy!=AbsPdf::kVirtual) {

    const Results *results(0);
    
    if (m_parallelReduction==NLL::kParallel && m_doCalculationBy==AbsPdf::kOpenMP) {
      m_sums.resize(OpenMP::GetMaxNumThreads());
      m_sums.assign(OpenMP::GetMaxNumThreads(),m_zeroValueAndError);
      m_logs.clear();
      m_logs.resize(OpenMP::GetMaxNumThreads());

    }
    
    // Use block splitting only in case of OpenMP
    if (m_doCalculationBy==AbsPdf::kOpenMP && m_nBlockEvents>0) {

      if (m_externalLoop) {
#pragma omp parallel
	{
	  if (OpenMP::GetRankThread()==0) // master
	    results = RunEvaluationBlockSplitting();
	  else
	    RunEvaluationBlockSplitting();
	}	
      }
      else {
	// Internal loop parallelization
	results = RunEvaluationBlockSplitting();
      }
    }
    else { // Without block splitting, TBB or Cilk

      // External loop parallelization (OpenMP case)
      if (m_doCalculationBy==AbsPdf::kOpenMP && m_externalLoop) {
#pragma omp parallel
	{
	  if (OpenMP::GetRankThread()==0) // master
	    results = RunEvaluation();
	  else
	    RunEvaluation();
	}	
      }
      // Case TBB (always with external loop parallelization) or
      // Cilk or
      // OpenMP with internal loop parallelization
      else if (m_doCalculationBy==AbsPdf::kTBB || 
	       m_doCalculationBy==AbsPdf::kCilk_for || 
	       m_doCalculationBy==AbsPdf::kCilk_spawn || 
	       !m_externalLoop) {
	results = RunEvaluation();
      }
    }

    if (results!=0 && MPISvc::SumReduce((int)results->GetSize())==m_data->GetEntries()) {

      // TBB and Cilk_spawn always run the reduction in parallel
      if (m_doCalculationBy!=AbsPdf::kTBB && m_doCalculationBy!=AbsPdf::kCilk_spawn) {

	if (m_parallelReduction==NLL::kParallel) {
	  // Cilk_for already runs the parallel reduction
	  if (m_doCalculationBy!=AbsPdf::kCilk_for) {
	    //final reduction
            for (unsigned int i=0; i!=m_logs.size(); ++i) m_sums[i].value=-0.693147182464599609375*m_logs[i].value();
	    m_result = TMath::DoubleDoubleAccumulation(m_sums);
	  }
	}
	else
	  // Sequential reduction (using Knuth)
	  m_result = results->NegReduction();
	
      }
      
      // MPI Reduction
      // Note that the reduction is locally done in sequential or parallel for OpenMP, 
      // but it is always executed in parallel with MPI
      // Based on Knuth summation
      m_result.value = MPISvc::SumReduce(m_result);

      doVirtualAlgo = kFALSE;
    }
    else {
      m_result.value = 0;
      std::cerr << "Error in the execution of the SIMD Algo for the PDF " << m_pdf->GetName()
		<< ". Do Virtual CPU Algo instead." << std::endl;
    }

  }

  if (doVirtualAlgo) {
    UInt_t nEvents = m_data->GetEntries();
    for (UInt_t i = 0; i<nEvents; ++i) {
      m_data->Get(i);
      m_result.value -= m_pdf->GetLogVal();
    }
  }
  
  if (m_pdf->IsExtended()) {
    m_result.value += m_pdf->ExtendedTerm(m_data->GetEntries());
  }

  return m_result.value;
}

const Results *NLL::RunEvaluationBlockSplitting()
{
  const Results *results(0);

  UInt_t iBlockStart(0);
  
  while (1) {
    results = &(m_pdf->GetValSIMD(iBlockStart,m_nBlockEvents));

    if (RunParallelReduction(results,iBlockStart))
      break;

    iBlockStart += m_nBlockEvents;
  }

  return results;
}

const Results *NLL::RunEvaluation()
{
  const Results *results(0);

  switch (m_doCalculationBy) {
  case AbsPdf::kOpenMP:
    results = &(m_pdf->GetLogValSIMD());
    RunParallelReduction(results);
    break;
  case AbsPdf::kTBB:
    m_result = TBB::reduce(RunTBB(m_pdf),results,m_nBlockEvents);
    // TBB already runs the parallel reduction
    break;
  case AbsPdf::kCilk_for:
      results = &(m_pdf->GetLogValSIMD());
      RunParallelReduction(results);
      break;
  case AbsPdf::kCilk_spawn:
    // Use blocks and external loops
    // Runs always parallel reduction
    if (m_externalLoop) {

      // Do the decomposition
      UInt_t nEvents = MPISvc::GetProcElements(m_data->GetEntries());
      // If the block size is zero, then use CILK_NWORKERS as number of workers
      // otherwise use m_nBlockEvents to determine the number of blocks
      // Minimum number of blocks is always CILK_NWORKERS
      Int_t nBlocks = (m_nBlockEvents==0) ? Cilk::GetNWorkers() : 
	std::max(Cilk::GetNWorkers(),Int_t(double(nEvents)/m_nBlockEvents+0.5));

      int iBlockStart(0), iBlockEnd(0);
      
      CilkSafeCall(
		   cilk::reducer_opadd<TMath::ValueAndError_t> result;
		   for (Int_t block = 0; block<nBlocks-1; block++) {
		     int nBlockEvents = Partitioner::GetElements(nBlocks,block,nEvents,iBlockStart,iBlockEnd);
		     _Cilk_spawn RunCILK(iBlockStart,nBlockEvents,result);

		   }

		   // Take care of the last block
		   results = RunCILK(iBlockEnd,nEvents-iBlockEnd,result);

		   _Cilk_sync;

		   m_result = result.get_value();

		   );
    }
    break;
  }

  return results;
}

Bool_t NLL::RunParallelReduction(const Results *results, UInt_t iStart)
{
  Bool_t doneAllElements(kTRUE);

  if (results->IsEmpty())
    return doneAllElements;

  if (m_parallelReduction==NLL::kSequential || 
      m_doCalculationBy==AbsPdf::kTBB) {
    UInt_t nEvents, iEnd;
    results->GetData(nEvents,iEnd,iStart,m_nBlockEvents);
    // Set the condition when all elements are analyzed
    return (iEnd==nEvents);
  }
  else if (m_doCalculationBy==AbsPdf::kCilk_for) {
    // Reduction in Cilk
    UInt_t nEvents, iEnd;
    const Double_t *pResults = results->GetData(nEvents,iEnd,iStart,m_nBlockEvents);
    m_result = Cilk::reduce<TMath::ValueAndError_t>(pResults,iStart,iEnd);
  }
  else if (m_doCalculationBy==AbsPdf::kOpenMP) {
    // Run using OpenMP
    // Active OpenMP, if not already there
    if (OpenMP::IsInParallel()) {
      doneAllElements = RunParallelReductionOpenMP(results,iStart);
    }
    else {
#pragma omp parallel reduction(&& : doneAllElements)
      {
	doneAllElements = RunParallelReductionOpenMP(results,iStart);
      }
    }
  }

  return doneAllElements;

}

Bool_t NLL::RunParallelReductionOpenMP(const Results *results, UInt_t iStart)
{
  // Parallel reduction (thread by thread) using Knuth
  UInt_t nEvents, iEnd;
  const Double_t *pResults = results->GetData(nEvents,iEnd,iStart,m_nBlockEvents);
  auto localValue = m_logs[OpenMP::GetRankThread()];  
  PartialNegReduction(localValue,pResults,iEnd,iStart);
  m_logs[OpenMP::GetRankThread()] = localValue;

  // Set the condition when all elements are analyzed
  return (iEnd==nEvents);

}
