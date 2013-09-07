
#include "AbsPdf.h"
#include "TRandom.h"

AbsPdf::AbsPdf(const Char_t* name, const Char_t* title) :
  Named(name,title), m_doCalculationBy(kVirtual),
  m_data(0),
  m_integral(0), m_cacheIntegral(kTRUE)
{
}

AbsPdf::~AbsPdf()
{
  AbsPdf::ClearResults();
}

void AbsPdf::ClearResults(Bool_t /*recursive*/)
{
  // Delete the data containers
  m_resultsCPU.ClearAll();
  m_logResultsCPU.ClearAll();

  // Reset the evaluation to the default method, i.e. kVirtual
  m_doCalculationBy = kVirtual;
}

void AbsPdf::RandomizeFloatParameters()
{
  List<Variable> pdfPars;
  GetParameters(pdfPars);
  pdfPars.Sort();
  pdfPars.ResetIterator();
  Variable *par(0);
  TRandom rand;
  while ((par = pdfPars.Next())!=0) {
    if (!par->IsConstant()) {
      std::cout << par->GetName() << " = " << par->GetVal();
      par->SetVal(rand.Uniform(par->getMin(),par->getMax()));
      std::cout << " --> " << par->GetVal() << std::endl;
    }
  }

}

void AbsPdf::Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy)
{
  m_data = &data;  
  m_doCalculationBy = doCalculationBy;

  switch (m_doCalculationBy) {
  case kOpenMP:
    // Use parallel container for OpenMP
    m_resultsCPU.ResizeParallelContainer();
    if (*doLog)
      m_logResultsCPU.ResizeParallelContainer();
    break;
  case kTBB:
  case kCilk_for:
  case kCilk_spawn:
    // Use a sequential container for TBB and Cilk
    m_resultsCPU.ResizeParallelContainer(kTRUE);

    // Set the dimension of the container
    m_resultsCPU.AllocateData(GetParallelNumElements(m_data->GetEntries()));

    if (*doLog) {
      // Use a sequential container for TBB and Cilk 
      m_logResultsCPU.ResizeParallelContainer(kTRUE);
      // Set the dimension of the container
      m_logResultsCPU.AllocateData(GetParallelNumElements(m_data->GetEntries()));
    }
    break;
  default:
    std::cerr << "Init not valid for PDF " << GetName() << "!!!" << std::endl; 
    exit(1);
  }

  // Only the first call will have the log results
  if (*doLog)
    *doLog = kFALSE;

}

const Results &AbsPdf::GetLogValSIMD(UInt_t iStart, UInt_t nPartialEvents)
{
  assert(m_doCalculationBy!=kVirtual && !GetValSIMD(iStart,nPartialEvents).IsEmpty());

  // Does the Log of the results
  if (m_doCalculationBy==kOpenMP) {
    // Run using OpenMP
    // Active OpenMP, if not already there
    if (OpenMP::IsInParallel()) {
      RunGetLogVal(iStart,nPartialEvents);
    }
    else {
#pragma omp parallel
      {
	RunGetLogVal(iStart,nPartialEvents);
      }
    }
  }
  else // case TBB or Cilk
    RunGetLogVal(iStart,nPartialEvents);

  return m_logResultsCPU; // GetLogVal returns always the CPU results

}

Int_t AbsPdf::RunGetLogVal(UInt_t iStart, UInt_t nPartialEvents)
{
  UInt_t nEventsThread, iEnd;
  Double_t* pResults = m_resultsCPU.GetData(nEventsThread,iEnd,iStart,nPartialEvents);
  Double_t* resultsCPU(0);
  if (m_doCalculationBy==kOpenMP)
    resultsCPU = m_logResultsCPU.AllocateData(nEventsThread);
  else // case TBB or Cilk
    resultsCPU = m_logResultsCPU.GetData(nEventsThread);

  if (m_doCalculationBy==kCilk_for) {
    CilkSafeCall(
#ifdef __CILK
#pragma cilk grainsize = Cilk::grainsize
#endif
#pragma ivdep
		 _Cilk_for (Int_t idx = (Int_t)iStart; idx<(Int_t)iEnd; idx++) {
		   resultsCPU[idx] = TMath::Log(pResults[idx]);
		 }
		 );
  }
  else {

#ifndef USE_CEAN
#pragma ivdep
    for (Int_t idx = (Int_t)iStart; idx<(Int_t)iEnd; idx++) {
      resultsCPU[idx] = TMath::Log(pResults[idx]);
    }
#else
    resultsCPU[iStart:iEnd-iStart] = TMath::Log(pResults[iStart:iEnd-iStart]);
#endif

  }

  return nEventsThread;

}

const Results &AbsPdf::GetValSIMD(UInt_t iStart, UInt_t nPartialEvents)
{
  assert(m_doCalculationBy!=kVirtual);

  Bool_t isOK(kFALSE);

  if (m_doCalculationBy==kOpenMP) {
    // Run using OpenMP
    // Active OpenMP, if not already there
    // Note: TBB always with external loop
    if (OpenMP::IsInParallel()) {
      isOK = evaluateSIMD(iStart,nPartialEvents,1./GetIntegral());
    }
    else {
      isOK = kTRUE;
#pragma omp parallel reduction(&& : isOK)
      {
	isOK = evaluateSIMD(iStart,nPartialEvents,1./GetIntegral());
      }
    }
  } 
  else // case TBB or Cilk
    isOK = evaluateSIMD(iStart,nPartialEvents,1./GetIntegral());

  assert(isOK);

  return m_resultsCPU;

}

