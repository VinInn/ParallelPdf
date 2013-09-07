
#include "PdfAdd.h"
#include "openmp.h"

#include <iostream>

PdfAdd::PdfAdd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2, Variable &fraction) :
  AbsPdf(name,title), m_isExtended(kFALSE)
{
  m_pdfs.AddElement(pdf1);
  m_pdfs.AddElement(pdf2);
  m_fractions.AddElement(fraction);

}

PdfAdd::PdfAdd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs, List<Variable> fractions) :
  AbsPdf(name,title), m_isExtended(kFALSE)
{
  if (pdfs.GetSize()!=fractions.GetSize() && pdfs.GetSize()!=fractions.GetSize()-1) {
    std::cerr << GetName() << ":: Wrong number of fractions!" << std::endl;
    assert(0);
  }
  
  if (pdfs.GetSize()==fractions.GetSize())
    m_isExtended = kTRUE;
  
  m_pdfs.AddElement(pdfs);
  m_fractions.AddElement(fractions);
  
}

void PdfAdd::GetParameters(List<Variable>& parameters)
{
  parameters.AddElement(m_fractions);
  AbsPdf *pdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  while ((pdf = iter_pdfs.Next())!=0) {
    pdf->GetParameters(parameters);
  }
}

Double_t PdfAdd::ExtendedTerm(UInt_t observed) const
{
  Double_t expected = ExpectedEvents();
  return expected-observed*TMath::Log(expected);

}

Double_t PdfAdd::ExpectedEvents() const
{
  Double_t nEvents(0);
  if (m_isExtended) {
    Variable *var(0);
    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    while ((var = iter_fractions.Next())!=0)
      nEvents += var->GetVal();
  }

  return nEvents;
}

void PdfAdd::ClearResults(Bool_t recursive)
{
  AbsPdf::ClearResults();
  if (recursive) {
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->ClearResults(recursive);
    }
  }
}

void PdfAdd::Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy)
{
  AbsPdf::Init(data,doLog,doCalculationBy);
  AbsPdf *pdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  while ((pdf = iter_pdfs.Next())!=0) {
    pdf->Init(data,doLog,doCalculationBy);
  }

}

void PdfAdd::CacheIntegral()
{
  AbsPdf::CacheIntegral();
  AbsPdf *pdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  while ((pdf = iter_pdfs.Next())!=0) {
    pdf->CacheIntegral();
  }

}

Double_t PdfAdd::evaluate() const
{
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
  
  AbsPdf *pdf = iter_pdfs.Next();
  Variable *var(0);
  Double_t lastFraction = 1.;
  Double_t ret(0);
  
  while ((var = iter_fractions.Next())!=0) {
    lastFraction -= var->GetVal();
    ret += var->GetVal()*pdf->GetVal();
    pdf = iter_pdfs.Next();
  }
  
  if (!m_isExtended)
    ret += lastFraction*pdf->GetVal();
  
  return ret;
}

const Results &PdfAdd::GetValSIMD(UInt_t iStart, UInt_t nPartialEvents)
{
  assert(m_doCalculationBy!=kVirtual && 
	 evaluateSIMD(iStart,nPartialEvents,1./GetIntegral()));

  return m_resultsCPU;
}


Bool_t PdfAdd::evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents, 
			    const Double_t invIntegral)
{
  List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
  Variable *var = iter_fractions.Next();
  Double_t coeff = var->GetVal();
  Double_t lastFraction = 1.;

  // Total number of events
  Int_t nEvents = (Int_t)m_data->GetEntries();
  Bool_t isOK(kTRUE);
  UInt_t iPdf(0);
  AbsPdf *pdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
 
  while ((pdf = iter_pdfs.Next()) && isOK) {
    const Results* pResultsPdf = &(pdf->GetValSIMD(iPartialStart,nPartialEvents));
    if (pResultsPdf->IsEmpty()) {
      isOK = kFALSE;
      break;
    }

    if (iPdf>0) {
      lastFraction -= coeff;
      if ((var = iter_fractions.Next())!=0)
	coeff = var->GetVal();
      else
	coeff = lastFraction;
    }

    // Take in account the normalization
    coeff *= invIntegral;

    if (m_doCalculationBy==AbsPdf::kTBB || m_doCalculationBy==AbsPdf::kCilk_for || 
	m_doCalculationBy==AbsPdf::kCilk_spawn || OpenMP::IsInParallel()) {
      isOK = RunEvaluate(nEvents,iPartialStart,nPartialEvents,
			 pResultsPdf,iPdf,coeff);
    }
    else {
#pragma omp parallel reduction(&& : isOK)
      {
	isOK = RunEvaluate(nEvents,iPartialStart,nPartialEvents,
			   pResultsPdf,iPdf,coeff);
      }
    }

    iPdf++;

  }

  return isOK;
  
}


Bool_t PdfAdd::RunEvaluate(Int_t nEvents, const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			   const Results* pResultsPdf, UInt_t iPdf, Double_t coeff)
{
  Bool_t isOK;
  Int_t idx, nEventsThread = GetParallelNumElements(nEvents);

  UInt_t nResults(0), iPartialEnd;
  const Double_t* resultsPdf = pResultsPdf->GetData(nResults,iPartialEnd,
						    iPartialStart,nPartialEvents);

  if ((isOK=(nEventsThread==nResults))) {
    
    Double_t* resultsCPU(0);
    if (m_doCalculationBy==AbsPdf::kOpenMP)
      resultsCPU = m_resultsCPU.AllocateData(nEventsThread);
    else // case TBB or Cilk
      resultsCPU = m_resultsCPU.GetData(nResults);

    const Int_t iE = iPartialEnd;

    if (0==iPdf) {

      if (m_doCalculationBy==kCilk_for) {
	CilkSafeCall(
#ifdef __CILK
#pragma cilk grainsize = Cilk::grainsize
#endif
#pragma ivdep
		     _Cilk_for (Int_t idx = (Int_t)iPartialStart; idx<iE; idx++) {
		       resultsCPU[idx] = coeff*resultsPdf[idx];
		     }
		     );
      }
      else {

#ifndef USE_CEAN
#pragma ivdep
	for (idx = (Int_t)iPartialStart; idx<iE; idx++) {
	  resultsCPU[idx] = coeff*resultsPdf[idx];
	}
#else
	resultsCPU[iPartialStart:iPartialEnd-iPartialStart] = coeff*resultsPdf[iPartialStart:iPartialEnd-iPartialStart];
#endif

      }
    }
    else {
     
      if (m_doCalculationBy==kCilk_for) {
	CilkSafeCall(
#ifdef __CILK
#pragma cilk grainsize = Cilk::grainsize
#endif
#pragma ivdep
		     _Cilk_for (Int_t idx = (Int_t)iPartialStart; idx<iE; idx++) {
		       resultsCPU[idx] += coeff*resultsPdf[idx];
		     }
		     );
      }
      else {

#ifndef USE_CEAN
#pragma ivdep
	for (idx = (Int_t)iPartialStart; idx<iE; idx++) {
	  resultsCPU[idx] += coeff*resultsPdf[idx];
	}
#else
	resultsCPU[iPartialStart:iPartialEnd-iPartialStart] += coeff*resultsPdf[iPartialStart:iPartialEnd-iPartialStart];
#endif
      }

    }
    
  }
  
  return isOK;

}
