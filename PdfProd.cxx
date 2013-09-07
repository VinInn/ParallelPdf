
#include "PdfProd.h"
#include "openmp.h"

#include <iostream>

PdfProd::PdfProd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2) :
  AbsPdf(name,title)
{
  m_pdfs.AddElement(pdf1);
  m_pdfs.AddElement(pdf2);

}

PdfProd::PdfProd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs) :
  AbsPdf(name,title)
{
  m_pdfs.AddElement(pdfs);

}

void PdfProd::GetParameters(List<Variable>& parameters)
{
  AbsPdf *pdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  while ((pdf = iter_pdfs.Next())!=0) {
    pdf->GetParameters(parameters);
  }
}

void PdfProd::ClearResults(Bool_t recursive)
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


void PdfProd::CacheIntegral()
{
  AbsPdf::CacheIntegral();
  AbsPdf *pdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  while ((pdf = iter_pdfs.Next())!=0) {
    pdf->CacheIntegral();
  }

}

void PdfProd::Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy)
{
  AbsPdf::Init(data,doLog,doCalculationBy);
  AbsPdf *pdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  while ((pdf = iter_pdfs.Next())!=0) {
    pdf->Init(data,doLog,doCalculationBy);
  }

}


Double_t PdfProd::evaluate() const
{
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  AbsPdf *pdf = iter_pdfs.Next();
  Double_t ret = pdf->GetVal();
  
  while ((pdf = iter_pdfs.Next())!=0) {
    ret *= pdf->GetVal();
  }

  return ret;
}

const Results &PdfProd::GetValSIMD(UInt_t iStart, UInt_t nPartialEvents)
{
  // No normalization is required
  assert(m_doCalculationBy!=kVirtual && 
         evaluateSIMD(iStart,nPartialEvents,1.));

  return m_resultsCPU;

}


Bool_t PdfProd::evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			     const Double_t /*invIntegral*/)
{
  Int_t nEvents = (Int_t)m_data->GetEntries();
  Bool_t isOK(kTRUE);
  UInt_t iPdf(0);
  List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
  AbsPdf *pdf(0);

  while ((pdf = iter_pdfs.Next()) && isOK) {
    const Results* pResultsPdf = &(pdf->GetValSIMD(iPartialStart,nPartialEvents));
    if (pResultsPdf->IsEmpty()) {
      isOK = kFALSE;
      break;
    }

    if (m_doCalculationBy==AbsPdf::kTBB || m_doCalculationBy==AbsPdf::kCilk_for || 
	m_doCalculationBy==AbsPdf::kCilk_spawn || OpenMP::IsInParallel()) {
      isOK = RunEvaluate(nEvents,iPartialStart,nPartialEvents,
			 pResultsPdf,iPdf);
    }
    else {
#pragma omp parallel reduction(&& : isOK)
      {
	isOK = RunEvaluate(nEvents,iPartialStart,nPartialEvents,
			   pResultsPdf,iPdf);
      }
    }

    iPdf++;
  }

  return isOK;

}


Bool_t PdfProd::RunEvaluate(Int_t nEvents, const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			    const Results* pResultsPdf, UInt_t iPdf)
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
		       resultsCPU[idx] = resultsPdf[idx];
		     }
		     );
      }
      else {

#ifndef USE_CEAN
#pragma ivdep
	for (idx = (Int_t)iPartialStart; idx<iE; idx++) {
	  resultsCPU[idx] = resultsPdf[idx];
	}
#else
	resultsCPU[iPartialStart:iPartialEnd-iPartialStart] = resultsPdf[iPartialStart:iPartialEnd-iPartialStart];
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
		       resultsCPU[idx] *= resultsPdf[idx];
		     }
		     );
      }
      else {

#ifndef USE_CEAN
#pragma ivdep
	for (idx = (Int_t)iPartialStart; idx<iE; idx++) {
	  resultsCPU[idx] *= resultsPdf[idx];
	}
#else
	resultsCPU[iPartialStart:iPartialEnd-iPartialStart] *= resultsPdf[iPartialStart:iPartialEnd-iPartialStart];
#endif

      }
    }
  }

  return isOK;
}
