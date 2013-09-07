
#include "PdfPolynomial.h"
#include "openmp.h"

#include <algorithm>
#include <iostream>

Int_t PdfPolynomial::m_blockSize = 1000;

PdfPolynomial::PdfPolynomial(const Char_t* name, const Char_t* title, const Variable &x) :
  AbsPdf(name,title), m_x(&x)
{
}

PdfPolynomial::PdfPolynomial(const Char_t* name, const Char_t* title, const Variable &x,
			     List<Variable> coeff) :
  AbsPdf(name,title), m_x(&x)
{
  m_coeff.AddElement(coeff);
}

PdfPolynomial::~PdfPolynomial()
{
}


Double_t PdfPolynomial::evaluate() const
{
  UInt_t size = m_coeff.GetSize();
  Double_t coeffCPU[size+1];
  loadCoeff(coeffCPU,size);
  return evaluateLocal(m_x->GetVal(),coeffCPU,size);
}

Double_t PdfPolynomial::integral() const
{
  UInt_t order = m_coeff.GetSize();
  Double_t xmaxprod = m_x->GetMax();
  Double_t xminprod = m_x->GetMin();
  Double_t sum = xmaxprod-xminprod;
  for (UInt_t i = 0; i < order; i++) {
      xmaxprod *= m_x->GetMax();
      xminprod *= m_x->GetMin();
      sum += m_coeff.GetElement(i)->GetVal()*(xmaxprod - xminprod)/(i+2);
  }
  return sum;
  
}

Bool_t PdfPolynomial::evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
				   const Double_t invIntegral)
{
  const Double_t *dataCPU = m_data->GetCPUData(*m_x);
  if (dataCPU==0)
    return kFALSE;
  
  Int_t size = m_coeff.GetSize();
  Double_t coeffCPU[size+1];
  loadCoeff(coeffCPU,size);

  // Apply normalization factor
  Int_t idx;
  for (idx=0; idx<=size; idx++)
    coeffCPU[idx] *= invIntegral;

  UInt_t iPartialEnd(0);
  Double_t* resultsCPU = GetDataResultsCPUThread(dataCPU,iPartialEnd,iPartialStart,nPartialEvents);

  if (m_doCalculationBy==kCilk_for) {
    // Cilk doesn't like variable-length array...
    Double_t *vCoeffCPU = coeffCPU;
    CilkSafeCall(
		 // No vectorization for this loop if using Cilk
#ifdef __CILK
#pragma cilk grainsize = Cilk::grainsize
#endif
#pragma ivdep
		 _Cilk_for (Int_t idx = iPartialStart; idx<iPartialEnd; idx++) {
		   resultsCPU[idx] = evaluateLocal(dataCPU[idx],vCoeffCPU,size);
		 }
		 );
  }
  else {

    // Apply Block-Splitting
    Int_t iPartialStartBlock(iPartialStart);
    Int_t iPartialEndBlock(iPartialStart);

    do {

      iPartialEndBlock += m_blockSize;
      if (m_blockSize<=0 || iPartialEndBlock>iPartialEnd) iPartialEndBlock = iPartialEnd;

      std::fill(resultsCPU+iPartialStartBlock,resultsCPU+iPartialEndBlock,coeffCPU[size]);

      for (Int_t order = size-1; order>=0; order--) {
	Double_t vcoeff = coeffCPU[order];

#ifndef USE_CEAN  
#pragma ivdep
	for (idx = iPartialStartBlock; idx<iPartialEndBlock; idx++) {
	  resultsCPU[idx] = evaluateLocalSingleCoeff(dataCPU[idx],vcoeff,resultsCPU[idx]); 
	}
#else
	resultsCPU[iPartialStartBlock:iPartialEndBlock-iPartialStartBlock] = 
	  evaluateLocalSingleCoeff(dataCPU[iPartialStartBlock:iPartialEndBlock-iPartialStartBlock],
				   vcoeff,
				   resultsCPU[iPartialStartBlock:iPartialEndBlock-iPartialStartBlock]); 
#endif
      }

      iPartialStartBlock = iPartialEndBlock;

    } while (iPartialEndBlock<iPartialEnd);

  }

  return kTRUE;
}

