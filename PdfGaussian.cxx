
#include "PdfGaussian.h"

#include <iostream>

PdfGaussian::PdfGaussian(const Char_t* name, const Char_t* title, Variable &x,
			 Variable &mu, Variable &sigma) :
  AbsPdf(name,title), m_x(&x), m_mu(&mu), m_sigma(&sigma)
{
  
}

Double_t PdfGaussian::evaluate() const
{
  return evaluateLocal(m_x->GetVal(),m_mu->GetVal(),m_sigma->GetVal());
}

Double_t PdfGaussian::integral() const
{
  static const Double_t root2 = TMath::Sqrt2() ;
  static const Double_t rootPiBy2 = TMath::Sqrt(TMath::PiOver2());
  Double_t invxscale = 1./(root2*m_sigma->GetVal());
  Double_t ret = rootPiBy2*m_sigma->GetVal()*
    (TMath::Erf((m_x->GetMax()-m_mu->GetVal())*invxscale)-
     TMath::Erf((m_x->GetMin()-m_mu->GetVal())*invxscale));

  return ret;
  
}

Bool_t PdfGaussian::evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
				 const Double_t invIntegral)
{
  const Double_t *dataCPU = m_data->GetCPUData(*m_x);
  if (dataCPU==0)
    return kFALSE;

  UInt_t iPartialEnd(0);
  Double_t* resultsCPU = GetDataResultsCPUThread(dataCPU,iPartialEnd,iPartialStart,nPartialEvents);

  if (m_doCalculationBy==kCilk_for) {
    CilkSafeCall(
#ifdef __CILK
#pragma cilk grainsize = Cilk::grainsize
#endif
#pragma ivdep
		 _Cilk_for (Int_t idx = (Int_t)iPartialStart; idx<(Int_t)iPartialEnd; idx++) {
		   resultsCPU[idx] = evaluateLocal(dataCPU[idx],m_mu->GetVal(),m_sigma->GetVal())*invIntegral;
		 }
		 );
  }
  else {

#ifndef USE_CEAN
#pragma ivdep
    for (Int_t idx = (Int_t)iPartialStart; idx<(Int_t)iPartialEnd; idx++) {
      resultsCPU[idx] = evaluateLocal(dataCPU[idx],m_mu->GetVal(),m_sigma->GetVal())*invIntegral;
    }
#else
    resultsCPU[iPartialStart:iPartialEnd-iPartialStart] = evaluateLocal(dataCPU[iPartialStart:iPartialEnd-iPartialStart],
									m_mu->GetVal(),m_sigma->GetVal())*invIntegral;
#endif
  }

  return kTRUE;
}

