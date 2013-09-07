
#include "PdfBreitWigner.h"

#include <iostream>

PdfBreitWigner::PdfBreitWigner(const Char_t* name, const Char_t* title, Variable &x,
			       Variable &mu, Variable &width) :
  AbsPdf(name,title), m_x(&x), m_mu(&mu), m_width(&width)
{
  
}

Double_t PdfBreitWigner::evaluate() const
{
  return evaluateLocal(m_x->GetVal(),m_mu->GetVal(),m_width->GetVal());
}

Double_t PdfBreitWigner::integral() const
{
  Double_t c = 2./m_width->GetVal();
  Double_t ret = c*(TMath::ATan(c*(m_x->GetMax()-m_mu->GetVal())) - TMath::ATan(c*(m_x->GetMin()-m_mu->GetVal())));
  
  return ret;
  
}

Bool_t PdfBreitWigner::evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
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
		   resultsCPU[idx] = evaluateLocal(dataCPU[idx],m_mu->GetVal(),m_width->GetVal())*invIntegral;
		 }
		 );
  }
  else {

#ifndef USE_CEAN
#pragma ivdep
    for (Int_t idx = (Int_t)iPartialStart; idx<(Int_t)iPartialEnd; idx++) {
      resultsCPU[idx] = evaluateLocal(dataCPU[idx],m_mu->GetVal(),m_width->GetVal())*invIntegral;
    }
#else
    resultsCPU[iPartialStart:iPartialEnd-iPartialStart] = evaluateLocal(dataCPU[iPartialStart:iPartialEnd-iPartialStart],
									m_mu->GetVal(),m_width->GetVal())*invIntegral;
#endif
  }

  return kTRUE;
}

