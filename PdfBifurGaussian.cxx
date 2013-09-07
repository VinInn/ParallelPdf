
#include "PdfBifurGaussian.h"

#include <iostream>

PdfBifurGaussian::PdfBifurGaussian(const Char_t* name, const Char_t* title, Variable &x,
				   Variable &mu, Variable &sigmaL, Variable &sigmaR) :
  AbsPdf(name,title), m_x(&x), m_mu(&mu), m_sigmaL(&sigmaL), m_sigmaR(&sigmaR)
{

}

Double_t PdfBifurGaussian::evaluate() const
{
  return evaluateLocal(m_x->GetVal(),m_mu->GetVal(),m_sigmaL->GetVal(),m_sigmaR->GetVal());
}

Double_t PdfBifurGaussian::integral() const
{
  static const Double_t root2 = TMath::Sqrt2() ;
  static const Double_t rootPiBy2 = TMath::Sqrt(TMath::PiOver2());
  Double_t invxscaleL = 1./(root2*m_sigmaL->GetVal());
  Double_t invxscaleR = 1./(root2*m_sigmaR->GetVal());

  Double_t integral = 0.0;
  if(m_x->GetMax() < m_mu->GetVal())	{
    integral = m_sigmaL->GetVal()*(TMath::Erf((m_x->GetMax()-m_mu->GetVal())*invxscaleL)-TMath::Erf((m_x->GetMin()-m_mu->GetVal())*invxscaleL));
  }
  else if (m_x->GetMin() > m_mu->GetVal()) {
    integral = m_sigmaR->GetVal()*(TMath::Erf((m_x->GetMax()-m_mu->GetVal())*invxscaleR)-TMath::Erf((m_x->GetMin()-m_mu->GetVal())*invxscaleR));
  }
  else {
    integral = m_sigmaR->GetVal()*TMath::Erf((m_x->GetMax()-m_mu->GetVal())*invxscaleR)-m_sigmaL->GetVal()*TMath::Erf((m_x->GetMin()-m_mu->GetVal())*invxscaleL);
  }

  return integral*rootPiBy2;

}

Bool_t PdfBifurGaussian::evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
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
		   resultsCPU[idx] = evaluateLocal(dataCPU[idx],m_mu->GetVal(),m_sigmaL->GetVal(),m_sigmaR->GetVal())*invIntegral;
		 }
		 );
  }
  else {

#ifndef USE_CEAN
#pragma ivdep
    for (Int_t idx = (Int_t)iPartialStart; idx<(Int_t)iPartialEnd; idx++) {
      resultsCPU[idx] = evaluateLocal(dataCPU[idx],m_mu->GetVal(),m_sigmaL->GetVal(),m_sigmaR->GetVal())*invIntegral;
    }
#else
    resultsCPU[iPartialStart:iPartialEnd-iPartialStart] = evaluateLocal(dataCPU[iPartialStart:iPartialEnd-iPartialStart],
									m_mu->GetVal(),m_sigmaL->GetVal(),m_sigmaR->GetVal())*invIntegral;
#endif
  }

  return kTRUE;
}

