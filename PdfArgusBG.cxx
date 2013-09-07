
#include "PdfArgusBG.h"
#include "openmp.h"

#include <iostream>

PdfArgusBG::PdfArgusBG(const Char_t* name, const Char_t* title, Variable &m,
		       Variable &m0, Variable &c) :
  AbsPdf(name,title), m_m(&m), m_m0(&m0), m_c(&c)
{

}

Double_t PdfArgusBG::evaluate() const
{
  return evaluateLocal(m_m->GetVal(),1./m_m0->GetVal(),m_c->GetVal());
}

Double_t PdfArgusBG::integral() const
{
  const Double_t rootPi = TMath::Sqrt(TMath::Pi());
  Double_t min = (m_m->GetMin() < m_m0->GetVal()) ? m_m->GetMin() : m_m0->GetVal();
  Double_t max = (m_m->GetMax() < m_m0->GetVal()) ? m_m->GetMax() : m_m0->GetVal();
  Double_t f1 = (1.-TMath::Power(min/m_m0->GetVal(),2));
  Double_t f2 = (1.-TMath::Power(max/m_m0->GetVal(),2));
  Double_t aLow, aHigh ;
  aLow  = -0.5*m_m0->GetVal()*m_m0->GetVal()*(TMath::Exp(m_c->GetVal()*f1)*TMath::Sqrt(f1)/m_c->GetVal() + 
					      0.5/TMath::Power(-m_c->GetVal(),1.5)*rootPi*TMath::Erf(TMath::Sqrt(-m_c->GetVal()*f1)));
  aHigh = -0.5*m_m0->GetVal()*m_m0->GetVal()*(TMath::Exp(m_c->GetVal()*f2)*TMath::Sqrt(f2)/m_c->GetVal() + 
					      0.5/TMath::Power(-m_c->GetVal(),1.5)*rootPi*TMath::Erf(TMath::Sqrt(-m_c->GetVal()*f2)));
  Double_t area = aHigh - aLow;
  return area;

}

Bool_t PdfArgusBG::evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
				const Double_t invIntegral)
{
  const Data::Value_t *dataCPU = m_data->GetCPUData(*m_m);
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
                 auto om0 = 1./m_m0->GetVal();
                 _Cilk_for (Int_t idx = (Int_t)iPartialStart; idx<(Int_t)iPartialEnd; idx++) {
		   resultsCPU[idx] = evaluateLocal(dataCPU[idx],om0,m_c->GetVal())*invIntegral;
		 }
		 );
  }
  else {

#ifndef USE_CEAN
#pragma ivdep
    auto om0 = 1./m_m0->GetVal();
    for (Int_t idx = (Int_t)iPartialStart; idx<(Int_t)iPartialEnd; idx++) {
      resultsCPU[idx] = evaluateLocal(dataCPU[idx],om0,m_c->GetVal())*invIntegral;
    }
#else
    resultsCPU[iPartialStart:iPartialEnd-iPartialStart] = evaluateLocal(dataCPU[iPartialStart:iPartialEnd-iPartialStart],
									1./m_m0->GetVal(),m_c->GetVal())*invIntegral;
#endif
  }  

  return kTRUE;
}

