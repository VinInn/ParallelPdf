
#include "PdfBifurGaussian.h"

#include <iostream>

PdfBifurGaussian::PdfBifurGaussian(const Char_t* name, const Char_t* title, Variable &x,
				   Variable &mu, Variable &sigmaL, Variable &sigmaR) :
  AbsPdf(name,title), m_x(&x), m_mu(&mu), m_sigmaL(&sigmaL), m_sigmaR(&sigmaR)
{

}


Double_t PdfBifurGaussian::integral() const
{
  const Double_t root2 = TMath::Sqrt2() ;
  const Double_t rootPiBy2 = TMath::Sqrt(TMath::PiOver2());
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
