
#include "PdfBifurGaussian.h"

#include <iostream>

PdfBifurGaussian::PdfBifurGaussian(const Char_t* name, const Char_t* title, Variable &x,
				   Variable &mu, Variable &sigmaL, Variable &sigmaR) :
  AbsPdf(name,title,&x,&mu,&sigmaL,&sigmaR), 
  m_x(&x), m_mu(&mu), m_sigmaL(&sigmaL), m_sigmaR(&sigmaR)
{

}


Double_t PdfBifurGaussian::integral(PdfState const & state) const
{
  constexpr Double_t root2 = TMath::Sqrt2() ;
  constexpr Double_t rootPiBy2 = TMath::Sqrt(TMath::PiOver2());
  Double_t invxscaleL = 1./(root2*m_sigmaL->value(state));
  Double_t invxscaleR = 1./(root2*m_sigmaR->value(state));

  Double_t integral = 0.0;
  if(m_x->GetMax() < m_mu->value(state))	{
    integral = m_sigmaL->value(state)*(TMath::Erf((m_x->GetMax()-m_mu->value(state))*invxscaleL)-TMath::Erf((m_x->GetMin()-m_mu->value(state))*invxscaleL));
  }
  else if (m_x->GetMin() > m_mu->value(state)) {
    integral = m_sigmaR->GetVal()*(TMath::Erf((m_x->GetMax()-m_mu->value(state))*invxscaleR)-TMath::Erf((m_x->GetMin()-m_mu->value(state))*invxscaleR));
  }
  else {
    integral = m_sigmaR->GetVal()*TMath::Erf((m_x->GetMax()-m_mu->value(state))*invxscaleR)-m_sigmaL->value(state)*TMath::Erf((m_x->GetMin()-m_mu->value(state))*invxscaleL);
  }

  return integral*rootPiBy2;

}
