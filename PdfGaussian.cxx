#include "PdfGaussian.h"

#include <iostream>

PdfGaussian::PdfGaussian(const Char_t* name, const Char_t* title, Variable &x,
			 Variable &mu, Variable &sigma) :
  AbsPdf(name,title), m_x(&x), m_mu(&mu), m_sigma(&sigma)
{
  
}

Double_t PdfGaussian::integral() const
{
  const Double_t root2 = TMath::Sqrt2() ;
  const Double_t rootPiBy2 = TMath::Sqrt(TMath::PiOver2());
  Double_t invxscale = 1./(root2*m_sigma->GetVal());
  Double_t ret = rootPiBy2*m_sigma->GetVal()*
    (TMath::Erf((m_x->GetMax()-m_mu->GetVal())*invxscale)-
     TMath::Erf((m_x->GetMin()-m_mu->GetVal())*invxscale));

  return ret;
  
}

