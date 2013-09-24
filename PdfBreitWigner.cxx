
#include "PdfBreitWigner.h"

#include <iostream>

PdfBreitWigner::PdfBreitWigner(const Char_t* name, const Char_t* title, Variable &x,
			       Variable &mu, Variable &width) :
  AbsPdf(name,title,&x,&mu,&width), m_x(&x), m_mu(&mu), m_width(&width)
{
  
}


double PdfBreitWigner::integral(PdfState const & state) const
{
  Double_t c = 2./m_width->value(state);
  Double_t ret = c*(TMath::ATan(c*(m_x->GetMax()-m_mu->value(state))) - TMath::ATan(c*(m_x->GetMin()-m_mu->value(state))));
  
  return ret;
  
}

