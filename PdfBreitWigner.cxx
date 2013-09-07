
#include "PdfBreitWigner.h"

#include <iostream>

PdfBreitWigner::PdfBreitWigner(const Char_t* name, const Char_t* title, Variable &x,
			       Variable &mu, Variable &width) :
  AbsPdf(name,title), m_x(&x), m_mu(&mu), m_width(&width)
{
  
}


Double_t PdfBreitWigner::integral() const
{
  Double_t c = 2./m_width->GetVal();
  Double_t ret = c*(TMath::ATan(c*(m_x->GetMax()-m_mu->GetVal())) - TMath::ATan(c*(m_x->GetMin()-m_mu->GetVal())));
  
  return ret;
  
}

