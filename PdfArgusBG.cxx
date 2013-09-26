#include "PdfArgusBG.h"

PdfArgusBG::PdfArgusBG(const Char_t* name, const Char_t* title, Variable &m,
		       Variable &m0, Variable &c) :
  AbsPdf(name,title,&m,&m0,&c), m_m(&m), m_m0(&m0), m_c(&c)
{

}

double PdfArgusBG::integral(PdfState const & state) const
{
  constexpr Double_t rootPi = TMath::Sqrt(TMath::Pi());
  Double_t min = std::min(m_m->GetMin() , m_m0->value(state));
  Double_t max = std::max(m_m->GetMax() , m_m0->value(state));
  Double_t f1 = (1.-TMath::Power(min/m_m0->value(state),2));
  Double_t f2 = (1.-TMath::Power(max/m_m0->value(state),2));
  Double_t aLow, aHigh ;
  aLow  = -0.5*m_m0->value(state)*m_m0->value(state)*(TMath::Exp(m_c->value(state)*f1)*TMath::Sqrt(f1)/m_c->value(state) + 
					      0.5/TMath::Power(-m_c->value(state),1.5)*rootPi*TMath::Erf(TMath::Sqrt(-m_c->value(state)*f1)));
  aHigh = -0.5*m_m0->value(state)*m_m0->value(state)*(TMath::Exp(m_c->value(state)*f2)*TMath::Sqrt(f2)/m_c->value(state) + 
					      0.5/TMath::Power(-m_c->value(state),1.5)*rootPi*TMath::Erf(TMath::Sqrt(-m_c->value(state)*f2)));
  Double_t area = aHigh - aLow;
  return area;

}

