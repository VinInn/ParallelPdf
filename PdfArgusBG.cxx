#include "PdfArgusBG.h"

PdfArgusBG::PdfArgusBG(const Char_t* name, const Char_t* title, Variable &m,
		       Variable &m0, Variable &c) :
  AbsPdf(name,title,&m,&m0,&c), m_m(&m), m_m0(&m0), m_c(&c)
{

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

