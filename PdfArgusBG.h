#ifndef PDFARGUSBG
#define PDFARGUSBG

#include "AbsPdf.h"
#include "Variable.h"

#define RooArgusBG PdfArgusBG

class PdfArgusBG : public AbsPdf {
public:
  PdfArgusBG(const Char_t* name, const Char_t* title, Variable &m,
	     Variable &m0, Variable &c);

  virtual ~PdfArgusBG() { }

  virtual void GetParameters(List<Variable>& parameters) { parameters.AddElement(*m_m0); parameters.AddElement(*m_c); }

private:

  virtual double integral(PdfState const & state) const;
  
  
  virtual void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);
    
    Data::Value_t const * __restrict__ ldata = data.GetData(*m_m, dataOffset);
    
    auto invInt = invIntegral(state);
 
    auto om0 = 1./m_m0->value(state);
    auto c = m_c->value(state);
    for (auto idx = 0U; idx!=bsize; ++idx) {
      auto x = ldata[idx];
      auto y = evaluateOne(x,om0,c)*invInt;
      res[idx] = y;
    }

  }

  static Double_t evaluateOne(const Double_t m, const Double_t om0,
			      const Double_t c) {
    
    Double_t t= m*om0;
    
    Double_t u=  1. - t*t;
    
    // return (t >= 1.) ? 0. :  m*TMath::Sqrt(u)*TMath::Exp(c*u) ;
    return m*TMath::Sqrt(u)*TMath::Exp(c*u) ;
  }
  
private:
  Variable *m_m;
  Variable *m_m0;
  Variable *m_c;
  
};

#endif
