#ifndef PDFGAUSSIAN
#define PDFGAUSSIAN

#include "AbsPdf.h"
#include "Variable.h"

#define RooGaussian PdfGaussian

class PdfGaussian : public AbsPdf {
public:
  PdfGaussian(const Char_t* name, const Char_t* title, Variable &x,
	      Variable &mu, Variable &sigma);
  virtual ~PdfGaussian() { }

  virtual void GetParameters(List<Variable>& parameters) { parameters.AddElement(*m_mu); parameters.AddElement(*m_sigma); }
  


  virtual double integral(PdfState const & state) const;
  
  virtual void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);
    
    Data::Value_t const * __restrict__ ldata = data.GetData(*m_x, dataOffset);
    
    auto invInt = invIntegral(state);
 
    auto coeff = -0.5/(m_sigma->value(state)*m_sigma->value(state));
    auto mu = m_mu->value(state);
    for (auto idx = 0U; idx!=bsize; ++idx) {
      auto x = ldata[idx];
      auto y = evaluateOne(x,mu,coeff)*invInt;
      res[idx] = y;
    }

  }  

  static Double_t evaluateOne(const Double_t x, const Double_t mu,
			      const Double_t coeff)  {
    auto arg = x-mu;
    return TMath::Exp(coeff*arg*arg);
  }

 
  
 private:
  Variable *m_x;
  Variable *m_mu;
  Variable *m_sigma;

};

#endif
