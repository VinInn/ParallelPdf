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
  

private:

  virtual Double_t integral() const;


  void GetVal(double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);

    Data::Value_t const * __restrict__ ldata = data.GetData(*m_x, dataOffset);
    
    auto invIntegral = GetInvIntegral();
 
    auto coeff = -0.5/(m_sigma->GetVal()*m_sigma->GetVal());
    for (auto idx = 0U; idx!=bsize; ++idx) {
      auto x = ldata[idx];
      auto y = evaluateOne(x,m_mu->GetVal(),coeff)*invIntegral;
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
