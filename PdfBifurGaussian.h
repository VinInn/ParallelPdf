
#ifndef PDFBIFURGAUSSIAN
#define PDFBIFURGAUSSIAN

#include "AbsPdf.h"
#include "Variable.h"

#define RooBifurGauss PdfBifurGaussian

class PdfBifurGaussian : public AbsPdf {
 public:
  PdfBifurGaussian(const Char_t* name, const Char_t* title, Variable &x,
		   Variable &mu, Variable &sigmaL, Variable &sigmaR);
  virtual ~PdfBifurGaussian() { }

  virtual void GetParameters(List<Variable>& parameters) { parameters.AddElement(*m_mu); 
    parameters.AddElement(*m_sigmaL); parameters.AddElement(*m_sigmaR); }

  virtual double integral(PdfState const & state) const;
  
    virtual void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);
    
    Data::Value_t const * __restrict__ ldata = data.GetData(*m_x, dataOffset);
    

    auto invInt = invIntegral(state);

    auto coeffL = -0.5/(m_sigmaL->value(state)*m_sigmaL->value(state));
    auto coeffR = -0.5/(m_sigmaR->value(state)*m_sigmaR->value(state));
    auto mu = m_mu->value(state);
    for (auto idx = 0U; idx!=bsize; ++idx) {
      auto x = ldata[idx];
      auto y = evaluateOne(x,mu,coeffL,coeffR)*invInt;
      res[idx] = y;
    }

  }  
 

  static Double_t evaluateOne(const Double_t x, const Double_t mu,
			      const Double_t coeffL, const Double_t coeffR) {

   Double_t arg = x - mu;
   Double_t coef = (arg<0) ? coeffL : coeffR;

   return TMath::Exp(coef*arg*arg);
  }



 private:
  Variable *m_x;
  Variable *m_mu;
  Variable *m_sigmaL;
  Variable *m_sigmaR;

};

#endif
