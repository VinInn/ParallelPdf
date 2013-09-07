
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
  
 protected:
  virtual Double_t evaluate() const;
  virtual Double_t integral() const;
  
  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			      const Double_t invIntegral);
  
 private:
  
  inline Double_t evaluateLocal(const Double_t x, const Double_t mu,
				const Double_t sigmaL, const Double_t sigmaR) const {

    Double_t arg = x - mu;
    Double_t coef = 0.0;
    
    if (arg < coef) {
      if (TMath::Abs(sigmaL)>1e-30) {
        coef = -((Double_t)0.5)/(sigmaL*sigmaL);
      }
    } else {
      if (TMath::Abs(sigmaR)>1e-30) {
        coef = -((Double_t)0.5)/(sigmaR*sigmaR);
      }
    }
    
    return TMath::Exp(coef*arg*arg);
  }

 private:
  Variable *m_x;
  Variable *m_mu;
  Variable *m_sigmaL;
  Variable *m_sigmaR;

};

#endif
