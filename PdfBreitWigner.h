
#ifndef PDFBREITWIGNER
#define PDFBREITWIGNER

#include "AbsPdf.h"
#include "Variable.h"

class PdfBreitWigner : public AbsPdf {
public:
  PdfBreitWigner(const Char_t* name, const Char_t* title, Variable &x,
		 Variable &mu, Variable &width);
  virtual ~PdfBreitWigner() { }

  virtual void GetParameters(List<Variable>& parameters) { parameters.AddElement(*m_mu); parameters.AddElement(*m_width); }
  
 protected:
  virtual Double_t evaluate() const;
  virtual Double_t integral() const;
  
  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			      const Double_t invIntegral);
  
 private:

  inline Double_t evaluateLocal(const Double_t x, const Double_t mu,
				const Double_t width) const {
    return ((Double_t)1.)/(TMath::Power(x-mu,2)+((Double_t)0.25)*width*width);
  }

 private:
  Variable *m_x;
  Variable *m_mu;
  Variable *m_width;

};

#endif
