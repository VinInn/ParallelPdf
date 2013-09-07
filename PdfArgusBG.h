
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

 protected:
  virtual Double_t evaluate() const;
  virtual Double_t integral() const;

  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			      const Double_t invIntegral);
  
 private:
  inline Double_t evaluateLocal(const Double_t m, const Double_t om0,
				const Double_t c) const {

    Double_t t= m*om0;

    Double_t u= ((Double_t)1.) - t*t;
    
    // return (t >= 1.) ? 0. :  m*TMath::Sqrt(u)*TMath::Exp(c*u) ;
    return m*TMath::Sqrt(u)*TMath::Exp(c*u) ;
  }

 private:
  Variable *m_m;
  Variable *m_m0;
  Variable *m_c;

};

#endif
