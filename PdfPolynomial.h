
#ifndef PDFPOLYNOMIAL
#define PDFPOLYNOMIAL

#include "AbsPdf.h"
#include "Variable.h"
#include "List.h"

#include "tbb.h"

#define RooPolynomial PdfPolynomial

class PdfPolynomial : public AbsPdf {
public:
  PdfPolynomial(const Char_t* name, const Char_t* title, const Variable &x);
  PdfPolynomial(const Char_t* name, const Char_t* title, const Variable &x,
		List<Variable> coeff);
  virtual ~PdfPolynomial();

  virtual void GetParameters(List<Variable>& parameters) { parameters.AddElement(m_coeff); }
  static void SetBlockSize(Int_t blockSize) { m_blockSize = blockSize; }

 protected:
  virtual Double_t evaluate() const;
  virtual Double_t integral() const;
  
  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			      const Double_t invIntegral);

 private:

  inline const Double_t *loadCoeff(Double_t *coeffCPU, UInt_t size) const {
    // the coeffCPU must have the correct dimension (size+1)
    coeffCPU[0] = 1.;
    for(UInt_t i = 0; i<size; i++) {
      coeffCPU[i+1] = m_coeff.GetElement(i)->GetVal();
    }
    
    return &coeffCPU[0];

  }

  inline Double_t evaluateLocal(const Double_t x, const Double_t *coeff, UInt_t order) const {

    double result = coeff[order];
    for (;order>0;--order)
      result = result*x+coeff[order-1];

    return result;

  }

  inline Double_t evaluateLocalSingleCoeff(const Double_t x, const Double_t coeff, const Double_t result) const {
    return (result*x+coeff);
  }

 private:
  const Variable *m_x;
  List<Variable> m_coeff;
  static Int_t m_blockSize;
};

#endif
