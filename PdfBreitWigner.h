
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
  
 private:

  virtual Double_t integral() const;

  void GetVal(double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    
    Data::Value_t const * __restrict__ ldata = data.GetData(*m_x, dataOffset);
    
    auto invIntegral = GetInvIntegral();
 
    for (auto idx = 0U; idx!=bsize; ++idx) {
      auto x = ldata[idx];
      auto y = evaluateOne(x,m_mu->GetVal(),m_width->GetVal())*invIntegral;
      res[idx] = y;
    }

  }  
  inline Double_t evaluateOne(const Double_t x, const Double_t mu,
				const Double_t width) const {
    auto arg = x-mu;
    return 1./(arg*arg+0.25*width*width);
  }

 private:
  Variable *m_x;
  Variable *m_mu;
  Variable *m_width;

};

#endif
