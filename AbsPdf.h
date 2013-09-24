#ifndef ABS_PDF
#define ABS_PDF

#include <iostream>
#include <cassert>

#include "Named.h"
#include "PdfReferenceState.h"

#include "TMath.h"
#include "Data.h"
#include <initializer_list>

#define RooAbsPdf AbsPdf

class AbsPdf : public Named {
 public:

  template<typename... Args> 
  AbsPdf(const Char_t* name, const Char_t* title, Args... args):
    Named(name,title, pdf){
    PdfReferenceState::registerPdf(this,{args...});
  }
  virtual ~AbsPdf() {}


  virtual void RandomizeFloatParameters();
  virtual void GetParameters(List<Variable>& parameters) { }

  void operator()(PdfState const & state, double * __restrict__ res, double * __restrict__ loc, unsigned int bsize, const Data & data, unsigned int dataOffset) const {
    state.pdfVal(num(),res, loc, bsize,data,dataOffset);
  }
  double invIntegral(PdfState const & state) const { return state.invIntegral(num()); }

  virtual double integral(PdfState const & state) const = 0;
  virtual void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const=0; 



  virtual Double_t ExtendedTerm(UInt_t observed) const { return .0; }
  virtual Bool_t IsExtended() const { return kFALSE; }
  virtual Double_t ExpectedEvents() const { return .0; }

  bool noCache() const { return m_nocache;}

protected:

  bool m_nocache=false;


private:

};



struct NoCacheAbsPdf : public AbsPdf {
  template<typename... Args> 
   NoCacheAbsPdf(Args... args): AbsPdf(args...){ m_nocache=true;}
};

#endif
