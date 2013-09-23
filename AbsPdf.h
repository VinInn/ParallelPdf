#ifndef ABS_PDF
#define ABS_PDF
#include <iostream>
#include <cassert>

#include "Named.h"
#include "PdfReferenceState.h"

#include "TMath.h"
#include "Data.h"
#include "omp.h"
#include <vector>
#include <initializer_list>

#define RooAbsPdf AbsPdf

class AbsPdf : public Named {
 public:

  template<typename... Args> 
  AbsPdf(const Char_t* name, const Char_t* title, Args... args):
    Named(name,title, pdf), m_InvIntegral(omp_get_max_threads()){
    PdfReferenceState::registerPdf(this,{args...});
  }
  virtual ~AbsPdf() {}

  // FIXME all this needs a clean up
  virtual void makeCache(unsigned int){}
  virtual int verifyCache(bool){ return 0;}
  virtual unsigned int cacheSize() const { return 0;}


  virtual void RandomizeFloatParameters();
  virtual void GetParameters(List<Variable>& parameters) { }


  virtual void GetVal(double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const=0; 
  
  // to be called outside parallel loop 
  virtual void CacheAllIntegral() {
    auto li = 1./integral();
    for (auto & iI :  m_InvIntegral) iI = li;
  }

  // to be called inside parallel loop 
  virtual void CacheIntegral(int lpar=-2) {
    m_InvIntegral[omp_get_thread_num()] = 1./integral();
  }

  
  Double_t GetInvIntegral() const { return  m_InvIntegral[omp_get_thread_num()];}

  virtual Double_t ExtendedTerm(UInt_t observed) const { return .0; }
  virtual Bool_t IsExtended() const { return kFALSE; }
  virtual Double_t ExpectedEvents() const { return .0; }

  bool noCache() const { return m_nocache;}

protected:

  bool m_nocache=false;


private:

  virtual Double_t integral() const = 0;

 
  std::vector<double> m_InvIntegral;


};



struct NoCacheAbsPdf : public AbsPdf {
  template<typename... Args> 
   NoCacheAbsPdf(Args... args): AbsPdf(args...){ m_nocache=true;}
};

#endif
