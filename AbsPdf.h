#ifndef ABS_PDF
#define ABS_PDF
#include <iostream>
#include <cassert>

#include "Named.h"
#include "TMath.h"
#include "Data.h"

#define RooAbsPdf AbsPdf

class AbsPdf : public Named {
 public:

 
  AbsPdf(const Char_t* name, const Char_t* title): Named(name,title){}
  virtual ~AbsPdf() {}


  virtual void makeCache(unsigned int){}
  virtual int verifyCache(){ return 0;}

  virtual void RandomizeFloatParameters();
  virtual void GetParameters(List<Variable>& parameters) { }


  virtual void GetVal(double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const=0; 
  
  // to be called outside parallel loop 
  virtual void CacheIntegral() {
    m_InvIntegral = 1./integral();
  }

  
  Double_t GetInvIntegral() const { return  m_InvIntegral;}

  virtual Double_t ExtendedTerm(UInt_t observed) const { return .0; }
  virtual Bool_t IsExtended() const { return kFALSE; }
  virtual Double_t ExpectedEvents() const { return .0; }



private:

  virtual Double_t integral() const = 0;

 
  Double_t m_InvIntegral;


};

#endif
