#ifndef ROO_MINIMIZER
#define ROO_MINIMIZER

#include "NLL.h"
#include "Named.h"
#include "TMath.h"

#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "RooMinimizerFcn.h"

class RooMinimizer {
public:

  RooMinimizer(NLL& nll) ;
  virtual ~RooMinimizer() ;

  void setStrategy(Int_t strat) ;
  void setErrorLevel(Double_t level) ;
  void setEps(Double_t eps) ;

  inline Int_t NumCallsFCN() const { return _fcn->NumCallsFCN(); }
  inline Int_t NumInvalidNLL() const { return _fcn->GetNumInvalidNLL(); }
  inline Double_t MinFCN() const { return _minimumResult!=0 ? _minimumResult->Fval() : -1; }
  inline Double_t Edm() const { return _minimumResult!=0 ? _minimumResult->Edm() : -1; }

  Int_t migrad() ;
  Int_t hesse() ;
  Int_t minos() ;

private:

  RooMinimizerFcn* _fcn;
  ROOT::Minuit2::FunctionMinimum *_minimumResult;
  Int_t _status;
  Int_t _strategy;
  Int_t _eps;
  Int_t _maxfcn;

} ;


#endif

