
#ifndef ROO_MINIMIZER_FCN
#define ROO_MINIMIZER_FCN

#include "NLL.h"
#include <iostream>
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Variable.h"
#include "List.h"

class RooMinimizerFcn : public ROOT::Minuit2::FCNBase {

 public:

  RooMinimizerFcn(NLL& nll);
  virtual ~RooMinimizerFcn();

  virtual double Up() const {return _errorDef; } 
  virtual double operator()(const std::vector<double>&) const;
  void SetErrorDef(Double_t errorDef) { _errorDef = errorDef; }
  ROOT::Minuit2::MnUserParameters &Parameters();
  List<Variable>& FloatPdfPars() const { return _floatPdfPars; }
  inline Int_t GetNPar() const { return _floatPdfPars.GetSize(); }

  Double_t& GetMaxFCN() { return _maxFCN; }
  Int_t GetNumInvalidNLL() { return _numBadNLL; }
  inline Int_t NumCallsFCN() const { return _nCalls; }

  void BackProp(const ROOT::Minuit2::MnUserParameters &results);  

private:
  
  NLL* m_nll;
  mutable Double_t _maxFCN;
  mutable Int_t _numBadNLL;
  Double_t _errorDef;
  mutable Int_t _nCalls;

  ROOT::Minuit2::MnUserParameters _pars;
  mutable List<Variable> _floatPdfPars;

};

#endif

