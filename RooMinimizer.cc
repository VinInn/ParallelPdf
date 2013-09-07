
#include "RooMinimizer.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnUserParameterState.h"

RooMinimizer::RooMinimizer(NLL& nll) :
  _fcn(new RooMinimizerFcn(nll)),
  _minimumResult(0),
  _status(-1),
  _strategy(1), _eps(1.0)
{
  _maxfcn = 500*_fcn->GetNPar();
  setErrorLevel(0.5) ;

  // Set print level
  ROOT::Minuit2::MnPrint::SetLevel(2);

}

RooMinimizer::~RooMinimizer()
{
  delete _fcn;
  delete _minimumResult;
}

void RooMinimizer::setStrategy(Int_t istrat)
{
  _strategy = istrat;
}

void RooMinimizer::setErrorLevel(Double_t level)
{
  _fcn->SetErrorDef(level);

}

void RooMinimizer::setEps(Double_t eps)
{
  _eps = eps;
}

Int_t RooMinimizer::migrad()
{
  std::cout << "RooMinimizer::migrad: Minimize with max-calls " << _maxfcn 
	    << " convergence for edm < " << _eps 
	    << " strategy " << _strategy << std::endl;

  ROOT::Minuit2::MnMigrad minimizer(*_fcn,_fcn->Parameters(),_strategy);
  delete _minimumResult;
  _minimumResult = new ROOT::Minuit2::FunctionMinimum(minimizer(_maxfcn,_eps));

  // check if Hesse needs to be run 
  if (_minimumResult->IsValid() && 
      _minimumResult->State().Error().Dcovar() != 0 ) {
    // run Hesse (Hesse will add results in the last state of fMinimum
    ROOT::Minuit2::MnHesse hesse(_strategy);
    hesse(*_fcn,*_minimumResult,_maxfcn);
  }

  _fcn->BackProp(_minimumResult->UserParameters());

  _status = (_minimumResult->IsValid() && _minimumResult->HasAccurateCovar() ? 0 : 1);

  return _status ;
}

Int_t RooMinimizer::hesse()
{

  if (_minimumResult==0 && _status!=0) {
    return (_status = -1);
  }
   
  ROOT::Minuit2::MnHesse hesse(_strategy);
  hesse(*_fcn,*_minimumResult,_maxfcn);
  _fcn->BackProp(_minimumResult->UserParameters());

  _status = (_minimumResult->IsValid() && _minimumResult->HasAccurateCovar() ? 0 : 1);

  return _status ;

}

Int_t RooMinimizer::minos()
{

  if (_minimumResult==0 && _status!=0) {
    return (_status = -1);
  }

  ROOT::Minuit2::MnMinos minos(*_fcn,*_minimumResult,_strategy);

  ROOT::Minuit2::MnUserParameters pars(_fcn->Parameters());
  List<Variable> floatPars;
  floatPars.AddElement(_fcn->FloatPdfPars());
  Variable *par(0);
  floatPars.ResetIterator();
  Int_t index(0);
  while ((par = floatPars.Next())!=0) {
    index = pars.Index(par->GetName());
    std::pair<Double_t,Double_t> asymErrors = minos(index,_maxfcn);
    par->SetAsymError(asymErrors);
  }
  _fcn->BackProp(_minimumResult->UserParameters());

  return _status ;
}


