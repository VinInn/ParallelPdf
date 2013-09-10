
#include "RooMinimizerFcn.h"

RooMinimizerFcn::RooMinimizerFcn(NLL& nll) :
  m_nll(&nll),
  _maxFCN(-1e30), _numBadNLL(0),  
  _errorDef(0.5),
  _nCalls(0)
{ 
  List<Variable> pdfPars;
  m_nll->GetPdf()->GetParameters(pdfPars);
  Variable *par(0);
  pdfPars.ResetIterator();
  while ((par = pdfPars.Next())!=0) {
    if (!par->IsConstant()) {
      _floatPdfPars.AddElement(*par);
      _pars.Add(par->GetName(),par->GetVal(),par->GetError(),par->GetMin(),par->GetMax());
    }
  }
  m_nll->makeCache();
  m_nll->GetVal(false); // init cache if needed
}

RooMinimizerFcn::~RooMinimizerFcn()
{

}

ROOT::Minuit2::MnUserParameters &RooMinimizerFcn::Parameters()
{
  Variable *par(0);
  _floatPdfPars.ResetIterator();
  Int_t index(0);
  while ((par = _floatPdfPars.Next())!=0) {
    index = _pars.Index(par->GetName());
    _pars.SetValue(index,par->GetVal());
    _pars.SetError(index,par->GetError());
    _pars.SetLimits(index,par->GetMin(),par->GetMax());
  }
  
  return _pars;

}

double RooMinimizerFcn::operator()(const std::vector<double>& parameters) const
{
  Variable *par(0);
  _floatPdfPars.ResetIterator();
  std::vector<double>::const_iterator iterPars = parameters.begin();
  while ((par = _floatPdfPars.Next())!=0) {
    par->SetAllVal(*(iterPars));
    iterPars++;
  }

  _nCalls++;
  double fvalue = m_nll->GetVal();
  if (isinf(fvalue) || isnan(fvalue)) {
    _numBadNLL++ ;
    return _maxFCN;
  }

  if (fvalue>_maxFCN) {
    _maxFCN = fvalue ;
  }

  return fvalue;

}

void RooMinimizerFcn::BackProp(const ROOT::Minuit2::MnUserParameters &results)
{
  Variable *par(0);
  _floatPdfPars.ResetIterator();
  Int_t index(0);
  while ((par = _floatPdfPars.Next())!=0) {
    index = results.Index(par->GetName());
    par->SetAllVal(results.Value(index));
    par->SetError(results.Error(index));
  }

}

