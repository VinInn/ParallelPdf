#include "AbsPdf.h"
#include "TRandom.h"
#include "Variable.h"

void AbsPdf::RandomizeFloatParameters()
{
  List<Variable> pdfPars;
  GetParameters(pdfPars);
  pdfPars.Sort();
  pdfPars.ResetIterator();
  TRandom rand;
  for(auto par : pdfPars() ) {
    if (!par->IsConstant()) {
      std::cout << par->GetName() << " = " << par->GetVal();
      par->SetAllVal(rand.Uniform(par->getMin(),par->getMax()));
      std::cout << " --> " << par->GetVal() << std::endl;
    }
  }
}

