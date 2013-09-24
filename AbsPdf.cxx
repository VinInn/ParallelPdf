#include "AbsPdf.h"
#include "TRandom.h"
#include "Variable.h"

void AbsPdf::RandomizeFloatParameters()
{
  List<Variable> pdfPars;
  GetParameters(pdfPars);
  pdfPars.Sort();
  pdfPars.ResetIterator();
  Variable *par(0);
  TRandom rand;
  while ((par = pdfPars.Next())!=0) {
    if (!par->IsConstant()) {
      std::cout << par->GetName() << " = " << par->GetVal();
      par->SetAllVal(rand.Uniform(par->getMin(),par->getMax()));
      std::cout << " --> " << par->GetVal() << std::endl;
    }
  }
}

