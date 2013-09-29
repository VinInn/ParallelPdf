
#include "PdfScheduler.h"


#include "PdfReferenceState.h"

#include "models/model.h"
// #include "models/gauss1.h"


#include "TRandom.h"
#include <cmath>

#include <iostream>

void fillRandom(List<Variable> &  variables, Data & data, unsigned int N) {
  std::cout << "Generate " << N << " events..." << std::endl;
  TRandom rand;
  for (auto i=0U; i<N; i++) {
    variables.ResetIterator();
    while (Variable *var = variables.Next()) {
      var->SetAllVal(rand.Uniform(var->GetMin(),var->GetMax()));
    }
    data.Push_back();
  }
}

AbsPdf *Model(Variable &x, Variable &y, Variable &z, const Int_t N)
{

  // Define the model
  //  AbsPdf *model = Extended1(x,y,z);
  AbsPdf *model = ModelEtapRGKs(x,y,z,N);
  //  AbsPdf *model = Gauss1(x);

  return model;
}


int main() {

  const unsigned int N = 10000;

  // Define the variables
  DataVariable x("x","",-0.2,0.2); // DE
  DataVariable y("y","",5.25,5.29); // mES
  DataVariable z("z","",-3,1.5); // Fisher
  List<Variable> variables(x,y,z);

  // Fill the data
  Data data("data","",N,variables);
  
  fillRandom(variables, data, N);

  auto model = Model(x,y,z,N);
  
  PdfReferenceState & refState = PdfReferenceState::me();

  refState.init(data);
  
  refState.print();


  delete model; // sic
  return 0;

}
