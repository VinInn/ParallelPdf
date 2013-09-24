#define private public
#include "PdfReferenceState.h"
#undef private
#include "Data.h"
#include "List.h"


// #include "models/extended1.h"
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


void refresh(int ivar , bool force) {
  std::vector<unsigned short> res; std::vector<unsigned short>  dep;
  PdfReferenceState::me().refresh(res,dep,ivar,force);
  assert(res.size()==dep.size());
  for (auto i=0U; i<res.size(); ++i)
    std::cout << res[i] <<"," << dep[i] <<" ";
  std::cout << std::endl;
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
  
  PdfReferenceState::me().init(data.GetEntries());
  
  PdfReferenceState::me().print();

  refresh(-1,true);

  // this is not the way how it will be done as it is not thread safe...

  auto & vars = PdfReferenceState::me().m_Params;

  for (auto i = 0U; i!=vars.size(); ++i) {
    if (vars[i]->isData() || vars[i]->IsConstant()) continue;
    auto v = vars[i]->GetVal();
    auto e = vars[i]->GetError();
    vars[i]->SetVal(v+e);
    std::cout << "var " << i << ":   ";
    refresh(-1,false);
  }
  std::cout << std::endl;
  std::cout << std::endl;

  for (auto i = 0U; i!=vars.size(); ++i) {
    if (vars[i]->isData() || vars[i]->IsConstant()) continue;
    std::cout << "var " << i << ":   ";
    refresh(i,true);
  }
  std::cout << std::endl;
  std::cout << std::endl;


  delete model; // sic

  return 0;

};
