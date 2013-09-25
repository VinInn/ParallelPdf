#define private public
#include "PdfReferenceState.h"
#undef private
#include "Data.h"
#include "List.h"
#include "TMath.h"

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


void refresh(const Data & data, int ivar , bool all) {
  auto mpdf = PdfReferenceState::me().pdfs().back();

  std::vector<unsigned short> pdfs; std::vector<unsigned short>  dep;
  PdfReferenceState const & state = PdfReferenceState::me();
  PdfReferenceState::me().refresh(pdfs,dep,ivar, true, all);
  assert(pdfs.size()==dep.size());
  for (auto i=0U; i<pdfs.size(); ++i)
    std::cout << pdfs[i] <<"," << dep[i] <<" ";
  std::cout << std::endl;
  if (!pdfs.empty()) assert(mpdf->num()==pdfs.back());

  // the order is correct...
  for (auto i: pdfs) state.cacheIntegral(i);
  
  auto tot = data.GetEntries();
  alignas(ALIGNMENT) double lres[256];
  double * res=0;
  TMath::IntLog localValue;
  for (auto ie=0U; ie<tot; ie+= 256) {
    auto offset = ie;
    auto bsize = std::min(256U,tot-ie);
    for (auto i: pdfs) state.cachePdf(i,bsize,data,offset);
    res = state.pdfVal(mpdf->num(), lres, bsize,data,offset);
    assert(res==&lres[0]);
    localValue = IntLogAccumulate(localValue, lres, bsize);
  }
  auto ret = -0.693147182464599609375*localValue.value();
  
  ret += mpdf->ExtendedTerm(state,tot);

  std::cout << "result " << ret << std::endl;

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

  refresh(data, -1,true);
  std::cout << std::endl;
  refresh(data,-1,false);
  std::cout << std::endl;
  std::cout << std::endl;


  // this is not the way how it will be done as it is not thread safe...

  auto & vars = PdfReferenceState::me().m_Params;
  for (auto ak=0; ak!=2; ++ak) {
    for (auto i = 0U; i!=vars.size(); ++i) {
      if (vars[i]->isData() || vars[i]->IsConstant()) continue;
      auto v = vars[i]->GetVal();
      auto e = vars[i]->GetError();
      vars[i]->SetVal(v+e);
      std::cout << "var " << i << ":   ";
      refresh(data,-1,false);
      if (ak==1) vars[i]->SetVal(v-e);
    }
    std::cout << std::endl;
  }
  refresh(data,-1,false);
  std::cout << std::endl;
  refresh(data,-1,false);
  std::cout << std::endl;
  std::cout << std::endl;
  for (auto i = 0U; i!=vars.size(); ++i) {
    if (vars[i]->isData() || vars[i]->IsConstant()) continue;
    std::cout << "var " << i << ":   ";
    refresh(data,i,false);
  }
  refresh(data,-1,false);
  std::cout << std::endl;
  refresh(data, -1,true);
  std::cout << std::endl;
  std::cout << std::endl;


  delete model; // sic

  return 0;

};
