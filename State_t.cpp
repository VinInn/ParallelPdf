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


double refresh(PdfState & state,  int ivar , bool all, bool docache, bool print=true) {
  auto const & pdfsV = PdfReferenceState::me().pdfs();
  auto mpdf = PdfReferenceState::me().pdfs().back();
  auto npdfs = PdfReferenceState::me().pdfs().size();
  std::vector<unsigned short> pdfs; std::vector<unsigned short>  dep;
  if (all) state.allDeps(pdfs,dep, docache);
  else state.deps(pdfs,dep, docache);

  assert(pdfs.size()==dep.size());
  if (print) {
    for (auto i=0U; i<pdfs.size(); ++i)
      std::cout << pdfs[i] <<"," << dep[i] <<" ";
    std::cout << std::endl;
  }
  if (!pdfs.empty()) assert(mpdf->num()==pdfs.back());

  for (auto i: pdfs) state.cacheIntegral(i);

  for (auto i=0U; i<npdfs; ++i)
    assert( pdfsV[i]->invIntegral(state) == 1./pdfsV[i]->integral(state));
  // std::cout << i<<':' << pdfsV[i]->invIntegral(state) << ',' << 1./pdfsV[i]->integral(state) << ' ';
  // std::cout << std::endl;


  auto tot = state.data().GetEntries();
  alignas(ALIGNMENT) double lres[256];
  double * res=0;
  auto npar = Data::inPart();
  auto ret = mpdf->ExtendedTerm(state,tot);

  size_t ltot=0;
  for (auto ip=0U; ip<npar; ++ip) {
    auto lp = state.data().startP(ip);
    auto ln = state.data().sizeP(ip);
    TMath::IntLog localValue;
    for (auto ie=0U; ie<ln; ie+= 256) {
      auto offset = lp+ie;
      auto bsize = std::min(256UL,ln-ie);
      ltot +=bsize;
      // the order is correct...
      if (ivar<0) for (auto i: pdfs) state.cachePdf(i,bsize,offset);
      res = state.value(lres, bsize,offset);
      assert(res==&lres[0]);
      localValue = IntLogAccumulate(localValue, res, bsize);
    }
    ret += -0.693147182464599609375*localValue.value();
  }
  assert(tot==ltot);

  if (print) std::cout << "result " << ltot << " " << ret << std::endl;
  return ret;
}


int doit (unsigned int N) {


  //  const unsigned int N = 10000;

  std::cout << "Partitions " << Data::inPart() << std::endl;

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


  PdfNoCacheState ncState(&refState);

  refresh(ncState, -1,true,false);
  std::cout << std::endl;


  refresh(refState, -1,true,true);
  std::cout << std::endl;
  refresh(refState,-1,false,true);
  std::cout << std::endl;
  std::cout << std::endl;


  // this is not the way how it will be done as it is not thread safe...

  auto & vars = PdfReferenceState::me().m_Params;
  for (auto ak=0; ak!=4; ++ak) {
    for (auto i = 0U; i!=vars.size(); ++i) {
      if (vars[i]->isData() || vars[i]->IsConstant()) continue;
      auto v = vars[i]->GetVal();
      auto e = vars[i]->GetError();
      vars[i]->SetVal(v+e);
      std::cout << "var " << i << ":   ";
      if (ak>1) 
	refresh(ncState, -1,true,false);
      else
	refresh(refState,-1,false, true);
      if (ak%2==1) vars[i]->SetVal(v-e);
    }
    std::cout << std::endl;
  }
  refresh(refState,-1,false, true);
  std::cout << std::endl;
  refresh(refState,-1,false,true);
  std::cout << std::endl;
  std::cout << std::endl;

  for (auto i = 0U; i!=vars.size(); ++i) {
    if (vars[i]->isData() || vars[i]->IsConstant()) continue;
    auto v = vars[i]->GetVal();
    auto e = vars[i]->GetError();
    PdfModifiedState mstate1(&refState,i, v+e);
    std::cout << "var " << i << ":   ";
    refresh(mstate1, i,false,false);
  }


  std::cout << "\nderivs:\n";

  for (auto ak=0; ak!=2; ++ak) {
    for (auto i = 0U; i!=vars.size(); ++i) {
      if (vars[i]->isData() || vars[i]->IsConstant()) continue;
      auto v = vars[i]->GetVal();
      auto e = vars[i]->GetError();
      vars[i]->SetVal(v+e);
      auto p =  (ak>0) ?
	refresh(ncState, -1,true,false,false) :
	refresh(refState,-1,false, true,false);
      vars[i]->SetVal(v-e);
      auto n =  (ak>0) ?
	refresh(ncState, -1,true,false,false) :
	refresh(refState,-1,false, true,false);
      vars[i]->SetVal(v);
      std::cout << (p-n)/(2*e) << ", ";
    }
    std::cout << std::endl;
  }

  refresh(refState,-1,true,true,false);


  for (auto i = 0U; i!=vars.size(); ++i) {
    if (vars[i]->isData() || vars[i]->IsConstant()) continue;
    auto v = vars[i]->GetVal();
    auto e = vars[i]->GetError();
    PdfModifiedState mstate1(&refState,i, v+e);
    auto p = refresh(mstate1, i,false,false,false);
    PdfModifiedState mstate2(&refState,i, v-e);
    auto n = refresh(mstate2, i,false,false,false);
    std::cout << (p-n)/(2*e) << ", ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  refresh(refState,-1,false,true);
  std::cout << std::endl;
  refresh(refState,-1,true,true);
  std::cout << std::endl;
  std::cout << std::endl;


  delete model; // sic

  return 0;

};


int main() {

  // Data::nPartions=2;
  const unsigned int N = 20000;
  
  doit(N);
  
  
  Data::nPartions=2;
  
  PdfReferenceState::me().reset();
  doit(N);

}
