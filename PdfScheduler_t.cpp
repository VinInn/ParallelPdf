#define  DOPRINT
#include "PdfScheduler.h"
#undef DOPRINT

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

void refresh(PdfState & state,  int ivar , bool all, bool docache) {
  auto const & pdfsV = PdfReferenceState::me().pdfs();
  auto mpdf = PdfReferenceState::me().pdfs().back();
  auto npdfs = PdfReferenceState::me().pdfs().size();
  std::vector<unsigned short> pdfs; std::vector<unsigned short>  dep;
  if (all) state.allDeps(pdfs,dep, docache);
  else state.deps(pdfs,dep, docache);

  assert(pdfs.size()==dep.size());
  for (auto i=0U; i<pdfs.size(); ++i)
    std::cout << pdfs[i] <<"," << dep[i] <<" ";
  std::cout << std::endl;
  if (!pdfs.empty()) assert(mpdf->num()==pdfs.back());

  for (auto i: pdfs) state.cacheIntegral(i);

  for (auto i=0U; i<npdfs; ++i)
    assert( pdfsV[i]->invIntegral(state) == 1./pdfsV[i]->integral(state));
  // std::cout << i<<':' << pdfsV[i]->invIntegral(state) << ',' << 1./pdfsV[i]->integral(state) << ' ';
  // std::cout << std::endl;


  auto tot = state.data().GetEntries();
  alignas(ALIGNMENT) double lres[256];
  double * res=0;
  TMath::IntLog localValue;
  for (auto ie=0U; ie<tot; ie+= 256) {
    auto offset = ie;
    auto bsize = std::min(256U,tot-ie);
    // the order is correct...
    for (auto i: pdfs) state.cachePdf(i,bsize,offset);
    res = state.value(lres, bsize,offset);
    assert(res==&lres[0]);
    localValue = IntLogAccumulate(localValue, res, bsize);
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
  
  PdfReferenceState & refState = PdfReferenceState::me();

  refState.init(data);
  
  refState.print();

  auto pdfPars = PdfReferenceState::me().variables();
  // first count 
  int nvar=0;
  for ( auto p :  pdfPars ) if (!p->IsConstant() && !p->isData()) ++nvar;
    
  int vars[nvar];
  int k=0; int ip=0;
  for ( auto p :  pdfPars ) { if (!p->IsConstant() && !p->isData()) vars[k++]= ip; ++ip;}
  assert(k==nvar);
 

  // evaluate and cache
  refresh(refState, -1,true,true);
  std::cout << std::endl;

  double steps[2*nvar];

  k=0;
  for (auto i : vars) {
    auto e = pdfPars[i]->GetError();
    steps[k++]=e;
    steps[k++]=-e;
  }
  assert(k==2*nvar);


  double deriv[nvar];
  differentiate(refState,nvar,vars,steps,deriv);

  for ( auto d : deriv) 
    std::cout << d << ' ';
  std::cout << std::endl;


  delete model; // sic
  return 0;

}
