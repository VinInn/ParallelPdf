#include "PdfReferenceState.h"
#include "AbsPdf.h"
#include "Variable.h"
#include "List"



void PdfModifiedState::pdfVal(size_t i, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const {
  auto k = findPdf(i);
  if (k<0) m_reference->pdfVal(i,res,bsize,data,dataOffset);
  else m_reference->pdf(i)->GetVal(res,bsize,data,dataOffset);
}

PdfReferenceState & PdfReferenceState::me() {
  static PdfReferenceState local;
  return local;
}


void PdfReferenceState::registerPdf(AbsPdf* pdf, std::initializer_list<Named *> pdfOrVar) {
  me().registerHere(pdf,pdfOrVar);
}


int PdfReferenceState::add(List<Variable> & vars) {
  auto n=0;
  for (auto v : vars()) n+=add(v);
  return n;
}

int PdfReferenceState::add(List<AbsPdf> & pdfs) {
  auto n=0;
  for (auto v : pdfs()) n+=add(v);
  return n;

}

int PdfReferenceState::add(Variable * lvar) {
  int k = find(m_Params.begin(),m_Params.end(),lvar)-m_Params.begin();
  if (k==int(m_Params.size())) {  lvar->setNum(m_Params.size()) ;m_Params.push_back(lvar);}
  m_PdfsPar.push_back(k);
  return 1;
}

  int PdfReferenceState::add(AbsPdf * lpdf) {
    auto k = find(m_pdfs.begin(),m_pdfs.end(),lpdf)-m_pdfs.begin();
    m_Dep.push_back(k);
    return 1;
  }


void PdfReferenceState::registerHere(AbsPdf* pdf, std::initializer_list<Named *> pdfOrVar) {

  assert(!initialized);

  pdf->setNum(m_pdfs.size());
  m_pdfs.push_back(pdf);


  auto nP=0; auto nV=0;
  for ( auto elem : pdfOrVar) {
    switch (elem->who()) {
    case Named::list :
      {
      auto pl = dynamic_cast<List<AbsPdf>*>(elem);
      if (pl) nP+=add(*pl);
      else {
	auto vl = dynamic_cast<List<Variable>*>(elem);
	if (vl) nV+=add(*vl);
      }
      break;
      }
    case Named::pdf :
      nP += add(reinterpret_cast<AbsPdf*>(elem));
      break;
    case Named::var :
      nV+=add(reinterpret_cast<Variable*>(elem));
      break;
    case Named::unknown :
      // error
      break;
    }

  }
  auto o1 = m_indexDep.back();
  auto o2 = m_indexPdf.back();

  m_indexDep.push_back(m_Dep.size());
  m_indexPdf.push_back(m_PdfsPar.size());

  assert( (m_indexDep.back()-o1) == nP);
  assert( (m_indexPdf.back()-o2) == nV);

}

namespace {
  void invert(size_t N, std::vector<unsigned short> & index,  std::vector<short> & list) {
    if (list.empty()) return;
    assert(index.size()>1);
    std::vector<std::vector<unsigned short> > direct(N);
    for (auto k=1U; k!=index.size(); ++k) {
      for (auto i=index[k-1]; i!=index[k]; ++i) {
	assert(list[i]<int(N));
	direct[list[i]].push_back(k-1);
      }
    }
    index.clear();
    list.clear();
    index.push_back(0);
    for (auto & v : direct) {
      index.push_back(index.back()+v.size());
      for ( auto i : v) list.push_back(i);
    }
    assert(index.size()==N+1);
    assert(list.size()==index.back());
    
  }
}

void PdfReferenceState::init(int size) {
  assert(!initialized);
  initialized=true;

  // invert dependency vectors....
  invert(m_pdfs.size(),m_indexDep, m_Dep);
  invert(m_Params.size(),m_indexPdf,m_PdfsPar);

  m_parCache.resize(m_Params.size());
  m_InvIntegrals.resize(m_pdfs.size());

  m_indexCache.resize(m_pdfs.size(),-1);


  auto k=0U; auto i=0U;
  for ( auto p : m_pdfs) { 
    if ( p->noCache() ) m_indexCache[i] = k++;
    ++i;
  }

  m_resCache = std::move(Data("","",size,k));
  for (auto i = 0U; i!=m_Params.size(); ++i)
    m_parCache[i]=m_Params[i]->GetVal();
 

}

#include<iostream>
void PdfReferenceState::print() const {
  std::cout << std::endl;

  for (auto k=0U; k!=m_pdfs.size(); ++k) {
    std::cout << m_pdfs[k]->num() <<"," <<  m_pdfs[k]->name() << ": ";
      for (auto i=m_indexDep[k]; i!=m_indexDep[k+1]; ++i)
	std::cout <<  m_pdfs[m_Dep[i]]->name() <<", ";
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for (auto k=0U; k!=m_Params.size(); ++k) {
    std::cout <<  m_Params[k]->num() <<"," << m_Params[k]->name() << ": ";
      for (auto i=m_indexPdf[k]; i!=m_indexPdf[k+1]; ++i)
	std::cout <<  m_pdfs[m_PdfsPar[i]]->name() <<", ";
    std::cout << std::endl;
  }
  std::cout << std::endl;


}
