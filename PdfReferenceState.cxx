#include "PdfReferenceState.h"
#include "AbsPdf.h"
#include "Variable.h"
#include "List.h"
#include <set>
#include <tuple>

#include<iostream>

double PdfNoCacheState::paramVal(size_t i) const { return  m_reference->variables()[i]->GetVal();}


double * PdfNoCacheState::pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize,  unsigned int dataOffset) const {
  m_reference->pdf(i)->values(*this,loc,bsize,data(),dataOffset); return loc;
}

double * PdfModifiedState::pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize, unsigned int dataOffset) const {
  auto k = findPdf(i);
  if (k<0) return m_reference->pdfVal(i,loc,bsize,dataOffset);
  else { m_reference->pdf(i)->values(*this,loc,bsize,data(),dataOffset); return loc;}
}

double * PdfReferenceState::pdfVal(size_t i,  double * __restrict__ loc, unsigned int bsize,  unsigned int dataOffset) const {
  auto k = m_indexCache[i];
  if (k<0) { /* std::cout << "eval " << i << std::endl; */ pdf(i)->values(*this,loc,bsize,data(),dataOffset); return loc; }
  else {  /* std::cout << "from chache " << i << std::endl; */ return m_resCache.GetData(k,dataOffset); }
}

void PdfModifiedState::cacheIntegral(size_t i) const {
  auto k = findPdf(i);
  if (k<0) m_reference->cacheIntegral(i);
  else m_InvIntegrals[k] = 1./m_reference->pdf(i)->integral(*this);
}


void PdfModifiedState::cacheYourIntegral(size_t i) const {
  m_InvIntegrals[i] = 1./m_reference->pdf(m_pdfs[i])->integral(*this);
}


void PdfReferenceState::cacheIntegral(size_t i) const {
  // std::cout << "cashing integral" << i << std::endl;
  m_InvIntegrals[i] = 1./pdf(i)->integral(*this);
}

void PdfNoCacheState::cacheIntegral(size_t i) const {
  // std::cout << "cashing integral" << i << std::endl;
  m_InvIntegrals[i] = 1./m_reference->pdf(i)->integral(*this);
}


void PdfModifiedState::cachePdf(size_t i, unsigned int bsize,  unsigned int dataOffset) const {
  auto k = findPdf(i);
  if (k<0) m_reference->cachePdf(i,bsize,dataOffset);
}

void PdfReferenceState::cachePdf(size_t i, unsigned int bsize, unsigned int dataOffset) const {
  auto k = m_indexCache[i];
  if (k<0) return;
  // std::cout << "cashing" << i << std::endl;
  auto res = m_resCache.GetData(k,dataOffset);
  pdf(i)->values(*this,res,bsize,data(),dataOffset);
}




void PdfReferenceState::refresh(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, int ivar, bool doCache, bool allPdf) const {

  struct ToRefresh {
    ToRefresh(){}
    ToRefresh(int i1, int i2) : ind(i1), deps(i2){}
    int ind=-1;
    mutable int deps=0;  // because of set...
    void incr() const { ++deps;}
    bool operator<(ToRefresh rh) const { return ind<rh.ind;}
  };
  std::set<ToRefresh> toRefresh;
  auto iter = toRefresh.end();
  bool ok=false;
  
  std::function<void(unsigned int)> walk;
  walk = [&](unsigned int kp) {
    for (auto d=m_indexDep[kp]; d!=m_indexDep[kp+1]; ++d) {
      auto kk = m_Dep[d];
      std::tie(iter, ok) = toRefresh.emplace(kk,2); // +1
      if (ok) walk(kk);
      else (*iter).incr();
    }
  };

  auto refreshI = [&](unsigned int i) {
    for (auto k=m_indexPdf[i]; k!=m_indexPdf[i+1]; ++k) {
      auto kp = m_PdfsPar[k];
      std::tie(iter, ok) = toRefresh.emplace(kp,1); // depends on integral
      if (ok) walk(kp);
    }
  };

  if (ivar>=0) refreshI(ivar);
  else if (allPdf) {
    if (doCache) for (auto i = 0U; i!=m_Params.size(); ++i) {
	if (m_Params[i]->isData()) continue;
	m_parCache[i]=m_Params[i]->GetVal();
      }
    for (auto k=0U; k!=m_pdfs.size(); ++k) {
      std::tie(iter, ok) = toRefresh.emplace(k,1); // depends on integral
      if (ok) walk(k);
    }
  } else {
    for (auto i = 0U; i!=m_Params.size(); ++i) {
      if (m_Params[i]->isData()) continue;
      if (m_parCache[i]!=m_Params[i]->GetVal()) {
	if (doCache) m_parCache[i]=m_Params[i]->GetVal();
	refreshI(i);
      }
    }  
  }
  for ( auto s : toRefresh) { res.push_back(s.ind); dep.push_back(s.deps);}

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

void PdfReferenceState::init(const Data & idata) {
  assert(!initialized);
  initialized=true;

  m_data = &idata;

  // invert dependency vectors....
  invert(m_pdfs.size(),m_indexDep, m_Dep);
  invert(m_Params.size(),m_indexPdf,m_PdfsPar);

  m_parCache.resize(m_Params.size());
  m_InvIntegrals.resize(m_pdfs.size());

  m_indexCache.resize(m_pdfs.size(),-1);


  auto k=0U; auto i=0U;
  for ( auto p : m_pdfs) { 
    if ( !p->noCache() ) m_indexCache[i] = k++;
    ++i;
  }

  m_resCache = std::move(Data("","",m_data->GetEntries(),k));

  for (auto i = 0U; i!=m_Params.size(); ++i)
    m_parCache[i]=m_Params[i]->GetVal();
 
}

#include<iostream>
void PdfReferenceState::print() const {
  std::cout << std::endl;

  for (auto k=0U; k!=m_pdfs.size(); ++k) {
    std::cout << m_pdfs[k]->num() <<"," <<  m_pdfs[k]->name() << ": ";
    std::cout << ( m_indexCache[k]<0 ? 'n' : 'c' ) << m_indexCache[k] << ' ';
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

  std::cout << "cache capacity " << m_resCache.capacity()  << std::endl;
}
