#ifndef PdfReferenceState_H
#define PdfReferenceState_H

#include "PdfState.h"

#include "Named.h"
#include <vector>
#include <initializer_list>


#include "List.h"


#include "Data.h"

class PdfReferenceState : public PdfState {

  friend class PdfScheduler;

public:

  std::vector<AbsPdf*> & pdfs() { return m_pdfs; }
  std::vector<AbsPdf*> const & pdfs() const { return m_pdfs; }

  std::vector<Variable*> & variables() { return m_Params; }
  std::vector<Variable*> const & variables() const { return m_Params; }

  size_t size() const final { return m_InvIntegrals.size();}
 
  // return values for the whole model...
  double * value(double * __restrict__ loc, unsigned int bsize,  unsigned int dataOffset) const final {
    return pdfVal(size()-1,loc, bsize ,dataOffset);
  }



  // return value for Paramer i;
  double paramVal(size_t i) const final { return m_parCache[i];}
  // return integral for pdf i;
  double invIntegral(size_t i) const final { return m_InvIntegrals[i];}
  // fill res for pdf i;
  double * pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize, unsigned int dataOffset) const final;

  void cacheIntegral(size_t i) const final;
  void cachePdf(size_t i, unsigned int bsize,  unsigned int dataOffset) const final;


  void cacheYourIntegral(size_t i) const final {
    return cacheIntegral(i);
  }



  void deps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, bool doCache) const {
    refresh(res,dep,-1, doCache, false);
  }

  void allDeps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, bool doCache) const final {
    refresh(res,dep,-1, doCache, true);
  }

  void depsI(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, int ivar, bool doCache) const {
    refresh(res,dep,ivar,doCache,false);
  }

  void refresh(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, int ivar, bool doCache, bool allPdf) const;

  AbsPdf * pdf(int i) {return m_pdfs[i];}
  AbsPdf const * pdf(int i) const {return m_pdfs[i];}

  PdfReferenceState() : m_indexDep(1,0), m_indexPdf(1,0), initialized(false){}

  void reset() {
    *this = PdfReferenceState();
  }

  void init(const Data & idata);

  static PdfReferenceState & me();
  static void registerPdf(AbsPdf * pdf, std::initializer_list<Named *> pdfOrVar);

  void print() const;

private:

  void registerHere(AbsPdf* pdf, std::initializer_list<Named *> pdfOrVar);
  
  int add(List<Variable> &);
  int add(List<AbsPdf>&);

  int add(Variable *);
  int add(AbsPdf *);



  std::vector<AbsPdf *> m_pdfs;
  std::vector<short> m_indexCache; // some are not in cache....

  std::vector<unsigned short> m_indexDep; // index in vector below
  std::vector<short> m_Dep; // direct dependencies

  std::vector<Variable*> m_Params;
  std::vector<unsigned short> m_indexPdf; // index in the vector below...
  std::vector<short> m_PdfsPar;  // pdf corresponding to a par... 

  mutable std::vector<double> m_parCache; // cache of param 
  mutable Data m_resCache;   // cache of pdfs results
  mutable std::vector<double> m_InvIntegrals; // cache of inverseIntegrals

  bool initialized;

};


class PdfModifiedState  : public PdfState {

public:

  PdfModifiedState(){}
  PdfModifiedState(PdfReferenceState const * ref, unsigned int ipar, double v) :
    m_reference(ref),  m_param(ipar), m_value(v){
    ref->depsI(m_pdfs,m_deps,ipar,false);
    m_InvIntegrals.resize(m_pdfs.size());
    m_data = &ref->data();
  }


  size_t size() const final { return m_InvIntegrals.size();}

  // return values for the whole model...
  double * value(double * __restrict__ loc, unsigned int bsize, unsigned int dataOffset) const final {
    return pdfVal(m_reference->size()-1,loc, bsize, dataOffset);
  }


  unsigned short param() const { return m_param; }


  // return value for Paramer i;
  double paramVal(size_t i) const final { return i== m_param ?  m_value : m_reference->paramVal(i);}
  // return integral for pdf i;
  double invIntegral(size_t i) const final { auto k = findPdf(i); return k>=0 ? m_InvIntegrals[k] : m_reference->invIntegral(i); }
  // fill res for pdf i;
  double * pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize,  unsigned int dataOffset) const final;

  void cacheIntegral(size_t i) const final;
  void cachePdf(size_t i, unsigned int bsize,  unsigned int dataOffset) const final;

  void cacheYourIntegral(size_t i) const final;

  void deps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, bool) const final {
    res = m_pdfs; dep = m_deps;
  }

  void allDeps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, bool) const final {
    m_reference->allDeps(res,dep, true);
  }

private:

  int findPdf(size_t i) const {
    int k = std::find(m_pdfs.begin(),m_pdfs.end(),i)-m_pdfs.begin();
    return (k==int(m_pdfs.size())) ? -1 : k;
  }

  PdfReferenceState const * m_reference=nullptr;

  std::vector<unsigned short> m_pdfs;
  std::vector<unsigned short> m_deps;

  
  unsigned short m_param;

  double m_value;

  mutable std::vector<double> m_InvIntegrals;

};


class PdfNoCacheState  : public PdfState {

public:

  explicit PdfNoCacheState(PdfReferenceState const * ref) :
  m_reference(ref),m_InvIntegrals(ref->pdfs().size()){m_data = &ref->data();}

  size_t size() const final { return m_InvIntegrals.size();}

  // return values for the whole model...
  double * value(double * __restrict__ loc, unsigned int bsize,  unsigned int dataOffset) const final {
    return pdfVal(m_reference->size()-1,loc, bsize,dataOffset);
  }

  // return value for Paramer i;
  double paramVal(size_t i) const final;
  // return integral for pdf i;
  double invIntegral(size_t i) const final { return m_InvIntegrals[i]; }
  // fill res for pdf i;
  double * pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize,  unsigned int dataOffset) const final;

  void cacheIntegral(size_t i) const final;
  void cachePdf(size_t, unsigned int,  unsigned int) const final {}

  void cacheYourIntegral(size_t i) const final {
    return cacheIntegral(i);
  }

 
  void deps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, bool) const final {
    m_reference->allDeps(res, dep, false);
  }

  void allDeps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, bool) const final {
    m_reference->allDeps(res,dep, false);
  }



private:

  PdfReferenceState const * m_reference;
  mutable std::vector<double> m_InvIntegrals;

};

#endif
