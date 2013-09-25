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


  // return value for Paramer i;
  double paramVal(size_t i) const final { return m_parCache[i];}
  // return integral for pdf i;
  double invIntegral(size_t i) const final { return m_InvIntegrals[i];}
  // fill res for pdf i;
  double * pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize, const Data & data, unsigned int dataOffset) const final;

  void cacheIntegral(size_t i) const final;
  void cachePdf(size_t i, unsigned int bsize, const Data & data, unsigned int dataOffset) const final;


  void allDeps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, bool doCache) {
    refresh(res,dep,-1, doCache, true);
  }

  void deps(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, int ivar, bool doCache) {
    refresh(res,dep,ivar,doCache,false);
  }

  void refresh(std::vector<unsigned short> & res, std::vector<unsigned short> & dep, int ivar, bool doCache, bool allPdf);

  AbsPdf * pdf(int i) {return m_pdfs[i];}
  AbsPdf const * pdf(int i) const {return m_pdfs[i];}

  PdfReferenceState() : m_indexDep(1,0), m_indexPdf(1,0), initialized(false){}

  void init(int size);

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

  std::vector<double> m_parCache; // cache of param 
  mutable Data m_resCache;   // cache of pdfs results
  mutable std::vector<double> m_InvIntegrals; // cache of inverseIntegrals

  bool initialized;

};


class PdfModifiedState  : public PdfState {

public:

  PdfModifiedState(PdfReferenceState const * ref, unsigned int ipar, double v, std::vector<unsigned short> const & ipdfs) :
    m_reference(ref),  m_param(ipar), m_pdfs(ipdfs), m_value(v){}


  // return value for Paramer i;
  double paramVal(size_t i) const final { return i== m_param ?  m_value : m_reference->paramVal(i);}
  // return integral for pdf i;
  double invIntegral(size_t i) const final { auto k = findPdf(i); return k>=0 ? m_InvIntegrals[k] : m_reference->invIntegral(i); }
  // fill res for pdf i;
  double * pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize, const Data & data, unsigned int dataOffset) const final;

  void cacheIntegral(size_t i) const final;
  void cachePdf(size_t i, unsigned int bsize, const Data & data, unsigned int dataOffset) const final;


private:

  int findPdf(size_t i) const {
    int k = std::find(m_pdfs.begin(),m_pdfs.end(),i)-m_pdfs.begin();
    return (k==int(m_pdfs.size())) ? -1 : k;
  }

  PdfReferenceState const * m_reference;
  
  unsigned short m_param;

  std::vector<unsigned short> m_pdfs;

  double m_value;

  mutable std::vector<double> m_InvIntegrals;

};


class PdfNoCacheState  : public PdfState {

public:

  explicit PdfNoCacheState(PdfReferenceState const * ref) :
    m_reference(ref){}


  // return value for Paramer i;
  double paramVal(size_t i) const final { return  m_reference->paramVal(i);}
  // return integral for pdf i;
  double invIntegral(size_t i) const final { return  m_reference->invIntegral(i); }
  // fill res for pdf i;
  double * pdfVal(size_t i, double * __restrict__ loc, unsigned int bsize, const Data & data, unsigned int dataOffset) const final;

  void cacheIntegral(size_t i) const final { m_reference->cacheIntegral(i); }
  void cachePdf(size_t, unsigned int, const Data &, unsigned int) const final {}

private:

  PdfReferenceState const * m_reference;

};

#endif
