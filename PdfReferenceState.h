#ifndef PdfReferenceState_H
#define PdfReferenceState_H

#include "Named.h"
#include <vector>
#include <initializer_list>
#include "Data.h"


class AbsPdf;
class Variable;
#include "List.h"



/*
 * state of the pdf: it could be either reference or modified
 * a reference state has all chaches initialized
 * a modified state 
 *
 */
class PdfState {
public:
  virtual ~PdfState(){}

  // return value for Paramer i;
  virtual double value(size_t i) const =0;
  // return integral for pdf i;
  virtual double invIntegral(size_t i) const =0;
  // fill res for pdf i;
  virtual void pdfVal(size_t i, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const =0;


};


class PdfReferenceState : public PdfState {

public:

  // return value for Paramer i;
  double value(size_t i) const final { return m_parCache[i];}
  // return integral for pdf i;
  double invIntegral(size_t i) const final { return m_InvIntegrals[i];}
  // fill res for pdf i;
  void pdfVal(size_t i, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const final {
    res = m_resCache.GetData(i,dataOffset);
  }

  AbsPdf * pdf(int i) {return m_pdfs[i];}
  AbsPdf const * pdf(int i) const {return m_pdfs[i];}

  PdfReferenceState() : m_indexDep(1,0), m_indexPdf(1,0), initialized(false){}

  void init();

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
  std::vector<unsigned short> m_indexDep; // index in vector below
  std::vector<short> m_Dep; // direct dependencies

  std::vector<Variable*> m_Params;
  std::vector<unsigned short> m_indexPdf; // index in the vector below...
  std::vector<short> m_PdfsPar;  // pdf corresponding to a par... (-1 is a sum  normalization) 

  std::vector<double> m_parCache; // cache of param 
  mutable Data m_resCache;   // cache of pdfs results
  std::vector<double> m_InvIntegrals; // cache of inverseIntegrals

  bool initialized;

};


class PdfModifiedState  : public PdfState {

public:

  // return value for Paramer i;
  double value(size_t i) const final { return i== m_param ?  m_value : m_reference->value(i);}
  // return integral for pdf i;
  double invIntegral(size_t i) const final { auto k = findPdf(i); return k>=0 ? m_InvIntegrals[k] : m_reference->invIntegral(i); }
  // fill res for pdf i;
  void pdfVal(size_t i, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const final;


private:

  int findPdf(size_t i) const {
    int k = std::find(m_pdfs.begin(),m_pdfs.end(),i)-m_pdfs.begin();
    return (k==int(m_pdfs.size())) ? -1 : k;
  }

  PdfReferenceState const * m_reference;
  
  unsigned short m_param;

  std::vector<unsigned short> m_pdfs;

  double m_value;

  std::vector<double> m_InvIntegrals;

};


#endif
