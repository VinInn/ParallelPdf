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
class PdfReferenceState {

public:

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


class PdfModifiedState {



private:

  PdfReferenceState const * m_reference;
  
  unsigned short m_param;

  std::vector<unsigned short> m_pdfs;

  double m_value;

  std::vector<double> m_InvIntegrals;

};


#endif
