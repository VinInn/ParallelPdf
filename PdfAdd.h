#ifndef PDFADD
#define PDFADD

#include "AbsPdf.h"
#include "List.h"
#include "TMath.h"

#define RooAddPdf PdfAdd


template<typename T, int N>
struct Add {

  T operator()(T const * v, T const * c, int i, int stride ) const { return (*c)*v[i] + add(v+stride,c+1,i,stride); }
  T operator()(std::initializer_list<T> il) const { return *il.begin()+add(std::begin(il)+1);}
  Add<T,N-1> add;
  
};

template<typename T>
struct Add<T,2> {
  
  T operator()(T x, T y) const { return x+y;}
  T operator()(T const * v, T const * c, int i, int stride ) const { return (*c)*v[i] + (*(c+1))*(v+stride)[i]; }
  T operator()(std::initializer_list<T> il) const { return *il.begin() + *(std::begin(il)+1);}
};


template<int N>
class PdfAdd : public AbsPdf {
public:
  PdfAdd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2, Variable &fraction) :
    AbsPdf(name,title), m_isExtended(kFALSE)
  {
    m_pdfs.AddElement(pdf1);
    m_pdfs.AddElement(pdf2);
    m_fractions.AddElement(fraction);
    
  }
  PdfAdd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs, List<Variable> fractions) :
    AbsPdf(name,title), m_isExtended(kFALSE)
  {
    if (pdfs.GetSize()!=fractions.GetSize() && pdfs.GetSize()!=fractions.GetSize()-1) {
      std::cerr << GetName() << ":: Wrong number of fractions!" << std::endl;
      assert(0);
    }
    
    if (pdfs.GetSize()==fractions.GetSize())
      m_isExtended = kTRUE;
    
    m_pdfs.AddElement(pdfs);
    m_fractions.AddElement(fractions);
    
  }
  
  virtual ~PdfAdd() { }
  
   
  virtual void GetParameters(List<Variable>& parameters) 
  {
    parameters.AddElement(m_fractions);
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->GetParameters(parameters);
    }
  }
  
  virtual void CacheIntegral() {
    AbsPdf::CacheIntegral();
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->CacheIntegral();
    }
    
  }
  

  
private:



  virtual Double_t integral() const { return m_isExtended ? ExpectedEvents() : 1.; }
  

  virtual void GetVal(double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);

    auto strid = stride(bsize);
    alignas(ALIGNMENT) double lres[N][strid];
    double coeff[N];

    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    
    Variable *var(0);
    AbsPdf *pdf = iter_pdfs.Next();
    Double_t lastFraction = 1.;
 
   
    int k=0;
    while ((var = iter_fractions.Next())!=0) {
      lastFraction -= var->GetVal();
      coeff[k]=  var->GetVal();
      pdf->GetVal(lres[k], bsize, data, dataOffset);
      pdf = iter_pdfs.Next();
      ++k;
    }

    if (!m_isExtended) {
      coeff[k]=lastFraction;
      pdf->GetVal(lres[k], bsize, data, dataOffset);
       assert(N==k+1);
    } else 
      assert(N==k);
    
    Add<double,N> add;
    double const *  kres = lres[0];
    auto invIntegral = GetInvIntegral();
    for (auto idx = 0; idx!=bsize; ++idx) {
      res[idx] = add(kres,coeff,idx, strid)*invIntegral;
    }

  }

  virtual Double_t ExtendedTerm(UInt_t observed) const {
    Double_t expected = ExpectedEvents();
    return expected-observed*TMath::Log(expected);
  }

  virtual Bool_t IsExtended() const { return m_isExtended; }

  virtual Double_t ExpectedEvents() const {
    Double_t nEvents(0);
    if (m_isExtended) {
      Variable *var(0);
      List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
      while ((var = iter_fractions.Next())!=0)
	nEvents += var->GetVal();
    }
    
  return nEvents;
}
  
private:
  
  mutable List<AbsPdf> m_pdfs;
  mutable List<Variable> m_fractions;
  
  Bool_t m_isExtended;

 
};

#endif

