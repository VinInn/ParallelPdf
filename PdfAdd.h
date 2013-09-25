#ifndef PDFADD
#define PDFADD

#include "AbsPdf.h"
#include "List.h"
#include "TMath.h"

#define RooAddPdf PdfAdd


template<typename T, int N>
struct Add {
  T operator()(T const *  __restrict__ const * v, T const * __restrict__ c, int i ) const { return (*c)*(*v)[i] + add(v+1,c+1,i); }
  T operator()(T const * v, T const * c, int i, int stride ) const { return (*c)*v[i] + add(v+stride,c+1,i,stride); }
  T operator()(std::initializer_list<T> il) const { return *il.begin()+add(std::begin(il)+1);}
  Add<T,N-1> add;
  
};

template<typename T>
struct Add<T,2> {
  
  T operator()(T x, T y) const { return x+y;}
  T operator()(T const *  __restrict__ const * v, T const * __restrict__ c, int i ) const { return (*c)*(*v)[i] + (*(c+1))*(*(v+1))[i]; }
  T operator()(T const * v, T const * c, int i, int stride ) const { return (*c)*v[i] + (*(c+1))*(v+stride)[i]; }
  T operator()(std::initializer_list<T> il) const { return *il.begin() + *(std::begin(il)+1);}
};


template<int N>
class PdfAdd : public NoCacheAbsPdf {
public:
  PdfAdd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2, Variable &fraction) :
    NoCacheAbsPdf(name,title,&pdf1,&pdf2,&fraction), m_isExtended(kFALSE)
  {
    m_pdfs.AddElement(pdf1);
    m_pdfs.AddElement(pdf2);
    m_fractions.AddElement(fraction);
    
  }
  PdfAdd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs, List<Variable> fractions) :
    NoCacheAbsPdf(name,title,&pdfs,&fractions), m_isExtended(kFALSE)
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
  
 
  
  virtual double integral(PdfState const & state) const { return m_isExtended ? ExpectedEvents(state) : 1.; }
  
  
  virtual void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);

    auto strid = stride(bsize);
    double * __restrict__ pres[N];
    alignas(ALIGNMENT) double lres[N][strid];
    double coeff[N];

    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    
    Variable *var(0);
    Double_t lastFraction = 1.;
 
    int k=0;
    while ((var = iter_fractions.Next())!=0) {
      lastFraction -= var->value(state);
      coeff[k++]=  var->value(state);
    }

    if (!m_isExtended) {
      coeff[k]=lastFraction;
      assert(N==k+1);
    } else 
      assert(N==k);
 
    for (int l=0; l!=N; ++l) {
      auto pdf = m_pdfs()[l];
      (*pdf)(state,pres[l], &(lres[l][0]), bsize, data, dataOffset);
    }
  
    
    Add<double,N> add;
   auto invInt = invIntegral(state);
    double const * __restrict__  const *  kres = pres;
    // double const * kres = lres[0];
#pragma omp simd
    for (auto idx = 0; idx<bsize; ++idx) {
      // res[idx] = add(kres,coeff,idx,strid)*invIntegral;
      res[idx] = add(kres,coeff,idx)*invInt;
    }

  }


  virtual Double_t ExtendedTerm(PdfState const & state, UInt_t observed) const {
    Double_t expected = ExpectedEvents(state);
    return expected-observed*TMath::Log(expected);
  }

  virtual Bool_t IsExtended() const { return m_isExtended; }

  virtual Double_t ExpectedEvents(PdfState const & state) const {
    Double_t nEvents(0);
    if (m_isExtended) {
      Variable *var(0);
      List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
      while ((var = iter_fractions.Next())!=0)
	nEvents += var->value(state);
    }
    
  return nEvents;
}
 

  
  mutable List<AbsPdf> m_pdfs;
  mutable List<Variable> m_fractions;
  
  Bool_t m_isExtended;

 
};

#endif

