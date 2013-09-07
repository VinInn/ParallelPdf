#ifndef PDFADD3Prod
#define PDFADD3Prod

#include "AbsPdf.h"
#include "List.h"
#include "TMath.h"

#define RooAddPdf PdfAdd


template<typename T, int N>
struct Add3Prod {

  T operator()(T const *  __restrict__ v, T const * __restrict__ c, int i, int stride ) const { return (*c)*(v)[i]*(v+stride)[i]*(v+2*stride)[i] + add(v+3*stride,c+1,i,stride); }
  Add3Prod<T,N-1> add;
  
};

template<typename T>
struct Add3Prod<T,2> {
  
  T operator()(T const * __restrict__ v, T const * __restrict__ c, int i, int stride ) const { return (*c)*(v)[i]*(v+stride)[i]*(v+2*stride)[i] + (*(c+1))*(v+3*stride)[i]*(v+4*stride)[i]*(v+5*stride)[i]; }
};


// very ad hoc...
template<int N>
class PdfAdd3Prod : public AbsPdf {
public:

  PdfAdd3Prod (const Char_t* name, const Char_t* title, List<AbsPdf> pdfs, List<Variable> fractions) :
    AbsPdf(name,title), m_isExtended(kFALSE)
  {
    
    if (pdfs.GetSize()!=3*fractions.GetSize() && pdfs.GetSize()!=3*fractions.GetSize()-1) {
      std::cerr << GetName() << ":: Wrong number of fractions!" << std::endl;
      assert(0);
    }
    
    if (pdfs.GetSize()==3*fractions.GetSize())
      m_isExtended = kTRUE;
    
    m_pdfs.AddElement(pdfs);
    m_fractions.AddElement(fractions);
    
  }
  
  virtual ~PdfAdd3Prod () { }
    
  
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
    alignas(ALIGNMENT) double lres[3*N][strid];
    double coeff[N];

    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    
    Variable *var(0);
    AbsPdf *pdf = iter_pdfs.Next();
    Double_t lastFraction = 1.;
    
    int k=0; int l=0;
    while ((var = iter_fractions.Next())!=0) {
      lastFraction -= var->GetVal();
      coeff[k++]=  var->GetVal();
      for (int j=0; j!=3; ++j) {
	pdf->GetVal(lres[l], bsize, data, dataOffset);
	pdf = iter_pdfs.Next();
	++l;
      }
    }

    if (!m_isExtended) {
      coeff[k]=lastFraction;
      for (int j=0; j!=3; ++j) {
	pdf->GetVal(lres[l], bsize, data, dataOffset);
	++l;
      }
      assert(N==k+1);
   } else 
      assert(N==k);
    assert(3*N==l);    


    Add3Prod<double,N> add;
    auto invIntegral = GetInvIntegral();
    double const * kres = lres[0];
    for (auto idx = 0; idx!=bsize; ++idx) {
      res[idx] = add(kres,coeff,idx,strid)*invIntegral;
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

