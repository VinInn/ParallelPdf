#ifndef PDFADD3Prod
#define PDFADD3Prod

#include "AbsPdf.h"
#include "List.h"
#include "TMath.h"

#define RooAddPdf PdfAdd


template<typename T, int N>
struct Add3Prod {
  T operator()(T const *  __restrict__ const * v, T const * __restrict__ c, int i ) const { return (*c)*(*v)[i]*(*(v+1))[i]*(*(v+2))[i] + add(v+3,c+1,i); }
  T operator()(T const *  __restrict__ v, T const * __restrict__ c, int i, int stride ) const { return (*c)*(v)[i]*(v+stride)[i]*(v+2*stride)[i] + add(v+3*stride,c+1,i,stride); }
  Add3Prod<T,N-1> add;
  
};

template<typename T>
struct Add3Prod<T,2> {
  T operator()(T const *  __restrict__ const * v, T const * __restrict__ c, int i ) const { return (*c)*(*v)[i]*(*(v+1))[i]*(*(v+2))[i] + (*(c+1))*(*(v+3))[i]*(*(v+4))[i]*(*(v+5))[i]; }
  T operator()(T const * __restrict__ v, T const * __restrict__ c, int i, int stride ) const { return (*c)*(v)[i]*(v+stride)[i]*(v+2*stride)[i] + (*(c+1))*(v+3*stride)[i]*(v+4*stride)[i]*(v+5*stride)[i]; }
};


// very ad hoc...
template<int N>
class PdfAdd3Prod final : public NoCacheAbsPdf {
public:

  PdfAdd3Prod (const Char_t* name, const Char_t* title, List<AbsPdf> pdfs, List<Variable> fractions) :
    NoCacheAbsPdf(name,title, &pdfs, &fractions),
    m_isExtended(kFALSE)
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

  

  virtual double integral(PdfState const & state) const { return m_isExtended ? ExpectedEvents(state) : 1.; }
  

  virtual void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);

    auto strid = stride(bsize);
    double * __restrict__ pres[3*N];
    alignas(ALIGNMENT) double lres[3*N][strid];  // to store not cached one....
    double coeff[N];


    Double_t lastFraction = 1.;
    
    int k=0;
    for (auto var : m_fractions()) {
      lastFraction -= var->value(state);
      coeff[k++]=  var->value(state);
    }

    // this is extended...
    if (!m_isExtended) {
      coeff[k]=lastFraction;
      assert(N==k+1);
    } else 
      assert(N==k);
 
    for (int l=0; l!=3*N; ++l) {
      auto pdf = m_pdfs()[l];
      pres[l] = (*pdf)(state, &(lres[l][0]), bsize, dataOffset);
    }



    Add3Prod<double,N> add;
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
      for (auto var : m_fractions())
	nEvents += var->value(state);
    }
    
  return nEvents;
}
  
private:
  
  mutable List<AbsPdf> m_pdfs;
  mutable List<Variable> m_fractions;
  
 
  Bool_t m_isExtended;
 
};

#endif

