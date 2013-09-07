#ifndef PDFPROD
#define PDFPROD

#include "AbsPdf.h"
#include "List.h"
#include "openmp.h"

template<typename T, int N>
struct Mult {

  T operator()(T const * const * v,  int i ) const { return (*v)[i] * mult(v+1,i); }
  T operator()(std::initializer_list<T> il) const { return *il.begin() * mult(std::begin(il)+1);}
  Mult<T,N-1> mult;
  
};

template<typename T>
struct Mult<T,2> {
  
  T operator()(T x, T y) const { return x*y;}
  T operator()(T const * const * v, int i ) const { return (*v)[i] * (*(v+1))[i]; }
  T operator()(std::initializer_list<T> il) const { return *il.begin() * *(std::begin(il)+1);}
};



#define RooProdPdf PdfProd

template<int N>
class PdfProd : public AbsPdf {
public:
  PdfProd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2) :
    AbsPdf(name,title)
  {
    m_pdfs.AddElement(pdf1);
    m_pdfs.AddElement(pdf2);
    
  }
  PdfProd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs):
    AbsPdf(name,title)
  {
    m_pdfs.AddElement(pdfs);
  }

  virtual ~PdfProd() { }
  
  virtual const Results &GetValSIMD(UInt_t iStart = 0, UInt_t nPartialEvents = 0){
    // No normalization is required
    assert(m_doCalculationBy!=kVirtual && 
	   evaluateSIMD(iStart,nPartialEvents,1.));
    
    return m_resultsCPU;
  }
  
  
  
  virtual void CacheIntegral(){
    AbsPdf::CacheIntegral();
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->CacheIntegral();
    }
    
  }
  
  
  virtual void GetParameters(List<Variable>& parameters){
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->GetParameters(parameters);
    }
  }
  
  
  virtual void Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy = kOpenMP){
    AbsPdf::Init(data,doLog,doCalculationBy);
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->Init(data,doLog,doCalculationBy);
    }
    
  }
  
  
  virtual void ClearResults(Bool_t recursive = kFALSE){
    AbsPdf::ClearResults();
    if (recursive) {
      AbsPdf *pdf(0);
      List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
      while ((pdf = iter_pdfs.Next())!=0) {
	pdf->ClearResults(recursive);
      }
    }
  }
  
protected:

  virtual Double_t evaluate() const {
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    AbsPdf *pdf = iter_pdfs.Next();
    Double_t ret = pdf->GetVal();
    
    while ((pdf = iter_pdfs.Next())!=0) {
      ret *= pdf->GetVal();
    }
    
    return ret;
  }
  
  virtual Double_t integral() const { return 1.; }
  
  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			      const Double_t /* invIntegral */) {
    
    auto nEvents = m_data->GetEntries();
    auto nEventsThread = GetParallelNumElements(nEvents);
    
    UInt_t nResults(0), iPartialEnd;
    
    double const * res[N];
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    
    int k=0;
    AbsPdf *pdf(0);
    while ((pdf = iter_pdfs.Next())!=0) {
      const Results* pResultsPdf = &(pdf->GetValSIMD(iPartialStart,nPartialEvents));
      const Double_t* __restrict__ resultsPdf = pResultsPdf->GetData(nResults,iPartialEnd,
								     iPartialStart,nPartialEvents);
      res[k++] = (Double_t const *)__builtin_assume_aligned (resultsPdf, 32, 0);
      
    }
    assert(k==N);
    
    double * resultsCPU = m_resultsCPU.AllocateData(nEventsThread);
    
    Mult<double,N> mult;
    const Int_t iE = iPartialEnd;
    for (auto idx = (Int_t)iPartialStart; idx<iE; idx++) {
      resultsCPU[idx] = mult(res,idx);
    }
    
    return true;
    
  }
  
private:
  
  mutable List<AbsPdf> m_pdfs;
  
  
};

#endif

