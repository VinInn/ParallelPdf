#ifndef PDFPROD
#define PDFPROD

#include "AbsPdf.h"
#include "List.h"
#include "openmp.h"

template<typename T, int N>
struct Mult {

  T operator()(T const * v,  int i,  int stride ) const { return (v)[i] * mult(v+stride,i,stride); }
  T operator()(std::initializer_list<T> il) const { return *il.begin() * mult(std::begin(il)+1);}
  Mult<T,N-1> mult;
  
};

template<typename T>
struct Mult<T,2> {
  
  T operator()(T x, T y) const { return x*y;}
  T operator()(T const * v, int i,  int stride ) const { return (v)[i] * (v+stride)[i]; }
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
  
    
  virtual void CacheIntegral(){
    AbsPdf::CacheIntegral();
    for (auto pdf : m_pdfs()) pdf->CacheIntegral();
  }
 

  virtual void CacheAllIntegral(){
    AbsPdf::CacheAllIntegral();
    for (auto pdf : m_pdfs()) pdf->CacheAllIntegral();
  }     
    
  
  
  virtual void GetParameters(List<Variable>& parameters){
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->GetParameters(parameters);
    }
  }
  
  
private:  

  virtual Double_t integral() const { return 1.; }
  
  virtual void GetVal(double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 

    auto strid = stride(bsize);
    alignas(ALIGNMENT) double lres[N][strid];

    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    
    int k=0;
    AbsPdf *pdf(0);
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->GetVal(lres[k++], bsize, data, dataOffset);
    }
    assert(k==N);
    
    
    Mult<double,N> mult;
    double const *  kres = lres[0];
    for (auto idx = 0; idx!=bsize; ++idx) {
      res[idx] = mult(kres,idx, strid);
    }
  }
  
private:
  
  mutable List<AbsPdf> m_pdfs;
  
  
};

#endif

