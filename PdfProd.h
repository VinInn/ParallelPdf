#ifndef PDFPROD
#define PDFPROD

#include "AbsPdf.h"
#include "List.h"
#include "openmp.h"

template<typename T, int N>
struct Mult {
  T operator()(T const *  __restrict__ const * v,  int i ) const { return (*v)[i] * mult(v+1,i); }
  T operator()(T const * v,  int i,  int stride ) const { return (v)[i] * mult(v+stride,i,stride); }
  T operator()(std::initializer_list<T> il) const { return *il.begin() * mult(std::begin(il)+1);}
  Mult<T,N-1> mult;
  
};

template<typename T>
struct Mult<T,2> {
  
  T operator()(T x, T y) const { return x*y;}
  T operator()(T const *  __restrict__ const * v,  int i ) const { return (*v)[i] * (*(v+1))[i]; }
  T operator()(T const * v, int i,  int stride ) const { return (v)[i] * (v+stride)[i]; }
  T operator()(std::initializer_list<T> il) const { return *il.begin() * *(std::begin(il)+1);}
};



#define RooProdPdf PdfProd

template<int N>
class PdfProd : public NoCacheAbsPdf {
public:
  PdfProd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2) :
    NoCacheAbsPdf(name,title,&pdf1,&pdf2)
  {
    m_pdfs.AddElement(pdf1);
    m_pdfs.AddElement(pdf2);
    
  }
  PdfProd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs):
    NoCacheAbsPdf(name,title,&pdfs)
  {
    m_pdfs.AddElement(pdfs);
  }

  virtual ~PdfProd() { }
  
    
  virtual void GetParameters(List<Variable>& parameters){
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->GetParameters(parameters);
    }
  }
  
  
private:  

  virtual double integral(PdfState const & state) const { return  1.; }
  
  
  virtual void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);
    
    auto strid = stride(bsize);
    double * __restrict__ pres[N];
    alignas(ALIGNMENT) double lres[N][strid];

    for (int l=0; l!=N; ++l) {
      auto pdf = m_pdfs()[l];
      (*pdf)(state,pres[l], &(lres[l][0]), bsize, data, dataOffset);
    }

    
    Mult<double,N> mult;
    double const * __restrict__  const *  kres = pres;
#pragma omp simd
    for (auto idx = 0; idx<bsize; ++idx) {
      res[idx] = mult(kres,idx);
    }
  }
  
private:
  
  mutable List<AbsPdf> m_pdfs;
  
  
};

#endif

