#ifndef PDFPOLYNOMIAL
#define PDFPOLYNOMIAL

#include "AbsPdf.h"
#include "Variable.h"
#include "List.h"


#include<algorithm>

template<typename T, int N>
class HornerPoly {
public:
  HornerPoly(){}
  HornerPoly(std::initializer_list<T> il) : p(std::begin(il)+1), c0(*il.begin()){}
  HornerPoly(T const coeff[N+1]) : p(coeff+1), c0(*(coeff)){};
  T operator()(T x) const { return c0 + x*p(x); }
private:
  HornerPoly<T,N-1> p;
  T c0; 
};

template<typename T>
class HornerPoly<T,0> {
public:
  HornerPoly(){}
  HornerPoly(T coeff) : c0(coeff){};
  HornerPoly(T const * coeff) : c0(*coeff){};
  T operator()(T) const { return c0; }
private:
  T c0; 
};


#define RooPolynomial PdfPolynomial


template<int N>
class PdfPolynomial : public AbsPdf {
public:
  using Poly = HornerPoly<double, N>;  // N is the number of coefficient, not the order, still below one adds one...
  
  PdfPolynomial(const Char_t* name, const Char_t* title, Variable &x)  : AbsPdf(name,title,&x), m_x(&x) {   m_nocache=true; }
  
  PdfPolynomial(const Char_t* name, const Char_t* title, Variable &x,
		List<Variable> coeff) :
    AbsPdf(name,title, &x,&coeff), m_x(&x)
  {
    m_coeff.AddElement(coeff);
    assert(m_coeff.GetSize()==N);
  }
  
  virtual ~PdfPolynomial(){}
  
  virtual void GetParameters(List<Variable>& parameters) { parameters.AddElement(m_coeff); }
  


  
  virtual double integral(PdfState const & state) const
  {
    UInt_t order = m_coeff.GetSize();
    Double_t xmaxprod = m_x->GetMax();
    Double_t xminprod = m_x->GetMin();
    Double_t sum = xmaxprod-xminprod;
    for (UInt_t i = 0; i < order; i++) {
      xmaxprod *= m_x->GetMax();
      xminprod *= m_x->GetMin();
      sum += m_coeff.GetElement(i)->value(state)*(xmaxprod - xminprod)/(i+2);
    }
    return sum;
    
  }
  
  
  
  void values(PdfState const & state, double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);

    Data::Value_t const * __restrict__ ldata = data.GetData(*m_x, dataOffset);
    
    auto invInt = invIntegral(state);
    
    Int_t size = m_coeff.GetSize();
    Double_t coeffCPU[size+1];
    loadCoeff(state,coeffCPU,size);
    assert(m_coeff.GetSize()==N);
    
    Poly poly(coeffCPU);
    
#pragma omp simd
    for (auto idx = 0; idx<bsize; ++idx) {
      auto x = ldata[idx];
      auto y = poly(x)*invInt;
      res[idx] = y;
    }

}


 private:

  void loadCoeff(PdfState const & state, Double_t *coeffCPU, UInt_t size) const {
    // the coeffCPU must have the correct dimension (size+1)
    coeffCPU[0] = 1.;
    for(UInt_t i = 0; i<size; i++) {
      coeffCPU[i+1] = m_coeff.GetElement(i)->value(state);
    }
  }

 

 private:
  const Variable *m_x;
  List<Variable> m_coeff;
 
};

#endif

