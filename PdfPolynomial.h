#ifndef PDFPOLYNOMIAL
#define PDFPOLYNOMIAL

#include "AbsPdf.h"
#include "Variable.h"
#include "List.h"

#include "tbb.h"


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
  
  PdfPolynomial(const Char_t* name, const Char_t* title, const Variable &x)  : AbsPdf(name,title), m_x(&x) {}
  
  PdfPolynomial(const Char_t* name, const Char_t* title, const Variable &x,
		List<Variable> coeff) :
    AbsPdf(name,title), m_x(&x)
  {
    m_coeff.AddElement(coeff);
    assert(m_coeff.GetSize()==N);
  }
  
  virtual ~PdfPolynomial(){}
  
  virtual void GetParameters(List<Variable>& parameters) { parameters.AddElement(m_coeff); }
  static void SetBlockSize(Int_t blockSize) {  }
  
protected:
  virtual Double_t evaluate() const
  {
    UInt_t size = m_coeff.GetSize();
    Double_t coeffCPU[size+1];
    loadCoeff(coeffCPU,size);
    return evaluateLocal(m_x->GetVal(),coeffCPU,size);
  }
  
  virtual Double_t integral() const
  {
    UInt_t order = m_coeff.GetSize();
    Double_t xmaxprod = m_x->GetMax();
    Double_t xminprod = m_x->GetMin();
    Double_t sum = xmaxprod-xminprod;
    for (UInt_t i = 0; i < order; i++) {
      xmaxprod *= m_x->GetMax();
      xminprod *= m_x->GetMin();
      sum += m_coeff.GetElement(i)->GetVal()*(xmaxprod - xminprod)/(i+2);
    }
    return sum;
    
  }
  
  
  
  
  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			      const Double_t invIntegral) {
    
    const Data::Value_t  *__restrict__ dataCPU = m_data->GetCPUData(*m_x);
    dataCPU = (const Data::Value_t *)__builtin_assume_aligned (dataCPU, 32, 0);

    if (dataCPU==0)
      return kFALSE;
    
    Int_t size = m_coeff.GetSize();
    Double_t coeffCPU[size+1];
    loadCoeff(coeffCPU,size);
    assert(m_coeff.GetSize()==N);
    
   Poly poly(coeffCPU);
  
   UInt_t iPartialEnd(0);
   Double_t* __restrict__ resultsCPU = GetDataResultsCPUThread(dataCPU,iPartialEnd,iPartialStart,nPartialEvents);
   resultsCPU = (Double_t *)__builtin_assume_aligned (resultsCPU, 32, 0);
  

  int is = iPartialStart; int ie = iPartialEnd;
  double iI = invIntegral; 
  for (auto idx = is; idx<ie; ++idx) {
    auto x = dataCPU[idx];
    auto y = poly(x)*iI;
    resultsCPU[idx] = y;

  }
  
  
  return kTRUE;
}

  
  


 private:

  inline const Double_t *loadCoeff(Double_t *coeffCPU, UInt_t size) const {
    // the coeffCPU must have the correct dimension (size+1)
    coeffCPU[0] = 1.;
    for(UInt_t i = 0; i<size; i++) {
      coeffCPU[i+1] = m_coeff.GetElement(i)->GetVal();
    }
    
    return &coeffCPU[0];

  }

  inline Double_t evaluateLocal(const Double_t x, const Double_t *coeff, UInt_t order) const {

    double result = coeff[order];
    for (;order>0;--order)
      result = result*x+coeff[order-1];

    return result;

  }

  inline Double_t evaluateLocalSingleCoeff(const Double_t x, const Double_t coeff, const Double_t result) const {
    return (result*x+coeff);
  }

 private:
  const Variable *m_x;
  List<Variable> m_coeff;
  // static Int_t m_blockSize;
};

#endif

