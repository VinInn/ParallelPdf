#ifndef PDFADD
#define PDFADD

#include "AbsPdf.h"
#include "List.h"
#include "TMath.h"
#include "Results.h"
#include "openmp.h"

#define RooAddPdf PdfAdd


template<typename T, int N>
struct Add {

  T operator()(T const * const * v, T const * c, int i ) const { return (*c)*(*v)[i] + add(v+1,c+1,i); }
  T operator()(std::initializer_list<T> il) const { return *il.begin()+add(std::begin(il)+1);}
  Add<T,N-1> add;
  
};

template<typename T>
struct Add<T,2> {
  
  T operator()(T x, T y) const { return x+y;}
  T operator()(T const * const * v, T const * c, int i ) const { return (*c)*(*v)[i] + (*(c+1))*(*(v+1))[i]; }
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
  
  virtual const Results & GetValSIMD(UInt_t iStart = 0, UInt_t nPartialEvents = 0) {
    assert(m_doCalculationBy!=kVirtual && 
	   evaluateSIMD(iStart,nPartialEvents,1./GetIntegral()));
    
    return m_resultsCPU;
  }
  
  
  
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
  
  virtual void Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy = kOpenMP)
  {
    AbsPdf::Init(data,doLog,doCalculationBy);
    AbsPdf *pdf(0);
    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    while ((pdf = iter_pdfs.Next())!=0) {
      pdf->Init(data,doLog,doCalculationBy);
    }
    
  }
  
  virtual void ClearResults(Bool_t recursive = kFALSE) {
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
    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    
    AbsPdf *pdf = iter_pdfs.Next();
    Variable *var(0);
    Double_t lastFraction = 1.;
    Double_t ret(0);
    
    while ((var = iter_fractions.Next())!=0) {
      lastFraction -= var->GetVal();
      ret += var->GetVal()*pdf->GetVal();
      pdf = iter_pdfs.Next();
    }
    
    if (!m_isExtended)
      ret += lastFraction*pdf->GetVal();
    
    return ret;
  }

  virtual Double_t integral() const { return m_isExtended ? ExpectedEvents() : 1.; }
  
  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents, 
			      const Double_t invIntegral) {

    auto nEvents = m_data->GetEntries();
    auto nEventsThread = GetParallelNumElements(nEvents);

    const double * res[N];
    double coeff[N];

    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    
    Variable *var(0);
    AbsPdf *pdf = iter_pdfs.Next();
    Double_t lastFraction = 1.;
 
    UInt_t nResults(0), iPartialEnd;
    
    int k=0;
    while ((var = iter_fractions.Next())!=0) {
      lastFraction -= var->GetVal();
      coeff[k]=  var->GetVal();
      const Results* pResultsPdf = &(pdf->GetValSIMD(iPartialStart,nPartialEvents));
      res[k] = pResultsPdf->GetData(nResults,iPartialEnd,
				    iPartialStart,nPartialEvents);
      pdf = iter_pdfs.Next();
      ++k;
    }

    if (!m_isExtended) {
      coeff[k]=lastFraction;
      const Results* pResultsPdf = &(pdf->GetValSIMD(iPartialStart,nPartialEvents));
      res[k] = pResultsPdf->GetData(nResults,iPartialEnd,
				  iPartialStart,nPartialEvents);
      assert(N==k+1);
   } else 
      assert(N==k);
    
    double * resultsCPU = m_resultsCPU.AllocateData(nEventsThread);

    Add<double,N> add;
    const Int_t iE = iPartialEnd;
    for (auto idx = (Int_t)iPartialStart; idx<iE; idx++) {
      resultsCPU[idx] = add(res,coeff,idx)*invIntegral;
    }
    
    return true;

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

