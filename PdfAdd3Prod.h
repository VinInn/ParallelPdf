#ifndef PDFADD3Prod
#define PDFADD3Prod

#include "AbsPdf.h"
#include "List.h"
#include "TMath.h"
#include "Results.h"
#include "openmp.h"

#define RooAddPdf PdfAdd


template<typename T, int N>
struct Add3Prod {

  T operator()(T const * const * __restrict__ v, T const * __restrict__ c, int i ) const { return (*c)*(*v)[i]*(*(v+1))[i]*(*(v+2))[i] + add(v+3,c+1,i); }
  Add3Prod<T,N-1> add;
  
};

template<typename T>
struct Add3Prod<T,2> {
  
  T operator()(T const * const * __restrict__ v, T const * __restrict__ c, int i ) const { return (*c)*(*v)[i]*(*(v+1))[i]*(*(v+2))[i] + (*(c+1))*(*(v+3))[i]*(*(v+4))[i]*(*(v+5))[i]; }
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
      ret += var->GetVal()*pdf->GetVal();
      pdf = iter_pdfs.Next();
      ret += var->GetVal()*pdf->GetVal();
      pdf = iter_pdfs.Next();
    }
    
    if (!m_isExtended){
      ret += lastFraction*pdf->GetVal();
      pdf = iter_pdfs.Next();
      ret += var->GetVal()*pdf->GetVal();
      pdf = iter_pdfs.Next();
      ret += var->GetVal()*pdf->GetVal();
    }
    return ret;
  }

  virtual Double_t integral() const { return m_isExtended ? ExpectedEvents() : 1.; }
  
  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents, 
			      const Double_t invIntegral) {

    auto nEvents = m_data->GetEntries();
    auto nEventsThread = GetParallelNumElements(nEvents);

    const double * res[3*N];
    double coeff[N];

    List<AbsPdf>::Iterator iter_pdfs(m_pdfs.GetIterator());
    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    
    Variable *var(0);
    AbsPdf *pdf = iter_pdfs.Next();
    Double_t lastFraction = 1.;
 
    UInt_t nResults(0), iPartialEnd;
    
    int k=0; int l=0;
    while ((var = iter_fractions.Next())!=0) {
      lastFraction -= var->GetVal();
      coeff[k++]=  var->GetVal();
      for (int j=0; j!=3; ++j) {
	const Results* pResultsPdf = &(pdf->GetValSIMD(iPartialStart,nPartialEvents));
	res[l] = pResultsPdf->GetData(nResults,iPartialEnd,
				    iPartialStart,nPartialEvents);
	pdf = iter_pdfs.Next();
	++l;
      }
    }

    if (!m_isExtended) {
      coeff[k]=lastFraction;
      for (int j=0; j!=3; ++j) {
	const Results* pResultsPdf = &(pdf->GetValSIMD(iPartialStart,nPartialEvents));
	res[l] = pResultsPdf->GetData(nResults,iPartialEnd,
				      iPartialStart,nPartialEvents);
	++l;
      }
      assert(N==k+1);
   } else 
      assert(N==k);
    assert(3*N==l);    

    double * __restrict__ resultsCPU = m_resultsCPU.AllocateData(nEventsThread);
    resultsCPU = (double *)__builtin_assume_aligned (resultsCPU, 32, 0);

    Add3Prod<double,N> add;
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

