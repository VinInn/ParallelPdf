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
class PdfAdd3Prod final : public AbsPdf {
public:

  PdfAdd3Prod (const Char_t* name, const Char_t* title, List<AbsPdf> pdfs, List<Variable> fractions) :
    AbsPdf(name,title, &pdfs, &fractions), 
    m_lmodPdf(omp_get_max_threads(),-2),
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
    makeParameterCache();
    
  }
  
  void makeCache(unsigned int size) {
    m_resCache = std::move(Data("","",size,m_pdfs.GetSize()));
    for (auto i = 0U; i!=m_AllParams.size(); ++i)
      m_parCache[i]=m_AllParams[i]->GetVal();
  } 

  void makeParameterCache() {
    
    // better this loop to be always the same...
    // first me
    for ( auto p : m_fractions())  { 
      m_AllParams.push_back(p);
      m_PdfsPar.push_back(-1);
    }
    int k=0;
    for ( auto pdf : m_pdfs() ) { 
      m_modPdfs.push_back(true);
      m_parPdfs.push_back(m_AllParams.size());
      List<Variable> parameters;
      pdf->GetParameters(parameters);
      for (auto p : parameters()) {
	m_AllParams.push_back(p);
	m_PdfsPar.push_back(k);
      }
      ++k;
    }
    assert(k==m_pdfs().size());
    m_parPdfs.push_back(m_AllParams.size());
    assert(m_parPdfs.size()==m_pdfs().size()+1);
    m_parCache.resize(m_AllParams.size());
  }

  int verifyCache(bool init) {
    for(auto & p : m_lmodPdf) p=-2;
    int n=0;
    if (init) {
      // actually just tell to initiale chaches
      for ( auto k=0; k!=m_pdfs().size(); ++k ) m_modPdfs[k]=true;
      doNotCache = false;
      n = m_AllParams.size();
    } else {
      // assumption: minuit either change one param or all...
      // here we should verify fractions....
      
      // m_parPdfs[0] start of param first pdf
      for ( auto k=0; k!=m_pdfs().size(); ++k ) {
	m_modPdfs[k]=false;
	for (int i = m_parPdfs[k]; i!=m_parPdfs[k+1]; ++i)
	  if (m_AllParams[i]->GetVal()!=m_parCache[i]) { m_modPdfs[k]=true; ++n; break;}
      }
      doNotCache= (n==1);
    }
    if (n>1) {  // resetCache
       // m_parPdfs[0] start of param first pdf
      for (int i = m_parPdfs.front(); i!=m_parPdfs.back(); ++i)
	m_parCache[i]=m_AllParams[i]->GetVal();
    }
    return n;
  }

  virtual ~PdfAdd3Prod () { }
    
  
  virtual void GetParameters(List<Variable>& parameters) 
  {
    parameters()=m_AllParams;
 
  }
  
  virtual void CacheIntegral(int lpar=-2) final {
    doNotCache=true;
    // thread local...
    AbsPdf::CacheIntegral();
    assert(lpar<int(m_PdfsPar.size()));
    auto k = m_PdfsPar[lpar];
    assert(k<int(m_pdfs().size()));
    if(k>=0) m_pdfs()[k]->CacheIntegral(lpar);
    // nasty trick to avoid to add a new interface....
    lmodPdf() = k;
  }

  virtual void CacheAllIntegral() {
    // not thread safe... need to be cached as used for the whole loop
    AbsPdf::CacheAllIntegral();
    for ( auto k=0; k!=m_pdfs().size(); ++k )
      if(m_modPdfs[k]) m_pdfs()[k]->CacheAllIntegral();
  }
  
  unsigned int cacheSize() const final { 
    return m_resCache.capacity();
  }

private:

  virtual Double_t integral() const { return m_isExtended ? ExpectedEvents() : 1.; }
  

  virtual void GetVal(double * __restrict__ res, unsigned int bsize, const Data & data, unsigned int dataOffset) const { 
    res = (double * __restrict__)__builtin_assume_aligned(res,ALIGNMENT);

    auto strid = stride(bsize);
    double * __restrict__ pres[3*N];
    alignas(ALIGNMENT) double lres[3*N][strid];
    double coeff[N];


    List<Variable>::Iterator iter_fractions(m_fractions.GetIterator());
    
    Variable *var(0);
    Double_t lastFraction = 1.;
    
    int k=0;
    while ((var = iter_fractions.Next())!=0) {
      lastFraction -= var->GetVal();
      coeff[k++]=  var->GetVal();
    }
    // this is extended...
    if (!m_isExtended) {
      coeff[k]=lastFraction;
      assert(N==k+1);
    } else 
      assert(N==k);
 
    // form cache
    auto lnocache = doNotCache;
    if (m_resCache.empty()) lnocache=true;
    else
      for (int l=0; l!=15; ++l) { pres[l] = m_resCache.GetData(l,dataOffset);}
 
    auto lp = lmodPdf();
    if (lp>=0) {
      m_pdfs()[lp]->GetVal(lres[lp], bsize, data, dataOffset);
      pres[lp] = &(lres[lp][0]);
      
    } else if(lp!=-1)
      for (int l=0; l!=15; ++l) {
	auto pdf = m_pdfs()[l];
	if (m_modPdfs[l]) { 
	  if (lnocache) {
	    pdf->GetVal(lres[l], bsize, data, dataOffset);
	    pres[l] = &(lres[l][0]);
	  } else {
	    pdf->GetVal(pres[l], bsize, data, dataOffset);
	  }
	}
      }



    Add3Prod<double,N> add;
    auto invIntegral = GetInvIntegral();
    double const * __restrict__  const *  kres = pres;
    // double const * kres = lres[0];
#pragma omp simd
    for (auto idx = 0; idx<bsize; ++idx) {
      // res[idx] = add(kres,coeff,idx,strid)*invIntegral;
      res[idx] = add(kres,coeff,idx)*invIntegral;
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
  
  std::vector<bool> m_modPdfs; // pdf to be called
  std::vector<unsigned short> m_parPdfs; // index in vectors below
  std::vector<Variable*> m_AllParams;
  std::vector<short> m_PdfsPar;  // pdf corresponding to this par... (-1 is this)

  std::vector<double> m_parCache; // cache of param (from  previous call)


  int & lmodPdf() { return m_lmodPdf[omp_get_thread_num()];}
  int  lmodPdf() const { return m_lmodPdf[omp_get_thread_num()];}
  std::vector<int> m_lmodPdf;

  mutable Data m_resCache;

  Bool_t m_isExtended;
  bool doNotCache=true;
 
};

#endif

