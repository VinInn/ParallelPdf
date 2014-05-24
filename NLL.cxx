#include "NLL.h"
#include "openmp.h"
#include <mutex>
#include <atomic>
typedef std::mutex Mutex;
// typedef std::lock_guard<std::mutex> Lock;
typedef std::unique_lock<std::mutex> Lock;

namespace global {
  // control cout....
  Mutex coutLock;
}

#include <iostream>

NLL::NLL(const Char_t* name, const Char_t* title, Data &data, AbsPdf &pdf,
	 int dyn, bool idocache): Named(name,title), m_data(&data), m_pdf(&pdf),  
				   minLoop(OpenMP::GetMaxNumThreads(),1000000), 
				   maxLoop(OpenMP::GetMaxNumThreads(),0), 
				   aveLoop(OpenMP::GetMaxNumThreads(),0), 
				   m_ngroups(dyn), docache(idocache) {}

NLL::~NLL() {

  if(m_ngroups>0) {
    std::cout << "min dyn sched "; 
    for (auto l : minLoop) 
      std::cout << l << " ";
    std::cout << std::endl;
    std::cout << "max dyn sched "; 
    for (auto l : maxLoop) 
      std::cout << l << " ";
    std::cout << std::endl;
     std::cout << "ave dyn sched "; 
    for (auto l : aveLoop) 
      std::cout << double(l)/double(m_nLoops) << " ";
    std::cout << std::endl;
  }

  // std::cout << "\nCache size " << m_pdf->cacheSize() << std::endl;

}




double  NLL::GetVal(PdfState& state) {
  // assume to be in a thread..


  auto ntot = m_data->GetEntries();
  auto par = 0U;
  auto start =0U;
  auto size   =  ntot;

  // a hack
  if ( OpenMP::IsInParallel()) {
    par = Data::partition();
    start =  m_data->startP();
    size   =  m_data->sizeP();
  }

  std::vector<unsigned short> pdfs; std::vector<unsigned short>  dep;
  state.deps(pdfs,dep, false);

  for (auto i: pdfs) state.cacheIntegral(i);


  alignas(ALIGNMENT) double lres[m_nBlockEvents];
  double * res=0;
  TMath::IntLog localValue;
  for (auto ie=start; ie<start+size; ie+= m_nBlockEvents) {
    auto offset = ie;
    auto bsize = std::min(m_nBlockEvents,(start+size)-ie);
    res = (*m_pdf)(state, lres, bsize, offset);
    assert(res==&lres[0]);
    PartialNegReduction(localValue,res,bsize);
  }
  auto ret = -0.693147182464599609375*localValue.value();

  if (par==0 && m_pdf->IsExtended())
    ret += m_pdf->ExtendedTerm(state,ntot);

  return ret;
}


Double_t NLL::GetVal(bool verify)
{

  std::vector<unsigned short> pdfs; std::vector<unsigned short>  dep;

  PdfReferenceState& refState = PdfReferenceState::me();
  PdfNoCacheState ncState(&refState);

  PdfState * pState = docache ? (PdfState *)(&refState) : (PdfState *)(&ncState);
  PdfState& state = *pState;

  if (verify) state.deps(pdfs,dep, docache);
  else state.allDeps(pdfs,dep, docache);

  static bool first=true;
  if (first) {
    first=false;
    std::cout << "max threads " << OpenMP::GetMaxNumThreads() << std::endl;
  }
  m_nLoops++;

  for (auto i: pdfs) state.cacheIntegral(i);

  
  m_logs.clear();
  m_logs.resize(OpenMP::GetMaxNumThreads());

  
  if (m_ngroups>0) {

    int nloops[OpenMP::GetMaxNumThreads()];
    for (int i=0; i<OpenMP::GetMaxNumThreads(); ++i) nloops[i]=0;

    std::atomic<int> istart[m_ngroups];
    int iend[m_ngroups];
    for (int ig=0; ig!=m_ngroups; ++ig) {
      int ls=0;
      Partitioner::GetElements(m_ngroups,ig,m_data->GetEntries(),ls,iend[ig]);
      istart[ig]=ls;
    }


  //int isOk=0;
#pragma omp parallel 
  // reduction(+ : isOk)
    {
      
      // isOk = 
      nloops[omp_get_thread_num()]  = RunEvaluationBlockSplittingDynamic(state, pdfs, istart, iend);
      
    }
    int k=0;
    for (auto l : nloops) {
      minLoop[k]=std::min(minLoop[k],l);
      maxLoop[k]=std::max(maxLoop[k],l);
      aveLoop[k]+=l;
      k++;
    }
    /*
    std::cout << "dyn sched "; 
    for (auto l : nloops) 
      std::cout << l << " ";
    std::cout << std::endl;
    */
  } else {
    
    
    //int isOk=0;
#pragma omp parallel 
    // reduction(+ : isOk)
    {
      
      // isOk = 
      RunEvaluationBlockSplittingStatic(state, pdfs);
      
    }

  }

  /*
  std::cout << "tot done " << isOk << std::endl;

  if(omp_in_parallel()) std::cout << "in parallel" << std::endl;
  std::cout << "thread " << OpenMP::GetRankThread() << " of " << OpenMP::GetNumThreads() << std::endl;
  */
  
  //final reduction
  __float128 ss=0.;
  for (int i=0; i!=OpenMP::GetMaxNumThreads(); ++i)
    ss+=  __float128(-0.693147182464599609375*m_logs[i].value());

  if (m_pdf->IsExtended()) {
    ss += m_pdf->ExtendedTerm(state,m_data->GetEntries());
  }


  return ss;
}

int NLL::RunEvaluationBlockSplittingStatic(PdfState const & state, std::vector<unsigned short> const & pdfs) {
  
  int iStart=0, iEnd=0;
  unsigned int ntot = OpenMP::GetThreadElements(m_data->GetEntries(),iStart,iEnd);
 
  /*
  static __thread int first(true);
  if (first) 
  {
    first=false;
    Lock l(global::coutLock);
    if(omp_in_parallel()) std::cout << "in parallel" << std::endl;
    
    std::cout << "thread " << OpenMP::GetRankThread() << " of " << OpenMP::GetNumThreads() << std::endl;
    std::cout <<  m_data->GetEntries() << " " << iStart << " " << iEnd << " " << ntot << std::endl;
  }
  */

  alignas(ALIGNMENT) double lres[m_nBlockEvents];
  double * res = lres;
  auto localValue = m_logs[OpenMP::GetRankThread()];  
  for (UInt_t ie=0; ie<ntot; ie+= m_nBlockEvents) {
    auto offset = iStart+ie;
    auto bsize = std::min(m_nBlockEvents,ntot-ie);
    for (auto i: pdfs) state.cachePdf(i,bsize, offset);
    res = (*m_pdf)(state, lres, bsize, offset);
    assert(res==&lres[0]);
    PartialNegReduction(localValue,res,bsize);
  }
  m_logs[OpenMP::GetRankThread()] = localValue;

  return 1;
}

int NLL::RunEvaluationBlockSplittingDynamic(PdfState const & state, std::vector<unsigned short> const & pdfs, std::atomic<int> * istart, int const * iend) {

  int chunk = 4*m_nBlockEvents;
  int endgame = omp_get_num_threads()*chunk;

  int k = omp_get_thread_num();
  int ig =  k/(omp_get_num_threads()/m_ngroups);
  assert(ig<m_ngroups);
  std::atomic<int> & start = istart[ig];
  auto end = iend[ig];

  alignas(ALIGNMENT) double lres[m_nBlockEvents];
  double * res=0;
  auto localValue = m_logs[k];  

  int lp=0;
  while (true) {
    int ls = start; 
    if (ls>=end) break;
    if ( (end-ls)<endgame) chunk = 2*m_nBlockEvents;

    while (ls<end && !std::atomic_compare_exchange_weak(&start,&ls,ls+chunk)); 
    auto ln = std::min(chunk,end-ls);
    if (ln<=0) break;
    lp++;
    for (int ie=0; ie<ln; ie+= m_nBlockEvents) {
      auto offset = ls+ie;
      auto bsize = std::min(int(m_nBlockEvents),ln-ie);
      for (auto i: pdfs) state.cachePdf(i,bsize,offset);
      res = (*m_pdf)(state, lres, bsize, offset);
      assert(res==&lres[0]);
      PartialNegReduction(localValue,lres,bsize);
    }

  }
  m_logs[k] = localValue;

  return lp;

}
