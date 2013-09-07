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
	 bool dyn): Named(name,title), m_data(&data), m_pdf(&pdf),  
		    minLoop(OpenMP::GetMaxNumThreads(),1000000), 
		    maxLoop(OpenMP::GetMaxNumThreads(),0), 
		    aveLoop(OpenMP::GetMaxNumThreads(),0), 
		    dynamic(dyn) {}

NLL::~NLL() {

  if(dynamic) {
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

}

Double_t NLL::GetVal()
{

  static bool first=true;
  if (first) {
    first=false;
    std::cout << "max threads " << OpenMP::GetMaxNumThreads() << std::endl;
  }
  m_nLoops++;

  m_pdf->CacheIntegral();
  
  m_logs.clear();
  m_logs.resize(OpenMP::GetMaxNumThreads());
  
  if (dynamic) {

    int nloops[OpenMP::GetMaxNumThreads()]={0,};
    std::atomic<int> start(0);
  //int isOk=0;
#pragma omp parallel 
  // reduction(+ : isOk)
    {
      
      // isOk = 
      nloops[omp_get_thread_num()]  = RunEvaluationBlockSplittingDynamic(start);
      
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
      RunEvaluationBlockSplittingStatic();
      
    }

  }

  /*
  std::cout << "tot done " << isOk << std::endl;

  if(omp_in_parallel()) std::cout << "in parallel" << std::endl;
  std::cout << "thread " << OpenMP::GetRankThread() << " of " << OpenMP::GetNumThreads() << std::endl;
  */
  
  //final reduction
  __float128 ss=0.;
  for (unsigned int i=0; i!=OpenMP::GetMaxNumThreads(); ++i)
    ss+=  __float128(-0.693147182464599609375*m_logs[i].value());

  if (m_pdf->IsExtended()) {
    ss += m_pdf->ExtendedTerm(m_data->GetEntries());
  }


  return ss;
}

int NLL::RunEvaluationBlockSplittingStatic() {
  
  int iStart=0, iEnd=0;
  unsigned int ntot = OpenMP::GetThreadElements(m_data->GetEntries(),iStart,iEnd);
  
  static __thread int first(true);
  if (first) 
  {
    first=false;
    Lock l(global::coutLock);
    if(omp_in_parallel()) std::cout << "in parallel" << std::endl;
    
    std::cout << "thread " << OpenMP::GetRankThread() << " of " << OpenMP::GetNumThreads() << std::endl;
    std::cout <<  m_data->GetEntries() << " " << iStart << " " << iEnd << " " << ntot << std::endl;
  }
  

  alignas(ALIGNMENT) double res[m_nBlockEvents];
  auto localValue = m_logs[OpenMP::GetRankThread()];  
  for (UInt_t ie=0; ie<ntot; ie+= m_nBlockEvents) {
    auto offset = iStart+ie;
    auto bsize = std::min(m_nBlockEvents,ntot-ie);
    m_pdf->GetVal(res, bsize, *m_data, offset);  
    PartialNegReduction(localValue,res,bsize);
  }
  m_logs[OpenMP::GetRankThread()] = localValue;

  return 1;
}

int NLL::RunEvaluationBlockSplittingDynamic(std::atomic<int> & start) {

  int ntot = m_data->GetEntries();
  int chunk = 4*m_nBlockEvents;
  int endgame = ntot -  omp_get_num_threads()*chunk;

  int k = omp_get_thread_num();
  alignas(ALIGNMENT) double res[m_nBlockEvents];
  auto localValue = m_logs[k];  

  int lp=0;
  while (true) {
    int ls = start; 
    if (ls>=ntot) break;
    if (ls<endgame) chunk = 2*m_nBlockEvents;

    while (ls<ntot && !std::atomic_compare_exchange_weak(&start,&ls,ls+chunk)); 
    auto ln = std::min(chunk,ntot-ls);
    if (ln<=0) break;
    lp++;
    for (int ie=0; ie<ln; ie+= m_nBlockEvents) {
      auto offset = ls+ie;
      auto bsize = std::min(int(m_nBlockEvents),ln-ie);
      m_pdf->GetVal(res, bsize, *m_data, offset);  
      PartialNegReduction(localValue,res,bsize);
    }

  }
  m_logs[k] = localValue;

  return lp;

}
