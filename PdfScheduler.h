#ifndef PdfScheduler_H
#define PdfScheduler_H

#include "PdfReferenceState.h"
#include "AbsPdf.h"
#include "Variable.h"
#include "List.h"

#include<atomic>


#include<omp.h>

#include "CircularBuffer.h"


#include <thread>
#include <mutex>
typedef std::mutex Mutex;
// typedef std::lock_guard<std::mutex> Lock;
typedef std::unique_lock<std::mutex> Lock;
typedef std::unique_lock<std::mutex> Guard;
//typedef std::condition_variable Condition;


namespace global {
  // control cout....
  Mutex coutLock;
}




class PdfScheduler {
  
public:

  enum What { integral, chunk, reduction};
 

  PdfScheduler(size_t inevals, PdfModifiedState const * imstates, double * ires, size_t bsize) :
    m_nBlockEvents(bsize), nevals(inevals), mstates(imstates), res(ires), istate(0),
    integToDo(nevals/2), integDone(nevals/2), stateReady(nevals/2) {
    for ( auto k=0U; k!=integDone.size(); ++k)  { integToDo[k]=integDone[k]=mstates[k].size(); }
  }
  

  ~PdfScheduler();

  size_t nBlockEvents() const { return m_nBlockEvents;}

  void doTasks() noexcept;
  

  void chunkResult(size_t i, TMath::IntLog value){}

private:
  
  size_t m_nBlockEvents;
  
  What todo=integral;
  
  size_t nevals;
  PdfModifiedState const * mstates;
  double * res;
  
  std::atomic<int> istate;
  
  std::vector<std::atomic<int> > integToDo;
  std::vector<std::atomic<int> > integDone;
  std::vector<std::atomic<int> > stateReady;

  
};


PdfScheduler::~PdfScheduler() {
#ifdef DOPRINT

  for ( auto k=0U; k!=integDone.size(); ++k)
    std::cout << k << ':' << integToDo[k] <<','<< integDone[k] << ' ';
  std::cout << std::endl;


#endif
}


/*
// check depedency, if 0 schedule pdf...
void ready(int i) {
auto k = dep[i];
while (k>0 && !std::atomic_compare_exchange_weak(&dep[i],&k,k-1));
  if (0==k) {
  auto c = qsize;
  while (!std::atomic_compare_exchange_weak(&qsize,&c,c+1));
  queue[c]=i;
  }
  }
*/


/*
void computeChunk() const {
  unsigned int block =  scheduler().nBlockEvents();
  alignas(ALIGNMENT) double lres[block];
    double * res=0;
    TMath::IntLog localValue;
    for (auto ie=0U; ie<ln; ie+= block) {
      auto offset = ls+ie;
      auto bsize = std::min(block,ln-ie);
      res = state().value(lres, bsize,offset);
      assert(res==&lres[0]);
      localValue = IntLogAccumulate(localValue, res, bsize);
    }
    scheduler().chunkResult(stateId,localValue);
}
*/


void compute(PdfReferenceState & refState, size_t nevals, PdfModifiedState const * mstates, double * res) {
  
  
  PdfScheduler scheduler(nevals, mstates, res, 512);
  
  
  
#pragma omp parallel
  {
    try {
      scheduler.doTasks();
      
    } catch(...) {}
  }    
  
  
  
}





void differentiate(PdfReferenceState & refState, unsigned int nvar, int const * vars, double const * steps, double *  res) {
  auto nevals = 2*nvar;
  
  PdfModifiedState mstate[nevals];
  for (auto il=0U; il!=nevals; ++il) {
    auto ik = il/2; // hope optmize in >1
    auto v = refState.paramVal(vars[ik]);
    auto nv = v + steps[il];
    mstate[il] = PdfModifiedState(&refState,vars[ik], nv);
  }
  double d[nevals];
  compute(refState,nevals, mstate, d);
  for (auto il=0U; il!=nevals; il+=2) {
    auto ik = il/2; // hope optmize in >1
    res[ik] = (d[il+1]-d[il])/(steps[il+1]-steps[il]);
  }
}


// inside a single thred...
void PdfScheduler::doTasks() noexcept {
  
  // auto allN = omp_get_num_threads();
  // auto meN = omp_get_thread_num();
 
  int ls = istate;
  switch (todo) {
   case integral:
    // integrals
    while(true) {
      if (ls==int(nevals)) break;
      auto & aw = integToDo[ls]; 
      while(true) {
	int is =  aw;
	while (is>0 && !std::atomic_compare_exchange_weak(&aw,&is,is-1));
	if (is<=0) break;
	mstates[ls].cacheIntegral(is-1);
	--integDone[ls];
      }
      ls = istate;
      while (ls<int(nevals) && !std::atomic_compare_exchange_weak(&istate,&ls,ls+1));
      if (ls==int(nevals)) break;
    }
    todo=chunk;
  case chunk:
    if (istate==int(nevals)) { 
      todo=reduction;
    }
  case reduction:
    //    { Guard g(global::coutLock);  std::cout << "pushing "<< omp_get_thread_num() << " : " << "reduction" << " " << buff.size() << std::endl; }
    break; 
  }
}

/*
  // ok now events chunks...
  
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
    computeChunk(..);

  }
  m_logs[k] = localValue;

    */



  /*
  // cache what needed
  while (true) {
    unsigned int ls = curr; 
    if (ls>=qsize) break;
    while (ls<qsize && !std::atomic_compare_exchange_weak(&curr,&ls,ls+1));
    if (ls>=qsize) break;
    auto ipdf = queue[ls];
    state().cachePdf(ipdf,bsize,data,offset);
  */





#endif
