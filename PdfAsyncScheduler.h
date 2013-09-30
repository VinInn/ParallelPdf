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


class PdfScheduler;

// can be either integral on chunk, final reduction 
class Task {
public:
  enum What { integral, chunk, reduction};
  
  Task(){}
  Task(PdfScheduler * is, What iw, PdfModifiedState const & ms, unsigned int sid, unsigned int st, unsigned int en=0) :
    what(iw), ls(st), ln(en), stateId(sid), m_state(&ms),  m_scheduler(is){}
  
  void evaluate() const {
    assert(done==false);
    done=true;
    switch (what) {
    case integral:
      doIntegral();
      break;
    case chunk:
      computeChunk();
      break;
    case reduction:
      reduce();
    }
    
  }
  
  
  What what;
  
  // for integral ls is pdf
  // chunk ls, ln are start and end event
  unsigned int ls, ln;

  unsigned int stateId;
  PdfModifiedState const * m_state=nullptr;
  mutable PdfScheduler * m_scheduler=nullptr;

  mutable bool done=false;
  
  PdfState const & state() const { return *m_state;}
  PdfScheduler & scheduler() const { return *m_scheduler; }
  
  void doIntegral() const;
  void computeChunk()  const;
  
  
  void reduce() const {} 
};




class PdfScheduler {
  
public:
  PdfScheduler(size_t inevals, PdfModifiedState const * imstates, double * ires, size_t bsize, CircularBuffer<Task> & itasks) :
    tasks(itasks),
    m_nBlockEvents(bsize), nevals(inevals), mstates(imstates), res(ires),
    integDone(nevals/2), stateReady(nevals/2) {
    for ( auto k=0U; k!=integDone.size(); ++k) integDone[k]=mstates[k].size();
  }
  
  size_t nBlockEvents() const { return m_nBlockEvents;}

  void pushTasks(CircularBuffer<Task> & buff) noexcept;
  
  void integralDone(size_t i) {  
    --integDone[i]; 
  }
  void chunkResult(size_t i, TMath::IntLog value){}

private:
  CircularBuffer<Task> & tasks;
  
  size_t m_nBlockEvents;
  
  Task::What todo=Task::integral;
  
  size_t nevals;
  PdfModifiedState const * mstates;
  double * res;
  
  size_t istate=0;
  size_t idoing=0;
  
  std::vector<std::atomic<int> > integDone;
  std::vector<std::atomic<int> > stateReady;
  
};

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


void Task::doIntegral() const {
  state().cacheIntegral(ls);
  scheduler().integralDone(stateId);    
}

void Task::computeChunk() const {
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



void compute(PdfReferenceState & refState, size_t nevals, PdfModifiedState const * mstates, double * res) {
  
  
  auto busize = std::min(size_t(4*omp_get_max_threads()),nevals);
  CircularBuffer<Task> tasks(busize);
  
  PdfScheduler scheduler(nevals, mstates, res, 512, tasks );
  
  int pushing[omp_get_max_threads()]={0,};
  
  
#pragma omp parallel
  {
    try {
      auto me = omp_get_thread_num();
      Task curr;
      while (true) {
	// first try to push
	if ( (!tasks.draining()) && tasks.halfEmpty() && tasks.tryLock()) {
	  scheduler.pushTasks(tasks);
	  tasks.releaseLock();
	  ++pushing[me];
	}
	int k=busize/4+1;
	while( (tasks.draining() || (--k))  && tasks.pop(curr, false)) {
	  curr.evaluate();
	}
	if (tasks.draining()) break;
      }
      // needed???
      while( tasks.pop(curr)) {
	curr.evaluate();
      }
    } catch(...) {}
  }    
  
  std::cout << "pushers: ";
  for (int i=0; i<omp_get_max_threads(); ++i) std::cout << pushing[i] << ", ";
  std::cout << std::endl;std::cout << std::endl;
  
  
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
void PdfScheduler::pushTasks(CircularBuffer<Task> & buff) noexcept {
  
  // auto allN = omp_get_num_threads();
  // auto meN = omp_get_thread_num();
  
  switch (todo) {
  case Task::integral:
    // integrals
    for (;istate<nevals; ++istate) {
      auto Npdfs = mstates[istate].size();
      for (;idoing<Npdfs; ++idoing) {
	//	{ Guard g(global::coutLock); std::cout << "pushing "<< omp_get_thread_num() << " : integral " << idoing << " " << buff.size() << std::endl; }
	if (!buff.push(Task(this,Task::integral,mstates[istate],istate, idoing),false)) break;
      }
      if (idoing==Npdfs) idoing=0;
      else break;
    }
    if (istate==nevals) { 
      istate=0;
      todo=Task::chunk;
    }
    else break;
  case Task::chunk:
    
    for (;istate<nevals; ++istate) {
    }
    if (istate==nevals) { 
      istate=0;
      todo=Task::reduction;
    }
    else break;
  case Task::reduction:
    //    { Guard g(global::coutLock);  std::cout << "pushing "<< omp_get_thread_num() << " : " << "reduction" << " " << buff.size() << std::endl; }
    buff.drain();
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
