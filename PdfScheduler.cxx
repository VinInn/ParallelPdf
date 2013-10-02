#include "PdfScheduler.h"


namespace global {
  // control cout....
  extern Mutex coutLock;
}

void PdfScheduler::chunkResult(size_t i, TMath::IntLog value) {
  Guard g(global::coutLock); // FIXME
  values[i].reduce(value);
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



void PdfScheduler::computeChunk(unsigned int ist, unsigned int icu) {
  // { Guard g(global::coutLock);  std::cout << nPdfToEval << " chunk "<< omp_get_thread_num() << " : " <<  ist << " " << icu << std::endl; }

  // stupid spinlock waiting for integrals...
  while (ist>=nPdfToEval) std::this_thread::yield();
  auto k = pdfToEval[ist];

  auto block =  nBlockEvents();
  auto chunk = 4*block;  // we know
  auto const & data = mstates[k].data();

  auto ls = data.startP()+icu*chunk;
  auto ln = std::min(data.sizeP()-icu*chunk,chunk);

  alignas(ALIGNMENT) double lres[block];
  double * res=0;
  TMath::IntLog localValue;
  for (auto ie=0U; ie<ln; ie+= block) {
    auto offset = ls+ie;
    auto bsize = std::min(block,ln-ie);
    res = mstates[k].value(lres, bsize,offset);
    assert(res==&lres[0]);
    localValue = IntLogAccumulate(localValue, res, bsize);
  }
  chunkResult(k,localValue);
  
  /*
  { Guard g(global::coutLock);  
    std::cout << nPdfToEval << " chunk "<< omp_get_thread_num() << " " 
	      << mstates[k].param() <<',' << mstates[k].paramVal(mstates[k].param())
	      << " " << ls<<',' <<ln
	      << " : " <<  ist << ' ' << k << " " << icu << ' ' << localValue.value() << std::endl; 
	      }
  */
}



namespace {
  void compute(PdfReferenceState & refState, size_t nevals, PdfModifiedState const * mstates, double * res) {
    
    TMath::IntLog  values[nevals];
    
    PdfScheduler scheduler(nevals, mstates, values, 512);
  
    
    
#pragma omp parallel
    {
      try {
	scheduler.doTasks();
	
      } catch(...) {}
    }    
    
    auto m_pdf = refState.pdfs().back();
    auto ntot = refState.data().size();

    for (auto k=0U; k<nevals; ++k) {
      res[k] = -0.693147182464599609375*values[k].value();
      if (m_pdf->IsExtended())
      res[k] += m_pdf->ExtendedTerm(mstates[k],ntot);
    }
    
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
  
  //auto const meT = omp_get_thread_num();

  auto const meG = Data::partition();



  int ls = istate;
  
  int lc = nChunks[meG];
  int start=iChunks[meG];

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
	mstates[ls].cacheYourIntegral(is-1);

	// pseudo queue
	is = integDone[ls];
	assert(is>0);
	while (!std::atomic_compare_exchange_weak(&integDone[ls],&is,is-1));
	assert(is>0);
	if (1==is) {
	  //push ls
	  unsigned int k = nPdfToEval;
	  while (!std::atomic_compare_exchange_weak(&nPdfToEval,&k,k+1));
	  assert(k<nevals);
	  pdfToEval[k]=ls;
	} 
      }
      ls = istate;
      while (ls<int(nevals) && !std::atomic_compare_exchange_weak(&istate,&ls,ls+1));
      if (ls==int(nevals)) break;
    }
    todo=chunk;
    // break;
  case chunk:
    while(true) {
      if (lc<=0) break;
      auto k = start+lc-1;
      auto & aw = pdfToDo[k]; 
      while(true) {
	int ip =  aw;
	while (ip<int(nevals) && !std::atomic_compare_exchange_weak(&aw,&ip,ip+1));
	if (ip==int(nevals)) break;
	computeChunk(ip,lc-1);
	--pdfDone[lc-1];
      }
      lc = nChunks[meG];
      while (lc>0 && !std::atomic_compare_exchange_weak(&nChunks[meG],&lc,lc-1));
      if (lc<=0) break;
    }
    todo=reduction;
    
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




