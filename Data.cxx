#include "Data.h"
#include "Partitioner.h"
#include "Variable.h"
#ifndef __APPLE__
#define NUMACTL
#endif
#ifdef NUMACTL
#include <numa.h>
#endif

int Data::nPartions=0;

Data::Data(const Char_t* name, const Char_t* title, UInt_t size,
	   List<Variable> &vars) :
  Named(name,title), m_vars(vars()),  m_data(inPart()),  m_stride(inPart(),0), m_size(inPart(),0), m_capacity(inPart(),0), m_start(inPart(),0)
{
  allocate(size,vars.GetSize());
}

Data::Data(const Char_t* name, const Char_t* title, UInt_t size, UInt_t nvars) :
    Named(name,title), m_data(inPart()),  m_stride(inPart(),0), m_size(inPart(),0), m_capacity(inPart(),0), m_start(inPart(),0){

  allocate(size, nvars);

}

void Data::allocate(UInt_t size, UInt_t nvars) {

  if ( 0==nPartions ) {
      auto  nev = size;
      auto me = 0U;
      m_stride[me] = stride(nev); 
      m_capacity[me]= nvars*m_stride[me]; 
      m_data[me]= (Value_t*)memalign(ALIGNMENT,m_capacity[me]*sizeof(Value_t));
      m_start[me]=0;
      return;
  }

  if ( 0>nPartions ) {
    for (auto me = 0U; me!=inPart(); ++me) {
      int ls=0; int le=0;
      auto nev = Partitioner::GetElements(inPart(),me,size,ls,le);
      m_stride[me] = stride(nev);
      m_capacity[me]= nvars*m_stride[me];
      m_data[me]= (Value_t*)memalign(ALIGNMENT,m_capacity[me]*sizeof(Value_t));
      // force the OS to allocate physical memory for the region
      memset(m_data[me], -1, m_capacity[me]*sizeof(Value_t));
      m_start[me]=ls;
      // mess up the interleave (DOES NOT)
      //if (me==0U) {
      // auto pp = (Value_t*)memalign(ALIGNMENT,m_capacity[me]*sizeof(Value_t));
      // memset(pp, -1, m_capacity[me]*sizeof(Value_t));
      //}
    }
    return;
  }

#pragma omp parallel
  {
    // assume each thread will allocate in its own NUMA side
    // select one for each partion
    bool t0 = 0== omp_get_thread_num()%(OpenMP::GetNumThreads()/inPart());
    if (t0) {
#ifdef NUMACTL
//      numa_setlocal_memory();
#endif
      auto me = partition();
      int ls=0; int le=0;
      auto nev = Partitioner::GetElements(inPart(),me,size,ls,le);
      m_stride[me] = stride(nev); 
      m_capacity[me]= nvars*m_stride[me];
#ifdef NUMACTL
      m_data[me]= (Value_t*)numa_alloc_local(m_capacity[me]*sizeof(Value_t));
      memset(m_data[me], -1, m_capacity[me]*sizeof(Value_t));  // required???
#else
      m_data[me]= (Value_t*)memalign(ALIGNMENT,m_capacity[me]*sizeof(Value_t));
      // force the OS to allocate physical memory for the region
      memset(m_data[me], -1, m_capacity[me]*sizeof(Value_t));
#endif
      assert(m_data[me]!=nullptr);
      assert(0==((size_t)(m_data[me])&(size_t)(ALIGNMENT-1)));
      m_start[me]=ls;
    }

  }
  for (auto me = 0U; me!=inPart(); ++me) assert(m_data[me]!=nullptr);
}


#ifdef NUMACTL
Data::~Data() { int k=0; for(auto d:m_data) numa_free(d,m_capacity[k++]*sizeof(Value_t) ); }
#else
Data::~Data() { for(auto d:m_data) free(d); }
#endif


void Data::Push_back()
{
  auto me = partition(m_totSize);
  assert(me<inPart());
  assert(m_size[me]<m_stride[me]);
  auto iter = m_data[me]+m_size[me];
  for (auto var : m_vars) {
    (*iter) = var->GetVal();
    iter+=m_stride[me];
  }
  ++m_size[me];
  ++m_totSize;
}

Bool_t Data::Get(UInt_t iEvent)
{
  if (iEvent>m_totSize) return kFALSE;
  auto me = partition(iEvent);
  iEvent -=m_start[me];
  assert(iEvent>=0);
  auto iter = m_data[me]+iEvent;
  for (auto var : m_vars) {
    var->SetAllVal(*iter);
    iter+=m_stride[me];
  }

  return kTRUE;

}

