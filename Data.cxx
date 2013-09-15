#include "Data.h"
#include "Partitioner.h"


size_t Data::nPartions=1;

Data::Data(const Char_t* name, const Char_t* title, UInt_t size,
	   List<Variable> &vars) :
  Named(name,title), m_vars(vars()),  m_data(nPartions),  m_stride(nPartions,0), m_size(nPartions,0), m_capacity(nPartions,0), m_start(nPartions,0)
{
  allocate(size,vars.GetSize());
}

Data::Data(const Char_t* name, const Char_t* title, UInt_t size, UInt_t nvars) :
    Named(name,title), m_data(nPartions),  m_stride(nPartions,0), m_size(nPartions,0), m_capacity(nPartions,0), m_start(nPartions,0){

  allocate(size, nvars);

}

void Data::allocate(UInt_t size, UInt_t nvars) {
/*
  if ( 1==nPartions ) {
      auto  nev = size;
      auto me = 0U;
      m_stride[me] = stride(nev); 
      m_capacity[me]= nvars*m_stride[me]; 
      m_data[me]= (Value_t*)memalign(ALIGNMENT,m_capacity[me]*sizeof(Value_t));
      m_start[me]=0;
  }
*/
#pragma omp parallel
  {
    // assume each thread will allocate in its own NUMA side
    // select one for each partion
    bool t0 = 0== omp_get_thread_num()%(omp_get_num_threads()/nPartions);
    if (t0) {
      auto me = partition();
      int ls=0; int le=0;
      auto nev = Partitioner::GetElements(nPartions,me,size,ls,le);
      m_stride[me] = stride(nev); 
      m_capacity[me]= nvars*m_stride[me]; 
      m_data[me]= (Value_t*)memalign(ALIGNMENT,m_capacity[me]*sizeof(Value_t));
      m_start[me]=ls;
    }

  }

}


Data::~Data() { for(auto d:m_data) free(d); }



void Data::Push_back()
{
  auto me = partition(m_totSize);
  assert(me<nPartions);
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

