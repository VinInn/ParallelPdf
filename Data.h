#ifndef DATA
#define DATA

#include "Named.h"
#include "TMath.h"

#include "openmp.h"


class Variable;

#include "List.h"

#include <vector>
#include <iostream>
#include <cassert>

class Data final : public Named {
 public:
  using Value_t = double;
  // using Value_t = float;

  static int nPartions;
  static size_t inPart() { return std::max(size_t(1),size_t(std::abs(nPartions)));}

  Data(){}
  Data(const Data&) = delete;
  Data & operator=(const Data&) = delete;
  
  Data(Data&& rh) : Named(rh.GetName(),rh.GetTitle()), 
    m_vars(std::move(rh.m_vars)), m_data(std::move(rh.m_data)), 
    m_stride(std::move(rh.m_stride)), m_size(std::move(rh.m_size)),
    m_capacity(std::move(rh.m_capacity)),m_start(std::move(rh.m_start)),
    m_totSize(rh.m_totSize) {}

  Data & operator=(Data&& rh) {
    Named::operator=(rh);
    std::swap(m_vars,rh.m_vars);
    std::swap(m_data,rh.m_data);
    std::swap(m_stride,rh.m_stride);
    std::swap(m_size,rh.m_size);
    std::swap(m_capacity, rh.m_capacity);
    std::swap(m_start,rh.m_start);
    std::swap(m_totSize,rh.m_totSize);
    return *this;
  } 
  

  Data(const Char_t* name, const Char_t* title, UInt_t size, List<Variable> &vars);
  Data(const Char_t* name, const Char_t* title, UInt_t size, UInt_t nvars);

  void allocate(UInt_t size, UInt_t nvars);

  virtual ~Data();

  unsigned int size() const { return m_totSize; }
  unsigned int capacity() const { auto k=0U; for (auto c:m_capacity) k+=c; return k;}

  void Push_back();
  inline UInt_t GetEntries() const { return m_totSize; }

  bool empty() const { return m_data.empty();}

  Value_t * data() { return m_data[partition()];}
  Value_t const * data() const { return m_data[partition()];}

  static size_t partition() {
    auto ig =  omp_get_thread_num()/(OpenMP::GetNumThreads()/inPart());
    assert(ig<inPart());
    return ig;
  }

  size_t startP() const {
    return m_start[partition()];
  }

  size_t sizeP() const {
    return m_size[partition()];
  }


 size_t startP(int i) const {
    return m_start[i];
  }
  size_t sizeP(int i) const {
    return m_size[i];
  }


  // find in which partition is this event
  size_t partition(unsigned int nev) const {
    for ( auto i=inPart(); i!=0; --i)
      if (nev>=m_start[i-1]) return i-1;
    return inPart();
  }



  Bool_t Get(UInt_t iEvent);

  Value_t const * GetData(const Variable &var, unsigned int dataOffset) const {
    auto index = std::find(m_vars.begin(),m_vars.end(),&var)-m_vars.begin();
    return GetData(index, dataOffset);
  }


  Value_t const * GetData(unsigned int index, unsigned int dataOffset) const {
    auto me = partition(dataOffset); auto off = dataOffset-m_start[me];
    return (Value_t const *)__builtin_assume_aligned(m_data[me]+off+index*m_stride[me],ALIGNMENT);
  }


  Value_t * GetData(unsigned int index, unsigned int dataOffset) {
    auto me = partition(dataOffset); auto off = dataOffset-m_start[me];
    return (Value_t *)__builtin_assume_aligned(m_data[me]+off+index*m_stride[me],ALIGNMENT);

  }

  
 private:
  std::vector<Variable*> m_vars;
  
  std::vector<Value_t *> m_data;

  std::vector<size_t>  m_stride;
  std::vector<size_t>  m_size;
  std::vector<size_t>  m_capacity;
  std::vector<size_t>  m_start;

  size_t m_totSize=0;

};

#endif
