#ifndef DATA
#define DATA

#include "Named.h"
#include "TMath.h"

#include "Variable.h"
#include "List.h"

#include <vector>
#include <iostream>
#include <cassert>

class Data : public Named {
 public:
  using Value_t = double;
  // using Value_t = float;

  Data(){}
  Data(const Data&) = delete;
  Data & operator=(const Data&) = delete;
  
  Data(Data&& rh) : Named(rh.GetName(),rh.GetTitle()), 
		    m_vars(std::move(rh.m_vars)), m_data(rh.m_data), 
		    m_stride(rh.m_stride), m_size(rh.m_size),
		    m_capacity(rh.m_capacity) { rh.m_data=nullptr;}
  Data & operator=(Data&& rh) {
    Named::operator=(rh);
    std::swap(m_vars,rh.m_vars);
    std::swap(m_data,rh.m_data);
    std::swap(m_stride,rh.m_stride);
    std::swap(m_size,rh.m_size);
    std::swap(m_capacity, rh.m_capacity);
    return *this;
  } 
  

  Data(const Char_t* name, const Char_t* title, UInt_t size, List<Variable> &vars);
  Data(const Char_t* name, const Char_t* title, UInt_t size, UInt_t nvars);

  virtual ~Data();

  unsigned int size() const { return m_size; }
  unsigned int capacity() const { return m_capacity;}

  void Push_back();
  inline UInt_t GetEntries() const { return m_size; }

  bool empty() const { return m_data==nullptr;}

  Bool_t Get(UInt_t iEvent);

  Value_t const * GetData(const Variable &var, unsigned int dataOffset) const {
    auto index = std::find(m_vars.begin(),m_vars.end(),&var)-m_vars.begin();
    return (Value_t const *)__builtin_assume_aligned(m_data+dataOffset+index*m_stride,ALIGNMENT);
  }


  Value_t const * GetData(unsigned int index, unsigned int dataOffset) const {
    return (Value_t const *)__builtin_assume_aligned(m_data+dataOffset+index*m_stride,ALIGNMENT);
  }


  Value_t * GetData(unsigned int index, unsigned int dataOffset) {
    return (Value_t *)__builtin_assume_aligned(m_data+dataOffset+index*m_stride,ALIGNMENT);
  }

  
 private:
  std::vector<Variable*> m_vars;
  
  Value_t * m_data=nullptr; //!

  unsigned int m_stride=0;
  unsigned int m_size=0;
  unsigned int m_capacity=0;

};

#endif
