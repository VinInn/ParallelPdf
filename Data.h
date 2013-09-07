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

  Data(const Char_t* name, const Char_t* title, UInt_t size, List<Variable> &vars);
  virtual ~Data();

  void Push_back();
  inline UInt_t GetEntries() const { return m_size; }

  Bool_t Get(UInt_t iEvent);

  Value_t const * GetData(const Variable &var, unsigned int dataOffset) const {
    auto index = m_vars.Index(var);
    return (Value_t const *)__builtin_assume_aligned(m_data+dataOffset+index*m_stride,ALIGNMENT);
  }

  
 private:
  List<Variable> m_vars;
  
  Value_t * m_data=nullptr; //!

  unsigned int m_stride=0;
  unsigned int m_size=0;

};

#endif
