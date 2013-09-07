#ifndef DATA
#define DATA

#include "Named.h"
#include "TMath.h"

#include "Variable.h"
#include "List.h"

#include <vector>
#include <iostream>

#include "tbb.h"

class Data : public Named {
 public:
  using Value_t = double;
  // using Value_t = float;

  Data(const Char_t* name, const Char_t* title, UInt_t size, Variable &var1);
  Data(const Char_t* name, const Char_t* title, UInt_t size, Variable &var1, Variable &var2);
  Data(const Char_t* name, const Char_t* title, UInt_t size, Variable &var1, Variable &var2, Variable &var3);
  Data(const Char_t* name, const Char_t* title, UInt_t size, List<Variable> &vars);
  virtual ~Data();

  void Push_back();
  inline UInt_t GetEntries() const { return m_data.size()/m_vars.GetSize(); }

  Bool_t Get(UInt_t iEvent);

  Bool_t DoVectors(bool force = false);
  inline Bool_t IsVectorized() const { return m_dataCPU!=nullptr; }
  const Value_t *GetCPUData(const Variable &var) const {

    auto index = m_vars.Index(var);
    /*
    if (index<0) {
      std::cerr << "Data for variable " << var.GetName() << " are not in the data sample " << GetName() << "!!!" << std::endl;
      return 0;
    }  

   if (!IsVectorized()) {
      std::cerr << "Data for variable " << var.GetName() << " of data sample " << GetName() << " are vectorized!!!" << std::endl;
     return 0;
   }
   */
   return (Value_t *)__builtin_assume_aligned(&m_dataCPU[index*GetEntries()],32); //  need to be protected for GetEntries()%8!=0
   // return &(m_dataCPU[index*GetEntries()]);

  }

  
 private:
  List<Variable> m_vars;
  
  VectorSTD(Value_t) m_data; // matrix container

  Value_t * m_dataCPU=nullptr; //!

};

#endif
