
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
  Data(const Char_t* name, const Char_t* title, UInt_t size, Variable &var1);
  Data(const Char_t* name, const Char_t* title, UInt_t size, Variable &var1, Variable &var2);
  Data(const Char_t* name, const Char_t* title, UInt_t size, Variable &var1, Variable &var2, Variable &var3);
  Data(const Char_t* name, const Char_t* title, UInt_t size, List<Variable> &vars);
  virtual ~Data();

  void Push_back();
  inline UInt_t GetEntries() const { return m_data.size()/m_vars.GetSize(); }

  Bool_t Get(UInt_t iEvent);

  Bool_t DoVectors(bool force = false);
  inline Bool_t IsVectorized() const { return m_dataCPU.size()>0; }
  const Double_t *GetCPUData(const Variable &var) const;
  
 private:
  List<Variable> m_vars;
  
  VectorSTD(Double_t) m_data; // matrix container

  VectorSTD(Double_t) m_dataCPU; //!

};

#endif
