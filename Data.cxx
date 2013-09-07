#include "Data.h"
#include<malloc.h>

Data::Data(const Char_t* name, const Char_t* title, UInt_t size,
	   Variable &var1) :
  Named(name,title)
{
  m_vars.AddElement(var1);
  m_data.reserve(size);
}

Data::Data(const Char_t* name, const Char_t* title, UInt_t size,
	   Variable &var1, Variable &var2) :
  Named(name,title)
{
  m_vars.AddElement(var1);
  m_vars.AddElement(var2);
  m_data.reserve(size*2);
}

Data::Data(const Char_t* name, const Char_t* title, UInt_t size,
	   Variable &var1, Variable &var2, Variable &var3) :
  Named(name,title)
{
  m_vars.AddElement(var1);
  m_vars.AddElement(var2);
  m_vars.AddElement(var3);
  m_data.reserve(size*3);
}

Data::Data(const Char_t* name, const Char_t* title, UInt_t size,
	   List<Variable> &vars) :
  Named(name,title)
{
  m_vars.AddElement(vars);
  m_data.reserve(size*vars.GetSize());
}


void Data::Push_back()
{
  m_vars.ResetIterator();
  while (Variable *var = m_vars.Next()) {
    m_data.push_back(var->GetVal());
  }

}

Bool_t Data::Get(UInt_t iEvent)
{
  if ((iEvent+1)*m_vars.GetSize()>m_data.size())
    return kFALSE;

  VectorSTD(Value_t)::iterator iter = m_data.begin()+iEvent*m_vars.GetSize();
  m_vars.ResetIterator();
  while (Variable *var = m_vars.Next()) {
    var->SetVal(*iter);
    ++iter;
  }

  return kTRUE;

}

Bool_t Data::DoVectors(bool force)
{
  // Something to do (possible cache misses)


  // if (m_dataCPU && !force)
  //  return kFALSE;
  
  VectorSTD(Value_t)::iterator iter = m_data.begin();
  m_dataCPU = (Value_t *)__builtin_assume_aligned(memalign(32,m_data.size()*sizeof(Value_t)),32);
  Value_t * iterCPU=m_dataCPU;  

  UInt_t iVar(0);
  UInt_t nVars = m_vars.GetSize();
  UInt_t nEvents = GetEntries();
  while (iter!=m_data.end()) {
    
    for (iVar = 0; iVar<nVars && iter!=m_data.end(); iVar++) {
      (*(iterCPU+iVar*nEvents)) = (*iter);
      ++iter;
    }
    ++iterCPU;
  }

  
  return kTRUE;

}

Data::~Data()
{ free(m_dataCPU); 
}

