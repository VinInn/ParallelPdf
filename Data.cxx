#include "Data.h"


Data::Data(const Char_t* name, const Char_t* title, UInt_t size,
	   List<Variable> &vars) :
  Named(name,title), m_stride(stride(size))
{
  m_vars.AddElement(vars);
  m_data = (Value_t*)memalign(ALIGNMENT,vars.GetSize()*m_stride*sizeof(Value_t));
}


Data::~Data()
{ free(m_data); 
}



void Data::Push_back()
{
  assert(m_size<m_stride);
  m_vars.ResetIterator();
  auto iter = m_data+m_size;
  while (Variable *var = m_vars.Next()) {
    (*iter) = var->GetVal();
    iter+=m_stride;
  }
  ++m_size;
}

Bool_t Data::Get(UInt_t iEvent)
{
  if (iEvent>m_size) return kFALSE;

  auto iter = m_data+iEvent;
  m_vars.ResetIterator();
  while (Variable *var = m_vars.Next()) {
    var->SetVal(*iter);
    iter+=m_stride;
  }

  return kTRUE;

}

