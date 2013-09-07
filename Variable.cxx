
#include "Variable.h"

#include <iostream>

Variable::Variable(const Char_t* name, const Char_t* title, Double_t value) :
  Named(name,title), m_value(value), m_error(0), m_errorLo(0), m_errorHi(0),
  m_min(value), m_max(value), m_isConstant(kTRUE)
{
  
}

Variable::Variable(const Char_t* name, const Char_t* title, Double_t min, Double_t max) :
  Named(name,title), m_value((min+max)/2.), m_error(0), m_errorLo(0), m_errorHi(0),
  m_min(std::min(min,max)), m_max(std::max(min,max)),
  m_isConstant(kFALSE)
{

}

Variable::Variable(const Char_t* name, const Char_t* title, Double_t value, Double_t min, Double_t max) :
  Named(name,title), m_value(value), m_error(0), m_errorLo(0), m_errorHi(0),
  m_min(std::min(min,max)), m_max(std::max(min,max)),
  m_isConstant(kFALSE)
{
  if (m_value<m_min) m_value = m_min;
  if (m_value>m_max) m_value = m_max;

}

Variable::~Variable()
{

}

Double_t Variable::SetVal(Double_t value)
{
  m_changed = (value!=m_value);
  if (!m_changed) return value;
  if (value>=m_min && value<=m_max)
    m_value = value;
  else {
    std::cerr << "Value " << value << " out of range for variable " << GetName()
	      << " [" << m_min << ", " << m_max << "]" << std::endl;
    if (value<m_min)
      value = m_min;
    else
      value = m_max;
  }
  
  return m_value;

}

void Variable::Print() {
  std::cout << GetName() << " = " << GetVal();
  if (!IsConstant()) {
    std::cout << " +/- " << m_error;
    if (m_errorLo!=0 || m_errorHi!=0)
      std::cout << " [" << m_errorLo << ", " << m_errorHi << "]";
    std::cout << " {" << GetMin() << ", " << GetMax() << "}";
  }
  std::cout << std::endl;

}
