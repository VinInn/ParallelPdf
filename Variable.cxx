
#include "Variable.h"

#include <iostream>

Variable::Variable(const Char_t* name, const Char_t* title, Double_t value) :
  Named(name,title), m_value(omp_get_max_threads(),value), m_error(0), m_errorLo(0), m_errorHi(0),
  m_min(value), m_max(value), m_isConstant(kTRUE)
{
  
}

Variable::Variable(const Char_t* name, const Char_t* title, Double_t min, Double_t max) :
  Named(name,title), m_value(omp_get_max_threads(),(min+max)/2.), m_error(0), m_errorLo(0), m_errorHi(0),
  m_min(std::min(min,max)), m_max(std::max(min,max)),
  m_isConstant(kFALSE)
{

}

Variable::Variable(const Char_t* name, const Char_t* title, Double_t value, Double_t min, Double_t max) :
  Named(name,title), m_value(omp_get_max_threads(), value), m_error(0), m_errorLo(0), m_errorHi(0),
  m_min(std::min(min,max)), m_max(std::max(min,max)),
  m_isConstant(kFALSE)
{}

Variable::~Variable()
{

}

Double_t Variable::SetAllVal(Double_t value)
{
  if (value<m_min || value>m_max)  {
    std::cerr << "Value " << value << " out of range for variable " << GetName()
	      << " [" << m_min << ", " << m_max << "]" << std::endl;
    if (value<m_min)
      value = m_min;
    else
      value = m_max;
  }
  
  for ( auto & v : m_value) v=value;
  return value;

}


Double_t Variable::SetVal(Double_t value)
{
  return m_value[omp_get_thread_num()] =value;
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
