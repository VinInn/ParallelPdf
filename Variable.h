#ifndef VARIABLE_H
#define VARIABLE_H


#include "PdfState.h"


#include "Named.h"
#include "TMath.h"


#include <omp.h>
#include <vector>

#define RooRealVar Variable

class Variable : public Named {
public:
  Variable(const Char_t* name, const Char_t* title, Double_t value);
  Variable(const Char_t* name, const Char_t* title, Double_t min, Double_t max);
  Variable(const Char_t* name, const Char_t* title, Double_t value, Double_t min, Double_t max);
  virtual ~Variable();

  inline Double_t GetVal() const { return m_value[omp_get_thread_num()]; }
  
  double value(PdfState const & state) const {
    return state.paramVal(num());
  }

  inline Double_t GetError() const { return m_error; }
  inline Double_t GetMin() const { return m_min; }
  inline Double_t GetMax() const { return m_max; }
  inline Double_t getMin() const { return m_min; }
  inline Double_t getMax() const { return m_max; }
  
  Double_t SetAllVal(Double_t value);
  Double_t SetVal(Double_t value);
  
  inline void SetError(Double_t error) { m_error = error; }
  inline void SetAsymError(std::pair<Double_t,Double_t> asymError) { m_errorLo = asymError.first; m_errorHi = asymError.second; }
  inline void setError(Double_t error) { m_error = error; }
  
  inline bool IsConstant() const { return m_isConstant; }
  inline void setConstant(Bool_t isConstant = kTRUE) { m_isConstant = isConstant; }
  bool isData() { return m_isData;}

  void Print();
  
protected:
  std::vector<double> m_value;
  
  Double_t m_error;
  Double_t m_errorLo;
  Double_t m_errorHi;
  Double_t m_min;
  Double_t m_max;
  
  bool m_isConstant;
  bool m_isData=false;
};

struct DataVariable : public Variable {
  template<typename... Args> 
  DataVariable(Args... args): Variable(args...){ m_isData=true;}
};

#endif
