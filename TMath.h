
#ifndef TMATH_HH
#define TMATH_HH

#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include <iostream>

// Basic type used by ROOT
typedef double Double_t;
//typedef float Double_t;
typedef char Char_t;
typedef unsigned int UInt_t;
typedef int Int_t;
typedef bool Bool_t;
typedef float Float_t;

// If we're not compiling with GNU C++, elide __attribute__
#if !defined(__GNUC__) || defined(__INTEL_COMPILER)
#define  __attribute__(x)  /*NOTHING*/
#endif


namespace TMath {

  inline Double_t Pi()       { return 3.14159265358979323846; }
  inline Double_t TwoPi()    { return 2.0 * Pi(); }
  inline Double_t PiOver2()  { return Pi() / 2.0; }
  inline Double_t PiOver4()  { return Pi() / 4.0; }
  inline Double_t InvPi()    { return 1.0 / Pi(); }
  inline Double_t RadToDeg() { return 180.0 / Pi(); }
  inline Double_t DegToRad() { return Pi() / 180.0; }
  inline Double_t Sqrt2()    { return 1.4142135623730950488016887242097; }
  inline Double_t LnGamma(Double_t z) { return lgamma(z); }
  inline Double_t Floor(Double_t x) { return std::floor(x); }


  inline Double_t Tan(Double_t x) { return std::tan(x); }
  inline Double_t Cos(Double_t x) { return std::cos(x); }
  inline Double_t Sin(Double_t x) { return std::sin(x); }
  inline Double_t Log(Double_t x) { return std::log(x); }
  inline Double_t Abs(Double_t x) { return std::abs(x); }
  inline Double_t Exp(Double_t x) { return std::exp(x); }
  inline Double_t Sqrt(Double_t x) { return std::sqrt(x); }
  inline Double_t Power(Double_t x, Double_t y) { return std::pow(x,y); }
  inline Double_t Erf(Double_t x) { return erf(x); }
  inline Double_t ATan(Double_t x) { return std::atan(x); }
  inline Double_t ATan2(Double_t y, Double_t x)
  { if (x != 0) return  atan2(y, x);
    if (y == 0) return  0;
    if (y >  0) return  Pi()/2;
    else        return -Pi()/2;
  }

  inline Double_t Min(Double_t a, Double_t b) { return std::min(a,b); }
  inline Int_t Min(Int_t a, Int_t b) { return std::min(a,b); }
  inline UInt_t Min(UInt_t a, UInt_t b) { return std::min(a,b); }

  // Does the sum between a and b, return the sum and the error due to rounding
  // Use the Kahan algorithm (Fast2Sum algorithm)
  // NOTE: the algorithm works only if abs(a)>=abs(b)
  // See: P. Kornerup at al, "On the Computation of Correctly-Rounded Sums" 
  // pp.155-160, 2009 19th IEEE Symposium on Computer Arithmetic, 2009
  //#pragma intel optimization_level 0
  //#pragma GCC optimize ("O0")

  //  inline Double_t KahanSummation(const Double_t a, const Double_t b, Double_t& error) 
  //    __attribute__((optimize("O0")));

  inline Double_t KahanSummation(const Double_t a, const Double_t b, Double_t& error) {
#ifdef __INTEL_COMPILER
#pragma float_control(precise,on)
#endif
    {
      Double_t sum(a + b);
      Double_t z(sum-a);
      error = b - z;
      return sum;
    }
  }

  // Does the sum between a and b, return the sum and the error due to rounding
  // Uses the Knuth algorithm (2Sum algorithm)
  // See: P. Kornerup at al, "On the Computation of Correctly-Rounded Sums" 
  // pp.155-160, 2009 19th IEEE Symposium on Computer Arithmetic, 2009
  //#pragma intel optimization_level 0
  //#pragma GCC optimize ("O0")
  //  inline Double_t KnuthSummation(const Double_t a, const Double_t b, Double_t& error)
  //    __attribute__((optimize("O0")));

  inline Double_t KnuthSummation(const Double_t a, const Double_t b, Double_t& error) {
#ifdef __INTEL_COMPILER
#pragma float_control(precise,on)
#endif
    {
      Double_t sum(a + b);
      Double_t bprime(sum - a);
      Double_t aprime(sum - bprime);
      Double_t deltaB(b - bprime);
      Double_t deltaA(a - aprime);
      error = deltaA + deltaB;
      return sum;
    }
  }


  // Does the accumulation of sum += value and calculate the error due to rounding
  // Use Kahan algorithm
  inline void KahanAccumulationAdd(double& sum, double& error, const double value) {
    sum = KahanSummation(sum,error+value,error);
  }

  // Does the accumulation of sum += value and calculate the error due to rounding
  // Use Knuth algorithm
  inline void KnuthAccumulationAdd(double& sum, double& error, const double value) {
    sum = KnuthSummation(sum,error+value,error);
  }

  // Does the accumulation of sum -= value and calculate the error due to rounding
  // Use Kahan algorithm
  inline void KahanAccumulationSub(double& sum, double& error, const double value) {
    sum = KahanSummation(sum,error-value,error);
  }

  // Does the accumulation of sum -= value and calculate the error due to rounding
  // Use Knuth algorithm
  inline void KnuthAccumulationSub(double& sum, double& error, const double value) {
#ifdef __INTEL_COMPILER
#pragma float_control(precise,on)
#endif
    sum = KnuthSummation(sum,error-value,error);
  }

  // Does the accumulation using Double-Double Knuth algorithm
  //#pragma intel optimization_level 0
  //#pragma GCC optimize ("O0")

  // Value and Error used in the DoubleDouble accumulation
  struct ValueAndError_t {
    ValueAndError_t() { value = 0.; error = 0.; }
    ValueAndError_t(double invalue) { value = invalue; error = 0.; }

    // Note: ignore the error on other
    // Use += to take in account this error
    ValueAndError_t &operator-=(const ValueAndError_t &other) {
#ifdef __INTEL_COMPILER
#pragma float_control(precise,on)
#endif
      {
	value = KnuthSummation(value,error-other.value,error);

	return *this;
      }
    }

    ValueAndError_t &operator+=(const ValueAndError_t &other) {
#ifdef __INTEL_COMPILER
#pragma float_control(precise,on)
#endif
      {
	Double_t t2(0); // t
	Double_t t1 = KnuthSummation(value,other.value,t2);
	t2 += (error+other.error);
	  
	value = t1+t2;
	error = t2-(value-t1);

	return *this;
      }
    }

    Double_t value;
    Double_t error;
  };

  //  inline ValueAndError_t DoubleDoubleAccumulation(const std::vector<ValueAndError_t> &values)
  //    __attribute__((optimize("O2")));

  inline ValueAndError_t DoubleDoubleAccumulation(const std::vector<ValueAndError_t> &values) {
#ifdef __INTEL_COMPILER
#pragma float_control(precise,on)
#endif
    {
      ValueAndError_t result;
      std::vector<ValueAndError_t>::const_iterator iterValues = values.begin();
      
      if (iterValues!=values.end()) {
	result = *(iterValues);
	Double_t t1(0); // s
	Double_t t2(0); // t
	
	for (iterValues++; iterValues!=values.end(); ++iterValues) {
	  t1 = KnuthSummation(result.value,iterValues->value,t2);
	  t2 += (result.error+iterValues->error);
	  
	  result.value = t1+t2;
	  result.error = t2-(result.value-t1);
	}
      }
      
      return result;
    }

  }

}

namespace ROOT {
  namespace Math {
    double landau_quantile(double z, double xi);

  }
}

#endif
