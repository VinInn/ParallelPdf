#ifndef TMATH_HH
#define TMATH_HH

#include<cstdlib>

#ifdef __APPLE__
inline void * memalign(size_t a, size_t n) { void * p=0;  posix_memalign(&p,a,n); return p;}
#else
#include <malloc.h>
#endif
#ifndef ALIGNMENT
#define ALIGNMENT 32
#endif

// #define ASSUME_ALIGNED(x,i) __builtin_assume_aligned(x,i)
// #define ASSUME_ALIGNED(x,i) x


inline size_t stride(size_t size) {  return (ALIGNMENT/sizeof(double))*(size/(ALIGNMENT/sizeof(double))+1);} // ok ok I am wasting ALIGNMENT bytes!

#include <vdt/vdtMath.h>
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

  constexpr Double_t Pi()       { return 3.14159265358979323846; }
  constexpr Double_t TwoPi()    { return 2.0 * Pi(); }
  constexpr Double_t PiOver2()  { return Pi() / 2.0; }
  constexpr Double_t PiOver4()  { return Pi() / 4.0; }
  constexpr Double_t InvPi()    { return 1.0 / Pi(); }
  constexpr Double_t RadToDeg() { return 180.0 / Pi(); }
  constexpr Double_t DegToRad() { return Pi() / 180.0; }
  constexpr Double_t Sqrt2()    { return 1.4142135623730950488016887242097; }

  inline Double_t LnGamma(Double_t z) { return std::lgamma(z); }
  inline Double_t Floor(Double_t x) { return std::floor(x); }


#ifdef STD_MATH
  inline Double_t Tan(Double_t x) { return std::tan(x); }
  inline Double_t Cos(Double_t x) { return std::cos(x); }
  inline Double_t Sin(Double_t x) { return std::sin(x); }
  inline Double_t Log(Double_t x) { return std::log(x); }
  inline Double_t Exp(Double_t x) { return std::exp(x); }
  inline Double_t ATan(Double_t x) { return std::atan(x); }
  inline Double_t ATan2(Double_t y, Double_t x) { return std::atan2(y, x);  }
#else
  inline Double_t Tan(Double_t x) { return vdt::fast_tan(x); }
  inline Double_t Cos(Double_t x) { return vdt::fast_cos(x); }
  inline Double_t Sin(Double_t x) { return vdt::fast_sin(x); }
  inline Double_t Log(Double_t x) { return vdt::fast_log(x); }
  inline Double_t Exp(Double_t x) { return vdt::fast_exp(x); }
  inline Double_t ATan(Double_t x) { return vdt::fast_atan(x); }
  inline Double_t ATan2(Double_t y, Double_t x) { return vdt::fast_atan2(y, x);  }
#endif
  constexpr Double_t Sqrt(Double_t x) { return std::sqrt(x); }
  inline Double_t Abs(Double_t x) { return std::abs(x); }
  inline Double_t Power(Double_t x, Double_t y) { return std::pow(x,y); }
  inline Double_t Erf(Double_t x) { return erf(x); }

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


  ValueAndError_t DoubleDoubleAccumulation(const std::vector<ValueAndError_t> &values);
//      __attribute__((optimize("O2")));  // does not work with inline....


  union binary64 {
    binary64() : ui64(0) {};
    binary64(double ff) : f(ff) {};
    binary64(int64_t ii) : i64(ii){}
    binary64(uint64_t ui) : ui64(ui){}
    
    uint64_t ui64; /* unsigned int */                
    int64_t i64; /* Signed int */                
    double f;
  };
  
  inline
  void frex(double x, int & er, double & mr) {
    
    binary64 xx,m;
    xx.f = x;
    
    // as many integer computations as possible, most are 1-cycle only, and lots of ILP.
    int e=  int( ( xx.ui64 >> 52) & 0x7FF) -1023; // extract exponent
    m.ui64 = (xx.ui64 & 0x800FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL; // extract mantissa as an FP number
    
    long long adjust = (xx.ui64>>51)&1; // first bit of the mantissa, tells us if 1.m > 1.5
    m.i64 -= adjust << 52; // if so, divide 1.m by 2 (exact operation, no rounding)
    e += adjust;           // and update exponent so we still have x=2^E*y
    
    er = e;
    // now back to floating-point
    mr = m.f; // 
    // all the computations so far were free of rounding errors...
  }


  struct alignas(128) IntLog {
    IntLog() : iexp(0), mantissa(1.){}
    IntLog(int i, double m) : iexp(i), mantissa(m){}
    int iexp=0;
    double mantissa=0;
    double value() const { return iexp+std::log2(mantissa); }
  };

inline  IntLog IntLogAccumulate(IntLog start, double const * value, int nval) {
    int sf=start.iexp;; double mf=start.mantissa; 
    int NN = nval;
    for (int i=0;i<NN; i+=128) {
      double mi=1.f; 
      for (auto k=i; k!=std::min(i+128,NN); ++k) {
	double x = value[k];
	int er=0; double mr=0;  frex(x,er,mr);
	sf+=er; mi*=mr;
      }
      int ei=0; frex(mf*mi,ei,mf); sf+=ei;
    }
    return IntLog(sf,mf);
  }



}
namespace ROOT {
  namespace Math {
    double landau_quantile(double z, double xi);

  }
}

#endif
