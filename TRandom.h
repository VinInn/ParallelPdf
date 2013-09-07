#ifndef ROOT_TRandom
#define ROOT_TRandom

#include "Named.h"
#include "TMath.h"

class TRandom : public Named {

protected:
   UInt_t   fSeed;  //Random number generator seed

public:
   TRandom(UInt_t seed=65539);
   virtual ~TRandom();
   virtual  Int_t    Binomial(Int_t ntot, Double_t prob);
   virtual  Double_t BreitWigner(Double_t mean=0, Double_t gamma=1);
   virtual  void     Circle(Double_t &x, Double_t &y, Double_t r);
   virtual  Double_t Exp(Double_t tau);
   virtual  Double_t Gaus(Double_t mean=0, Double_t sigma=1);
   virtual  UInt_t   GetSeed() const {return fSeed;}
   virtual  UInt_t   Integer(UInt_t imax);
   virtual  Double_t Landau(Double_t mean=0, Double_t sigma=1);
   virtual  Int_t    Poisson(Double_t mean);
   virtual  Double_t PoissonD(Double_t mean);
   virtual  void     Rannor(Float_t &a, Float_t &b);
   virtual  void     Rannor(Double_t &a, Double_t &b);
   //   virtual  void     ReadRandom(const char *filename);
   virtual  void     SetSeed(UInt_t seed=0);
   virtual  Double_t Rndm(Int_t i=0);
   virtual  void     RndmArray(Int_t n, Float_t *array);
   virtual  void     RndmArray(Int_t n, Double_t *array);
   virtual  void     Sphere(Double_t &x, Double_t &y, Double_t &z, Double_t r);
   virtual  Double_t Uniform(Double_t x1=1);
   virtual  Double_t Uniform(Double_t x1, Double_t x2);
   //   virtual  void     WriteRandom(const char *filename);

};

#endif
