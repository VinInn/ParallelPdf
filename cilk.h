
#ifndef CILK_H
#define CILK_H

#include "TMath.h"
#include <cstdlib>

#ifdef USE_CILK

#include <cilk/reducer_opadd.h>

#define CilkSafeCall(p)      p

#ifdef __cilk
#define __CILK
#endif

#else

#define CilkSafeCall(p)      std::cerr << "Not compiled with Cilk Plus!!!" << std::endl; exit(1);

#endif

namespace Cilk {
  
  extern int grainsize;

  template<typename Type, typename Value> 
    Type reduce(Value values, Int_t iStart, Int_t iEnd) {

    CilkSafeCall(
                 cilk::reducer_opadd<Type> result;
                 _Cilk_for (Int_t idx = iStart; idx<iEnd; idx++) {
                   result -= Type(values[idx]);
                 }
                 return result.get_value();
                 );

    return Type();

  }

  int GetNWorkers();

}

#endif
