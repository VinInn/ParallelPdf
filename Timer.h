
#ifndef TIMER_H
#define TIMER_H

#include "mpisvc.h"
#include "openmp.h"

namespace Timer {

  enum Timer_type {
    TIMER_LOCAL=0,
    TIMER_GLOBAL
  };

  inline double Wtime(Timer_type type = TIMER_GLOBAL) { 
    // Local case
    if (type==TIMER_LOCAL || !MPISvc::IsInitialized())
      return OpenMP::GetTime();
    
    // Global case
    return MPISvc::GetTime();

  }

}

#endif
