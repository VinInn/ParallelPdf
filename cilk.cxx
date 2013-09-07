
#include "cilk.h"

namespace Cilk {

  // Default value for the grainsize
  int grainsize = 0;

  int GetNWorkers()
  { 
    CilkSafeCall(
		 return __cilkrts_get_nworkers();
		 );
    return 0;
  }

}
