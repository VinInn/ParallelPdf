#include "Results.h"

void Results::ClearAll()
{
  m_resultsCPU.clear();
  for (Int_t i=0; i<m_parallelResultsCPU.size(); i++) {
    free(m_parallelResultsCPU[i]);
    m_parallelResultsCPU[i] = 0;
  }

  m_parallelResultsCPU.clear();
  m_parallelSize.clear();
  m_parallelCapacity.clear();

}


