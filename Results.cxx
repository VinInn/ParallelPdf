
#include "Results.h"

Results::Results() :
  m_resultsCPU(0), m_parallelResultsCPU(0)
{
}

void Results::ClearAll()
{
  m_resultsCPU.clear();
  for (Int_t i=0; i<m_parallelResultsCPU.size(); i++) {
    if (m_parallelResultsCPU[i]==0)
      continue;
    m_parallelResultsCPU[i]->clear();
    delete m_parallelResultsCPU[i];
    m_parallelResultsCPU[i] = 0;
  }

  m_parallelResultsCPU.clear();

}

