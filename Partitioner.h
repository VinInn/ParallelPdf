
#ifndef PARTITIONER
#define PARTITIONER

namespace Partitioner {

  // determine the number of elements per each process/thread (statically allocated)
  // determine also the first element
  inline int GetElements(int num, int rank, int nEvents, int& iStart, int& iEnd) {

  // allign to 4 double boudaries!
  int nEvents4 = nEvents/4;
  
    int numEventsIn = nEvents4/num;
    
    int numEventsOut = nEvents4%num;

    iStart = 4*rank*numEventsIn; 
    if (rank<numEventsOut) {
      iStart += 4*rank;
      iEnd = iStart + 4*(numEventsIn + 1);
    }
    else {
      iStart += 4*numEventsOut;
      iEnd = iStart + 4*numEventsIn;
    }
    if (rank==num-1) iEnd+= (nEvents-4*nEvents4);
    



    return iEnd-iStart;
  }

}

#endif
