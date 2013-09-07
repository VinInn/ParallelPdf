
#ifndef PARTITIONER
#define PARTITIONER

namespace Partitioner {

  // determine the number of elements per each process/thread (statically allocated)
  // determine also the first element
  inline int GetElements(int num, int rank, int nEvents, int& iStart, int& iEnd) {

    int numEventsIn = nEvents/num;
    int numEventsOut = nEvents%num;

    iStart = rank*numEventsIn; 
    if (rank<numEventsOut) {
      iStart += rank;
      iEnd = iStart + numEventsIn + 1;
    }
    else {
      iStart += numEventsOut;
      iEnd = iStart + numEventsIn;
    }

    return iEnd-iStart;
  }

}

#endif
