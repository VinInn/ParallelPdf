
#ifndef MPISERVICE_H
#define MPISERVICE_H

// the purpose of this macro is shutting down MPI
#ifdef USE_MPI
#define MPISafeCall(p)      p
#define MPIElse(p)
#else
#define MPISafeCall(p)
#define MPIElse(p)          p
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "Partitioner.h"
#include "TMath.h"

class MPISvc {
 public:

  static inline void Init() {
    MPISafeCall(
		if (!IsInitialized())
		  m_instance = new MPISvc;
		);
  }

  static inline void Finalize() {
    if (IsInitialized()) {
      delete m_instance;
      m_instance = 0;
    }
  }
  
  static inline bool IsInitialized() { return m_instance!=0; }

  static inline bool UseMPI() { 
    MPISafeCall(
		return true; 
		);
    MPIElse(
	    return false;
	    );
  }

  static inline int GetRankProc() {
    MPISafeCall(
		Init();
		);
    return m_rank;
  }

  static inline int GetNumProcs() {
    MPISafeCall(
		Init();
		);
    return m_size;
  }

  static inline void Barrier() {
    MPISafeCall(
		Init();
		if (m_size>1)
		  MPI::COMM_WORLD.Barrier();
		); 
  }

  static inline double SumReduce(TMath::ValueAndError_t partialSum) {
    MPISafeCall(
		Init();
		if (m_size==1)
		  return partialSum.value;
		std::vector<TMath::ValueAndError_t> results(m_size);
		MPI::COMM_WORLD.Allgather(&partialSum,2,MPI::DOUBLE,
					  &results[0],2,MPI::DOUBLE);
		return TMath::DoubleDoubleAccumulation(results).value;
		//		return Reduce(partialSum.value,MPI::DOUBLE,MPI::SUM);
		);
    MPIElse(
	    return partialSum.value;
	    );
  }

  static inline double SumReduce(double partialSum) {
    MPISafeCall(
		return Reduce(partialSum,MPI::DOUBLE,MPI::SUM);
		);
    MPIElse(
	    return partialSum;
	    );
  }

  static inline int SumReduce(int partialSum) {
    MPISafeCall(
		return Reduce(partialSum,MPI::INT,MPI::SUM);
		);
    MPIElse(
	    return partialSum;
	    );
  }

  static inline int GetProcElements(int nEvents, int& iStart, int& iEnd) {
    if (m_size>1)
      return Partitioner::GetElements(GetNumProcs(),GetRankProc(),nEvents,iStart,iEnd);
    else {
      iStart = 0;
      iEnd = nEvents;
      return nEvents;
    }
  }

  static inline int GetProcElements(int nEvents) {
    int iStart(0), iEnd(0);
    return GetProcElements(nEvents,iStart,iEnd);
  }

  static inline double GetTime() {
    MPISafeCall(
		Init();
		// Sync timing
		if (m_size>1)
		  MPI::COMM_WORLD.Barrier();
		return MPI::Wtime();
		);
    MPIElse(
	    return 0.;
	    );
  }

 protected:
  MPISvc();
  virtual ~MPISvc();
  
  MPISafeCall(
	      template<typename T> static T Reduce(T partialSum, MPI::Datatype datatype, MPI::Op op) {
		T total(partialSum);
		MPISafeCall(
			    Init();
			    if (m_size>1)
			      MPI::COMM_WORLD.Allreduce(&partialSum,&total,1,datatype,op);
			    );
		
		return total;
	      }
	      );

 private:
  static MPISvc* m_instance;
  static unsigned int m_size;
  static unsigned int m_rank; 

};


#endif
