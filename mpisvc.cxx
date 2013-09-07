#include "mpisvc.h"
#include <iostream>

MPISvc* MPISvc::m_instance = 0;
unsigned int MPISvc::m_size = 1;
unsigned int MPISvc::m_rank = 0;


MPISvc::MPISvc()
{
  MPISafeCall(
	      MPI::Init();
	      m_size = MPI::COMM_WORLD.Get_size();
	      m_rank = MPI::COMM_WORLD.Get_rank();
	      printf("Info --> MPISvc::ctor Start MPI on #%i/%i processor\n",
		     m_rank,m_size);
	      MPI::COMM_WORLD.Barrier();
	      );
}

MPISvc::~MPISvc()
{
  MPISafeCall(
	      MPI::COMM_WORLD.Barrier();
	      printf("Info --> MPISvc::ctor Finalize MPI on #%i/%i processor\n",
		     m_rank,m_size);
	      m_rank = 0; m_size = 1;
	      MPI::Finalize();
	      );
}
