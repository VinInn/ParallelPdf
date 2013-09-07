#include "mpisvc.h"
#include "MsgService.h"
#include <iostream>
#include <cstdio>
#include "openmp.h"

parallel_ostream *parallel_ostream::m_instance[] = {0};
parallel_ostream::ostream_mask parallel_ostream::m_mask = parallel_ostream::ostream_mask(0,0);
parallel_ostream::ostream_mask parallel_ostream::m_mask_previous = parallel_ostream::ostream_mask(0,0);

parallel_ostream& parallel_ostream::init(ostream_type type, ostream_mask mask, bool silent)
{
  setvbuf(stdout,NULL,_IONBF,0);

  static bool initSingle(true);
  if (initSingle)
    setmask(mask,silent);

  if (type==kALL) { // initialize all streams
    initSingle = false;
    for (int itype = 0; itype<type; itype++)
      init(static_cast<ostream_type>(itype));
    initSingle = true;
    return *(Instance(static_cast<ostream_type>(0))); // return the first instance
  }
  
  if (is_initialized(type))
    return *(Instance(type));
  
  switch (type) {
  case kCLOG:
    return *(Instance(type) = new parallel_ostream(std::clog));
  default:
    return *(Instance(type) = new parallel_ostream(std::cout));
  }

}

parallel_ostream& parallel_ostream::init(ostream_type type, int omp_mask, int mpi_mask, bool silent)
{
  return init(type,ostream_mask(omp_mask,mpi_mask),silent);
}

void parallel_ostream::cleanup(ostream_type type)
{
  if (type==kALL) { // clean all streams
    for (int itype = 0; itype<type; itype++)
      cleanup(static_cast<ostream_type>(itype));
    return;
  }
      
  if (!is_initialized(type))
    return;
  
  delete Instance(type);
  Instance(type) = 0;

}

void parallel_ostream::setmask(ostream_mask mask, bool silent) 
{
  m_mask_previous = m_mask; // store the previou mask
  m_mask = mask;

  if (!silent)
    std::printf("******************** SET STREAM MASK ********************\n");
  unsigned int nthread(1);
  if (m_mask.omp<0) {
    nthread = -m_mask.omp;
    if (!silent)
      printf("OpenMP: all threads can output. Set the maximum number of ostream to %i.\n",nthread);
  }
  else {
    if (!silent)
      printf("OpenMP: only thread with id=%i can output.\n",m_mask.omp);
  }

  for (int itype = 0; itype<kALL; itype++) {
    if (!is_initialized(static_cast<ostream_type>(itype)))
      continue;
    Instance(static_cast<ostream_type>(itype))->SetNStreams(nthread);
  }

  if (!silent) {
    if (MPISvc::UseMPI()) {
      if (m_mask.mpi<0)
	printf("   MPI: all processes can output.\n");
      else
	printf("   MPI: only proc with rank=%i can output.\n",m_mask.mpi);
    }

    printf("*********************************************************\n");
  }

}

parallel_ostream::parallel_ostream(std::ostream& stream) :
  m_stream(stream)
{
  SetNStreams(m_mask.omp>=0 ? 1 : -m_mask.omp);
  m_stdBuffer = m_stream.rdbuf(this);
}

parallel_ostream::~parallel_ostream()
{
  m_stream.rdbuf(m_stdBuffer);
}

void parallel_ostream::SetNStreams(unsigned int nstreams)
{
  if (m_buffer.size()==nstreams)
    return;

  m_buffer.resize(nstreams);
  for (unsigned int istream=0; istream<nstreams; istream++)
    m_buffer[istream].reserve(4096);
}

parallel_ostream::int_type parallel_ostream::overflow(parallel_ostream::int_type nChar)
{
  if (traits_type::not_eof(nChar) && 
      (m_mask.mpi<0 || m_mask.mpi==MPISvc::GetRankProc())) {
    if (m_mask.omp>=0 && m_mask.omp==OpenMP::GetRankThread()) {
      m_buffer[0].push_back(nChar);
    }
    else if (m_mask.omp<0 && OpenMP::GetRankThread()<m_buffer.size()) {
      m_buffer[OpenMP::GetRankThread()].push_back(nChar);
    }
  }

  return nChar;
}

int parallel_ostream::sync()
{
  if (m_mask.mpi<0 || m_mask.mpi==MPISvc::GetRankProc()) {
    if (m_mask.omp>=0 && m_mask.omp==OpenMP::GetRankThread()) {
      m_stdBuffer->sputn(m_buffer[0].c_str(),m_buffer[0].size());
      m_buffer[0].clear();
    }
    else if (m_mask.omp<0 && OpenMP::GetRankThread()<m_buffer.size()) {
#pragma omp critical
      {
	m_label << "(MPI=" << MPISvc::GetRankProc() 
		<< ", OMP=" << OpenMP::GetRankThread() 
		<< ") \t" << m_buffer[OpenMP::GetRankThread()];
	m_stdBuffer->sputn(m_label.str().c_str(),m_label.str().size());
	m_label.str(" ");
      }
      m_buffer[OpenMP::GetRankThread()].clear();
    }
  }

  return 0;
}

