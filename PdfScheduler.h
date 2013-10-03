#ifndef PdfScheduler_H
#define PdfScheduler_H

#include "PdfReferenceState.h"
#include "AbsPdf.h"
#include "Variable.h"
#include "List.h"

#include <thread>
#include<atomic>


#include<omp.h>

#include "CircularBuffer.h"


#include <thread>
#include <mutex>
typedef std::mutex Mutex;
// typedef std::lock_guard<std::mutex> Lock;
typedef std::unique_lock<std::mutex> Lock;
typedef std::unique_lock<std::mutex> Guard;
//typedef std::condition_variable Condition;



void differentiate(PdfReferenceState & refState, unsigned int nvar, int const * vars, double const * steps, double *  res);


class PdfScheduler {
  
public:

  enum What { integral, chunk, reduction};
 

  PdfScheduler(size_t inevals, PdfModifiedState const * imstates, TMath::IntLog * ivalues, size_t bsize) :
    m_nBlockEvents(bsize), nevals(inevals), mstates(imstates), values(ivalues),locks(nevals),
    istate(0), integToDo(inevals), integDone(inevals), 
    nPdfToEval(0), pdfToEval(inevals,-1),stateReady(inevals) {
    for (auto & l : locks) l=false;
    for ( auto k=0U; k!=integDone.size(); ++k)  { integToDo[k]=integDone[k]=mstates[k].size(); }
    setChunks();
  }
  
  inline void setChunks();

  inline ~PdfScheduler();

  size_t nBlockEvents() const { return m_nBlockEvents;}

  void doTasks() noexcept;
  
  void computeChunk(unsigned  ist, unsigned  icu) ;

  void chunkResult(size_t i, TMath::IntLog value);

private:
  
  size_t m_nBlockEvents;
  
  What todo=integral;
  
  size_t nevals;
  PdfModifiedState const * mstates;
  TMath::IntLog * values;
  std::vector<std::atomic<bool> > locks;

  std::atomic<int> istate;
  std::vector<std::atomic<int> > integToDo;
  std::vector<std::atomic<int> > integDone;

  std::atomic<unsigned int> nPdfToEval;
  std::vector<int> pdfToEval;

  std::vector<std::atomic<int> > nChunks;
  std::vector<int> iChunks;

  std::vector<std::atomic<int> >  pdfToDo;
  std::vector<std::atomic<int> >  pdfDone;


  std::vector<std::atomic<int> > stateReady;
  
};

inline
PdfScheduler::~PdfScheduler() {
#ifdef DOPRINT

  for ( auto k=0U; k!=integDone.size(); ++k)
    std::cout << k << ':' << integToDo[k] <<','<< integDone[k] << ' ';
  std::cout << std::endl;

  std::cout << nPdfToEval << " ";
  for ( auto k=0U; k!=pdfToEval.size(); ++k)
    std::cout << k << ':' << pdfToEval[k]<< ' ';;
  std::cout << std::endl;

  auto npar = Data::inPart();
  std::cout << "partions " << npar << ": ";
  for (auto i=0U; i<npar; ++i) std::cout << nChunks[i] <<','<< iChunks[i] <<' ';
  std::cout << "\n tot chunks " << iChunks.back() << ": ";
  for (auto i=0; i<iChunks.back(); ++i) std::cout << pdfToDo[i] <<','<< pdfDone[i] <<' ';
  std::cout << std::endl;


#endif
}

inline
void PdfScheduler::setChunks() {
  
  // auto allN = omp_get_num_threads();
  // auto meN = omp_get_thread_num();

  auto const & data = mstates[0].data();

  auto npar = Data::inPart();
  // auto ntot = data.size();

  auto chunk = 4*m_nBlockEvents;

  nChunks=std::vector<std::atomic<int> >(npar);
  iChunks.resize(npar+1,0);
  for (auto i=0U; i<npar; ++i) {
    nChunks[i] = data.sizeP(i)/chunk;
    if (0!= data.sizeP(i)%chunk) ++nChunks[i];
    iChunks[i+1] = iChunks[i]+nChunks[i];
  }
  pdfToDo=std::vector<std::atomic<int> >(iChunks.back());
  pdfDone=std::vector<std::atomic<int> >(iChunks.back());
  for( auto & a : pdfToDo) { a=0;}
  for( auto & a : pdfDone) { a=nevals;}


#ifdef DOPRINT
  std::cout << "partions " << npar << ": ";
  for (auto i=0U; i<npar; ++i) std::cout << nChunks[i] <<','<< iChunks[i] <<' ';
  std::cout << "\n tot chunks " << iChunks.back() << ": ";
  for (auto i=0; i<iChunks.back(); ++i) std::cout << pdfToDo[i] <<','<< pdfDone[i] <<' ';
  std::cout << std::endl;

#endif

}

#endif
