#ifndef NLL_H
#define NLL_H

#include "AbsPdf.h"

#include "Named.h"
#include "TMath.h"

#include "Data.h"
#include <vector>
#include<atomic>


// template<unsigned int BLSIZE>
class NLL : public Named {
 public:
  NLL(const Char_t* name, const Char_t* title, Data &data, AbsPdf &pdf,
      int dyn, bool idocache);

  ~NLL();  

  // parallel esecution
  double GetVal(bool verify=true);

  //  sequential
  double GetVal(PdfState& state);


  inline AbsPdf *GetPdf() { return m_pdf; }
  void makeCache();


  inline void SetBlockEventsSize(UInt_t nBlockEvents) {
    m_nBlockEvents = nBlockEvents; 
  }

  unsigned int dataSize() const { return m_data->size();}

private:

  int RunEvaluationBlockSplittingStatic(PdfState const & state, std::vector<unsigned short> const & pdfs);
  int RunEvaluationBlockSplittingDynamic(PdfState const & state, std::vector<unsigned short> const & pdfs, std::atomic<int> * istart, int const * iend);


  // Sequential vectorized reduction using IntLog
  static void PartialNegReduction(TMath::IntLog &value,
				  const Double_t *  __restrict__ pResults, unsigned int bsize) {
    pResults = (double * __restrict__)__builtin_assume_aligned(pResults,ALIGNMENT);
    value = IntLogAccumulate(value, pResults, bsize);
  }

 private:

  Data *m_data;
  AbsPdf *m_pdf;
 
 private:


  std::vector<TMath::IntLog> m_logs; // result the partial parallel sums (still as exp,mant)

  double m_nLoops=0;
  std::vector<int> minLoop,maxLoop,aveLoop;


  UInt_t m_nBlockEvents;
  int m_ngroups=0;
  bool docache=true;
 
 
};

#endif
