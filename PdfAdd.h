
#ifndef PDFADD
#define PDFADD

#include "AbsPdf.h"
#include "List.h"
#include "TMath.h"
#include "Results.h"

#define RooAddPdf PdfAdd

class PdfAdd : public AbsPdf {
public:
  PdfAdd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2, Variable &fraction);
  PdfAdd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs, List<Variable> fractions);
  virtual ~PdfAdd() { }

  virtual const Results &GetValSIMD(UInt_t iStart = 0, UInt_t nPartialEvents = 0);

  virtual void GetParameters(List<Variable>& parameters);
  virtual void CacheIntegral();
  virtual void Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy = kOpenMP);

  virtual void ClearResults(Bool_t recursive = kFALSE);

 protected:
  virtual Double_t evaluate() const;
  virtual Double_t integral() const { return m_isExtended ? ExpectedEvents() : 1.; }

  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents, 
			      const Double_t invIntegral);

  virtual Double_t ExtendedTerm(UInt_t observed) const;
  virtual Bool_t IsExtended() const { return m_isExtended; }
  virtual Double_t ExpectedEvents() const;

 private:

  mutable List<AbsPdf> m_pdfs;
  mutable List<Variable> m_fractions;

  Bool_t m_isExtended;

  Bool_t RunEvaluate(Int_t nEvents, const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
		     const Results* pResultsPdf, UInt_t iPdf, Double_t coeff);

};

#endif
