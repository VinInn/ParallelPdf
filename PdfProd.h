
#ifndef PDFPROD
#define PDFPROD

#include "AbsPdf.h"
#include "List.h"

#define RooProdPdf PdfProd

class PdfProd : public AbsPdf {
public:
	PdfProd(const Char_t* name, const Char_t* title, AbsPdf &pdf1, AbsPdf &pdf2);
	PdfProd(const Char_t* name, const Char_t* title, List<AbsPdf> pdfs);
  virtual ~PdfProd() { }
  
  virtual const Results &GetValSIMD(UInt_t iStart = 0, UInt_t nPartialEvents = 0);

  virtual void CacheIntegral();
  virtual void GetParameters(List<Variable>& parameters);
  virtual void Init(const Data& data, Bool_t *doLog, DoCalculationBy doCalculationBy = kOpenMP);

  virtual void ClearResults(Bool_t recursive = kFALSE);

 protected:
  virtual Double_t evaluate() const;
  virtual Double_t integral() const { return 1.; }

  virtual Bool_t evaluateSIMD(const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
			      const Double_t invIntegral);

 private:

  mutable List<AbsPdf> m_pdfs;

  Bool_t RunEvaluate(Int_t nEvents, const UInt_t& iPartialStart, const UInt_t& nPartialEvents,
		     const Results* pResultsPdf, UInt_t iPdf);


};

#endif
