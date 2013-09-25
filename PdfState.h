#ifndef PdfState_H
#define PdfState_H

#include <cstddef>

class AbsPdf;
class Variable;
class Data;



/*
 * state of the pdf: it could be either reference or modified
 * a reference state has all chaches initialized
 * a modified state 
 *
 */
class PdfState {
public:
  virtual ~PdfState(){}

  // return value for Paramer i;
  virtual double paramVal(size_t i) const =0;
  // return integral for pdf i;
  virtual double invIntegral(size_t i) const =0;
  // fill res for pdf i;  res will contain the result, loc is a incoming "cache" to store values. if values stored in loc, res will point to loc 
  virtual double * pdfVal(size_t i,  double * __restrict__ loc, unsigned int bsize, const Data & data, unsigned int dataOffset) const =0;

  virtual void cacheIntegral(size_t i) const =0;
  virtual void cachePdf(size_t i, unsigned int bsize, const Data & data, unsigned int dataOffset)  const =0;


};


#endif
