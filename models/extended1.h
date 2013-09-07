
// This model does (n_s P_s + n_b P_b), where P_s = P_b = (G1a+G1b) G2 G3 
// (G is a Gaussian) 

#include "Variable.h"
#include "PdfGaussian.h"
#include "PdfProd.h"
#include "PdfAdd.h"

AbsPdf *Extended1(Variable &x, Variable &y, Variable &z)
{
  // Define the model
  Variable *mu1a = new Variable("mu1a","",0.5);
  Variable *sigma1a = new Variable("sigma1a","",2);
  AbsPdf *gauss1a = new PdfGaussian("gauss1a","",x,*mu1a,*sigma1a);
  
  Variable *mu1b = new Variable("mu1b","",0.5);
  Variable *sigma1b = new Variable("sigma1b","",2);
  AbsPdf *gauss1b = new PdfGaussian("gauss1b","",y,*mu1b,*sigma1b);
  Variable *frac1 = new Variable("frac1","",0.5);
  AbsPdf *add1 = new PdfAdd("add1","",*gauss1a,*gauss1b,*frac1);

  Variable *mu2 = new Variable("mu2","",0.5);
  Variable *sigma2 = new Variable("sigma2","",2);
  AbsPdf *gauss2 = new PdfGaussian("gauss2","",y,*mu2,*sigma2);
  
  Variable *mu3 = new Variable("mu3","",0.1);
  Variable *sigma3 = new Variable("sigma3","",8);
  AbsPdf *gauss3 = new PdfGaussian("gauss3","",z,*mu3,*sigma3);

  AbsPdf *prod123 = new PdfProd("prod123","",List<AbsPdf>(*add1,*gauss2,*gauss3));

  Variable *n1 = new Variable("n1","",100);
  Variable *n2 = new Variable("n2","",1000);

  return new PdfAdd("extended","",List<AbsPdf>(*prod123,*prod123),List<Variable>(*n1,*n2));

}


