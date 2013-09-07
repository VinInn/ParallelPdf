
// This model does (n_s P_s + n_b P_b), where P_s = P_b = (G1a+G1b) G2 G3 
// (G is a Gaussian) 

#include "Variable.h"
#include "PdfGaussian.h"

AbsPdf *Gauss1(Variable &x)
{
  // Define the model
  Variable *mu = new Variable("mu","",0.5,-2,2);
  Variable *sigma = new Variable("sigma","",2);
  return new PdfGaussian("gauss","",x,*mu,*sigma);

}


