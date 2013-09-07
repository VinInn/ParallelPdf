
// Data analysis of branching fraction measurement of eta'(rg) K0S (PhysRevD.80.112002)

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <sstream>

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooBifurGauss.h"
#include "RooNLLVar.h"
#include "RooPolynomial.h"
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "TRandom.h"

#include <omp.h>

// Define the model
RooAbsPdf *ModelEtapRGKs(RooRealVar &x, RooRealVar &y, RooRealVar &z,
			 const Int_t N)
{
  // Var x
  
  RooRealVar *muA1x = new RooRealVar("muA1x"," ",-0.0018,-0.01,0.01); muA1x->setError(0.001);
  RooRealVar *sigmaA1x = new RooRealVar("sigmaA1x"," ",0.01512,0,0.1); sigmaA1x->setError(0.001);
  RooAbsPdf *gaussA1x = new RooGaussian("gaussA1x"," ",x,*muA1x,*sigmaA1x);

  RooRealVar *coeff1B1x = new RooRealVar("coeff1B1x"," ",-0.3156,-1,1); coeff1B1x->setError(0.001);
  RooArgList coefficientsB1x(*coeff1B1x);
  RooAbsPdf *polyB1x = new RooPolynomial("polyB1x"," ",x,coefficientsB1x);

  RooRealVar *coeff1C1x = new RooRealVar("coeff1C1x"," ",-0.7728);
  RooRealVar *coeff2C1x = new RooRealVar("coeff2C1x"," ",0.0067);
  RooRealVar *coeff3C1x = new RooRealVar("coeff3C1x"," ",0.1047); 
  RooRealVar *coeff4C1x = new RooRealVar("coeff4C1x"," ",-0.1120);
  RooArgList coefficientsC1x;
  coefficientsC1x.add(*coeff1C1x);
  coefficientsC1x.add(*coeff2C1x);
  coefficientsC1x.add(*coeff3C1x);
  coefficientsC1x.add(*coeff4C1x);
  RooAbsPdf *polyC1x = new RooPolynomial("polyC1x"," ",x,coefficientsC1x);

  RooRealVar *muD1x = new RooRealVar("muD1x"," ",0.1115);
  RooRealVar *sigmaD1x = new RooRealVar("sigmaD1x"," ",0.0464);
  RooAbsPdf *gaussD1x = new RooGaussian("gaussD1x"," ",x,*muD1x,*sigmaD1x);
  RooRealVar *coeff1D2x = new RooRealVar("coeff1D2x"," ",0.4146);
  RooArgList coefficientsD2x(*coeff1D2x);
  RooAbsPdf *polyD2x = new RooPolynomial("polyD2x"," ",x,coefficientsD2x);
  RooRealVar *fracDx = new RooRealVar("fracDx"," ",0.3821);
  RooAbsPdf *addDx = new RooAddPdf("addDx"," ",*gaussD1x,*polyD2x,*fracDx);

  RooRealVar *coeff1E1x = new RooRealVar("coeff1E1x"," ",-0.9392);
  RooRealVar *coeff2E1x = new RooRealVar("coeff2E1x"," ",-0.0793);
  RooRealVar *coeff3E1x = new RooRealVar("coeff3E1x"," ",0.2838); 
  RooRealVar *coeff4E1x = new RooRealVar("coeff4E1x"," ",-0.1428);
  RooArgList coefficientsE1x;
  coefficientsE1x.add(*coeff1E1x);
  coefficientsE1x.add(*coeff2E1x);
  coefficientsE1x.add(*coeff3E1x);
  coefficientsE1x.add(*coeff4E1x);
  RooAbsPdf *polyE1x = new RooPolynomial("polyE1x"," ",x,coefficientsE1x);

  // Var y

  RooRealVar *muA1y = new RooRealVar("muA1y"," ",5.2798,5.27,5.29); muA1y->setError(0.001);
  RooRealVar *sigmaA1y = new RooRealVar("sigmaA1y"," ",0.002640,0,0.01); sigmaA1y->setError(0.001);
  RooAbsPdf *gaussA1y = new RooGaussian("gaussA1y"," ",y,*muA1y,*sigmaA1y);
  
  RooRealVar *mB1y = new RooRealVar("mB1y"," ",y.getMax());
  RooRealVar *cB1y = new RooRealVar("cB1y"," ",-27.8171,-40,-10); cB1y->setError(0.001);
  RooAbsPdf *argusB1y = new RooArgusBG("argusB1y"," ",y,*mB1y,*cB1y);

  RooRealVar *mC1y = new RooRealVar("mC1y"," ",y.getMax());
  RooRealVar *cC1y = new RooRealVar("cC1y"," ",-65.2194);
  RooAbsPdf *argusC1y = new RooArgusBG("argusC1y"," ",y,*mC1y,*cC1y);
  RooRealVar *muC2y = new RooRealVar("muC2y"," ",5.2808); 
  RooRealVar *sigmaC2y = new RooRealVar("sigmaC2y"," ",0.0041); 
  RooAbsPdf *gaussC2y = new RooGaussian("gaussC2y"," ",y,*muC2y,*sigmaC2y);
  RooRealVar *fracCy = new RooRealVar("fracCy"," ",0.8576);
  RooAbsPdf *addCy = new RooAddPdf("addCy"," ",*argusC1y,*gaussC2y,*fracCy);

  RooRealVar *muD1y = new RooRealVar("muD1y"," ",5.2785);
  RooRealVar *sigmaD1y = new RooRealVar("sigmaD1y"," ",0.0054);
  RooAbsPdf *gaussD1y = new RooGaussian("gaussD1y"," ",y,*muD1y,*sigmaD1y);

  RooRealVar *mE1y = new RooRealVar("mE1y"," ",y.getMax());
  RooRealVar *cE1y = new RooRealVar("cE1y"," ",-61.2961);
  RooAbsPdf *argusE1y = new RooArgusBG("argusE1y"," ",y,*mE1y,*cE1y);
  RooRealVar *muE2y = new RooRealVar("muE2y"," ",5.2784);
  RooRealVar *sigmaE2y = new RooRealVar("sigmaE2y"," ",0.0050);
  RooAbsPdf *gaussE2y = new RooGaussian("gaussE2y"," ",y,*muE2y,*sigmaE2y);
  RooRealVar *fracEy = new RooRealVar("fracEy"," ",0.6793);
  RooAbsPdf *addEy = new RooAddPdf("addEy"," ",*argusE1y,*gaussE2y,*fracEy);
  
  // Var z

  RooRealVar *muA1z = new RooRealVar("muA1z"," ",-0.5518,-1,1); muA1z->setError(0.001);
  RooRealVar *sigmaA1z = new RooRealVar("sigmaA1z"," ",0.3314,0,1); sigmaA1z->setError(0.001);
  RooAbsPdf *gaussA1z = new RooGaussian("gaussA1z"," ",z,*muA1z,*sigmaA1z);

  RooRealVar *muB1z = new RooRealVar("muB1z"," ",-1.1352,-1.5,-0.5); muB1z->setError(0.001);
  RooRealVar *sigmaLB1z = new RooRealVar("sigmaLB1z"," ",0.3321,0,1); sigmaLB1z->setError(0.001);
  RooRealVar *sigmaRB1z = new RooRealVar("sigmaRB1z"," ",0.4441,0,1); sigmaRB1z->setError(0.001);
  RooAbsPdf *bifurgaussB1z = new RooBifurGauss("bifurgaussB1z"," ",z,*muB1z,*sigmaLB1z,*sigmaRB1z);
  RooAbsPdf *polyB1z = new RooPolynomial("polyB1z"," ",z);
  RooRealVar *fracBz = new RooRealVar("fracBz"," ",0.99);
  RooAbsPdf *addBz = new RooAddPdf("addBz"," ",*bifurgaussB1z,*polyB1z,*fracBz);

  RooRealVar *muC1z = new RooRealVar("muC1z"," ",-0.6762); 
  RooRealVar *sigmaLC1z = new RooRealVar("sigmaLC1z"," ",0.3241);
  RooRealVar *sigmaRC1z = new RooRealVar("sigmaRC1z"," ",0.3477);
  RooAbsPdf *bifurgaussC1z = new RooBifurGauss("bifurgaussC1z"," ",z,*muC1z,*sigmaLC1z,*sigmaRC1z);

  RooRealVar *muD1z = new RooRealVar("muD1z"," ",-0.6529);
  RooRealVar *sigmaLD1z = new RooRealVar("sigmaLD1z"," ",0.3472); 
  RooRealVar *sigmaRD1z = new RooRealVar("sigmaRD1z"," ",0.3577);
  RooAbsPdf *bifurgaussD1z = new RooBifurGauss("bifurgaussD1z"," ",z,*muD1z,*sigmaLD1z,*sigmaRD1z);

  RooRealVar *muE1z = new RooRealVar("muE1z"," ",-0.6336);
  RooRealVar *sigmaLE1z = new RooRealVar("sigmaLE1z"," ",0.3440);
  RooRealVar *sigmaRE1z = new RooRealVar("sigmaRE1z"," ",0.3570);
  RooAbsPdf *bifurgaussE1z = new RooBifurGauss("bifurgaussE1z"," ",z,*muE1z,*sigmaLE1z,*sigmaRE1z);

  
  RooRealVar *nA = new RooRealVar("nA"," ",10,0,N); nA->setError(1);
  RooRealVar *nB = new RooRealVar("nB"," ",40,0,N); nB->setError(1);
  RooRealVar *nC = new RooRealVar("nC"," ",30,0,N); nC->setError(1);
  RooRealVar *nD = new RooRealVar("nD"," ",10,0,N); nD->setError(1);
  RooRealVar *nE = new RooRealVar("nE"," ",10,0,N); nE->setError(1);
  RooArgList nevents;
  nevents.add(*nA); nevents.add(*nB); nevents.add(*nC); nevents.add(*nD); nevents.add(*nE);

  RooAbsPdf *pdfA = new RooProdPdf("pdfA"," ",RooArgList(*gaussA1x,*gaussA1y,*gaussA1z));
  RooAbsPdf *pdfB = new RooProdPdf("pdfB"," ",RooArgList(*polyB1x,*argusB1y,*addBz));
  RooAbsPdf *pdfC = new RooProdPdf("pdfC"," ",RooArgList(*polyC1x,*addCy,*bifurgaussC1z));
  RooAbsPdf *pdfD = new RooProdPdf("pdfD"," ",RooArgList(*addDx,*gaussD1y,*bifurgaussD1z));
  RooAbsPdf *pdfE = new RooProdPdf("pdfE"," ",RooArgList(*polyE1x,*addEy,*bifurgaussE1z));
  RooArgList Pdfs;
  Pdfs.add(*pdfA); Pdfs.add(*pdfB); Pdfs.add(*pdfC); Pdfs.add(*pdfD); Pdfs.add(*pdfE);
  
  return new RooAddPdf("extended"," ",Pdfs,nevents);
  //  return pdfA;
  //  return gaussA1x;
  //  return addDx;

}

int FindOption(int argc, char **argv, const char *option);
int ReadIntOption(int argc, char **argv, const char *option, int default_value = 1000);
char *ReadStringOption(int argc, char **argv, const char *option, char* default_value = 0);

std::string outputStatus(int status)
{

  std::stringstream buffer;

  if (status==-999) 
    buffer << "---";
  else
    buffer << status;

  return buffer.str();
}


int main(int argc, char **argv)
{

  // list of available commands
  if (FindOption(argc,argv,"-h")>=0) {
    std::cout << "Options:\n";
    std::cout << "-h to see this help\n";
    std::cout << "-n <int> to set the number of events (100K by default)\n";
    std::cout << "-t <int> to set the number of threads (1 by default)\n";
    std::cout << "-m to run minos (false by default)\n";
    std::cout << "-g to generate events (read from file by default)\n";
    std::cout << "-f <string> to set the file to read (\"data1M.dat\" by default)\n";
    std::cout << std::endl;
    return 0;
  }

  const unsigned int N         = ReadIntOption(argc,argv,"-n",100000);
  const unsigned int NTHREADS  = ReadIntOption(argc,argv,"-t",1);
  const bool         runMinos  = (FindOption(argc,argv,"-m")>=0);
  const bool         genEvents = (FindOption(argc,argv,"-g")>=0);
  const char*        datafile  = ReadStringOption(argc,argv,"-f","data1M.dat");

  RooRealVar x("x"," ",-0.2,0.2); // DE
  RooRealVar y("y"," ",5.25,5.29); // mES
  RooRealVar z("z"," ",-3,1.5); // Fisher
  RooArgList vars(x,y,z);
  RooAbsPdf *model = ModelEtapRGKs(x,y,z,N);
  //  model->graphVizTree("model.dot","",true, true);
  //  return 0;

  RooDataSet *data(0);
  
  if (genEvents) {
    data = (RooDataSet*)model->generate(vars,N);
    data->SetName("data");
    data->Print();
    //  data->write("data1M.dat");
    //  return 0;
  }
  else {
    std::cout << "Read " << N << " events for file " << datafile << "..." << std::endl;
    std::ifstream filedat(datafile);
    if (!filedat.is_open()) {
      std::cout << "Cannot read " << datafile << "! Abort..." << std::endl;
      return 0;
    }

    data = new RooDataSet("data","data",vars);
    UInt_t i(0);
    RooLinkedListIter iterVars = vars.iterator();
    RooRealVar *var(0);
    for (i=0; i<N && !filedat.eof(); i++) {
      iterVars.Reset();
      while (var = (RooRealVar*)iterVars.Next()) {
        double value;
        filedat >> value;
        if (filedat.eof())
          break;
        var->setVal(value);
      }

      data->add(vars);

    }

  }

  data->Print();

  RooNLLVar nll("nll","",*model,*data,RooFit::Extended(),RooFit::NumCPU(NTHREADS));
  //  RooNLLVar nll("nll","",*model,*data);

  RooMinimizer minimizer(nll);
  minimizer.optimizeConst(true);
  minimizer.setStrategy(1);
  minimizer.setMinimizerType("Minuit2");

  RooArgSet *paramSet = model->getParameters(RooArgSet());
  RooArgList *floatParamList = (RooArgList*)paramSet->selectByAttrib("Constant",kFALSE);
  floatParamList->remove(vars);

  // Randomize the float parameters
  TRandom rand;
  for (int i=0; i<floatParamList->getSize(); i++) {
    RooRealVar *par = (RooRealVar*)floatParamList->at(i);
    std::cout << par->GetName() << " = " << par->getVal();
    par->setVal(rand.Uniform(par->getMin(),par->getMax()));
    std::cout << " --> " << par->getVal() << std::endl;
  }

  double start = omp_get_wtime();
  int statusMigrad = minimizer.migrad();
  int statusHesse = minimizer.hesse();
  int statusMinos = -999; 
  if (runMinos)
    statusMinos = minimizer.minos();
  double end = omp_get_wtime();

  double value = nll.getVal();
  floatParamList->Print("v");

  RooFitResult *results = minimizer.save();
  results->Print("v");

  delete paramSet;
  delete floatParamList;
  delete data;
  delete results;

  std::cout << "Result = " << std::setprecision(20) << value << std::endl;
  std::cout << "Real Time (s) = " << std::setprecision(5) << end-start << std::endl;
  //  std::cout << "# FCN calls = " << minimizer.getNFcnCalls() << std::endl;
  std::cout << "Status (Migrad/Hesse/Minos) = " 
  	    << statusMigrad << "/" << statusHesse << "/" << outputStatus(statusMinos) << std::endl;

  return 0;

}

//
//  command line option processing
//
int FindOption(int argc, char **argv, const char *option)
{
  for (int i = 1; i<argc; i++)
    if(std::strcmp(argv[i],option)==0)
      return i;
  return -1;
}

int ReadIntOption(int argc, char **argv, const char *option, int default_value)
{
  int iplace = FindOption(argc,argv,option);
  if (iplace>=0 && iplace<argc-1)
    return std::atoi(argv[iplace+1]);
  return default_value;
}

char *ReadStringOption(int argc, char **argv, const char *option, char *default_value)
{
  int iplace = FindOption(argc,argv,option);
  if (iplace>=0 && iplace<argc-1)
    return argv[iplace+1];
  return default_value;
}

