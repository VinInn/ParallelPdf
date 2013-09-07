
#include "mpisvc.h"

#include "Variable.h"
#include "List.h"
#include "Data.h"
#include "NLL.h"
#include "tbb.h"
#include "cilk.h"
#include "RooMinimizer.h"
#include "MsgService.h"

#include "TRandom.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include "common.h"

#include "models/extended1.h"
#include "models/model.h"
#include "models/gauss1.h"

#include "Timer.h"

std::string outputStatus(int status)
{
  std::stringstream buffer;

  if (status==-1) 
    buffer << "---";
  else
    buffer << status;

  return buffer.str();
}


double DoNLL(const unsigned int Iter, const unsigned int blockSize, Data &data, 
	     AbsPdf *model, bool runMinos,
	     bool externalLoop,
	     AbsPdf::DoCalculationBy doCalculationBy,
	     const char* label)
{
  // Do the calculation
  std::cout << "Runs the algorithm with `" << label << "'" << std::endl;
  std::cout << std::endl;

  NLL nll("nll","",data,*model,externalLoop,doCalculationBy);
  nll.SetBlockEventsSize(blockSize);
  RooMinimizer minimizer(nll);

  double value(0);
  int statusMigrad(-1), statusHesse(-1), statusMinos(-1);
  int callsAfterMigrad(0), callsAfterHesse(0), callsAfterMinos(0);
  double start = Timer::Wtime();

  if (Iter>0) {
    for (unsigned int i=0; i<Iter; i++)
      value += nll.GetVal();
  }
  else {
    model->RandomizeFloatParameters();
    statusMigrad = minimizer.migrad();
    callsAfterMigrad = minimizer.NumCallsFCN();
    statusHesse = minimizer.hesse();
    callsAfterHesse = minimizer.NumCallsFCN();
    if (runMinos) {
      statusMinos = minimizer.minos();
      callsAfterMinos = minimizer.NumCallsFCN();
    }
  }

  double end = Timer::Wtime();

  if (Iter==0) {
    //    value = minimizer.MinFCN();
    value = nll.GetVal();
    double edm = minimizer.Edm();

    std::cout << label << " # FCN Calls (After Migrad/Hesse/Minos) = " 
	      << callsAfterMigrad << "/" << callsAfterHesse << "/" << callsAfterMinos << std::endl;
    std::cout << "Status (Migrad/Hesse/Minos) = " 
	      << statusMigrad << "/" << statusHesse << "/" << outputStatus(statusMinos) << std::endl;
    std::cout << label << " # Invalid NLL = " << minimizer.NumInvalidNLL() << std::endl;
    std::cout << label << " Edm = " << std::setprecision(10) << edm << std::endl;

    List<Variable> pdfPars;
    nll.GetPdf()->GetParameters(pdfPars);
    pdfPars.Sort();
    pdfPars.Print(kTRUE);
    std::cout << std::endl;
  }

  std::cout << label << " Result = " << std::setprecision(10) << value << std::endl;
  std::cout << label << " Real Time (s) = " << std::setprecision(5) << end-start << std::endl;
  std::cout << std::endl;

  return value;

}

bool CheckResults(double value, double reference, const char* label)
{
  std::cout << "Check " << label << ": ";
  if ((value+reference)!=0. && (std::abs((value-reference)/(value+reference)))<1e-8) {
    std::cout << "OK" << std::endl;
    return kTRUE;
  }

  std::cout << "FAILED" << std::endl;
  return kFALSE;

}

AbsPdf *Model(Variable &x, Variable &y, Variable &z, const Int_t N)
{

  // Define the model
  //  AbsPdf *model = Extended1(x,y,z);
  AbsPdf *model = ModelEtapRGKs(x,y,z,N);
  //  AbsPdf *model = Gauss1(x);

  return model;
}


int main(int argc, char **argv)
{
  MPISvc::Init(); // Start MPI
  parallel_ostream::init(); // Only Master thread in master process can output

  // list of available commands
  if (FindOption(argc,argv,"-h")>=0) {
    std::cout << "Options:\n";
    std::cout << "-h to see this help\n";
    std::cout << "-n <int> to set the number of events (100K by default)\n";
    std::cout << "-m to run minos (false by default)\n";
    std::cout << "-e to run the external loop parallelization, only OpenMP (false by default)\n";
    std::cout << "-i <int> to set the number of iterations (100 by default)\n";
    std::cout << "-b <int> to set the number of events per block (default: 0 for OpenMP and Cilk, 1000 for TBB)\n";
    //    std::cout << "-k to run on MIC, offload mode (false by default)\n";
    std::cout << "-a <int> to choose the algorithm to run: 0=OpenMP (default), 1=TBB, 2=Cilk\n";
    std::cout << "-c to check the results against the virtual functions algorithm (false by default)\n";
    std::cout << std::endl;
    return 0;
  }
  const unsigned int N            = ReadIntOption(argc,argv,"-n",100000);
  const bool         runMinos     = (FindOption(argc,argv,"-m")>=0);
  const bool         externalLoop = (FindOption(argc,argv,"-e")>=0);
  const unsigned int Iter         = ReadIntOption(argc,argv,"-i",0);
  unsigned int blockSize          = ReadIntOption(argc,argv,"-b",0);
  //  const bool useMIC               = (FindOption(argc,argv,"-k")>=0);
  const unsigned int algo         = ReadIntOption(argc,argv,"-a",0);
  const bool useVirtualAlgo       = (FindOption(argc,argv,"-c")>=0);
  const char* datafile = "data1M.dat";

  // Define the variables
  Variable x("x","",-0.2,0.2); // DE
  Variable y("y","",5.25,5.29); // mES
  Variable z("z","",-3,1.5); // Fisher
  List<Variable> variables(x,y,z);

  // Fill the data
  Data data("data","",N,variables);

  if (Iter>0) {
    std::cout << "Generate " << N << " events..." << std::endl;
    TRandom rand;
    for (UInt_t i=0; i<N; i++) {
      variables.ResetIterator();
      while (Variable *var = variables.Next()) {
	var->SetVal(rand.Uniform(var->GetMin(),var->GetMax()));
      }
      data.Push_back();
    }
  } 
  else {
    std::cout << "Read " << N << " events for file " << datafile << "..." << std::endl;
    std::ifstream filedat(datafile);
    if (!filedat.is_open()) {
      std::cout << "Cannot read " << datafile << "! Abort..." << std::endl;
      return 0;
    }

    UInt_t i(0);
    for (i=0; i<N && !filedat.eof(); i++) {
      variables.ResetIterator();
      while (Variable *var = variables.Next()) {
	double value;
	filedat >> value;
	if (filedat.eof())
	  break;
        var->SetVal(value);
      }
      data.Push_back();
    }

    if (i<N) {
      std::cout << "Cannot read " << N << " events in " 
		<< datafile << "! Abort..." << std::endl;
      return 0;
    }

  }

  std::cout << std::endl 
	    << "# Events = " << std::fixed << std::setw(9) << data.GetEntries() 
	    << std::endl << std::endl;

  AbsPdf *model = Model(x,y,z,N);
  double valueAlgo(0);
  std::string label;

  if (algo==0) {
    // Do the calculation with OpenMP
    label = "OpenMP";
    valueAlgo = DoNLL(Iter,blockSize,data,model,runMinos,
		      externalLoop,
		      AbsPdf::kOpenMP,label.c_str());
  }
  else if (algo==1) {
    // Do the calculation with TBB
    label = "TBB";
    // Check if it was compiled with TBB
    TBBSafeCall(
		// Sequential case, for debugging only
		// tbb::task_scheduler_init init(1); // initiallize TBB task scheduler
                );
    if (blockSize<=0) 
      blockSize = 1000;
    valueAlgo = DoNLL(Iter,blockSize,data,model,runMinos,
		      externalLoop,
		      AbsPdf::kTBB,label.c_str());
  }
  else if (algo==2) {
    // Do the calculation with Cilk Plus
    label = "Cilk";
    // Check if it was compiled with Cilk
    CilkSafeCall(
		// Sequential case, for debugging only
                );
    // Cilk grainsize==0 means automatic splitting
    Cilk::grainsize = blockSize;
    if (externalLoop)
      valueAlgo = DoNLL(Iter,blockSize,data,model,runMinos,
			externalLoop,
			AbsPdf::kCilk_spawn,label.c_str());
    else
      valueAlgo = DoNLL(Iter,blockSize,data,model,runMinos,
			externalLoop,
			AbsPdf::kCilk_for,label.c_str());
  }

  if (useVirtualAlgo) {
    AbsPdf *modelVirtualAlgo = Model(x,y,z,N);
    // Do the calculation with Virtual functions
    double valueVirtual = DoNLL(Iter,blockSize,data,modelVirtualAlgo,runMinos,
				externalLoop,
				AbsPdf::kVirtual,"Virtual");

    std::cout << std::endl;

    label += "==Virtual";
    CheckResults(valueAlgo,valueVirtual,label.c_str());
  }

  std::cout << std::endl;

  model->ClearResults(kTRUE);
  delete model;
  
  parallel_ostream::cleanup();
  MPISvc::Finalize(); // Finalize MPI

  return 0;

}
