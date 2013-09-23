
#include "mpisvc.h"

#include "Variable.h"
#include "List.h"
#include "Data.h"
#include "NLL.h"
#ifdef DO_MINUIT
#include "RooMinimizer.h"
#endif
#include "MsgService.h"
#include "PdfReferenceState.h"

#include "TRandom.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include "common.h"

// #include "models/extended1.h"
#include "models/model.h"
// #include "models/gauss1.h"

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
	     AbsPdf & model, bool runMinos,
	     const char* label, int dynamic, bool docache, bool parderiv)
{
  // Do the calculation
  std::cout << "Runs the algorithm with `" << label << "'" << std::endl;
  std::cout << std::endl;

  assert(blockSize%8==0);
  // assert(blockSize<BLSIZE);

  NLL nll("nll","",data,model,dynamic,docache);
  nll.SetBlockEventsSize(blockSize);

#ifdef  DO_MINUIT
  RooMinimizer minimizer(nll);
#endif

  double value(0);
#ifdef  DO_MINUIT
  int statusMigrad(-1), statusHesse(-1), statusMinos(-1);
  int callsAfterMigrad(0), callsAfterHesse(0), callsAfterMinos(0);
#endif

  auto nEval=1U;

  double start = Timer::Wtime();

  if (Iter>0) {
#ifndef  DO_MINUIT
    // done in minimizer
    nll.makeCache();
#endif
    List<Variable> pdfPars;
    nll.GetPdf()->GetParameters(pdfPars);
    auto v1 = nll.GetVal(false); // init cache if needed
    auto v2 = nll.GetVal();

    // first count 
    int nvar=0;
    for ( auto ip=0U; ip!=pdfPars().size(); ++ip) if (!pdfPars()[ip]->IsConstant()) ++nvar;
    
    int var[nvar];
    int k=0;
    for ( auto ip=0U; ip!=pdfPars().size(); ++ip) if (!pdfPars()[ip]->IsConstant()) var[k++]= ip;
    assert(k==nvar);
   

    if (!parderiv) {
      {
	auto vr = pdfPars()[var[0]];
	auto v = vr->GetVal();
	auto e = vr->GetError();
	vr->SetAllVal(v-e);
	auto v3 = nll.GetVal();
	vr->SetAllVal(v+e);
	auto v4 = nll.GetVal();
	vr->SetAllVal(v);
	//	auto v5 = nll.GetVal();
	//auto v6 = nll.GetVal(false);
	std::cout << "init test " 
		  << v1 << " " 
		  << v2 << " " 
		  << v3 << " "  
		  << v4 << " "
	  //	  << v5 << " "
	  //	  << v6 << "\n"
		  << std::endl;
      }
      {
	auto vr = pdfPars()[var[nvar-1]];
	auto v = vr->GetVal();
	auto e = vr->GetError();
	vr->SetAllVal(v-e);
	auto v3 = nll.GetVal();
	vr->SetAllVal(v+e);
	auto v4 = nll.GetVal();
	vr->SetAllVal(v);
	auto v5 = nll.GetVal();
	auto v6 = nll.GetVal(false);
	std::cout << "init test " 
		  << v1 << " " 
		  << v2 << " " 
		  << v3 << " "  
		  << v4 << " "
		  << v5 << " "
		  << v6 << "\n"
		  << std::endl;
      }
    }else {
      auto v3 = nll.GetVal(var[0]);
      auto v4 = nll.GetVal(var[nvar-1]);
      std::cout << "init test " 
		<< v1 << " " 
		<< v2 << " " 
		<< v3 << " "  
		<< v4 << " "
		<< std::endl;
      auto vr = pdfPars()[var[0]];
      auto v = vr->GetVal();
      auto e = vr->GetError();
      vr->SetAllVal(v+e);
      v1 = nll.GetVal();
      vr->SetAllVal(v);
      v2 = nll.GetVal();
      vr->SetVal(v+e);
      v3 = nll.GetVal(var[0]);
      v4 = nll.GetVal(var[nvar-1]);
      vr->SetVal(v);
      std::cout << "init test " 
		<< v1 << " " 
		<< v2 << " " 
		<< v3 << " "  
		<< v4 << " "
		<< std::endl;
      vr = pdfPars()[var[nvar-1]];
      v = vr->GetVal();
      e = vr->GetError();
      vr->SetAllVal(v+e);
      v1 = nll.GetVal();
      vr->SetAllVal(v);
      v2 = nll.GetVal();
      vr->SetVal(v+e);
      v4 = nll.GetVal(var[nvar-1]);
      v3 = nll.GetVal(var[0]);
      vr->SetVal(v);
      auto v5 = nll.GetVal();
      std::cout << "init test " 
		<< v1 << " " 
		<< v2 << " " 
		<< v3 << " "  
		<< v4 << " "  
		<< v5 << "\n"
		<< std::endl;
      
      vr = pdfPars()[var[nvar-1]];
      v = vr->GetVal();
      e = vr->GetError();
      vr->SetVal(v+e);
      v3 = nll.GetVal(var[nvar-1]);
      vr->SetVal(v-e);
      v4 = nll.GetVal(var[nvar-1]);
      vr->SetVal(v);
      nll.GetPdf()->CacheIntegral(var[nvar-1]); // does not really matter as pdf cache is used...
      vr = pdfPars()[var[0]];
      v = vr->GetVal();
      e = vr->GetError();
      vr->SetVal(v+e);
      v1 = nll.GetVal(var[0]);
      vr->SetVal(v-e);
      v2 = nll.GetVal(var[0]);
      vr->SetVal(v);
      nll.GetPdf()->CacheIntegral(var[0]);
 
      v5 = nll.GetVal();
      std::cout << "init test " 
		<< v1 << " " 
		<< v2 << " " 
		<< v3 << " "  
		<< v4 << " "  
		<< v5 << "\n"
		<< std::endl;
      
    }
      





    
    for (unsigned int i=0; i<Iter; i++) {
      nll.GetVal(false);
      // fake computation of derivatives
      // double dvup[nvar*Data::ipar()]={0,}, vdown[nvar::ipar()]={0,};

      if (parderiv) {
	// outer parallel...
	auto nloops = 2*nvar;
	std::atomic<int> ok[nvar];for ( auto & o : ok) o=0;
	std::atomic<int> alP[Data::inPart()]; for ( auto & o : alP) o=0;
	double lval[OpenMP::GetMaxNumThreads()]={0.};
#pragma omp parallel
	{
	  auto par = Data::partition();
	  // auto start = data.startP();
	  // auto size   = data.sizeP();
	  std::atomic<int> & al = alP[par]; 
	  int il=0;
	  while(true) {
	    if (il>=nloops) break;
	    while (il<nloops && !std::atomic_compare_exchange_weak(&al,&il,il+1)); 
	    if (il>=nloops) break;
	    auto ik = il/2; // hope optmize in >1
	    bool pm = 0 == il%2;  // hope optmize in &1
	    auto vr = pdfPars()[var[ik]];
	    auto v = vr->GetVal();
	    auto e = vr->GetError();
	    auto nv = pm ? v+e : v-e;
	    vr->SetVal(nv);
	    if (pm) 
	      lval[omp_get_thread_num()] += nll.GetVal(var[ik]);
	    else 
	      lval[omp_get_thread_num()] -= nll.GetVal(var[ik]);

	    vr->SetVal(v);
	    // now the integral is wrong fir this thread...  (does not matter!)
	    nll.GetPdf()->CacheIntegral(var[ik]);
	    ok[ik]++;
	  }
	} // end parallel section
	for ( auto const & o : ok) assert(2*Data::inPart()==o);
	for (auto v: lval) value+=v;
      }
      else {
	for(int k=0; k!=nvar; ++k) {
	  auto vr = pdfPars()[var[k]];
	  auto v = vr->GetVal();
	  auto e = vr->GetError();
	  vr->SetAllVal(v+e);
	  value += nll.GetVal();
	  vr->SetAllVal(v-e);
	  value -= nll.GetVal();
	  vr->SetAllVal(v);
	}
      }

      nEval+=2*nvar+1;
    }
  }
  else {
    model.RandomizeFloatParameters();
#ifdef  DO_MINUIT
    statusMigrad = minimizer.migrad();
    callsAfterMigrad = minimizer.NumCallsFCN();
    statusHesse = minimizer.hesse();
    callsAfterHesse = minimizer.NumCallsFCN();
    if (runMinos) {
      statusMinos = minimizer.minos();
      callsAfterMinos = minimizer.NumCallsFCN();
    }
#endif
  }

  double end = Timer::Wtime();

#ifdef  DO_MINUIT
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
  }
#endif

  if (Iter>0) std::cout << "NUmber of Iterations " << Iter << "\nNUmber of Evaluation " << nEval << std::endl;

  List<Variable> pdfPars;
  
  nll.GetPdf()->GetParameters(pdfPars);
  pdfPars.Sort();
  pdfPars.Print(kTRUE);
  std::cout << std::endl;
 

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
  //  MPISvc::Init(); // Start MPI
  //parallel_ostream::init(); // Only Master thread in master process can output

  // list of available commands
  if (FindOption(argc,argv,"-h")>=0) {
    std::cout << "Options:\n";
    std::cout << "-h to see this help\n";
    std::cout << "-n <int> to set the number of events (100K by default)\n";
    std::cout << "-m to run minos (false by default)\n";
     std::cout << "-i <int> to set the number of iterations (0 by default)\n";
    std::cout << "-b <int> to set the number of events per block (default: 1024)\n";
    //    std::cout << "-k to run on MIC, offload mode (false by default)\n";
    std::cout << "-d dynamic scheduling (0 by default >1 is number of groups)\n";
    std::cout << "-c use cache (false by default)\n";
    std::cout << "-a numa affinity (0 by default; otherwise is number of partitions)\n";
    std::cout << "-p compute derivative in parallel (false by default)\n";

    std::cout << std::endl;
    return 0;
  }
  const unsigned int N            = ReadIntOption(argc,argv,"-n",100000);
  const bool         runMinos     = (FindOption(argc,argv,"-m")>=0);
  int Iter                        = ReadIntOption(argc,argv,"-i",0);
  unsigned int blockSize          = ReadIntOption(argc,argv,"-b",1024);
  //  const bool useMIC               = (FindOption(argc,argv,"-k")>=0);
  const int  dynamic              = ReadIntOption(argc,argv,"-d",0);
  bool docache              = (FindOption(argc,argv,"-c")>=0);
  const bool parderiv             = (FindOption(argc,argv,"-p")>=0);
  int numa                  = ReadIntOption(argc,argv,"-a",0);

  if (parderiv) docache = true;

  const char* datafile = "data1M.dat";

  // Define the variables
  DataVariable x("x","",-0.2,0.2); // DE
  DataVariable y("y","",5.25,5.29); // mES
  DataVariable z("z","",-3,1.5); // Fisher
  List<Variable> variables(x,y,z);

  Data::nPartions=numa;

  // Fill the data
  Data data("data","",N,variables);

  if (Iter>0) {
    std::cout << "Generate " << N << " events..." << std::endl;
    TRandom rand;
    for (UInt_t i=0; i<N; i++) {
      variables.ResetIterator();
      while (Variable *var = variables.Next()) {
	var->SetAllVal(rand.Uniform(var->GetMin(),var->GetMax()));
      }
      data.Push_back();
    }
  } 
  else {
    Iter = -Iter;
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
        var->SetAllVal(value);
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

  auto model = Model(x,y,z,N);
 
  PdfReferenceState::me().init();

  PdfReferenceState::me().print();


  std::string label;

  // Do the calculation with OpenMP
  label = "OpenMP";
  DoNLL(Iter,blockSize,data,*model,runMinos,label.c_str(),dynamic,docache,parderiv);

  delete model;

  std::cout << std::endl;
 
  // parallel_ostream::cleanup();
  // MPISvc::Finalize(); // Finalize MPI

  return 0;

}
