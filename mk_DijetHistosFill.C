// Purpose: Fill dijet analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: June 6, 2021

//#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "DijetHistosFill.h"

#include "TSystem.h"

#include <fstream>
#include <string>

#define GPU
//#define LOCAL

#ifdef LOCAL
// Compile these libraries into *.so first with root -l -b -q mk_CondFormats.C
// (works for 6.18.04?)
/*
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)
*/
//R__LOAD_LIBRARY(DijetHistosFill.C+g)
// As in jetphys/mk2_histosFill.C:
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc)

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty_cc)
//
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetResolutionObject_cc)
R__LOAD_LIBRARY(JetMETCorrections/Modules/src/JetResolution_cc)
//
R__LOAD_LIBRARY(DijetHistosFill_C)
#else
// (works for 6.26/10)
R__LOAD_LIBRARY(DijetHistosFill_C.so)
#endif

void mk_DijetHistosFill(string dataset = "X", string version = "vX") {

  // Get JMENANO from either location:
  // - lxplus:/eos/cms/store/group/phys_jetmet/JMENanoRun3/v2p1/JetMET
  // - Hefaistos:/media/DATA/JME_NANO_DATA
  
  if (!(dataset=="UL2016BCD" || dataset=="UL2016EF" || 
	dataset=="UL2016APVMG" ||
        dataset=="UL2016GH" || dataset=="UL2016MG" || dataset=="UL2016Flat" || 
	dataset=="UL2017B" || dataset=="UL2017C" || dataset=="UL2017D" ||
	dataset=="UL2017E" || dataset=="UL2017F" ||
	dataset=="UL2017MG" ||
	dataset=="UL2018A" || dataset=="UL2018B" || dataset=="UL2018C" ||
	//dataset=="UL2018D" ||
	dataset=="UL2018D1" || dataset=="UL2018D2" ||
	dataset=="UL2018MG" ||

	dataset=="UL2016BCD_ZB" || dataset=="UL2016EF_ZB" || dataset=="UL2016GH_ZB" ||
	dataset=="UL2017B_ZB" || dataset=="UL2017C_ZB" || dataset=="UL2017D_ZB" ||
	dataset=="UL2017E_ZB" || dataset=="UL2017F_ZB" ||
	dataset=="UL2018A_ZB" || dataset=="UL2018B_ZB" || dataset=="UL2018C_ZB" ||
	dataset=="UL2018D_ZB"
	)) {
    cout << "Dataset not supported" << endl << flush;
    cout << "Supported datasets are:" << endl
	 << "UL2016BCD, UL2016EF, UL2016APVMG" << endl
	 << "UL2016GH, UL2016MG, UL2016Flat," << endl
	 << "UL2017B, UL2017C, UL2017D, UL2017E, UL2017F" << endl
	 << "UL2017MG" << endl
      //<< "UL2018A, UL2018B, UL2018C, UL2018D" << endl
	 << "UL2018A, UL2018B, UL2018C, UL2018D1, UL2018D2" << endl
	 << "UL2018MG"
	 << endl
      	 << "UL2016BCD_ZB, UL2016EF_ZB, UL2016GH_ZB" << endl
	 << "UL2017B_ZB, UL2017C_ZB, UL2017D_ZB, UL2017E_ZB, UL2017F_ZB" << endl
      	 << "UL2018A_ZB, UL2018B_ZB, UL2018C_ZB, UL2018D_ZB" << endl
	 << endl;
  }
  
  // Settings
  bool addData =
    (dataset=="UL2016BCD" || dataset=="UL2016EF" ||
     dataset=="UL2016GH" ||
     dataset=="UL2017B" || dataset=="UL2017C" || dataset=="UL2017D" ||
     dataset=="UL2017E" || dataset=="UL2017F" ||
     dataset=="UL2018A" || dataset=="UL2018B" || dataset=="UL2018C" ||
     //dataset=="UL2018D"
     dataset=="UL2018D1" || dataset=="UL2018D2" ||

     dataset=="UL2016BCD_ZB" || dataset=="UL2016EF_ZB" || dataset=="UL2016GH_ZB" ||
     dataset=="UL2017B_ZB" || dataset=="UL2017C_ZB" || dataset=="UL2017D_ZB" ||
     dataset=="UL2017E_ZB" || dataset=="UL2017F_ZB" ||
     dataset=="UL2018A_ZB" || dataset=="UL2018B_ZB" || dataset=="UL2018C_ZB" ||
     dataset=="UL2018D_ZB"
     );
  bool addMC =
    (dataset=="UL2016APVMG" ||
     dataset=="UL2016MG"  || dataset=="UL2016Flat" ||
     dataset=="UL2017MG"  || dataset=="UL2018MG"
     ); 

  //cout << "Clean old shared objects and link files" << endl << flush;
  //gSystem->Exec("rm *.d");
  //gSystem->Exec("rm *.so");
  //gSystem->Exec("rm *.pcm");	

  string path = gSystem->pwd();

  gSystem->AddIncludePath(Form("-I%s",path.c_str()));
  gSystem->AddIncludePath(Form("-I%s/CondFormats/JetMETObjects/interface",path.c_str()));

#ifdef GPU
  // Compile these libraries into *.so first with root -l -b -q mk_CondFormats.C

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetResolutionObject.cc+");
  gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+");

  cout << "Load library in GPU mode" << endl << flush;
  gROOT->ProcessLine(".L DijetHistosFill.C+g");
#endif


  TChain *c = new TChain("Events");
  
  // Automatically figure out where we are running the job
  bool runGPU = (path=="/media/storage/dijet");
  bool runLocal = (path=="/Users/voutila/Dropbox/Cern/dijet" ||
		   path=="/Users/manvouti/Dropbox/Cern/dijet");
  if (!runLocal) assert(runGPU);
  
  if (addData) {
    ifstream fin(runLocal ? Form("dataFiles_local_%s.txt",dataset.c_str()) : 
		 Form("dataFiles_%s.txt",dataset.c_str()), ios::in);
    string filename;
    cout << "Chaining data files:" << endl << flush;
    int nFiles(0), nFilesMax(9999);
    while (fin >> filename && nFiles<nFilesMax) {
      ++nFiles;
      c->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;

    //bool isZB = (dataset=="UL2017B_ZB" || dataset=="UL2017C_ZB" || dataset=="UL2017D_ZB" ||
    //		 dataset=="UL2017E_ZB" || dataset=="UL2017F_ZB");
    // => decide internally from dataset.Contains("_ZB")
    
    DijetHistosFill filler(c,0,dataset,version);
    filler.Loop();
  }
  
  if (addMC) {
    ifstream fin(runLocal ? Form("mcFiles_local_%s.txt",dataset.c_str()) :
		 Form("mcFiles_%s.txt",dataset.c_str()), ios::in);
    string filename;
    cout << "Chaining MC files:" << endl << flush;
    int nFiles(0), nFilesMax(9999);
    while (fin >> filename && nFiles<nFilesMax) {
      ++nFiles;
      c->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;

    bool isMG = (dataset=="UL2016APVMG" || dataset=="UL2016MG" ||
		 dataset=="UL2017MG" || dataset=="UL2018MG");
    
    DijetHistosFill filler(c, isMG ? 2 : 1, dataset,version);
    filler.Loop();
  }

}
