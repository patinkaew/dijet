// Purpose: Fill gamma+jet analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: June 6, 2021
/*
//#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
*/
#include "DijetHistosFill.h"

#include "TSystem.h"

#include <fstream>
#include <string>

#define GPU

#ifdef LOCAL
// Compile these libraries into *.so first with root -l -b -q mk_CondFormats.C
/*
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)
*/
R__LOAD_LIBRARY(DijetHistosFill.C+g)
#else
R__LOAD_LIBRARY(DijetHistosFill_C.so)
#endif

void mk_DijetHistosFill(string dataset = "RunC") {

  // Settings
  bool addData = (dataset=="RunCearly" || dataset=="RunC");
  bool addMC = (dataset=="FlatQCD");

  //cout << "Clean old shared objects and link files" << endl << flush;
  //gSystem->Exec("rm *.d");
  //gSystem->Exec("rm *.so");
  //gSystem->Exec("rm *.pcm");	

  string path = gSystem->pwd();
  /*
  gSystem->AddIncludePath(Form("-I%s",path.c_str()));
  gSystem->AddIncludePath(Form("-I%s/CondFormats/JetMETObjects/interface",path.c_str()));
  */
#ifdef GPU
  // Compile these libraries into *.so first with root -l -b -q mk_CondFormats.C
  /*
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");
  */
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
    
    DijetHistosFill filler(c,0,dataset);
    filler.Loop();
  }
  
  if (addMC) {
    ifstream fin(runLocal ? "mcFiles_local.txt" :
		 Form("mcFiles_%s.txt",dataset.c_str()), ios::in);
    string filename;
    cout << "Chaining MC files:" << endl << flush;
    int nFiles(0), nFilesMax(100);
    while (fin >> filename && nFiles<nFilesMax) {
      ++nFiles;
      c->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;
  
    DijetHistosFill filler(c,1,dataset);
    filler.Loop();
  }

}
