// Compare 19Dec2023 vs 22Sep2023
// Only load a few branches (event ID, jet pt, eta, phi)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include "ROOT/RVec.hxx"
//#include "ROOT/RVec.hxx"
//#include "c++/v1/vector"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStopwatch.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "tools.h"
#include "tdrstyle_mod22.C"

#include <iostream>
#include <map>

using namespace std;
using namespace tools;

//bool patchJESA = true;

class evtid {
private:
  UInt_t run_, ls_;
  ULong64_t evt_;
public :
  evtid() : run_(0), ls_(0), evt_(0) {}
  evtid(UInt_t run, UInt_t ls, ULong64_t evt) : run_(run), ls_(ls), evt_(evt) {}
  bool operator()(evtid const& a, evtid const& b) const {
    if (a.run_ < b.run_) return true;
    if (a.run_ > b.run_) return false;
    if (a.ls_  < b.ls_)  return true;
    if (a.ls_  > b.ls_)  return false;
    return (a.evt_ < b.evt_);
  }
  UInt_t run() const { return run_; }
  UInt_t lbn() const { return ls_; }
  ULong64_t evt() const { return evt_; }
};

std::map<int, std::map<int, int> > _json;
bool LoadJSON(string json) {
  cout << "Processing LoadJSON() with " + json << endl << flush;
  ifstream file(json, ios::in);
  if (!file.is_open()) { assert(false); return false; }
  char c;
  string s, s2, s3;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  if (c!='{') return false;
  while (file >> s and sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    //if (_gh_debug) PrintInfo(Form("\"%d\": ",rn),true);

    while (file.get(c) and c==' ') {};
    //if (_gh_debug) { PrintInfo(Form("%c",c),true); assert(c=='['); }
    ++nrun;

    bool endrun = false;
    while (!endrun and file >> s >> s2 and (sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3 or (file >> s3 and sscanf((s+s2+s3).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3))) {

      s2 = s1;
      if (s2=="]") { file >> s3; s2 += s3; }

      //if (_gh_debug) PrintInfo(Form("[%d,%d,'%s']",ls1,ls2,s1),true);

      for (int ls = ls1; ls != ls2+1; ++ls) {
        _json[rn][ls] = 1;
        ++nls;
      }

      endrun = (s2=="]," || s2=="]}");
      //if (_gh_debug and !endrun and s2!=",") { PrintInfo(string("s1: ")+s2,true); assert(s2==","); }
    } // while ls
    //if (_gh_debug) PrintInfo("",true);

    if (s2=="]}") continue;
    //else if (_gh_debug and s2!="],") PrintInfo(string("s2: ")+s2,true);
    assert(s2=="],");
  } // while run
  //if (s2!="]}") { PrintInfo(string("s3: ")+s2,true); return false; }
  if (s2!="]}") { return false; }

  cout << "Called LoadJSON() with " << json << endl;
  cout << Form("Loaded %d good runs and %d good lumi sections\n",nrun,nls);
  return true;
} // LoadJSON


void compareLite(string run="2023D") {

  TDirectory *curdir = gDirectory;
  setTDRStyle();
  const char *crun = run.c_str();

  cout << endl << "Processing Run" << run << endl << flush;
  cout <<         "======================" << endl << flush;

  // Book tree tA (19Dec2023)
  // Also set JSON filter and jesAfix
  double jesAfix(1.);
  TChain *c_tA = new TChain("Events");
  cout << "A is 19Dec2023" << endl;
  //if (run=="2023BCv123" || run=="2023CCv4" || run=="2023D") {
  {
    //c_tA->AddFile("../data/dijet/data/JetMET0_Run2023D-19Dec2023-v1-NanoV12.root");
    LoadJSON("rootfiles/Cert_Collisions2022_355100_362760_Golden.json");
    LoadJSON("rootfiles/Cert_Collisions2023_366442_370790_Golden.json");

    string filename = Form("input_files/dataFiles_%s.txt.19Dec2023.v12",crun);
    ifstream fin(filename.c_str(), ios::in);

    cout << "Chaining data files for A: " << filename << endl << flush;
    int nFiles(0);
    while (fin >> filename) {
      ++nFiles;
      c_tA->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;
  }
  
  
  // Set branches to sort events
  TBranch *b_run_tA, *b_lbn_tA, *b_evt_tA;
  UInt_t run_tA, lbn_tA;
  ULong64_t evt_tA;
  c_tA->SetBranchAddress("run",&run_tA,&b_run_tA);
  c_tA->SetBranchAddress("luminosityBlock",&lbn_tA,&b_lbn_tA);
  c_tA->SetBranchAddress("event",&evt_tA,&b_evt_tA);

  cout << "Sort TA entries" << endl << flush;
  map<evtid, pair<Long64_t, Long64_t>, evtid> mtA;
  Long64_t nentries = c_tA->GetEntries();//Fast();
  cout << "..Processing " << nentries << " entries" << endl << flush;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = c_tA->LoadTree(jentry);
    if (ientry < 0) break;
    //nb = c_tA->GetEntry(jentry);   nbytes += nb;
    b_run_tA->GetEntry(ientry);
    b_lbn_tA->GetEntry(ientry);
    b_evt_tA->GetEntry(ientry);
    assert(run_tA);
    assert(lbn_tA);
    assert(evt_tA);

    mtA[evtid(run_tA, lbn_tA, evt_tA)]
      = pair<Long64_t, Long64_t>(jentry, ientry);

    if (jentry%1000000==0) cout << "." << flush;
  }
  cout << endl;
  cout << "Found " << mtA.size() << " unique entries" << endl;
  cout << endl;


  // Book TB tree (22Sep2023)
  TChain *c_tB = new TChain("Events");
  cout << "B is 22Sep2023" << endl;
  //if (run=="2023D") {
  {
    //c_tB->AddFile("../data/dijet/data/JetMET0_Run2023D-19Dec2023-v1-NanoV12.root"); // test tree A duplicate
    //c_tB->Add("../data/dijet/data/JetMET0_Run2023D-22Sep2023-v1-NanoV12.root");
    //LoadJSON("rootfiles/Cert_Collisions2023_366442_370790_Golden.json");

    string filename = Form("input_files/dataFiles_%s.txt.22Sep2023.v12",crun);
    ifstream fin(Form("input_files/dataFiles_%s.txt",crun), ios::in);
    
    cout << "Chaining data files for B: " << filename << endl << flush;
    int nFiles(0);
    while (fin >> filename) {
      ++nFiles;
      c_tB->AddFile(filename.c_str());
    }
    cout << "Chained " << nFiles <<  " files" << endl << flush;
  }

  // Set branches to sort events
  TBranch *b_run_tB, *b_lbn_tB, *b_evt_tB;
  UInt_t run_tB, lbn_tB;
  ULong64_t evt_tB;
  c_tB->SetBranchAddress("run",&run_tB,&b_run_tB);
  c_tB->SetBranchAddress("luminosityBlock",&lbn_tB,&b_lbn_tB);
  c_tB->SetBranchAddress("event",&evt_tB,&b_evt_tB);

  cout << "Sort TB entries" << endl << flush;
  map<evtid, pair<Long64_t, Long64_t>, evtid> mtB;
  nentries = c_tB->GetEntries();//Fast();
  cout << "..Processing " << nentries << " entries" << endl << flush;

  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = c_tB->LoadTree(jentry);
    if (ientry < 0) break;
    //nb = c_tB->GetEntry(jentry);   nbytes += nb;
    b_run_tB->GetEntry(ientry);
    b_lbn_tB->GetEntry(ientry);
    b_evt_tB->GetEntry(ientry);
    assert(run_tB);
    assert(lbn_tB);
    assert(evt_tB);

    mtB[evtid(run_tB, lbn_tB, evt_tB)]
      = pair<Long64_t, Long64_t>(jentry, ientry);

    if (jentry%1000000==0) cout << "." << flush;
  }
  cout << endl;
  cout << "Found " << mtB.size() << " unique entries" << endl;
  cout << endl;

  // Use sorted events to find matching ones
  cout << "Match TA to TB events" << endl;
  map<Long64_t, Long64_t> mAtoB;
  Long64_t ntotA = mtA.size();
  int failsAtoB(0), nAtoB(0);
  typedef map<evtid, pair<Long64_t, Long64_t> > IT;
  for (IT::const_iterator it = mtA.begin(); it != mtA.end(); ++it) {

    if (mtB.find(it->first)!=mtB.end()) {
      assert(mAtoB.find(it->second.first)==mAtoB.end());
      mAtoB[it->second.first] = mtB[it->first].first;
    }
    else if (++failsAtoB<10) {
      
      cout << "For tA("<<it->first.run()<<","<<it->first.lbn()<<","
	   <<it->first.evt()<<"), did not find a matching entry in tB" << endl;
    }
    if (++nAtoB%1000000==0) cout << "." << flush;
  }
  cout << "\nFound " << mAtoB.size() << " matching entries" << endl;

  /*
  map<Long64_t, Long64_t> mBtoA;
  Long64_t ntotB = mtB.size();
  int failsBtoA(0), int nBtoA(0);
  for (IT::const_iterator it = mtB.begin(); it != mtB.end(); ++it) {

    if (mtA.find(it->first)!=mtA.end()) {
      assert(mBtoA.find(it->second.first)==mBtoA.end());
      mBtoA[it->second.first] = mtA[it->first].first;
    }
    else if (++failsBtoA<10) {
      
      cout << "For tB("<<it->first.run()<<","<<it->first.lbn()<<","
	   <<it->first.evt()<<"), did not find a matching entry in tA" << endl;
    }
    if (++nBtoA%1000000==0) cout << "." << flush;
  }
  cout << "\nFound " << mBtoA.size() << " matching entries" << endl;
  */
  
  // Set unique branches needed from tree A
  //TBranch *b_flag_tA;//, *b_flag_tB;
  //TBranch *b_jttrg_tA;//, *b_jttrg_tB;
  //TBranch *b_jttrg500_tA, *b_jttrg450_tA, *b_jttrg400_tA, *b_jttrg320_tA;
  //TBranch *b_jttrg260_tA, *b_jttrg200_tA, *b_jttrg140_tA, *b_jttrg140_tA;
  //TBranch *b_jttrg80_tA, *b_jttrg60_tA, *b_jttrg40_tA;

  // Jet triggers
  TBranch        *b_HLT_PFJet40;
  TBranch        *b_HLT_PFJet60;
  TBranch        *b_HLT_PFJet80;
  TBranch        *b_HLT_PFJet140;
  TBranch        *b_HLT_PFJet200;
  TBranch        *b_HLT_PFJet260;
  TBranch        *b_HLT_PFJet320;
  TBranch        *b_HLT_PFJet400;
  TBranch        *b_HLT_PFJet450;
  TBranch        *b_HLT_PFJet500;
  TBranch        *b_HLT_PFJet550;

  // Set branches for MET filters
  TBranch          *b_Flag_goodVertices;
  TBranch          *b_Flag_globalSuperTightHalo2016Filter;
  TBranch          *b_Flag_EcalDeadCellTriggerPrimitiveFilter;
  TBranch          *b_Flag_BadPFMuonFilter;
  TBranch          *b_Flag_BadPFMuonDzFilter;
  TBranch          *b_Flag_hfNoisyHitsFilter;
  TBranch          *b_Flag_eeBadScFilter;
  TBranch          *b_Flag_ecalBadCalibFilter;
  
  // Set common branches needed from both trees
  TBranch *b_njt_tA, *b_njt_tB;
  TBranch *b_jtpt_tA, *b_jtpt_tB;
  TBranch *b_jteta_tA, *b_jteta_tB;
  TBranch *b_jtphi_tA, *b_jtphi_tB;
  TBranch *b_jtid_tA, *b_jtid_tB;
  TBranch *b_jtjes_tA, *b_jtjes_tB;

  // Set unique branches needed from tree A
  //Bool_t flag_tA;//, flag_tB;
  //Bool_t jttrg_tA;//, jttrg_tB;
  //Bool_t jttrg500_tA, jttrg450_tA, jttrg400_tA, jttrg320_tA;
  //Bool_t jttrg260_tA, jttrg200_tA, jttrg140_tA, jttrg140_tA;
  //Bool_t jttrg80_tA, jttrg60_tA, jttrg40_tA;

  // Jet triggers
  Bool_t          HLT_PFJet40;
  Bool_t          HLT_PFJet60;
  Bool_t          HLT_PFJet80;
  Bool_t          HLT_PFJet140;
  Bool_t          HLT_PFJet200;
  Bool_t          HLT_PFJet260;
  Bool_t          HLT_PFJet320;
  Bool_t          HLT_PFJet400;
  Bool_t          HLT_PFJet450;
  Bool_t          HLT_PFJet500;
  Bool_t          HLT_PFJet550;

  // Set MET filters
  Bool_t          Flag_goodVertices;
  Bool_t          Flag_globalSuperTightHalo2016Filter;
  Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
  Bool_t          Flag_BadPFMuonFilter;
  Bool_t          Flag_BadPFMuonDzFilter;
  Bool_t          Flag_hfNoisyHitsFilter;
  Bool_t          Flag_eeBadScFilter;
  Bool_t          Flag_ecalBadCalibFilter;

  // Set common branches needed from both trees
  Int_t njt_tA, njt_tB;
  //Int_t npv_tA(0), npv_tB(0);
  //
  const int njt = 100;
  Float_t jtpt_tA[njt], jtpt_tB[njt];
  Float_t jteta_tA[njt], jteta_tB[njt];
  Float_t jtphi_tA[njt], jtphi_tB[njt];
  UChar_t jtid_tA[njt], jtid_tB[njt];
  Float_t jtjes_tA[njt], jtjes_tB[njt];

  // Book tree A common branches
  c_tA->SetBranchAddress("nJet",&njt_tA,&b_njt_tA);
  c_tA->SetBranchAddress("Jet_pt",jtpt_tA,&b_jtpt_tA);
  c_tA->SetBranchAddress("Jet_eta",jteta_tA,&b_jteta_tA);
  c_tA->SetBranchAddress("Jet_phi",jtphi_tA,&b_jtphi_tA);
  c_tA->SetBranchAddress("Jet_jetId",jtid_tA,&b_jtid_tA);
  c_tA->SetBranchAddress("Jet_rawFactor",jtjes_tA,&b_jtjes_tA);

  // Jet triggers
  c_tA->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
  c_tA->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
  c_tA->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
  c_tA->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
  c_tA->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
  c_tA->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
  c_tA->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
  c_tA->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
  c_tA->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
  c_tA->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
  
  // MET filters
  c_tA->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
  c_tA->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
  c_tA->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
  c_tA->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
  c_tA->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
  c_tA->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
  c_tA->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
  c_tA->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
  
  // Book tree B unique branches
  c_tB->SetBranchAddress("nJet",&njt_tB,&b_njt_tB);
  c_tB->SetBranchAddress("Jet_pt",jtpt_tB,&b_jtpt_tB);
  c_tB->SetBranchAddress("Jet_eta",jteta_tB,&b_jteta_tB);
  c_tB->SetBranchAddress("Jet_phi",jtphi_tB,&b_jtphi_tB);
  c_tB->SetBranchAddress("Jet_jetId",jtid_tB,&b_jtid_tB);
  c_tB->SetBranchAddress("Jet_rawFactor",jtjes_tB,&b_jtjes_tB);
  
  const int nsample = 1;//100; // 0.5h
  const float frac = 1;//0.01;
  cout << "Pairing TA and TB" << endl << flush;
  cout << "Sampling 1/"<<nsample<<" of events" << endl << flush;
  cout << "Keeping first "<<frac*100<<"% of events" << endl << flush;

  // Speed up processing by selecting only used branches for tree A
  c_tA->SetBranchStatus("*",0);  // disable all branches

  // Activate tree A unique branches
  //c_tA->SetBranchStatus("Flag_Run3",1);
  c_tA->SetBranchStatus("run",1);
  c_tA->SetBranchStatus("luminosityBlock",1);
  c_tA->SetBranchStatus("event",1);
  //c_tA->SetBranchStatus("HLT_PFJet500",1);

  // Jet triggers
  c_tA->SetBranchStatus("HLT_PFJet40",1);
  c_tA->SetBranchStatus("HLT_PFJet60",1);
  c_tA->SetBranchStatus("HLT_PFJet80",1);
  c_tA->SetBranchStatus("HLT_PFJet140",1);
  c_tA->SetBranchStatus("HLT_PFJet200",1);
  c_tA->SetBranchStatus("HLT_PFJet260",1);
  c_tA->SetBranchStatus("HLT_PFJet320",1);
  c_tA->SetBranchStatus("HLT_PFJet400",1);
  c_tA->SetBranchStatus("HLT_PFJet450",1);
  c_tA->SetBranchStatus("HLT_PFJet500",1);
     
  // MET Filters
  c_tA->SetBranchStatus("Flag_goodVertices", 1);
  c_tA->SetBranchStatus("Flag_globalSuperTightHalo2016Filter", 1);
  c_tA->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
  c_tA->SetBranchStatus("Flag_BadPFMuonFilter", 1);
  c_tA->SetBranchStatus("Flag_BadPFMuonDzFilter", 1);
  c_tA->SetBranchStatus("Flag_hfNoisyHitsFilter", 1);
  c_tA->SetBranchStatus("Flag_eeBadScFilter", 1);
  c_tA->SetBranchStatus("Flag_ecalBadCalibFilter", 1);
  
  // Activate tree A common branches
  c_tA->SetBranchStatus("nJet",1);
  c_tA->SetBranchStatus("Jet_pt",1);
  c_tA->SetBranchStatus("Jet_eta",1);
  c_tA->SetBranchStatus("Jet_phi",1);
  c_tA->SetBranchStatus("Jet_jetId",1);
  c_tA->SetBranchStatus("Jet_rawFactor",1);

  // Speed up processing by selecting only used branches for tree B
  c_tB->SetBranchStatus("*",0);  // disable all branches
  //c_tB->SetBranchStatus("Flag_Run3",1);
  //c_tB->SetBranchStatus("run",1);
  //c_tB->SetBranchStatus("luminosityBlock",1);
  //c_tB->SetBranchStatus("event",1);
  //c_tB->SetBranchStatus("HLT_PFJet500",1);

  // Activate tree B common branches
  c_tB->SetBranchStatus("nJet",1);
  c_tB->SetBranchStatus("Jet_pt",1);
  c_tB->SetBranchStatus("Jet_eta",1);
  c_tB->SetBranchStatus("Jet_phi",1);
  c_tB->SetBranchStatus("Jet_jetId",1);
  c_tB->SetBranchStatus("Jet_rawFactor",1);

  // Results of interest
  // pT binning from JEC?
  /*
  double vx[] =
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, //638, 686, 737, 790, 846, 905, 967,
     //1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450};
     638, 790, 967, 1172, 1410, 1684, 2000, 2500, 3000, 3500,
     4000, 4500, 5000, 6000};
  int nx = sizeof(vx)/sizeof(vx[0])-1;
  */
  /*
  double vxw[] =
      {10, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 395, 468,
     548, 638, 790, 967, 1172, 1410, 1684, 2000, 2500, 3450};
  int nxw = sizeof(vxw)/sizeof(vxw[0])-1;
  */
  double vx[] =
    {10, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 395, 468,
     548, 638, 790, 967, 1172, 1410, 1684, 2000, 2500, 3000, 3500,
     4000, 4500, 5000, 6000};
  int nx = sizeof(vx)/sizeof(vx[0])-1;
  
  // Open file for outputting results
  TFile *f = new TFile(Form("rootfiles/compareLite_%s.root",run.c_str()),
		       "RECREATE");

  // Tag-and-probe method
  TProfile *pta_tp = new TProfile("pta_tp",";p_{T,tag};p_{T,A}",nx,vx);
  TProfile *pa_tp = new TProfile("pa_tp",";p_{T,tag};p_{T,A}/p_{T,tag}",nx,vx);
  TProfile *pb_tp = new TProfile("pb_tp",";p_{T,tag};p_{T,B}/p_{T,tag}",nx,vx);
  TProfile *pd_tp = new TProfile("pd_tp",";p_{T,tag};"
				 "(p_{T,B}-p_{T,A})/p_{T,tag}",nx,vx);
  TH2D *h2a_tp =new TH2D("h2a_tp",";p_{T,tag};p_{T,A}/p_{T,tag}",nx,vx,400,-1,3);
  TH2D *h2b_tp =new TH2D("h2b_tp",";p_{T,tag};p_{T,B}/p_{T,tag}",nx,vx,400,-1,3);
  TH2D *h2d_tp =new TH2D("h2d_tp",";p_{T,tag};"
			 "0.5*(p_{T,B}-p_{T,A})/p_{T,tag}",nx,vx,600,-3,3);
  
  // Direct match method in two variants
  TProfile *pjesa_dm = new TProfile("pjesa_dm",";p_{T,A};JES(A)",nx,vx);
  TProfile *pjesb_dm = new TProfile("pjesb_dm",";p_{T,B};JES(B)",nx,vx);
  TProfile *pta_dm = new TProfile("pta_dm",";p_{T,A};p_{T,A}",nx,vx);
  TProfile *pa_dm = new TProfile("pa_dm",";p_{T,A};p_{T,B}/p_{T,A}",nx,vx);
  TProfile *ptb_dm = new TProfile("ptb_dm",";p_{T,B};p_{T,A}",nx,vx);
  TProfile *pb_dm = new TProfile("pb_dm",";p_{T,B};p_{T,A}/p_{T,B}",nx,vx);
  TProfile *ptd_dm = new TProfile("ptd_dm",";(p_{T,B}+p_{T,A})/2;"
				  "p_{T,A}",nx,vx);
  TProfile *pd_dm = new TProfile("pd_dm",";(p_{T,B}+p_{T,A})/2;"
				 "(p_{T,B}-p_{T,A})/(p_{T,B}+p_{T,A})",nx,vx);
  TH2D *h2a_dm = new TH2D("h2a_dm",";p_{T,A};p_{T,B}/p_{T,A}",nx,vx,400,-1,3);
  TH2D *h2b_dm = new TH2D("h2b_dm",";p_{T,B};p_{T,A}/p_{T,B}",nx,vx,400,-1,3);
  TH2D *h2d_dm = new TH2D("h2d_dm",";(p_{T,B}+p_{T,A})/2;"
			  "(p_{T,B}-p_{T,A})/(p_{T,B}+p_{T,A})",nx,vx,600,-3,3);

  curdir->cd();
  
  TStopwatch t;
  t.Start();

  int nev = 0;
  int ngood = 0;
  int nmatch = 0;
  //int nj = 0;
  for (map<Long64_t,Long64_t>::const_iterator it = mAtoB.begin();
       it != mAtoB.end(); ++it) {

    if (++nev%10000==0) cout << "." << flush;
    if (nev%nsample!=0) continue;
    if (nev>frac*ntotA) continue;

    if (nev==1 || nev==10000 || nev==100000 || nev==1000000 || nev==5000000){
      cout << endl
	   << Form("Processed %ld events (%1.1f%%) in %1.0f sec. ETA:",
		   (long int)nev, 100.*nev/ntotA,
		   t.RealTime()) << endl;
      TDatime now; now.Set(now.Convert()+t.RealTime()*ntotA/nev);
      now.Print();
      t.Continue();
    }

    // Load matching entries
    Long64_t jentrytA = it->first;
    Long64_t jentrytB = it->second;
    
    if (jentrytA<0 || jentrytA>=ntotA) continue;
    Long64_t ientrytA = c_tA->LoadTree(jentrytA);
    if (ientrytA < 0) break;
    c_tA->GetEntry(jentrytA);
    if (jentrytB<0 || jentrytB>=ntotB) continue;
    Long64_t ientrytB = c_tB->LoadTree(jentrytB);
    if (ientrytA < 0) break;
    c_tB->GetEntry(jentrytB);

    // Temporary patch for missing MET filters
    //flag_tA = true;
    //flag_tB = true;

    // MET filters
    bool flag_tA =
      (Flag_goodVertices ||
       Flag_globalSuperTightHalo2016Filter ||
       Flag_EcalDeadCellTriggerPrimitiveFilter ||
       Flag_BadPFMuonFilter ||
       Flag_BadPFMuonDzFilter ||
       Flag_hfNoisyHitsFilter ||
       Flag_eeBadScFilter ||
       Flag_ecalBadCalibFilter);
    
    // Does the run/LS pass the latest JSON selection?
    if (_json[run_tA][lbn_tA]==0) {
      continue;
    }
    else
      ++ngood;
    
    // Loop over two leading jets to find probe pairs
    for (int i = 0; i != min(2,int(njt_tA)); ++i) {
      
      double pt = jtpt_tA[i];
      double jes = (1-jtjes_tA[i]);
      double eta = jteta_tA[i];
      double phi = jtphi_tA[i];
      bool idtightA = (jtid_tA[i]>=4);

      bool hasmatch = false;
      for (int j = 0; j != min(2,int(njt_tB)) && !hasmatch; ++j) {
	
	double ptB = jtpt_tB[j];
	double jesB = (1-jtjes_tB[j]);
	double etaB = jteta_tB[j];
	double phiB = jtphi_tB[j];
	bool idtightB = (jtid_tB[j]>=4);
	
	// Match probe jets with deltaR<R/cone2
	double dr = tools::oplus(delta_eta(eta,etaB),delta_phi(phi,phiB));
	if (dr < 0.20) {

	  ++nmatch;
	  hasmatch = true;
	  double ptave = 0.5 * (ptB + pt);

	  // Tag selection
	  bool istp(false);
	  double pttag(0);
	  if (i<2 && j<2 && njt_tA>1 && njt_tB>1) {

	    // Tag is the other one of the two leading jets
	    int k = (i==0 ? 1 : 0);
	    int l = (j==0 ? 1 : 0);
	    
	    double pttagA = jtpt_tA[k];
	    double pttagB = jtpt_tB[l];
	    pttag = 0.5 * (pttagA+pttagB);

	    double phiTA = jtphi_tA[k];
	    double phiTB = jtphi_tB[l];
	    double dphiTA = delta_phi(phiTA,phi);
	    double dphiTB = delta_phi(phiTB,phiB);
	    double dphitag = 0.5*(dphiTA + dphiTB);

	    double etaTA = jteta_tA[k];
	    double etaTB = jteta_tB[l];
	    double etatag = 0.5*(etaTA + etaTB);

	    double pt3A = (njt_tA>2 ? jtpt_tA[2] : 0);
	    double pt3B = (njt_tB>2 ? jtpt_tB[2] : 0);
	    double alphaTA = pt3A / pttagA;
	    double alphaTB = pt3B / pttagB;
	    double alphatag = (alphaTA>0 && alphaTB>0 ?
			       0.5*(alphaTA+alphaTB) : max(alphaTA,alphaTB));

	    //bool trigger = true; // pttag>600
	    bool trigger =
	      ((HLT_PFJet500 && pttag>638) ||
	       (HLT_PFJet450 && pttag>548) ||
	       (HLT_PFJet400 && pttag>468) ||
	       (HLT_PFJet320 && pttag>395) ||
	       (HLT_PFJet260 && pttag>300) ||
	       (HLT_PFJet200 && pttag>245) ||
	       (HLT_PFJet140 && pttag>196) ||
	       (HLT_PFJet80  && pttag>114) ||
	       (HLT_PFJet60  && pttag>84) ||
	       (HLT_PFJet40  && pttag>37));
	    //{10,15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 395, 468,
	    //548, 638, 790, 967, 1172, 1410, 1684, 2000, 2500, 3000, 3500,
	    //4000, 4500, 5000, 6000};
	    //
	    //{10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
	    //97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395,
	    //430, 468, 507, 548, 592,
	    //638, 790, 967, 1172, 1410, 1684, 2000, 2500, 3000, 3500,
	    //4000, 4500, 5000, 6000};
	    istp = (dphitag>2.7 && alphatag<0.3 && fabs(etatag)<1.3 &&
		    pt/pttag<1.35 && ptB/pttag<1.35 &&
		    pt/pttag>0.45 && ptB/pttag>0.45 && trigger);
	  } // tag selection

	  // Direct match probe trigger selections
	  bool trigD =
	    ((HLT_PFJet500 && ptave>638) ||
	     (HLT_PFJet450 && ptave>548) ||
	     (HLT_PFJet400 && ptave>468) ||
	     (HLT_PFJet320 && ptave>395) ||
	     (HLT_PFJet260 && ptave>300) ||
	     (HLT_PFJet200 && ptave>245) ||
	     (HLT_PFJet140 && ptave>196) ||
	     (HLT_PFJet80  && ptave>114) ||
	     (HLT_PFJet60  && ptave>84) ||
	     (HLT_PFJet40  && ptave>37));
	  bool trigB =
	    ((HLT_PFJet500 && ptB>638) ||
	     (HLT_PFJet450 && ptB>548) ||
	     (HLT_PFJet400 && ptB>468) ||
	     (HLT_PFJet320 && ptB>395) ||
	     (HLT_PFJet260 && ptB>300) ||
	     (HLT_PFJet200 && ptB>245) ||
	     (HLT_PFJet140 && ptB>196) ||
	     (HLT_PFJet80  && ptB>114) ||
	     (HLT_PFJet60  && ptB>84) ||
	     (HLT_PFJet40  && ptB>37));
	  bool trigA =
	    ((HLT_PFJet500 && pt>638) ||
	     (HLT_PFJet450 && pt>548) ||
	     (HLT_PFJet400 && pt>468) ||
	     (HLT_PFJet320 && pt>395) ||
	     (HLT_PFJet260 && pt>300) ||
	     (HLT_PFJet200 && pt>245) ||
	     (HLT_PFJet140 && pt>196) ||
	     (HLT_PFJet80  && pt>114) ||
	     (HLT_PFJet60  && pt>84) ||
	     (HLT_PFJet40  && pt>37));

	  // Look at good probe jets in barrel in good events
	  if (fabs(eta)<1.3 && fabs(etaB)<1.3 &&
	      idtightA && idtightB &&
	      flag_tA && //flag_tB &&
	      pt > 0.5*ptB && ptB > 0.5*pt
	    //dr < 0.10) {
	    ) {

	    // Tag-and-probe method
	    if (istp) {
	      pta_tp->Fill(pttag, pt);
	      pa_tp->Fill(pttag, pt / pttag);
	      pb_tp->Fill(pttag, ptB / pttag);
	      pd_tp->Fill(pttag, 0.5*(ptB-pt) / pttag);
	      
	      h2a_tp->Fill(pttag, pt / pttag);
	      h2b_tp->Fill(pttag, ptB / pttag);
	      h2d_tp->Fill(pttag, 0.5*(ptB-pt) / pttag);
	    } //istp

	    // Direct matching method
	    if (trigA) pjesa_dm->Fill(pt, jes);
	    if (trigB) pjesb_dm->Fill(ptB, jesB);
	      
	    if (trigA) pta_dm->Fill(pt, pt);
	    if (trigA) pa_dm->Fill(pt, ptB / pt);
	    if (trigB) ptb_dm->Fill(ptB, pt);
	    if (trigB) pb_dm->Fill(ptB, pt / ptB);
	    if (trigD) ptd_dm->Fill(ptave, pt);
	    if (trigD) pd_dm->Fill(ptave, 0.5*(ptB-pt) / ptave);
	    
	    if (trigA) h2a_dm->Fill(pt, ptB / pt);
	    if (trigB) h2b_dm->Fill(ptB, pt / ptB);
	    if (trigD) h2d_dm->Fill(ptave, 0.5*(ptB-pt) / ptave);
	  } // good barrel probe
	} // dr match
      } // for j
    } // for i
  } // for events
  cout << endl << "Found " << nev << " matching events"// << endl;
    //<< " of which " << nj << " had same number of jets" << endl;
       << " of which " << ngood << " passed JSON selection" << endl;
  cout << "Found " << nmatch << " matching jets in these events" << endl;
    
  cout << "Output stored to " << f->GetName() << endl;
  f->Write();
  f->Close();

  t.Stop();
  cout << "Processing used " << t.CpuTime() << "s CPU time ("
       << t.CpuTime()/3600. << "h)" << endl;
  cout << "Processing used " << t.RealTime() << "s real time ("
       << t.RealTime()/3600. << "h)" << endl;
  cout << endl << endl;

}
