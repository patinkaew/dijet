#define DijetHistosFill_cxx
#include "DijetHistosFill.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <string>

// Fill multijet histograms
bool doJetveto = true;
bool doIncjet = true;
bool doDijet = true;
bool doDijetOrig = false;
bool doMultijet = true;
bool debug = false; // general debug
bool debugevent = false; // per-event debug

// Maximum asymmetry of 2/3 corresponds to x2 ratio of tag and probe
// Permit ~0.7 extra scaling to allow for HF L3Res
const double maxa = 10; // no cut with 10


//double DELTAPHI(double a, double b) {
//double phi1 = max(a,b);
//double phi2 = min(a,b);
//double d = phi1-phi2;
//if (d>TMath::Pi()) d -= TMath::TwoPi();
//return fabs(d);
//}
double DELTAPHI(double phi1, double phi2) {
  double dphi = fabs(phi1-phi2);
  return (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}
double DELTAR(double phi1, double phi2, double eta1, double eta2) {
  return sqrt(pow(DELTAPHI(phi1,phi2),2) + pow(eta1-eta2,2));
}

// Hardcoded pT, eta thresholds for each trigger
// used in e.g. jetvetoHistos
struct range {
  double ptmin;
  double ptmax;
  double absetamin;
  double absetamax;
};
std::map<std::string, struct range> mt;

class jetvetoHistos {
public:

  // Basic information about the trigger
  string trg;
  int trgpt;
  double ptmin, ptmax, absetamin, absetamax;
  
  TH1D *hpt;       // jet/event counts without veto
  TH1D *heta;      // jet/event counts without veto 
  TH1D *hpt_veto;  // jet/event counts with veto
  TH1D *heta_veto; // jet/event counts with veto
  TH1D *heta_pretrg; // jet/event counts without veto and trigger eta

  TH2D *h2etaphi; // jet counts
  TProfile2D *p2asymm;     // balancing

  // optional composition plots
  TProfile2D *p2chf, *p2nef, *p2nhf;
  TProfile2D *p2chftp, *p2neftp, *p2nhftp;
};

class incjetHistos {
public:

  // Basic information about the trigger
  string trg;
  int trgpt;
  double ptmin, ptmax, absetamin, absetamax;

  static const int ny = 10;
  TH1D *hall;
  TH1D *hsel;
  TH1D *hpt13;
  TH1D* vpt[ny];
  TH2D *h2pt;
};

class dijetHistos {
public:

  // Basic information about the trigger
  string trg;
  int trgpt;
  double ptmin, ptmax, absetamin, absetamax;

  TH2D *h2all, *h2sel;
  TProfile2D *p2jes2; // JEC L2L3Res for undoing
  TProfile2D *p2m0, *p2m0x, *p2m2, *p2m2x; // JER MPFX, DBX methods
  TProfile2D *p2m0ab, *p2m2ab, *p2mnab, *p2muab; // pT,ave (bisector)
  TProfile2D *p2m0tc, *p2m2tc, *p2mntc, *p2mutc; // pT,tag (central)
  TProfile2D *p2m0pf, *p2m2pf, *p2mnpf, *p2mupf; // pt,probe (forward)
};

class multijetHistos {
public:

  // Basic information about the trigger
  string trg;
  int trgpt;

  TH1D *hna, *hnl, *hnr;
  TProfile *presa, *presl, *presr;
  TProfile *pm0a, *pm2a, *pmna, *pmua, *pmoa;
  TProfile *pm0l, *pm2l, *pmnl, *pmul, *pmol;
  TProfile *pm0r, *pm2r, *pmnr, *pmur, *pmor;

  // controls
  TH2D *h2m0a;
  TH2D *h2m2a;
  TH1D *hcosdphi;
};

class dijetHistosOrig {
public:

  // Basic information about the trigger
  string trg;
  int trgpt;
  double ptmin;
  double ptmax;

  // Control pT distributions from the trigger
  TH1D *hpta;  // pTave
  TH1D *hptt;  // pTtag
  TH1D *hptp;  // pTprobe
  TH1D *hptaf; // pTave, forward eta (|eta|>2.8)
  TH1D *hpttf; // pTtag, forward eta (|eta|>2.8)
  TH1D *hptpf; // pTprobe, forward eta (|eta|>2.8)

  // Number of events vs eta for e.g. determining Fwd, HFJEC thresholds
  TH1D *hna;
  TH1D *hnt;
  TH1D *hnp;

  // Profiles for dijet balance vs eta using various schemes:
  // ab==pTave, bisector axis
  TProfile *pm0ab;
  TProfile *pm2ab;
  TProfile *pmnab;
  TProfile *pmuab;
  TProfile *pmoab;

  // tb==pTtag, bisector axis
  TProfile *pm0tb;
  TProfile *pm2tb;
  TProfile *pmntb;
  TProfile *pmutb;
  TProfile *pmotb;

  // pb==pTprobe, bisector axis
  TProfile *pm0pb;
  TProfile *pm2pb;
  TProfile *pmnpb;
  TProfile *pmupb;
  TProfile *pmopb;

  // ac==pTave, central jet axis
  TProfile *pm0ac;
  TProfile *pm2ac;
  TProfile *pmnac;
  TProfile *pmuac;
  TProfile *pmoac;

  // tc==pTtag, central jet axis
  TProfile *pm0tc;
  TProfile *pm2tc;
  TProfile *pmntc;
  TProfile *pmutc;
  TProfile *pmotc;

  // af==pTave, forward jet axis
  TProfile *pm0af;
  TProfile *pm2af;
  TProfile *pmnaf;
  TProfile *pmuaf;
  TProfile *pmoaf;

  TProfile *pm0pf;
  TProfile *pm2pf;
  TProfile *pmnpf;
  TProfile *pmupf;
  TProfile *pmopf;

  // 2D distributions for MPF-MPFX method
  // binning as in above TProfiles (ab, tb, pb, tc, pf)
  // m0=sum from all jets, m2=sum from 2 leading jets only
  TH2D *h2m0ab;
  TH2D *h2m0abx;
  TH2D *h2m0tb;
  TH2D *h2m0tbx;
  TH2D *h2m0pb;
  TH2D *h2m0pbx;
  //
  TH2D *h2m2ab;
  TH2D *h2m2abx;
  TH2D *h2m2tb;
  TH2D *h2m2tbx;
  TH2D *h2m2pb;
  TH2D *h2m2pbx;

  TH2D *h2m0tc;
  TH2D *h2m0tcx;
  TH2D *h2m0pf;
  TH2D *h2m0pfx;
   //
  TH2D *h2m2tc;
  TH2D *h2m2tcx;
  TH2D *h2m2pf;
  TH2D *h2m2pfx;

  // 2D distributions for jet veto maps (pTave, pTtag, pTprobe bins)
  TProfile2D *p2etaphia, *p2etaphit, *p2etaphip;

  // Mass profiles for 
  //TH1D *hmjji, *hmjji13, *hmjj1, *hmjj113;
  TH1D *hmjj2, *hmjj213;
  //TH1D *hmjj, *hmjj113;

  // PF composition plots (+JES)
  TProfile *pjes, *pres, *pdjes, *pjes13, *pres13, *pdjes13;
  TProfile *pchf, *pnhf, *pnef, *pcef, *pmuf, *phef, *phhf;
  TProfile *pchf13, *pnhf13, *pnef13, *pcef13, *pmuf13;

  // HLT-offline trigger matching plots for HLT JES
  TProfile *phashlt, *pishlt, *pismax, *pisnear, *pisclose;
  TProfile *pjeshlt, *pjeshltmax, *pjeshltnear, *pjeshltclose;
  TProfile *ppthlt, *pptoff, *ppttag;
  TH2D *h2jeshlt;
};

// Helper function to retrieve FactorizedJetCorrector 
FactorizedJetCorrector *getFJC(string l1="", string l2="", string res="",
			       string path="") {

  // Set default jet algo
  if (l1!="" && !(TString(l1.c_str()).Contains("_AK")))
    l1 += "_AK4PFPuppi";
  if (l2!="" && !(TString(l2.c_str()).Contains("_AK")))
    l2 += "_AK4PFPuppi";
  if (res!="" && !(TString(res.c_str()).Contains("_AK")))
    res += "_AK4PFPuppi";

  // Set default path
  if (path=="") path = "CondFormats/JetMETObjects/data";
  const char *cd = path.c_str();
  const char *cl1 = l1.c_str();
  const char *cl2 = l2.c_str();
  const char *cres = res.c_str();
  string s("");

  vector<JetCorrectorParameters> v;
  if (l1!=""){
    s = Form("%s/%s.txt",cd,cl1);
    cout << s << endl << flush;
    JetCorrectorParameters *pl1 = new JetCorrectorParameters(s);
    v.push_back(*pl1);
  }
  if (l2!="") {
    s = Form("%s/%s.txt",cd,cl2);
    cout << s << endl << flush;
    JetCorrectorParameters *pl2 = new JetCorrectorParameters(s);
    v.push_back(*pl2);
  }
  if (res!="") {
    s = Form("%s/%s.txt",cd,cres);
    cout << s << endl << flush;
    JetCorrectorParameters *pres = new JetCorrectorParameters(s);
    v.push_back(*pres);
  }
  FactorizedJetCorrector *jec = new FactorizedJetCorrector(v);

  return jec;
} // getJF

void DijetHistosFill::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DijetHistosFill.C
//      root> DijetHistosFill t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);

   if (debug) cout << "Setting branch status for "
		   << (isMC ? "MC" : "DATA") << endl << flush;
   
   if (isMC) fChain->SetBranchStatus("genWeight",1); // v2

   fChain->SetBranchStatus("run",1);
   fChain->SetBranchStatus("luminosityBlock",1);
   fChain->SetBranchStatus("event",1);
   //fChain->SetBranchStatus("Rho_fixedGridRhoAll",1);
   if (isRun2) fChain->SetBranchStatus("fixedGridRhoFastjetAll",1);
   
   //fChain->SetBranchStatus("HLT_DiPFJetAve40",1);
   //fChain->SetBranchStatus("HLT_PFJet40",1); // no events
   //fChain->SetBranchStatus("HLT_AK8PFJet40",1);
   //fChain->SetBranchStatus("HLT_PFJetFwd40",1);
   //fChain->SetBranchStatus("HLT_AK8PFJetFwd40",1);

   // Listing of available triggers
   vector<string> vtrg;

   vtrg.push_back("HLT_ZeroBias");
   
   vtrg.push_back("HLT_PFJet40");
   vtrg.push_back("HLT_PFJet60");
   vtrg.push_back("HLT_PFJet80");
   //vtrg.push_back("HLT_PFJet110");
   vtrg.push_back("HLT_PFJet140");
   vtrg.push_back("HLT_PFJet200");
   vtrg.push_back("HLT_PFJet260");
   vtrg.push_back("HLT_PFJet320");
   vtrg.push_back("HLT_PFJet400"); // v14
   vtrg.push_back("HLT_PFJet450");
   vtrg.push_back("HLT_PFJet500");
   if (isRun2>2) vtrg.push_back("HLT_PFJet550");

   vtrg.push_back("HLT_DiPFJetAve40");
   vtrg.push_back("HLT_DiPFJetAve60");
   vtrg.push_back("HLT_DiPFJetAve80");
   vtrg.push_back("HLT_DiPFJetAve140");
   vtrg.push_back("HLT_DiPFJetAve200");
   vtrg.push_back("HLT_DiPFJetAve260");
   vtrg.push_back("HLT_DiPFJetAve320");
   vtrg.push_back("HLT_DiPFJetAve400");
   vtrg.push_back("HLT_DiPFJetAve500");

   vtrg.push_back("HLT_DiPFJetAve60_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve80_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve100_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve160_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve220_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve300_HFJEC");

   //vtrg.push_back("HLT_PFJetFwd15");
   //vtrg.push_back("HLT_PFJetFwd25");
   if (isRun2>2) {
     vtrg.push_back("HLT_PFJetFwd40");
     vtrg.push_back("HLT_PFJetFwd60");
     vtrg.push_back("HLT_PFJetFwd80");
     vtrg.push_back("HLT_PFJetFwd140");
     vtrg.push_back("HLT_PFJetFwd200");
     vtrg.push_back("HLT_PFJetFwd260");
     vtrg.push_back("HLT_PFJetFwd320");
     vtrg.push_back("HLT_PFJetFwd400");
     vtrg.push_back("HLT_PFJetFwd450");
     vtrg.push_back("HLT_PFJetFwd500");
   }
   int ntrg = vtrg.size();

   for (int i = 0; i != ntrg; ++i) {
     fChain->SetBranchStatus(vtrg[i].c_str(),1);
     if (mtrg[vtrg[i]]==0) {
       cout << "Missing branch info for " << vtrg[i] << endl << flush;
     }
     assert(mtrg[vtrg[i]]!=0);
   }

   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_mass",1);

   fChain->SetBranchStatus("Jet_rawFactor",1);
   if (isRun2) fChain->SetBranchStatus("Jet_area",1);
   
   bool doPFComposition = true;
   if (doPFComposition) {
     fChain->SetBranchStatus("Jet_chHEF",1);  // h+
     fChain->SetBranchStatus("Jet_neHEF",1);  // h0
     fChain->SetBranchStatus("Jet_neEmEF",1); // gamma
     fChain->SetBranchStatus("Jet_chEmEF",1); // e
     fChain->SetBranchStatus("Jet_muEF",1);   // mu
     fChain->SetBranchStatus("Jet_hfEmEF",1); // HFe
     fChain->SetBranchStatus("Jet_hfHEF",1);  // HFh
   }

   if (isRun2) {
     fChain->SetBranchStatus("MET_pt",1);
     fChain->SetBranchStatus("MET_phi",1);
   }
   else {
     fChain->SetBranchStatus("PuppiMET_pt",1);
     fChain->SetBranchStatus("PuppiMET_phi",1);
   }

   fChain->SetBranchStatus("Flag_METFilters",1);

   // Trigger studies => TrigObjAK4 later (fixed now)
   bool doTriggerMatch = false;
   nTrigObjJMEAK4 = 0; // turn off
   if (doTriggerMatch) {
     /*
     fChain->SetBranchStatus("nTrigObj",1);
     fChain->SetBranchStatus("TrigObj_pt",1);
     fChain->SetBranchStatus("TrigObj_eta",1);
     fChain->SetBranchStatus("TrigObj_phi",1);
     fChain->SetBranchStatus("TrigObj_id",1); // Jet==1, FatJet==6?
     */
     // https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_8/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L136-L180

     fChain->SetBranchStatus("nTrigObjJMEAK4",1);
     fChain->SetBranchStatus("TrigObjJMEAK4_pt",1);
     fChain->SetBranchStatus("TrigObjJMEAK4_eta",1);
     fChain->SetBranchStatus("TrigObjJMEAK4_phi",1);
   }

   // List reference pT and abseta thresholds for triggers
   mt["HLT_ZeroBias"]  = range{15,  3000,  0, 5.2};
   
   mt["HLT_DiPFJetAve40"]  = range{40,  85,  0, 5.2};
   mt["HLT_DiPFJetAve60"]  = range{85,  100, 0, 5.2};
   mt["HLT_DiPFJetAve80"]  = range{100, 155, 0, 5.2};
   mt["HLT_DiPFJetAve140"] = range{155, 210, 0, 5.2};
   mt["HLT_DiPFJetAve200"] = range{210, 300, 0, 5.2};
   mt["HLT_DiPFJetAve260"] = range{300, 400, 0, 5.2};
   mt["HLT_DiPFJetAve320"] = range{400, 500, 0, 5.2};
   mt["HLT_DiPFJetAve400"] = range{500, 600, 0, 5.2};
   mt["HLT_DiPFJetAve500"] = range{600,3000, 0, 5.2};
   
   //2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
   double fwdeta = 3.139; // was 2.853. 80% (100%) on negative (positive) side
   double fwdeta0 = 2.964;//2.853; // 40 and 260 up
   mt["HLT_DiPFJetAve60_HFJEC"]  = range{85,  100, fwdeta, 5.2};
   mt["HLT_DiPFJetAve80_HFJEC"]  = range{100, 125, fwdeta, 5.2};
   mt["HLT_DiPFJetAve100_HFJEC"] = range{125, 180, fwdeta, 5.2};
   mt["HLT_DiPFJetAve160_HFJEC"] = range{180, 250, fwdeta, 5.2};
   mt["HLT_DiPFJetAve220_HFJEC"] = range{250, 350, fwdeta0, 5.2};
   mt["HLT_DiPFJetAve300_HFJEC"] = range{350,3000, fwdeta0, 5.2};
   
   mt["HLT_PFJet40"]  = range{40,  85,  0, 5.2};
   mt["HLT_PFJet60"]  = range{85,  100, 0, 5.2};
   mt["HLT_PFJet80"]  = range{100, 155, 0, 5.2};
   mt["HLT_PFJet140"] = range{155, 210, 0, 5.2};
   mt["HLT_PFJet200"] = range{210, 300, 0, 5.2};
   mt["HLT_PFJet260"] = range{300, 400, 0, 5.2};
   mt["HLT_PFJet320"] = range{400, 500, 0, 5.2};
   mt["HLT_PFJet400"] = range{500, 600, 0, 5.2};
   mt["HLT_PFJet450"] = range{500, 600, 0, 5.2};
   mt["HLT_PFJet500"] = range{600, 700, 0, 5.2};
   mt["HLT_PFJet550"] = range{700,3000, 0, 5.2};
   
   mt["HLT_PFJetFwd40"] = range{40,  85,  fwdeta0, 5.2};
   mt["HLT_PFJetFwd60"] = range{85,  100, fwdeta, 5.2};
   mt["HLT_PFJetFwd80"] = range{100, 155, fwdeta, 5.2};
   mt["HLT_PFJetFwd140"] = range{155, 210, fwdeta, 5.2};
   mt["HLT_PFJetFwd200"] = range{210, 300, fwdeta0, 5.2};
   mt["HLT_PFJetFwd260"] = range{300, 400, fwdeta0, 5.2};
   mt["HLT_PFJetFwd320"] = range{400, 500, fwdeta0, 5.2};
   mt["HLT_PFJetFwd400"] = range{500, 600, fwdeta0, 5.2};
   mt["HLT_PFJetFwd450"] = range{500, 600, fwdeta0, 5.2}; // x
   mt["HLT_PFJetFwd500"] = range{600,3000, fwdeta0, 5.2};
   
   
   if (debug) cout << "Setting up JEC corrector" << endl << flush;

   // Redo JEC
   FactorizedJetCorrector *jec(0);
   //jec = getFJC("","Winter22Run3_V1_MC_L2Relative","","");
   if (isRun2==0) {
     jec = getFJC("","Winter22Run3_V1_MC_L2Relative",
		  isMC ? "":"Winter22Run3_RunC_V2_DATA_L2L3Residual_AK4PFPuppi");
   }
   if (isRun2==1) {
     exit(0);
   }
   if (isRun2==2) {
     if (isMC) 
       jec = getFJC("Summer19UL16_V7_MC_L1FastJet_AK4PFchs",
		    "Summer19UL16_V7_MC_L2Relative_AK4PFchs","");
     else
       jec = getFJC("Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFchs",
		    "Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs",
		    "Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs");
   }
   if (isRun2==3) {
     exit(0);
   }
   if (isRun2==4) {
     exit(0);
   }
   
   
   TLorentzVector p4met, p4dj;
   //TLorentzVector p4, p4s, p4mht, p4mht2, p4mhtc, p4mhtc3, p4t, p4p;
   TLorentzVector p4, p4s, p4t, p4p;
   TLorentzVector p4lead, p4recoil, p4other;
   TLorentzVector p4leadRES, p4recoilRES;
   TLorentzVector p4b, p4bx, p4c, p4cx, p4f, p4fx, p4l, p4r;
   TLorentzVector p4m0, p4m2, p4mn, p4mu, p4mo;
   TLorentzVector p4m3, p4mo3;
   TFile *fout = new TFile(Form("rootfiles/jmenano_%s_out.root",
				isMC ? "mc" : "data"), "RECREATE");
   
   // Monitor trigger rates
   TH1D *htrg = new TH1D("htrg","Triggers;Trigger;N_{events}",
			 vtrg.size(),0,vtrg.size());
   for (int i = 1; i != htrg->GetNbinsX()+1; ++i) {
     htrg->GetXaxis()->SetBinLabel(i,vtrg[i-1].c_str());
   }

   if (debug) cout << "Setting up histograms" << endl << flush;   

   // Inclusive jets pT binning
   double vpti[] = 
     {1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
      97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
      507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
      1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
      2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
      4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
   double npti = sizeof(vpti)/sizeof(vpti[0])-1;

   // Regular L2Relative and L2Res eta binning
   double vx[] =
     {-5.191,
      -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
      -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
      -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
      -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, 
      -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
      0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
      1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
      2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
      4.363, 4.538, 4.716, 4.889, 5.191};
   const int nx = sizeof(vx)/sizeof(vx[0])-1;
   const int ny = 800;
   double vy[ny+1];
   for (int i = 0; i != ny+1; ++i) vy[i] = -3. + (+5.+3.)/ny*i;
   const int nz = 400;
   double vz[nz+1];
   for (int i = 0; i != nz+1; ++i) vz[i] = -1. + (+3.+1.)/nz*i;
   
   // Extras for JMENANO
   TH2D *h2mhtvsmet = new TH2D("h2mhtvsmet","MHT vs MET;MET;MHT",
			       500,0,500,500,0,500);
   TH2D *h2dphi = new TH2D("h2dphi","#Delta#phi vs #eta;#eta;#Delta#phi",
			   nx,vx,126,-TMath::TwoPi(),+TMath::TwoPi());

   // L2Res profiles for HDM method
   // coding: m0=MPF, m2=DB, mn=n-jet, mu=uncl. (observable)
   //         a=PtAve, t=PtTag, p=PtProbe (binning)
   //         b=bisector, c=central(tag), f=forward(probe) (projection axis)
   fout->mkdir("Refs");
   fout->cd("Refs");
   TH1D *hna = new TH1D("hna","PtAve;#eta;N_{events}",nx,vx);
   TH1D *hnt = new TH1D("hnt","PtTag;#eta;N_{events}",nx,vx);
   TH1D *hnp = new TH1D("hnp","PtProbe;#eta;N_{events}",nx,vx);

   TProfile *pm0ab = new TProfile("pm0ab","PtAve Bisector;#eta;MPF",nx,vx);
   TProfile *pm2ab = new TProfile("pm2ab","PtAve Bisector;#eta;MPF2",nx,vx);
   TProfile *pmnab = new TProfile("pmnab","PtAve Bisector;#eta;MPFn",nx,vx);
   TProfile *pmuab = new TProfile("pmuab","PtAve Bisector;#eta;MPFu",nx,vx);
   TProfile *pmoab = new TProfile("pmoab","PtAve Bisector;#eta;MPFo",nx,vx);

   TProfile *pm0tb = new TProfile("pm0tb","PtTag Bisector;#eta;MPF",nx,vx);
   TProfile *pm2tb = new TProfile("pm2tb","PtTag Bisector;#eta;MPF2",nx,vx);
   TProfile *pmntb = new TProfile("pmntb","PtTag Bisector;#eta;MPFn",nx,vx);
   TProfile *pmutb = new TProfile("pmutb","PtTag Bisector;#eta;MPFu",nx,vx);
   TProfile *pmotb = new TProfile("pmotb","PtTag Bisector;#eta;MPFo",nx,vx);

   TProfile *pm0pb = new TProfile("pm0pb","PtProbe Bisector;#eta;MPF",nx,vx);
   TProfile *pm2pb = new TProfile("pm2pb","PtProbe Bisector;#eta;MPF2",nx,vx);
   TProfile *pmnpb = new TProfile("pmnpb","PtProbe Bisector;#eta;MPFn",nx,vx);
   TProfile *pmupb = new TProfile("pmupb","PtProbe Bisector;#eta;MPFu",nx,vx);
   TProfile *pmopb = new TProfile("pmopb","PtProbe Bisector;#eta;MPFo",nx,vx);
 
   TProfile *pm0ac = new TProfile("pm0ac","PtAve Central;#eta;MPF",nx,vx);
   TProfile *pm2ac = new TProfile("pm2ac","PtAve Central;#eta;MPF2",nx,vx);
   TProfile *pmnac = new TProfile("pmnac","PtAve Central;#eta;MPFn",nx,vx);
   TProfile *pmuac = new TProfile("pmuac","PtAve Central;#eta;MPFu",nx,vx);
   TProfile *pmoac = new TProfile("pmoac","PtAve Central;#eta;MPFo",nx,vx);

   TProfile *pm0tc = new TProfile("pm0tc","PtTag Central;#eta;MPF",nx,vx);
   TProfile *pm2tc = new TProfile("pm2tc","PtTag Central;#eta;MPF2",nx,vx);
   TProfile *pmntc = new TProfile("pmntc","PtTag Central;#eta;MPFn",nx,vx);
   TProfile *pmutc = new TProfile("pmutc","PtTag Central;#eta;MPFu",nx,vx);
   TProfile *pmotc = new TProfile("pmotc","PtTag Central;#eta;MPFo",nx,vx);

   TProfile *pm0af = new TProfile("pm0af","PtAve Forward;#eta;MPF",nx,vx);
   TProfile *pm2af = new TProfile("pm2af","PtAve Forward;#eta;MPF2",nx,vx);
   TProfile *pmnaf = new TProfile("pmnaf","PtAve Forward;#eta;MPFn",nx,vx);
   TProfile *pmuaf = new TProfile("pmuaf","PtAve Forward;#eta;MPFu",nx,vx);
   TProfile *pmoaf = new TProfile("pmoaf","PtAve Forward;#eta;MPFo",nx,vx);

   TProfile *pm0pf = new TProfile("pm0pf","PtProbe Forward;#eta;MPF",nx,vx);
   TProfile *pm2pf = new TProfile("pm2pf","PtProbe Forward;#eta;MPF2",nx,vx);
   TProfile *pmnpf = new TProfile("pmnpf","PtProbe Forward;#eta;MPFn",nx,vx);
   TProfile *pmupf = new TProfile("pmupf","PtProbe Forward;#eta;MPFu",nx,vx);
   TProfile *pmopf = new TProfile("pmopf","PtProbe Forward;#eta;MPFo",nx,vx);

   // 2D distributions for MPF-MPFX method
   TH2D *h2m0ab = new TH2D("h2m0ab","PtAve Bisector;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0abx = new TH2D("h2m0abx","PtAve Bisector;#eta;MPFX",nx,vx,ny,vy);
   TH2D *h2m0tb = new TH2D("h2m0tb","PtTag Bisector;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0tbx = new TH2D("h2m0tbx","PtTag Bisector;#eta;MPFX",nx,vx,ny,vy);
   TH2D *h2m0pb  = new TH2D("h2m0pb","PtProbe Bisector;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0pbx = new TH2D("h2m0pbx","PtProbe Bisector;#eta;MPFX",nx,vx,ny,vy);
   //
   TH2D *h2m2ab = new TH2D("h2m2ab","PtAve Bisector;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2abx = new TH2D("h2m2abx","PtAve Bisector;#eta;MPF2X",nx,vx,nz,vz);
   TH2D *h2m2tb = new TH2D("h2m2tb","PtTag Bisector;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2tbx = new TH2D("h2m2tbx","PtTag Bisector;#eta;MPF2X",nx,vx,nz,vz);
   TH2D *h2m2pb = new TH2D("h2m2pb","PtProb Bisector;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2pbx =new TH2D("h2m2pbx","PtProb Bisector;#eta;MPF2X",nx,vx,nz,vz);

   TH2D *h2m0tc = new TH2D("h2m0tc","PtTag Central;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0tcx = new TH2D("h2m0tcx","PtTag Central;#eta;MPFX",nx,vx,ny,vy);
   TH2D *h2m0pf  = new TH2D("h2m0pf","PtProbe Forward;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0pfx = new TH2D("h2m0pfx","PtProbe Forward;#eta;MPFX",nx,vx,ny,vy);
   //
   TH2D *h2m2tc = new TH2D("h2m2tc","PtTag Central;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2tcx = new TH2D("h2m2tcx","PtTag Central;#eta;MPF2X",nx,vx,nz,vz);
   TH2D *h2m2pf = new TH2D("h2m2pf","PtProb Forward;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2pfx =new TH2D("h2m2pfx","PtProb Forward;#eta;MPF2X",nx,vx,nz,vz);

   // HLT-offline matching using tag-and-probe (trigger-matched tag)
   TProfile *phashlt = new TProfile("phashlt",";#eta;Has HLT match",nx,vx);
   TProfile *pishlt = new TProfile("pishlt",";#eta;Is good HLT match",nx,vx);
   TProfile *pismax = new TProfile("pismax",";#eta;Hardest HLT match",nx,vx);
   TProfile *pisnear = new TProfile("pisnear",";#eta;Nearest HLT match",nx,vx);
   TProfile *pisclose = new TProfile("pisclose",";#eta;Closest in pT HLT match",nx,vx);
   TProfile *pjeshlt = new TProfile("pjeshlt",";#eta;HLT/Off",nx,vx);
   TProfile *pjeshltmax = new TProfile("pjeshltmax",";#eta;HLT/Off",nx,vx);
   TProfile *pjeshltnear = new TProfile("pjeshltnear",";#eta;HLT/Off",nx,vx);
   TProfile *pjeshltclose = new TProfile("pjeshltclose",";#eta;HLT/Off",nx,vx);
   TProfile *ppthlt = new TProfile("ppthlt",";#eta;PtProbe(HLT)",nx,vx);
   TProfile *pptoff = new TProfile("pptoff",";#eta;PtProbe(off)",nx,vx);
   TProfile *ppttag = new TProfile("ppttag",";#eta;PtTag",nx,vx);
   TH2D *h2jeshlt = new TH2D("h2jeshlt",";#eta;HLT/Off",nx,vx,ny,vy);
   
   // 2D eta-phi histograms and profiles for jet veto maps
   TH2D *h2etaphi = new TH2D("h2etaphi","All;#eta;#phi",nx,vx,
			     72,-TMath::Pi(),TMath::Pi());
   TH2D *h2etaphi40 = new TH2D("h2etaphi40","All;#eta;#phi",nx,vx,
			       72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphia = new TProfile2D("p2etaphia",
					 "Tag-and-probe PtAve;#eta;#phi;MPF2",
					  nx,vx,
					  72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphit = new TProfile2D("p2etaphit",
					 "Tag-and-probe PtTag;#eta;#phi;MPF2",
					  nx,vx,
					  72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphip = new TProfile2D("p2etaphip",
					 "Tag-and-probe PtProbe;#eta;#phi;MPF2",
					  nx,vx,
					  72,-TMath::Pi(),TMath::Pi());
   
   // Controls
   TH1D *hpt = new TH1D("hpt","hpt",500,0,500);
   TH1D *hpt30 = new TH1D("hpt30","hpt30",500,0,500);
   TH1D *heta = new TH1D("heta","heta",104,-5.2,5.2);
   TH1D *heta40 = new TH1D("heta40","heta40",104,-5.2,5.2);
   TH1D *hnjet = new TH1D("hnjet","hnjet",500,0,500);
   //TH1D *hnjet30 = new TH1D("hnjet30","hnjet30",500,0,500);

   TH1D *hpta = new TH1D("hpta","PtAve",500,0,500);
   TH1D *hptt = new TH1D("hptt","PtTag",500,0,500);
   TH1D *hptp = new TH1D("hptp","PtProbe",500,0,500);

   const double twopi = TMath::TwoPi();
   TH1D *hdphi = new TH1D("hdphi","#Delta#phi",126,-twopi,+twopi);
   TH1D *hdphib = new TH1D("hdphib","#Delta#phi Barrel",126,-twopi,+twopi);
   TH1D *hdr = new TH1D("hdr","#DeltaR",220,0,11);//126,0,twopi);
   TH1D *hdrb = new TH1D("hdrb","#DeltaR Barrel",220,0,11);//126,0,twopi);
   TH1D *hdri = new TH1D("hdri","#DeltaR Inclusive",220,0,11);//126,0,twopi);
   TH1D *hdrib = new TH1D("hdrib","#DeltaR Inclusive Barrel",220,0,11);
   TH1D *hdribb = new TH1D("hdribb","#DeltaR Inclusive Barrel^{2}",220,0,11);
   TH1D *hdri2 = new TH1D("hdri2","#DeltaR Leads",220,0,11);//126,0,twopi);
   TH1D *hdri2b = new TH1D("hdri2b","#DeltaR Leads Barrel",220,0,11);
   TH1D *hdri2bb = new TH1D("hdri2bb","#DeltaR Leads Barrel^{2}",220,0,11);
   
   TH1D *hav = new TH1D("hav","Asymmetry",400,-2,2);
   TH1D *hat = new TH1D("hat","Pt,probe / Pt,tag",400,0,4);
   TH1D *hap = new TH1D("hap","Pt,tag / Pt,probe",400,0,4);

   TH1D *havb = new TH1D("havb","Asymmetry Barrel",400,-2,2);
   TH1D *hatb = new TH1D("hatb","Pt,probe / Pt,tag Barrel",400,0,4);
   TH1D *hapb = new TH1D("hapb","Pt,tag / Pt,probe Barrel",400,0,4);
   
   //TH1D *hmv = new TH1D("hmv","MPF PtAve",800,-2,6);
   //TH1D *hmt = new TH1D("hmt","MPF PtTag",800,-2,6);
   //TH1D *hmp = new TH1D("hmp","MPF PtProbe",800,-2,6);

   //TH1D *hmva = new TH1D("hmva","MPFA PtAve",800,-2,6);
   //TH1D *hmvc = new TH1D("hmvc","MPFC PtAve",800,-2,6);
   //TH1D *hmvc3 = new TH1D("hmvc3","MPFC3 PtAve",800,-2,6);
   //TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",800,-2,6);
   //TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",400,-2,2);

   // PF composition plots
   // Copy L2Res histograms for multiple pT bins
   const double vpt[] = {40, 50, 60, 70, 85, 100, 125, 155, 180, 210, 250, 300,
			 350, 400, 500, 600, 800, 1000, 1200, 1500,
			 1800, 2100, 2400, 2700, 3000};
   const int npt = sizeof(vpt)/sizeof(vpt[0])-1;
   TH1D *hptbins = new TH1D("hptbins",";p_{T} (GeV);N_{events}",npt,vpt);

   TProfile *pjes = new TProfile("pjes","JES",nx,vx);
   TProfile *pres = new TProfile("pres","RES",nx,vx);
   TProfile *pdjes = new TProfile("pdjes","#DELTAJES",nx,vx);
   TProfile *pjes13 = new TProfile("pjes13","JES",npt,vpt);
   TProfile *pres13 = new TProfile("pres13","RES",npt,vpt);
   TProfile *pdjes13 = new TProfile("pdjes13","#DELTAJES",npt,vpt);

   TProfile *pchf = new TProfile("pchf","CHF",nx,vx);
   TProfile *pnhf = new TProfile("pnhf","NHF",nx,vx);
   TProfile *pnef = new TProfile("pnef","NEF",nx,vx);
   TProfile *pcef = new TProfile("pcef","CEF",nx,vx);
   TProfile *pmuf = new TProfile("pmuf","MUF",nx,vx);
   TProfile *phef = new TProfile("phfhf","HEF",nx,vx);
   TProfile *phhf = new TProfile("phfef","HHF",nx,vx);

   TProfile *pchf13 = new TProfile("pchf13","CHF",npt,vpt);
   TProfile *pnhf13 = new TProfile("pnhf13","NHF",npt,vpt);
   TProfile *pnef13 = new TProfile("pnef13","NEF",npt,vpt);
   TProfile *pcef13 = new TProfile("pcef13","CEF",npt,vpt);
   TProfile *pmuf13 = new TProfile("pmuf13","MUF",npt,vpt);
   
   TH1D *hmjj2 = new TH1D("hmjj2","Dijet mass, 2-lead",1000,0,1000);
   TH1D *hmjj213 = new TH1D("hmjj213","Dijet mass, |#DeltaEta|<1.3, 2-lead",
   			    1000,0,1000);
   fout->cd();
   
   if (debug) cout << "Load JSON (or not)" << endl << flush;
   
   bool doJSON = (true && !isMC);
   if (doJSON) {
     if (!LoadJSON()) {
       cout << "Issues loading the JSON file; aborting..." << endl << flush;
       exit(0);
     }
   }
   int _nbadevts_json(0);

   bool doTrigger = (true && !isMC);
   int _nbadevts_trg(0);
   int _nbadevts_fwdtrg(0);
   int _ngoodevts(0);
   
   // Listing of runs and LS
   int nrun(0), nls(0), nevt(0);
   map<int, map<int, int> > mrunls;

   if (debug) cout << "Setup pT bins (triggers)" << endl << flush;
   
   //vector<basicHistos*> vh(npt);
   map<string, jetvetoHistos*> mhjv;
   map<string, incjetHistos*> mhij;
   map<string, dijetHistos*> mhdj;
   map<string, multijetHistos*> mhmj;
   map<string, vector<dijetHistosOrig*> > mhdjo;
   const bool doPtBins = true;
   if (doPtBins) {

     for (int itrg = 0; itrg != ntrg; ++itrg) {

       if (debug) cout << "Trigger " << vtrg[itrg] << endl << flush;
       
       //vector<dijetHistos*> vh(npt);
       fout->mkdir(vtrg[itrg].c_str());
       fout->cd(vtrg[itrg].c_str());
       TDirectory *dout = gDirectory;
       vector<dijetHistosOrig*> &vh = mhdjo[vtrg[itrg]];
       vh.resize(npt);

       // Figure out trigger pT threshold from the name
       int trgpt(-1), nfound(0);
       if (nfound!=1) {
	 nfound = (vtrg[itrg]=="HLT_ZeroBias" ? 1 : 0);
	 trgpt = 0;
       }
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_PFJet%d",&trgpt);
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_PFJetFwd%d",&trgpt);
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_DiPFJetAve%d",&trgpt);
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_DiPFJetAve%d_HFJEC",&trgpt);
       if (!(nfound==1 && trgpt!=-1)) {
	 cout << "trigger " << vtrg[itrg] << ": nfound="<<nfound
	      << ", trgpt = " << trgpt << endl << flush;
       }
       assert(nfound==1);
       assert(trgpt!=-1);

       // Jet veto per trigger
       if (doJetveto) {
	 if (debug) cout << "Setup doJetveto " << trgpt << endl << flush;
	 
	 dout->mkdir("Jetveto");
	 dout->cd("Jetveto");

	 jetvetoHistos *h = new jetvetoHistos();

	 string &t = vtrg[itrg];
	 mhjv[t] = h;
	 h->trg = t;
	 h->trgpt = trgpt;

	 struct range &r  = mt[t];
	 h->ptmin = r.ptmin;
	 h->ptmax = r.ptmax;
	 h->absetamin = r.absetamin;
	 h->absetamax = r.absetamax;

	 h->hpt = new TH1D("hpt",";p_{T,jet} (GeV);N_{jet} (no veto)",
			   npt,vpt);
	 h->hpt_veto = new TH1D("hpt_veto",";p_{T,jet} (GeV);N_{jet} (veto)",
				npt,vpt);
	 h->heta = new TH1D("heta",";p_{T,jet} (GeV);N_{jet} (no veto)",
			    nx,vx);
	 h->heta_veto = new TH1D("heta_veto",";p_{T,jet} (GeV);N_{jet} (veto)",
				 nx,vx);
	 h->heta_pretrg = new TH1D("heta_pretrg",";p_{T,jet} (GeV);"
				   "N_{jet} (veto)",nx,vx);

	 h->h2etaphi = new TH2D("h2pt",";#eta;#phi;N_{jet}",
				nx,vx, 72,-TMath::Pi(),+TMath::Pi());
	 h->p2asymm = new TProfile2D("p2asymm",";#eta;#phi;Asymmetry",
				     nx,vx, 72,-TMath::Pi(),+TMath::Pi());

	 if (doPFComposition) {
	   h->p2chf = new TProfile2D("p2chf",";#eta;#phi;CHF (DM)",
				     nx,vx, 72,-TMath::Pi(),+TMath::Pi());
	   h->p2nef = new TProfile2D("p2nef",";#eta;#phi;NEF (DM)",
				     nx,vx, 72,-TMath::Pi(),+TMath::Pi());
	   h->p2nhf = new TProfile2D("p2nhf",";#eta;#phi;NHF (DM)",
				     nx,vx, 72,-TMath::Pi(),+TMath::Pi());
	   
	   h->p2chftp = new TProfile2D("p2chftp",";#eta;#phi;CHF (TP)",
				       nx,vx, 72,-TMath::Pi(),+TMath::Pi());
	   h->p2neftp = new TProfile2D("p2neftp",";#eta;#phi;NEF (TP)",
				       nx,vx, 72,-TMath::Pi(),+TMath::Pi());
	   h->p2nhftp = new TProfile2D("p2nhftp",";#eta;#phi;NHF (TP)",
				       nx,vx, 72,-TMath::Pi(),+TMath::Pi());
	 }
       }

       // Inclusive jet per trigger
       if (doIncjet) {
	 if (debug) cout << "Setup doIncjet " << trgpt << endl << flush;
	 
	 dout->mkdir("Incjet");
	 dout->cd("Incjet");

	 incjetHistos *h = new incjetHistos();

	 string &t = vtrg[itrg];
	 mhij[t] = h;
	 h->trg = t;
	 h->trgpt = trgpt;

	 struct range &r  = mt[t];
	 h->ptmin = r.ptmin;
	 h->ptmax = r.ptmax;
	 h->absetamin = r.absetamin;
	 h->absetamax = r.absetamax;
	 
	 h->hall = new TH1D("hall",";p_{T,jet} (GeV)",npti,vpti);
	 h->hsel = new TH1D("hsel",";p_{T,jet} (GeV)",npti,vpti);
	 h->hpt13 = new TH1D("hpt13",";p_{T,jet} (GeV)",npti,vpti);
	 h->h2pt = new TH2D("h2pt",";#eta;p_{T} (GeV);N_{jet}",
			    nx,vx,npti,vpti);
	 for (int iy = 0; iy != h->ny; ++iy) {
	   h->vpt[iy] = new TH1D(Form("hpt%d",iy),";p_{T} (GeV);N_{jet}",
				 npti,vpti);
	 } // for iy
       } // incjet

       // Dijet per trigger
       if (doDijet) {
	 if (debug) cout << "Setup doDijet " << trgpt << endl << flush;
	 
	 dout->mkdir("Dijet");
	 dout->cd("Dijet");

	 dijetHistos *h = new dijetHistos();

	 string &t = vtrg[itrg];
	 mhdj[t] = h;
	 h->trg = t;
	 h->trgpt = trgpt;

	 struct range &r  = mt[t];
	 h->ptmin = r.ptmin;
	 h->ptmax = r.ptmax;
	 h->absetamin = r.absetamin;
	 h->absetamax = r.absetamax;

	 // Counting of events
	 h->h2all = new TH2D("h2all",";#eta;p_{T,ave} (GeV);N_{events}",
			     nx,vx, npt, vpt);
	 h->h2sel = new TH2D("h2sel",";#eta;p_{T,ave} (GeV);N_{events}",
			     nx,vx, npt, vpt);
	 
	 // JEC L2L3Res for undoing
	 h->p2jes2 = new TProfile2D("p2jes2",";#eta;p_{T,ave} (GeV);"
				    "JES(probe)/JES(tag)",
				    nx,vx, npt, vpt);
	 
	 // Basic profiles with RMS as error ("S") for JER studies
	 h->p2m0 = new TProfile2D("p2m0",";#eta;p_{T,ave} (GeV);MPF0 (MPF)",
				  nx,vx, npt, vpt, "S");
	 h->p2m0x = new TProfile2D("p2m0x",";#eta;p_{T,ave} (GeV);MPFX0 (MPFX)",
				   nx,vx, npt, vpt, "S");
	 h->p2m2 = new TProfile2D("p2m2",";#eta;p_{T,ave} (GeV);MPF2 (DB)",
				  nx,vx, npt, vpt, "S");
	 h->p2m2x = new TProfile2D("p2m2x",";#eta;p_{T,ave} (GeV);MPF2 (DBX)",
				   nx,vx, npt, vpt, "S");

	 // Variants with different binnings and with error on the mean
	 h->p2m0ab = new TProfile2D("p2m0ab",";#eta;p_{T,ave} (GeV);MPF0",
				    nx,vx, npt, vpt);
	 h->p2m2ab = new TProfile2D("p2m2ab",";#eta;p_{T,ave} (GeV);MPF2",
				    nx,vx, npt, vpt);
	 h->p2mnab = new TProfile2D("p2mnab",";#eta;p_{T,ave} (GeV);MPFn",
				    nx,vx, npt, vpt);
	 h->p2muab = new TProfile2D("p2muab",";#eta;p_{T,ave} (GeV);MPFu",
				    nx,vx, npt, vpt);

	 h->p2m0tc = new TProfile2D("p2m0tc",";#eta;p_{T,tag} (GeV);MPF0",
				    nx,vx, npt, vpt);
	 h->p2m2tc = new TProfile2D("p2m2tc",";#eta;p_{T,tag} (GeV);MPF2",
				    nx,vx, npt, vpt);
	 h->p2mntc = new TProfile2D("p2mntc",";#eta;p_{T,tag} (GeV);MPFn",
				    nx,vx, npt, vpt);
	 h->p2mutc = new TProfile2D("p2mutc",";#eta;p_{T,tag} (GeV);MPFu",
				    nx,vx, npt, vpt);

	 h->p2m0pf = new TProfile2D("p2m0pf",";#eta;p_{T,probe} (GeV);MPF0",
				    nx,vx, npt, vpt);
	 h->p2m2pf = new TProfile2D("p2m2pf",";#eta;p_{T,probe} (GeV);MPF2",
				    nx,vx, npt, vpt);
	 h->p2mnpf = new TProfile2D("p2mnpf",";#eta;p_{T,probe} (GeV);MPFn",
				    nx,vx, npt, vpt);
	 h->p2mupf = new TProfile2D("p2mupf",";#eta;p_{T,probe} (GeV);MPFu",
				    nx,vx, npt, vpt);

       } // doDijet

       // Multijet per trigger
       if (doMultijet) {

	 if (debug) cout << "Setup doMultijet " << trgpt << endl << flush;
	 
	 dout->mkdir("Multijet");
	 dout->cd("Multijet");
	 multijetHistos *h = new multijetHistos();
	 mhmj[vtrg[itrg]] = h;
	 h->trg = vtrg[itrg];
	 h->trgpt = trgpt;
	 //h->ptmin = vpt[ipt];
	 //h->ptmax = vpt[ipt+1];

	 h->hna = new TH1D("hna","",npti,vpti);
	 h->hnl = new TH1D("hnl","",npti,vpti);
	 h->hnr = new TH1D("hnr","",npti,vpti);

	 h->presa = new TProfile("presa","",npti,vpti);
	 h->presl = new TProfile("presl","",npti,vpti);
	 h->presr = new TProfile("presr","",npti,vpti);

	 h->pm0a = new TProfile("pm0a","",npti,vpti);
	 h->pm2a = new TProfile("pm2a","",npti,vpti);
	 h->pmna = new TProfile("pmna","",npti,vpti);
	 h->pmua = new TProfile("pmua","",npti,vpti);
	 h->pmoa = new TProfile("pmoa","",npti,vpti);

	 h->pm0l = new TProfile("pm0l","",npti,vpti);
	 h->pm2l = new TProfile("pm2l","",npti,vpti);
	 h->pmnl = new TProfile("pmnl","",npti,vpti);
	 h->pmul = new TProfile("pmul","",npti,vpti);
	 h->pmol = new TProfile("pmol","",npti,vpti);

	 h->pm0r = new TProfile("pm0r","",npti,vpti);
	 h->pm2r = new TProfile("pm2r","",npti,vpti);
	 h->pmnr = new TProfile("pmnr","",npti,vpti);
	 h->pmur = new TProfile("pmur","",npti,vpti);
	 h->pmor = new TProfile("pmor","",npti,vpti);

	 h->h2m0a = new TH2D("h2m0a","",npti,vpti,200,-1,3);
	 h->h2m2a = new TH2D("h2m2a","",npti,vpti,200,-1,3);
	 h->hcosdphi = new TH1D("hcosdphi","",102,-1.01,1.01);
       } // doMultijet


       if (debug) cout << "Setup doDijetOrig pT bins" << endl << flush;
       
       // Dijet in pT bins
       if (doDijetOrig) {
       for (int ipt = 0; ipt != npt; ++ipt) {

	 if (!dout->FindObject("DijetOrig")) dout->mkdir("DijetOrig");
	 dout->mkdir(Form("DijetOrig/Pt_%d_%d",int(vpt[ipt]),int(vpt[ipt+1])));
	 dout->cd(Form("DijetOrig/Pt_%d_%d",int(vpt[ipt]),int(vpt[ipt+1])));
	 dijetHistosOrig *h = new dijetHistosOrig();
	 vh[ipt] = h;
	 h->trg = vtrg[itrg];
	 h->trgpt = trgpt;
	 h->ptmin = vpt[ipt];
	 h->ptmax = vpt[ipt+1];
	 
	 h->hpta = new TH1D("hpta",";p_{T} (GeV);N_{evts}",1000,0,1000);
	 h->hptt = new TH1D("hptt",";p_{T} (GeV);N_{evts}",1000,0,1000);
	 h->hptp = new TH1D("hptp",";p_{T} (GeV);N_{evts}",1000,0,1000);
	 h->hptaf = new TH1D("hptaf",";p_{T} (GeV);N_{evts}",1000,0,1000);
	 h->hpttf = new TH1D("hpttf",";p_{T} (GeV);N_{evts}",1000,0,1000);
	 h->hptpf = new TH1D("hptpf",";p_{T} (GeV);N_{evts}",1000,0,1000);

	 //if (ptave >= 40 && ptave <50) {
	 h->hna = (TH1D*)hna->Clone();
	 h->hnt = (TH1D*)hnt->Clone();
	 h->hnp = (TH1D*)hnp->Clone();

	 h->pm0ab = (TProfile*)pm0ab->Clone();
	 h->pm2ab = (TProfile*)pm2ab->Clone();
	 h->pmnab = (TProfile*)pmnab->Clone();
	 h->pmuab = (TProfile*)pmuab->Clone();
	 h->pmoab = (TProfile*)pmoab->Clone();
	 
	 h->pm0ac = (TProfile*)pm0ac->Clone();
	 h->pm2ac = (TProfile*)pm2ac->Clone();
	 h->pmnac = (TProfile*)pmnac->Clone();
	 h->pmuac = (TProfile*)pmuac->Clone();
	 h->pmoac = (TProfile*)pmoac->Clone();
	 
	 h->pm0af = (TProfile*)pm0af->Clone();
	 h->pm2af = (TProfile*)pm2af->Clone();
	 h->pmnaf = (TProfile*)pmnaf->Clone();
	 h->pmuaf = (TProfile*)pmuaf->Clone();
	 h->pmoaf = (TProfile*)pmoaf->Clone();
	 
	 h->h2m0ab = (TH2D*)h2m0ab->Clone();
	 h->h2m0abx = (TH2D*)h2m0abx->Clone();
	 h->h2m2ab = (TH2D*)h2m2ab->Clone();
	 h->h2m2abx = (TH2D*)h2m2abx->Clone();
	 
	 h->p2etaphia = (TProfile2D*)p2etaphia->Clone();
	 //}
	 //if (pttag >= 40 && pttag <50) {
	 h->pm0tb = (TProfile*)pm0tb->Clone();
	 h->pm2tb = (TProfile*)pm2tb->Clone();
	 h->pmntb = (TProfile*)pmntb->Clone();
	 h->pmutb = (TProfile*)pmutb->Clone();
	 h->pmotb = (TProfile*)pmotb->Clone();
	 
	 h->pm0tc = (TProfile*)pm0tc->Clone();
	 h->pm2tc = (TProfile*)pm2tc->Clone();
	 h->pmntc = (TProfile*)pmntc->Clone();
	 h->pmutc = (TProfile*)pmutc->Clone();
	 h->pmotc = (TProfile*)pmotc->Clone();
	 
	 h->h2m0tb = (TH2D*)h2m0tb->Clone();
	 h->h2m0tbx = (TH2D*)h2m0tbx->Clone();
	 h->h2m2tb = (TH2D*)h2m2tb->Clone();
	 h->h2m2tbx = (TH2D*)h2m2tbx->Clone();
	 
	 h->h2m0tc = (TH2D*)h2m0tc->Clone();
	 h->h2m0tcx = (TH2D*)h2m0tcx->Clone();
	 h->h2m2tc = (TH2D*)h2m2tc->Clone();
	 h->h2m2tcx = (TH2D*)h2m2tcx->Clone();
	 
	 h->p2etaphit = (TProfile2D*)p2etaphit->Clone();
	 //}
	 //if (ptprobe >= 40 && ptprobe <50) {
	 h->pm0pb = (TProfile*)pm0pb->Clone();
	 h->pm2pb = (TProfile*)pm2pb->Clone();
	 h->pmnpb = (TProfile*)pmnpb->Clone();
	 h->pmupb = (TProfile*)pmupb->Clone();
	 h->pmopb = (TProfile*)pmopb->Clone();
	 
	 h->pm0pf = (TProfile*)pm0pf->Clone(); // fix v9
	 h->pm2pf = (TProfile*)pm2pf->Clone();
	 h->pmnpf = (TProfile*)pmnpf->Clone();
	 h->pmupf = (TProfile*)pmupf->Clone();
	 h->pmopf = (TProfile*)pmopf->Clone();
	 
	 h->h2m0pb = (TH2D*)h2m0pb->Clone();
	 h->h2m0pbx = (TH2D*)h2m0pbx->Clone();
	 h->h2m2pb = (TH2D*)h2m2pb->Clone();
	 h->h2m2pbx = (TH2D*)h2m2pbx->Clone();
	 
	 h->h2m0pf = (TH2D*)h2m0pf->Clone();
	 h->h2m0pfx = (TH2D*)h2m0pfx->Clone();
	 h->h2m2pf = (TH2D*)h2m2pf->Clone();
	 h->h2m2pfx = (TH2D*)h2m2pfx->Clone();
	 
	 h->p2etaphip = (TProfile2D*)p2etaphip->Clone();

	 h->hmjj2 = (TH1D*)hmjj2->Clone();
	 h->hmjj213 = (TH1D*)hmjj213->Clone();

	 h->pjes = (TProfile*)pjes->Clone();
	 h->pres = (TProfile*)pres->Clone();
	 h->pdjes = (TProfile*)pdjes->Clone();
	 h->pjes13 = (TProfile*)pjes13->Clone();
	 h->pres13 = (TProfile*)pres13->Clone();
	 h->pdjes13 = (TProfile*)pdjes13->Clone();

	 if (doPFComposition) {
	   h->pchf = (TProfile*)pchf->Clone();
	   h->pnhf = (TProfile*)pnhf->Clone();
	   h->pnef = (TProfile*)pnef->Clone();
	   h->pcef = (TProfile*)pcef->Clone();
	   h->pmuf = (TProfile*)pmuf->Clone();
	   h->phef = (TProfile*)phef->Clone();
	   h->phhf = (TProfile*)phhf->Clone();
	   
	   h->pchf13 = (TProfile*)pchf13->Clone();
	   h->pnhf13 = (TProfile*)pnhf13->Clone();
	   h->pnef13 = (TProfile*)pnef13->Clone();
	   h->pcef13 = (TProfile*)pcef13->Clone();
	   h->pmuf13 = (TProfile*)pmuf13->Clone();
	 }
	   
	 if (doTriggerMatch) {
	   h->phashlt = (TProfile*)phashlt->Clone();
	   h->pishlt = (TProfile*)pishlt->Clone();
	   h->pismax = (TProfile*)pismax->Clone();
	   h->pisnear = (TProfile*)pisnear->Clone();
	   h->pisclose = (TProfile*)pisclose->Clone();
	   h->pjeshlt = (TProfile*)pjeshlt->Clone();
	   h->pjeshltmax = (TProfile*)pjeshltmax->Clone();
	   h->pjeshltnear = (TProfile*)pjeshltnear->Clone();
	   h->pjeshltclose = (TProfile*)pjeshltclose->Clone();
	   h->ppthlt = (TProfile*)ppthlt->Clone();
	   h->pptoff = (TProfile*)pptoff->Clone();
	   h->ppttag = (TProfile*)ppttag->Clone();
	   h->h2jeshlt = (TH2D*)h2jeshlt->Clone();
	 }
       } // for ipt
       } // doDijetOrig
     } // for itrg
   } // doPtBins

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   cout << "Loaded " << nentries << " entries" << endl << flush;

   // For trigger matching studies
   const int kMaxTrigJet = 3;
   Float_t Jet_hltPt[kMaxTrigJet];
   Float_t Jet_hltPtClose[kMaxTrigJet];
   Float_t Jet_hltPtNear[kMaxTrigJet];
   Float_t Jet_hltPtMax[kMaxTrigJet];
   Float_t Jet_RES[nJetMax];
   Float_t Jet_deltaJES[nJetMax];
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      //if (jentry%1000==0) cout << "." << flush;
      if (jentry%100000==0) cout << "." << flush;
      if (jentry%5000000==0) cout << "n="<<jentry<<endl<<flush;

      double w = (isMC ? genWeight : 1.);

      if (doJSON) {

	if (debugevent) cout << "doJSON: Read in branches" << endl << flush;

	//b_branchname->GetEntry(ientry); //read only this branch
	b_run->GetEntry(ientry); //read only this branch
	b_luminosityBlock->GetEntry(ientry); //read only this branch
	
	// Does the run/LS pass the latest JSON selection?
	if (_json[run][luminosityBlock]==0) {
	  //_badjson.insert(pair<int, int>(run, lbn));
	  ++_nbadevts_json;
	  //return false;
	  continue;
	}
      } // doJSON

      if (debugevent) cout << "Read in entry" << endl << flush;

      // Do this now before trigger
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (debugevent) cout << "Check run, LS in JSON" << endl << flush;
      
      if (mrunls.find(run)==mrunls.end()) ++nrun;
      if (mrunls[run].find(luminosityBlock)==mrunls[run].end()) ++nls;
      ++nevt;
      mrunls[run][luminosityBlock] = 1;

      if (doTrigger) {

	if (debugevent) cout << "Check trigger" << endl << flush;

	bool fired = false;
	for (int i = 0; i != ntrg; ++i) {
	  fired = (fired || (*mtrg[vtrg[i]]));
	  if (*mtrg[vtrg[i]]) htrg->Fill(i);
	}
	if (!fired) {
	  ++_nbadevts_trg;
	  continue;
	}
      }
      ++_ngoodevts;

      if (debugevent) cout << "Redo JEC" << endl << flush;
      
      // Redo JEC right after event cuts but before anything else
      // Do not re-sort
      int njet = nJet;
      for (int i  = 0; i != njet; ++i) {
	double rawJetPt = Jet_pt[i] * (1.0 - Jet_rawFactor[i]);
	double rawJetMass = Jet_mass[i] * (1.0 - Jet_rawFactor[i]);
	jec->setJetPt(rawJetPt);
	jec->setJetEta(Jet_eta[i]);
	if (isRun2) {
	  jec->setJetA(Jet_area[i]);
	  jec->setRho(Rho_fixedGridRhoFastjetAll);
	}
	//double corr = jec->getCorrection();
	vector<float> v = jec->getSubCorrections();
	double corr = v.back();
	double res = (v.size()>1 ? v[v.size()-1]/v[v.size()-2] : 1.);
	Jet_RES[i] = 1./res;
	Jet_deltaJES[i] = (1./corr) / (1.0 - Jet_rawFactor[i]);
	Jet_pt[i] = corr * rawJetPt;
	Jet_mass[i] = corr * rawJetMass;
	Jet_rawFactor[i] = (1.0 - 1.0/corr);
      } // for njet

      // Match leading jets to HLT objects
      for (int i  = 0; i != min(njet,kMaxTrigJet); ++i) {
	Jet_hltPt[i] = 0;
	Jet_hltPtMax[i] = 0;
	Jet_hltPtNear[i] = 0;
	Jet_hltPtClose[i] = 0;
	double dptmin(9999);
	double drmin(999);
	//double drmin2(999);
	for (unsigned int j  = 0; j != nTrigObjJMEAK4; ++j) {
	  //if (TrigObjJMEAK4_id[j]==1) {
	  if (true) {
	    double dphi = DELTAPHI(Jet_phi[i], TrigObjJMEAK4_phi[j]);
	    double deta = fabs(Jet_eta[i] - TrigObjJMEAK4_eta[j]);
	    double dr  = sqrt(dphi*dphi + deta*deta);
	    // Closest in Pt HLT object within jet radius
	    if (dr < 0.4 && fabs(Jet_pt[i] - TrigObjJMEAK4_pt[j]) < dptmin) {
	      Jet_hltPtClose[i] = TrigObjJMEAK4_pt[j];
	      dptmin = fabs(Jet_pt[i] - TrigObjJMEAK4_pt[j]);
	    }
	    // Hardest HLT object within jet radius
	    if (dr < 0.4 && Jet_hltPtMax[i] < TrigObjJMEAK4_pt[j]) {
	      Jet_hltPtMax[i] = TrigObjJMEAK4_pt[j];
	    }
	    // Nearest HLT object within jet radius
	    if (dr < 0.4 && dr < drmin) {
	      Jet_hltPtNear[i] = TrigObjJMEAK4_pt[j];
	      drmin = dr;
	    }
	    // Nearest HLT object within half of jet radius
	    //if (dr < 0.2 && dr < drmin2) {
	    //Jet_hltPt[i] = TrigObjJMEAK4_pt[j];
	    //drmin2 = dr;
	    //}
	  } // TrigObjJMEAK4==Jet
	} // for nTrigObjJMEAK4
	// Require best match within dR<R/2 and also closest in pT at dR<R
	if (drmin<0.2 && Jet_hltPtClose[i]==Jet_hltPtNear[i]) {
	  Jet_hltPt[i] = Jet_hltPtClose[i];
	}
      } // for njet

      if (debugevent) cout << "Sum four-vectors for MHT" << endl << flush;
      
      //int njet3 = 0;
      int njetn = 0;

      if (isRun2)
	p4met.SetPtEtaPhiM(MET_pt,0,MET_phi,0);
      else
	p4met.SetPtEtaPhiM(PuppiMET_pt,0,PuppiMET_phi,0);
      //p4mht.SetPtEtaPhiM(0,0,0,0);
      //p4mht2.SetPtEtaPhiM(0,0,0,0);
      //p4mhtc.SetPtEtaPhiM(0,0,0,0);
      //p4mhtc3.SetPtEtaPhiM(0,0,0,0);
      p4m0.SetPtEtaPhiM(0,0,0,0);
      p4m2.SetPtEtaPhiM(0,0,0,0);
      p4mn.SetPtEtaPhiM(0,0,0,0);
      p4mu.SetPtEtaPhiM(0,0,0,0);
      p4mo.SetPtEtaPhiM(0,0,0,0);

      bool ismultijet = true; // multijet pre-setting
      //Jet_pt[1] < 0.6*Jet_pt[0] && Jet_pt[1]>30. && fabs(Jet_eta[1])<2.5 &&
      //Jet_pt[2] < 0.6*Jet_pt[0] && Jet_pt[2]>30. && fabs(Jet_eta[2])<2.5);
	 //Jet_pt[1] < 0.7*Jet_pt[0] && Jet_pt[1]>30. && fabs(Jet_eta[1])<2.5 &&
	 //Jet_pt[2] < 0.7*Jet_pt[0] && Jet_pt[2]>30. && fabs(Jet_eta[2])<2.5);
      // add dphi veto and dphijet (recoil,lead) later
      // check both Jet[1] and Jet[2] incase big JEC changes. Add |eta|<2.5
	   
      p4lead.SetPtEtaPhiM(0,0,0,0);
      p4recoil.SetPtEtaPhiM(0,0,0,0);
      p4other.SetPtEtaPhiM(0,0,0,0);
      p4m3.SetPtEtaPhiM(0,0,0,0);
      p4mo3.SetPtEtaPhiM(0,0,0,0);
      p4leadRES.SetPtEtaPhiM(0,0,0,0);
      p4recoilRES.SetPtEtaPhiM(0,0,0,0);
      for (int i = 0; i != njet; ++i) {

	p4.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

	for (int j = i+1; j != njet; ++j) {
	  p4s.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	  double dr = p4.DeltaR(p4s);
	  hdri->Fill(dr,w);
	  if (i<2 && j<2) hdri2->Fill(dr,w);
	  if (fabs(p4.Eta())<1.3 || fabs(p4s.Eta())<1.3) {
	    hdrib->Fill(dr,w);
	    if (i<2 && j<2) hdri2b->Fill(dr,w);
	    if (fabs(p4.Eta())<1.3 && fabs(p4s.Eta())<1.3) {
	      hdribb->Fill(dr,w);
	      if (i<2 && j<2) hdri2bb->Fill(dr,w);
	    }
	  }
	} // for j

	hpt->Fill(p4.Pt(),w);
	if (fabs(p4.Eta())<3.0)
	  hpt30->Fill(p4.Pt(),w);
	heta->Fill(p4.Eta(),w);
	h2etaphi->Fill(p4.Eta(),p4.Phi(),w);
	if (p4.Pt()>40.) {
	  heta40->Fill(p4.Eta(),w);
	  h2etaphi40->Fill(p4.Eta(),p4.Phi(),w);
	}
	
	// Jet veto maps
	if (doJetveto) {

	  for (int itrg = 0; itrg != ntrg; ++itrg) {
	    string &trg = vtrg[itrg];
	    if (!(*mtrg[trg])) continue;
	    
	    jetvetoHistos *h = mhjv[trg]; assert(h);
	    double abseta = fabs(p4.Eta());
	    double pt = p4.Pt();
	    double passveto = true; // add jet veto maps
	    if (pt >= h->ptmin && pt < h->ptmax) {
	      h->heta_pretrg->Fill(p4.Eta(), w);
	    }
	    if (pt >= h->ptmin && pt < h->ptmax &&
		abseta >= h->absetamin && abseta < h->absetamax) {
	      h->hpt->Fill(pt, w);
	      h->heta->Fill(p4.Eta(), w);
	      h->h2etaphi->Fill(p4.Eta(), p4.Phi(), w);
	      if (passveto) {
		h->hpt_veto->Fill(pt, w);
		h->heta_veto->Fill(p4.Eta(), w);
	      }

	      if (doPFComposition) {
		h->p2chf->Fill(p4.Eta(), p4.Phi(), Jet_chHEF[i], w);
		h->p2nef->Fill(p4.Eta(), p4.Phi(), Jet_neEmEF[i], w);
		h->p2nhf->Fill(p4.Eta(), p4.Phi(), Jet_neHEF[i], w);
	      }
	    }
	  }
	} // doJetveto

	// Inclusive jets
	if (doIncjet) {
	  for (int itrg = 0; itrg != ntrg; ++itrg) {
	    
	    string &trg = vtrg[itrg];
	    if (!(*mtrg[trg])) continue;
	    
	    incjetHistos *h = mhij[trg];

	    h->hall->Fill(p4.Pt(), w);
	    if (p4.Pt() >= h->ptmin && p4.Pt() < h->ptmax &&
		fabs(p4.Rapidity()) > h->absetamin &&
		fabs(p4.Rapidity()) < h->absetamax)
	      h->hsel->Fill(p4.Pt(), w);
	    if (fabs(p4.Rapidity())<1.3)
	      h->hpt13->Fill(p4.Pt(), w);
	    h->h2pt->Fill(p4.Eta(), p4.Pt(), w);
	    int iy = int(fabs(p4.Rapidity()) / 0.5);
	    if (iy<h->ny) h->vpt[iy]->Fill(p4.Pt(), w);
	  } // for itrg 
        }
		 

	//p4mht -= p4;
	//if (i<2) p4mht2 -= p4;
	//if (fabs(p4.Eta())<3.0 || i<2) {
	//p4mhtc -= p4;
	//if (njet3<3) {
	//  p4mhtc3 -= p4;
	//}
	//++njet3;
	//}

	// All jets as MET proxy
	p4m0 -= p4; // all jets (any pT)

	// L2Res HDM (dijet)
	if (i<2) {  // two leading jets
	  p4m2 -= p4;
	}
	else if (fabs(p4.Eta())<3.0 && njetn<3) { // soft jets
	  p4mn -= p4;
	  p4mo -= p4;
	  ++njetn;
	}
	else { // other "unclustered" jets
	  p4mu -= p4;
	  p4mo -= p4;
	}

	// L3Res HDM (multijet)
	if (i==0 && p4.Pt()>30.) { // leading jet
	  p4lead += p4;
	  p4m3 -= p4;
	  p4leadRES += Jet_RES[i]*p4;
	}
	else if (i>0 && p4.Pt()>15. && fabs(p4.Eta())<2.5 &&
		 DELTAPHI(p4.Phi(),p4lead.Phi())>1.0) {
	  //else if (i>0 && p4.Pt()>30. &&
	  //	 DELTAR(p4.Phi(),p4lead.Phi(),p4.Eta(),p4lead.Eta())>1.0) {
	  // recoil jets
	  p4recoil += p4;
	  p4m3 -= p4;
	  p4recoilRES += Jet_RES[i]*p4;
	}
	else { // all other "unclustered" jets
	  p4other += p4;
	  p4mo3 -= p4;
	}
	// veto nearby jets for multijet topology
	if (i>0 && p4.Pt()>30. && fabs(p4.Eta())<2.5 &&
	    DELTAPHI(p4.Phi(),p4lead.Phi())<=1.0)
	  //if (i>0 && p4.Pt()>30. &&
	  //DELTAR(p4.Phi(),p4lead.Phi(),p4.Eta(),p4lead.Eta())<=1.0)
	  ismultijet = false;//(ismultijet && false);
	
      } // for i in njet
      hnjet->Fill(njet,w);
      //hnjet30->Fill(njet3,w);

      // also check recoil phi for multijet selection
      double ptrecoil = p4recoil.Pt();
      double dphirecoil = DELTAPHI(p4lead.Phi(), p4recoil.Phi());
      //ismultijet = (ismultijet && dphirecoil>2.7);
      //ismultijet = (ismultijet && fabs(dphirecoil-TMath::Pi())<0.3);
      ismultijet =
	(ismultijet && fabs(dphirecoil-TMath::Pi())<0.3 &&
	 Jet_pt[0]>30. && fabs(Jet_eta[0])<1.3 &&
	 Jet_pt[1] < 0.6*ptrecoil && Jet_pt[1]>30. && fabs(Jet_eta[1])<2.5 &&
	 Jet_pt[2] < 0.6*ptrecoil && Jet_pt[2]>30. && fabs(Jet_eta[2])<2.5);

      
      // dijet pre-selection
      if (njet>=2) {

	if (debugevent) cout << "Dijet analysis" << endl << flush;
      
	// both leading jets act as tag and probe in turn
	for (int itag = 0; itag != 2; ++itag) {

	  // tag and probe jet selection
	  int iprobe = (itag == 0 ? 1 : 0);
	  p4t.SetPtEtaPhiM(Jet_pt[itag], Jet_eta[itag], Jet_phi[itag],
			   Jet_mass[itag]);
	  p4p.SetPtEtaPhiM(Jet_pt[iprobe], Jet_eta[iprobe], Jet_phi[iprobe],
			   Jet_mass[iprobe]);

	  // dijet observables
	  double eta = p4p.Eta();
	  double pttag = p4t.Pt();
	  double ptprobe = p4p.Pt();
	  double ptave = 0.5*(pttag+ptprobe);
	  double asymm = (ptprobe - pttag) / ptave;

	  double dphi = DELTAPHI(p4t.Phi(),p4p.Phi());
	  double dr = p4t.DeltaR(p4p);
	  //double mpf = 1 + (p4mht.Vect().Dot(p4b.Vect()))/ptave;
	  //double mpf2 = 1 + (p4mht2.Vect().Dot(p4b.Vect()))/ptave;
	  //double mpfc = 1 + (p4mhtc.Vect().Dot(p4b.Vect()))/ptave;
	  //double mpfc3 = 1 + (p4mhtc3.Vect().Dot(p4b.Vect()))/ptave;

	  // bisector axis
	  p4b.SetPtEtaPhiM(0,0,0,0);
	  p4b += p4t;
	  p4b -= p4p;
	  p4b.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi(),0.);
	  p4b *= 1./p4b.Pt();
	  p4bx.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi()+0.5*TMath::Pi(),0.);

	  double m0b = 1 + (p4m0.Vect().Dot(p4b.Vect()))/ptave;
	  double m2b = 1 + (p4m2.Vect().Dot(p4b.Vect()))/ptave;
	  double mnb = 0 + (p4mn.Vect().Dot(p4b.Vect()))/ptave;
	  double mub = 0 + (p4mu.Vect().Dot(p4b.Vect()))/ptave;
	  double mob = 0 + (p4mo.Vect().Dot(p4b.Vect()))/ptave;

	  double m0bx = 1 + (p4m0.Vect().Dot(p4bx.Vect()))/ptave;
	  double m2bx = 1 + (p4m2.Vect().Dot(p4bx.Vect()))/ptave;

	  // central axis + pttag norm&binning
	  p4c.SetPtEtaPhiM(0,0,0,0);
	  p4c += p4t;
	  p4c.SetPtEtaPhiM(p4c.Pt(),0.,p4c.Phi(),0.);
	  p4c *= 1./p4c.Pt();
	  p4cx.SetPtEtaPhiM(p4c.Pt(),0.,p4c.Phi()+0.5*TMath::Pi(),0.);

	  double m0c = 1 + (p4m0.Vect().Dot(p4c.Vect()))/pttag;
	  double m2c = 1 + (p4m2.Vect().Dot(p4c.Vect()))/pttag;
	  double mnc = 0 + (p4mn.Vect().Dot(p4c.Vect()))/pttag;
	  double muc = 0 + (p4mu.Vect().Dot(p4c.Vect()))/pttag;
	  double moc = 0 + (p4mo.Vect().Dot(p4c.Vect()))/pttag;

	  double m0cx = 1 + (p4m0.Vect().Dot(p4cx.Vect()))/pttag;
	  double m2cx = 1 + (p4m2.Vect().Dot(p4cx.Vect()))/pttag;

	  // forward axis ("backwards" to tag hemisphere) + ptprobe norm&binning
	  p4f.SetPtEtaPhiM(0,0,0,0);
	  p4f -= p4p;
	  p4f.SetPtEtaPhiM(p4f.Pt(),0.,p4f.Phi(),0.);
	  p4f *= 1./p4f.Pt();
	  p4fx.SetPtEtaPhiM(p4f.Pt(),0.,p4f.Phi()+0.5*TMath::Pi(),0.);

	  double m0f = 1 + (p4m0.Vect().Dot(p4f.Vect()))/ptprobe;
	  double m2f = 1 + (p4m2.Vect().Dot(p4f.Vect()))/ptprobe;
	  double mnf = 0 + (p4mn.Vect().Dot(p4f.Vect()))/ptprobe;
	  double muf = 0 + (p4mu.Vect().Dot(p4f.Vect()))/ptprobe;
	  double mof = 0 + (p4mo.Vect().Dot(p4f.Vect()))/ptprobe;

	  double m0fx = 1 + (p4m0.Vect().Dot(p4fx.Vect()))/ptprobe;
	  double m2fx = 1 + (p4m2.Vect().Dot(p4fx.Vect()))/ptprobe;

	  // Dijet mass
	  p4dj = p4t; p4dj += p4p;
	  double mjj = p4dj.M();
	  double deta = fabs(p4p.Eta()-p4t.Eta());
	  
	  dijetHistos *h(0);
	  bool isdijet = fabs(p4t.Eta())<1.3 && dphi>2.7 && fabs(asymm)<maxa;
	  if (isdijet) {

	    for (int itrg = 0; itrg != ntrg; ++itrg) {
	      
	      if (debugevent) cout << "Check trigger #"<<itrg<<" for dijet "
				   << endl << flush;
	      
	      string &trg = vtrg[itrg];
	      if (!(*mtrg[trg])) continue;

	      if (doJetveto) {
		jetvetoHistos *h = mhjv[trg];
		if (ptave >= h->ptmin && ptave < h->ptmax &&
		    fabs(eta) >= h->absetamin && fabs(eta) < h->absetamax) {
		  h->p2asymm->Fill(eta, p4p.Phi(), asymm, w);

		  if (doPFComposition) {
		    h->p2chftp->Fill(eta, p4p.Phi(), Jet_chHEF[iprobe], w);
		    h->p2neftp->Fill(eta, p4p.Phi(), Jet_neEmEF[iprobe], w);
		    h->p2nhftp->Fill(eta, p4p.Phi(), Jet_neHEF[iprobe], w);
		  }
		}
	      } // doJetveto

	      if (doDijet) {
		dijetHistos *h = mhdj[trg];
		h->h2all->Fill(eta, ptave, w);
		if (ptave >= h->ptmin && ptave < h->ptmax &&
		    fabs(eta) >= h->absetamin && fabs(eta) < h->absetamax) {
		  h->h2sel->Fill(eta, ptave, w);
		}
		//if (fabs(eta) >= h->absetamin && fabs(eta) < h->absetamax) {
		{
		  double jes2 = Jet_RES[iprobe] / Jet_RES[itag];
		  h->p2jes2->Fill(eta, ptave, jes2, w);

		  h->p2m0->Fill(eta, ptave, m0b, w);
		  h->p2m0x->Fill(eta, ptave, m0bx, w);
		  h->p2m2->Fill(eta, ptave, m2b, w);
		  h->p2m2x->Fill(eta, ptave, m2bx, w);

		  h->p2m0ab->Fill(eta, ptave, m0b, w);
		  h->p2m2ab->Fill(eta, ptave, m2b, w);
		  h->p2mnab->Fill(eta, ptave, mnb, w);
		  h->p2muab->Fill(eta, ptave, mub, w);
		}
		//if (pttag >= h->ptmin && pttag < h->ptmax &&
		//if (fabs(eta) >= h->absetamin && fabs(eta) < h->absetamax) {
		{
		  h->p2m0tc->Fill(eta, pttag, m0c, w);
		  h->p2m2tc->Fill(eta, pttag, m2c, w);
		  h->p2mntc->Fill(eta, pttag, mnc, w);
		  h->p2mutc->Fill(eta, pttag, muc, w);
		}
		//if (ptprobe >= h->ptmin && ptprobe < h->ptmax &&
		//if (fabs(eta) >= h->absetamin && fabs(eta) < h->absetamax) {
		{
		  h->p2m0pf->Fill(eta, ptprobe, m0f, w);
		  h->p2m2pf->Fill(eta, ptprobe, m2f, w);
		  h->p2mnpf->Fill(eta, ptprobe, mnf, w);
		  h->p2mupf->Fill(eta, ptprobe, muf, w);
		}
	      } // doDijet
	    
	      if (doDijetOrig) {
	      // pTave binning
	      int iptave = hptbins->FindBin(ptave)-1;
	      if (iptave>=0 && iptave<npt) {
		dijetHistosOrig *h = mhdjo[trg][iptave];

		h->hpta->Fill(ptave, w);
		if (fabs(eta)>2.8) h->hptaf->Fill(ptave, w);

		h->hmjj2->Fill(mjj);
		if (deta<1.3) h->hmjj213->Fill(mjj);
		
		h->hna->Fill(eta, w);
		
		h->pm0ab->Fill(eta, m0b, w);
		h->pm2ab->Fill(eta, m2b, w);
		h->pmnab->Fill(eta, mnb, w);
		h->pmuab->Fill(eta, mub, w);
		h->pmoab->Fill(eta, mob, w);
		
		h->pm0ac->Fill(eta, m0c, w);
		h->pm2ac->Fill(eta, m2c, w);
		h->pmnac->Fill(eta, mnc, w);
		h->pmuac->Fill(eta, muc, w);
		h->pmoac->Fill(eta, moc, w);
		
		h->pm0af->Fill(eta, m0f, w);
		h->pm2af->Fill(eta, m2f, w);
		h->pmnaf->Fill(eta, mnf, w);
		h->pmuaf->Fill(eta, muf, w);
		h->pmoaf->Fill(eta, mof, w);
		
		h->h2m0ab->Fill(eta, m0b, w);
		h->h2m0abx->Fill(eta, m0bx, w);
		h->h2m2ab->Fill(eta, m2b, w);
		h->h2m2abx->Fill(eta, m2bx, w);
		
		h->p2etaphia->Fill(eta, p4p.Phi(), asymm, w);

		h->pjes->Fill(eta, (1.0-Jet_rawFactor[iprobe]), w);
		h->pres->Fill(eta, Jet_RES[iprobe], w);
		h->pdjes->Fill(eta, Jet_deltaJES[iprobe], w);
		if (fabs(eta)<1.3) {
		  h->pjes13->Fill(ptave, (1.0-Jet_rawFactor[iprobe]), w);
		  h->pres13->Fill(ptave, Jet_RES[iprobe], w);
		  h->pdjes13->Fill(ptave, Jet_deltaJES[iprobe], w);
		}
		  
		if (doPFComposition) {
		  h->pchf->Fill(eta, Jet_chHEF[iprobe], w);
		  h->pnhf->Fill(eta, Jet_neHEF[iprobe], w);
		  h->pnef->Fill(eta, Jet_neEmEF[iprobe], w);
		  h->pcef->Fill(eta, Jet_chEmEF[iprobe], w);
		  h->pmuf->Fill(eta, Jet_muEF[iprobe], w);
		  h->phef->Fill(eta, Jet_hfEmEF[iprobe], w);
		  h->phhf->Fill(eta, Jet_hfHEF[iprobe], w);

		  if (fabs(eta)<1.3) {
		    h->pchf13->Fill(ptave, Jet_chHEF[iprobe], w);
		    h->pnhf13->Fill(ptave, Jet_neHEF[iprobe], w);
		    h->pnef13->Fill(ptave, Jet_neEmEF[iprobe], w);
		    h->pcef13->Fill(ptave, Jet_chEmEF[iprobe], w);
		    h->pmuf13->Fill(ptave, Jet_muEF[iprobe], w);
		  }
		} // doPFcomposition
	      } // pTave binning

	      // pTtag binning
	      int ipttag = hptbins->FindBin(pttag)-1;
	      if (ipttag>=0 && ipttag<npt) {
		dijetHistosOrig *h = mhdjo[trg][ipttag];

		h->hptt->Fill(pttag, w);
		if (fabs(eta)>2.8) h->hpttf->Fill(pttag, w);

		hmjj2->Fill(mjj);
		if (deta<1.3) hmjj213->Fill(mjj);
		
		h->hnt->Fill(eta, w);

		h->pm0tb->Fill(eta, m0b, w);
		h->pm2tb->Fill(eta, m2b, w);
		h->pmntb->Fill(eta, mnb, w);
		h->pmutb->Fill(eta, mub, w);
		h->pmotb->Fill(eta, mob, w);
		
		h->pm0tc->Fill(eta, m0c, w);
		h->pm2tc->Fill(eta, m2c, w);
		h->pmntc->Fill(eta, mnc, w);
		h->pmutc->Fill(eta, muc, w);
		h->pmotc->Fill(eta, moc, w);
		
		h->h2m0tb->Fill(eta, m0b, w);
		h->h2m0tbx->Fill(eta, m0bx, w);
		h->h2m2tb->Fill(eta, m2b, w);
		h->h2m2tbx->Fill(eta, m2bx, w);
		
		h->h2m0tc->Fill(eta, m0c, w);
		h->h2m0tcx->Fill(eta, m0cx, w);
		h->h2m2tc->Fill(eta, m2c, w);
		h->h2m2tcx->Fill(eta, m2cx, w);
		
		h->p2etaphit->Fill(eta, p4p.Phi(), p4p.Pt()/pttag, w);

		// Trigger matching studies
		if (doTriggerMatch && Jet_hltPtMax[itag] > h->trgpt) {
		  phashlt->Fill(eta, Jet_hltPtMax[iprobe]>0 ? 1 : 0, w);
		  h->phashlt->Fill(eta, Jet_hltPtMax[iprobe]>0 ? 1 : 0, w);
		  pishlt->Fill(eta, Jet_hltPt[iprobe]>0 ? 1 : 0, w);
		  h->pishlt->Fill(eta, Jet_hltPt[iprobe]>0 ? 1 : 0, w);
		  if (Jet_hltPt[iprobe]>0) {
		    h->pismax->Fill(eta, Jet_hltPt[iprobe]==Jet_hltPtMax[iprobe] ? 1 : 0,w);
		    h->pjeshlt->Fill(eta, Jet_hltPt[iprobe]/Jet_pt[iprobe], w);
		    h->pjeshltmax->Fill(eta, Jet_hltPtMax[iprobe]/Jet_pt[iprobe], w);
		    h->h2jeshlt->Fill(eta,Jet_hltPt[iprobe]/Jet_pt[iprobe], w);
		    h->ppthlt->Fill(eta, Jet_hltPt[iprobe], w);
		    h->pptoff->Fill(eta, Jet_pt[iprobe], w);
		    h->ppttag->Fill(eta, Jet_pt[itag], w);

		    pismax->Fill(eta, Jet_hltPt[iprobe]==Jet_hltPtMax[iprobe] ? 1 : 0,w);
		    pjeshlt->Fill(eta, Jet_hltPt[iprobe]/Jet_pt[iprobe], w);
		    pjeshltmax->Fill(eta,Jet_hltPtMax[iprobe]/Jet_pt[iprobe],w);
		    h2jeshlt->Fill(eta,Jet_hltPt[iprobe]/Jet_pt[iprobe], w);
		    ppthlt->Fill(eta, Jet_hltPt[iprobe], w);
		    pptoff->Fill(eta, Jet_pt[iprobe], w);
		    ppttag->Fill(eta, Jet_pt[itag], w);
		  } // good HLT match
		  if (Jet_hltPtNear[iprobe]>0) {
		    // Check if nearest is also closest in pT
		    h->pisclose->Fill(eta, Jet_hltPtNear[iprobe]==Jet_hltPtClose[iprobe] ? 1 : 0,w);
		    pisclose->Fill(eta, Jet_hltPtNear[iprobe]==Jet_hltPtClose[iprobe] ? 1 : 0,w);
		    h->pjeshltnear->Fill(eta, Jet_hltPtNear[iprobe]/Jet_pt[iprobe], w);
		    pjeshltnear->Fill(eta,Jet_hltPtNear[iprobe]/Jet_pt[iprobe], w);
		  } // nearest HLT match
		  if (Jet_hltPtClose[iprobe]>0) {
		    // Check if closest in pT is also nearest
		    h->pisnear->Fill(eta, Jet_hltPtClose[iprobe]==Jet_hltPtNear[iprobe] ? 1 : 0,w);
		    pisnear->Fill(eta, Jet_hltPtClose[iprobe]==Jet_hltPtNear[iprobe] ? 1 : 0,w);
		    h->pjeshltclose->Fill(eta, Jet_hltPtClose[iprobe]/Jet_pt[iprobe], w);
		    pjeshltclose->Fill(eta, Jet_hltPtClose[iprobe]/Jet_pt[iprobe], w);
		  } // closest HLT match
		  
		} // doTriggerMatch && tag has HLT match above threshold
	      } // pTtag binning

	      int iptprobe = hptbins->FindBin(ptprobe)-1;
	      if (iptprobe>=0 && iptprobe<npt) {
		dijetHistosOrig *h = mhdjo[trg][iptprobe];

		h->hptp->Fill(ptprobe, w);
		if (fabs(eta)>2.8) h->hptpf->Fill(ptprobe, w);

		hmjj2->Fill(mjj);
		if (deta<1.3) hmjj213->Fill(mjj);
		
		h->hnp->Fill(eta, w);

		h->pm0pb->Fill(eta, m0b, w);
		h->pm2pb->Fill(eta, m2b, w);
		h->pmnpb->Fill(eta, mnb, w);
		h->pmupb->Fill(eta, mub, w);
		h->pmopb->Fill(eta, mob, w);
		
		h->pm0pf->Fill(eta, m0f, w);
		h->pm2pf->Fill(eta, m2f, w);
		h->pmnpf->Fill(eta, mnf, w);
		h->pmupf->Fill(eta, muf, w);
		h->pmopf->Fill(eta, mof, w);
		
		h->h2m0pb->Fill(eta, m0b, w);
		h->h2m0pbx->Fill(eta, m0bx, w);
		h->h2m2pb->Fill(eta, m2b, w);
		h->h2m2pbx->Fill(eta, m2bx, w);
		
		h->h2m0pf->Fill(eta, m0f, w);
		h->h2m0pfx->Fill(eta, m0fx, w);
		h->h2m2pf->Fill(eta, m2f, w);
		h->h2m2pfx->Fill(eta, m2fx, w);
		
		h->p2etaphip->Fill(eta, p4p.Phi(), p4t.Pt()/ptprobe, w);
	      } // pTprobe binning
	      } // doDijetOrig
	      
	    } // for itrg
	  } // dijet tag-and-probe selection
	  
	  // dijet without deltaphi cut
	  if (fabs(p4t.Eta()<1.3) && fabs(asymm)<maxa) {

	    if (ptave>=40) {
	      hdphi->Fill(dphi, w);
	      h2dphi->Fill(p4p.Eta(),dphi, w);
	      hdr->Fill(dr);
	      if (fabs(p4p.Eta())<1.3) {
		hdphib->Fill(dphi, w);
		hdrb->Fill(dr, w);
	      }
	    }
	  }

	  // dijet without asymmetry cut
	  if (fabs(p4t.Eta())<1.3 && dphi>2.7) {

	    hpta->Fill(ptave, w);
	    hptt->Fill(p4t.Pt(), w);
	    hptp->Fill(p4p.Pt(), w);
	    
	    if (ptave>=40) {
	      hav->Fill(asymm, w);
	      if (fabs(p4p.Eta())<1.3) {
		havb->Fill(asymm, w);
	      }

	      //hmv->Fill(mpf, w);
	      //hmva->Fill(mpf, w);
	      //hmvc->Fill(mpfc, w);
	      //hmvc3->Fill(mpfc3, w);
	      //hmv2->Fill(mpf2, w);
	    }
	    if (p4t.Pt()>=40) {
	      hat->Fill(p4p.Pt() / p4t.Pt(), w);
	      if (fabs(p4p.Eta())<1.3) {
		hatb->Fill(p4p.Pt() / p4t.Pt(), w);
	      }
	      //hmt->Fill(mpf, w);
	    }
	    if (p4p.Pt()>=40) {
	      hap->Fill(p4t.Pt() / p4p.Pt(), w);
	      if (fabs(p4p.Eta())<1.3) {
		hapb->Fill(p4t.Pt() / p4p.Pt(), w);
	      }
	      //hmp->Fill(mpf, w);
	    }

	  } // dijet without asymmetry cut
	} // for itag
      } // dijet

      // Multijet selection
      if (ismultijet && doMultijet) {


	if (debugevent) cout << "Analyze multijet" << endl << flush;
	  
	// pTave binning
	double ptlead = p4lead.Pt();
	double ptrecoil = p4recoil.Pt();
	double ptave = 0.5*(ptlead+ptrecoil);
	
	// Bisector axis
	p4b.SetPtEtaPhiM(0,0,0,0);
	p4b -= p4lead;
	p4b += p4recoil;
	p4b.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi(),0.);
	p4b *= 1./p4b.Pt();
	//p4bx.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi()+0.5*TMath::Pi(),0.);

	// Projection to transverse plane (is this necessary?)
	p4m0.SetPtEtaPhiM(p4m0.Pt(),0.,p4m0.Phi(),0.);
	p4m3.SetPtEtaPhiM(p4m3.Pt(),0.,p4m3.Phi(),0.);
	p4mo3.SetPtEtaPhiM(p4mo3.Pt(),0.,p4mo3.Phi(),0.);
	
	double m0b = 1 + (p4m0.Vect().Dot(p4b.Vect()))/ptave;
	double m3b = 1 + (p4m3.Vect().Dot(p4b.Vect()))/ptave;
	double mob = 0 + (p4mo3.Vect().Dot(p4b.Vect()))/ptave;

	p4l.SetPtEtaPhiM(0,0,0,0);
	p4l -= p4lead;
	p4l.SetPtEtaPhiM(p4l.Pt(),0.,p4l.Phi(),0.);
	p4l *= 1./p4l.Pt();

	double m0l = 1 + (p4m0.Vect().Dot(p4l.Vect()))/ptlead;
	double m3l = 1 + (p4m3.Vect().Dot(p4l.Vect()))/ptlead;
	double mol = 0 + (p4mo3.Vect().Dot(p4l.Vect()))/ptlead;

	p4r.SetPtEtaPhiM(0,0,0,0);
	p4r += p4recoil;
	p4r.SetPtEtaPhiM(p4r.Pt(),0.,p4r.Phi(),0.);
	p4r *= 1./p4r.Pt();

	double m0r = 1 + (p4m0.Vect().Dot(p4r.Vect()))/ptrecoil;
	double m3r = 1 + (p4m3.Vect().Dot(p4r.Vect()))/ptrecoil;
	double mor = 0 + (p4mo3.Vect().Dot(p4r.Vect()))/ptrecoil;

	for (int itrg = 0; itrg != ntrg; ++itrg) {
	  
	  string &trg = vtrg[itrg];
	  if (!(*mtrg[trg])) continue;

	  multijetHistos *h = mhmj[trg];
	  
	  h->hna->Fill(ptave, w);
	  h->hnl->Fill(ptlead, w);
	  h->hnr->Fill(ptrecoil, w);

	  //h->pres->Fill(ptave, Jet_RES[0], w);
	  double res = (p4leadRES.Pt()/p4recoilRES.Pt()) /
	    (p4lead.Pt()/p4recoil.Pt());
	  h->presa->Fill(ptave, res, w);
	  h->presl->Fill(ptlead, res, w);
	  h->presr->Fill(ptrecoil, res, w);

	  h->pm0a->Fill(ptave, m0b, w);
	  h->pm2a->Fill(ptave, m3b, w);
	  h->pmna->Fill(ptave, mob, w);
	  h->pmua->Fill(ptave, 0., w);
	  h->pmoa->Fill(ptave, mob+0., w);

	  h->pm0l->Fill(ptlead, m0l, w);
	  h->pm2l->Fill(ptlead, m3l, w);
	  h->pmnl->Fill(ptlead, mol, w);
	  h->pmul->Fill(ptlead, 0., w);
	  h->pmol->Fill(ptlead, mol+0., w);

	  h->pm0r->Fill(ptrecoil, m0r, w);
	  h->pm2r->Fill(ptrecoil, m3r, w);
	  h->pmnr->Fill(ptrecoil, mor, w);
	  h->pmur->Fill(ptrecoil, 0., w);
	  h->pmor->Fill(ptrecoil, mor+0., w);

	  h->h2m0a->Fill(ptave, m0b, w);
	  h->h2m2a->Fill(ptave, m3b, w);
	  if (ptave>1.25*h->trgpt)
	    h->hcosdphi->Fill(cos(DELTAPHI(p4lead.Phi(),p4recoil.Phi())), w);
	} // for itrg
      } // ismultijet
      
      //h2mhtvsmet->Fill(p4met.Pt(), p4mht.Pt(), w);
      h2mhtvsmet->Fill(p4met.Pt(), p4m0.Pt(), w);
   } // for jentry
   cout << endl << flush;

   cout << "Finished looping over " << nevt << " of which " << _ngoodevts
	<< " passed trigger. Start writing file." << endl << flush;
   fout->Write();
   fout->Close();
   cout << "File written and closed." << endl << flush;

   if (doJSON)
     cout << Form("Found %d bad events according to new JSON (events cut)",_nbadevts_json) << endl;
   if (doTrigger) {
     cout << Form("Found %d bad events according to trigger bits (events cut)",_nbadevts_trg) << endl;
     cout << Form("Found %d bad events not in fwd trigger phase space (events cut)",_nbadevts_fwdtrg) << endl;
   }

  cout << "Processed " << nrun << " runs, "
       << nls << " luminosity blocks and " << nevt << " events" << endl;
  cout << "Saving these to file rootfiles/jmenano.json for brilcalc" << endl;

   ofstream fjson("rootfiles/jmenano.json");
   fjson << "{" << endl;
   for (map<int, map<int,int> >::iterator it = mrunls.begin();
	it != mrunls.end(); ++it) {

     int run = it->first;
     int ls_prev(0), ls(0), ok(0);
     bool firstblock(true), newblock(true);

     if (it!=mrunls.begin()) fjson << "," << endl;
     fjson << "  \"" << run << "\": [";
     for (map<int,int>::iterator jt = it->second.begin();
	  jt != it->second.end(); ++jt) {

       ok = jt->second;
       if (ok) {
	 ls_prev = ls;
	 ls = jt->first;
       }
       else
	 continue;

       newblock = (firstblock || (ls - ls_prev > 1));
       if (newblock && ok) {
	 if (firstblock) {
	   fjson << "[" << ls << ", ";
	   firstblock = false;
	 }
	 else {
	   fjson << ls_prev << "], [" << ls << ", ";
	 }
       }
     } // for jt
     fjson << ls << "]]";
   } // for it
   fjson << "}" << endl;
   
   cout << Form("Analyzed %d events",_ngoodevts) << endl;
   cout << "Saving these to " << fout->GetName() << " for drawJMENANO.C" << endl;

   //h2mhtvsmet->Draw("COLZ");
}


bool DijetHistosFill::LoadJSON()
{
  // Get the JSON files from here:
  // - /eos/user/c/cmsdqm/www/CAF/certification/Collisions22/
  
  // Golden 1.44/fb
  // string json = "rootfiles/Cert_Collisions2022_355100_356615_Golden.json";
  // Golden JSON, 4.86/fb
  //string json = "rootfiles/Cert_Collisions2022_355100_357550_Golden..json";
  // Golden JSON, 7.67/fb
  //string json = "rootfiles/Cert_Collisions2022_355100_357900_Golden.json";
  // Golden JSON BCDEF, 9.71/fb
  //string json = "rootfiles/Cert_Collisions2022_355100_359812_Golden.json";
  // Golden JSON BCDEF, 14.6/fb
    string json = "rootfiles/Cert_Collisions2022_355100_360491_Golden.json";
// Golden JSON RunB, 0.0846/fb
//string json = "rootfiles/Cert_Collisions2022_eraB_355100_355769_Golden.json";
// Golden JSON RunC, 4.84/fb
//string json = "rootfiles/Cert_Collisions2022_eraC_355862_357482_Golden.json"
// Golden JSON RunD, 2.74/fb
//string json = "rootfiles/Cert_Collisions2022_eraD_357538_357900_Golden.json";
    if (isRun2==1 || isRun2==2)
      json="rootfiles/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt";
    if (isRun2==3)
      json="rootfiles/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt";
    if (isRun2==4)
    json="rootfiles/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt";
  
  cout << "Processing LoadJSON() with " + json + " ..." << flush;
  ifstream file(json, ios::in);
  if (!file.is_open()) return false;
  char c;
  string s, s2;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  if (c!='{') return false;
  while (file >> s and sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (debug) cout << Form("\"%d\": ",rn);
    
    while (file.get(c) and c==' ') {};
    if (debug) { cout << Form("%c",c) << flush; assert(c=='['); }
    ++nrun;

    bool endrun = false;
    while (!endrun and file >> s >> s2 and sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) {
      if (debug) cout << Form("[%d,%d,%s]",ls1,ls2,s1);

      for (int ls = ls1; ls != ls2+1; ++ls) {
        _json[rn][ls] = 1;
        ++nls;
      }

      s2 = s1;
      endrun = (s2=="]," || s2=="]}");
      if (debug and !endrun and s2!=",") { cout << string("s1: ")+s2 << endl << flush; assert(s2==","); }
    } // while ls
    if (debug) cout << endl;

    if (s2=="]}") continue;
    else if (debug and s2!="],") cout << string("s2: ")+s2 << endl << flush;
    assert(s2=="],");
  } // while run
  //if (s2!="]}") { PrintInfo(string("s3: ")+s2,true); return false; }
  if (s2!="]}") { cout <<  string("s3: ")+s2 << endl; return false; }

  cout << string("Called LoadJSON() with ") + json + ":" << endl;
  cout << Form("Loaded %d good runs and %d good lumi sections",nrun,nls) << endl;
  return true;
} // LoadJSON
