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

double DELTAPHI(double a, double b) {
  double phi1 = max(a,b);
  double phi2 = min(a,b);
  double d = phi1-phi2;
  if (d>TMath::Pi()) d -= TMath::TwoPi();
  return fabs(d);
}

class basicHistos {
public:
  string trg;
  int trgpt;
  double ptmin;
  double ptmax;

  TH1D *hpta;
  TH1D *hptt;
  TH1D *hptp;
  TH1D *hptaf;
  TH1D *hpttf;
  TH1D *hptpf;

  TProfile *pm0ab;
  TProfile *pm2ab;
  TProfile *pmnab;
  TProfile *pmuab;
  TProfile *pmoab;

  TProfile *pm0tb;
  TProfile *pm2tb;
  TProfile *pmntb;
  TProfile *pmutb;
  TProfile *pmotb;

  TProfile *pm0pb;
  TProfile *pm2pb;
  TProfile *pmnpb;
  TProfile *pmupb;
  TProfile *pmopb;
 
  TProfile *pm0ac;
  TProfile *pm2ac;
  TProfile *pmnac;
  TProfile *pmuac;
  TProfile *pmoac;

  TProfile *pm0tc;
  TProfile *pm2tc;
  TProfile *pmntc;
  TProfile *pmutc;
  TProfile *pmotc;

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

  TProfile2D *p2etaphia, *p2etaphit, *p2etaphip;

  //TH1D *hmjji, *hmjji13, *hmjj1, *hmjj113;
  TH1D *hmjj2, *hmjj213;
  //TH1D *hmjj, *hmjj113;

  TProfile *pjes, *pdjes, *pjes13, *pdjes13;
  TProfile *pchf, *pnhf, *pnef, *pcef, *pmuf, *phef, *phhf;
  TProfile *pchf13, *pnhf13, *pnef13, *pcef13, *pmuf13;

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

   if (isMC) fChain->SetBranchStatus("genWeight",1); // v2

   fChain->SetBranchStatus("run",1);
   fChain->SetBranchStatus("luminosityBlock",1);
   fChain->SetBranchStatus("event",1);

   //fChain->SetBranchStatus("HLT_DiPFJetAve40",1);
   //fChain->SetBranchStatus("HLT_PFJet40",1); // no events
   //fChain->SetBranchStatus("HLT_AK8PFJet40",1);
   //fChain->SetBranchStatus("HLT_PFJetFwd40",1);
   //fChain->SetBranchStatus("HLT_AK8PFJetFwd40",1);

   // Listing of available triggers
   vector<string> vtrg;
   vtrg.push_back("HLT_DiPFJetAve40");
   vtrg.push_back("HLT_DiPFJetAve60");
   vtrg.push_back("HLT_DiPFJetAve80");
   vtrg.push_back("HLT_DiPFJetAve140");
   vtrg.push_back("HLT_DiPFJetAve200");
   vtrg.push_back("HLT_DiPFJetAve260");
   vtrg.push_back("HLT_DiPFJetAve320");
   vtrg.push_back("HLT_DiPFJetAve400");
   vtrg.push_back("HLT_DiPFJetAve500");

   vtrg.push_back("HLT_PFJet40");
   vtrg.push_back("HLT_PFJet60");
   vtrg.push_back("HLT_PFJet80");
   //vtrg.push_back("HLT_PFJet110");
   vtrg.push_back("HLT_PFJet140");
   vtrg.push_back("HLT_PFJet200");
   vtrg.push_back("HLT_PFJet260");
   vtrg.push_back("HLT_PFJet320");
   vtrg.push_back("HLT_PFJet450");
   vtrg.push_back("HLT_PFJet500");
   vtrg.push_back("HLT_PFJet550");

   vtrg.push_back("HLT_DiPFJetAve60_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve80_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve100_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve160_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve220_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve300_HFJEC");

   vtrg.push_back("HLT_PFJetFwd15");
   vtrg.push_back("HLT_PFJetFwd25");
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
   int ntrg = vtrg.size();

   for (int i = 0; i != ntrg; ++i) {
     fChain->SetBranchStatus(vtrg[i].c_str(),1);
   }

   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_mass",1);

   fChain->SetBranchStatus("Jet_rawFactor",1);

   fChain->SetBranchStatus("Jet_chHEF",1);  // h+
   fChain->SetBranchStatus("Jet_neHEF",1);  // h0
   fChain->SetBranchStatus("Jet_neEmEF",1); // gamma
   fChain->SetBranchStatus("Jet_chEmEF",1); // e
   fChain->SetBranchStatus("Jet_muEF",1);   // mu
   fChain->SetBranchStatus("Jet_hfEmEF",1); // HFe
   fChain->SetBranchStatus("Jet_hfHEF",1);  // HFh

   fChain->SetBranchStatus("PuppiMET_pt",1);
   fChain->SetBranchStatus("PuppiMET_phi",1);

   fChain->SetBranchStatus("Flag_METFilters",1);

   // Trigger studies
   fChain->SetBranchStatus("nTrigObj",1);
   fChain->SetBranchStatus("TrigObj_pt",1);
   fChain->SetBranchStatus("TrigObj_eta",1);
   fChain->SetBranchStatus("TrigObj_phi",1);
   fChain->SetBranchStatus("TrigObj_id",1); // Jet==1, FatJet==6?
   // https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_8/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L136-L180
   
   // Select appropriate L1RC for type-I MET L1L2L3-RC calculation
   //FactorizedJetCorrector *jecl1rc(0);
   //string& ds = dataset;
   //if (ds=="2016B") jecl1rc = getFJC("Summer19UL16APV_RunBCD_V7_DATA_L1RC");
   // Redo JEC
   FactorizedJetCorrector *jec(0);
   jec = getFJC("","Winter22Run3_V1_MC_L2Relative","","");

   
   TLorentzVector p4met, p4dj;
   TLorentzVector p4, p4s, p4mht, p4mht2, p4mhtc, p4mhtc3, p4t, p4p;
   TLorentzVector p4b, p4bx, p4c, p4cx, p4f, p4fx;
   TLorentzVector p4m0, p4m2, p4mn, p4mu, p4mo;
   TFile *fout = new TFile(Form("rootfiles/jmenano_%s_out.root",
				isMC ? "mc" : "data"), "RECREATE");

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
   TH1D *hnjet30 = new TH1D("hnjet30","hnjet30",500,0,500);
   //TH1D *hndup = new TH1D("hndup","hndup",100,0,100);

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
   
   /*
   TH2D *h2etaphidr0 = new TH2D("h2etaphidr0","#DeltaR=0;#eta;#phi",nx,vx,
				72,-TMath::Pi(),TMath::Pi());
   TH2D *h2ijdr0 = new TH2D("hijdr0","#DeltaR=0;i;j",50,0,50,50,0,50);
   TH1D *havdr0 = new TH1D("havddr0","Asymmetry #DeltaR=0",400,-2,2);
   TH1D *hpt1dr0 = new TH1D("hpt1dr0","Pt1 #DeltaR=0",500,0,500);
   TH1D *hpt2dr0 = new TH1D("hpt2dr0","Pt2 #DeltaR=0",500,0,500);
   */

   TH1D *hav = new TH1D("hav","Asymmetry",400,-2,2);
   TH1D *hat = new TH1D("hat","Pt,probe / Pt,tag",400,0,4);
   TH1D *hap = new TH1D("hap","Pt,tag / Pt,probe",400,0,4);

   TH1D *havb = new TH1D("havb","Asymmetry Barrel",400,-2,2);
   TH1D *hatb = new TH1D("hatb","Pt,probe / Pt,tag Barrel",400,0,4);
   TH1D *hapb = new TH1D("hapb","Pt,tag / Pt,probe Barrel",400,0,4);
   
   TH1D *hmv = new TH1D("hmv","MPF PtAve",800,-2,6);
   TH1D *hmt = new TH1D("hmt","MPF PtTag",800,-2,6);
   TH1D *hmp = new TH1D("hmp","MPF PtProbe",800,-2,6);

   TH1D *hmva = new TH1D("hmva","MPFA PtAve",800,-2,6);
   TH1D *hmvc = new TH1D("hmvc","MPFC PtAve",800,-2,6);
   TH1D *hmvc3 = new TH1D("hmvc3","MPFC3 PtAve",800,-2,6);
   TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",800,-2,6);
   //TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",400,-2,2);

   // PF composition plots
   // Copy L2Res histograms for multiple pT bins
   const double vpt[] = {40, 50, 60, 70, 85, 100, 125, 155, 180, 210, 250, 300,
			 350, 400, 500, 600, 800, 1000, 1200, 1500,
			 1800, 2100, 2400, 2700, 3000};
   const int npt = sizeof(vpt)/sizeof(vpt[0])-1;
   TH1D *hptbins = new TH1D("hptbins",";p_{T} (GeV);N_{events}",npt,vpt);

   TProfile *pjes = new TProfile("pjes","JES",nx,vx);
   TProfile *pdjes = new TProfile("pdjes","#DELTAJES",nx,vx);
   TProfile *pchf = new TProfile("pchf","CHF",nx,vx);
   TProfile *pnhf = new TProfile("pnhf","NHF",nx,vx);
   TProfile *pnef = new TProfile("pnef","NEF",nx,vx);
   TProfile *pcef = new TProfile("pcef","CEF",nx,vx);
   TProfile *pmuf = new TProfile("pmuf","MUF",nx,vx);
   TProfile *phef = new TProfile("phfhf","HEF",nx,vx);
   TProfile *phhf = new TProfile("phfef","HHF",nx,vx);

   TProfile *pjes13 = new TProfile("pjes13","JES",npt,vpt);
   TProfile *pdjes13 = new TProfile("pdjes13","#DELTAJES",npt,vpt);
   TProfile *pchf13 = new TProfile("pchf13","CHF",npt,vpt);
   TProfile *pnhf13 = new TProfile("pnhf13","NHF",npt,vpt);
   TProfile *pnef13 = new TProfile("pnef13","NEF",npt,vpt);
   TProfile *pcef13 = new TProfile("pcef13","CEF",npt,vpt);
   TProfile *pmuf13 = new TProfile("pmuf13","MUF",npt,vpt);
   
   // Extra fun
   // Run vs BX => need BX
   //TH1D *hmjji = new TH1D("hmjji","Dijet mass, inclusive",1000,0,1000);
   //TH1D *hmjji13 = new TH1D("hmjji13","Dijet mass, |#DeltaEta|<1.3, incl.",
   //			    1000,0,1000);
   //TH1D *hmjj1 = new TH1D("hmjj1","Dijet mass, 1-lead",1000,0,1000);
   //TH1D *hmjj113 = new TH1D("hmjj113","Dijet mass, |#DeltaEta|<1.3, 1-lead",
   //			    1000,0,1000);
   TH1D *hmjj2 = new TH1D("hmjj2","Dijet mass, 2-lead",1000,0,1000);
   TH1D *hmjj213 = new TH1D("hmjj213","Dijet mass, |#DeltaEta|<1.3, 2-lead",
   			    1000,0,1000);
   //TH1D *hmjj = new TH1D("hmjj","Dijet mass, tag-and-probe",1000,0,1000);
   //TH1D *hmjj13 = new TH1D("hmjj13","Dijet mass, |#DeltaEta|<1.3, TnP",
   //	                    1000,0,1000);


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

   //vector<basicHistos*> vh(npt);
   map<string, vector<basicHistos*> > mh;
   const bool doPtBins = true;
   if (doPtBins) {

     for (int itrg = 0; itrg != ntrg; ++itrg) {

       //vector<basicHistos*> vh(npt);
       fout->mkdir(vtrg[itrg].c_str());
       fout->cd(vtrg[itrg].c_str());
       TDirectory *dout = gDirectory;
       vector<basicHistos*> &vh = mh[vtrg[itrg]];
       vh.resize(npt);

       // Figure out trigger pT threshold from the name
       int trgpt(0), nfound(0);
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_PFJet%d",&trgpt);
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_PFJetFwd%d",&trgpt);
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_DiPFJetAve%d",&trgpt);
       if (nfound!=1)
	 nfound = sscanf(vtrg[itrg].c_str(),"HLT_DiPFJetAve%d_HFJEC",&trgpt);
       assert(nfound==1);
       assert(trgpt!=0);
       
       for (int ipt = 0; ipt != npt; ++ipt) {
       
	 dout->mkdir(Form("Pt_%d_%d",int(vpt[ipt]),int(vpt[ipt+1])));
	 dout->cd(Form("Pt_%d_%d",int(vpt[ipt]),int(vpt[ipt+1])));
	 basicHistos *h = new basicHistos();
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
	 h->pdjes = (TProfile*)pdjes->Clone();
	 h->pchf = (TProfile*)pchf->Clone();
	 h->pnhf = (TProfile*)pnhf->Clone();
	 h->pnef = (TProfile*)pnef->Clone();
	 h->pcef = (TProfile*)pcef->Clone();
	 h->pmuf = (TProfile*)pmuf->Clone();
	 h->phef = (TProfile*)phef->Clone();
	 h->phhf = (TProfile*)phhf->Clone();

	 h->pjes13 = (TProfile*)pjes13->Clone();
	 h->pdjes13 = (TProfile*)pdjes13->Clone();
	 h->pchf13 = (TProfile*)pchf13->Clone();
	 h->pnhf13 = (TProfile*)pnhf13->Clone();
	 h->pnef13 = (TProfile*)pnef13->Clone();
	 h->pcef13 = (TProfile*)pcef13->Clone();
	 h->pmuf13 = (TProfile*)pmuf13->Clone();

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
       } // for ipt
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

      // Do this now before trigger
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (mrunls.find(run)==mrunls.end()) ++nrun;
      if (mrunls[run].find(luminosityBlock)==mrunls[run].end()) ++nls;
      ++nevt;
      mrunls[run][luminosityBlock] = 1;

      if (doTrigger) {
	//b_HLT_DiPFJetAve40->GetEntry(ientry);
	//b_HLT_PFJet40->GetEntry(ientry); // no events?
	//b_HLT_AK8PFJet40->GetEntry(ientry);
	//b_HLT_PFJetFwd40->GetEntry(ientry);
	//b_HLT_AK8PFJetFwd40->GetEntry(ientry);
	//if (HLT_DiPFJetAve40==false && HLT_PFJet40==false &&
	//  HLT_AK8PFJet40==false &&
	//  HLT_PFJetFwd40==false && HLT_AK8PFJetFwd40==false) {
	//++_nbadevts_trg;
	//continue;
	//}
	//if (HLT_DiPFJetAve40==false HLT_PFJetFwd40==false) {
	//if (HLT_DiPFJetAve40==false) {
	//if (HLT_PFJet40==false) { // no events?
	bool fired = false;
	for (int i = 0; i != ntrg; ++i) {
	  fired = (fired || (*mtrg[vtrg[i]]));
	}
	//if (HLT_PFJetFwd40==false) {
	if (!fired) {
	  ++_nbadevts_trg;
	  continue;
	}
      }
      ++_ngoodevts;

      //nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      // Extra event cuts for trigger eta
      //if (doTrigger) {
	//if (!(HLT_DiPFJetAve40==true || HLT_PFJet40==true ||
	//    HLT_AK8PFJet40==true)) {
	//if (HLT_DiPFJetAve40==false) {
	//if (HLT_PFJet40==false) {
	//if (HLT_PFJetFwd40==true) {
      //if (!(nJet>1 && (fabs(Jet_eta[0])>2.853 || fabs(Jet_eta[1])>2.853))) {
      //    ++_nbadevts_fwdtrg;
      //    continue;
      //  }
      //} // fwd-only trigger
      //}

      // Redo JEC right after event cuts but before anything else
      // Do not re-sort
      int njet = nJet;
      for (int i  = 0; i != njet; ++i) {
	jec->setJetPt(Jet_pt[i] * (1.0 - Jet_rawFactor[i])); 
	jec->setJetEta(Jet_eta[i]);
	double corr = jec->getCorrection();
	Jet_deltaJES[i] = (1./corr) / (1.0 - Jet_rawFactor[i]);
	Jet_pt[i] = corr * Jet_pt[i] * (1.0 - Jet_rawFactor[i]); 
	Jet_mass[i] = corr * Jet_mass[i] * (1.0 - Jet_rawFactor[i]);
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
	for (unsigned int j  = 0; j != nTrigObj; ++j) {
	  if (TrigObj_id[j]==1) {
	    double dphi = DELTAPHI(Jet_phi[i], TrigObj_phi[j]);
	    double deta = fabs(Jet_eta[i] - TrigObj_eta[j]);
	    double dr  = sqrt(dphi*dphi + deta*deta);
	    // Closest in Pt HLT object within jet radius
	    if (dr < 0.4 && fabs(Jet_pt[i] - TrigObj_pt[j]) < dptmin) {
	      Jet_hltPtClose[i] = TrigObj_pt[j];
	      dptmin = fabs(Jet_pt[i] - TrigObj_pt[j]);
	    }
	    // Hardest HLT object within jet radius
	    if (dr < 0.4 && Jet_hltPtMax[i] < TrigObj_pt[j]) {
	      Jet_hltPtMax[i] = TrigObj_pt[j];
	    }
	    // Nearest HLT object within jet radius
	    if (dr < 0.4 && dr < drmin) {
	      Jet_hltPtNear[i] = TrigObj_pt[j];
	      drmin = dr;
	    }
	    // Nearest HLT object within half of jet radius
	    //if (dr < 0.2 && dr < drmin2) {
	    //Jet_hltPt[i] = TrigObj_pt[j];
	    //drmin2 = dr;
	    //}
	  } // TrigObj==Jet
	} // for nTrigObj
	// Require best match within dR<R/2 and also closest in pT at dR<R
	if (drmin<0.2 && Jet_hltPtClose[i]==Jet_hltPtNear[i]) {
	  Jet_hltPt[i] = Jet_hltPtClose[i];
	}
      } // for njet

      
      int njet3 = 0;
      int njetn = 0;

      p4met.SetPtEtaPhiM(PuppiMET_pt,0,PuppiMET_phi,0);
      p4mht.SetPtEtaPhiM(0,0,0,0);
      p4mht2.SetPtEtaPhiM(0,0,0,0);
      p4mhtc.SetPtEtaPhiM(0,0,0,0);
      p4mhtc3.SetPtEtaPhiM(0,0,0,0);
      p4m0.SetPtEtaPhiM(0,0,0,0);
      p4m2.SetPtEtaPhiM(0,0,0,0);
      p4mn.SetPtEtaPhiM(0,0,0,0);
      p4mu.SetPtEtaPhiM(0,0,0,0);
      p4mo.SetPtEtaPhiM(0,0,0,0);
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
	p4mht -= p4;
	if (i<2) p4mht2 -= p4;
	if (fabs(p4.Eta())<3.0 || i<2) {
	  p4mhtc -= p4;
	  if (njet3<3) {
	    p4mhtc3 -= p4;
	  }
	  ++njet3;
	}
	
	// L2Res HDM
	p4m0 -= p4;
	if (i<2) { // leading jet
	  p4m2 -= p4;
	}
	else if (fabs(p4.Eta())<3.0 && njetn<3) { // soft jets
	  p4mn -= p4;
	  p4mo -= p4;
	  ++njetn;
	}
	else { // other "unclustered"
	  p4mu -= p4;
	  p4mo -= p4;
	}
      } // for i in njet
      hnjet->Fill(njet,w);
      hnjet30->Fill(njet3,w);

      // dijet selection
      if (njet>=2) {

	// both leading jets act as tag and probe in turn
	for (int itag = 0; itag != 2; ++itag) {

	  int iprobe = (itag == 0 ? 1 : 0);
	  p4t.SetPtEtaPhiM(Jet_pt[itag], Jet_eta[itag], Jet_phi[itag],
			   Jet_mass[itag]);
	  p4p.SetPtEtaPhiM(Jet_pt[iprobe], Jet_eta[iprobe], Jet_phi[iprobe],
			   Jet_mass[iprobe]);
	  // bisector axis
	  p4b.SetPtEtaPhiM(0,0,0,0);
	  p4b += p4t;
	  p4b -= p4p;
	  p4b.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi(),0.);
	  p4b *= 1./p4b.Pt();
	  p4bx.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi()+0.5*TMath::Pi(),0.);

	  double eta = p4p.Eta();
	  double pttag = p4t.Pt();
	  double ptprobe = p4p.Pt();
	  double ptave = 0.5*(pttag+ptprobe);
	  double asymm = (ptprobe - pttag) / ptave;
	  // Maximum asymmetry of 2/3 corresponds to x2 ratio of tag and probe
	  // Permit ~0.7 extra scaling to allow for HF L3Res
	  const double maxa = 10;//0.963;//2./3;
	  //double deltaphi = fabs(p4t.DeltaPhi(p4p));
	  double deltaphi = DELTAPHI(p4t.Phi(),p4p.Phi());
	  double dr = p4t.DeltaR(p4p);
	  double mpf = 1 + (p4mht.Vect().Dot(p4b.Vect()))/ptave;
	  double mpf2 = 1 + (p4mht2.Vect().Dot(p4b.Vect()))/ptave;
	  double mpfc = 1 + (p4mhtc.Vect().Dot(p4b.Vect()))/ptave;
	  double mpfc3 = 1 + (p4mhtc3.Vect().Dot(p4b.Vect()))/ptave;

	  // bisector axis
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
	  //hmjj2->Fill(mjj);
	  //if (deta<1.3) hmjj213->Fill(mjj);
	  
	  basicHistos *h(0);
	  if (fabs(p4t.Eta())<1.3 && deltaphi>2.7 && fabs(asymm)<maxa) {

	    for (int itrg = 0; itrg != ntrg; ++itrg) {
	      
	      string &trg = vtrg[itrg];
	      if (!(*mtrg[trg])) continue;

	      //double eta = p4p.Eta();
	      //if (ptave >= 40 && ptave <50) {
	      int iptave = hptbins->FindBin(ptave)-1;
	      if (iptave>=0 && iptave<npt) {
		//h = vh[iptave];
		h = mh[trg][iptave];

		h->hpta->Fill(ptave, w);
		if (fabs(eta)>2.8) h->hptaf->Fill(ptave, w);

		h->hmjj2->Fill(mjj);
		if (deta<1.3) h->hmjj213->Fill(mjj);
		
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
		h->pdjes->Fill(eta, Jet_deltaJES[iprobe], w);
		h->pchf->Fill(eta, Jet_chHEF[iprobe], w);
		h->pnhf->Fill(eta, Jet_neHEF[iprobe], w);
		h->pnef->Fill(eta, Jet_neEmEF[iprobe], w);
		h->pcef->Fill(eta, Jet_chEmEF[iprobe], w);
		h->pmuf->Fill(eta, Jet_muEF[iprobe], w);
		h->phef->Fill(eta, Jet_hfEmEF[iprobe], w);
		h->phhf->Fill(eta, Jet_hfHEF[iprobe], w);

		if (fabs(eta)<1.3) {
		  h->pjes13->Fill(ptave, (1.0-Jet_rawFactor[iprobe]), w);
		  h->pdjes13->Fill(eta, Jet_deltaJES[iprobe], w);
		  h->pchf13->Fill(ptave, Jet_chHEF[iprobe], w);
		  h->pnhf13->Fill(ptave, Jet_neHEF[iprobe], w);
		  h->pnef13->Fill(ptave, Jet_neEmEF[iprobe], w);
		  h->pcef13->Fill(ptave, Jet_chEmEF[iprobe], w);
		  h->pmuf13->Fill(ptave, Jet_muEF[iprobe], w);
		}
	      }
	      //double pttag = p4t.Pt();
	      //if (pttag >= 40 && pttag <50) {
	      int ipttag = hptbins->FindBin(pttag)-1;
	      if (ipttag>=0 && ipttag<npt) {
		//h = vh[ipttag];
		h = mh[trg][ipttag];

		h->hptt->Fill(pttag, w);
		if (fabs(eta)>2.8) h->hpttf->Fill(pttag, w);

		hmjj2->Fill(mjj);
		if (deta<1.3) hmjj213->Fill(mjj);
		
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
		if (Jet_hltPtMax[itag] > h->trgpt) {
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
		  
		} // tag has HLT match above threshold		
	      }
	      //double ptprobe = p4p.Pt();
	      //if (ptprobe >= 40 && ptprobe <50) {
	      int iptprobe = hptbins->FindBin(ptprobe)-1;
	      if (iptprobe>=0 && iptprobe<npt) {
		//h = vh[iptprobe];
		h = mh[trg][iptprobe];

		h->hptp->Fill(ptprobe, w);
		if (fabs(eta)>2.8) h->hptpf->Fill(ptprobe, w);

		hmjj2->Fill(mjj);
		if (deta<1.3) hmjj213->Fill(mjj);
		
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
	      }
	    } // for itrg
	  } // tag-and-probe selection

	  if (fabs(p4t.Eta()<1.3) && fabs(asymm)<maxa) {

	    if (ptave>=40) {
	      hdphi->Fill(deltaphi, w);
	      h2dphi->Fill(p4p.Eta(),deltaphi, w);
	      hdr->Fill(dr);
	      if (fabs(p4p.Eta())<1.3) {
		hdphib->Fill(deltaphi, w);
		hdrb->Fill(dr, w);
	      }
	    }
	  }

	  if (fabs(p4t.Eta())<1.3 && deltaphi>2.7) {

	    hpta->Fill(ptave, w);
	    hptt->Fill(p4t.Pt(), w);
	    hptp->Fill(p4p.Pt(), w);
	    
	    if (ptave>=40) {
	      hav->Fill(asymm, w);
	      if (fabs(p4p.Eta())<1.3) {
		havb->Fill(asymm, w);
	      }

	      /*
	      cout << "itag="<<itag<<" hav="<<asymm<<" hmv2="<<mpf2
		   << " ptt="<<p4t.Pt()<<" ptp="<<p4p.Pt()<<endl
		   << " phit="<<p4t.Phi()<<" phip="<<p4p.Phi()
		   << " ptave="<<ptave<<endl
		   << " phib="<<p4b.Phi()
		   << " dphi="<<deltaphi 
		   << " magb="<<p4b.Mag()<<endl
		   << " ptmht2="<<p4mht2.Pt()
		   << " phimht2="<<p4mht2.Phi()
		   << endl;
	      */
	      hmv->Fill(mpf, w);
	      hmva->Fill(mpf, w);
	      hmvc->Fill(mpfc, w);
	      hmvc3->Fill(mpfc3, w);
	      hmv2->Fill(mpf2, w);
	      //cout << p4mht.Pt() << " / " << ptave << " vs a=" << mpf << endl;
	      //cout << p4mhtc.Pt() << " / " << ptave << " vs c=" << mpfc << endl;
	      //cout << p4mhtc3.Pt() << " / "<< ptave << " vs c3="<< mpfc3<< endl;
	      //if (fabs(p4p.Eta())<1.3)
	      //hmvb->Fill(mpf);
	    }
	    if (p4t.Pt()>=40) {
	      hat->Fill(p4p.Pt() / p4t.Pt(), w);
	      if (fabs(p4p.Eta())<1.3) {
		hatb->Fill(p4p.Pt() / p4t.Pt(), w);
	      }
	      hmt->Fill(mpf, w);
	    }
	    if (p4p.Pt()>=40) {
	      hap->Fill(p4t.Pt() / p4p.Pt(), w);
	      if (fabs(p4p.Eta())<1.3) {
		hapb->Fill(p4t.Pt() / p4p.Pt(), w);
	      }
	      hmp->Fill(mpf, w);
	    }

	  } // tag and deltaphi
	} // for itag
      } // dijet

      h2mhtvsmet->Fill(p4met.Pt(), p4mht.Pt(), w);
   } // for jentry
   cout << endl << flush;

   fout->Write();
   fout->Close();

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
  // Golden 1.44/fb
  // string json = "rootfiles/Cert_Collisions2022_355100_356615_Golden.json";
  // Golden JSON, 4.86/fb
  //string json = "rootfiles/Cert_Collisions2022_355100_357550_Golden..json";
  // Golden JSON, 7.67/fb
  string json = "rootfiles/Cert_Collisions2022_355100_357900_Golden.json";
// Golden JSON RunB, 0.0846/fb
//string json = "rootfiles/Cert_Collisions2022_eraB_355100_355769_Golden.json";
// Golden JSON RunC, 4.84/fb
//string json = "rootfiles/Cert_Collisions2022_eraC_355862_357482_Golden.json"
// Golden JSON RunD, 2.74/fb
//string json = "rootfiles/Cert_Collisions2022_eraD_357538_357900_Golden.json";
  if (isRun2)
    json="rootfiles/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt";
  
  bool debug = true;
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
