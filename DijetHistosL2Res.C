// Purpose: calculate JEC L2Res using MPF (update to HMD) method
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "../jecsys2020/tdrstyle_mod15.C"

TH1D *getHDM(TProfile2D* p2, TProfile2D *p22, TProfile2D *p2n, TProfile2D *p2u,
	     double eta1, double eta2);

void DijetHistosL2Ress(string file, string dir, string bin="");

// Process several directories in a uniform way
void DijetHistosL2Res() {

  // Data before final L2Res
  /*
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2016APV_v26c.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2016GH_v26c.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2017_v26.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2018_v26c.root","Dijet2");
  */
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_Run2_v26c.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_Run2_v26c.root","Dijet2","pf");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_Run2_v26c.root","Dijet2","tc");
  //DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2017_v26.root","Dijet2");

  // MC before JER SF
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v26.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v26.root","Dijet2","pf");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v26.root","Dijet2","tc");
  //DijetHistosL2Ress("rootfiles/jmenano_mc_cmb_UL2017MG_v26.root","Dijet2");

  // MC after JER SF
  /*
  DijetHistosL2Ress("rootfiles/jmenano_mc_cmb_UL2016APVMG_v27.root","Dijet2");
  DijetHistosL2Ress("rootfiles/jmenano_mc_cmb_UL2016MG_v27.root","Dijet2");
  DijetHistosL2Ress("rootfiles/jmenano_mc_cmb_UL2017MG_v27.root","Dijet2");
  DijetHistosL2Ress("rootfiles/jmenano_mc_cmb_UL2018MG_v27.root","Dijet2");
  */
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v27.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v27.root","Dijet2","pf");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v27.root","Dijet2","tc");
  //DijetHistosL2Ress("rootfiles/jmenano_mc_cmb_UL2017MG_v27.root","Dijet2");

}

// Update cmb.root to add L2Res results
void DijetHistosL2Ress(string file, string dir, string bin) {

  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile(file.c_str(),"UPDATE");
  assert(f && !f->IsZombie());

  // Determine if MC
  TString s(file.c_str());
  bool isMC = (s.Contains("_mc"));
  if (isMC) cout << "DijetHistosL2Ress(\"" << file << "\") (MC)" << endl;
  else      cout << "DijetHistosL2Ress(\"" << file << "\") (Data)" << endl;
  
  f->cd(dir.c_str());
  TDirectory *d = gDirectory;

  curdir->cd();
  
  const char *cb = bin.c_str();
  TProfile2D *p2m0 = (TProfile2D*)d->Get(Form("p2m0%s",cb)); assert(p2m0);
  TProfile2D *p2m2 = (TProfile2D*)d->Get(Form("p2m2%s",cb)); assert(p2m2);
  TProfile2D *p2mn = (TProfile2D*)d->Get(Form("p2mn%s",cb)); assert(p2mn);
  TProfile2D *p2mu = (TProfile2D*)d->Get(Form("p2mu%s",cb)); assert(p2mu);
  
  // Create 2D histograms for storing HDM JEC
  TH2D *h2jec = p2m0->ProjectionXY(Form("h2jec%s",cb)); h2jec->Reset();

  for (int i = 1; i != h2jec->GetNbinsX()+1; ++i) {

    // Calculate HDM (ratio of tag and probe)
    double eta = h2jec->GetXaxis()->GetBinCenter(i);
    TH1D *h1jec(0);
    h1jec = getHDM(p2m0,p2m2,p2mn,p2mu,eta,eta);

    for (int j = 1; j != h2jec->GetNbinsY()+1; ++j) {

      double jec = h1jec->GetBinContent(j);
      double err = h1jec->GetBinError(j);
      h2jec->SetBinContent(i, j, jec);
      h2jec->SetBinError(i, j, err);
    } // for j
    
    delete h1jec;
  } // for i

  d->cd();

  h2jec->Write(Form("h2jec%s",cb),TObject::kOverwrite);
  
  f->Write(nullptr,TObject::kOverwrite);
  curdir->cd();  

  delete h2jec;

} // DijetHistosL2Ress

// Calculate HDM JEC for one given eta slice at a time
TH1D *getHDM(TProfile2D* p2, TProfile2D *p22, TProfile2D *p2n, TProfile2D *p2u,
	     double eta1, double eta2) {
  
  string s = Form("_%s",p2->GetName());
  const char *c = s.c_str();

  int i1 = p2->GetXaxis()->FindBin(eta1);
  int i2 = p2->GetXaxis()->FindBin(eta2);

  TProfile *p1 = p2->ProfileY(Form("p1%s",c),i1,i2,"");
  TProfile *p12 = p22->ProfileY(Form("p12%s",c),i1,i2,"");
  TProfile *p1n = p2n->ProfileY(Form("p1n%s",c),i1,i2,"");
  TProfile *p1u = p2u->ProfileY(Form("p1u%s",c),i1,i2,"");

  // Calculate HDM JEC
  TH1D *h1jec = p1->ProjectionX(Form("h1jec_%s",c));
  for (int i = 1; i != h1jec->GetNbinsX()+1; ++i) {

    double mpf = p1->GetBinContent(i);
    double err = p1->GetBinError(i);

    // Add proper hdm calculation here, now just test with MPF
    double hdm = mpf;
    double jec = hdm;
    
    h1jec->SetBinContent(i, jec);
    h1jec->SetBinError(i, err);
  }

  //cout << h1 << ", " << (*h1) << endl;

  delete p1;
  delete p12;
  delete p1n;
  delete p1u;
  
  return h1jec;
} // getHDM
