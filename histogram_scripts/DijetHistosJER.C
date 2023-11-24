// Purpose: calculate JER SF using MPFX method
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "../tdrstyle_mod22.C"

string version = "v35";

TH1D *getJER(TProfile2D* p2, TProfile2D *p2x,
	     double eta1, double eta2, TH1D **h1 = 0, TH1D **h1x = 0);
TH1D *getJERZ(TH2D *h2, TH2D *h2x, const char *c = "",
	      TH1D **h1 = 0, TH1D **h1x = 0);

void DijetHistosJERs(string file, string dir);
void drawDijetHistosJER(string sd="",string sm="",string era="");
void drawDijetHistosJERtest();

bool scaleJER = false;

// Process several directories in a uniform way
void DijetHistosJER(string rootdir="../rootfiles", string hadddir="../haddfiles") {


  // Run3 files (v29->v30)
  scaleJER = false;
  DijetHistosJERs(rootdir + "/jmenano_data_cmb_2022CD_JME_"+version+".root","Dijet2");
  DijetHistosJERs(rootdir + "/jmenano_data_cmb_2022E_JME_"+version+".root","Dijet2");
  DijetHistosJERs(rootdir + "/jmenano_data_cmb_2022FG_JME_"+version+".root","Dijet2");
  DijetHistosJERs(rootdir + "/jmenano_data_cmb_2023BCv123_JME_"+version+".root","Dijet2");
  DijetHistosJERs(rootdir + "/jmenano_data_cmb_2023Cv4_JME_v"+version+"root","Dijet2");
  DijetHistosJERs(rootdir + "/jmenano_data_cmb_2023D_JME_"+version+".root","Dijet2");
  //
  scaleJER = true; // something weird with MC RMS
  DijetHistosJERs(rootdir + "/jmenano_mc_cmb_Summer22MG_"+version+".root","Dijet2");
  //DijetHistosJERs("rootfiles/jmenano_mc_cmb_Summer22EEMG_v30.root","Dijet2");
  scaleJER = false;
  
  string mc22 = rootdir + "/jmenano_mc_cmb_Summer22MG_"+version+".root";
  drawDijetHistosJER(rootdir + "/jmenano_data_cmb_2022CD_JME_"+version+".root",mc22,
		     "2022CD_vs_Summer22_"+version+"");
  drawDijetHistosJER(rootdir + "/jmenano_data_cmb_2022E_JME_"+version+".root",mc22,
		     "2022E_vs_Summer22_"+version+"");
  drawDijetHistosJER(rootdir + "/jmenano_data_cmb_2022FG_JME_"+version+".root",mc22,
		     "2022FG_vs_Summer22_"+version+"");
  drawDijetHistosJER(rootdir + "/jmenano_data_cmb_2023BCv123_JME_"+version+".root",mc22,
		     "2023BCv123_vs_Summer22_"+version+"");
  drawDijetHistosJER(rootdir + "/jmenano_data_cmb_2023Cv4_JME_"+version+".root",mc22,
		     "2023Cv4_vs_Summer22_"+version+"");
  drawDijetHistosJER(rootdir + "/jmenano_data_cmb_2023D_JME_"+version+".root",mc22,
		     "2023D_vs_Summer22_"+version+"");


  //Regarding the pt dependent JER SFs. Here is the pre approval of the DP Note, where we kept a slide for the SFs https://indico.cern.ch/event/1302471/#8-dp-notes-for-run3
  // MC pThat / qScale (LHE_HT)
  // qScale is in JMENANO, XPOG decides what stays in NANO?
  // => find two files where pthatmax is not is, and ones where it is
  
  /*
  DijetHistosJERs("rootfiles/jmenano_data_cmb_v21ul16.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_data_cmb_v21ul16.root","Dijet/JER");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","Dijet/JER");
  */
  /*
  DijetHistosJERs("rootfiles/jmenano_data_cmb_v22ul16.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v22ul16mg.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","Dijet2");
  */
  /*
  DijetHistosJERs("rootfiles/jmenano_data_cmb_v22ul16.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v23ul16mg.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v23ul16flat.root","Dijet2");
  */
  //drawDijetHistosJER();
  //drawDijetHistosJERtest();

  // Before MC JER SF
  /*
  DijetHistosJERs("haddfiles/jmenano_data_cmb_UL2016APV_v26c.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2016APVMG_v26.root","Dijet2");
  DijetHistosJERs("haddfiles/jmenano_data_cmb_UL2016GH_v26c.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2016MG_v26.root","Dijet2");
  DijetHistosJERs("haddfiles/jmenano_data_cmb_UL2017_v26.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2017MG_v26.root","Dijet2");
  DijetHistosJERs("haddfiles/jmenano_data_cmb_UL2018_v26c.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2018MG_v26.root","Dijet2");
  DijetHistosJERs("haddfiles/jmenano_data_cmb_Run2_v26c.root","Dijet2");
  DijetHistosJERs("haddfiles/jmenano_mc_cmb_Run2_v26.root","Dijet2");

  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2016APV_v26c.root",
  		     "rootfiles/jmenano_mc_cmb_UL2016APVMG_v26.root",
  		     "UL2016APV_ZB_v26c");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2016GH_v26c.root",
  		     "rootfiles/jmenano_mc_cmb_UL2016MG_v26.root",
  		     "UL2016GH_ZB_v26c");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2017_v26.root",
  		     "rootfiles/jmenano_mc_cmb_UL2017MG_v26.root",
  		     "UL2017_ZB_v26c");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2018_v26c.root",
  		     "rootfiles/jmenano_mc_cmb_UL2018MG_v26.root",
  		     "UL2018_ZB_v26c");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_Run2_v26c.root",
  		     "haddfiles/jmenano_mc_cmb_Run2_v26.root",
  		     "Run2_ZB_v26c");
  */

  /*
  // After MC JER SF
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2016APVMG_v27.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2016MG_v27.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2017MG_v27.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2018MG_v27.root","Dijet2");
  DijetHistosJERs("haddfiles/jmenano_mc_cmb_Run2_v27.root","Dijet2");

  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2016APV_v26c.root",
  		     "rootfiles/jmenano_mc_cmb_UL2016APVMG_v27.root",
  		     "UL2016APV_ZB_v27");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2016GH_v26c.root",
  		     "rootfiles/jmenano_mc_cmb_UL2016MG_v27.root",
  		     "UL2016GH_ZB_v27");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2017_v26.root",
  		     "rootfiles/jmenano_mc_cmb_UL2017MG_v27.root",
  		     "UL2017_ZB_v27");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_UL2018_v26c.root",
  		     "rootfiles/jmenano_mc_cmb_UL2018MG_v27.root",
  		     "UL2018_ZB_v27");
  drawDijetHistosJER("haddfiles/jmenano_data_cmb_Run2_v26c.root",
  		     "haddfiles/jmenano_mc_cmb_Run2_v27.root",
  		     "Run2_ZB_v27");
  */

  // New Run3 files from Iita and Mikael
  // Iita_20230814/*_v1.root -> Iita_20230824_jetveto/*_JME_v1.root
  // 2022
  /*
  scaleJER = false;
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022C_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022D_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022E_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022F_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022G_JME_v1.root","Dijet2");
  // 2023
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023B_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023Cv123_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023BCv123_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023Cv4_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023D_JME_v1.root","Dijet2");
  // Main combos (after checking stability for L2Res)
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022CD_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022FG_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022EFG_JME_v1.root","Dijet2");
  //DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2022BCv123_JME_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023Cv4D_JME_v1.root","Dijet2");
  */
  
  // Plot results vs UL2018 for now
  //DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2018MG_v27.root","Dijet2"); // with JER SF
  //DijetHistosJERs("rootfiles/jmenano_mc_cmb_UL2018MG_v26.root","Dijet2"); // no JER SF
  /*
  scaleJER = true;
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_mc_cmb_Summer22_v1.root","Dijet2");
  DijetHistosJERs("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_mc_cmb_Summer22EE_v1.root","Dijet2");
  scaleJER = false;
  */
  
  // mc_cmb_UL2018MG_v26 -> mc_cmb_Summer22_v1
  /*
  const char *mc22 = "../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_mc_cmb_Summer22_v1.root";
  //const char *mc22 = "rootfiles/jmenano_mc_cmb_UL2018MG_v26.root";
  //const char *mc22ee = "../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_mc_cmb_Summer22EE_v1.root";
  const char *mc22ee = "rootfiles/jmenano_mc_cmb_UL2018MG_v26.root";
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022C_JME_v1.root",mc22,"2022C_v1_vs_Summer22_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022D_JME_v1.root",mc22,"2022D_v1_vs_Summer22_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022E_JME_v1.root",mc22ee,"2022E_v1_vs_Summer22EE_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022F_JME_v1.root",mc22ee,"2022F_v1_vs_Summer22EE_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022G_JME_v1.root",mc22ee,"2022G_v1_vs_Summer22EE_v1");
  //
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023B_JME_v1.root",mc22,"2023B_v1_vs_Summer22_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023Cv123_JME_v1.root",mc22,"2023Cv123_v1_vs_Summer22_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023BCv123_JME_v1.root",mc22,"2023BCv123_v1_vs_Summer22_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023Cv4_JME_v1.root",mc22,"2023Cv4_v1_vs_Summer22_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023D_JME_v1.root",mc22,"2023D_v1_vs_Summer22_v1");
  //
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022CD_JME_v1.root",mc22ee,"2022CD_v1_vs_Summer22EE_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022FG_JME_v1.root",mc22ee,"2022FG_v1_vs_Summer22EE_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/jmenano_data_cmb_2022EFG_JME_v1.root",mc22ee,"2022EFG_v1_vs_Summer22EE_v1");
  drawDijetHistosJER("../jecsys3/rootfiles/Iita_20230824_jetveto/nano_data_cmb_2023Cv4D_JME_v1.root",mc22,"2023Cv4D_v1_vs_Summer22_v1");
  */
}

// Update cmb.root to add JER results
void DijetHistosJERs(string file, string dir) {

  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile(file.c_str(),"UPDATE");
  assert(f && !f->IsZombie());

  // Determine if MC
  TString s(file.c_str());
  bool isMC = (s.Contains("_mc"));
  if (isMC) cout << "DijetHistosJERs(\"" << file << "\") (MC)" << endl;
  else      cout << "DijetHistosJERs(\"" << file << "\") (Data)" << endl;
  
  f->cd(dir.c_str());
  TDirectory *d = gDirectory;

  curdir->cd();
  
  TProfile2D *p2m0 = (TProfile2D*)d->Get("p2m0"); assert(p2m0);
  TProfile2D *p2m0x = (TProfile2D*)d->Get("p2m0x"); assert(p2m0x);
  
  // Calculate MPFX JER (without PU+UE components yet)
  TH1D *h1jer13(0), *h1m0s13(0), *h1m0xs13(0);
  h1jer13 = getJER(p2m0,p2m0x,-1.3,1.3,&h1m0s13,&h1m0xs13);
  h1jer13->SetName("h1jer13");
  h1m0s13->SetName("h1m0s13");
  h1m0xs13->SetName("h1m0xs13");

  // Create 2D histograms for storing JER2, JER, RMS(MPF) and RMS(MPFX)
  TH2D *h2jer2 = p2m0->ProjectionXY("h2jer2"); h2jer2->Reset();
  TH2D *h2jer = (TH2D*)h2jer2->Clone("h2jer");
  TH2D *h2m0s = (TH2D*)h2jer2->Clone("h2m0s");
  TH2D *h2m0xs = (TH2D*)h2jer2->Clone("h2m0xs");

  // Alternative longer calculation starting from MPF2X
  TProfile2D *p2m2 = (TProfile2D*)d->Get("p2m2"); assert(p2m2);
  TProfile2D *p2m2x = (TProfile2D*)d->Get("p2m2x"); assert(p2m2x);
  TProfile2D *p2mo = (TProfile2D*)d->Get("p2mnu"); assert(p2mo);
  TProfile2D *p2mox = (TProfile2D*)d->Get("p2mnux"); assert(p2mox);
  //
  TH1D *h1bal13(0), *h1m2s13(0), *h1m2xs13(0);
  h1bal13 = getJER(p2m2,p2m2x,-1.3,1.3,&h1m2s13,&h1m2xs13);
  h1bal13->SetName("h1bal13");
  h1m2s13->SetName("h1m2s13");
  h1m2xs13->SetName("h1m2xs13");
  TH2D *h2bal2 = p2m2->ProjectionXY("h2bal2"); h2bal2->Reset();
  TH2D *h2bal = (TH2D*)h2bal2->Clone("h2bal");
  TH2D *h2m2s = (TH2D*)h2bal2->Clone("h2m2s");
  TH2D *h2m2xs = (TH2D*)h2bal2->Clone("h2m2xs");
  //
  TH1D *h1fsr13(0), *h1mos13(0), *h1moxs13(0);
  h1fsr13 = getJER(p2mo,p2mox,-1.3,1.3,&h1mos13,&h1moxs13);
  h1fsr13->SetName("h1fsr13");
  h1mos13->SetName("h1mos13");
  h1moxs13->SetName("h1moxs13");
  TH2D *h2fsr2 = p2mo->ProjectionXY("h2fsr2"); h2fsr2->Reset();
  TH2D *h2fsr = (TH2D*)h2fsr2->Clone("h2fsr");
  TH2D *h2mos = (TH2D*)h2fsr2->Clone("h2mos");
  TH2D *h2moxs = (TH2D*)h2fsr2->Clone("h2moxs");
  //
  TH1D *h1jrb13 = (TH1D*)h1bal13->Clone("h1jrb13"); h1jrb13->Reset();
  TH2D *h2jrb2 = p2mo->ProjectionXY("h2jrb2"); h2jrb2->Reset();
  TH2D *h2jrb = (TH2D*)h2jrb2->Clone("h2jrb");
  
  for (int i = 1; i != h1jrb13->GetNbinsX()+1; ++i) {
    double bal13 = h1bal13->GetBinContent(i);
    double bal13e = h1bal13->GetBinError(i);
    double fsr13 = h1fsr13->GetBinContent(i);
    double fsr13e = h1fsr13->GetBinError(i);
    double jrb13 = sqrt(max(bal13*bal13e - fsr13*fsr13,0.));
    double jrb13e = (jrb13 ? sqrt(pow(bal13*bal13e,2)+pow(fsr13*fsr13e,2))/jrb13 : 0);
    h1jrb13->SetBinContent(i, jrb13);
    h1jrb13->SetBinError(i, jrb13e);    
  } // for i
  
  for (int i = 1; i != h2jer->GetNbinsX()+1; ++i) {

    // Calculate JER2 (average of tag and probe)
    double eta = h2jer->GetXaxis()->GetBinCenter(i);
    TH1D *h1jer2(0), *h1m0s(0), *h1m0xs(0);
    h1jer2 = getJER(p2m0,p2m0x,eta,eta,&h1m0s,&h1m0xs);

    // Alternative longer calculation
    TH1D *h1bal2(0), *h1m2s(0), *h1m2xs(0);
    h1bal2 = getJER(p2m2,p2m2x,eta,eta,&h1m2s,&h1m2xs);
    TH1D *h1fsr2(0), *h1mos(0), *h1moxs(0);
    h1fsr2 = getJER(p2mo,p2mox,eta,eta,&h1mos,&h1moxs);
    
    for (int j = 1; j != h2jer->GetNbinsY()+1; ++j) {

      double jer2 = h1jer2->GetBinContent(j);
      double err2 = h1jer2->GetBinError(j);
      h2jer2->SetBinContent(i, j, jer2);
      h2jer2->SetBinError(i, j, err2);
      
      h2m0s->SetBinContent(i, j, h1m0s->GetBinContent(j));
      h2m0s->SetBinError(i, j, h1m0s->GetBinError(j));
      h2m0xs->SetBinContent(i, j, h1m0xs->GetBinContent(j));
      h2m0xs->SetBinError(i, j, h1m0xs->GetBinError(j));

      // Calculate probe JER by subtracting tag JER (FE, forward extension)
      double jer13 = h1jer13->GetBinContent(j);
      double err13 = h1jer13->GetBinError(j);
      double jer = sqrt(max(2.*jer2*jer2 - jer13*jer13,0.));
      double err = (jer ? sqrt(2.*pow(jer2*err2,2)+pow(jer13*err13,2))/jer : 0);
      h2jer->SetBinContent(i, j, jer);
      h2jer->SetBinError(i, j, err);

      //////////////////////////////////////////////////////////////////////
      
      // Alternative longer calculation
      double bal2 = h1bal2->GetBinContent(j);
      double bal2e = h1bal2->GetBinError(j);
      h2bal2->SetBinContent(i, j, bal2);
      h2bal2->SetBinError(i, j, bal2e);

      h2m2s->SetBinContent(i, j, h1m2s->GetBinContent(j));
      h2m2s->SetBinError(i, j, h1m2s->GetBinError(j));
      h2m2xs->SetBinContent(i, j, h1m2xs->GetBinContent(j));
      h2m2xs->SetBinError(i, j, h1m2xs->GetBinError(j));

      double fsr2 = h1fsr2->GetBinContent(j);
      double fsr2e = h1fsr2->GetBinError(j);
      h2fsr2->SetBinContent(i, j, fsr2);
      h2fsr2->SetBinError(i, j, fsr2e);

      h2mos->SetBinContent(i, j, h1mos->GetBinContent(j));
      h2mos->SetBinError(i, j, h1mos->GetBinError(j));
      h2moxs->SetBinContent(i, j, h1moxs->GetBinContent(j));
      h2moxs->SetBinError(i, j, h1moxs->GetBinError(j));

      //double jrb13 = h1jrb13->GetBinContent(j);
      //double jrbe13 = h1jrb13->GetBinError(j);
      double jrb2 = sqrt(max(bal2*bal2 - fsr2*fsr2,0.));
      double jrb2e = (jrb2 ? sqrt(pow(bal2*bal2e,2)+pow(fsr2*fsr2e,2))/jrb2 : 0);
      h2jrb2->SetBinContent(i, j, jrb2);
      h2jrb2->SetBinError(i, j, jrb2e);

      // Calculate B, F and JERB with forward extension
      double bal13 = h1bal13->GetBinContent(j);
      double bal13e = h1bal13->GetBinError(j);
      double bal = sqrt(max(2.*bal2*bal2 - bal13*bal13,0.));
      double bale = (bal ? sqrt(2.*pow(bal2*bal2e,2)+pow(bal13*bal13e,2))/bal : 0);
      h2bal->SetBinContent(i, j, bal);
      h2bal->SetBinError(i, j, bale);
      
      double fsr13 = h1fsr13->GetBinContent(j);
      double fsr13e = h1fsr13->GetBinError(j);
      double fsr = sqrt(max(2.*fsr2*fsr2 - fsr13*fsr13,0.));
      double fsre = (fsr ? sqrt(2.*pow(fsr2*fsr2e,2)+pow(fsr13*fsr13e,2))/fsr : 0);
      h2fsr->SetBinContent(i, j, fsr);
      h2fsr->SetBinError(i, j, fsre);
      
      //double jrb13 = h1jrb13->GetBinContent(j);
      //double jrb13e = h1jrb13->GetBinError(j);
      //double jrb = sqrt(max(1.*jrb2*jrb2 - jrb13*jrb13,0.));
      //double jrbe = (jrb ? sqrt(2.*pow(jrb2*jrb2e,2)+pow(jrb13*jrb13e,2))/jrb : 0);
      double jrb = sqrt(max(bal*bal - fsr*fsr,0.));
      double jrbe = (jrb ? sqrt(pow(bal*bale,2)+pow(fsr*fsre,2))/jrb : 0); 
      h2jrb->SetBinContent(i, j, jrb);
      h2jrb->SetBinError(i, j, jrbe);

    } // for j
    
    delete h1jer2;
    delete h1m0s;
    delete h1m0xs;

    delete h1bal2;
    delete h1m2s;
    delete h1m2xs;

    delete h1fsr2;
    delete h1mos;
    delete h1moxs;
  } // for i

  d->cd();

  // Store core results for MPFX only
  TH1D *h1mpfx13 = (TH1D*)h1jer13->Clone("h1mpfx13");
  TH2D *h2mpfx = (TH2D*)h2jer->Clone("h2mpfx");
  h1mpfx13->Write("h1mpfx13",TObject::kOverwrite);
  h2mpfx->Write("h2mpfx",TObject::kOverwrite);
  
  // Extras for MPFX only
  h1bal13->Write("h1bal13",TObject::kOverwrite);
  h1fsr13->Write("h1fsr13",TObject::kOverwrite);
  h1jrb13->Write("h1jrb13",TObject::kOverwrite);

  h2bal->Write("h2bal",TObject::kOverwrite);
  h2fsr->Write("h2fsr",TObject::kOverwrite);
  h2jrb->Write("h2jrb",TObject::kOverwrite);
  
  h2jer2->Write("h2mpfx2",TObject::kOverwrite);
  h2bal2->Write("h2bal2",TObject::kOverwrite);
  h2fsr2->Write("h2fsr2",TObject::kOverwrite);
  h2jrb2->Write("h2jrb2",TObject::kOverwrite);

  h1m0s13->Write(h1m0s13->GetName(),TObject::kOverwrite);
  h1m0xs13->Write(h1m0xs13->GetName(),TObject::kOverwrite);
  h1m2s13->Write(h1m2s13->GetName(),TObject::kOverwrite);
  h1m2xs13->Write(h1m2xs13->GetName(),TObject::kOverwrite);
  h1mos13->Write(h1mos13->GetName(),TObject::kOverwrite);
  h1moxs13->Write(h1moxs13->GetName(),TObject::kOverwrite);

  h2m0s->Write(h2m0s->GetName(),TObject::kOverwrite);
  h2m0xs->Write(h2m0xs->GetName(),TObject::kOverwrite);
  h2m2s->Write(h2m2s->GetName(),TObject::kOverwrite);
  h2m2xs->Write(h2m2xs->GetName(),TObject::kOverwrite);
  h2mos->Write(h2mos->GetName(),TObject::kOverwrite);
  h2moxs->Write(h2moxs->GetName(),TObject::kOverwrite);

  //f->Write()
  delete h2jer2;
  delete h2bal2;
  delete h2fsr2;
  delete h2jrb2;
  curdir->cd();

  // Continue to add PU noise from random cone measurements
  TH2D *h2n(0);
  if (true) { // adding RC noise

    // Only 2016 file for now, need to add others later
    TFile *frc(0);
    if (s.Contains("UL2016BCDEF") || s.Contains("UL2016APV")) {
      frc = new TFile("../JERCProtoLab/Summer20UL16APV/JER_noise/RC_noise_UL16APV.root","READ");
    }
    if (s.Contains("UL2016GH") || s.Contains("UL2016MG")) {
      frc = new TFile("../JERCProtoLab/Summer20UL16/JER_noise/RC_noise_UL16nonAPV.root","READ");
    }
    if (s.Contains("UL2017")) {
      frc = new TFile("../JERCProtoLab/Summer20UL17/JER_noise/RC_noise_UL17.root","READ");
    }
    if (s.Contains("UL2018")) {
      frc = new TFile("../JERCProtoLab/Summer20UL18/JER_noise/RC_noise_UL18.root","READ");
    }
    if (s.Contains("Run2")) {
      // Use UL2018 until proper Run2 combo available
      frc = new TFile("../JERCProtoLab/Summer20UL18/JER_noise/RC_noise_UL18.root","READ");
    }
    
    // Run3 2022+2023
    if (!frc && (s.Contains("202") || s.Contains("Summer22"))) {
      frc = new TFile("../JERCProtoLab/Summer20UL18/JER_noise/RC_noise_UL18.root","READ");
    }


    assert(frc && !frc->IsZombie());
    curdir->cd();

    TGraphAsymmErrors *g(0);
    if (isMC) g = (TGraphAsymmErrors*)frc->Get("rc_noiseterm_jer_MC_nominal");
    else      g = (TGraphAsymmErrors*)frc->Get("rc_noiseterm_jer_Data_nominal");
    assert(g);

    // Add RC for |eta|<1.3 region
    double sumrms2(0), sumeta(0);
    for (int i = 0; i != 5; ++i) {
      double deta = g->GetErrorXlow(i) + g->GetErrorXhigh(i);
      sumrms2 += deta*pow(g->GetY()[i],2);
      sumeta += deta;
    }
    double rms13 = (sumeta ? sqrt(sumrms2 / sumeta) : 0);

    // Update JER for |eta|<1.3
    TH1D *h1n13 = (TH1D*)h1jer13->Clone("h1n13"); h1n13->Reset();
    for (int i = 1; i != h1jer13->GetNbinsX()+1; ++i) {

      double pt = h1jer13->GetBinCenter(i);
      double n = rms13 / pt;
      double sc = h1jer13->GetBinContent(i);
      double esc = h1jer13->GetBinError(i);
      double jer = sqrt(n*n + sc*sc);
      double ejer = (jer!=0 ? sc * esc / jer : 0); // en=0
      h1n13->SetBinContent(i, n); 
      if (sc!=0 && esc!=0) {
	h1jer13->SetBinContent(i, jer);
	h1jer13->SetBinError(i, ejer);
      }
    }

    // Update JER for other eta bins
    h2n = (TH2D*)h2jer->Clone("h2n"); h2n->Reset();
    for (int i = 1; i != h2jer->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2jer->GetNbinsY()+1; ++j) {

	double etamin = h2jer->GetXaxis()->GetBinLowEdge(i);
	double etamax = h2jer->GetXaxis()->GetBinLowEdge(i+1);
	double absetamin = min(fabs(etamin),fabs(etamax));
	double absetamax = max(fabs(etamin),fabs(etamax));
	double sumrms2(0), sumeta(0);
	for (int k = 0; k != g->GetN(); ++k) {
	  if (absetamin <= g->GetX()[k] && g->GetX()[k] <= absetamax) {
	    double deta = g->GetErrorXlow(k) + g->GetErrorXhigh(k);
	    sumrms2 += deta*pow(g->GetY()[k],2);
	    sumeta += deta;
	  }
	}
	double rms = (sumeta ? sqrt(sumrms2 / sumeta) : 0.);

	double pt = h2jer->GetYaxis()->GetBinCenter(j);
	double n = rms / pt;
	double sc = h2jer->GetBinContent(i, j);
	double esc = h2jer->GetBinError(i, j);
	double jer = sqrt(n*n + sc*sc);
	double ejer = (jer!=0 ? sc * esc / jer : 0); // en=0
	h2n->SetBinContent(i, j, n); 
	if (sc!=0 && esc!=0) {
	  h2jer->SetBinContent(i, j, jer);
	  h2jer->SetBinError(i, j, ejer);
	}
      } // for j
    } // for i
  } // add RC noise

  d->cd();

  // Store core results for full JER
  h1jer13->Write(h1jer13->GetName(),TObject::kOverwrite);
  h2jer->Write(h2jer->GetName(),TObject::kOverwrite);
  if (h2n) h2n->Write(h2n->GetName(),TObject::kOverwrite);
  
  f->Write(nullptr,TObject::kOverwrite);
  curdir->cd();  

} // DijetHistosJERs


void drawDijetHistosJER(string sd, string sm, string era) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/jmenano_data_cmb_v21ul16.root","READ");
  //TFile *f = new TFile("rootfiles/jmenano_data_cmb_v22ul16.root","READ");
  TFile *f = new TFile(sd.c_str(),"READ");
  assert(f && !f->IsZombie());

  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","READ");
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v22ul16mg.root","READ");
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v23ul16mg.root","READ");
  TFile *fm = new TFile(sm.c_str(),"READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  TString tera = era.c_str();
  if (tera.Contains("UL2016APV")) lumi_13TeV = "2016early, 19.7 fb^{-1}";
  if (tera.Contains("UL2016BCDEF")) lumi_13TeV = "2016BCDEF, 19.7 fb^{-1}";
  if (tera.Contains("UL2016GH")) lumi_13TeV = "2016late, 16.8 fb^{-1}";
  if (tera.Contains("UL2017")) lumi_13TeV = "2017, 41.5 fb^{-1}";
  if (tera.Contains("UL2018")) lumi_13TeV = "2018, 59.9 fb^{-1}";
  if (tera.Contains("Run2")) lumi_13TeV = "2018, 137.9 fb^{-1}";

  const char *c = "_Dijet2";
  const char *cera = era.c_str();
  TCanvas *c1 = new TCanvas(Form("c1_%s_%s",c,cera),Form("c1_%s_%s",c,cera),
			    1200,600);
  c1->Divide(6,3,0,0);
  
  TCanvas *c2 = new TCanvas(Form("c2_%s_%s",c,cera),Form("c2_%s_%s",c,cera),
			    1200,600);
  c2->Divide(6,3,0,0);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  
  TH1D *h1jer13 =  (TH1D*)f->Get("Dijet2/h1jer13"); assert(h1jer13);
  TH1D *h1jer13m =  (TH1D*)fm->Get("Dijet2/h1jer13"); assert(h1jer13m);
  
  TH2D *h2jer =  (TH2D*)f->Get("Dijet2/h2jer"); assert(h2jer);
  TH2D *h2jerm =  (TH2D*)fm->Get("Dijet2/h2jer"); assert(h2jerm);

  TH2D *h2m0s =  (TH2D*)f->Get("Dijet2/h2m0s"); assert(h2m0s);
  TH2D *h2m0sm =  (TH2D*)fm->Get("Dijet2/h2m0s"); assert(h2m0sm);
  TH2D *h2m0xs =  (TH2D*)f->Get("Dijet2/h2m0xs"); assert(h2m0xs);
  TH2D *h2m0xsm =  (TH2D*)fm->Get("Dijet2/h2m0xs"); assert(h2m0xsm);

  for (int i = 1; i != h2jer->GetNbinsX()+1; ++i) {

    c1->cd(min(i,18));
    gPad->SetLogx();

    TH1D *h = tdrHist(Form("h%d_%s_%s",i,c,cera),"JER",0,0.3);
    tdrDraw(h,"",kNone);

    // RMS(MPFX) as low pT reference (sensitive to PU)
    TH1D *h1xs = h2m0xs->ProjectionY(Form("h1xs%d_%s_%s",i,c,cera),i,i);
    TH1D *h1xsm = h2m0xsm->ProjectionY(Form("h1xsm%d_%s_%s",i,c,cera),i,i);
    // Calculate effective parallel jet areas for MET (nmet):
    // int_0_2pi 1./(2pi)*cos^2(theta)dtheta = 1/2 to scale phi direction area
    double nmet = ((1./2.)*2.*5.*TMath::TwoPi()) / (TMath::Pi()*0.4*0.4);
    h1xs->Scale(1./sqrt(nmet));
    h1xsm->Scale(1./sqrt(nmet));
    h1xs->SetLineWidth(2);
    h1xs->GetXaxis()->SetRangeUser(15,200);
    //tdrDraw(h1xs,"HIST",kNone,kRed,kSolid,-1,kNone);
    h1xsm->GetXaxis()->SetRangeUser(15,200);
    //tdrDraw(h1xsm,"HIST",kNone,kRed,kSolid,-1,kNone);
    
    //h1jer13->SetLineWidth(2);
    //tdrDraw(h1jer13,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
    //tdrDraw(h1jer13m,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(h1jer13,"HIST",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(h1jer13m,"HIST",kNone,kBlue-9,kSolid,-1,kNone);

    TH1D *h1jer = h2jer->ProjectionY(Form("h1jer%d_%s_%s",i,c,cera),i,i);
    tdrDraw(h1jer,"Pz",kFullCircle,kGreen+2);
    h1jer->SetMarkerSize(0.7);
    TH1D *h1jerm = h2jerm->ProjectionY(Form("h1jerm%d_%s_%s",i,c,cera),i,i);
    tdrDraw(h1jerm,"Pz",kOpenCircle,kGreen+2);
    h1jerm->SetMarkerSize(0.7);

    double eta1 = h2jer->GetXaxis()->GetBinLowEdge(i);
    double eta2 = h2jer->GetXaxis()->GetBinLowEdge(i+1);
    tex->DrawLatex(0.60,0.88,Form("%1.3f < |#eta| < %1.3f",eta1,eta2));

    double ptmax = min(2000.,6500.*0.7/cosh(eta1));
    TF1 *f1 = new TF1(Form("f1%d_%s_%s",i,c,cera),"sqrt([0]*fabs([0])/(x*x)+"
		      "[1]*[1]*pow(x,[3])+[2]*[2])",30,ptmax);
    //f1->SetParameters(-1,1,0.04,0.);
    f1->SetParameters(0,1,0.04,-1);
    //f1->SetParLimits(0,-5.,5);
    f1->FixParameter(0,0);
    f1->SetParLimits(1,0.5,1.5);
    f1->SetParLimits(2,0.02,0.15);
    //f1->SetParLimits(3,-1.5,-0.5);
    f1->SetParLimits(3,-1.25,-0.75);
    if (min(fabs(eta1),fabs(eta2))>2.6) {
      f1->FixParameter(3,-1);
      //f1->SetParameter(0,0.);
      //f1->SetParLimits(0,0.,10.);
    }
    h1jerm->Fit(f1,"QRN");
    f1->SetRange(15,3500);
    f1->SetLineColor(kBlack);//kGreen+1);
    f1->Draw("SAME");

    // PATCH Run3
    h1jer->Scale(0.5);
    h1jer13->Scale(0.5);
    // \PATCH Run3
    
    if (i%6==0 || i==1) {
      //TLegend *leg = tdrLeg(0.70,0.85-4*0.05,1.00,0.85);
      TLegend *leg = tdrLeg(0.70,0.85-3*0.05,1.00,0.85);
      leg->AddEntry(h1jer,"DATA","PLE");
      leg->AddEntry(h1jerm,"MC","PLE");
      leg->AddEntry(h1jer13,"|#eta|<1.3","L");
      //leg->AddEntry(h1xs,"Est. PU","L");
    }
    
    c2->cd(min(i,18));
    gPad->SetLogx();
    
    TH1D *h2 = tdrHist(Form("h2_%d_%s_%s",i,c,cera),"Data/MC",0.8,2.0);
    tdrDraw(h2,"",kNone);
    l->DrawLine(15,1,3500,1);

    // RMS(MPF) as reference (less statistical uncertainty)
    //TH1D *h1m0s = h2m0s->ProjectionY(Form("h1m0s%d%s",i,c),i,i);
    //TH1D *h1m0sm = h2m0sm->ProjectionY(Form("h1m0sm%d%s",i,c),i,i);
    //TH1D *h1m0sr = (TH1D*)h1m0s->Clone(Form("h1m0sr%d%s",i,c));
    //h1m0sr->Divide(h1m0sm);
    //tdrDraw(h1m0sr,"HIST",kNone,kBlue,kSolid,-1,kNone);

    // RMS(MPFX) as low pT reference (sensitive to PU)
    TH1D *h1m0xs = h2m0xs->ProjectionY(Form("h1m0xs%d_%s_%s",i,c,cera),i,i);
    TH1D *h1m0xsm = h2m0xsm->ProjectionY(Form("h1m0xsm%d_%s_%s",i,c,cera),i,i);
    TH1D *h1m0xsr = (TH1D*)h1m0xs->Clone(Form("h1m0xsr%d_%s_%s",i,c,cera));
    h1m0xsr->Divide(h1m0xsm);
    h1m0xsr->GetXaxis()->SetRangeUser(15,200);
    //tdrDraw(h1m0xsr,"HIST",kNone,kRed,kSolid,-1,kNone);

    
    TH1D *h1jer13r = (TH1D*)h1jer13->Clone("h1jer13r");
    h1jer13r->Divide(h1jer13m);
    //h1jer13r->SetLineWidth(2);
    //tdrDraw(h1jer13r,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(h1jer13r,"HIST",kNone,kBlue,kSolid,-1,kNone);
    
    TH1D *h1jerr = (TH1D*)h1jer->Clone("h1jerr");
    h1jerr->Divide(h1jerm);
    tdrDraw(h1jerr,"Pz",kFullCircle,kGreen+2);
    
    tex->DrawLatex(0.60,0.88,Form("%1.3f < |#eta| < %1.3f",eta1,eta2));

    
    TF1 *f2 = new TF1(Form("f2%d_%s_%s",i,c,cera),
		      //"sqrt([0]*[4]*fabs([0]*[4])/(x*x)+"
		      //"[1]*[5]*[1]*[5]*pow(x,[3]*[7])+[2]*[6]*[2]*[6])"
		      "sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])"
		      "/sqrt([4]*fabs([4])/(x*x)+[5]*[5]*pow(x,[7])+[6]*[6])",
		      30,ptmax);
    //f2->SetParameters(0,1,0.04,-1,f1->GetParameter(0),f1->GetParameter(1),
    //		      f1->GetParameter(2),f1->GetParameter(3));
    for (int ip = 0; ip != f1->GetNpar(); ++ip) {
      //f1->SetParameter(ip, 1.);//f1->GetParameter(ip));
      f2->SetParameter(ip, f1->GetParameter(ip));
      f2->SetParameter(4+ip, f1->GetParameter(ip));
      f2->FixParameter(4+ip, f1->GetParameter(ip));
    }
    //f2->SetParLimits(0,-5.,5);
    //f2->FixParameter(0,1.);//f1->GetParameter(0));
    f2->FixParameter(0,f1->GetParameter(0));
    //f2->SetParLimits(1,0.5,1.5);
    //f2->SetParLimits(1,1,1.5);//f1->GetParameter(1),1.5);
    f2->SetParLimits(1,f1->GetParameter(1),1.5);
    //f2->SetParLimits(2,0.02,0.15);
    //f2->SetParLimits(2,1,1.5);//f1->GetParameter(2),0.15);
    f2->SetParLimits(2,f1->GetParameter(2),0.25);
    //f2->SetParLimits(3,-1.5,-0.5);
    //f2->FixParameter(3,1.);//f1->GetParameter(3));
    f2->FixParameter(3,f1->GetParameter(3));
    if (min(fabs(eta1),fabs(eta2))>2.6) {
      //f2->FixParameter(3,1.);//-1);
      f2->FixParameter(3,-1);
      //f2->SetParameter(0,0.);
      //f2->SetParLimits(0,0.,10.);
    }
    h1jerr->Fit(f2,"QRN");
    //f2->SetRange(15,3500);
    f2->SetLineWidth(2);
    f2->SetLineColor(kBlack);//Green+1);
    f2->Draw("SAME");

    // Add HCAL prefire pT^2 shape to JER as well
    TF1 *f2x2 = new TF1(Form("f2x2%d_%s_%s",i,c,cera),
			"sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2]+"
			//"[4]*[4]*x*x)"
			"[4]*[4]*pow(x,2))"
			"/sqrt([5]*fabs([5])/(x*x)+[6]*[6]*pow(x,[8])+[7]*[7])",
			30,ptmax);
    assert(f1->GetNpar()==4);
    f2x2->SetParameter(4, 0.);
    for (int ip = 0; ip != f1->GetNpar(); ++ip) {
      f2x2->SetParameter(ip, f1->GetParameter(ip));
      f2x2->SetParameter(5+ip, f1->GetParameter(ip));
      f2x2->FixParameter(5+ip, f1->GetParameter(ip));
    }
    f2x2->FixParameter(0,f1->GetParameter(0));
    f2x2->SetParLimits(1,f1->GetParameter(1),1.5);
    f2x2->SetParLimits(2,f1->GetParameter(2),0.25);
    f2x2->FixParameter(3,f1->GetParameter(3));
    if (min(fabs(eta1),fabs(eta2))>1.3) {
      f2x2->FixParameter(4, 0.);
    }
    if (min(fabs(eta1),fabs(eta2))>2.6) {
      f2x2->FixParameter(3,-1);
    }
    h1jerr->Fit(f2x2,"QRN");
    f2x2->SetLineWidth(1);
    f2x2->SetLineColor(kRed);
    f2x2->Draw("SAME");
    
    if (i%6==0 || i==1) {
      //TLegend *leg = tdrLeg(0.70,0.85-4*0.05,1.00,0.85);
      TLegend *leg = tdrLeg(0.70,0.85-3*0.05,1.00,0.85);
      leg->AddEntry(h1jer,"DATA","PLE");
      leg->AddEntry(h1jerm,"MC","PLE");
      leg->AddEntry(h1jer13,"|#eta|<1.3","L");
      //leg->AddEntry(h1m0xsr,"Est. PU","L");
    }
  } // for i

  c1->SaveAs(Form("pdf/DijetHistosJER_JERMC_%s.pdf",era.c_str()));
  c2->SaveAs(Form("pdf/DijetHistosJER_JERSF_%s.pdf",era.c_str()));
} // drawDijetHistosJER

void drawDijetHistosJERtest() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/jmenano_data_cmb_v22ul16.root","READ");
  //TFile *f = new TFile("rootfiles/jmenano_data_cmb_v21ul16.root","READ");
  //TFile *f = new TFile("rootfiles/jmenano_data_cmb.root","READ");
  assert(f && !f->IsZombie());

  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v23ul16mg.root","READ"); // with sF
  TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v22ul16mg.root","READ");// no SF
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","READ");
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","READ");
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb.root","READ");
  assert(fm && !fm->IsZombie());

  //TFile *fp = new TFile("rootfiles/jmenano_mc_cmb_v23ul16flat.root","READ"); // with SF
  TFile *fp = new TFile("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","READ"); // no SF
  assert(fp && !fp->IsZombie());


  TProfile2D *p2m0 = (TProfile2D*)f->Get("Dijet2/p2m0"); assert(p2m0);
  TProfile2D *p2m0x = (TProfile2D*)f->Get("Dijet2/p2m0x"); assert(p2m0x);
  TProfile2D *p2m0m = (TProfile2D*)fm->Get("Dijet2/p2m0"); assert(p2m0m);
  TProfile2D *p2m0xm = (TProfile2D*)fm->Get("Dijet2/p2m0x"); assert(p2m0xm);
  TProfile2D *p2m0p = (TProfile2D*)fp->Get("Dijet2/p2m0"); assert(p2m0p);
  TProfile2D *p2m0xp = (TProfile2D*)fp->Get("Dijet2/p2m0x"); assert(p2m0xp);

  TProfile2D *p2m2 = (TProfile2D*)f->Get("Dijet2/p2m2"); assert(p2m2);
  TProfile2D *p2m2x = (TProfile2D*)f->Get("Dijet2/p2m2x"); assert(p2m2x);
  TProfile2D *p2mnu = (TProfile2D*)f->Get("Dijet2/p2mnu"); assert(p2mnu);
  TProfile2D *p2mnux = (TProfile2D*)f->Get("Dijet2/p2mnux"); assert(p2mnux);
  
  TProfile2D *p2m2m = (TProfile2D*)fm->Get("Dijet2/p2m2"); assert(p2m2m);
  TProfile2D *p2m2xm = (TProfile2D*)fm->Get("Dijet2/p2m2x"); assert(p2m2xm);
  TProfile2D *p2mnum = (TProfile2D*)fm->Get("Dijet2/p2mnu"); assert(p2mnum);
  TProfile2D *p2mnuxm = (TProfile2D*)fm->Get("Dijet2/p2mnux"); assert(p2mnuxm);

  TProfile2D *p2m2p = (TProfile2D*)fp->Get("Dijet2/p2m2"); assert(p2m2p);
  TProfile2D *p2m2xp = (TProfile2D*)fp->Get("Dijet2/p2m2x"); assert(p2m2xp);
  TProfile2D *p2mnup = (TProfile2D*)fp->Get("Dijet2/p2mnu"); assert(p2mnup);
  TProfile2D *p2mnuxp = (TProfile2D*)fp->Get("Dijet2/p2mnux"); assert(p2mnuxp);
  curdir->cd();

  lumi_13TeV = "2016GH, 16.8 fb^{-1}";
  TH1D *h = tdrHist("h","JER",0,0.85,"p_{T,avb} (GeV)");
  TH1D *hd = tdrHist("hd","Data/MC",0.85,1.3,"p_{T,avb} (GeV)");
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,4,11);//,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  
  c1->cd(1);
  gPad->SetLogx();

  //double etamin(0.), etamax(1.3);
  double etamin(0.261), etamax(0.522);
  TH1D *h1jer(0), *h1m0s(0), *h1m0xs(0);
  h1jer = getJER(p2m0,p2m0x,etamin,etamax,&h1m0s,&h1m0xs);
  TH1D *h1jerm(0), *h1m0sm(0), *h1m0xsm(0);
  p2m0m->SetName("p2m0m");
  h1jerm = getJER(p2m0m,p2m0xm,etamin,etamax,&h1m0sm,&h1m0xsm);
  TH1D *h1jerp(0), *h1m0sp(0), *h1m0xsp(0);
  p2m0p->SetName("p2m0p");
  h1jerp = getJER(p2m0p,p2m0xp,etamin,etamax,&h1m0sp,&h1m0xsp);

  TH1D *h1bal(0), *h1m2s(0), *h1m2xs(0);
  h1bal = getJER(p2m2,p2m2x,etamin,etamax,&h1m2s,&h1m2xs);
  TH1D *h1balm(0), *h1m2sm(0), *h1m2xsm(0);
  p2m2m->SetName("p2m2m");
  h1balm = getJER(p2m2m,p2m2xm,etamin,etamax,&h1m2sm,&h1m2xsm);
  TH1D *h1balp(0), *h1m2sp(0), *h1m2xsp(0);
  p2m2p->SetName("p2m2p");
  h1balp = getJER(p2m2p,p2m2xp,etamin,etamax,&h1m2sp,&h1m2xsp);

  TH1D *h1fsr(0), *h1mnus(0), *h1mnuxs(0);
  h1fsr = getJER(p2mnu,p2mnux,etamin,etamax,&h1mnus,&h1mnuxs);
  TH1D *h1fsrm(0), *h1mnusm(0), *h1mnuxsm(0);
  p2mnum->SetName("p2mnum");
  h1fsrm = getJER(p2mnum,p2mnuxm,etamin,etamax,&h1mnusm,&h1mnuxsm);
  TH1D *h1fsrp(0), *h1mnusp(0), *h1mnuxsp(0);
  p2mnup->SetName("p2mnup");
  h1fsrp = getJER(p2mnup,p2mnuxp,etamin,etamax,&h1mnusp,&h1mnuxsp);

  // Extract JER from MPF2 by subtracting FSR = MPFNU-MPFNUX
  TH1D *h1jer2 = (TH1D*)h1m2s->Clone("h1jer2");
  TH1D *h1jer2m = (TH1D*)h1m2sm->Clone("h1jer2m");
  TH1D *h1jer2p = (TH1D*)h1m2sp->Clone("h1jer2p");
  for (int i = 1; i != h1jer2->GetNbinsX()+1; ++i) {
    //double bal = h1m2sm->GetBinContent(i);
    double bal = h1bal->GetBinContent(i);
    double fsr = h1fsr->GetBinContent(i);
    h1jer2->SetBinContent(i, sqrt(max(bal*bal - fsr*fsr,0.)));
    double balm = h1balm->GetBinContent(i);
    double fsrm = h1fsrm->GetBinContent(i);
    h1jer2m->SetBinContent(i, sqrt(max(balm*balm - fsrm*fsrm,0.)));
    double balp = h1balp->GetBinContent(i);
    double fsrp = h1fsrp->GetBinContent(i);
    h1jer2p->SetBinContent(i, sqrt(max(balp*balp - fsrp*fsrp,0.)));
  }
  
  TF1 *fjer = new TF1("fjer","sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2])",
		      30,3000);
  TF1 *fjerm = new TF1("fjerm","sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2])",
		       30,3000);
  TF1 *fjerp = new TF1("fjerp","sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2])",
		       30,3000);

  fjer->SetParameters(-1,1,0.05);
  h1jer->Fit(fjer,"QRN");
  fjer->SetRange(15,3500);

  fjerm->SetParameters(fjer->GetParameter(0),fjer->GetParameter(1),
		       fjer->GetParameter(2));
  h1jerm->Fit(fjerm,"QRN");
  fjerm->SetRange(15,3500);

  fjerp->SetParameters(fjer->GetParameter(0),fjer->GetParameter(1),
		       fjer->GetParameter(2));
  h1jerp->Fit(fjerm,"QRN");
  fjerp->SetRange(15,3500);

  fjer->SetLineWidth(2);
  fjer->SetLineColor(kGreen+2);
  fjer->Draw("SAME");
  fjerm->SetLineColor(kGreen+2);
  fjerm->Draw("SAME");
  fjerp->SetLineColor(kGreen-9);
  fjerp->Draw("SAME");

  tdrDraw(h1m0sp,"PE",kOpenCircle,kBlue-9);
  tdrDraw(h1m0xsp,"PE",kOpenDiamond,kRed-9);
  tdrDraw(h1jerp,"PE",kOpenCircle,kGreen-9);

  tdrDraw(h1m0sm,"PE",kOpenCircle,kBlue);
  tdrDraw(h1m0xsm,"PE",kOpenDiamond,kRed);
  tdrDraw(h1jerm,"PE",kOpenCircle,kGreen+2);
  
  tdrDraw(h1m0s,"PE",kFullCircle,kBlue);
  tdrDraw(h1m0xs,"PE",kFullDiamond,kRed);
  tdrDraw(h1jer,"PE",kFullCircle,kGreen+2);

  tdrDraw(h1balp,"PE",kOpenDiamond,kMagenta-9);
  tdrDraw(h1fsrp,"PE",kOpenStar,kOrange-9);
  tdrDraw(h1jer2p,"PE",kOpenSquare,kGreen-9);

  //tdrDraw(h1m2sm,"PE",kOpenDiamond,kMagenta+2);
  tdrDraw(h1balm,"PE",kOpenDiamond,kMagenta+2);
  //tdrDraw(h1mnusm,"PE",kOpenStar,kOrange+2);
  tdrDraw(h1fsrm,"PE",kOpenStar,kOrange+2);
  tdrDraw(h1jer2m,"PE",kOpenSquare,kGreen+3);

  //tdrDraw(h1mnuxsm,"PE",kFullDiamond,kBlack);

  tdrDraw(h1bal,"PE",kFullDiamond,kMagenta+2);
  tdrDraw(h1fsr,"PE",kFullStar,kOrange+2);
  tdrDraw(h1jer2,"PE",kFullSquare,kGreen+3); h1jer2->SetMarkerSize(0.7);

  
  TLegend *leg = tdrLeg(0.60,0.88-7*0.05,0.85,0.88);
  leg->SetHeader("  Data");

  //leg->AddEntry(h1m0s,"RMS(MPF)","PLE");
  //leg->AddEntry(h1m0xs,"RMS(MPFX)","PLE");
  leg->AddEntry(h1m0s,"MPF","PLE");
  leg->AddEntry(h1m0xs,"MPFX","PLE");
  /*
  // MPF2X can effectively only remove ISR
  leg->AddEntry(h1m0s,"RMS(MPF2)","PLE");
  leg->AddEntry(h1m0xs,"RMS(MPF2X)","PLE");
  */
  //leg->AddEntry(h1jer,"JER","PLE");
  leg->AddEntry(h1jer,"JER = M #oplus -MX","PLE");

  TLegend *legm = tdrLeg(0.53,0.88-7*0.05,0.78,0.88);
  legm->SetHeader("MG");
  legm->AddEntry(h1m0sm,"","PLE");
  legm->AddEntry(h1m0xsm,"","PLE");
  legm->AddEntry(h1jerm,"","PLE");
  TLegend *legp = tdrLeg(0.47,0.88-7*0.05,0.72,0.88);
  legp->SetHeader("P8");
  legp->AddEntry(h1m0sp,"","PLE");
  legp->AddEntry(h1m0xsp,"","PLE");
  legp->AddEntry(h1jerp,"","PLE");

  leg->AddEntry(h1bal,"B = M2 #oplus -M2X","PLE");
  legm->AddEntry(h1balm,"","PLE");
  legp->AddEntry(h1balp,"","PLE");
  leg->AddEntry(h1fsr,"F = MNU #oplus -MNUX","PLE");
  legm->AddEntry(h1fsrm,"","PLE");
  legp->AddEntry(h1fsrp,"","PLE");
  leg->AddEntry(h1jer2,"JER = B #oplus -F","PLE");
  legm->AddEntry(h1jer2m,"","PLE");
  legp->AddEntry(h1jer2p,"","PLE");

  //tex->DrawLatex(0.53,0.60,"|#eta| < 1.3");
  tex->DrawLatex(0.53,0.45,"|#eta| < 1.3");

  c1->cd(2);
  gPad->SetLogx();
  
  l->DrawLine(15,1,3500,1);

  TH1D *h1m0spr = (TH1D*)h1m0s->Clone("h1m0spr"); h1m0spr->Divide(h1m0sp);
  TH1D *h1m0xspr = (TH1D*)h1m0xs->Clone("h1m0xspr"); h1m0xspr->Divide(h1m0xsp);
  TH1D *h1jerpr = (TH1D*)h1jer->Clone("hjerpr"); h1jerpr->Divide(h1jerp);

  TH1D *h1m0sr = (TH1D*)h1m0s->Clone("h1m0sr"); h1m0sr->Divide(h1m0sm);
  TH1D *h1m0xsr = (TH1D*)h1m0xs->Clone("h1m0xsr"); h1m0xsr->Divide(h1m0xsm);
  TH1D *h1jerr = (TH1D*)h1jer->Clone("hjerr"); h1jerr->Divide(h1jerm);

  TH1D *h1fsrpr = (TH1D*)h1fsr->Clone("h1fsrpr"); h1fsrpr->Divide(h1fsrp);
  TH1D *h1jer2pr = (TH1D*)h1jer2->Clone("h1jer2pr"); h1jer2pr->Divide(h1jer2p);

  TH1D *h1fsrr = (TH1D*)h1fsr->Clone("h1fsrr"); h1fsrr->Divide(h1fsrm);
  TH1D *h1jer2r = (TH1D*)h1jer2->Clone("h1jer2r"); h1jer2r->Divide(h1jer2m);

  TF1 *fjerpr = new TF1("fjerpr","sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2]) /"
			"sqrt([3]*fabs([3])/(x*x)+[4]*[4]/x+[5]*[5])",15,3000);
  fjerpr->SetParameters(fjer->GetParameter(0),fjer->GetParameter(1),
			fjer->GetParameter(2),
			fjerp->GetParameter(0),fjerp->GetParameter(1),
			fjerp->GetParameter(2));
  fjerpr->SetLineWidth(1);
  fjerpr->SetLineColor(kGreen-9);
  fjerpr->Draw("SAME");
  
  TF1 *fjerr = new TF1("fjerr","sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2]) /"
		       "sqrt([3]*fabs([3])/(x*x)+[4]*[4]/x+[5]*[5])",15,3000);
  fjerr->SetParameters(fjer->GetParameter(0),fjer->GetParameter(1),
		       fjer->GetParameter(2),
		       fjerm->GetParameter(0),fjerm->GetParameter(1),
		       fjerm->GetParameter(2));
  fjerr->SetLineWidth(2);
  fjerr->SetLineColor(kGreen+2);
  fjerr->Draw("SAME");

  tdrDraw(h1m0spr,"PE",kFullCircle,kBlue-9);
  tdrDraw(h1m0xspr,"PE",kFullDiamond,kRed-9);
  tdrDraw(h1jerpr,"PE",kFullCircle,kGreen-9);

  tdrDraw(h1fsrpr,"PE",kFullStar,kOrange-9);
  tdrDraw(h1jer2pr,"PE",kFullSquare,kGreen-9); 
  
  tdrDraw(h1m0sr,"PE",kFullCircle,kBlue);
  tdrDraw(h1m0xsr,"PE",kFullDiamond,kRed);
  tdrDraw(h1jerr,"PE",kFullCircle,kGreen+2);

  tdrDraw(h1fsrr,"PE",kFullStar,kOrange+2);
  tdrDraw(h1jer2r,"PE",kFullSquare,kGreen+3); 

  c1->SaveAs("pdf/DijetHistosJER_JER13.pdf");

  TH1D *_h1jern(0), *_h1jernm(0), *_h1jernr(0);
  { // adding RC noise

    // This is probably 2018 file. Would need to ask other years as well
    //TFile *frc = new TFile("../jecsys2020/rootfiles/jerCombo/RC.root","READ");
    TFile *frc = new TFile("../JERCProtoLab/Summer20UL16/JER_noise/RC_noise_UL16nonAPV.root","READ");
    assert(frc && !frc->IsZombie());
    curdir->cd();

    //TGraphAsymmErrors *gd = (TGraphAsymmErrors*)frc->Get("Data/RMS");
    TGraphAsymmErrors *gd = (TGraphAsymmErrors*)frc->Get("rc_noiseterm_jer_Data_nominal");
    assert(gd);
    //TGraphAsymmErrors *gm = (TGraphAsymmErrors*)frc->Get("MC/RMS");
    TGraphAsymmErrors *gm = (TGraphAsymmErrors*)frc->Get("rc_noiseterm_jer_MC_nominal");
    assert(gm);

    double sumrmsdt2(0), sumrmsmc2(0), sumeta(0);
    for (int i = 0; i != 5; ++i) {
      double deta = 2.* gm->GetErrorXlow(i);
      //sumrmsdt2 += deta*pow(gd->GetY()[i],2);
      // PATCH Run3
      sumrmsdt2 += 0;
      // \PATCH Run3
      sumrmsmc2 += deta*pow(gm->GetY()[i],2);
      sumeta += deta;
    }
    double rmsdt = sqrt(sumrmsdt2 / sumeta);
    double rmsmc = sqrt(sumrmsmc2 / sumeta);

    TH1D *h1n = (TH1D*)h1jer->Clone("h1n"); h1n->Reset();
    TH1D *h1nm = (TH1D*)h1jer->Clone("h1nm"); h1nm->Reset();
    TH1D *h1nr = (TH1D*)h1jer->Clone("h1nr"); h1nr->Reset();
    TH1D *h1jern = (TH1D*)h1jer->Clone("h1jern"); h1jern->Reset();
    TH1D *h1jernm = (TH1D*)h1jer->Clone("h1jernm"); h1jernm->Reset();
    TH1D *h1jernr = (TH1D*)h1jer->Clone("h1jernr"); h1jernr->Reset();
    for (int i = 1; i != h1n->GetNbinsX()+1; ++i) {
      double n = rmsdt / h1n->GetBinCenter(i);
      double sc = h1jer->GetBinContent(i);
      double esc = h1jer->GetBinError(i);
      double jern = sqrt(n*n + sc*sc);
      double ejern = (jern ? sc * esc / jern : 0); // en=0
      h1n->SetBinContent(i, n); 
      h1jern->SetBinContent(i, jern);
      h1jern->SetBinError(i, ejern);

      double nm = rmsmc / h1nm->GetBinCenter(i);
      double scm = h1jerm->GetBinContent(i);
      double escm = h1jerm->GetBinError(i);
      double jernm = sqrt(nm*nm + scm*scm);
      double ejernm = (jernm ? scm * escm / jernm : 0); // enm=0
      h1nm->SetBinContent(i, nm); 
      h1jernm->SetBinContent(i, jernm);
      h1jernm->SetBinError(i, ejernm);

      h1nr->SetBinContent(i, nm ? n/nm : 0.); 
      h1jernr->SetBinContent(i, jernm ? jern/jernm : 0.);
      if (jern!=0 && jernm!=0) 
	h1jernr->SetBinError(i, sqrt(pow(ejern/jern,2) + pow(ejernm/jernm,2))*
			     jern / jernm);
    }

    //double neff = (2.*5.*TMath::TwoPi()) / (TMath::Pi()*0.4*0.4) / sqrt(2.);
    //TH1D *h1m0xsn = (TH1D*)h1m0xs->Clone("h1m0xsn");
    // undo 1./sqrt(2.) earlier, then scale detector resolution of energy
    // fluctuations (58%?) to local energy fluctuations
    //h1m0xsn->Scale(sqrt(2.)/sqrt(neff) / 0.58); 

    // Separate out PF low PU behavior with N<0 and extra from PU
    // TBD: add noise from UE, which is not there in ZeroBias and MPFX
    double ptmax = min(2000.,6500.*0.7/cosh(0.));
    TF1 *f1m = new TF1("f13m","sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
		       "[1]*[1]*pow(x,[3])+[2]*[2])",30,ptmax);
    f1m->SetParameters(0,1,0.04,-1,rmsmc);
    f1m->SetParLimits(0,-5.,0.);
    f1m->FixParameter(4,rmsmc);
    f1m->SetParLimits(1,0.5,1.5);
    f1m->SetParLimits(2,0.02,0.25);
    f1m->SetParLimits(3,-1.25,-0.75);
    h1jernm->Fit(f1m,"QRN");

    TF1 *f1 = new TF1("f13","sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
		      "[1]*[1]*pow(x,[3])+[2]*[2])",30,ptmax);
    f1->SetParameters(f1m->GetParameter(0),f1m->GetParameter(1),
		      f1m->GetParameter(2),f1m->GetParameter(3),rmsdt);
    f1->FixParameter(0,f1m->GetParameter(0));
    f1->FixParameter(4,rmsdt);
    f1->SetParLimits(1,f1m->GetParameter(1),1.5);
    f1->SetParLimits(2,f1m->GetParameter(2),0.25);
    f1->FixParameter(3,f1m->GetParameter(3));
    h1jern->Fit(f1,"QRN");

    TF1 *f1r = new TF1("f13r",
		       "sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
		       "[1]*[1]*pow(x,[3])+[2]*[2])/"
		       "sqrt([5]*fabs([5])/(x*x)+[9]*[9]/(x*x)+"
		       "[6]*[6]*pow(x,[8])+[7]*[7])",
		       30,ptmax);
    f1r->SetParameters(f1->GetParameter(0),f1->GetParameter(1),
		       f1->GetParameter(2),f1->GetParameter(3),
		       f1->GetParameter(4),
		       f1m->GetParameter(0),f1m->GetParameter(1),
		       f1m->GetParameter(2),f1m->GetParameter(3),
		       f1m->GetParameter(4));
    
    TH1D *h = tdrHist("h2","RMS",0,0.35,"p_{T,avb} (GeV)");
    TH1D *hd = tdrHist("hd2","Data/MC",0.90,1.25,"p_{T,avb} (GeV)");
    TCanvas *c2 = tdrDiCanvas("c2",h,hd,4,11);
    c2->cd(1);
    gPad->SetLogx();

    //tdrDraw(h1m0xsn,"Pz",kFullDiamond,kRed-9,kSolid,-1,kNone);
    tdrDraw(h1n,"Pz",kFullDiamond,kRed,kSolid,-1,kNone);
    tdrDraw(h1jer,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone);
    tdrDraw(h1jern,"Pz",kFullStar,kBlack,kSolid,-1,kNone);
    _h1jern = h1jern;
    _h1jernm = h1jernm;
    
    f1->SetRange(15,3500);
    f1->SetLineColor(kBlack);
    f1->Draw("SAME");

    TLegend *leg = tdrLeg(0.60,0.89-3*0.06,.90,0.89);
    leg->AddEntry(h1jern,"Full JER","PLE");
    leg->AddEntry(h1jer,"MPFX (N_{0}SCd)","PLE");
    leg->AddEntry(h1n,"RC (N_{PU})","PLE");

    tex->DrawLatex(0.45,0.60,Form("#chi^{2} / ndf (Data) = %1.1f / %d",
				  f1->GetChisquare(), f1->GetNDF()));
    tex->DrawLatex(0.45,0.54,Form("#chi^{2} / ndf (MC) = %1.1f / %d",
				  f1m->GetChisquare(), f1m->GetNDF()));

    f1r->SetRange(15,3500);
    double chi2(0); int ndf(-2);
    for (int i = 1; i != h1jernr->GetNbinsX()+1; ++i) {
      if (h1jernr->GetBinError(i)!=0 && h1jernr->GetBinContent(i)!=0) {
	double pt = h1jernr->GetBinCenter(i);
	chi2 += pow(h1jernr->GetBinContent(i)-f1r->Eval(pt),2) /
	  pow(h1jernr->GetBinError(i),2);
	++ndf;
      }
    }
    tex->DrawLatex(0.45,0.49,Form("#chi^{2} / ndf (Ratio) = %1.1f / %d",
				  chi2,ndf));
    tex->DrawLatex(0.42,0.80,"|#eta| < 1.3");

    c2->cd(2);
    gPad->SetLogx();

    tdrDraw(h1nr,"Pz",kFullDiamond,kRed,kSolid,-1,kNone);
    tdrDraw(h1jerr,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone);
    tdrDraw(h1jernr,"Pz",kFullStar,kBlack,kSolid,-1,kNone);
    _h1jernr = h1jernr;
    
    f1r->SetLineColor(kBlack);
    f1r->Draw("SAME");


    TF1 *fc = new TF1("fc","[0]",15,3500);
    h1jerr->Fit(fc,"QRN");
    fc->Draw("SAME");
    tex->DrawLatex(0.45,0.49,Form("#chi^{2} / ndf (const) = %1.1f / %d",
				  fc->GetChisquare(), fc->GetNDF()));
		      
    c2->SaveAs("pdf/DijetHistosJER_JER13_JERwithRC.pdf");
  } // RC noise
  
  // test with Z+jet and dijet traditional
  if (false) {
    TFile *fdj = new TFile("rootfiles/dijet_balance_UL18_Summer19UL18_V5_AK4CHS.root","READ");
    assert(fdj && !fdj->IsZombie());
    // older dijet files
    //TH1F *h1dj = (TH1F*)fdj->Get("Data_jer_dijet_SM_1_nominal"); assert(h1dj);
    //TH1F *h1djm = (TH1F*)fdj->Get("MC_jer_dijet_SM_1_nominal"); assert(h1djm);
    //TH1F *h1djr = (TH1F*)h1dj->Clone("h1djr");
    //h1djr->Divide(h1djm);
    TGraphAsymmErrors *gdj(0), *gdjm(0);
    gdj = (TGraphAsymmErrors*)fdj->Get("dijet_balance_jer_Data_0p261_0p522_SM_nominal"); assert(gdj);
    gdjm = (TGraphAsymmErrors*)fdj->Get("dijet_balance_jer_MC_0p261_0p522_SM_nominal"); assert(gdjm);
    TGraphErrors *gdjr = new TGraphErrors(gdj->GetN());
    for (int i = 0; i != gdjr->GetN(); ++i) {
      gdjr->SetPoint(i, gdj->GetX()[i], gdj->GetY()[i]/gdjm->GetY()[i]);
      gdjr->SetPointError(i, gdj->GetErrorXlow(i),
			  sqrt(pow(gdj->GetErrorYlow(i)/gdj->GetY()[i],2)+
			       pow(gdjm->GetErrorYlow(i)/gdjm->GetY()[i],2))
			  * gdjr->GetY()[i]);
    }

    TFile *fzj = new TFile("rootfiles/zjet_balance_UL2018_jetpt_nominal_small.root","READ");
    assert(fzj && !fzj->IsZombie());
    TGraphErrors *gz(0), *gzm(0), *gzr(0);
    gz = (TGraphErrors*)fzj->Get("data"); assert(gz);
    gzm = (TGraphErrors*)fzj->Get("mc"); assert(gzm);
    gzr = (TGraphErrors*)fzj->Get("ratio"); assert(gzr);
    
    
    //TFile *fz = new TFile("../jecsys2020/rootfiles/jme_ZplusJet_Muon_Run2016FGH_v8.root","READ");
    TFile *fz = new TFile("../jecsys2020/rootfiles/jme_bplusZ_2016FGH_Zmm_sync_v11.root","READ");
    assert(fz && !fz->IsZombie());

    TH2D *h2z(0), *h2zx(0), *h2zm(0), *h2zxm(0);
    //const char *cd = "Run2016F-H/2016/data/eta_00_13/";
    //const char *cm = "Run2016F-H/2016/mc/eta_00_13/";
    //h2z = (TH2D*)fz->Get(Form("%s/zpt_mpf_zmmjet_a100",cd)); assert(h2z);
    //h2zx = (TH2D*)fz->Get(Form("%s/zpt_mpfx_zmmjet_a100",cd)); assert(h2zx);
    //h2zm = (TH2D*)fz->Get(Form("%s/zpt_mpf_zmmjet_a100",cm)); assert(h2zm);
    //h2zxm = (TH2D*)fz->Get(Form("%s/zpt_mpfx_zmmjet_a100",cm)); assert(h2zxm);
    const char *cd = "data/eta_00_13/";
    const char *cm = "mc/eta_00_13/";
    h2z = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPF_alpha100",cd)); assert(h2z);
    h2zx = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPFx_alpha100",cd)); assert(h2zx);
    h2zm = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPF_alpha100",cm)); assert(h2zm);
    h2zxm = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPFx_alpha100",cm)); assert(h2zxm);
    //
    TH2D *h2z2(0), *h2z2x(0), *h2z2m(0), *h2z2xm(0);
    // NB: some bug with RMPFJet1 cutting of at 0.5? pT1/pTZ>0.5?
    h2z2 = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPFjet1_alpha100",cd)); assert(h2z2);
    h2z2x = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPFjet1x_alpha100",cd)); assert(h2z2x);
    h2z2m = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPFjet1_alpha100",cm)); assert(h2z2m);
    h2z2xm = (TH2D*)fz->Get(Form("%s/h_Zpt_RMPFjet1x_alpha100",cm)); assert(h2z2xm);
    curdir->cd();
    
    //TProfile *pz = h2z->ProfileX("pz",1,-1,"S");
    //TProfile *pzx = h2zx->ProfileX("pzx",1,-1,"S");
    //TProfile *pzm = h2zm->ProfileX("pz",1,-1,"S");
    //TProfile *pzxm = h2zxm->ProfileX("pzx",1,-1,"S");

    TH1D *h1jerz(0), *h1m0zs(0), *h1m0xzs(0);
    h1jerz = getJERZ(h2z,h2zx,"data",&h1m0zs,&h1m0xzs);
    TH1D *h1jerzm(0), *h1m0zsm(0), *h1m0xzsm(0);
    h1jerzm = getJERZ(h2zm,h2zxm,"mc",&h1m0zsm,&h1m0xzsm);

    TH1D *h1balz(0), *h1m2zs(0), *h1m2xzs(0);
    h1balz = getJERZ(h2z2,h2z2x,"data",&h1m2zs,&h1m2xzs);
    TH1D *h1balzm(0), *h1m2zsm(0), *h1m2xzsm(0);
    h1balzm = getJERZ(h2z2m,h2z2xm,"mc",&h1m2zsm,&h1m2xzsm);

    c1->cd(1);

    tdrDraw(h1m0xzsm,"PE",kOpenDiamond,kRed-9);
    tdrDraw(h1m0zsm,"PE",kOpenDiamond,kBlue-9);
    tdrDraw(h1m0xzs,"PE",kFullDiamond,kRed-9);
    tdrDraw(h1m0zs,"PE",kFullDiamond,kBlue-9);
    
    tdrDraw(h1jerz,"PE",kFullStar,kBlack);
    tdrDraw(h1jerzm,"PE",kOpenStar,kBlack);

    c1->cd(2);

    TH1D *h1m0zsr = (TH1D*)h1m0zs->Clone("h1m0zsr");
    h1m0zsr->Divide(h1m0zsm);
    TH1D *h1m0xzsr = (TH1D*)h1m0xzs->Clone("h1m0xzsr");
    h1m0xzsr->Divide(h1m0xzsm);
    TH1D *h1jerzr = (TH1D*)h1jerz->Clone("h1jerzr");
    h1jerzr->Divide(h1jerzm);
    //
    TH1D *h1m2zsr = (TH1D*)h1m2zs->Clone("h1m2zsr");
    h1m2zsr->Divide(h1m2zsm);
    TH1D *h1m2xzsr = (TH1D*)h1m2xzs->Clone("h1m2xzsr");
    h1m2xzsr->Divide(h1m2xzsm);
    TH1D *h1balzr = (TH1D*)h1balz->Clone("h1balzr");
    h1balzr->Divide(h1balzm);
    //
    TH1D *h1m2sr = (TH1D*)h1m2s->Clone("h1m2sr");
    h1m2sr->Divide(h1m2sm);
    TH1D *h1m2xsr = (TH1D*)h1m2xs->Clone("h1m2xsr");
    h1m2xsr->Divide(h1m2xsm);
    TH1D *h1balr = (TH1D*)h1bal->Clone("h1balr");
    h1balr->Divide(h1balm);

    tdrDraw(h1m0xzsr,"PE",kOpenDiamond,kRed-9);
    tdrDraw(h1m0zsr,"PE",kOpenDiamond,kBlue-9);
    tdrDraw(h1jerzr,"PE",kFullStar,kBlack);

    c1->SaveAs("pdf/DijetHistosJER_JER13_DijetAndZJet.pdf");
    c1->SaveAs("pdf/DijetHistosJER_JER13_DijetAndZJet.root");


    TH1D *hz = tdrHist("hz","RMS",0,0.35);//0.8);
    TH1D *hdz = tdrHist("hdz","Data/MC",0.85,1.5);
    TCanvas *c2z = tdrDiCanvas("c2z",hz,hdz,4,11);

    TH1D *h1balb = (TH1D*)h1m2s->Clone("h1balb");
    h1balb->Scale(1./sqrt(2.));
    TH1D *h1isrb = (TH1D*)h1m2xs->Clone("h1isrb");
    h1isrb->Scale(1./sqrt(2.));
    TH1D *h1jerb = (TH1D*)h1bal->Clone("h1jerb");
    h1jerb->Scale(1./sqrt(2.));
    
    c2z->cd(1);
    gPad->SetLogx();
    //tdrDraw(h1m2s,"PE",kFullCircle,kRed);
    //    tdrDraw(h1balb,"PE",kFullCircle,kRed);
    //tdrDraw(h1m2xs,"PE",kFullCircle,kBlue);
    //    tdrDraw(h1isrb,"PE",kFullCircle,kBlue);
    //tdrDraw(h1bal,"PE",kFullCircle,kMagenta+2);
    //    tdrDraw(h1jerb,"PE",kFullCircle,kMagenta+2);
    //tdrDraw(h1m2zs,"PE",kFullDiamond,kRed-9);
    //tdrDraw(h1m2xzs,"PE",kFullDiamond,kBlue-9);
    //tdrDraw(h1balz,"PE",kFullDiamond,kMagenta-7);//9);
    //tdrDraw(h1balzm,"PE",kOpenDiamond,kMagenta-7);//9);
    tdrDraw(h1jerz,"PE",kFullDiamond,kMagenta+1);
    h1jerz->GetXaxis()->SetRangeUser(15,400);
    //tdrDraw(h1jerzm,"PE",kOpenDiamond,kMagenta+1);
    
    tdrDraw(gz,"PEz",kFullStar,kRed);
    tdrDraw(gzm,"PEz",kOpenStar,kRed);
    //tdrDraw(h1dj,"PE",kFullStar,kBlack);
    tdrDraw(gdj,"PEz",kFullStar,kBlack);
    tdrDraw(gdjm,"PEz",kOpenStar,kBlack);
    //tdrDraw(h1jer,"PE",kFullCircle,kGreen+2);
    tdrDraw(_h1jern,"PE",kFullDiamond,kGreen+2);//3);
    tdrDraw(_h1jernm,"PE",kOpenDiamond,kGreen+2);//3);
    
    c2z->cd(2);
    gPad->SetLogx();
    l->DrawLine(15,1,3500,1);

    //    tdrDraw(h1m2sr,"PE",kFullCircle,kRed);
    //    tdrDraw(h1m2xsr,"PE",kFullCircle,kBlue);
    //    tdrDraw(h1balr,"PE",kFullCircle,kMagenta+2);
    //tdrDraw(h1m2zsr,"PE",kFullDiamond,kRed-9);
    //tdrDraw(h1m2xzsr,"PE",kFullDiamond,kBlue-9);
    //tdrDraw(h1balzr,"PE",kFullDiamond,kMagenta-7);//-9);
    tdrDraw(h1jerzr,"PE",kFullDiamond,kMagenta+1);
    h1jerzr->Scale(1./0.929586); // h1jer13 ratio for v22/v23
    h1jerzr->GetXaxis()->SetRangeUser(15,400);
    
    tdrDraw(gzr,"PEz",kFullStar,kRed);
    //tdrDraw(h1djr,"PE",kFullStar,kBlack);
    tdrDraw(gdjr,"PEz",kFullStar,kBlack);
    //tdrDraw(h1jerr,"PE",kFullCircle,kGreen+2);
    tdrDraw(_h1jernr,"PE",kFullDiamond,kGreen+2);//3);
    //tdrDraw(h1jerzr,"PE",kFullDiamond,kGreen+2);

    //c2z->SaveAs("pdf/DijetHistosJER_Bal13_DijetAndZJet.pdf");
    c2z->SaveAs("pdf/DijetHistosJER_Bal13_vsDijetJER.pdf");
  } // Z+jet


} // drawDijetHistosJERtest

// Calculate JER for one given eta slice at a time
TH1D *getJER(TProfile2D* p2, TProfile2D *p2x,
	     double eta1, double eta2, TH1D **h1, TH1D **h1x) {
  
  string s = Form("_%s",p2->GetName());
  const char *c = s.c_str();

  int i1 = p2->GetXaxis()->FindBin(eta1);
  int i2 = p2->GetXaxis()->FindBin(eta2);

  TProfile *p1 = p2->ProfileY(Form("p1%s",c),i1,i2,"");
  TProfile *p1s = p2->ProfileY(Form("p1s%s",c),i1,i2,"S");
  TH1D *h1s = p1s->ProjectionX(Form("h1%s",c));
  TProfile *p1x = p2x->ProfileY(Form("p1x%s",c),i1,i2,"");
  TProfile *p1xs = p2x->ProfileY(Form("p1xs%s",c),i1,i2,"S");
  TH1D *h1xs = p1xs->ProjectionX(Form("h1xs%s",c));

  // Extract JER (RMS)
  TH1D *h1jer = p1xs->ProjectionX(Form("h1jer_%s",c));
  for (int i = 1; i != h1s->GetNbinsX()+1; ++i) {
    double k = (scaleJER ? 1./2. : 1.);
    double rmsy = p1s->GetBinError(i) * k / sqrt(2.);
    double erry = p1->GetBinError(i) * k / sqrt(2.);
    h1s->SetBinContent(i, rmsy);
    h1s->SetBinError(i, rmsy ? erry / rmsy : 0);
    double rmsx = p1xs->GetBinError(i) * k / sqrt(2.);
    double errx = p1x->GetBinError(i) * k / sqrt(2.);
    h1xs->SetBinContent(i, rmsx);
    h1xs->SetBinError(i, rmsx ? errx / rmsx : 0);

    double jer = sqrt(max(rmsy*rmsy - rmsx*rmsx,0.));
    double err = (jer ? sqrt(pow(rmsy*erry/jer,2) + pow(rmsx*errx/jer,2)
			     + pow(0.001,2)) : 0);
    h1jer->SetBinContent(i, jer);
    h1jer->SetBinError(i, err);
  }

  //cout << h1 << ", " << (*h1) << endl;

  assert(h1);
  if (h1) *h1 = h1s;
  else delete h1s;
  assert(h1x);
  if (h1x) *h1x = h1xs;
  else delete h1xs;

  delete p1;
  delete p1s;
  delete p1x;
  delete p1xs;
  
  return h1jer;
} // getJER

TH1D *getJERZ(TH2D *h2, TH2D *h2x, const char *c,
	      TH1D **h1, TH1D **h1x) {

  TProfile *p = h2->ProfileX(Form("pz_%s",c),1,-1,"");
  TProfile *ps = h2->ProfileX(Form("pzs_%s",c),1,-1,"S");
  TH1D *h1s = ps->ProjectionX(Form("h1zs_%s",c));
  TProfile *px = h2x->ProfileX(Form("pzx_%s",c),1,-1,"");
  TProfile *pxs = h2x->ProfileX(Form("pzxs_%s",c),1,-1,"S");
  TH1D *h1xs = pxs->ProjectionX(Form("h1zxs_%s",c));
  
  // Extract JER (RMS)
  double k = sqrt(2.); // 1./sqrt(2) to compare MPFX to dijet, 1 for JER
  TH1D *h1jer = p->ProjectionX(Form("h1jer_%s",c));
  for (int i = 1; i != p->GetNbinsX()+1; ++i) {
    double rmsy = ps->GetBinError(i) /k;
    double erry = p->GetBinError(i) /k;
    h1s->SetBinContent(i, rmsy);
    h1s->SetBinError(i, rmsy ? erry / rmsy : 0);
    double rmsx = pxs->GetBinError(i) /k;
    double errx = px->GetBinError(i) /k;
    h1xs->SetBinContent(i, rmsx);
    h1xs->SetBinError(i, rmsx ? errx / rmsx : 0);

    double jer = sqrt(max(rmsy*rmsy - rmsx*rmsx,0.));
    double err = (jer ? sqrt(pow(rmsy*erry/jer,2) + pow(rmsx*errx/jer,2)
			     + pow(0.001,2)) : 0);
    h1jer->SetBinContent(i, jer *k);
    h1jer->SetBinError(i, err *k);
  } // for i

  assert(h1);
  if (h1) *h1 = h1s;
  else delete h1s;
  assert(h1x);
  if (h1x) *h1x = h1xs;
  else delete h1xs;

  delete p;
  delete ps;
  delete px;
  delete pxs;
  
  return h1jer;
} // getJERZ
