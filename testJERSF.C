// Purpose: test JER SF text file produced by minitools/jerSF
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TMath.h"

#include "tdrstyle_mod22.C"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

const bool debug = true;

void testJERSF(string filename, string era="X") {

  setTDRStyle();

  //filename = "../JECDatabase/textFiles/Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK4PFchs.txt"; // L2Res example, runs
  //filename = "../JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt"; // JER example, doesn't run
  cout << "Open " << filename << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(filename.c_str());
  vector<JetCorrectorParameters> v;
  v.push_back(*p);
  FactorizedJetCorrector *jer = new FactorizedJetCorrector(v);

  if (debug) {
    jer->setJetEta(0.);
    jer->setJetPt(80.);
    jer->setRho(20.85); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=80,rho=20.85)="<<jer->getCorrection()<<endl;
    jer->setJetEta(0.);
    jer->setJetPt(80.);
    jer->setRho(5.); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=80,rho=5)="<<jer->getCorrection()<<endl;
  }

  TString tera(era.c_str());
  string fileref;
  if (tera.Contains("Run2") || tera.Contains("UL2018"))
    fileref="../JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt";
  if (tera.Contains("UL2017"))
    fileref="../JRDatabase/textFiles/Summer19UL17_JRV3_MC/Summer19UL17_JRV3_MC_SF_AK4PFchs.txt";
  if (tera.Contains("UL2016GH"))
    fileref="../JRDatabase/textFiles/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt";
  if (tera.Contains("UL2016APV") || tera.Contains("UL2016BCDEF"))
    fileref="../JRDatabase/textFiles/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt";

  ifstream fin(fileref.c_str());
  char c[512]; fin.getline(c, 512);
  cout << c << endl;
  string s;
  double x1, x2, sf, sfdw, sfup;
  int n;
  vector<double> vx; vx.push_back(0);
  vector<double> vy, vym,vye;
  while (fin >> x1 >> x2 >> n >> sf >> sfdw >> sfup) {
    if (x2>0) {
      vx.push_back(x2);
      vy.push_back(sf);
      vym.push_back(0.5*(sfup+sfdw));
      vye.push_back(0.5*(sfup-sfdw));
    }
  }
  TH1D *hrefjerband = new TH1D("hrefjerband",";|#eta|",vx.size()-1,&vx[0]);
  for (unsigned int i = 0; i != vym.size(); ++i) {
    hrefjerband->SetBinContent(i+1, vym[i]);
    hrefjerband->SetBinError(i+1, vye[i]);
  }

  TH1D *h = tdrHist("h","JER SF",0.95,1.75,"#eta_{jet}",0,5.2);
  //lumi_13TeV = "UL2018, 59.9 fb^{-1}";
  lumi_13TeV = "Run2, 137.9 fb^{-1}";
  if (tera.Contains("UL2016APV")) lumi_13TeV = "2016early, 19.7 fb^{-1}";
  if (tera.Contains("UL2016BCDEF")) lumi_13TeV = "2016BCDEF, 19.7 fb^{-1}";
  if (tera.Contains("UL2016GH")) lumi_13TeV = "2016late, 16.8 fb^{-1}";
  if (tera.Contains("UL2017")) lumi_13TeV = "2017, 41.5 fb^{-1}";
  if (tera.Contains("UL2018")) lumi_13TeV = "2018, 59.9 fb^{-1}";
  if (tera.Contains("Run2")) lumi_13TeV = "Run2, 137.9 fb^{-1}";

  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  //tdrDraw(h8,"HPz",kOpenDiamond,kRed,kSolid,-1,kNone);
  //tdrDraw(h80,"HPz",kFullCircle,kBlack,kSolid,-1,kNone);
  //tdrDraw(hxmax,"HPz",kOpenDiamond,kBlue,kSolid,-1,kNone);
  
  tdrDraw(hrefjerband,"E0",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);

  const int ndiv = 1;//5;
  const int nbins = 60*ndiv;
  const double deta = TMath::TwoPi()/72./ndiv;
  const double maxeta = nbins*deta; // 5.236
  const double rho = 20.85;
  const double vpt[] = {8, 15, 30, 60, 120, 240, 480, 1000, 2000, 4000};
  const int npt = sizeof(vpt)/sizeof(vpt[0]);
  const int color[npt] = {kGray, kBlack, kMagenta+1, kBlue, kCyan+2, kGreen+2,
			  kYellow+1, kOrange+1, kRed, kBlack};

  TLegend *leg = tdrLeg(0.65,0.89-(npt+1)*0.035,0.90,0.89);
  leg->SetTextSize(0.04);
  leg->AddEntry(hrefjerband,Form("%s ref.",era.c_str()),"PLE");

  for (int ipt = 0; ipt != npt; ++ipt) {

    double pt = vpt[ipt];

    TH1D *hjer = new TH1D(Form("hjer%d",ipt),";|#eta_{jet}|;JER SF",nbins,0.,maxeta);

    for (int i = 0; i != hjer->GetNbinsX()+1; ++i) {

      double eta = hjer->GetBinCenter(i);
      double etamin = hjer->GetBinLowEdge(i);
      double ptmax = 6500./cosh(etamin);

      jer->setJetEta(eta);
      jer->setJetPt(pt);
      jer->setRho(rho);
      double cfp = jer->getCorrection();
      jer->setJetEta(-eta);
      jer->setJetPt(pt);
      jer->setRho(rho);
      double cfm = jer->getCorrection();
      if (pt < ptmax)
	hjer->SetBinContent(i, 0.5*(cfp+cfm));
    } // for i

    tdrDraw(hjer,"HIST][",kNone /*marker[ipt]*/,color[ipt],kSolid,-1,kNone);

    leg->AddEntry(hjer,Form("%1.0f GeV",pt),"L");
  } // for ipt

  c1->SaveAs(Form("pdf/jersf/testJERSF_%s.pdf",era.c_str()));
} // testJERSF
