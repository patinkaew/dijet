#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "tdrstyle_mod15.C"
#include "TLine.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TLatex.h"
#include "TFitter.h"
#include "Fit/FitUtil.h"
#include "TStopWatch.h"

#include <map>

const int kRightHatch = 3654;
const int kLeftHatch = 3645;
const int kVertHatch = 3699;
const int kHorizHatch = 3600;
const int kCrossHatch = 3644;

using namespace std;

bool quiet = false; // silence printout

// return values from toy fits
struct ensemble {
  Int_t chi2_tot;
  Int_t chi2_mjj;
  Int_t ndf_tot;
  Int_t ndf_mjj;
  Double_t xsec;
  Double_t xsecerr;
  Double_t xsectrue;
} ;

// hard-coded settings
bool sigmatot = false;   // sigma = sigma_stat (oplus) sigma_fit
bool plotchi2 = true;    // Display fit chi2/NDF
bool plotchi2gbl = true; // Display global fit chi2/NDF
bool plotkfactor = true; // k-factor of data normalization
bool plotparams = true;  // plot 2+4 fit parameters
bool datanorm = true;    // apply data normalization (as opposed to lumi norm)
bool logx  = false;      // Draw plot with log(pT)
bool statband = false;   // Draw fit stat uncertainty band
bool jecbands = false;   // Draw fit JEC uncertainty band

bool combofit = true;    // Combined fit to mass and turn-on
bool fithisto = false;   // fit wide-binned histogram
bool fitgraph = true;    // fit LW-centered graph (better)
bool fitgevbin = false;  // log-likelihood fit of 1 GeV binned histogra (best)

bool drawsig = true;     // Disable signal histogram
bool drawmon = false;    // Draw scouting-monitoring data
bool drawmc = false;     // Disable QCD MC
bool drawfitunc = true;  // Draw fit uncertainty on (data-fit)/sigma pad
bool drawref = true;     // Draw refence fit on (data-fit)/sigma pad
bool dicanvas = true;    // Draw both cross section and significance
bool drawlogx = true;    // Use logarithmic x-axis
bool savecov = false;    // Save covariance matrix plot

bool scouting = true;    // Writes scouting text to plot
bool blind = true;       // Blind 693 < mjj < 838 for scouting analysis
bool trigfit = true; // If false, use trigger ineff as prescale and normal fit
bool addtrigerr = false; // Add trigger uncertainty (if also trigfit=true)

// Just as reminder, pT binning of mjj
// 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328,
const int ntrg = 10;
Double_t xtrg[] = {453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890};

double blindlow = 649;
double blindhigh = 838;

// ****************************************************
// *** The blinded region is mjj > 693 && mjj < 838 ***
// ****************************************************

// FIXME Check more accurate HLT-lumi for scouting
double lumi = 1800;      // Integrated luminosity for normalization
double lumimon = 1.8;    // Scouting monitoring luminosity
double sigxs = 15e0;     // Signal cross section
double sigacc = 0.939;   // Mjj acceptance for Mjj<1118 GeV
double xmin = 419;       // Starting bin of plot
double xmax = 6328;      // Last bin of histogram
double fitmin = 565;//453;   // First bin for scouting fit
double fitmax = 2037;//6328; //Last bin for scouting fit
double trgfitmin = 453;  // First bin for trigger turn-on fit
double trgfitmax = 890;  // Last bin for trigger turn-on fit
double sigfitmax = 1118; // Last bin of signal shape
double miny = 0.2e-6;    // Plot minimum y-axis range
double maxy = 10000e0;   // Plot maximum y-axis range
double resymax = 3.5;    // Suboad maximum y-axis range
double qcd_min = 1200;   // Count MC and data xsec from here

bool fixSepTrg(false); bool fixSepSim(false);
const char *ctag = "_sigout"; bool sigin = false;
bool fitSig(false);

TF1 *_fitError_func(0);
TMatrixD *_fitError_emat(0);
Double_t fitError(Double_t *xx, Double_t *p);

// Combined fit to dijet mass and turn-on implemented as TFitter
TObject *_fitObj, *_fitObj2;
TF1 *_fitFunc, *_fitFunc2;
void comboFitter(Int_t &npar, Double_t *grad, Double_t &chi2, Double_t *par,
		Int_t flag);

TH1D *_hS(0);
const int nFitBtS = 7;
Double_t fitBtS(Double_t *xx, Double_t *p) {

  if (_hS==0) {
    TDirectory *curdir = gDirectory;
    TFile *fs = new TFile("h_gg_750.root","READ");
    TH1D *tmp = (TH1D*)fs->Get("h_gg_750cut"); assert(tmp); // Mjj*1.5 cut
    curdir->cd();

    _hS = (TH1D*)tmp->Clone("hFitBtS");
    _hS->Scale(1./tmp->Integral()); // is probility per bin
    fs->Close();
    assert(_hS);
  }

  double x = *xx;
  double B = p[0] * pow(1-x/13000.,p[1]) / 
    pow(x/13000., p[2] + p[3]*log(x/13000.));
  double t = (1./2. * (1 + TMath::Erf((x-p[4])/p[5])));
  assert(_hS);
  int k = _hS->FindBin(x);
  double S = sigacc*p[6]*_hS->GetBinContent(k)/_hS->GetBinWidth(k); // pb / GeV

  return ((B+S)*t);
}


ensemble dijetmass_scouting_sig_pars(const char *pname,
				     const double pfitmin, const double pfitmax,
				     const bool pSepFit, const bool pfitSig,
				     const string toy="");

const bool kSepFit(true);
const bool kSimFit(false);
void dijetmass_scouting_sig() {

  // zzz (short cut to get here)

  // Blinding steps
  bool sigfit = false;
  dijetmass_scouting_sig_pars("_blind", 565, 2037, kSimFit, sigfit);

  // No signal fit, just background
  /*
  blind = true;
  bool sigfit = false;
  dijetmass_scouting_sig_pars("_seps", 565, 2037, kSepFit, sigfit);
  dijetmass_scouting_sig_pars("_sims", 565, 2037, kSimFit, sigfit);
  dijetmass_scouting_sig_pars("_simn", 453, 2037, kSimFit, sigfit);
  dijetmass_scouting_sig_pars("_simm", 565, 6328, kSimFit, sigfit);
  dijetmass_scouting_sig_pars("_siml", 453, 6328, kSimFit, sigfit);
  */

  // Studies on a single toyMC sample
  /*
  string toy = "sbtoy1";
  blind = false;
  dijetmass_scouting_sig_pars("_seps_toy", 565, 2037, kSepFit,true,toy);
  dijetmass_scouting_sig_pars("_sims_toy", 565, 2037, kSimFit,false,toy);
  dijetmass_scouting_sig_pars("_siml_toy", 453, 6328, kSimFit,true,toy);
  */

  // Ensemble studies on toys
  /*
  TDirectory *curdir = gDirectory;
  TFile *fout = new TFile("toyhist_sims_x5acc.root","RECREATE");

  TH1D *bbchi2 = new TH1D("bbchi2","",1000,0,100);
  TH1D *bschi2 = new TH1D("bschi2","",1000,0,100);
  TH1D *sschi2 = new TH1D("sschi2","",1000,0,100);
  TH1D *bspull = new TH1D("bspull","",1000,-10,10);
  TH1D *sspull = new TH1D("sspull","",1000,-10,10);
  //
  TH1D *bsxsec   = new TH1D("bsxsec","",2000,-100,100);
  TH1D *bsxsecerr = new TH1D("bsxsecerr","",1000,0,100);
  TH1D *bsxsecsig = new TH1D("bsxsecsig","",200,-10,10);
  TH1D *ssxsec   = new TH1D("ssxsec","",2000,-100,100);
  TH1D *ssxsecerr = new TH1D("ssxsecerr","",1000,0,100);
  TH1D *ssxsecsig = new TH1D("ssxsecsig","",200,-10,10);
  
  const int nxs = 5;
  TH1D* vsschi2[nxs];
  TH1D* vsspull[nxs];
  TH1D* vssxsec[nxs];
  TH1D* vssxsecerr[nxs];
  TH1D* vssxsecsig[nxs];
  for (int i = 0 ; i != nxs; ++i) {
    vsschi2[i] = new TH1D(Form("sschi2_%d",i),"",1000,0,100);
    vsspull[i] = new TH1D(Form("sspull_%d",i),"",1000,-10,10);
    vssxsec[i]   = new TH1D(Form("ssxsec_%d",i),"",2000,-100,100);
    vssxsecerr[i] = new TH1D(Form("ssxsecerr_%d",i),"",1000,0,100);
    vssxsecsig[i] = new TH1D(Form("ssxsecsig_%d",i),"",200,-10,10);
  }
  curdir->cd();

  // Basic toy studies: B+B, B+SB, SB+SB
  TStopwatch stop; double time(0);
  stop.Start();
  ensemble tbb, tbs, tss;
  blind = false;
  const int ntoy = 1000;
  for (int i = 0; i != ntoy; ++i) {

    string s = Form("toy%d",i+1);
    // Regular short range fits
    tbb = dijetmass_scouting_sig_pars("_bb_xtoy",565,2037,kSimFit,false,"b"+s);
    tbs = dijetmass_scouting_sig_pars("_bs_xtoy",565,2037,kSimFit,true, "b"+s);
    tss = dijetmass_scouting_sig_pars("_ss_xtoy",565,2037,kSimFit,true,"sb"+s);

    // Complementary long range fits
    //tbb =dijetmass_scouting_sig_pars("_bb_xtoy",453,6328,kSimFit,false,"b"+s);
    //tbs =dijetmass_scouting_sig_pars("_bs_xtoy",453,6328,kSimFit,true, "b"+s);
    //tss =dijetmass_scouting_sig_pars("_ss_xtoy",453,6328,kSimFit,true,"sb"+s);

    // Separate fits, to check uncertainty
    //tbb =dijetmass_scouting_sig_pars("_bb_xtoy",565,2037,kSepFit,false,"b"+s);
    //tbs =dijetmass_scouting_sig_pars("_bs_xtoy",565,2037,kSepFit,true, "b"+s);
    //tss =dijetmass_scouting_sig_pars("_ss_xtoy",565,2037,kSepFit,true,"sb"+s);

    // Fill information on chi2 and pulls
    bbchi2->Fill(tbb.chi2_tot);
    bschi2->Fill(tbs.chi2_tot);
    sschi2->Fill(tss.chi2_tot);
    bspull->Fill( (tbs.xsec - tbs.xsectrue) / tbs.xsecerr);
    sspull->Fill( (tss.xsec - tss.xsectrue) / tss.xsecerr);

    // Extracted cross section and it's uncertainty
    bsxsec->Fill( tbs.xsec );
    bsxsecerr->Fill( tbs.xsecerr );
    bsxsecsig->Fill( tbs.xsec / tbs.xsecerr);
    ssxsec->Fill( tss.xsec );
    ssxsecerr->Fill( tss.xsecerr );
    ssxsecsig->Fill( tss.xsec / tss.xsecerr);

    cout << "Processed " << s << " with "
	 << "xsec = "<<tss.xsectrue<<" ("<<tbs.xsectrue<<")"<<endl << flush;
    quiet = true; // only make plots and numbers for first toy 

    // Analyse pulls and significances for difference cross sections
    for (int j = 0; j != nxs; ++j) {

      tss = dijetmass_scouting_sig_pars("_ss_xtoy",565,2037,kSimFit,true,
					Form("sb%s_%d",s.c_str(),j));

      vsschi2[j]->Fill( tss.chi2_tot );
      vsspull[j]->Fill( (tss.xsec - tss.xsectrue) / tss.xsecerr);
      vssxsec[j]->Fill( tss.xsec );
      vssxsecerr[j]->Fill( tss.xsecerr );
      vssxsecsig[j]->Fill( tss.xsec / tss.xsecerr);
    } // for j

    cout << "Generation took " << stop.RealTime()-time << " s" << endl;
    time = stop.RealTime();
    stop.Continue();
  } // for i
  fout->Write();
  fout->Close();
  */
  // false/true
}

TFile *f(0), *fq(0), *fs(0), *ft(0);
ensemble dijetmass_scouting_sig_pars(const char *pname,
				     const double pfitmin, const double pfitmax,
				     const bool pSepFit, const bool pfitSig,
				     const string toy) {

  // Set tunable parameters
  ctag = pname;
  fitmin = pfitmin;
  fitmax = pfitmax;
  fixSepTrg = pSepFit;
  sigin = false; // no signal injections
  fitSig = pfitSig;
  combofit = (!fixSepTrg && !fixSepSim);

  bool isToy = (toy!="");
  const char *ctoy = toy.c_str();
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();
  gStyle->SetOptStat(0);

  // Input data sample
  TFile *f = new TFile(isToy ? "toymc.root" :
		       "rawhists/rawhistV7_Run2015D_scoutingPFHT_BLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root","READ");
  assert(f && !f->IsZombie());

  // QCD MC for comparisons
  TFile *fq = new TFile("mjj_qcd_300toinf_reweight_v2.root", "READ");
  assert(fq && !fq->IsZombie());

  // Signal template
  TFile *fs = new TFile("h_gg_750.root","READ");
  assert(fs && !fs->IsZombie());

  // Trigger turn-on data
  TFile *ft = new TFile(isToy ? "toymc.root" :
			"triggerEfficiency_L1HTT150seed_HT450_DetaJJLess1p3_HLTv7Corr_output.root","READ");
  assert(ft && !ft->IsZombie());
  
  // Read in data histograms
  TH1D *hmjj, *hmjj1;
  hmjj  = (TH1D*)f->Get(isToy ? Form("mjj_%s",ctoy) : "mjj_mjjcor");
  hmjj1 = (TH1D*)f->Get(isToy ? Form("mjj_%s",ctoy) : "mjj_mjjcor_gev");
  assert(hmjj);
  assert(hmjj1);

  // Load QCD background MC
  TH1D *hqcd = (TH1D*)fq->Get("mjj");
  assert(hqcd);

  // Load signal template
  TH1D *hsig0 = (TH1D*)fs->Get("h_gg_750cut"); // Mjj*1.5 cut signal
  assert(hsig0);
  TH1D *hsig = (TH1D*)hsig0->Clone("hsig");
  hsig->GetXaxis()->SetRangeUser(trgfitmin,sigfitmax);
  // http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/EXO-14-005/CMS-PAS-EXO-14-005_Figure_003.png
  // => 8 TeV limit about 2pb, so 10 pb at 13 TeV (gg scaling should be x5)
  double sigmaGG750 = sigxs;
  double integral = hsig0->Integral();
  if (!quiet) cout << "Normalize gg750 by " << sigmaGG750 / integral << endl;
  for (int i = 1; i != hsig->GetNbinsX()+1; ++i) {

    double n = hsig->GetBinContent(i);
    double m = hsig->GetBinCenter(i);
    double dm = hsig->GetBinWidth(i);
    
    double y = n / dm * sigmaGG750*sigacc/integral;

    if (m>xmin && m<xmax) {

      hsig->SetBinContent(i, y);
    }
  } // for i

  TH1D *hS;
  hS = (TH1D*)hsig->Clone("hS");
  hS->Scale(1./sigmaGG750);

  // Print re-binned signal shape for generateMC
  vector<double> xs(0);
  vector<double> ys(0);
  for (int i = 1; i != hS->GetNbinsX()+1; ++i) {
    double x = hS->GetBinCenter(i);
    if (x>=419 && x<6328) {
      if (xs.size()==0) xs.push_back(hS->GetBinLowEdge(i));
      xs.push_back(hS->GetBinLowEdge(i+1));
      ys.push_back(hS->GetBinContent(i));
    }
  }
  if (!quiet) {
    cout << "_hSx["<<xs.size()<<"] = {";
    for (int i = 0; i != xs.size(); ++i) cout << xs[i] << ", "; cout<<"}"<<endl;
    cout << "_hSy["<<ys.size()<<"] = {";
    for (int i = 0; i != ys.size(); ++i) cout << ys[i] << ", "; cout<<"}"<<endl;
  }    

  // Background fit without signal injection for blind region placeholder
  TF1 *fsigout = new TF1("fsigout",fitBtS,fitmin,fitmax,nFitBtS);
  //fsigout->SetParameters(3.68542e-07, 5.1268, 7.38103, 0.395473, 495.253, 93.9084, 0); // values used for toyMC generation
  fsigout->SetParameters(3.74058e-07, 5.20434, 7.38254, 0.397342, 495.772, 92.4408, 0); // blind fit on April 19

  // Load trigger turn-on data
  TEfficiency *efftrg = (TEfficiency*)ft->Get(isToy ?
					      Form("efficiency_%s",ctoy) :
					      "efficiency");
  assert(efftrg);

  // Rebin trigger turn-on and scouting monitoring data to match analysis bins
  TH1D *hpass = (TH1D*)efftrg->GetCopyPassedHisto();
  TH1D *htot = (TH1D*)efftrg->GetCopyTotalHisto();
  TH1D *hpass_rebin = (TH1D*)hpass->Rebin(ntrg,"hpass_rebin",xtrg);
  TH1D *hmon = (TH1D*)htot->Rebin(ntrg,"hmon",xtrg);
  TH1D *htrgeff = (TH1D*)hpass_rebin->Clone("htrgeff");
  htrgeff->Sumw2(); // need this, even if gives warning?
  htrgeff->Divide(hpass_rebin,hmon,1,1,"B");
  
  // Scale signal for trigger efficiency (used in plotting only)
  TH1D *hsigmon = (TH1D*)hsig->Clone("hsigmon");
  for (int i = 0; i != htrgeff->GetNbinsX()+1 && trigfit; ++i) {
    int j = hsig->FindBin(htrgeff->GetBinCenter(i));
    hsig->SetBinContent(j, hsig->GetBinContent(j)
			* htrgeff->GetBinContent(i));
  }
  
  curdir->cd();


  // Re-fit turn-on for adding trigger uncertainty to data for separate fits
  TF1 *ftrg0 = new TF1("ftrg0","0.5 * (1 + TMath::Erf((x-[0])/[1]))",
		       trgfitmin, trgfitmax);
  ftrg0->SetParameters(500,100);
  efftrg->Fit(ftrg0, quiet ? "QRN" : "RN");
  if (!quiet) cout << "trg0: " << ftrg0->GetParameter(0) << ", "
		   << ftrg0->GetParameter(1) << endl;

  TMatrixD emat0(ftrg0->GetNpar(),ftrg0->GetNpar());
  gMinuit->mnemat(emat0.GetMatrixArray(), ftrg0->GetNpar());

  // Statistical uncertaity band used in mass normalization loop below
  TF1 *fet = new TF1("trgFitError", fitError, trgfitmin, trgfitmax, 1);
  fet->SetParameter(0,+1);
  _fitError_func = ftrg0;
  _fitError_emat = &emat0;
  TH1D *hmjjet = (TH1D*)hmjj->Clone("hmjjet");
  for (int i = 1; i != hmjj->GetNbinsX()+1; ++i) {
    double x = hmjj->GetBinCenter(i);
    hmjjet->SetBinContent(i, addtrigerr ? fabs(fet->Eval(x)-ftrg0->Eval(x)) : 0);
  }

  // Normalize by bin width and luminosity,
  // plus calculate asymmetric statistical uncertainty band
  // See: https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/scripts/doTestBackgoundFit.py#L109-L152

  TGraphAsymmErrors *gmjj = new TGraphAsymmErrors(0);
  TGraphAsymmErrors *gmjj_fit = new TGraphAsymmErrors(0);
  map<int, int> htog;
  const double alpha = 1 - 0.6827; // where from?

  double dataxsec(0);
  for (int i = 1; i != hmjj->GetNbinsX()+1; ++i) {
    double m = hmjj->GetBinCenter(i);
    if (m<qcd_min) continue;
    if (m>xmin && m<xmax) {
      dataxsec += hmjj->GetBinContent(i) / lumi;
    }
  }
  double mcxsec(0);
  for (int i = 1; i != hqcd->GetNbinsX()+1; ++i) {
    double m = hqcd->GetBinCenter(i);
    if (m<qcd_min) continue;
    
    if (m>xmin && m<xmax) {
      mcxsec += hqcd->GetBinContent(i);
    }
  }
  double kfactor = (dataxsec / mcxsec);
  if (!quiet) {
    cout << "The k-factor is " << kfactor << endl;
    cout << "We apply " << (datanorm ? "data normalization" :
			    "luminosity normalization") << endl;
  }    

  for (int i = 1; i != hmjj->GetNbinsX()+1; ++i) {

    double n = hmjj->GetBinContent(i);
    double m = hmjj->GetBinCenter(i);
    double dm = hmjj->GetBinWidth(i);

    if (sigin) { // signal injection test
      int j = hsig->FindBin(m); //assert(i==j);
      n += hsig->GetBinContent(j) * dm * lumi;
    }

    // zero the errors of blinded bins to help fitting
    if (blind && m > blindlow && m < blindhigh){
      hmjj->SetBinContent(i,0);
      hmjj->SetBinError(i,0);
      if (!quiet) cout << "Zeroed bin error with bin center " << m << endl;
    }

    double y = n / (dm * lumi);
    double l = 0.5*TMath::ChisquareQuantile(alpha/2, 2*n);
    double h = 0.5*TMath::ChisquareQuantile(1-alpha/2, 2*(n+1));
    double eyl = (n-l) / (dm * lumi);
    double eyh = (h-n) / (dm * lumi);

    if (!trigfit) {
      // Correct trigger efficiency bin-by-bin and add uncertainty
      int j = htrgeff->FindBin(m);
      double eff = htrgeff->GetBinContent(j);
      double erreff = htrgeff->GetBinError(j);

      y *= 1./eff;
      eyl = sqrt(eyl*eyl/(eff*eff) + y*y*erreff*erreff);
      eyh = sqrt(eyh*eyh/(eff*eff) + y*y*erreff*erreff);
    }
    if (trigfit && !combofit && addtrigerr) {
      // Add trigger fit uncertainty to separate fits
      // assume bin uncertainties uncorrelated
      // (not exactly true, but hopefully close enough)
      double erreff = hmjjet->GetBinContent(i);
      if (!quiet) cout << "m="<<m<<" erreff="<<erreff<<endl;
      eyl = sqrt(eyl*eyl + y*y*erreff*erreff);
      eyh = sqrt(eyh*eyh + y*y*erreff*erreff);
    }

    if (y==0) eyl = eyh; // patch for ey=0.5*(eyl+eyh) when y==0

    hmjj->SetBinContent(i, y);
    hmjj->SetBinError(i, 0.5*(eyl+eyh));

    // Graph with asymmetric Poisson errors
    if (m>xmin && m<xmax) {

      dataxsec += n / lumi;

      int j = gmjj->GetN();
      double x = m;
      double ex = 0.5*dm;

      // Only use points outside the blinded region
      if (!(blind && m>=blindlow && m<=blindhigh)) {
	
	gmjj->SetPoint(j, x, y);
	gmjj->SetPointError(j, ex, ex, eyl, eyh);
	
	gmjj_fit->SetPoint(j, x, y);
	gmjj_fit->SetPointError(j, 0, 0, eyl, eyh);
      }

      htog[i] = j;
    }
  }

  // Repeat above to add points from scouting monitoring trigger
  TGraphAsymmErrors *gmon = new TGraphAsymmErrors(0);
  TGraphAsymmErrors *gmon_fit = new TGraphAsymmErrors(0);
  map<int, int> htog_trg;
  for (int i = 1; i != hmon->GetNbinsX()+1; ++i) {

    double n = hmon->GetBinContent(i);
    double m = hmon->GetBinCenter(i);
    double dm = hmon->GetBinWidth(i);
    
    if (sigin) { // signal injection test
      int j = hsig->FindBin(m); //assert(i==j);
      n += hsig->GetBinContent(j) * dm * lumimon;
    }

    double y = n / (dm * lumimon);
    double l = 0.5*TMath::ChisquareQuantile(alpha/2, 2*n);
    double h = 0.5*TMath::ChisquareQuantile(1-alpha/2, 2*(n+1));
    double eyl = (n-l) / (dm * lumimon);
    double eyh = (h-n) / (dm * lumimon);

    if (y==0) eyl = eyh; // patch for ey=0.5*(eyl+eyh) when y==0

    hmon->SetBinContent(i, y);
    hmon->SetBinError(i, 0.5*(eyl+eyh));

    // Graph with asymmetric Poisson errors
    if (m>trgfitmin && m<trgfitmax) {

      int j = gmon->GetN();
      double x = m;
      double ex = 0.5*dm;

      // Use points everywhere, also in blinded region
      if (true) {
	gmon->SetPoint(j, x, y);
	gmon->SetPointError(j, ex, ex, eyl, eyh);
	
	gmon_fit->SetPoint(j, x, y);
	gmon_fit->SetPointError(j, 0, 0, eyl, eyh);
      }

      htog_trg[i] = j;
    }
  }


  TGraph *gqcd = new TGraph(0);
  for (int i = 1; i != hqcd->GetNbinsX()+1; ++i) {

    double n = hqcd->GetBinContent(i);
    double m = hqcd->GetBinCenter(i);
    double dm = hqcd->GetBinWidth(i);
    
    double y = n / dm;
    if (datanorm) y *= kfactor;

    if (m>xmin && m<xmax) {

      int j = gqcd->GetN();
      gqcd->SetPoint(j, m, y);
    }
  } // for i


  TF1 *f1;
  if (!trigfit){
    assert(false); // this block has not been maintained
    f1 = new TF1("f1","[0]*pow(1-x/13000.,[1])/pow(x/13000.,[2]"
		 "+[3]*log(x/13000.))", fitmin, fitmax); // 4-par
    f1->SetParameters(2.853e-07, 4.604, 7.358, 0.364);
  }
  else{

     f1 = new TF1("f1",fitBtS,fitmin,fitmax,nFitBtS);
     f1->SetParameters(3.68542e-07, 5.1268, 7.38103, 0.395473, 495.253, 93.9084, 0); // values used for toyMC generation

     // Set signal xsec to zero if not fitting it
     if (!fitSig) f1->FixParameter(nFitBtS-1, 0);

     const int np = f1->GetNpar();
     if (fixSepTrg) {
       // Trigger parameters from separate turnon fit (TEfficiency)
       f1->FixParameter(np-3, isToy ? 495.253 : 495.7); // MjjCorr
       f1->FixParameter(np-2, isToy ? 93.9084 : 96.0);  // MjjCorr
     }
     if (fixSepSim) {
       // Trigger parameters from simultaneous fit
       f1->FixParameter(np-3, 494.8); // HLTV7
       f1->FixParameter(np-2, 93.0); // HLTV7
     }

   }
  const int np = f1->GetNpar();

  _fitObj = 0; // set later
  efftrg->Draw(); gPad->Update();
  _fitObj2 = efftrg->GetPaintedGraph();
  _fitFunc = f1;
  _fitFunc2 = new TF1("f2","(1./2.*(1+TMath::Erf((x-[0])/[1])))",
		      trgfitmin,trgfitmax);
  _fitFunc2->SetParameters(f1->GetParameter(np-3),f1->GetParameter(np-2));

  // Fitting schemes in order of preference (best as last)
  // C: Simple chi2 fit to wide-binned graph
  // (updated custom fitter uses integrals, so pretty good now as well)
  if (fithisto) {

    if (!quiet) cout << "Histogram fit (fithisto: hmjj)" << endl;

    TVirtualFitter::Fitter(hmjj)->SetFCN(comboFitter);
    _fitObj = hmjj;
    hmjj->Fit(f1, quiet ? "QURN" : "URN");
  }
  // B: Chi2 fit to graph with proper bin centering (iterate until convergence)
  if (fitgraph){

    if (!quiet) cout << "Graph fit (fitgraph: gmjj_fit)" << endl;

    TVirtualFitter::Fitter(gmjj_fit)->SetFCN(comboFitter);
    _fitObj = gmjj_fit;
    gmjj_fit->Fit(f1, quiet ? "QURN" : "URN"); // initial fit
    
    // Iterate graph fit with properly centered points until convergence
    // or maximum number of iterations reached
    const int nmax_iter = 10;
    int itcount(0);
    bool converged(false);
    while (!converged && itcount<nmax_iter) {
      
      // Apply proper bin centering to graph after doing initial fit
      for (int i = 0; i != gmjj->GetN(); ++i) {
	
	double x = gmjj->GetX()[i];
	double y = gmjj->GetY()[i];
	double xl = x - gmjj->GetErrorXlow(i);
	double xh = x + gmjj->GetErrorXhigh(i);
	double eyl = gmjj->GetErrorYlow(i);
	double eyh = gmjj->GetErrorYhigh(i);
	
	// A 2-liner to implement the Lafferty-Wyatt x_LW. For details, see:
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/StatComWideBins?
	double yavg = f1->Integral(xl, xh) / (xh-xl);
	double x_LW = f1->GetX(yavg, xl, xh);
	
	double exl = x_LW - xl;
	double exh = xh - x_LW;
	
	gmjj->SetPoint(i, x_LW, y);
	gmjj_fit->SetPoint(i, x_LW, y);
	gmjj->SetPointError(i, exl, exh, eyl, eyh);
	gmjj_fit->SetPointError(i, 0, 0, eyl, eyh);
      } // for i

      // Iterate additional fits to properly centered points
      // (question: are exl, exh considered in the fit by ROOT?)
      if (!quiet) cout << Form("Graph fit (fitgraph iteration #%d: gmjj_fit)",
			       itcount+2) << endl;
      
      _fitObj = gmjj_fit;
      TVirtualFitter::Fitter(gmjj_fit)->SetFCN(comboFitter);
      
      // Details on how to fit and retrieve results:
      // https://root.cern.ch/doc/master/classTGraph.html#aa978c8ee0162e661eae795f6f3a35589
      TFitResultPtr fr = gmjj_fit->Fit(f1, quiet ? "QURN" : "URN");
      Int_t fitstatus = fr;
      if (fitstatus==0) converged = true;
      ++itcount;

    } // while !converged

    // Save centered graph for later interactive analysis (optimizing fit shape)
    if (!quiet) {
      TDirectory *curdir = gDirectory;
      TFile *fout = new TFile("datagraph.root","RECREATE");
      gmjj_fit->Write("gmjj_fit");
      fout->Close();
      curdir->cd();
    }
  } // fitgraph
  // A: Log-likelihood fit to 1 GeV binned histogram
  if (fitgevbin) {
    
    assert(false); // currently does not support combined fit
    if (!quiet) cout << "1 GeV histogram fit (hmjj1)" << endl;

    // Problem: how can we do simultaneous log-likelihood fit,
    // or add turn-on parameters as nuisances in LL?
    hmjj1->Fit(f1, quiet ? "QLRN" : "LRN");
  }

  if (!quiet) {
    cout << "Final fit params: " << endl;
    for (int i=0; i!=f1->GetNpar();++i){
      cout << Form("%1.4g",f1->GetParameter(i));
      if (i!=f1->GetNpar()-1) cout << ", ";
      else cout << endl;
    } 
  }

  // Calculate chi2/NDF
  double chi2(0);
  int ndf = 0;
  for (int i = 0; i != gmjj->GetN(); ++i) {
    double y = gmjj->GetY()[i];
    double x = gmjj->GetX()[i];
    
    double ibin = hmjj->FindBin(x);
    double xl = hmjj->GetBinLowEdge(ibin);
    double xh = hmjj->GetBinLowEdge(ibin+1);
    double yfit = f1->Integral(xl,xh)/(xh-xl);
    
    double dy = yfit - y;
    double ey = (dy>=0 ? gmjj->GetErrorYhigh(i) : gmjj->GetErrorYlow(i));
    assert(ey>0); // can't assert in blind analysis!
    if (x>fitmin && x<fitmax && !(blind && x>blindlow && x<blindhigh)) {
      chi2 += dy*dy / (ey*ey);
      ++ndf;
    }
  }

  // Semiblind
  
  if (!quiet) cout << "f1 has " << f1->GetNpar() << " parameters, "
		   << f1->GetNumberFreeParameters() << " free parameters and "
		   << f1->GetNDF() << " degrees of freedom (np="
		   << ndf << ")" << endl;
  //ndf -= f1->GetNumberFreeParameters();
  // above does not seem to count degrees of freedom correctly
  ndf -= f1->GetNpar();
  if (!fitSig) ++ndf; // p[6] fixed
  if (!combofit) ndf += 2; // p[4] and p[5] fixed


  // Graphics objects
  TLine *l = (quiet ? 0 : new TLine());
  TLatex *tex = (quiet ? 0 : new TLatex());
  if (tex) tex->SetTextSize(0.045);
  if (tex) tex->SetNDC();
  TCanvas *c1(0);

  TH1D *he2r(0), *hscstat(0), *hsigoy(0);
  if (!quiet) { // Mjj spectrum

  const int Npar = f1->GetNpar();
  TMatrixD emat(Npar, Npar);
  gMinuit->mnemat(emat.GetMatrixArray(), Npar);

  // Statistical uncertaity band
  TF1 *fke = new TF1("fitError", fitError, xmin, xmax, 1);
  _fitError_func = f1;
  _fitError_emat = &emat;

  // Systematic uncertainty band
  TF1 *f1s = new TF1("f1s","[0]*pow(1-(x*[3])/13000.,[1])/"
		     "pow((x*[3])/13000.,[2])"
		     // normalize to 1 at 1181 GeV (change from 1118 21.10.2015)
		     "* pow((13000-1181.)/(13000.-1181*[3]), [1])"
		     "* pow([3],[2])", 
		     xmin, xmax);
  f1s->SetParameters(f1->GetParameter(0), f1->GetParameter(1),
		     f1->GetParameter(2), 1.05);

  // Statistical uncertainty band for the fit
  TH1D *herr = (TH1D*)hmjj->Clone("herr"); herr->Clear();
  fke->SetParameter(0,+1);
  for (int i = 1; i != herr->GetNbinsX()+1; ++i) {
    double x = herr->GetBinCenter(i);
    if (x>xmin && x < xmax) {
      herr->SetBinContent(i, f1->Eval(x));
      herr->SetBinError(i, fabs(fke->Eval(x)-f1->Eval(x)));
    }
  }

  // Fit uncertainty band relative to statistical uncertainty for data
  TH1D *he2 = (TH1D*)hmjj->Clone("he2"); he2->Clear();
  he2r = (TH1D*)hmjj->Clone("he2r"); he2r->Clear();
  fke->SetParameter(0,+1);
  for (int i = 1; i != he2->GetNbinsX()+1; ++i) {
    double x = he2->GetBinCenter(i);
    if (x>xmin && x < xmax) {
      double k = he2->GetBinWidth(i)*lumi;
      double nfit = f1->Eval(x) * k;
      he2->SetBinContent(i, 0.);
      he2->SetBinError(i, herr->GetBinError(i)*k / sqrt(nfit));
      he2r->SetBinContent(i, 0.);
      he2r->SetBinError(i, 100.*herr->GetBinError(i)*k / nfit);
    }
  }


  // Stat+JEC uncertainty band
  TH1D *herr2 = (TH1D*)hmjj->Clone("herr2"); herr2->Clear();
  for (int i = 1; i != herr2->GetNbinsX()+1; ++i) {
    double x = herr2->GetBinCenter(i);
    if (x>xmin && x < xmax) {
      herr2->SetBinContent(i, f1->Eval(x));
      //f1s->SetParameter(3, jecshift[i]); // bin-by-bin
      herr2->SetBinError(i, sqrt( pow(fabs(f1s->Eval(x)-f1->Eval(x)),2)+
				  pow(herr->GetBinError(i),2)));
    }
  }

  // (Data-Fit) / ((Total) uncertainty)
  TH1D *hratio = (TH1D*)hmjj->Clone("hratio");
  TH1D *hsignalratio = (TH1D*)hmjj->Clone("hsignalratio");
  hscstat = (TH1D*)hmjj->Clone("hscstat");
  hsigoy = (TH1D*)hmjj->Clone("hsigoy");
  TH1D *hfitratio = (TH1D*)hmjj->Clone("hfitratio");
  TH1D *hfitratio_inv = (TH1D*)hmjj->Clone("hfitratio_inv");
  for (int i = 1; i != hratio->GetNbinsX()+1; ++i) {

    double x = hmjj->GetBinCenter(i);
    double xl = hmjj->GetBinLowEdge(i);
    double xh = hmjj->GetBinLowEdge(i+1);

    // sigma as statistical uncertainty only
    if (hmjj->GetBinError(i)!=0 && x>xmin && x<xmax
	&& x>fitmin && x<fitmax 
	&& !(blind && x>blindlow && x<blindhigh) ) {

      int j = htog[i];
      double y = gmjj->GetY()[j];

      double y_fit = f1->Integral(xl,xh)/(xh-xl);
      double y_fit0 = fsigout->Integral(xl,xh)/(xh-xl);
      double eyl = gmjj->GetErrorYlow(j);
      double eyh = gmjj->GetErrorYhigh(j);
      double sigma = (y_fit>y ? eyh : eyl);
      double sigma0 = sqrt(y_fit/((xh-xl)*lumi));
      assert(sigma!=0);
    
      if (sigmatot){
	sigma = sqrt(pow(hmjj->GetBinError(i),2) + 
		     pow(herr2->GetBinError(i),2));
      }
      hratio->SetBinContent(i, (hmjj->GetBinContent(i)
				- y_fit) / sigma);
      if (sigin)
	hsignalratio->SetBinContent(i, (hsig->GetBinContent(hsig->FindBin(x))
					- (y_fit-y_fit0)) / sigma);
      if (!sigin) 
	hsignalratio->SetBinContent(i, hsig->GetBinContent(hsig->FindBin(x))
				    / sigma);
      
      hfitratio->SetBinContent(i, (y_fit0 - y_fit) / sigma);
      hfitratio_inv->SetBinContent(i, (y_fit - y_fit0) / sigma);
      hscstat->SetBinError(i, 100. * sigma / y_fit);
      hscstat->SetBinContent(i, 0.);
      hsigoy->SetBinContent(i, 100. * hsig->GetBinContent(i) / y_fit);
      hsigoy->SetBinError(i, 0.);
    }
    else { // blinded region or outside fit range

      hratio->SetBinContent(i, 0);
      if (blind && x>blindlow && x<blindhigh) {
	double y_fit0 = fsigout->Integral(xl,xh)/(xh-xl);
	double y_fit = f1->Integral(xl,xh)/(xh-xl);
	double sigma0 = sqrt(y_fit0/((xh-xl)*lumi));
	if (trigfit && !combofit) {
	  sigma0 = sqrt(sigma0*sigma0+pow(y_fit0*hmjjet->GetBinContent(i),2.));
	}
	if (sigin) 
	  hsignalratio->SetBinContent(i, (hsig->GetBinContent(hsig->FindBin(x))
					  - (y_fit - y_fit0)) / sigma0);
	if (!sigin)
	  hsignalratio->SetBinContent(i, hsig->GetBinContent(hsig->FindBin(x))
				      / sigma0);
	hratio->SetBinContent(i, 0);

	hfitratio->SetBinContent(i, (y_fit0 - y_fit) / sigma0);
	hfitratio_inv->SetBinContent(i, (y_fit - y_fit0) / sigma0);
	hscstat->SetBinError(i, 100. * sigma0 / y_fit0);
	hscstat->SetBinContent(i, 0.);
	hsigoy->SetBinContent(i, 100. * hsig->GetBinContent(i) / y_fit);
	hsigoy->SetBinError(i, 0.);
      }
      else {
	hsignalratio->SetBinContent(i, 0);
	hfitratio->SetBinContent(i, 0);
	hfitratio_inv->SetBinContent(i, 0);
	hscstat->SetBinError(i, 100.*hmjj->GetBinError(i)
			     / hmjj->GetBinContent(i));
	hscstat->SetBinContent(i, 0.);
	hsigoy->SetBinContent(i, 0.);
	hsigoy->SetBinError(i, 0.);
      }
    }
    hratio->SetBinError(i, 0);
    hsignalratio->SetBinError(i, 0);
    hfitratio->SetBinError(i, 0);
    hfitratio_inv->SetBinError(i, 0);
  }

   
  TH1D *hup = new TH1D("hup",";Dijet mass M_{jj} (GeV);"
		       "d#sigma / dM_{jj} (pb / GeV)",
		       100,xmin,xmax);
  hup->SetMinimum(miny);
  hup->SetMaximum(maxy); // SetLinx
  if (logx) hup->SetMaximum(2e1); // SetLogx


  TH1D *hdw = new TH1D("hdw",Form(";M_{jj} (GeV);"
				  "(Data-Fit)/#sigma%s",
				  sigmatot ? "_{tot}" : ""),
		       100,xmin,xmax);
  hdw->SetMinimum(-resymax);
  hdw->SetMaximum(resymax);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  // Fewer x-axis labels for less clutter
  hup->GetXaxis()->SetNdivisions(505);
  hdw->GetXaxis()->SetNdivisions(505);

  lumi_13TeV = Form("%1.1f fb^{-1}",lumi/1000.);
  c1 = (dicanvas ? tdrDiCanvas("cmjj", hup, hdw, 4, 11) :
	tdrCanvas("cmjj", hup, 4, 11, kSquare));

  // Move titles a bit to not cut it off at the top or bottom
  hup->GetYaxis()->SetTitleOffset(hup->GetYaxis()->GetTitleOffset()-0.06);
  
  if (dicanvas) c1->cd(1);
  gPad->SetLogy();
  if (logx) gPad->SetLogx();

  hmjj->GetXaxis()->SetRangeUser(xmin,xmax);

  if (jecbands) tdrDraw(herr2,"E3",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);
  if (statband) tdrDraw(herr,"E3",kNone,kYellow+1,kSolid,-1,1001,kYellow+1);
  f1->SetLineWidth(2);
  if (drawmc) tdrDraw(gqcd,"L",kNone,kBlue,kDashed,kBlue);
  gqcd->SetLineWidth(2);
  if (drawsig && drawmon) tdrDraw(hsigmon,"HIST ][",kDashed,kCyan+1,kDashed,kCyan+1,kNone);
  hsigmon->SetLineWidth(2);
  if (drawsig) tdrDraw(hsig,"HIST ][",kDashed,kBlue,kDashed,kBlue,kNone,kBlue);
  hsig->SetLineWidth(2);

  // Draw points from scouting monitoring trigger
  if (drawmon) tdrDraw(gmon,"Pz",kOpenCircle);
  TF1 *fmon = new TF1("fmon","([0]*pow(1-x/13000.,[1])/pow(x/13000.,[2]"
		      "+[3]*log(x/13000.) ))",
		      trgfitmin, trgfitmax); // 4-par
  fmon->SetParameters(f1->GetParameter(0), f1->GetParameter(1),
		      f1->GetParameter(2), f1->GetParameter(3));
  fmon->SetLineColor(kMagenta);
  fmon->SetLineWidth(2);
  if (drawmon) fmon->Draw("SAME");

  tdrDraw(gmjj,"Pz",kFullCircle);
  f1->SetRange(453,xmax);
  TF1 *f1c = (TF1*)f1->DrawClone("SAME");
  f1c->SetLineStyle(kDotted);
  f1->SetRange(fitmin,fitmax);
  f1->Draw("SAME"); // Draw after to get fit ontop points

  gPad->RedrawAxis();

  tex->SetTextSize(0.045);
  tex->SetNDC();
 
  if (dicanvas){  
     c1->cd(2);
     if (logx) gPad->SetLogx();
     
     if (drawfitunc) tdrDraw(he2,"E3",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);

     l->SetLineStyle(kDashed);
     l->DrawLine(xmin,+2,xmax,+2);
     l->DrawLine(xmin,-2,xmax,-2);
     l->SetLineStyle(kDotted);
     l->DrawLine(xmin,+1,xmax,+1);
     l->DrawLine(xmin,-1,xmax,-1);

     tex->SetTextSize(0.045);
     tex->SetTextColor(kRed);
     tex->DrawLatex(0.40,50,"Blind region interpolated");

     if (blind){

       // Down panel       
       l->SetLineStyle(kDashed);
       l->SetLineWidth(2);
       l->SetLineColor(kRed);
       l->DrawLine(blindlow,-3,blindlow,3);
       l->DrawLine(blindhigh,-3,blindhigh,3);
       
       c1->cd(1);
       // Top panel
       l->DrawLine(blindlow,0.75,blindlow,150);
       l->DrawLine(blindhigh,0.75,blindhigh,150);
       c1->cd(2);

     }

     tdrDraw(hratio,"H");

     // Draw change in fit after signal injection
     hfitratio->SetLineWidth(2);
     hfitratio->GetXaxis()->SetRangeUser(xmin, xmax);
     if (drawref) tdrDraw(hfitratio,"H][",kNone,kCyan+2,kSolid,kCyan+2,kNone);
     hfitratio_inv->SetLineWidth(2);
     hfitratio_inv->GetXaxis()->SetRangeUser(xmin, xmax);
     if (drawref) tdrDraw(hfitratio_inv,"H][",kNone,kCyan+2,kDashed,kCyan+2,kNone);

     // Draw signal or signal injection
     hsignalratio->SetLineWidth(2);
     hsignalratio->GetXaxis()->SetRangeUser(fitmin,sigfitmax);
     tdrDraw(hsignalratio,"H][",kNone,kBlue,kSolid,kBlue,kNone);

     double xt = hdw->GetYaxis()->GetTitleSize()
       / hup->GetYaxis()->GetTitleSize();
     if (drawfitunc) {
       TLegend *leg2a = tdrLeg(0.73, 0.66, 0.98, 0.90);
       leg2a->SetTextSize(leg2a->GetTextSize()*xt);
       leg2a->AddEntry(hratio,"Data","F");
       if (drawfitunc) leg2a->AddEntry(he2,"Fit unc.","F");
       leg2a->AddEntry(hsignalratio,"Signal","L");
     }
     if (drawref) {
       TLegend *leg2b = tdrLeg(0.73,0.34,0.98,0.50);
       leg2b->SetTextSize(leg2b->GetTextSize()*xt);
       leg2b->AddEntry(hfitratio,"Ref. fit","L");
       leg2b->AddEntry(hfitratio_inv,"mirrored","L");
     }

     gPad->RedrawAxis();

     gPad->SetLogx();
     c1->cd(1);
     gPad->SetLogx();
  }

  tex->SetTextColor(kBlack);
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.22,"|#eta| < 2.5, |#Delta#eta| < 1.3");
  tex->DrawLatex(0.20,0.16,Form("M_{jj} > %1.1f TeV",xmin/1000.));
  tex->DrawLatex(0.20,0.10,"Wide jets");

  // chi2 may not display properly for TGraphAsymmErrors fit
  bool bands = (statband || jecbands);
  if (plotchi2 || bands || sigmatot)
    tex->DrawLatex(0.20,0.28,Form("#chi^{2} / NDF = %1.2f / %d",
				  chi2, ndf));
  if (plotkfactor)
    tex->DrawLatex(0.47, 0.60, Form("Data / MC = %1.2f (>%1.1f TeV)", kfactor, qcd_min/1000));

  if (scouting){
    tex->SetTextColor(kRed);
    tex->SetTextSize(0.06);
    tex->DrawLatex(0.56, 0.53, isToy ? "Toy MC" : "Scouting data");
    tex->SetTextSize(0.045);
    if (isToy) tex->DrawLatex(0.74, 0.53, Form("(%s)",ctoy));
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.67, 0.48, "Fit range:");
    tex->DrawLatex(0.67, 0.43, Form("%1.0f-%1.0f GeV",fitmin,fitmax));
  }
  tex->SetTextColor(kBlue);
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20, 0.47,
		 sigin ? Form("Signal inject.: %1.0f pb",sigmaGG750) :
		 Form("Signal: %1.0f pb",sigmaGG750));

  if (fitSig) {
    tex->SetTextColor(kBlue);
    tex->DrawLatex(0.50,0.10,Form("#sigma_{gg} = %1.1f #pm %1.1f pb",
				  f1->GetParameter(6), f1->GetParError(6)));
  }
  if (plotparams) {
    tex->SetTextColor(kGray+1);
    tex->SetTextSize(0.020);
    tex->DrawLatex(0.20,0.07,f1->GetExpFormula());
    tex->SetTextSize(0.030);
    tex->DrawLatex(0.20,0.04,Form("{%1.5g, %1.5g, %1.5g, %1.6g; %1.4g, %1.4g"
				  "; %1.4g}",
				  f1->GetParameter(0), f1->GetParameter(1),
				  f1->GetParameter(2), f1->GetParameter(3),
				  f1->GetParameter(4), f1->GetParameter(5),
				  f1->GetParameter(6)));
  }
  if (!quiet) {
    cout << "\nFit parameters (second line blind interpolation):" << endl;
    for (int i = 0; i != f1->GetNpar(); ++i) {
      cout << f1->GetParameter(i) << ", ";
    }
    cout << endl;
    for (int i = 0; i != fsigout->GetNpar(); ++i) {
      cout << fsigout->GetParameter(i) << ", ";
    }
    cout << endl << endl;
  }

  if (drawmon) {
    TLegend *leg2 = tdrLeg(0.48, 0.78, 0.73, 0.90);
    leg2->AddEntry(gmon," ","PL");
    leg2->AddEntry(fmon," ","L");
  }

  TLegend *leg = tdrLeg(0.55, 0.66, 0.80, 0.90);
  leg->AddEntry(gmjj,sigin ? "Data+signal" : "Data","PL");
  leg->AddEntry(f1,sigin ? "Fit to data+signal" : "Fit to data","L");
  if (!drawmc) leg->AddEntry(f1c,"Extrap. fit","L");
  if (drawmc) leg->AddEntry(gqcd,"QCD MC","L");
  if (drawsig) leg->AddEntry(hsig,"gg (750 GeV)","L");//"LF");
  if (statband) leg->AddEntry(herr,"Fit uncertainty","F");
  if (jecbands)  leg->AddEntry(herr2,"Fit+JEC uncert.","F");


  TH1D *hpull(0);
  if (dicanvas){
     c1->cd(2);
     hpull = new TH1D("hpull",";(Data-Fit)/#sigma_{tot};",12,-3,3);
     double sigmalocal(0), prevsign(0);
     int nbins(0), zerox(0);
     for (int i = 1; i != hratio->GetNbinsX()+1; ++i) {

       double x = hratio->GetBinCenter(i);
       if (x>xmin && x<xmax) {

         ++nbins;
         hpull->Fill(hratio->GetBinContent(i));
         if (fabs(hratio->GetBinContent(i))>sigmalocal)
	   sigmalocal = fabs(hratio->GetBinContent(i));

         // zero crossings
         if (prevsign * hratio->GetBinContent(i) <= 0) {
	   if (prevsign != 0) ++zerox;
	   prevsign = hratio->GetBinContent(i);
         }
       } // fit range
     } // for i

     if (!quiet) cout << "Found " << zerox << " zero crossings" << endl;

     // TMath::NormQuantile to turn probability to sigma
     // TMath::Erfc to turn sigma into probability
     // p_global = 1 - pow(1-p_local, nbins)
     // where instead of nbins should use zero crossings?
     double p_local = TMath::Erfc(sigmalocal/sqrt(2));
     double p_global = 1 - pow(1-p_local, nbins);
     double sigmaglobal = -TMath::NormQuantile(p_global*0.5);

     tex->SetTextSize(0.045*1.7);
     if (sigmatot)
       tex->DrawLatex(0.20,0.83,Form("Max local (global) significance"
				     " %1.1f#sigma (%1.1f#sigma)",
				     sigmalocal, sigmaglobal));

  }

  ////////////////////////////////
  // Draw fit covariance matrix //
  ////////////////////////////////

  TH1D *hc = new TH1D("hc",";Parameter;Parameter",np,0,np);
  hc->SetMinimum(0);
  hc->SetMaximum(np);

  TCanvas *cc = tdrCanvas("cc",hc,4,0,kSquare);
  cc->SetRightMargin(0.15);

  // Add safety checks for number of fit parameters changing
  TH2D *h2c = new TH2D("h2c",";Parameter;Parameter",np,0,np,np,0,np);
  for (int i = 1; i != h2c->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2c->GetNbinsY()+1; ++j) {
      h2c->SetBinContent(i,j, emat[i-1][i-1]*emat[j-1][j-1]!=0 ?
			 emat[i-1][j-1] / sqrt(emat[i-1][i-1]*emat[j-1][j-1])
			 : 0);
    }
  }
  
  cc->cd();
  gStyle->SetOptStat(0);
  h2c->Draw("SAME COLZ");
  h2c->GetZaxis()->SetRangeUser(-1,+1);
  gPad->RedrawAxis();

  if (savecov) cc->SaveAs(Form("pdf/dijetmass_cov%s.pdf",ctag));
  } // covariance - quiet / Mjj spectrum quiet

  //////////////////////////
  // Draw trigger turn-on //
  //////////////////////////

  // Chi2 calculation outside quiet region
  double chi2_3(0);

  TF1 *fpass = new TF1("fpass","([0]*pow(1-x/13000.,[1])/pow(x/13000.,[2]"
		       "+[3]*log(x/13000.) ))*"
		       "(1./2.*(1+TMath::Erf((x-[4])/[5])))",
		       fitmin, fitmax); // 4-par + '2-par' trigger
  for (int i = 0; i != fpass->GetNpar(); ++i)
    fpass->SetParameter(i, f1->GetParameter(i));

  TF1 *ftot = new TF1("ftot","([0]*pow(1-x/13000.,[1])/pow(x/13000.,[2]"
		      "+[3]*log(x/13000.) ))",
		      fitmin, fitmax); // 4-par + '2-par' trigger
  for (int i = 0; i != ftot->GetNpar(); ++i)
    ftot->SetParameter(i, f1->GetParameter(i));

  TH1D *htrgeff_fit3 = (TH1D*)htrgeff->Clone("htrgeff_fit3");
  for (int i = 1; i != htrgeff->GetNbinsX()+1; ++i) {
    double mjjmin = htrgeff->GetBinLowEdge(i);
    double mjjmax = htrgeff->GetBinLowEdge(i+1);
    double pass = fpass->Integral(mjjmin, mjjmax);
    double tot = ftot->Integral(mjjmin, mjjmax);
    htrgeff_fit3->SetBinContent(i, pass / tot);
  }

  TH1D *htrgratio3 = (TH1D*)htrgeff->Clone("htrgratio3");
  for (int i = 1; i != htrgratio3->GetNbinsX()+1; ++i) {
    htrgratio3->SetBinContent(i, htrgeff->GetBinError(i)!=0 ? 
			     (htrgeff->GetBinContent(i)
			      - htrgeff_fit3->GetBinContent(i))
			     / (htrgeff->GetBinError(i)) : 0);
    htrgratio3->SetBinError(i, 0);
    chi2_3 += pow(htrgratio3->GetBinContent(i),2);
  }


  TH1D *htrgstat(0), *hefferr(0);
  if (!quiet) { // trigger

  // https://root.cern.ch/doc/master/classTEfficiency.html

  TH1D *h1 = new TH1D("h1",";Dijet Mass [GeV];Trigger efficiency",
		      int(890-453),453,890);
  h1->SetMaximum(1.2);
  h1->SetMinimum(0.2+1e-4);
  TH1D *h2 = new TH1D("h2",";Dijet Mass [GeV];(Data-Fit)/#sigma",
		      int(890-453),453,890);
  h2->SetMaximum(+3.8);
  h2->SetMinimum(-3.8);

  // Re-fit turn-on
  TF1 *ftrg = new TF1("ftrg","([0]/2) * (1 + TMath::Erf((x-[1])/[2]))",
		      453, 890);
  ftrg->SetParameters(500,1,500,100);
  ftrg->FixParameter(0,1);
  TF1 *ftrg2 = new TF1("ftrg2",
		       "(1-[4])*(1./2.) * (1 + TMath::Erf((x-[0])/[1])) +"
		       "[4]*(1./2.) * (1 + TMath::Erf((x-[2])/[3]))",
		       453, 890);
  ftrg2->SetParameters(500,100, 500,100,0);
  ftrg2->FixParameter(4,0);
  ftrg->SetLineWidth(3);
  TF1 *ftrg3 = new TF1("ftrg3","(1./2.) * (1 + TMath::Erf((x-[0])/[1]))",
  		       453, 890);
  ftrg3->SetParameters(500,100);
  ftrg->SetLineWidth(2);

  gStyle->SetOptFit(0);
  efftrg->Fit(ftrg, quiet ? "QRN" : "RN");

  // Store efficiency in TGraphAsymmErrors and pulls to TH1D for plotting
  TH1D *htrgs = (TH1D*)h1->Clone("htrgs");
  TGraph *gtrgs = new TGraph(htrgs->GetNbinsX());
  TGraph *gtrgs3 = new TGraph(htrgs->GetNbinsX());
  TGraphAsymmErrors *gtrg = new TGraphAsymmErrors(htrgs->GetNbinsX());

  for (int i = 1; i != htrgs->GetNbinsX()+1; ++i) {

    double x = htrgs->GetBinCenter(i);
    double ex = htrgs->GetBinWidth(i);
    int j = efftrg->FindFixBin(x);
    double y = efftrg->GetEfficiency(j);
    double eyl = efftrg->GetEfficiencyErrorLow(j);
    double eyh = efftrg->GetEfficiencyErrorUp(j);
    double fy = ftrg->Eval(x);
    double ey = (fy<y ? eyl : eyh);
    htrgs->SetBinContent(i, (y - fy) / ey);
    htrgs->SetBinError(i, 0);
    int n = gtrg->GetN();
    gtrg->SetPoint(n, x, y);
    gtrg->SetPointError(n, ex, ex, eyl, eyh);
    int m = gtrgs->GetN();
    gtrgs->SetPoint(m, x, (y-fy) / ey);
  }

  // Optional, shorter way to retrieve the TGraphAsymmErrors
  //efftrg->Draw(); gPad->Update();
  //gtrg = efftrg->GetPaintedGraph();

  gtrg->Fit(ftrg2, quiet ? "QRN" : "RN");
  ftrg3->SetParameters(f1->GetParameter(np-3),f1->GetParameter(np-2));
  double tchi2(0); int tndf(0);
  for (int i = 0; i != gtrg->GetN(); ++i) {
    double x = gtrg->GetX()[i];
    double y = gtrg->GetY()[i];
    double fy = ftrg3->Eval(x);
    double ey = (fy>y ? gtrg->GetErrorYhigh(i) : gtrg->GetErrorYlow(i));
    if (x>trgfitmin && x<trgfitmax) {
      tchi2 += pow((gtrg->GetY()[i] - fy) / ey, 2);
      ++tndf;
    }
    int m = gtrgs3->GetN();
    gtrgs3->SetPoint(m, x, (y-fy) / ey);
  }
  

  const int npeff = ftrg2->GetNpar();
  TMatrixD emateff(npeff, npeff);
  gMinuit->mnemat(emateff.GetMatrixArray(), npeff);

  // Statistical uncertaity band
  TF1 *fkeff = new TF1("fitErrorEff", fitError, trgfitmin, trgfitmax, 1);
  _fitError_func = ftrg2;
  _fitError_emat = &emateff;

  // Statistical uncertainty band for the trigger turnon fit
  hefferr = (TH1D*)hmjj->Clone("hefferr"); hefferr->Clear();
  fkeff->SetParameter(0,+1);
  for (int i = 1; i != hefferr->GetNbinsX()+1; ++i) {
    double x = hefferr->GetBinCenter(i);
    hefferr->SetBinContent(i, 0.);
    hefferr->SetBinError(i, 100.*fabs(fkeff->Eval(x)-ftrg2->Eval(x)));
  }


  TCanvas *ctrg = (dicanvas ? tdrDiCanvas("ctrg", h1, h2, 4, 11) :
		   tdrCanvas("ctrg", h1, 4, 11, kSquare));
  
  if (blind){
    
    if (dicanvas) {
      ctrg->cd(2);
      // Down panel       
      l->SetLineStyle(kDashed);
      l->SetLineWidth(2);
      l->SetLineColor(kRed);
      l->DrawLine(blindlow,-2.5,blindlow,2.5);
      l->DrawLine(blindhigh,-2.5,blindhigh,2.5);
    }

    ctrg->cd(1);
    // Top panel
    l->DrawLine(blindlow,0.90,blindlow,1.1);
    l->DrawLine(blindhigh,0.90,blindhigh,1.1);
    
  }

  if (dicanvas) ctrg->cd(1);
  tdrDraw(gtrg,"Pz",kFullCircle,kGray);

  // Calculate fitted efficiency over wide bins
  // taking properly into account the exponential spectrum
  TF1 *fpass2 = new TF1("fpass2","([0]*pow(1-x/13000.,[1])/pow(x/13000.,[2]"
			"+[3]*log(x/13000.) ))*"
			"((1-[8])*(1./2.*(1+TMath::Erf((x-[4])/[5]))) + "
			"[8]*(1./2.*(1+TMath::Erf((x-[6])/[7]))))",
			fitmin, fitmax); // 4-par + '2-par' trigger
  for (int i = 0; i != fpass2->GetNpar(); ++i)
    fpass2->SetParameter(i, (i<f1->GetNpar() ? f1->GetParameter(i) : 0.));

  TH1D *htrgeff_fit2 = (TH1D*)htrgeff->Clone("htrgeff_fit2");
  fpass2->SetParameter(4, ftrg2->GetParameter(0));
  fpass2->SetParameter(5, ftrg2->GetParameter(1));
  fpass2->SetParameter(6, ftrg2->GetParameter(2));
  fpass2->SetParameter(7, ftrg2->GetParameter(3));
  fpass2->SetParameter(8, ftrg2->GetParameter(4));
  for (int i = 1; i != htrgeff->GetNbinsX()+1; ++i) {
    double mjjmin = htrgeff->GetBinLowEdge(i);
    double mjjmax = htrgeff->GetBinLowEdge(i+1);
    double pass2 = fpass2->Integral(mjjmin, mjjmax);
    double tot = ftot->Integral(mjjmin, mjjmax);
    htrgeff_fit2->SetBinContent(i, pass2 / tot);
  }
  //
  TH1D *htrgeff_fit1 = (TH1D*)htrgeff->Clone("htrgeff_fit1");
  fpass->SetParameter(4, ftrg->GetParameter(1));
  fpass->SetParameter(5, ftrg->GetParameter(2));
  for (int i = 1; i != htrgeff->GetNbinsX()+1; ++i) {
    double mjjmin = htrgeff->GetBinLowEdge(i);
    double mjjmax = htrgeff->GetBinLowEdge(i+1);
    double pass = fpass->Integral(mjjmin, mjjmax);
    double tot = ftot->Integral(mjjmin, mjjmax);
    htrgeff_fit1->SetBinContent(i, pass / tot);
  }

  // Draw rebinned efficiency as well
  htrgeff->SetLineWidth(2);
  tdrDraw(htrgeff,"Pz",kDot,kBlack);

  // Draw fitted efficiency
  htrgeff_fit1->SetLineWidth(1);
  tdrDraw(htrgeff_fit1,"HE",kNone,kRed,kSolid,kRed,kNone);
  htrgeff_fit2->SetLineWidth(1);
  tdrDraw(htrgeff_fit2,"HE",kNone,kGreen+1,kSolid,kGreen+1,kNone);
  htrgeff_fit3->SetLineWidth(1);
  tdrDraw(htrgeff_fit3,"HE",kNone,kBlue,kSolid,kBlue,kNone);

  ftrg->SetLineColor(kRed);
  ftrg->SetLineWidth(2);
  ftrg->Draw("SAME");
  ftrg2->SetLineColor(kGreen+1);
  ftrg2->SetLineWidth(2);
  ftrg2->Draw("SAME");
  ftrg3->SetLineColor(kBlue);
  ftrg3->SetLineWidth(2);
  ftrg3->Draw("SAME");

  gPad->RedrawAxis();
  

  if(dicanvas) ctrg->cd(2);
  tdrDraw(gtrgs,"P",kFullCircle,kGray);
  tdrDraw(gtrgs3,"P",kFullCircle,kBlue-9);
  gtrgs3->SetMarkerSize(0.5);

  int ndfr(0); double chi2_1(0), chi2_2(0);
  // Draw ratio of fitted and rebinned efficiencies
  TH1D *htrgratio1 = (TH1D*)htrgeff->Clone("htrgratio1");
  htrgstat = (TH1D*)htrgeff->Clone("htrgstat");
  for (int i = 1; i != htrgratio1->GetNbinsX()+1; ++i) {
    htrgratio1->SetBinContent(i, htrgeff->GetBinError(i)!=0 ? 
			     (htrgeff->GetBinContent(i)
			      - htrgeff_fit1->GetBinContent(i))
			     / (htrgeff->GetBinError(i)) : 0);
    htrgratio1->SetBinError(i, 0);
    chi2_1 += pow(htrgratio1->GetBinContent(i),2);
    if (htrgratio1->GetBinContent(i)!=0) ++ndfr;

    htrgstat->SetBinError(i, htrgeff->GetBinError(i)*100.);
    htrgstat->SetBinContent(i, 0);
    int j = hmjj->FindBin(htrgeff->GetBinCenter(i));
  }
  tdrDraw(htrgratio1,"HE",kNone,kRed,kSolid,kRed,kNone);
  //
  TH1D *htrgratio2 = (TH1D*)htrgeff->Clone("htrgratio2");
  for (int i = 1; i != htrgratio2->GetNbinsX()+1; ++i) {
    htrgratio2->SetBinContent(i, htrgeff->GetBinError(i)!=0 ? 
			     (htrgeff->GetBinContent(i)
			      - htrgeff_fit2->GetBinContent(i))
			     / (htrgeff->GetBinError(i)) : 0);
    htrgratio2->SetBinError(i, 0);
    chi2_2 += pow(htrgratio2->GetBinContent(i),2);
  }
  tdrDraw(htrgratio2,"HE",kNone,kGreen+1,kSolid,kGreen+1,kNone);
  tdrDraw(htrgratio3,"HE",kNone,kBlue,kSolid,kBlue,kNone);

  gPad->RedrawAxis();

  tex->SetTextColor(kRed);
  tex->DrawLatex(0.60,0.81,Form("#chi^{2} / NDF = %1.1f / %d", chi2_1, ndfr));
  tex->SetTextColor(kGreen+1);
  tex->DrawLatex(0.60,0.73,Form("#chi^{2} / NDF = %1.1f / %d", chi2_2, ndfr));
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.60,0.65,Form("#chi^{2} / NDF = %1.1f / %d", chi2_3, ndfr));
				

  ctrg->cd(1);
  tex->SetNDC(); tex->SetTextSize(0.035);

  tex->SetTextColor(kRed);
  tex->DrawLatex(0.50,0.60,"TEfficiency::Fit()");
  tex->DrawLatex(0.50,0.56,Form("#chi^{2} / NDF = %1.1f / %d",
				ftrg->GetChisquare(), ftrg->GetNDF()));
  tex->DrawLatex(0.50,0.52,Form("p_{0} = %1.4f #pm %1.4f",
				ftrg->GetParameter(0), ftrg->GetParError(0)));
  tex->DrawLatex(0.50,0.48,Form("p_{1} = %1.1f #pm %1.1f",
				ftrg->GetParameter(1), ftrg->GetParError(1)));
  tex->DrawLatex(0.50,0.44,Form("p_{2} = %1.1f #pm %1.1f",
  				ftrg->GetParameter(2), ftrg->GetParError(2)));

  tex->SetTextColor(kGreen+1);
  tex->DrawLatex(0.50,0.38,"TGraphAsymmErrors::Fit()");
  tex->DrawLatex(0.50,0.34,Form("#chi^{2} / NDF = %1.1f / %d",
				ftrg2->GetChisquare(), ftrg2->GetNDF()));
  tex->DrawLatex(0.50,0.30,Form("p_{0} = %1.1f #pm %1.4f",
				ftrg2->GetParameter(0), ftrg2->GetParError(0)));
  tex->DrawLatex(0.50,0.26,Form("p_{1} = %1.1f #pm %1.1f",
				ftrg2->GetParameter(1), ftrg2->GetParError(1)));

  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.50,0.20, combofit ? "Simultaneous dijet mass fit" :
		 "Separate dijet mass fit");
  tex->DrawLatex(0.50,0.16,Form("#chi^{2} / NDF = %1.1f / %d",tchi2,tndf));
  tex->DrawLatex(0.50,0.12,Form("p_{0} = %1.1f #pm %1.1f",
				f1->GetParameter(np-3), f1->GetParError(np-3)));
  tex->DrawLatex(0.50,0.08,Form("p_{1} = %1.1f #pm %1.1f",
				f1->GetParameter(np-2), f1->GetParError(np-2)));
  // Range
  tex->DrawLatex(0.50,0.04,Form("Fit %1.0f < p_{T} < %1.0f GeV",
				fitmin, fitmax));

  tex->SetTextColor(kGray+1);
  tex->DrawLatex(0.28,0.16, "Reference fit:");
  tex->DrawLatex(0.28,0.12,Form("p_{0} = %1.1f",fsigout->GetParameter(np-3)));
  tex->DrawLatex(0.28,0.08,Form("p_{1} = %1.1f",fsigout->GetParameter(np-2)));

  ctrg->cd(0);
  tex->SetTextColor(kRed); tex->SetTextSize(0.045);
  if (isToy) tex->DrawLatex(0.50, 0.87, Form("Toy MC (%s)",ctoy));


  ctrg->SaveAs(Form("pdf/dijetmass_trg%s.pdf",ctag));
  } // trigger - quiet


  // Global chi2 and ndf in analysis bins
  Double_t chi2_gbl = chi2 + chi2_3;
  int ndf_gbl = (fitSig ? -nFitBtS : -nFitBtS+1); // parameters (6+1)
  int ndf_mjj = ndf_gbl;
  for (int i = 1; i != hmjj->GetNbinsX()+1; ++i) {
    double x = hmjj->GetBinCenter(i);
    if (x>fitmin && x<fitmax &&
	!(blind && x>blindlow && x<blindhigh) ) {
      ++ndf_gbl; // spectrum (23)
      ++ndf_mjj;
    }
    if (x>trgfitmin && x<trgfitmax) ++ndf_gbl; // trigger turnon (7)
  }
  if (!quiet) {

  cout << "Global chi^2/ndf = " << chi2_gbl<<" / "<<ndf_gbl << endl;

  c1->cd(1);

  tex->SetTextColor(kBlue);
  if (plotchi2gbl)
    tex->DrawLatex(0.56,0.28,Form("(%1.2f / %d)", chi2_gbl, ndf_gbl));

  if (!quiet) c1->SaveAs(Form("pdf/dijetmass_mjj%s.pdf",ctag));
  } // global chi2 - quiet
  

  if (!quiet) { // statistical uncertainty

  TH1D *h3 = new TH1D("h3",";Dijet Mass [GeV];Statistical uncertainty (%)",
		      int(890-453),453,890);
  h3->SetMaximum(+1.5);
  h3->SetMinimum(-1);
  TCanvas *cstat = tdrCanvas("cstat", h3, 4, 11, kSquare);

  TF1 *ftrigratio = new TF1("ftrgratio", "100.*((1+TMath::Erf((x-[0])/[1])) / "
			    "(1+TMath::Erf((x-[2])/[3])) - 1)",
			    trgfitmin,trgfitmax);
  ftrigratio->SetParameters(f1->GetParameter(np-3),f1->GetParameter(np-2),
			    fsigout->GetParameter(np-3),
			    fsigout->GetParameter(np-2));
  TF1 *ffitratio = new TF1("ffitratio",
			   "100.*(([0]*pow(1-x/13000.,[1])/pow(x/13000.,[2]"
			   "+[3]*log(x/13000.) ))*"
			   "(1./2.*(1+TMath::Erf((x-[4])/[5]))) /"
			   "(([6]*pow(1-x/13000.,[7])/pow(x/13000.,[8]"
			   "+[9]*log(x/13000.) ))*"
			   "(1./2.*(1+TMath::Erf((x-[10])/[11])))) - 1)",
			   trgfitmin, trgfitmax);
  for (int i = 0; i != f1->GetNpar()-1; ++i) 
    ffitratio->SetParameter(i, f1->GetParameter(i));
  for (int i = 0; i != fsigout->GetNpar()-1; ++i)
    ffitratio->SetParameter(f1->GetNpar()-1+i, fsigout->GetParameter(i));

  // Blind region
  if (blind) {
    l->SetLineStyle(kDashed);
    l->SetLineWidth(2);
    l->SetLineColor(kRed);
    l->DrawLine(blindlow,-0.5,blindlow,0.8);
    l->DrawLine(blindhigh,-0.5,blindhigh,0.8);
  }

  tdrDraw(htrgstat,"E2",kNone,kYellow+1,kSolid,-1,1001,kYellow+1);
  tdrDraw(hefferr,"E2",kNone,kOrange+3,kSolid,-1,1001,kOrange+1);
  tdrDraw(hscstat,"E2",kNone,kGreen+1,kSolid,-1,kNone,kGreen+1);
  tdrDraw(he2r,"E2",kNone,kGreen+3,kSolid,-1,kNone,kGreen+3);

  ftrigratio->SetLineWidth(2);
  ftrigratio->SetLineColor(kCyan+1);
  ftrigratio->Draw("SAME");

  ffitratio->SetLineWidth(2);
  ffitratio->SetLineStyle(kDashed);
  ffitratio->SetLineColor(kCyan+3);
  ffitratio->Draw("SAME");

  hsigoy->SetLineWidth(2);
  tdrDraw(hsigoy,"H][",kNone,kBlue,kSolid,kBlue,kNone);

  TLegend *legs = tdrLeg(0.55,0.55,0.85,0.90);
  legs->SetTextSize(0.035);
  legs->AddEntry(htrgstat,"Trg. eff. (bin-by-bin)","FL");
  legs->AddEntry(hefferr,"#epsilon (sep. fit)","FL");
  legs->AddEntry(hscstat,"Data (bin-by-bin)","FL");
  legs->AddEntry(he2r,Form("B (%s fit %1.0f-%1.0f GeV)",
			   combofit ? "sim." : "sep.",
			   fitmin, fitmax), "FL");
  legs->AddEntry(ftrigratio,"Trig. eff. (fit/ref.)","L");
  legs->AddEntry(ffitratio,"B#times#epsilon (fit/ref.)","L");

  gPad->RedrawAxis();

  tex->SetTextColor(kRed); tex->SetTextSize(0.045);
  if (isToy) tex->DrawLatex(0.50, 0.20, Form("Toy MC (%s)",ctoy));

  if (!quiet) cstat->SaveAs(Form("pdf/dijetmass_stat%s.pdf",ctag));
  } // statistical uncertainty - quiet


  if (quiet) { // cleanup for loops

    f->Close();
    fq->Close();
    fs->Close();
    ft->Close();
  }

  ensemble t;
  t.chi2_tot = chi2_gbl;
  t.chi2_mjj = chi2;
  t.ndf_tot = ndf_gbl;
  t.ndf_mjj = ndf_mjj;  
  t.xsec = (fitSig ? f1->GetParameter(np-1) : 0);
  t.xsecerr = (fitSig ? f1->GetParError(np-1) : 0);
  t.xsectrue = (TString(ctoy).Contains("sbtoy") ? 15. : 0.);
  
  return t;
}


Double_t fitError(Double_t *xx, Double_t *pp) {

  assert(_fitError_func);
  assert(_fitError_emat);
  double x = *xx;
  double k = pp[0];
  TF1 *f = _fitError_func;
  int n = f->GetNpar();
  TMatrixD &emat = (*_fitError_emat);
  assert(emat.GetNrows()==n);
  assert(emat.GetNcols()==n);
  
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {

    double p = f->GetParameter(i);
    double dp = 0.1*sqrt(emat[i][i]);
    f->SetParameter(i, p+dp);
    double fup = f->Eval(x);
    f->SetParameter(i, p-dp);
    double fdw = f->Eval(x);
    f->SetParameter(i, p);
    df[i] = (dp ? (fup - fdw) / (2.*dp) : 0);
  }

  double sumerr2(0);
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      sumerr2 += emat[i][j]*df[i]*df[j];
    }
  }

  double err = sqrt(sumerr2);

  return (f->Eval(x) + k*err);
}

void comboFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
		 Int_t flag) {

  assert(_fitObj);   // dijet mass spectrum
  assert(_fitObj2);  // dijet mass spectrum
  assert(_fitFunc);  // dijet mass fit function
  assert(_fitFunc2); // turn-on fit function
  assert(_fitFunc->GetNpar()==npar);

  if (flag) {

    int cnt(0);
    chi2 = 0;

    // TH1D fit
    if (_fitObj->InheritsFrom("TH1")) {
      
      TH1D *h = (TH1D*)_fitObj;
      TF1 *f = _fitFunc;
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	
	// Retrieve central value and uncertainty for this point
	double mass_dw = h->GetBinLowEdge(i);
	double mass_up = h->GetBinLowEdge(i+1);
	double mass = 0.5*(mass_dw+mass_up);
	double xsec = h->GetBinContent(i);
	double xsec_err = h->GetBinError(i);

	if (mass<fitmin || mass>fitmax ||
	    (blind && mass>blindlow && mass<blindhigh)) continue;
	
	// Calculate fit value at this point
	for (int ipar = 0; ipar != f->GetNpar(); ++ipar) {
	  f->SetParameter(ipar, par[ipar]);
	}
	double fit = f->Integral(mass_dw,mass_up)/(mass_up-mass_dw);
	
	// Add chi2 from residual
	double chi = (xsec - fit) / xsec_err;
	chi2 += chi * chi;
      } // for i

    } // fit histo
    
    // TGraphAsymmErrors fit
    if (_fitObj->InheritsFrom("TGraphAsymmErrors")) {
      
      TGraphAsymmErrors *g = (TGraphAsymmErrors*)_fitObj;
      TF1 *f = _fitFunc;
      for (int i = 0; i != g->GetN(); ++i) {
	
	// Retrieve central value and uncertainty for this point
	double mass = g->GetX()[i];
	double xsec = g->GetY()[i];
	double xsec_eup = g->GetErrorYhigh(i);
	double xsec_edw = g->GetErrorYlow(i);

	if (mass<fitmin || mass>fitmax ||
	    (blind && mass>blindlow && mass<blindhigh)) continue;
	
	// Calculate fit value at this point
	for (int ipar = 0; ipar != f->GetNpar(); ++ipar) {
	  f->SetParameter(ipar, par[ipar]);
	}
	double fit = f->Eval(mass);//f->EvalPar(&mass,par);
	double xsec_err = (fit>xsec ? xsec_eup : xsec_edw);
	
	// Add chi2 from residual
	double chi = (xsec - fit) / xsec_err;
	chi2 += chi * chi;
      } // for i
    } // fit mass graph

    // Add turn-on to the chi2 for simultaneous fit
    if (combofit) {

      TGraphAsymmErrors *gt = (TGraphAsymmErrors*)_fitObj2;
      TF1 *ft = _fitFunc2;
      for (int i = 0; i != gt->GetN(); ++i) {
	
	// Retrieve central value and uncertainty for this point
	double mass = gt->GetX()[i];
	double xsec = gt->GetY()[i];
	double xsec_eup = gt->GetErrorYhigh(i);
	double xsec_edw = gt->GetErrorYlow(i);
	
	if (mass<trgfitmin || mass>trgfitmax) continue;
	
	// Calculate fit value at this point
	//for (int ipar = 0; ipar != ft->GetNpar(); ++ipar) {
	  //ft->SetParameter(ipar, par[npar-2+ipar]);
	for (int ipar = 0; ipar != ft->GetNpar()-1; ++ipar) {
	  ft->SetParameter(ipar, par[npar-3+ipar]);
	}
	double fit = ft->Eval(mass);//f->EvalPar(&mass,par);
	double xsec_err = (fit>xsec ? xsec_eup : xsec_edw);
	
	// Add chi2 from residual
	double chi = (xsec - fit) / xsec_err;
	chi2 += chi * chi;
      } // for i
    } // combofit
    
    // Give some feedback on progress in case loop gets stuck
    if ((++cnt)%1000==0) cout << "." << flush;
    
  } // if flag
  else {
    if (grad) {}; // suppress warning
    return;
  }

} // comboFitter
