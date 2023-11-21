// Purpose: calculate JEC L2Res using MPF (update to HMD) method
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "../jecsys2020/tdrstyle_mod15.C"

TH1D *getHDM(TProfile2D* p2, TProfile2D *p22, TProfile2D *p2n, TProfile2D *p2u,
	     double eta1, double eta2, string sb,
	     TH1D **h1mpf, TH1D **h1db);

void DijetHistosL2Ress(string file, string dir, string bin="");

// Process several directories in a uniform way
void DijetHistosL2Res() {

  /*
  // Data before final L2Res
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2016APV_v26c.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2016GH_v26c.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2017_v26.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2018_v26c.root","Dijet2");
  //
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_Run2_v26c.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_Run2_v26c.root","Dijet2","pf");
  DijetHistosL2Ress("haddfiles/jmenano_data_cmb_Run2_v26c.root","Dijet2","tc");
  //DijetHistosL2Ress("haddfiles/jmenano_data_cmb_UL2017_v26.root","Dijet2");

  // MC before JER SFs
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v26.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v26.root","Dijet2","pf");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v26.root","Dijet2","tc");
  //DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_UL2017MG_v26.root","Dijet2");

  // MC after JER SF
  DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_UL2016APVMG_v27.root","Dijet2");
  DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_UL2016MG_v27.root","Dijet2");
  DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_UL2017MG_v27.root","Dijet2");
  DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_UL2018MG_v27.root","Dijet2");
  //
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v27.root","Dijet2");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v27.root","Dijet2","pf");
  DijetHistosL2Ress("haddfiles/jmenano_mc_cmb_Run2_v27.root","Dijet2","tc");
  //DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_UL2017MG_v27.root","Dijet2");
  */

  // Process all binning variants (pTavp, pTprobe, pTtag)
  for (int i = 0; i != 3; ++i) {
    
    const char *c = (i==0 ? "" : (i==1 ? "pf" : "tc"));

    // Run3 (v29->v30)
    DijetHistosL2Ress("../rootfiles/jmenano_data_cmb_2022CD_JME_v30.root","Dijet2",c);
    DijetHistosL2Ress("../rootfiles/jmenano_data_cmb_2022E_JME_v30.root","Dijet2",c);
    DijetHistosL2Ress("../rootfiles/jmenano_data_cmb_2022FG_JME_v30.root","Dijet2",c);
    // 2023
    DijetHistosL2Ress("../rootfiles/jmenano_data_cmb_2023BCv123_JME_v30.root","Dijet2",c);
    DijetHistosL2Ress("../rootfiles/jmenano_data_cmb_2023Cv4_JME_v30.root","Dijet2",c);
    DijetHistosL2Ress("../rootfiles/jmenano_data_cmb_2023D_JME_v30.root","Dijet2",c);
    // MC
    DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_Summer22MG_v30.root","Dijet2",c);

    /*
    // New Run3 files from Iita and Mikael
    // 2022
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/jmenano_data_cmb_2022C_v1.root","Dijet2",c);
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/jmenano_data_cmb_2022D_v1.root","Dijet2",c);
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/jmenano_data_cmb_2022E_v1.root","Dijet2",c);
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/jmenano_data_cmb_2022F_v1.root","Dijet2",c);
    //DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/jmenano_data_cmb_2022G_v1.root","Dijet2",c);
    // 2023
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/nano_data_cmb_2023B_v1.root","Dijet2",c);
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/nano_data_cmb_2023Cv123_v1.root","Dijet2",c);
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/nano_data_cmb_2023Cv4_v1.root","Dijet2",c);
    DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/nano_data_cmb_2023D_v1.root","Dijet2",c);
    // Main combos (after checking stability for L2Res)
    //DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/jmenano_data_cmb_2022D_v1.root","Dijet2",c);
    //DijetHistosL2Ress("../jecsys3/rootfiles/Iita_20230814/nano_data_cmb_2022BCv123_v1.root","Dijet2",c);
    */
    // UL2018 reference
    DijetHistosL2Ress("../rootfiles/jmenano_mc_cmb_UL2018MG_v26.root","Dijet2",c);
  }
  
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
  
  // Create 2D histograms for storing HDM JES
  TH2D *h2jes = p2m0->ProjectionXY(Form("h2jes%s",cb)); h2jes->Reset();
  TH2D *h2mpf = p2m0->ProjectionXY(Form("h2mpf%s",cb)); h2mpf->Reset();
  TH2D *h2db = p2m0->ProjectionXY(Form("h2db%s",cb)); h2db->Reset();

  for (int i = 1; i != h2jes->GetNbinsX()+1; ++i) {

    // Calculate HDM (ratio of tag and probe)
    double eta = h2jes->GetXaxis()->GetBinCenter(i);
    TH1D *h1jec(0), *h1mpf(0), *h1db(0);
    h1jec = getHDM(p2m0,p2m2,p2mn,p2mu,eta,eta,cb,&h1mpf,&h1db);

    for (int j = 1; j != h2jes->GetNbinsY()+1; ++j) {

      double jec = h1jec->GetBinContent(j);
      double err = h1jec->GetBinError(j);
      h2jes->SetBinContent(i, j, jec);
      h2jes->SetBinError(i, j, err);

      if (h1mpf) {
	h2mpf->SetBinContent(i, j, h1mpf->GetBinContent(j));
	h2mpf->SetBinError(i, j, h1mpf->GetBinError(j));
      }
      if (h1db) {
	h2db->SetBinContent(i, j, h1db->GetBinContent(j));
	h2db->SetBinError(i, j, h1db->GetBinError(j));
      }
    } // for j
    
    delete h1jec;
  } // for i

  d->cd();

  h2jes->Write(Form("h2jes%s",cb),TObject::kOverwrite);
  h2mpf->Write(Form("h2mpf%s",cb),TObject::kOverwrite);
  h2db->Write(Form("h2db%s",cb),TObject::kOverwrite);
  
  f->Write(nullptr,TObject::kOverwrite);
  curdir->cd();  

  delete h2jes;

} // DijetHistosL2Ress

// Calculate HDM JEC for one given eta slice at a time
TF1 *_f1mpf(0), *_f1db(0);
TF1 *_f1mpf_ptavp(0), *_f1db_ptavp(0);
TF1 *_f1mpf_pttag(0), *_f1db_pttag(0);
TF1 *_f1mpf_ptprobe(0), *_f1db_ptprobe(0);
TH1D *getHDM(TProfile2D* p2, TProfile2D *p22, TProfile2D *p2n, TProfile2D *p2u,
	     double eta1, double eta2, string sb,
	     TH1D **h1mpf, TH1D **h1db) {
  
  string s = Form("_%s",p2->GetName());
  const char *c = s.c_str();

  int i1 = p2->GetXaxis()->FindBin(eta1);
  int i2 = p2->GetXaxis()->FindBin(eta2);
  double etamin = p2->GetXaxis()->GetBinLowEdge(i1);
  double etamax = p2->GetXaxis()->GetBinLowEdge(i2+1);

  TProfile *p1 = p2->ProfileY(Form("p1%s",c),i1,i2,"");
  TProfile *p12 = p22->ProfileY(Form("p12%s",c),i1,i2,"");
  TProfile *p1n = p2n->ProfileY(Form("p1n%s",c),i1,i2,"");
  TProfile *p1u = p2u->ProfileY(Form("p1u%s",c),i1,i2,"");

  // Calculate HDM JEC
  TH1D *h1jec = p1->ProjectionX(Form("h1jec_%s",c));
  if (h1mpf) *h1mpf = p1->ProjectionX(Form("h1mpf_%s",c));
  if (h1db)  *h1db  = p1->ProjectionX(Form("h1db_%s",c));
  for (int i = 1; i != h1jec->GetNbinsX()+1; ++i) {

    double mpf = p1->GetBinContent(i);
    double err = p1->GetBinError(i);

    double mpf0 = p1->GetBinContent(i);
    double mpf2 = p12->GetBinContent(i);
    double mpfn = p1n->GetBinContent(i);
    double mpfu = p1u->GetBinContent(i);

    // Sanity check on HDM inputs, it's critical this holds
    if (fabs(mpf0 - (mpf2+mpfn+mpfu)) > 1e-4) {
      cout << Form("(mpf0=%1.4g) != (mpf2+mpfn+mpfu=%1.4g) "
		   "at eta=[%1.3f,%1.3f], where\n"
		   "mpf2=%1.4g, mpfn=%1.4g, mpfu=%1.4g\n",
		   mpf0,  mpf2+mpfn+mpfu,
		   etamin,etamax, mpf2,mpfn,mpfu) << flush;
    }

    // MPF-21-001_temp.pdf Eq.(21),(22):
    // (R1-R2)/(R1+R2) = r0 + ((R1+R2)/(2Rn)-1)*rn + ((R1+R2)/(2Ru)-1)*ru, or
    // (R1-R2)/(R1+R2) = r2 + ((R1+R2)/(2Rn))*rn + ((R1+R2)/(2Ru))*ru
    // We want to solve for R1 (probe) over R2(tag), R1/R2, thus R=R1/R2:
    // (x-1)/(x+1) = r0 + ((x+1)/(2*(Rn/R2))-1)*rn + ((x+1)/(2*(Ru/R2))-1)*ru,
    // (x-1)/(x+1) = r2 + ((x+1)/(2*(Rn/R2)))*rn + ((x+1)/(2*(Ru/R2)))*ru
    //
    // Revising this formula now, starting from
    // 0 = pt + pp + pn +pu,
    // 2 = pt - pp,   (or could we re-normalize 2a=pt-pp?)
    // rX = Rx*px/a
    // a = (Rt*pt-Rp*pp)/2 ~ (Rt+Rp) + O((Rt-Rp)*(pt+pp)
    // => 2(Rt-Rp)/(Rt+Rp) = (MPF-1) -((Rt+Rp)/(2Rn)-1)*rn -((Rt+Rp)/(2Ru)-1)*ru
    double R2(1), Rn(1), Ru(0.92), eps(1e-5);
    if (_f1mpf_ptavp==0 && sb=="") {
      // [0]=MPF, [1]=DB, [2]=MPFn, [3]=MPFu, [4]=R2, [5]=Rn, [6]=Ru
      _f1mpf_ptavp = new TF1("_f1mpf_ptavp",
			     //"(x-1)/(x+1) - [0]"
			     //"- ((x+1)/(2*[5]/[4])-1)*[2]"
			     //"- ((x+1)/(2*[6]/[4])-1)*[3]",
			     //"2*(x-1)/(x+1) - ([0]-1)"
			     "(x-1)/(x+1) - ([0]-1)"
			     "+ ((x+1)/(2*[5]/[4])-1)*[2]"
			     "+ ((x+1)/(2*[6]/[4])-1)*[3]",
			     0.,100.);
    }
    // MPF-21-001_temp.pdf Eq. (14),(15) with explicit R2
    // R1/Rz = r0 + (R1/Rn-1)*rn + (R1/Ru-1)*ru
    // R1/Rz = r2 + (R1/Rn-0)*rn + (R1/Ru-0)*ru
    if (_f1mpf_pttag==0 && sb=="tc") {
      // [0]=MPF, [1]=DB, [2]=MPFn, [3]=MPFu, [4]=R2, [5]=Rn, [6]=Ru
      _f1mpf_pttag = new TF1("_f1mpf_pttag",
			     "x/[4] - [0]"
			     "- (x/[5]-1)*[2]"
			     "- (x/[6]-1)*[3]",
			     0.,100.);
    }
    // Guesswork, derive formula properly
    if (_f1mpf_ptprobe==0 && sb=="pf") {
      // [0]=MPF, [1]=DB, [2]=MPFn, [3]=MPFu, [4]=R2, [5]=Rn, [6]=Ru
      _f1mpf_ptprobe = new TF1("_f1mpf_ptprobe",
			       "(1-[4]/x) - [0]"
			       "- (([4]/x)/[5]-1)*[2]"
			       "- (([4]/x)/[6]-1)*[3]",
			       0.001,100.);
    }
    // Extension to pTtag (not double-checked yet):
    // R1/R2 = r0 + (R2/Rn-1)*rn + (R2/Ru-1)*ru, i.e.
    // x = r0 + (R2/Rn-1)*rn + (R2/Ru-1)*ru
    if (sb=="") { _f1mpf = _f1mpf_ptavp; }
    if (sb=="tc") { _f1mpf = _f1mpf_pttag; }
    if (sb=="pf") { _f1mpf = _f1mpf_ptprobe; }
    assert(_f1mpf);
    //_f1mpf->SetParameters(mpf0-1,mpf2-1,mpfn,mpfu, R2,Rn,Ru);
    _f1mpf->SetParameters(mpf0,mpf2,mpfn,mpfu, R2,Rn,Ru);
    double hdm_mpf = (mpf0!=0 ? _f1mpf->GetX(0.,0.,100.,eps) : 0);
    if (_f1db_ptavp==0 && sb=="") {
      // [0]=MPF, [1]=DB, [2]=Rn, [3]=Ru, [4]=R2, [5]=Rn, [6]=Ru
      _f1db_ptavp = new TF1("_f1db_ptavp",
			    //"(x-1)/(x+1) - [1]"
			    //"- ((x+1)/(2*[5]/[4]-0))*[2]"
			    //"- ((x+1)/(2*[6]/[4]-0))*[3]",
			    //"2.*(x-1)/(x+1) - ([1]-1)"
			     "(x-1)/(x+1) - ([1]-1)"
			     "+ ((x+1)/(2*[5]/[4])-0)*[2]"
			     "+ ((x+1)/(2*[6]/[4])-0)*[3]",
			    0.,100.);
    }
    if (_f1db_pttag==0 && sb=="tc") {
      // [0]=MPF, [1]=DB, [2]=MPFn, [3]=MPFu, [4]=R2, [5]=Rn, [6]=Ru
      _f1db_pttag = new TF1("_f1db_pttag",
			    "x/[4] - [0]"
			    "- (x/[5]-0)*[2]"
			    "- (x/[6]-0)*[3]",
			    0.,100.);
    }
    // Guesswork, derive formula properly
    if (_f1db_ptprobe==0 && sb=="pf") {
      // [0]=MPF, [1]=DB, [2]=MPFn, [3]=MPFu, [4]=R2, [5]=Rn, [6]=Ru
      _f1db_ptprobe = new TF1("_f1db_ptprobe",
			      "(1-[4]/x) - [0]"
			      "- (([4]/x)/[5]-1)*[2]"
			      "- (([4]/x)/[6]-1)*[3]",
			      0.001,100.);
    }
    if (sb=="") { _f1db = _f1db_ptavp; }
    if (sb=="tc") { _f1db = _f1db_pttag; }
    if (sb=="pf") { _f1db = _f1db_ptprobe; }
    if (!_f1db) {
      cout << "Missing _f1db for sb="<<sb<<endl<<flush;
      assert(_f1db);
    }
    //_f1db->SetParameters(mpf0-1,mpf2-1,mpfn,mpfu, R2,Rn,Ru);
    _f1db->SetParameters(mpf0,mpf2,mpfn,mpfu, R2,Rn,Ru);
    double hdm_db = (mpf2!=0 ? _f1db->GetX(0,0.,100.,eps) : 0);

    // Sanity check on HDM outputs, these should agree fully
    if (fabs(hdm_mpf - hdm_db) > 1e-4) {
      cout << Form("(hdm_mpf0=%1.4g) != (hdm_db=%1.4g) "
		   "at eta=[%1.3f,%1.3f], where\n"
		   "mpf0=%1.4g, mpf2=%1.4g, mpfn=%1.4g, mpfu=%1.4g\n",
		   hdm_mpf,hdm_db,etamin,etamax,mpf0,mpf2,mpfn,mpfu) << flush;
    }

    // Add proper hdm calculation here, now just test with MPF
    double hdm = (hdm_mpf>0 && hdm_mpf<100 ? hdm_mpf : 0);
    double jec = hdm;
    //double jec = mpf;
    
    h1jec->SetBinContent(i, jec);
    h1jec->SetBinError(i, err);

    if (*h1mpf) {
      (*h1mpf)->SetBinContent(i, mpf0);
      (*h1mpf)->SetBinError(i, p1->GetBinError(i));
    }
    if (*h1db) {
      (*h1db)->SetBinContent(i, mpf2);
      (*h1db)->SetBinError(i, p12->GetBinError(i));
    }
  }

  //cout << h1 << ", " << (*h1) << endl;

  delete p1;
  delete p12;
  delete p1n;
  delete p1u;
  
  return h1jec;
} // getHDM
