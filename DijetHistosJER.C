// Purpose: calculate JER SF using MPFX method
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TLine.h"

#include "../jecsys2020/tdrstyle_mod15.C"

TH1D *getJER(TProfile2D* p2, TProfile2D *p2x,
	     double eta1, double eta2, TH1D **h1 = 0, TH1D **h1x = 0);

void DijetHistosJERs(string file, string dir);
void drawDijetHistosJER();
void drawDijetHistosJERtest();

// Process several directories in a uniform way
void DijetHistosJER() {
  /*
  DijetHistosJERs("rootfiles/jmenano_data_cmb_v21ul16.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_data_cmb_v21ul16.root","Dijet/JER");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","Dijet/JER");
  */
  DijetHistosJERs("rootfiles/jmenano_data_cmb_v22ul16.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v22ul16mg.root","Dijet2");
  DijetHistosJERs("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","Dijet2");
  //drawDijetHistosJER();
  drawDijetHistosJERtest();
}

// Update cmb.root to add JER results
void DijetHistosJERs(string file, string dir) {

  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile(file.c_str(),"UPDATE");
  assert(f && !f->IsZombie());

  f->cd(dir.c_str());
  TDirectory *d = gDirectory;

  curdir->cd();
  
  TProfile2D *p2m0 = (TProfile2D*)d->Get("p2m0"); assert(p2m0);
  TProfile2D *p2m0x = (TProfile2D*)d->Get("p2m0x"); assert(p2m0x);

  // Calculate reference JER
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

  for (int i = 1; i != h2jer->GetNbinsX()+1; ++i) {

    // Calculate JER2 (average of tag and probe)
    double eta = h2jer->GetXaxis()->GetBinCenter(i);
    TH1D *h1jer2(0), *h1m0s(0), *h1m0xs(0);
    h1jer2 = getJER(p2m0,p2m0x,eta,eta,&h1m0s,&h1m0xs);

    for (int j = 1; j != h2jer->GetNbinsY()+1; ++j) {

      double jer2 = h1jer2->GetBinContent(j);
      double err2 = h1jer2->GetBinError(j);
      h2jer2->SetBinContent(i, j, jer2);
      h2jer2->SetBinError(i, j, err2);
      
      h2m0s->SetBinContent(i, j, h1m0s->GetBinContent(j));
      h2m0s->SetBinError(i, j, h1m0s->GetBinError(j));
      h2m0xs->SetBinContent(i, j, h1m0xs->GetBinContent(j));
      h2m0xs->SetBinError(i, j, h1m0xs->GetBinError(j));

      // Calculate probe JER by subtracting tag JER
      double jer13 = h1jer13->GetBinContent(j);
      double err13 = h1jer13->GetBinError(j);
      double jer = sqrt(max(2.*jer2*jer2 - jer13*jer13,0.));
      double err = (jer ? sqrt(2.*pow(jer2*err2,2)+pow(jer13*err13,2))/jer : 0);
      h2jer->SetBinContent(i, j, jer);
      h2jer->SetBinError(i, j, err);
    } // for j
    
    delete h1jer2;
    delete h1m0s;
    delete h1m0xs;
  } // for i

  d->cd();

  // Core results
  h1jer13->Write(h1jer13->GetName(),TObject::kOverwrite);
  h2jer->Write(h2jer->GetName(),TObject::kOverwrite);
  
  // Extras
  h2jer2->Write(h2jer2->GetName(),TObject::kOverwrite);

  h1m0s13->Write(h1m0s13->GetName(),TObject::kOverwrite);
  h1m0xs13->Write(h1m0xs13->GetName(),TObject::kOverwrite);
  h2m0s->Write(h2m0s->GetName(),TObject::kOverwrite);
  h2m0xs->Write(h2m0xs->GetName(),TObject::kOverwrite);

  f->Write();
  curdir->cd();
} // DijetHistosJERs


void drawDijetHistosJER() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/jmenano_data_cmb_v21ul16.root","READ");
  assert(f && !f->IsZombie());

  TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  lumi_13TeV = "2016GH, 16.8 fb^{-1}";
  const char *c = "_Dijet2";
  TCanvas *c1 = new TCanvas(Form("c1%s",c),Form("c1%s",c),1200,600);
  c1->Divide(6,3,0,0);
  
  TCanvas *c2 = new TCanvas(Form("c2%s",c),Form("c2%s",c),1200,600);
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

    TH1D *h = tdrHist(Form("h%d%s",i,c),"JER",0,0.3);
    tdrDraw(h,"",kNone);

    // RMS(MPFX) as low pT reference (sensitive to PU)
    TH1D *h1xs = h2m0xs->ProjectionY(Form("h1xs%d%s",i,c),i,i);
    TH1D *h1xsm = h2m0xsm->ProjectionY(Form("h1xsm%d%s",i,c),i,i);
    // Calculate effective parallel jet areas for MET (nmet):
    // int_0_2pi 1./(2pi)*cos^2(theta)dtheta = 1/2 to scale phi direction area
    double nmet = ((1./2.)*2.*5.*TMath::TwoPi()) / (TMath::Pi()*0.4*0.4);
    h1xs->Scale(1./sqrt(nmet));
    h1xsm->Scale(1./sqrt(nmet));
    h1xs->SetLineWidth(2);
    h1xs->GetXaxis()->SetRangeUser(15,200);
    tdrDraw(h1xs,"HIST",kNone,kRed,kSolid,-1,kNone);
    h1xsm->GetXaxis()->SetRangeUser(15,200);
    tdrDraw(h1xsm,"HIST",kNone,kRed,kSolid,-1,kNone);
    
    h1jer13->SetLineWidth(2);
    tdrDraw(h1jer13,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(h1jer13m,"HIST",kNone,kGreen+2,kSolid,-1,kNone);

    TH1D *h1jer = h2jer->ProjectionY(Form("h1jer%d%s",i,c),i,i);
    tdrDraw(h1jer,"Pz",kFullCircle,kGreen+2);
    h1jer->SetMarkerSize(0.7);
    TH1D *h1jerm = h2jerm->ProjectionY(Form("h1jerm%d%s",i,c),i,i);
    tdrDraw(h1jerm,"Pz",kOpenCircle,kGreen+2);
    h1jerm->SetMarkerSize(0.7);

    double eta1 = h2jer->GetXaxis()->GetBinLowEdge(i);
    double eta2 = h2jer->GetXaxis()->GetBinLowEdge(i+1);
    tex->DrawLatex(0.60,0.88,Form("%1.3f < |#eta| < %1.3f",eta1,eta2));

    double ptmax = min(2000.,6500.*0.7/cosh(eta1));
    TF1 *f1 = new TF1(Form("f1%d%s",i,c),"sqrt([0]*fabs([0])/(x*x)+"
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

    if (i%6==0 || i==1) {
      TLegend *leg = tdrLeg(0.70,0.85-4*0.05,1.00,0.85);
      leg->AddEntry(h1jer,"DATA","PLE");
      leg->AddEntry(h1jerm,"MC","PLE");
      leg->AddEntry(h1jer13,"|#eta|<1.3","L");
      leg->AddEntry(h1xs,"Est. PU","L");
    }
    
    c2->cd(min(i,18));
    gPad->SetLogx();
    
    TH1D *h2 = tdrHist(Form("h2%d%s",i,c),"Data/MC",0.5,2.0);
    tdrDraw(h2,"",kNone);
    l->DrawLine(15,1,3500,1);

    // RMS(MPF) as reference (less statistical uncertainty)
    //TH1D *h1m0s = h2m0s->ProjectionY(Form("h1m0s%d%s",i,c),i,i);
    //TH1D *h1m0sm = h2m0sm->ProjectionY(Form("h1m0sm%d%s",i,c),i,i);
    //TH1D *h1m0sr = (TH1D*)h1m0s->Clone(Form("h1m0sr%d%s",i,c));
    //h1m0sr->Divide(h1m0sm);
    //tdrDraw(h1m0sr,"HIST",kNone,kBlue,kSolid,-1,kNone);

    // RMS(MPFX) as low pT reference (sensitive to PU)
    TH1D *h1m0xs = h2m0xs->ProjectionY(Form("h1m0xs%d%s",i,c),i,i);
    TH1D *h1m0xsm = h2m0xsm->ProjectionY(Form("h1m0xsm%d%s",i,c),i,i);
    TH1D *h1m0xsr = (TH1D*)h1m0xs->Clone(Form("h1m0xsr%d%s",i,c));
    h1m0xsr->Divide(h1m0xsm);
    h1m0xsr->GetXaxis()->SetRangeUser(15,200);
    tdrDraw(h1m0xsr,"HIST",kNone,kRed,kSolid,-1,kNone);

    
    TH1D *h1jer13r = (TH1D*)h1jer13->Clone("h1jer13r");
    h1jer13r->Divide(h1jer13m);
    h1jer13r->SetLineWidth(2);
    tdrDraw(h1jer13r,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
    
    TH1D *h1jerr = (TH1D*)h1jer->Clone("h1jerr");
    h1jerr->Divide(h1jerm);
    tdrDraw(h1jerr,"Pz",kFullCircle,kGreen+2);
    
    tex->DrawLatex(0.60,0.88,Form("%1.3f < |#eta| < %1.3f",eta1,eta2));

    
    TF1 *f2 = new TF1(Form("f2%d%s",i,c),
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

    if (i%6==0 || i==1) {
      TLegend *leg = tdrLeg(0.70,0.85-4*0.05,1.00,0.85);
      leg->AddEntry(h1jer,"DATA","PLE");
      leg->AddEntry(h1jerm,"MC","PLE");
      leg->AddEntry(h1jer13,"|#eta|<1.3","L");
      leg->AddEntry(h1m0xsr,"Est. PU","L");
    }
  } // for i

  c1->SaveAs("pdf/DijetHistosJER_JERMC.pdf");
  c2->SaveAs("pdf/DijetHistosJER_JERSF.pdf");
} // drawDijetHistosJER

void drawDijetHistosJERtest() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/jmenano_data_cmb_v22ul16.root","READ");
  //TFile *f = new TFile("rootfiles/jmenano_data_cmb_v21ul16.root","READ");
  //TFile *f = new TFile("rootfiles/jmenano_data_cmb.root","READ");
  assert(f && !f->IsZombie());

  TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v22ul16mg.root","READ");
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","READ");
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","READ");
  //TFile *fm = new TFile("rootfiles/jmenano_mc_cmb.root","READ");
  assert(fm && !fm->IsZombie());

  TFile *fp = new TFile("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","READ");
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

  TH1D *h1jer(0), *h1m0s(0), *h1m0xs(0);
  h1jer = getJER(p2m0,p2m0x,0,1.3,&h1m0s,&h1m0xs);
  TH1D *h1jerm(0), *h1m0sm(0), *h1m0xsm(0);
  p2m0m->SetName("p2m0m");
  h1jerm = getJER(p2m0m,p2m0xm,0,1.3,&h1m0sm,&h1m0xsm);
  TH1D *h1jerp(0), *h1m0sp(0), *h1m0xsp(0);
  p2m0p->SetName("p2m0p");
  h1jerp = getJER(p2m0p,p2m0xp,0,1.3,&h1m0sp,&h1m0xsp);

  TH1D *h1bal(0), *h1m2s(0), *h1m2xs(0);
  h1bal = getJER(p2m2,p2m2x,0,1.3,&h1m2s,&h1m2xs);
  TH1D *h1balm(0), *h1m2sm(0), *h1m2xsm(0);
  p2m2m->SetName("p2m2m");
  h1balm = getJER(p2m2m,p2m2xm,0,1.3,&h1m2sm,&h1m2xsm);
  TH1D *h1balp(0), *h1m2sp(0), *h1m2xsp(0);
  p2m2p->SetName("p2m2p");
  h1balp = getJER(p2m2p,p2m2xp,0,1.3,&h1m2sp,&h1m2xsp);

  TH1D *h1fsr(0), *h1mnus(0), *h1mnuxs(0);
  h1fsr = getJER(p2mnu,p2mnux,0,1.3,&h1mnus,&h1mnuxs);
  TH1D *h1fsrm(0), *h1mnusm(0), *h1mnuxsm(0);
  p2mnum->SetName("p2mnum");
  h1fsrm = getJER(p2mnum,p2mnuxm,0,1.3,&h1mnusm,&h1mnuxsm);
  TH1D *h1fsrp(0), *h1mnusp(0), *h1mnuxsp(0);
  p2mnup->SetName("p2mnup");
  h1fsrp = getJER(p2mnup,p2mnuxp,0,1.3,&h1mnusp,&h1mnuxsp);

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
}

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
    double rmsy = p1s->GetBinError(i) / sqrt(2.);
    double erry = p1->GetBinError(i) / sqrt(2.);
    h1s->SetBinContent(i, rmsy);
    h1s->SetBinError(i, rmsy ? erry / rmsy : 0);
    double rmsx = p1xs->GetBinError(i) / sqrt(2.);
    double errx = p1x->GetBinError(i) / sqrt(2.);
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
}
