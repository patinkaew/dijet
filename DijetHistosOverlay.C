// Purpose: overlay DijetHistos with DESY results
#include "TFile.h"
#include "../jecsys2020/tdrstyle_mod15.C"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLine.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include <string>

bool debug = false;

int findBin(TH2D *h2, double x, double *xnew = 0);
void addBins(TH1D *h1to, TH2D *h2from, int i1, int i2,
	     TH1D *h1count, TH2D *h2count);
void rebin(TH1D *h1, TH1D *h1n, TH1D *h1old, TH1D *h1nold);

//void DijetHistosOverlays(string obs, string data, string spt = "PtAVP",
//void DijetHistosOverlays(string obs, string data, string spt = "PtTag",
void DijetHistosOverlays(string obs, string data, string spt = "PtProbe",
//void DijetHistosOverlays(string obs, string data, string spt = "PtAve",
			 bool data3=false);
void DijetHistosOverlayPtBins(string obs);
void DijetHistosOverlayJER();

void DijetHistosOverlay() {

  /*
  DjetHistosOverlays("MPF","data");
  DijetHistosOverlays("mpf2","data");
  DijetHistosOverlays("mpfN","data");
  DijetHistosOverlays("mpfU","data");

  DijetHistosOverlays("MPF","mc");
  DijetHistosOverlays("mpf2","mc");
  DijetHistosOverlays("mpfN","mc");
  DijetHistosOverlays("mpfU","mc");
  */
  //DijetHistosOverlays("MPF","mc");
  //DijetHistosOverlays("mpf2","mc");
  //DijetHistosOverlays("MPF","data");
  //DijetHistosOverlays("mpf2","data");

  //DijetHistosOverlays("mpf2","data",true);
  //DijetHistosOverlays("MPF","data",true);
  //DijetHistosOverlays("mpfU","data",true);

  /*
  DijetHistosOverlayPtBins("mpfU");
  DijetHistosOverlayPtBins("mpfN");
  DijetHistosOverlayPtBins("mpf2");
  DijetHistosOverlayPtBins("MPF");
  */

  DijetHistosOverlayJER();
} // DijetHistosOverlay

void DijetHistosOverlays(string obs, string data, string spt,
			 bool data3) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);
  TLine *l = new TLine();

  string s = Form("_%s_%s",data.c_str(),obs.c_str());
  const char *c = s.c_str();
  const char *cpt = spt.c_str();

  TFile *f1(0), *f12(0), *f13(0);
  if (data=="data") f1 = new TFile("rootfiles/jmenano_data_cmb_v22ul16.root","READ");
  if (data=="mc") f1 = new TFile("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","READ");
  //if (data=="data") f1 = new TFile("rootfiles/jmenano_data_cmb_v21ul16.root","READ");
  //if (data=="mc") f1 = new TFile("rootfiles/jmenano_mc_cmb_v20ul16flatmc.root","READ");
  //if (data=="data") f1 = new TFile("rootfiles/jmenano_data_cmb.root","READ");
  //if (data=="mc") f1 = new TFile("rootfiles/jmenano_mc_cmb.root","READ");

  if (data3) {
    f1 = new TFile("rootfiles/dijet2_a_dijet_cmb.root","READ");
    f12 = new TFile("rootfiles/dijet2_b_asymm_cmb.root","READ");
    f13 = new TFile("rootfiles/dijet2_c_allgood_cmb.root","READ");
    assert(f12 && !f12->IsZombie());
    assert(f13 && !f13->IsZombie());
  }
  assert(f1 && !f1->IsZombie());

  TProfile2D *p2(0), *p22(0), *p23(0);
  TH2D *h2n(0), *h2n2(0), *h2n3(0);
  if (spt=="PtAVP" || spt=="PtAve") {
    /*
      if (obs=="MPF") p2 = (TProfile2D*)f1->Get("Dijet/p2m0ab");
      if (obs=="mpf2") p2 = (TProfile2D*)f1->Get("Dijet/p2m2ab");
      if (obs=="mpfN") p2 = (TProfile2D*)f1->Get("Dijet/p2mnab");
      if (obs=="mpfU") p2 = (TProfile2D*)f1->Get("Dijet/p2muab");
    */
    if (obs=="MPF") p2 = (TProfile2D*)f1->Get("Dijet2/p2m0");
    if (obs=="mpf2") p2 = (TProfile2D*)f1->Get("Dijet2/p2m2");
    if (obs=="mpfN") p2 = (TProfile2D*)f1->Get("Dijet2/p2mn");
    if (obs=="mpfU") p2 = (TProfile2D*)f1->Get("Dijet2/p2mu");
    h2n = (TH2D*)f1->Get("Dijet2/h2pteta"); assert(h2n); // Dijet2
  }
  if (spt=="PtTag") {
    if (obs=="MPF") p2 = (TProfile2D*)f1->Get("Dijet2/p2m0tc");
    if (obs=="mpf2") p2 = (TProfile2D*)f1->Get("Dijet2/p2m2tc");
    if (obs=="mpfN") p2 = (TProfile2D*)f1->Get("Dijet2/p2mntc");
    if (obs=="mpfU") p2 = (TProfile2D*)f1->Get("Dijet2/p2mutc");
    h2n = (TH2D*)f1->Get("Dijet2/h2ptetatc"); assert(h2n); // Dijet2
  }
  if (spt=="PtProbe") {
    if (obs=="MPF") p2 = (TProfile2D*)f1->Get("Dijet2/p2m0pf");
    if (obs=="mpf2") p2 = (TProfile2D*)f1->Get("Dijet2/p2m2pf");
    if (obs=="mpfN") p2 = (TProfile2D*)f1->Get("Dijet2/p2mnpf");
    if (obs=="mpfU") p2 = (TProfile2D*)f1->Get("Dijet2/p2mupf");
    h2n = (TH2D*)f1->Get("Dijet2/h2ptetatc"); assert(h2n); // Dijet2
  }      
  assert(p2);
  assert(h2n);
  //TH2D *h2 = p2->ProjectionXY(Form("h2%s",c));
  //TH2D *h2n = (TH2D*)f1->Get("Dijet/h2pteta_aball"); assert(h2n); // Dijet
  //h2n = (TH2D*)f1->Get("Dijet2/h2pteta"); assert(h2n); // Dijet2
  
  if (data3) {
    p22 = (TProfile2D*)f12->Get(Form("Dijet2/%s",p2->GetName())); assert(p22);
    p23 = (TProfile2D*)f13->Get(Form("Dijet2/%s",p2->GetName())); assert(p23);
    h2n2 = (TH2D*)f12->Get("Dijet2/h2pteta"); assert(h2n2);
    h2n3 = (TH2D*)f13->Get("Dijet2/h2pteta"); assert(h2n3);
  }
  
  
  // TFile *f2 = new TFile("../jecsys2020/rootfiles/CombinationFiles-Run2016FGH-3.root","READ");
  //TFile *f2 = new TFile(Form("rootfiles/CombinationFiles-Run2016FGH-%s.root",cpt),"READ");
  TFile *f2 = new TFile(Form("rootfiles/CombinationFiles-Run2016FGH-%s.root",spt=="PtAve" ? "PtAVP" : cpt),"READ");
  assert(f2 && !f2->IsZombie());

  f2->cd(data.c_str());//"data");
  TDirectory *d2 = gDirectory;

  curdir->cd();

  TCanvas *c1 = new TCanvas(Form("c1%s",c),Form("c1%s",c),1200,600);
  c1->Divide(6,3,0,0);

  TCanvas *c2 = new TCanvas(Form("c2%s",c),Form("c2%s",c),1200,600);
  c2->Divide(6,3,0,0);

  int neta(0);
  TIter next(d2->GetListOfKeys());
  while (TKey *key = (TKey*)next()) {

    // Recurse directory structure
    if (!(string(key->GetClassName())=="TDirectoryFile")) continue;
    if (debug) cout << key->GetName() << "->";
    TDirectory *d2s = (TDirectory*)key->ReadObj();
    d2s->cd();
    
    int ietamin, ietamax;
    sscanf(key->GetName(),"eta_%d_%d",&ietamin,&ietamax);

    c1->cd(++neta);
    gPad->SetLogx();

    TH1D *h = tdrHist(Form("h%s_%d",c,neta),Form("%s (%s)",obs.c_str(),data.c_str()),0.98,1.07,"p_{T,avp} (GeV)");
    if (spt=="PtAVP") {
      if (obs=="MPF" || obs=="mpf2") {
	if (neta>=1)  h->GetYaxis()->SetRangeUser(0.97,1.07);
	//if (neta>=7)  h->GetYaxis()->SetRangeUser(0.98,1.08);
	if (neta>=7)  h->GetYaxis()->SetRangeUser(0.90,1.10);
	if (neta>=13) h->GetYaxis()->SetRangeUser(0.80,1.25);
      }
      if (obs=="mpfN" || obs=="mpfU") {
	if (neta>=1) h->GetYaxis()->SetRangeUser(-0.03,0.05);
	if (neta>=7) h->GetYaxis()->SetRangeUser(-0.05,0.15);
	if (neta>=13) h->GetYaxis()->SetRangeUser(-0.05,0.25);
      }
    }
    if (spt=="PtTag" || spt=="PtProbe" || spt=="PtAve") {
      if (obs=="MPF" || obs=="mpf2") {
	if (neta>=1)  h->GetYaxis()->SetRangeUser(0.6,1.4);
	if (neta>=7)  h->GetYaxis()->SetRangeUser(0.6,1.4);
	if (neta>=13) h->GetYaxis()->SetRangeUser(0.6,1.4);
      }
      if (obs=="mpfN" || obs=="mpfU") {
	if (neta>=1) h->GetYaxis()->SetRangeUser(-0.3,0.3);
	if (neta>=7) h->GetYaxis()->SetRangeUser(-0.3,0.3);
	if (neta>=13) h->GetYaxis()->SetRangeUser(-0.3,0.3);
      }
      if (spt=="PtTag") h->SetXTitle("p_{T,tag} (GeV)");
      if (spt=="PtProbe") h->SetXTitle("p_{T,probe} (GeV)");
      if (spt=="PtAve") h->SetXTitle("p_{T,avp} (GeV)");
    }

    TH1D *h1c = (TH1D*)d2s->Get(Form("%s_RawNEvents_a100",data.c_str()));
    assert(h1c);
    TGraphErrors *g0(0);
    if (obs=="MPF") g0 = (TGraphErrors*)d2s->Get("mpfchs_dijet_a100");
    if (obs=="mpf2") g0 = (TGraphErrors*)d2s->Get("mpfchs2_dijet_a100");
    if (obs=="mpfN") g0 = (TGraphErrors*)d2s->Get("mpfchsN_dijet_a100");
    if (obs=="mpfU") g0 = (TGraphErrors*)d2s->Get("mpfchsU_dijet_a100");
    assert(g0);

    h->Draw();
    l->DrawLine(15,1,3500,1);
    //g0->SetLineWidth(2);
    g0->SetLineColor(kRed);
    g0->SetMarkerColor(kRed);
    g0->Draw("SAMEPz");

    double etamin;
    int i1m = findBin(p2, -0.1*ietamin);
    int i1p = findBin(p2, +0.1*ietamin, &etamin);

    double etamax;
    int i2m = findBin(p2, -0.1*ietamax);
    int i2p = findBin(p2, +0.1*ietamax, &etamax);

    //int nbins = (i2p-i1p) + (i1m-i2m); // Dijet
    int nbins = (i2p-i1p); // Dijet2

    //TH1D *h0m = p2->ProjectionY(Form("h0m_%d",neta),i2m,i1m-1);
    //TH1D *h0p = p2->ProjectionY(Form("h0p_%d",neta),i1p,i2p-1);
    //TH1D *h0 = (TH1D*)h0m->Clone(Form("h0_%d",neta));
    //h0->Add(h0p);
    //h0->Scale(1./nbins);
    //h0->Draw("SAME");
						
    TH1D *h0 = p2->ProjectionY(Form("h0m%s_%d",c,neta),-1,-1);
    TH1D *h1 = (TH1D*)h0->Clone(Form("h1%s_%d",c,neta)); h1->Reset();
    TH1D *h1n = (TH1D*)h0->Clone(Form("h1n%s_%d",c,neta)); h1n->Reset();
    //addBins(h1,p2,i2m,i1m-1,h1n,h2n); // Dijet only
    addBins(h1,p2,i1p,i2p-1,h1n,h2n); // Dijet2+Dijet
    h1->SetLineColor(kGreen+2);
    //h1->Draw("SAMEH");

    //doublx xr[] = {15, 21, 28, 37, 49,
    //		     59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575,
    //		     638, 737, 846, 967, 1101, 1248,
    //               1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103};
    //double xb[] = {59, 85, 104, 170, 236, 302, 370, 460, 575, 1000, 2000};
    double xb[] = {15,28, 59, 86, 110, 170, 236, 302, 373, 460, 575,
		   1000, 2000, 3103};
    double nb = sizeof(xb)/sizeof(xb[0])-1;
    //double xf[] = {59, 86, 110, 132, 204, 279, 373};
    double xf[] = {59, 86, 110, 132, 204, 279, 373,  460, 575,
		   1000, 2000, 3103};
    int nf = sizeof(xf)/sizeof(xf[0])-1;
    double *x = (etamin>2.8 ? xf : xb);
    int n = (etamin>2.8 ? nf : nb);
    TH1D *h1old = h1;
    TH1D *h1nold = h1n;
    h1 = new TH1D(Form("h1new%s_%d",c,neta),"",n,x);
    h1n = new TH1D(Form("h1nnew%s_%d",c,neta),"",n,x);
    rebin(h1,h1n,h1old,h1nold);
    h1->SetLineColor(kGreen+2);
    h1->SetMarkerStyle(kNone);
    h1->Draw("SAMEH");

    if (data3) {
      TH1D *h12 = p22->ProjectionY(Form("h12_%s_%d",c,neta),i1p,i2p-1);//-1,-1);
      TH1D *h13 = p23->ProjectionY(Form("h13_%s_%d",c,neta),i1p,i2p-1);//-1,-1);
      h12->SetLineColor(kBlue);
      h12->SetLineStyle(kDotted);
      h12->Draw("SAMEH");
      h13->SetLineColor(kRed);
      h13->Draw("SAMEH");
    }

    tex->DrawLatex(0.60,0.92,Form("eta_%02d_%02d (%d)",ietamin,ietamax,nbins));
    tex->DrawLatex(0.60,0.86,Form("%1.3f<|#eta|<%1.3f",etamin,etamax));

    if (neta==1 || neta%6==0) {
      TLegend *leg = tdrLeg(0.60,0.83-2*0.06,0.90,0.83);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(g0,"DESY mini","PLE");
      leg->AddEntry(h1,"JME nano","PLE");
    }	      

    // Jet counts
    c2->cd(neta);
    gPad->SetLogx();
    gPad->SetLogy();

    TH1D *h2 = tdrHist(Form("h2%s_%d",c,neta),Form("%s (%s)","Counts",data.c_str()),0.5,1e7,"p_{T,avp} (GeV)");    

    h2->Draw();
    h1c->SetMarkerColor(kRed);
    h1c->SetLineColor(kRed);
    if (etamin<1.3) h1c->Scale(2.);
    h1c->Draw("SAMEH");
    h1n->SetMarkerStyle(kNone);
    h1n->SetLineColor(kGreen+2);
    h1n->Draw("SAMEHE");

    tex->DrawLatex(0.60,0.92,Form("eta_%02d_%02d (%d)",ietamin,ietamax,nbins));
    tex->DrawLatex(0.60,0.86,Form("%1.3f<|#eta|<%1.3f",etamin,etamax));

    gPad->RedrawAxis();

  } // while tkey

  c1->SaveAs(Form("pdf/DijetHistosOverlay_%s_%s_%s_%s.pdf","2016GH",
		  data.c_str(),obs.c_str(),cpt));
  c2->SaveAs(Form("pdf/DijetHistosOverlay_%s_%s_%s_%s.pdf","2016GH",
		  data.c_str(),"counts",cpt));
} // DijetHistosOverlay


void DijetHistosOverlayPtBins(string obs) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);
  TLine *l = new TLine();

  const char *co = obs.c_str();

  TFile *f1(0), *f1m(0), *f1p(0);
  f1 = new TFile("rootfiles/jmenano_data_cmb_v22ul16.root","READ");
  assert(f1 && !f1->IsZombie());
  f1m = new TFile("rootfiles/jmenano_mc_cmb_v22ul16mg.root","READ");
  assert(f1m && !f1m->IsZombie());
  f1p = new TFile("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","READ");
  assert(f1p && !f1p->IsZombie());

  TProfile2D *p2a(0), *p2t(0), *p2p(0);
  TProfile2D *p2am(0), *p2tm(0), *p2pm(0);
  TProfile2D *p2ap(0), *p2tp(0), *p2pp(0);
  //TH2D *h2na(0), *h2nt(0), *h2np(0);

  map<string,string> mh;
  mh["MPF"] = "Dijet2/p2m0";
  mh["mpf2"] = "Dijet2/p2m2";
  mh["mpfN"] = "Dijet2/p2mn";
  mh["mpfU"] = "Dijet2/p2mu";

  const char *ch = mh[obs].c_str();
  p2a = (TProfile2D*)f1->Get(Form("%s",ch));   assert(p2a);
  p2t = (TProfile2D*)f1->Get(Form("%stc",ch)); assert(p2t);
  p2p = (TProfile2D*)f1->Get(Form("%spf",ch)); assert(p2p);

  p2am = (TProfile2D*)f1m->Get(Form("%s",ch));   assert(p2am);
  p2tm = (TProfile2D*)f1m->Get(Form("%stc",ch)); assert(p2tm);
  p2pm = (TProfile2D*)f1m->Get(Form("%spf",ch)); assert(p2pm);

  p2ap = (TProfile2D*)f1p->Get(Form("%s",ch));   assert(p2ap);
  p2tp = (TProfile2D*)f1p->Get(Form("%stc",ch)); assert(p2tp);
  p2pp = (TProfile2D*)f1p->Get(Form("%spf",ch)); assert(p2pp);

  curdir->cd();

  TCanvas *c1 = new TCanvas(Form("c1%s",co),Form("c1%s",co),1200,600);
  c1->Divide(6,3,0,0);

  TCanvas *c2 = new TCanvas(Form("c2%s",co),Form("c2%s",co),1200,600);
  c2->Divide(6,3,0,0);

  for (int ieta = 1; ieta != p2a->GetNbinsX()+1; ++ieta) {
    
    c1->cd(ieta);
    gPad->SetLogx();
    
    TH1D *h = tdrHist(Form("h%s_%d",co,ieta),co,0.6,1.5);
    /*
    if (obs=="MPF")  h->GetYaxis()->SetRangeUser(0.8,1.4);
    if (obs=="mpf2") h->GetYaxis()->SetRangeUser(0.5,1.6);
    if (obs=="mpfN") h->GetYaxis()->SetRangeUser(-0.25,0.5);
    if (obs=="mpfU") h->GetYaxis()->SetRangeUser(-0.1,0.2);
    */
    // Common ranges
    if (obs=="MPF" || obs=="mpf2")  h->GetYaxis()->SetRangeUser(0.5,1.6);
    //if (obs=="mpfN" || obs=="mpfU") h->GetYaxis()->SetRangeUser(-0.25,0.5);
    if (obs=="mpfN" || obs=="mpfU") h->GetYaxis()->SetRangeUser(-0.5,0.6);
    h->Draw();

    l->SetLineStyle(kDotted);
    l->SetLineStyle(kGray+1);
    double ymin = h->GetMaximum();//GetYaxis()->GetXmin();
    double ymax = h->GetMinimum();//GetYaxis()->GetXmax();
    l->DrawLine(59,ymin,59,ymax);
    double etamin = p2a->GetXaxis()->GetBinLowEdge(ieta);
    double etamax = p2a->GetXaxis()->GetBinLowEdge(ieta+1);
    double ptmax = 0.5*6500/cosh(etamin);
    l->DrawLine(ptmax,ymin,ptmax,ymax);

    l->SetLineStyle(kDotted);
    l->SetLineStyle(kGray);
    double ptmax2 = 0.7*6500/cosh(etamin);
    l->DrawLine(ptmax2,ymin,ptmax2,ymax);

    l->SetLineStyle(kSolid);
    l->SetLineColor(kBlack);
    l->DrawLine(15,1,3500,1);
    
    TH1D *ha = p2a->ProjectionY(Form("ha%s_%d",co,ieta),ieta,ieta);
    TH1D *ht = p2t->ProjectionY(Form("ht%s_%d",co,ieta),ieta,ieta);
    TH1D *hp = p2p->ProjectionY(Form("hp%s_%d",co,ieta),ieta,ieta);

    TH1D *ham = p2am->ProjectionY(Form("ham%s_%d",co,ieta),ieta,ieta);
    TH1D *htm = p2tm->ProjectionY(Form("htm%s_%d",co,ieta),ieta,ieta);
    TH1D *hpm = p2pm->ProjectionY(Form("hpm%s_%d",co,ieta),ieta,ieta);

    TH1D *hap = p2ap->ProjectionY(Form("hap%s_%d",co,ieta),ieta,ieta);
    TH1D *htp = p2tp->ProjectionY(Form("htp%s_%d",co,ieta),ieta,ieta);
    TH1D *hpp = p2pp->ProjectionY(Form("hpp%s_%d",co,ieta),ieta,ieta);
    
    tdrDraw(htp,"HIST",kNone,kBlue+1,kDotted,-1,kNone);
    tdrDraw(hpp,"HIST",kNone,kRed+1,kDotted,-1,kNone);
    tdrDraw(hap,"HIST",kNone,kGreen+3,kDotted,-1,kNone);

    tdrDraw(htm,"HE",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(hpm,"HE",kNone,kRed,kSolid,-1,kNone);
    tdrDraw(ham,"HE",kNone,kGreen+2,kSolid,-1,kNone);

    tdrDraw(ht,"Pz",kFullDiamond,kBlue,kNone,-1,kNone);
    tdrDraw(hp,"Pz",kFullDiamond,kRed,kNone,-1,kNone);
    tdrDraw(ha,"Pz",kFullDiamond,kGreen+2,kNone,-1,kNone);
    
    //ht->SetMarkerSize(0.7);
    //hp->SetMarkerSize(0.7);
    //ha->SetMarkerSize(0.7);

    tex->DrawLatex(0.55,0.92,Form("%1.3f<|#eta|<%1.3f",etamin,etamax));

    if (ieta==2 || ieta%6==0) {
      TLegend *leg = tdrLeg(0.55,0.90-3*0.075,0.85,0.90);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(ht,"PtTag","PLE");
      leg->AddEntry(ha,"PtAVP","PLE");
      leg->AddEntry(hp,"PtProbe","PLE");
    }
    if (ieta==1 || ieta%6==5) {
      TLegend *leg = tdrLeg(0.55,0.90-3*0.075,0.85,0.90);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(ha,"Data","PLE");
      leg->AddEntry(ham,"MadGraph","PLE");
      leg->AddEntry(hap,"Pythia8","PL");
    }


    c2->cd(ieta);
    gPad->SetLogx();

    TH1D *h2 = tdrHist(Form("h2%s_%d",co,ieta),Form("Data-MC (%s)",co),-0.1,0.1);
    /*
    if (obs=="MPF")  h2->GetYaxis()->SetRangeUser(-0.05,0.10);
    if (obs=="mpf2") h2->GetYaxis()->SetRangeUser(-0.15,0.15);
    if (obs=="mpfN") h2->GetYaxis()->SetRangeUser(-0.15,0.15);
    if (obs=="mpfU") h2->GetYaxis()->SetRangeUser(-0.05,0.05);
    */
    // Common range
    h2->GetYaxis()->SetRangeUser(-0.15,0.15);
    h2->Draw();

    double ymin2 = h2->GetMaximum();
    double ymax2 = h2->GetMinimum();
    l->SetLineStyle(kDotted);
    l->SetLineStyle(kGray+1);
    l->DrawLine(59,ymin2,59,ymax2);
    l->DrawLine(ptmax,ymin2,ptmax,ymax2);
    l->SetLineStyle(kGray);
    l->DrawLine(ptmax2,ymin2,ptmax2,ymax2);
    l->SetLineStyle(kSolid);
    l->SetLineColor(kBlack);
    //l->DrawLine(15,0,3500,0);

    TH1D *hapr = (TH1D*)ha->Clone(Form("hapr%s_%d",co,ieta));
    TH1D *htpr = (TH1D*)ht->Clone(Form("htpr%s_%d",co,ieta));
    TH1D *hppr = (TH1D*)hp->Clone(Form("hppr%s_%d",co,ieta));

    hapr->Add(hap,-1);
    htpr->Add(htp,-1);
    hppr->Add(hpp,-1);

    htpr->GetXaxis()->SetRangeUser(59,ptmax2);
    hppr->GetXaxis()->SetRangeUser(59,ptmax2);
    hapr->GetXaxis()->SetRangeUser(59,ptmax2);

    // HISTP loses the histograms, so drwa HIST once at the back first
    tdrDraw(htpr,"HIST",kOpenDiamond,kBlue+1,kDotted,-1,kNone);
    tdrDraw(hppr,"HIST",kOpenDiamond,kRed+1,kDotted,-1,kNone);
    tdrDraw(hapr,"HIST",kOpenDiamond,kGreen+3,kDotted,-1,kNone);

    tdrDraw(htpr,"HISTP",kOpenDiamond,kBlue+1,kDotted,-1,kNone);
    tdrDraw(hppr,"HISTP",kOpenDiamond,kRed+1,kDotted,-1,kNone);
    tdrDraw(hapr,"HISTP",kOpenDiamond,kGreen+3,kDotted,-1,kNone);
    
    TH1D *har = (TH1D*)ha->Clone(Form("har%s_%d",co,ieta));
    TH1D *htr = (TH1D*)ht->Clone(Form("htr%s_%d",co,ieta));
    TH1D *hpr = (TH1D*)hp->Clone(Form("hpr%s_%d",co,ieta));

    har->Add(ham,-1);
    htr->Add(htm,-1);
    hpr->Add(hpm,-1);

    htr->GetXaxis()->SetRangeUser(59,ptmax2);
    hpr->GetXaxis()->SetRangeUser(59,ptmax2);
    har->GetXaxis()->SetRangeUser(59,ptmax2);

    tdrDraw(htr,"HEPz",kFullDiamond,kBlue,kSolid,-1,kNone);
    tdrDraw(hpr,"HEPz",kFullDiamond,kRed,kSolid,-1,kNone);
    tdrDraw(har,"HEPz",kFullDiamond,kGreen+2,kSolid,-1,kNone);

    tex->DrawLatex(0.55,0.92,Form("%1.3f<|#eta|<%1.3f",etamin,etamax));

    if (ieta==2 || ieta%6==0) {
      TLegend *leg = tdrLeg(0.55,0.90-3*0.075,0.85,0.90);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(htr,"PtTag","PLE");
      leg->AddEntry(har,"PtAVP","PLE");
      leg->AddEntry(hpr,"PtProbe","PLE");
    }
    if (ieta==1 || ieta%6==5) {
      TLegend *leg = tdrLeg(0.55,0.90-2*0.075,0.85,0.90);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(har,"MadGraph","PLE");
      leg->AddEntry(hapr,"Pythia8","PL");
    }
  } // for ieta

  c1->SaveAs(Form("pdf/DijetHistosOverlayPtBins_%s_%s.pdf","2016GH",co));
  c2->SaveAs(Form("pdf/DijetHistosOverlayPtBins_%s_%s_ratio.pdf","2016GH",co));
} // DijetHistosOverLayPtBins

void DijetHistosOverlayJER() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);
  TLine *l = new TLine();

  const char *co = "";
  
  TFile *f1(0), *f1m(0), *f1p(0);
  f1 = new TFile("rootfiles/jmenano_data_cmb_v22ul16.root","READ");
  assert(f1 && !f1->IsZombie());
  f1m = new TFile("rootfiles/jmenano_mc_cmb_v22ul16mg.root","READ");
  assert(f1m && !f1m->IsZombie());  
  f1p = new TFile("rootfiles/jmenano_mc_cmb_v22ul16flatmc.root","READ");
  assert(f1p && !f1p->IsZombie());

  curdir->cd();

  //TH1D *h1jer13(0), *h1jer13m(0), *h1jer13p(0);
  TH2D *h2jer(0), *h2jerm(0), *h2jerp(0);
  //h1jer13 = (TH1D*)f1->Get("Dijet2/h1jer13"); assert(h1jer13);
  //h1jer13m = (TH1D*)f1m->Get("Dijet2/h1jer"); assert(h1jer13m);
  //h1jer13p = (TH1D*)f1p->Get("Dijet2/h1jer"); assert(h1jer13p);
  h2jer = (TH2D*)f1->Get("Dijet2/h2jer"); assert(h2jer);
  h2jerm = (TH2D*)f1m->Get("Dijet2/h2jer"); assert(h2jerm);
  h2jerp = (TH2D*)f1p->Get("Dijet2/h2jer"); assert(h2jerp);
  TH2D *h2n(0), *h2nm(0), *h2np(0);
  h2n = (TH2D*)f1->Get("Dijet2/h2n"); assert(h2n);
  h2nm = (TH2D*)f1m->Get("Dijet2/h2n"); assert(h2nm);
  h2np = (TH2D*)f1p->Get("Dijet2/h2n"); assert(h2np);

  
  TCanvas *c1 = new TCanvas(Form("c1%s",co),Form("c1%s",co),1200,600);
  c1->Divide(6,3,0,0);

  TCanvas *c2 = new TCanvas(Form("c2%s",co),Form("c2%s",co),1200,600);
  c2->Divide(6,3,0,0);

  const double vy[] = {15, 30, 60, 120, 240, 480, 1000,  2000, 4000, 6500};
  const int ny = sizeof(vy)/sizeof(vy[0])-1;
  vector<double> vx(h2jer->GetNbinsX()+1);
  for (int i = 1; i != h2jer->GetNbinsX()+2; ++i) {
    vx[i-1] = h2jer->GetXaxis()->GetBinLowEdge(i);
  }
  int nx = vx.size()-1;
  TH2D *h2jerf = new TH2D("h2jerf",";#eta;JER",nx,&vx[0],ny,vy);
  TH2D *h2jersf = new TH2D("h2jersf",";#eta;JER SF",nx,&vx[0],ny,vy);
  
  for (int ieta = 1; ieta != h2jer->GetNbinsX()+1; ++ieta) {
    
    c1->cd(ieta);
    gPad->SetLogx();

    TH1D *h = tdrHist(Form("h%s_%d",co,ieta),"",0.,0.5);
    h->Draw();

    l->SetLineStyle(kDotted);
    l->SetLineStyle(kGray+1);
    double ymin = h->GetMaximum();
    double ymax = h->GetMinimum();
    l->DrawLine(59,ymin,59,ymax);
    double etamin = h2jer->GetXaxis()->GetBinLowEdge(ieta);
    double etamax = h2jer->GetXaxis()->GetBinLowEdge(ieta+1);
    double ptmax = 0.5*6500/cosh(etamin);
    l->DrawLine(ptmax,ymin,ptmax,ymax);

    l->SetLineStyle(kDotted);
    l->SetLineStyle(kGray);
    double ptmax2 = 0.7*6500/cosh(etamin);
    l->DrawLine(ptmax2,ymin,ptmax2,ymax);

    l->SetLineStyle(kSolid);
    l->SetLineColor(kBlack);
    l->DrawLine(15,1,3500,1);
    
    TH1D *h1jer = h2jer->ProjectionY(Form("h1jer%s_%d",co,ieta),ieta,ieta);
    TH1D *h1jerm = h2jerm->ProjectionY(Form("h1jerm%s_%d",co,ieta),ieta,ieta);
    TH1D *h1jerp = h2jerp->ProjectionY(Form("h1jerp%s_%d",co,ieta),ieta,ieta);

    h1jer->SetMarkerSize(0.7);
    
    tdrDraw(h1jerp,"HIST",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(h1jerm,"HEPz",kOpenCircle,kBlue,kSolid,-1,kNone);
    tdrDraw(h1jer,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone);
    
    tex->DrawLatex(0.55,0.92,Form("%1.3f<|#eta|<%1.3f",etamin,etamax));

    if (ieta==1 || ieta%6==0) {
      TLegend *leg = tdrLeg(0.55,0.90-3*0.075,0.85,0.90);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(h1jer,"Data","PLE");
      leg->AddEntry(h1jerm,"MadGraph","PLE");
      leg->AddEntry(h1jerp,"Pythia8","PLE");
    }

    // Add fits on top
    double fitptmax = min(2000.,0.7*6500./cosh(etamin));
    if (fabs(etamin)>=2.9)
      fitptmax = 0.5*6500./cosh(etamin);
    /*
    // v1 without separating out PU and "UE" noise
    TF1 *fjer = new TF1(Form("fjer%d",ieta),
			"sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])",
			30,fitptmax);
    TF1 *fjer2 = new TF1(Form("fjer2%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])",
			 30,fitptmax);
    TF1 *fjerm = new TF1(Form("fjerm%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])",
			 30,fitptmax);
    TF1 *fjerp = new TF1(Form("fjerp%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])",
			 30,fitptmax);

    fjerm->SetParameters(0,1,0.04,-1);
    fjerm->FixParameter(0,0);
    fjerm->SetParLimits(1,0.5,1.5);
    fjerm->SetParLimits(2,0.02,0.25);
    fjerm->SetParLimits(3,-1.25,-0.75);
    if (fabs(etamin)>2.6) {
      fjerm->FixParameter(3,-1);
    }    
    h1jerm->Fit(fjerm,"QRN");
    fjerm->SetRange(15,3500);

    fjerp->SetParameters(0,1,0.04,-1);
    fjerp->FixParameter(0,0);
    fjerp->SetParLimits(1,0.5,1.5);
    fjerp->SetParLimits(2,0.02,0.25);
    fjerp->SetParLimits(3,-1.25,-0.75);
    if (fabs(etamin)>2.6) {
      fjerp->FixParameter(3,-1);
    }    
    h1jerp->Fit(fjerp,"QRN");
    fjerp->SetRange(15,3500);

    fjer->SetParameters(fjerm->GetParameter(0),fjerm->GetParameter(1),
			fjerm->GetParameter(2),fjerm->GetParameter(3));
    fjer->FixParameter(0,fjerm->GetParameter(0));
    fjer->SetParLimits(1,fjerm->GetParameter(1),1.5);
    fjer->SetParLimits(2,fjerm->GetParameter(2),0.25);
    fjer->FixParameter(3,fjerm->GetParameter(3));
    h1jer->Fit(fjer,"QRN");
    fjer->SetRange(15,3500);

    fjer2->SetParameters(fjerp->GetParameter(0),fjerp->GetParameter(1),
			 fjerp->GetParameter(2),fjerp->GetParameter(3));
    fjer2->FixParameter(0,fjerp->GetParameter(0));
    fjer2->SetParLimits(1,fjerp->GetParameter(1),1.5);
    fjer2->SetParLimits(2,fjerp->GetParameter(2),0.25);
    fjer2->FixParameter(3,fjerp->GetParameter(3));
    h1jer->Fit(fjer2,"QRN");
    fjer2->SetRange(15,3500);
    */

    bool fixN0to0 = false;
    double fitptmin = 15;
    TF1 *fjer = new TF1(Form("fjer%d",ieta),
			"sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
			"[1]*[1]*pow(x,[3])+[2]*[2])",
			fitptmin,fitptmax);
    TF1 *fjer2 = new TF1(Form("fjer2%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
			 "[1]*[1]*pow(x,[3])+[2]*[2])",
			 fitptmin,fitptmax);
    TF1 *fjerm = new TF1(Form("fjerm%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
			 "[1]*[1]*pow(x,[3])+[2]*[2])",
			 fitptmin,fitptmax);
    TF1 *fjerp = new TF1(Form("fjerp%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
			 "[1]*[1]*pow(x,[3])+[2]*[2])",
			 fitptmin,fitptmax);
    

    double ptrcm = h2nm->GetYaxis()->GetBinCenter(1);
    double rcm = h2nm->GetBinContent(ieta, 1) * ptrcm;
    fjerm->SetParameters(0,1,0.04,-1,rcm);
    //if (!fixN0to0) fjerm->SetParLimits(0,-5,0);
    if (!fixN0to0) fjerm->SetParLimits(0,-2,+2);
    if (fixN0to0)  fjerm->FixParameter(0,0);
    fjerm->SetParLimits(1,0.5,1.5);
    fjerm->SetParLimits(2,0.02,0.25);
    fjerm->SetParLimits(3,-1.25,-0.75);
    fjerm->FixParameter(4,rcm);
    if (fabs(etamin)>2.6) {
      fjerm->FixParameter(0,0.);
      fjerm->FixParameter(3,-1);
    }    
    h1jerm->Fit(fjerm,"QRN");
    fjerm->SetRange(15,3500);

    double ptrc = h2n->GetYaxis()->GetBinCenter(1);
    double rcdt = h2n->GetBinContent(ieta, 1) * ptrc;
    fjer->SetParameters(fjerm->GetParameter(0),fjerm->GetParameter(1),
			fjerm->GetParameter(2),fjerm->GetParameter(3),
			rcdt);
    fjer->FixParameter(0,fjerm->GetParameter(0));
    fjer->SetParLimits(1,fjerm->GetParameter(1),1.5);
    fjer->SetParLimits(2,fjerm->GetParameter(2),0.25);
    fjer->SetParLimits(3,fjerm->GetParameter(3),-0.75);
    fjer->FixParameter(4,rcdt);
    if (fabs(etamin)>2.6) {
      //fjer->FixParameter(0,0.);
      //fjer->FixParameter(3,-1);
      fjer->FixParameter(3,fjerm->GetParameter(3));
    }    
    h1jer->Fit(fjer,"QRN");
    fjer->SetRange(15,3500);

    // Repeat for Pythia
    fjerp->SetParameters(0,1,0.04,-1,rcm);
    //if (!fixN0to0) fjerp->SetParLimits(0,-5,0);
    if (!fixN0to0) fjerp->SetParLimits(0,-2,+2);
    if (fixN0to0)  fjerp->FixParameter(0,0);
    fjerp->SetParLimits(1,0.5,1.5);
    fjerp->SetParLimits(2,0.02,0.25);
    fjerp->SetParLimits(3,-1.25,-0.75);
    fjerp->FixParameter(4,rcm);
    if (fabs(etamin)>2.6) {
      fjerp->FixParameter(0,0.);
      fjerp->FixParameter(3,-1);
    }    
    h1jerp->Fit(fjerp,"QRN");
    fjerp->SetRange(15,3500);
    //
    fjer2->SetParameters(fjerp->GetParameter(0),fjerp->GetParameter(1),
			 fjerp->GetParameter(2),fjerp->GetParameter(3),
			 rcdt);
    fjer2->FixParameter(0,fjerp->GetParameter(0));
    fjer2->SetParLimits(1,fjerp->GetParameter(1),1.5);
    fjer2->SetParLimits(2,fjerp->GetParameter(2),0.25);
    fjer2->SetParLimits(3,fjerp->GetParameter(3),-0.75);
    fjer2->FixParameter(4,rcdt);
    if (fabs(etamin)>2.6) {
      //fjer2->FixParameter(0,0.);
      //fjer2->FixParameter(3,-1);
      fjer2->FixParameter(3,fjerp->GetParameter(3));
    }    
    h1jer->Fit(fjer2,"QRN");
    fjer2->SetRange(15,3500);
    
    fjerp->SetLineColor(kBlue-9);
    fjerp->Draw("SAME");
    fjerm->SetLineColor(kBlue);
    fjerm->Draw("SAME");
    fjer2->SetLineColor(kGreen-9);
    fjer2->Draw("SAME");
    fjer->SetLineColor(kGreen+2);
    fjer->Draw("SAME");
    
    c2->cd(ieta);
    gPad->SetLogx();

    TH1D *h2 = tdrHist(Form("h2%s_%d",co,ieta),"Data/MC",0.8,2.0);
    h2->Draw();

    double ymin2 = h2->GetMaximum();
    double ymax2 = h2->GetMinimum();
    l->SetLineStyle(kDotted);
    l->SetLineStyle(kGray+1);
    l->DrawLine(59,ymin2,59,ymax2);
    l->DrawLine(ptmax,ymin2,ptmax,ymax2);
    l->SetLineStyle(kGray);
    l->DrawLine(ptmax2,ymin2,ptmax2,ymax2);
    l->SetLineStyle(kSolid);
    l->SetLineColor(kBlack);
    l->DrawLine(15,1,3500,1);

    TH1D *h1jerr = (TH1D*)h1jer->Clone(Form("h1jerr%s_%d",co,ieta));
    TH1D *h1jerpr = (TH1D*)h1jer->Clone(Form("h1jerpr%s_%d",co,ieta));

    h1jerr->Divide(h1jerm);
    h1jerpr->Divide(h1jerp);

    h1jerr->GetXaxis()->SetRangeUser(59,ptmax2);
    h1jerpr->GetXaxis()->SetRangeUser(59,ptmax2);

    h1jerpr->SetMarkerSize(0.7);
    h1jerr->SetMarkerSize(0.7);
    
    // HISTP loses the histograms, so drwa HIST once at the back first
    tdrDraw(h1jerpr,"HIST",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(h1jerr,"HEPz",kFullCircle,kGreen+2,kSolid,-1,kNone);

    tdrDraw(h1jerpr,"HIST",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(h1jerr,"HEPz",kFullCircle,kGreen+2,kSolid,-1,kNone);

    tex->DrawLatex(0.55,0.92,Form("%1.3f<|#eta|<%1.3f",etamin,etamax));

    if (ieta==1 || ieta%6==0) {
      TLegend *leg = tdrLeg(0.55,0.90-2*0.075,0.85,0.90);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(h1jerr,"MadGraph","PLE");
      leg->AddEntry(h1jerpr,"Pythia8","PL");
    }

    /*
    TF1 *fjerpr= new TF1(Form("fjerpr%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])/"
			 "sqrt([4]*fabs([4])/(x*x)+[5]*[5]*pow(x,[7])+[6]*[6])",
			 15,3500);
    fjerpr->SetParameters(fjer2->GetParameter(0),fjer2->GetParameter(1),
			  fjer2->GetParameter(2),fjer2->GetParameter(3),
			  fjerp->GetParameter(0),fjerp->GetParameter(1),
			  fjerp->GetParameter(2),fjerp->GetParameter(3));
    fjerpr->SetLineColor(kBlue);
    fjerpr->Draw("SAME");

    TF1 *fjerr = new TF1(Form("fjerr%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])/"
			 "sqrt([4]*fabs([4])/(x*x)+[5]*[5]*pow(x,[7])+[6]*[6])",
			 15,3500);
    fjerr->SetParameters(fjer->GetParameter(0),fjer->GetParameter(1),
			 fjer->GetParameter(2),fjer->GetParameter(3),
			 fjerm->GetParameter(0),fjerm->GetParameter(1),
			 fjerm->GetParameter(2),fjerm->GetParameter(3));
    fjerr->SetLineColor(kGreen+2);
    fjerr->SetLineWidth(2);
    fjerr->Draw("SAME");
    */

    TF1 *fjerpr= new TF1(Form("fjerpr%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
			 "[1]*[1]*pow(x,[3])+[2]*[2])/"
			 "sqrt([5]*fabs([5])/(x*x)+[9]*[9]/(x*x)+"
			 "[6]*[6]*pow(x,[8])+[7]*[7])",
			 15,3500);
    fjerpr->SetParameters(fjer2->GetParameter(0),fjer2->GetParameter(1),
			  fjer2->GetParameter(2),fjer2->GetParameter(3),
			  fjer2->GetParameter(4),
			  fjerp->GetParameter(0),fjerp->GetParameter(1),
			  fjerp->GetParameter(2),fjerp->GetParameter(3),
			  fjerp->GetParameter(4));
    fjerpr->SetLineColor(kBlue);
    fjerpr->Draw("SAME");

    TF1 *fjerr = new TF1(Form("fjerr%d",ieta),
			 "sqrt([0]*fabs([0])/(x*x)+[4]*[4]/(x*x)+"
			 "[1]*[1]*pow(x,[3])+[2]*[2])/"
			 "sqrt([5]*fabs([5])/(x*x)+[9]*[9]/(x*x)+"
			 "[6]*[6]*pow(x,[8])+[7]*[7])",
			 15,3500);
    fjerr->SetParameters(fjer->GetParameter(0),fjer->GetParameter(1),
			 fjer->GetParameter(2),fjer->GetParameter(3),
			 fjer->GetParameter(4),
			 fjerm->GetParameter(0),fjerm->GetParameter(1),
			 fjerm->GetParameter(2),fjerm->GetParameter(3),
			 fjerm->GetParameter(4));
    fjerr->SetLineColor(kGreen+2);
    fjerr->SetLineWidth(2);
    fjerr->Draw("SAME");
    
    for (int ipt = 1; ipt != h2jerf->GetNbinsY()+1; ++ipt) {
      double etamin = h2jerf->GetXaxis()->GetBinLowEdge(ieta);
      double etamax = h2jerf->GetXaxis()->GetBinLowEdge(ieta+1);
      double absetamin = min(fabs(etamin),fabs(etamax));
      double pt = h2jerf->GetYaxis()->GetBinLowEdge(ipt);
      //if (pt<=240 || (absetamin<3.2 && pt<=480) ||
      //  (absetamin<2.5 && pt<=1000) || (absetamin<1.5 && pt<=2000)) {
      if (cosh(absetamin)*pt < 6500.) {
	h2jerf->SetBinContent(ieta, ipt, fjer->Eval(pt));
	h2jersf->SetBinContent(ieta, ipt, fjerr->Eval(pt));
      }
      
      //int ietam = h2jerf->FindBin(-h2jerf->GetBinCenter(ieta));
      //h2jerf->SetBinContent(ietam, ipt, fjer->Eval(pt));
      //h2jersf->SetBinContent(ietam, ipt, fjerr->Eval(pt));
    }
  } // for ieta

  c1->SaveAs(Form("pdf/DijetHistosOverlayJER_%s.pdf","2016GH"));
  c2->SaveAs(Form("pdf/DijetHistosOverlayJER_%s_ratio.pdf","2016GH"));

  if (true) { // JER and JERSF vs eta

    int color[ny] = {kBlack, kMagenta+1, kBlue, kCyan+1, kGreen+2,
		     kYellow+2, kOrange+1, kRed, kBlack};
      
    TH1D *h = tdrHist("h","JER in data",0,0.55,"|#eta|",0.0,5.191);
    TH1D *hd = tdrHist("hd","Scale factor",0.9,1.6,"|#eta|",0.0,5.191);
    lumi_13TeV = "RunGH, 16.8 fb^{-1}";
    TCanvas *c3 = tdrDiCanvas("c3",h,hd,4,11);

    c3->cd(1);

    l->SetLineStyle(kDotted);
    l->DrawLine(1.3,0,1.3,0.35);
    l->DrawLine(2.5,0,2.5,0.40);
    l->DrawLine(3.139,0,3.139,0.50);
    
    TLegend *leg = tdrLeg(0.70,0.89-ny*0.035,0.95,0.89);
    leg->SetTextSize(0.04);

    for (int ipt = 1; ipt != h2jerf->GetNbinsY()+1; ++ipt) {
      TH1D *h1jerf = h2jerf->ProjectionX(Form("h1jerf_%d",ipt),ipt,ipt);
      tdrDraw(h1jerf,"HIST][",kNone,color[ipt-1],kSolid,-1,kNone);

      int pt = h2jerf->GetYaxis()->GetBinLowEdge(ipt);
      leg->AddEntry(h1jerf,Form("%d GeV",pt),"L");
    }
    
    c3->cd(2);
    
    l->SetLineStyle(kDashed);
    l->DrawLine(0,1,5.191,1);
    l->SetLineStyle(kDotted);
    l->DrawLine(1.3,0.9,1.3,1.2);
    l->DrawLine(2.5,0.9,2.5,1.45);
    l->DrawLine(3.139,0.9,3.139,1.55);
    
    for (int ipt = 1; ipt != h2jersf->GetNbinsY()+1; ++ipt) {
      TH1D *h1jersf = h2jersf->ProjectionX(Form("h1jersf_%d",ipt),ipt,ipt);
      tdrDraw(h1jersf,"HIST][",kNone,color[ipt-1],kSolid,-1,kNone);
    }

    c3->SaveAs(Form("pdf/DijetHistosOverlayJER_%s_JERvsEta.pdf","2016GH"));
  } // JER+JERSF
			      
} // DijetHistosOverLayJER

int findBin(TH2D *h2, double x, double *xnew) {

    int i = h2->GetXaxis()->FindBin(x);
    double x2 = h2->GetXaxis()->GetBinLowEdge(i);
    if (fabs(h2->GetXaxis()->GetBinLowEdge(i-1)-x)<fabs(x2-x)) --i;
    if (fabs(h2->GetXaxis()->GetBinLowEdge(i+1)-x)<fabs(x2-x)) ++i;
    x2 = h2->GetXaxis()->GetBinLowEdge(i);
    if (xnew) (*xnew) = x2;
    
    return i;
} // findBin


void addBins(TH1D *h1to, TH2D *h2from, int i1, int i2,
	     TH1D *h1count, TH2D *h2count) {

  assert(h1to->GetNbinsX()==h2from->GetNbinsY());
  
  for (int j = 1; j != h1to->GetNbinsX()+1; ++j) {

    double sumw = h1count->GetBinContent(j);
    double sumwz = sumw * h1to->GetBinContent(j);
    double sumwe2 = sumw * pow(h1to->GetBinError(j),2);
    
    for (int i = i1; i != i2+1; ++i) {
      
      double w = h2count->GetBinContent(i,j);
      double wz = w * h2from->GetBinContent(i,j);
      double we2 = w * pow(h2from->GetBinError(i,j),2);
      sumw += w;
      sumwz += wz;
      sumwe2 += we2;
    } // for j
    
    h1to->SetBinContent(j, sumw ? sumwz / sumw : 0);
    h1to->SetBinError(j, sumw ? sqrt(sumwe2 / sumw) : 0);
    h1count->SetBinContent(j, sumw);
  } // for i
} // addBins

void rebin(TH1D *h1, TH1D *h1n, TH1D *h1old, TH1D *h1nold) {

  // Add bins with sum of weights
  for (int i = 1; i != h1old->GetNbinsX()+1; ++i) {
    int j = h1->FindBin(h1old->GetBinCenter(i));
    double n = h1nold->GetBinContent(i);
    h1->SetBinContent(j, h1->GetBinContent(j) + n*h1old->GetBinContent(i));
    h1->SetBinError(j, h1->GetBinError(j) + n*pow(h1old->GetBinError(i),2));
    h1n->SetBinContent(j, h1n->GetBinContent(j) + n);
    h1n->SetBinError(j, h1n->GetBinContent(j)+n*pow(h1nold->GetBinError(i),2));
  }
  for (int i = 1; i != h1->GetNbinsX()+1; ++i) {
    double n = h1n->GetBinContent(i);
    h1->SetBinContent(i, n ? h1->GetBinContent(i) / n : 0);
    h1->SetBinError(i, n ? sqrt(h1->GetBinError(i)/n) : 0);
    h1n->SetBinError(i, n ? sqrt(h1n->GetBinError(i)/n) : 0);
  }

} // rebin
