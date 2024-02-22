// Purpose: Draw Multijet results to do sanity checks
//          Compare various pT bins and MPF, DB (HDM), data vs MC
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include <map>
#include <string>

#include "../tdrstyle_mod22.C"

string version = "v37_Summer23MG_NoL2L3Res_ReWPU_OnJet";

void drawMultijets(string epoch="2022E", string version="v35a");

void drawMultijet() {
    drawMultijets("Summer23MG", version);
}

void drawMultijets(string epoch, string version) {
    setTDRStyle();
    TDirectory *curdir = gDirectory;

    // Load requested MC file
    const char *ce = epoch.c_str();
    const char *cv = version.c_str();
    TFile *fm = new TFile("/media/storage/nestorma/dijet/rootfiles/v37_Summer23MG_NoL2L3Res_ReWPU_OnJet/jmenano_mc_cmb_Summer23MG_new_v37_Summer23MG_NoL2L3Res_ReWPU_OnJet.root");
    assert(fm && !fm->IsZombie());

    // Print the MC file name
    cout << "MC: " << fm->GetName() << endl;

    // List results to be plotted
    string vd[] = {"pm0l","pm0a","pm0r", "pm2l","pm2a","pm2r"};
    const int nvd = sizeof(vd)/sizeof(vd[0]);

    // Set plotting styles (color, marker, line style) for results
    map<string,int> color;
    color["pm0l"] = kRed;
    color["pm2l"] = kRed;
    color["pm0a"] = kGreen+2;
    color["pm2a"] = kGreen+2;
    color["pm0r"] = kBlue;
    color["pm2r"] = kBlue;
    map<string,int> marker;
    marker["pm0l"] = kFullTriangleDown;
    marker["pm0a"] = kFullDiamond;
    marker["pm0r"] = kFullTriangleUp;
    marker["pm2l"] = kOpenTriangleDown;
    marker["pm2a"] = kOpenDiamond;
    marker["pm2r"] = kOpenTriangleUp;
    map<string,int> style;
    style["pm0l"] = kSolid;
    style["pm0a"] = kSolid;
    style["pm0r"] = kSolid;
    style["pm2l"] = kDotted;
    style["pm2a"] = kDotted;
    style["pm2r"] = kDotted;
    map<string,const char*> label;
    label["pm0l"] = "MPF lead";
    label["pm2l"] = "DB lead";
    label["pm0a"] = "MPF avg.";
    label["pm2a"] = "DB avg.";
    label["pm0r"] = "MPF recoil";
    label["pm2r"] = "DB recoil";

    // Create canvas for plots
    double ptmin = 114;
    double ptmax = 2100;
    TH1D *hu = tdrHist("hu","Multijet #LTp_{T,lead}#GT / #LTp_{T,recoil}#GT",
                       0.85,1.30,"#LTp_{T,lead}#GT (GeV)",ptmin,ptmax);
    TLine *l = new TLine();
    map<string,const char*> title;
    title["2022C"] = "2022C";
    title["2022CD"] = "2022CD";
    title["2022E"] = "2022E";
    title["2022FG"] = "2022FG";
    title["Summer23MG"] = "Summer23MG";
    title["2023Cv4"] = "2023Cv4";
    title["Summer23MGBPix"] = "Summer23MGBPix";
    lumi_136TeV = title[ce];
    TCanvas *c1 = tdrDiCanvas("c1",hu,hu,8,11);//,kSquare);
    c1->cd(1);
    gPad->SetLogx();
    l->DrawLine(ptmin,1,ptmax,1);

    TLegend *legd = tdrLeg(0.60,0.90-0.05*nvd,0.90,0.90);
    TLegend *legm = tdrLeg(0.60,0.88-0.05*nvd,0.90,0.88);

    c1->cd(2);
    gPad->SetLogx();
    l->DrawLine(ptmin,1,ptmax,1);

    // Plot results
    for (int i = 0; i != nvd; ++i) {
        const char *ch = vd[i].c_str();

        cout << "Drawing " << ch << endl;

        c1->cd(1);
        TProfile *pm = (TProfile*)fm->Get(Form("Multijet/%s",ch)); assert(pm);
        pm->GetXaxis()->SetRangeUser(ptmin,ptmax);
        cout << "Drawing first profile" << endl;
        tdrDraw(pm,"HIST",marker[ch],color[ch],style[ch],-1,kNone);

        cout << "Drawing ratio plot" << endl;
        c1->cd(2);
        TH1D *hm = (TH1D*)pm->ProjectionX("hm");
        tdrDraw(hm,"Pz",marker[ch],color[ch],kSolid);
    }

    c1->SaveAs(Form("/media/storage/nestorma/dijet/rootfiles/pdf/v37_Summer23MG_NoL2L3Res_ReWPU_OnJet/drawMultijet/drawMultijet_%s_%s.pdf",ce,cv));
} // drawMultijet

