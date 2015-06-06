#include "../Utils/MyTools.hh"
#include "../Utils/CPlot.hh"
#include "../Utils/MitStyleRemix.hh"
#include "../Utils/const.h"
#include <TFile.h>

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooChi2Var.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"

void Wpt_muPtComp(const TString outputDir)
{
  TH1D* makeDiffHist(TH1D* h1, TH1D* h2, const TString name);
  const TString format("png"); 
  Int_t ratioColor = kGray+2;

  TFile *fname_MC;
  TFile *fname_RD;

  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;

  fname_MC = new TFile("Muon2012LoPU/Muon_DYToMuMu_S8.root");
  fname_RD = new TFile("Muon2012LoPU/Muon_RD_LowPU.root");

  CPlot *plotLep1Pt;
  CPlot *plotLep2Pt;
  CPlot *plotLepPt_p;
  CPlot *plotLepPt_m;

  CPlot *DiffPlotLep1Pt;
  CPlot *DiffPlotLep2Pt;
  CPlot *DiffPlotLepPt_p;
  CPlot *DiffPlotLepPt_m;
  
  TH1D *hDiffLep1;
  TH1D *hDiffLep2;
  TH1D *hDiffLepP;
  TH1D *hDiffLepM;

  TH1D *hMCLep1Pt;
  TH1D *hMCLep2Pt;
  TH1D *hMCLepPt_p;
  TH1D *hMCLepPt_m;

  TH1D *hRDLep1Pt;
  TH1D *hRDLep2Pt;
  TH1D *hRDLepPt_p;
  TH1D *hRDLepPt_m;

  TCanvas *c;
  c = MakeCanvas("c","c",800,600);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);

  hMCLep1Pt =(TH1D*)fname_MC->Get("h1_ZLep1Pt") ->Clone("hMCLep1Pt");  hMCLep1Pt->Sumw2();
  hMCLep2Pt =(TH1D*)fname_MC->Get("h1_ZLep2Pt") ->Clone("hMCLep2Pt");  hMCLep2Pt->Sumw2();
  hMCLepPt_p=(TH1D*)fname_MC->Get("h1_ZLepPt_p")->Clone("hMCLepPt_p"); hMCLepPt_p->Sumw2();
  hMCLepPt_m=(TH1D*)fname_MC->Get("h1_ZLepPt_m")->Clone("hMCLepPt_m"); hMCLepPt_m->Sumw2();
  
  hRDLep1Pt =(TH1D*)fname_RD->Get("h1_ZLep1Pt") ->Clone("hRDLep1Pt");  hRDLep1Pt->Sumw2();
  hRDLep2Pt =(TH1D*)fname_RD->Get("h1_ZLep2Pt") ->Clone("hRDLep2Pt");  hRDLep2Pt->Sumw2();
  hRDLepPt_p=(TH1D*)fname_RD->Get("h1_ZLepPt_p")->Clone("hRDLepPt_p"); hRDLepPt_p->Sumw2();
  hRDLepPt_m=(TH1D*)fname_RD->Get("h1_ZLepPt_m")->Clone("hRDLepPt_m"); hRDLepPt_m->Sumw2();
  
  //============================
  //Draw first Lepton pT
  //============================
  plotLep1Pt= new CPlot("Lep1Pt","","","Events");
  plotLep1Pt->setOutDir(CPlot::sOutDir);

  plotLep1Pt->AddHist1D(hMCLep1Pt,"Z #rightarrow #mu#mu","HIST",kRed);
  plotLep1Pt->AddHist1D(hRDLep1Pt,"Data","HIST",kBlack);
  plotLep1Pt->SetYRange(0.,1.2*(hRDLep1Pt->GetMaximum()));
  plotLep1Pt->SetLegend(0.68,0.57,0.9,0.8);
  plotLep1Pt->Draw(c,kFALSE,format,1);
  
  hDiffLep1 = makeDiffHist(hMCLep1Pt,hRDLep1Pt,"hDiffLep1");
  hDiffLep1->SetMarkerStyle(kFullCircle);
  hDiffLep1->SetMarkerSize(0.1);

  DiffPlotLep1Pt=new CPlot("Lep1Pt","","1^{st} muon p_{T} [GeV]","#chi");
  DiffPlotLep1Pt->setOutDir(CPlot::sOutDir);
  DiffPlotLep1Pt->AddHist1D(hDiffLep1,"EX0",ratioColor);
  DiffPlotLep1Pt->SetYRange(-8,8);
  DiffPlotLep1Pt->AddLine(0, 0,150, 0,kBlack,1);
  DiffPlotLep1Pt->AddLine(0, 5,150, 5,kBlack,3);
  DiffPlotLep1Pt->AddLine(0,-5,150,-5,kBlack,3);
  DiffPlotLep1Pt->Draw(c,kTRUE,format,2);
  //============================
  
  //============================
  //Draw second Lepton pT
  //============================
  plotLep2Pt= new CPlot("Lep2Pt","","","Events");
  plotLep2Pt->setOutDir(CPlot::sOutDir);

  plotLep2Pt->AddHist1D(hMCLep2Pt,"Z #rightarrow #mu#mu","HIST",kRed);
  plotLep2Pt->AddHist1D(hRDLep2Pt,"Data","HIST",kBlack);
  plotLep2Pt->SetYRange(0.,1.2*(hRDLep2Pt->GetMaximum()));
  plotLep2Pt->SetLegend(0.68,0.57,0.9,0.8);
  plotLep2Pt->Draw(c,kFALSE,format,1);
  
  hDiffLep2 = makeDiffHist(hMCLep2Pt,hRDLep2Pt,"hDiffLep2");
  hDiffLep2->SetMarkerStyle(kFullCircle);
  hDiffLep2->SetMarkerSize(0.1);

  DiffPlotLep2Pt=new CPlot("Lep2Pt","","2^{nd} muon p_{T} [GeV]","#chi");
  DiffPlotLep2Pt->setOutDir(CPlot::sOutDir);
  DiffPlotLep2Pt->AddHist1D(hDiffLep2,"EX0",ratioColor);
  DiffPlotLep2Pt->SetYRange(-8,8);
  DiffPlotLep2Pt->AddLine(0, 0,150, 0,kBlack,1);
  DiffPlotLep2Pt->AddLine(0, 5,150, 5,kBlack,3);
  DiffPlotLep2Pt->AddLine(0,-5,150,-5,kBlack,3);
  DiffPlotLep2Pt->Draw(c,kTRUE,format,2);
  //============================
  
  //============================
  //Draw Lepton plus pT
  //============================
  plotLepPt_p= new CPlot("LepPt_p","","","Events");
  plotLepPt_p->setOutDir(CPlot::sOutDir);

  plotLepPt_p->AddHist1D(hMCLepPt_p,"Z #rightarrow #mu#mu","HIST",kRed);
  plotLepPt_p->AddHist1D(hRDLepPt_p,"Data","HIST",kBlack);
  plotLepPt_p->SetYRange(0.,1.2*(hRDLepPt_p->GetMaximum()));
  plotLepPt_p->SetLegend(0.68,0.57,0.9,0.8);
  plotLepPt_p->Draw(c,kFALSE,format,1);
  
  hDiffLepP = makeDiffHist(hMCLepPt_p,hRDLepPt_p,"hDiffLepP");
  hDiffLepP->SetMarkerStyle(kFullCircle);
  hDiffLepP->SetMarkerSize(0.1);

  DiffPlotLepPt_p=new CPlot("LepPt_p","","#mu^{+} p_{T} [GeV]","#chi");
  DiffPlotLepPt_p->setOutDir(CPlot::sOutDir);
  DiffPlotLepPt_p->AddHist1D(hDiffLepP,"EX0",ratioColor);
  DiffPlotLepPt_p->SetYRange(-8,8);
  DiffPlotLepPt_p->AddLine(0, 0,150, 0,kBlack,1);
  DiffPlotLepPt_p->AddLine(0, 5,150, 5,kBlack,3);
  DiffPlotLepPt_p->AddLine(0,-5,150,-5,kBlack,3);
  DiffPlotLepPt_p->Draw(c,kTRUE,format,2);
  //============================

  //============================
  //Draw Lepton minus pT
  //============================
  plotLepPt_m= new CPlot("LepPt_m","","","Events");
  plotLepPt_m->setOutDir(CPlot::sOutDir);

  plotLepPt_m->AddHist1D(hMCLepPt_m,"Z #rightarrow #mu#mu","HIST",kRed);
  plotLepPt_m->AddHist1D(hRDLepPt_m,"Data","HIST",kBlack);
  plotLepPt_m->SetYRange(0.,1.2*(hRDLepPt_m->GetMaximum()));
  plotLepPt_m->SetLegend(0.68,0.57,0.9,0.8);
  plotLepPt_m->Draw(c,kFALSE,format,1);
  
  hDiffLepM = makeDiffHist(hMCLepPt_m,hRDLepPt_m,"hDiffLepM");
  hDiffLepM->SetMarkerStyle(kFullCircle);
  hDiffLepM->SetMarkerSize(0.1);

  DiffPlotLepPt_m=new CPlot("LepPt_m","","#mu^{-} p_{T} [GeV]","#chi");
  DiffPlotLepPt_m->setOutDir(CPlot::sOutDir);
  DiffPlotLepPt_m->AddHist1D(hDiffLepM,"EX0",ratioColor);
  DiffPlotLepPt_m->SetYRange(-8,8);
  DiffPlotLepPt_m->AddLine(0, 0,150, 0,kBlack,1);
  DiffPlotLepPt_m->AddLine(0, 5,150, 5,kBlack,3);
  DiffPlotLepPt_m->AddLine(0,-5,150,-5,kBlack,3);
  DiffPlotLepPt_m->Draw(c,kTRUE,format,2);
  //============================
}

TH1D *makeDiffHist(TH1D* h1, TH1D* h2, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=h1->GetNbinsX(); ibin++) {

    Double_t diff = (h1->GetBinContent(ibin)-h2->GetBinContent(ibin));

    Double_t err = sqrt(h1->GetBinContent(ibin));
    if(err==0) err= sqrt(h2->GetBinContent(ibin));

    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
  }
  hDiff->GetYaxis()->SetTitleOffset(0.42);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  return hDiff;
}
