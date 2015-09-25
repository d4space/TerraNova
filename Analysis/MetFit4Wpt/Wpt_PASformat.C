//================================================
//
// Perform fit to extract W->munu signal
//
//  * outputs plots and fit results summary
// from ksung fitWm.C
//________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "Math/LorentzVector.h"           // 4-vector class

#include "../Utils/MyTools.hh"	          // various helper functions
#include "../Utils/CPlot.hh"	          // helper class for plots
#include "../Utils/MitStyleRemix.hh"      // style settings for drawing
#include "../Utils/WModels.hh"            // definitions of PDFs for fitting
#include "../Utils/RecoilCorrector.hh"    // class to handle recoil corrections for MET
#endif
//From WMuNeu.h
#include "../Utils/const.h"

//=== MAIN MACRO ================================================ 

void Wpt_PASformat(const TString  outputDir,   // output directory
           const TString  filetype,		// Select input root files for Nominal, Up, Down and Before Recoil Correction
           const Double_t lumi,        // integrated luminosity (/fb)
	   const Double_t nsigma=0     // vary MET corrections by n-sigmas (nsigma=0 means nominal correction)
)
{
  gBenchmark->Start("Wpt_PASformat");

  //==================   
  // Settings 
  //==================   
  
  const TString format("png"); 
  //const TString format("pdf"); 

  // input ntuple file names
  TString fname = "MuonFitResultsModRayleighSimultNominal/RstMuon/SigYields_Nominal.root";
  if (filetype == "Electron")
    fname = "ElectronFitResultsModRayleighSimultNominal/RstElectron/SigYields_Nominal.root";

  TFile* infile = new TFile(fname);

  //-------------------------
  // Main analysis code 
  //=========================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
  //
  // Declare MET histograms
  // 0 is for total

  TH1D *hdataWPpt;
  TH1D *hdataWMpt;
  TH1D *hSigWPpt;
  TH1D *hSigWMpt;
  TH1D *hQCDWPpt;
  TH1D *hQCDWMpt;
  TH1D *hDYToMuMuP;
  TH1D *hDYToMuMuM;
  TH1D *hWToTauNuP;
  TH1D *hWToTauNuM;
  TH1D *hTTJetsP;
  TH1D *hTTJetsM;
  TH1D *hDYToTauTauP;
  TH1D *hDYToTauTauM;
  
  char histName[30],histName_org[30];

  hdataWPpt    = (TH1D*)infile->Get("hdataWPpt")    -> Clone("hdataWPpt");
  hdataWMpt    = (TH1D*)infile->Get("hdataWMpt")    -> Clone("hdataWMpt");
  hSigWPpt     = (TH1D*)infile->Get("hSigWPpt")     -> Clone("hSigWPpt");
  hSigWMpt     = (TH1D*)infile->Get("hSigWMpt")     -> Clone("hSigWMpt");
  hQCDWPpt     = (TH1D*)infile->Get("hQCDWPpt")     -> Clone("hQCDWPpt");
  hQCDWMpt     = (TH1D*)infile->Get("hQCDWMpt")     -> Clone("hQCDWMpt");
  hDYToMuMuP   = (TH1D*)infile->Get("hDYToMuMuP")   -> Clone("hDYToMuMuP");
  hDYToMuMuM   = (TH1D*)infile->Get("hDYToMuMuM")   -> Clone("hDYToMuMuM");
  hWToTauNuP   = (TH1D*)infile->Get("hWToTauNuP")   -> Clone("hWToTauNuP");
  hWToTauNuM   = (TH1D*)infile->Get("hWToTauNuM")   -> Clone("hWToTauNuM");
  hTTJetsP     = (TH1D*)infile->Get("hTTJetsP")     -> Clone("hTTJetsP");
  hTTJetsM     = (TH1D*)infile->Get("hTTJetsM")     -> Clone("hTTJetsM");
  hDYToTauTauP = (TH1D*)infile->Get("hDYToTauTauP") -> Clone("hDYToTauTauP");
  hDYToTauTauM = (TH1D*)infile->Get("hDYToTauTauM") -> Clone("hDYToTauTauM");
  
  TCanvas *c;
  c = MakeCanvas("c","c",800,800);
  c->SetPad(0,0,1.0,1.0);
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.15);
  c->SetLeftMargin(0.15);  
  c->SetRightMargin(0.07);  
  c->SetTickx(1);
  c->SetTicky(1);  
  gStyle->SetTitleOffset(1.400,"Y");
  TGaxis::SetMaxDigits(3);
  char ylabel[100];  // string buffer for y-axis label
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1} (8 TeV)",lumi*1000.);
  else         sprintf(lumitext,"%.2f fb^{-1}  at  #sqrt{s} = 8 TeV",lumi);
  char CMStext[100];
  sprintf(CMStext,"#font[61]{CMS}");
  char Preliminarytext[100];
  sprintf(Preliminarytext,"#font[52]{Preliminary}");
  char binlabel[50];
  
  // plot colors
  Int_t linecolorW   = kOrange-3;
  Int_t fillcolorW   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorQCD = kViolet+2;
  Int_t fillcolorQCD = kViolet-5;

  //
  // Dummy histograms for TLegend
  // (I can't figure out how to properly pass RooFit objects...)
  //
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  
  TH1D *hDummyW = new TH1D("hDummyW","",0,0,10);
  hDummyW->SetLineColor(linecolorW);
  hDummyW->SetFillColor(fillcolorW);
  hDummyW->SetFillStyle(1001);
  
  TH1D *hDummyEWK = new TH1D("hDummyEWK","",0,0,10);
  hDummyEWK->SetLineColor(linecolorEWK);
  hDummyEWK->SetFillColor(fillcolorEWK);
  hDummyEWK->SetFillStyle(1001);
  
  TH1D *hDummyQCD = new TH1D("hDummyQCD","",0,0,10);
  hDummyQCD->SetLineColor(linecolorQCD);
  hDummyQCD->SetFillColor(fillcolorQCD);
  hDummyQCD->SetFillStyle(1001);

// Wpt distribution=========================
  CPlot* plotWPptLog;
  CPlot* plotWMptLog;

//W plus pt distribution
  TH1D* hWptMC_p = (TH1D*)hDYToTauTauP->Clone("hWptMC_p");
  hWptMC_p->Add(hTTJetsP);
  hWptMC_p->Add(hWToTauNuP);
  hWptMC_p->Add(hDYToMuMuP);
  hWptMC_p->Add(hQCDWPpt);
  hWptMC_p->Add(hSigWPpt);
  
  double WptBinsLog[14]={1,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};
  double x1,x2,x3,x4,x5,x6,x7,err;
  TH1D *hDYToTauTauLogP = new TH1D("hDYToTauTauLogP","",13,WptBinsLog);
  TH1D *hTTJetsLogP     = new TH1D("hTTJetsLogP","",13,WptBinsLog);
  TH1D *hWToTauNuLogP   = new TH1D("hWToTauNuLogP","",13,WptBinsLog);
  TH1D *hDYToMuMuLogP   = new TH1D("hDYToMuMuLogP","",13,WptBinsLog);
  TH1D *hQCDWptLogP     = new TH1D("hQCDWptLogP","",13,WptBinsLog);
  TH1D *hSigWptLogP     = new TH1D("hSigWptLogP","",13,WptBinsLog);
  TH1D *hdataWptLogP    = new TH1D("hdataWptLogP","",13,WptBinsLog);
  TAxis *xaxis = hSigWPpt->GetXaxis();
  for (int i=1; i<=hSigWPpt->GetNbinsX(); i++){
    x1=hDYToTauTauP->GetBinContent(i);
    x2=hTTJetsP->GetBinContent(i);
    x3=hWToTauNuP->GetBinContent(i);
    x4=hDYToMuMuP->GetBinContent(i);
    x5=hQCDWPpt->GetBinContent(i);
    x6=hSigWPpt->GetBinContent(i);
    x7=hdataWPpt->GetBinContent(i);
    err = hdataWPpt->GetBinError(i);
    hDYToTauTauLogP->Fill(xaxis->GetBinCenter(i),x1);
    hTTJetsLogP->Fill(xaxis->GetBinCenter(i),x2);
    hWToTauNuLogP->Fill(xaxis->GetBinCenter(i),x3);
    hDYToMuMuLogP->Fill(xaxis->GetBinCenter(i),x4);
    hQCDWptLogP->Fill(xaxis->GetBinCenter(i),x5);
    hSigWptLogP->Fill(xaxis->GetBinCenter(i),x6);
    hdataWptLogP->Fill(xaxis->GetBinCenter(i),x7);
    hdataWptLogP->SetBinError(i,err);
  }
  
  TH1D* hEwkWPptLog = (TH1D*)hDYToTauTauLogP->Clone("hEwkWPptLog");
  hEwkWPptLog->Add(hTTJetsLogP);
  hEwkWPptLog->Add(hWToTauNuLogP);
  hEwkWPptLog->Add(hDYToMuMuLogP);
  
  TString plotName = "FitWDistribution_MuonPLog";
  if (filetype == "Electron")
    plotName = "FitWDistribution_ElePLog";
  plotWPptLog=new CPlot(plotName,"","p_{T}^{W} [GeV]","Events");
  plotWPptLog->setOutDir(CPlot::sOutDir);
  plotWPptLog->AddToStack(hEwkWPptLog,"#font[42]{EWK+t#bar{t}}",fillcolorEWK,linecolorEWK);
  plotWPptLog->AddToStack(hQCDWptLogP,"#font[42]{QCD}",fillcolorQCD,linecolorQCD);
  if(filetype == "Muon")
    plotWPptLog->AddToStack(hSigWptLogP,"#font[42]{W^{+}#rightarrow #mu^{+}#nu}",fillcolorW,linecolorW);
  if(filetype == "Electron")
    plotWPptLog->AddToStack(hSigWptLogP,"#font[42]{W^{+}#rightarrow e^{+}#nu}",fillcolorW,linecolorW);
  plotWPptLog->AddHist1D(hdataWptLogP,"#font[42]{data}","E");
  plotWPptLog->SetLegend(0.75,0.78,.98,0.88);
  plotWPptLog->SetYRange(0.25,1.4*(hWptMC_p->GetMaximum()));
  plotWPptLog->SetLogx();
  plotWPptLog->SetLogy();
  plotWPptLog->AddTextBox(CMStext,0.14,0.90,0.24,0.97,0);
  plotWPptLog->AddTextBox(lumitext,0.65,0.90,0.95,0.97,0);
  plotWPptLog->Draw(c,kFALSE,format);
  gPad->RedrawAxis();
  plotName = "Wpt_plots/FitWDistribution_MuonPLog.pdf";
  if (filetype == "Electron")
    plotName = "Wpt_plots/FitWDistribution_ElePLog.pdf";
  c->SaveAs(plotName);

//W minus pt distribution
  TH1D* hWptMC_m = (TH1D*)hDYToTauTauM->Clone("hWptMC_m");
  hWptMC_m->Add(hTTJetsM);
  hWptMC_m->Add(hWToTauNuM);
  hWptMC_m->Add(hDYToMuMuM);
  hWptMC_m->Add(hQCDWMpt);
  hWptMC_m->Add(hSigWMpt);
  
  TH1D *hDYToTauTauLogM = new TH1D("hDYToTauTauLogM","",13,WptBinsLog);
  TH1D *hTTJetsLogM = new TH1D("hTTJetsLogM","",13,WptBinsLog);
  TH1D *hWToTauNuLogM = new TH1D("hWToTauNuLogM","",13,WptBinsLog);
  TH1D *hDYToMuMuLogM = new TH1D("hDYToMuMuLogM","",13,WptBinsLog);
  TH1D *hQCDWptLogM = new TH1D("hQCDWptLogM","",13,WptBinsLog);
  TH1D *hSigWptLogM = new TH1D("hSigWptLogM","",13,WptBinsLog);
  TH1D *hdataWptLogM = new TH1D("hdataWptLogM","",13,WptBinsLog);
  for (int i=1; i<=hSigWMpt->GetNbinsX(); i++){
    x1=hDYToTauTauM->GetBinContent(i);
    x2=hTTJetsM->GetBinContent(i);
    x3=hWToTauNuM->GetBinContent(i);
    x4=hDYToMuMuM->GetBinContent(i);
    x5=hQCDWMpt->GetBinContent(i);
    x6=hSigWMpt->GetBinContent(i);
    x7=hdataWMpt->GetBinContent(i);
    err = hdataWMpt->GetBinError(i);
    hDYToTauTauLogM->Fill(xaxis->GetBinCenter(i),x1);
    hTTJetsLogM->Fill(xaxis->GetBinCenter(i),x2);
    hWToTauNuLogM->Fill(xaxis->GetBinCenter(i),x3);
    hDYToMuMuLogM->Fill(xaxis->GetBinCenter(i),x4);
    hQCDWptLogM->Fill(xaxis->GetBinCenter(i),x5);
    hSigWptLogM->Fill(xaxis->GetBinCenter(i),x6);
    hdataWptLogM->Fill(xaxis->GetBinCenter(i),x7);
    hdataWptLogM->SetBinError(i,err);
  }
  
  TH1D* hEwkWMptLog = (TH1D*)hDYToTauTauLogM->Clone("hEwkWMptLog");
  hEwkWMptLog->Add(hTTJetsLogM);
  hEwkWMptLog->Add(hWToTauNuLogM);
  hEwkWMptLog->Add(hDYToMuMuLogM);
  
  plotName = "FitWDistribution_MuonMLog";
  if (filetype == "Electron")
    plotName = "FitWDistribution_EleMLog";
  plotWMptLog=new CPlot(plotName,"","p_{T}^{W} [GeV]","Events");
  plotWMptLog->setOutDir(CPlot::sOutDir);
  plotWMptLog->AddToStack(hEwkWMptLog,"#font[42]{EWK+t#bar{t}}",fillcolorEWK,linecolorEWK);
  plotWMptLog->AddToStack(hQCDWptLogM,"#font[42]{QCD}",fillcolorQCD,linecolorQCD);
  if(filetype == "Muon")
    plotWMptLog->AddToStack(hSigWptLogM,"#font[42]{W^{-}#rightarrow #mu^{-} #bar{#nu}}",fillcolorW,linecolorW);
  if(filetype == "Electron")
    plotWMptLog->AddToStack(hSigWptLogM,"#font[42]{W^{-}#rightarrow e^{-} #bar{#nu}}",fillcolorW,linecolorW);
  plotWMptLog->AddHist1D(hdataWptLogM,"#font[42]{data}","E");
  plotWMptLog->SetLegend(0.75,0.78,.98,0.88);
  plotWMptLog->SetYRange(0.25,1.4*(hWptMC_m->GetMaximum()));
  plotWMptLog->SetLogx();
  plotWMptLog->SetLogy();
  plotWMptLog->AddTextBox(CMStext,0.14,0.90,0.24,0.97,0);
  plotWMptLog->AddTextBox(lumitext,0.65,0.90,0.95,0.97,0);
  plotWMptLog->Draw(c,kFALSE,format);
  gPad->RedrawAxis();
  plotName = "Wpt_plots/FitWDistribution_MuonMLog.pdf";
  if (filetype == "Electron")
    plotName = "Wpt_plots/FitWDistribution_EleMLog.pdf";
  c->SaveAs(plotName);

  gBenchmark->Show("Wpt_PASformat");
}
