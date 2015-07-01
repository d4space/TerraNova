#include <set>
#include <time.h>

#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<THStack.h>
#include<TGraphErrors.h>
#include<TGraphAsymmErrors.h>
#include<TCanvas.h>
#include<TFrame.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TChain.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>
using namespace std;
#include "../Utils/tdrstyle.C"
#include "../Utils/CMS_lumi.C"
#include "selections.h"
#include "_unfoldData.h"
#include "systematics.h"

TH1F* ConvertGraphToHisto(TGraph *pGraph);

const TString format("pdf");
//const TString format("png");

int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
int iPos = 0;

//
// Simple example of macro: plot with CMS name and lumi text
//  (this script does not pretend to work in all configurations)
// iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
// For instance:
//               iPeriod = 3 means: 7 TeV + 8 TeV
//               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV
// Initiated by: Gautier Hamel de Monchenault (Saclay)
//
int W = 800;
int H = 800;

int H_ref = 600;
int W_ref = 800;

// references for T, B, L, R
float T = 0.08*H_ref;
float B = 0.12*H_ref;
float L = 0.12*W_ref;
float R = 0.04*W_ref;

void PrintIt(TPad *pad, TString title);

void fillHistos(TTree* tree,
                TH1F *hPt,
                double mLow,
                double mHigh,
                int mode);


TGraphAsymmErrors* Powheg();
TGraphAsymmErrors* Resbos();
TGraphAsymmErrors* Madgraph();
TGraphAsymmErrors* FEWZ();


// x-section for the bkg processes (in pb)
double xsec_ttbarjets  = 225.197;
double xsec_WW         = 57.1097;
double xsec_WZ         = 32.3161;
double xsec_ZZ         = 8.25561;
double xsec_ztautau    = 5745.25/3.;
double xsec_WJetsToLNu = 30400; 
double xsec_QCD        = 134680;


// number of original MC events
// this you will get it from crab
int nMCEvTTbarJets   = 6736135;
int nMCEvZtautau     = 1737776;
int nMCEvWW          = 2519385;
int nMCEvWZ          = 2182479;
int nMCEvZZ          = 1890152;
int nMCEvWJetsToLNu  = 2053500;
int nMCEvQCD         = 7529312;

 
/* 
void DrawWithRatio(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, 
                   TGraphAsymmErrors* gDen, TGraph* band, TLine* line,
                   TLegend* leg=NULL, TLatex* tex=NULL, TLatex* tex1=NULL);

void DrawWithRatioTotal(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, 
                   TGraphAsymmErrors* gDen, 
                   TGraphAsymmErrors* gDen2, 
                   TGraphAsymmErrors* gDen3, 
		   TGraph* band, 
		   TGraph* band2, 
		   TGraph* band3, 
		   TLine* line,
		   TLine* line2,
		   TLine* line3,

		   TLegend* leg=NULL, 
		   TLegend* leg2=NULL, 
		   TLegend* leg3=NULL, 
		   
		   TLatex* tex=NULL, TLatex* tex1=NULL, TLatex* tex2=NULL);
*/
void DrawWithRatioTotal(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, 
                   TGraphAsymmErrors* gDen, 
                   TGraphAsymmErrors* gDen2, 
                   TGraphAsymmErrors* gDen3, 
                   TGraphAsymmErrors* gDen4, 
		   TGraph* band, 
		   TGraph* band2, 
		   TGraph* band3, 

		   TLegend* leg=NULL, 
		   TLegend* leg2=NULL, 
		   TLegend* leg3=NULL, 
		   
		   TLatex* tex=NULL);

void axis1F(TH1F  *histo,
            TAxis *xaxis,
            TAxis *yaxis,
            char  *xtitle,
            char  *ytitle);


//void diffXSecResAcc_resbos(double intLumi_invpb  =   18.424,
void diffXSecResAcc_Powheg_Madgraph_ResBos_FEWZ(double intLumi_invpb  =   18.424,
                    int    mBins          =   30,  // Z mass
                    double mLow           =   60,  // Z mass
                    double mHigh          =  120,  // Z mass
                    bool   isSave         = !true) {

  gROOT->LoadMacro("../Utils/tdrstyle.C");
  setTDRStyle();
  gROOT->LoadMacro("../Utils/CMS_lumi.C");
  writeExtraText = "true";
  extraText = "Preliminary";
  lumi_8TeV = "18.4 pb^{-1}";

  
  gROOT->Clear();
  gStyle->SetOptStat(0);  
 
  // ---------------------------------------------------------------------------
  // general variables
  char cTitle[200];
  //sprintf(cTitle, "CMS Preliminary");
  sprintf(cTitle, "");

  char cTitleExtd[200];
  //sprintf(cTitleExtd, "CMS Preliminary, 18.4 pb^{-1} at #sqrt{s} = 8 TeV");

  //TString png      = "png/diffXSec/";
  //TString rootPlot = "rootPlot/diffXSec/";
  // ---------------------------------------------------------------------------

  
  const int nptBins = 18;
  Double_t xbins_pt[nptBins+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};


  // this is useful only if we stick the data in center of the bins
  //    Double_t error_xbins[nptBins];
  //for (int i=0;i<nptBins;i++) error_xbins[i] = (xbins_pt[i+1] - xbins_pt[i]) / 2;
  
  Double_t error_xbinsl[nptBins];
  Double_t error_xbinsh[nptBins];
  
  int i;
  for (i=0;i<nptBins;i++) {
    
    error_xbinsl[i] = (point[i] - xbins_pt[i]);
    error_xbinsh[i] = (xbins_pt[i+1]-point[i]);
   
  }


  TH1F* hXSec          = new TH1F("hXSec"      , "", nptBins, xbins_pt);
  TH1F* hXSecFSR       = new TH1F("hXSecFSR"      , "", nptBins, xbins_pt);

  // Fill the data
  //  TFile* fileData = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_5_ecalpatch1/Data/ntuplesWITHdz/SingleMuRun2012A_190949_191090_191367_193112_193116.root");
  //TFile* fileData = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_5_ecalpatch1/Data/ntuplesWITHdz/fullselection/data2012A_LowPU.root");
  TFile* fileData = new TFile("data2012A_LowPU.root");
  TTree* treeData = (TTree*) fileData -> Get("tree");

  // Fill the TTbarJets MC  
  // TFile* fileMC_ttbarjets = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/TTJets/mergeFile/TTJets_merged.root");
  TFile* fileMC_ttbarjets = new TFile("TTJets.root");
  TTree* treeMC_ttbarjets = (TTree*) fileMC_ttbarjets -> Get("tree");
 
  // Fill the ZTauTau MC  
  //  TFile* fileMC_ztautau = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/DYToTauTau/mergeFile/DYToTauTau_merged.root");
  TFile* fileMC_ztautau = new TFile("DYToTauTau.root");
  TTree* treeMC_ztautau = (TTree*) fileMC_ztautau -> Get("tree");
 
  // Fill the WZ MC  
  //  TFile* fileMC_WZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WZ/mergeFile/WZ_merged.root");
  TFile* fileMC_WZ = new TFile("WZ.root");
  TTree* treeMC_WZ = (TTree*) fileMC_WZ -> Get("tree");

  // Fill the ZZ MC  
  //  TFile* fileMC_ZZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/ZZ/mergeFile/ZZ_merged.root");
  TFile* fileMC_ZZ = new TFile("ZZ.root");
  TTree* treeMC_ZZ = (TTree*) fileMC_ZZ -> Get("tree");

 
  // Fill the WW MC  
  // TFile* fileMC_WW = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WW/mergeFile/WW_merged.root");
  TFile* fileMC_WW = new TFile("WW.root");
  TTree* treeMC_WW = (TTree*) fileMC_WW -> Get("tree");

  // Fill the WJetsToLNu MC  
  //  TFile* fileMC_WJetsToLNu = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WJetsToLNu/mergeFile/WJetsToLNu_merged.root");
  TFile* fileMC_WJetsToLNu = new TFile("WJetsToLNu.root");
  TTree* treeMC_WJetsToLNu = (TTree*) fileMC_WJetsToLNu -> Get("tree");

  //Fill QCD MC
  // TFile* fileMC_QCD = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/QCD/mergeFile/QCD_Pt_20_MuEnrichedPt_15_merged.root");
  // TFile* fileMC_QCD = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/QCD.root");
  //TTree* treeMC_QCD = (TTree*) fileMC_QCD -> Get("tree");

  // ------------------------------------------------------------------------------

  // ------------------------------------------------------------------------------
  // Pt

  TH1F* hPtData           = new TH1F("hPtData",         "", nptBins, xbins_pt);
  
  TH1F* hPtMC_ttbarjets   = new TH1F("hPtMC_ttbarjets", "", nptBins, xbins_pt);
  TH1F* hPtMC_ztautau     = new TH1F("hPtMC_ztautau",   "", nptBins, xbins_pt);
  TH1F* hPtMC_WZ          = new TH1F("hPtMC_WZ",        "", nptBins, xbins_pt);
  TH1F* hPtMC_ZZ          = new TH1F("hPtMC_ZZ",        "", nptBins, xbins_pt);
  TH1F* hPtMC_WW          = new TH1F("hPtMC_WW",        "", nptBins, xbins_pt);
  TH1F* hPtMC_WJetsToLNu  = new TH1F("hPtMC_WJetsToLNu","", nptBins, xbins_pt);
  TH1F* hPtMC_QCD         = new TH1F("hPtMC_QCD",       "", nptBins, xbins_pt);

  // Filling events
  std::cout << "Filling Collisions Data Events\n";
  fillHistos(treeData, hPtData, mLow, mHigh, 0);
   
  std::cout << "filling TTbar+Jets events\n";
  fillHistos(treeMC_ttbarjets, hPtMC_ttbarjets, mLow, mHigh, 0);

  std::cout << "filling ZToTauTau events\n";
  fillHistos(treeMC_ztautau, hPtMC_ztautau, mLow, mHigh, 0);
     
  std::cout << "filling WZ events\n";
  fillHistos(treeMC_WZ, hPtMC_WZ, mLow, mHigh, 0);

  std::cout << "filling ZZ events\n";
  fillHistos(treeMC_ZZ, hPtMC_ZZ, mLow, mHigh, 0);

  std::cout << "filling WW events\n";
  fillHistos(treeMC_WW, hPtMC_WW, mLow, mHigh, 0);
   
  std::cout << "filling WJetsToLNu events\n";
  fillHistos(treeMC_WJetsToLNu, hPtMC_WJetsToLNu, mLow, mHigh, 0);

  //  std::cout << "filling QCD events\n";
  //fillHistos(treeMC_QCD, hPtMC_QCD, mLow, mHigh, 0);
  //  The QCD comes from data-driven estimation
     
  hPtMC_QCD -> SetBinContent( 1, 0.0314441  );
  hPtMC_QCD -> SetBinContent( 2, 0.0431494  );
  hPtMC_QCD -> SetBinContent( 3, 0.075017   );
  hPtMC_QCD -> SetBinContent( 4, 0.0446895  );
  hPtMC_QCD -> SetBinContent( 5, 0.026708   );
  hPtMC_QCD -> SetBinContent( 6, 0.0457004  );
  hPtMC_QCD -> SetBinContent( 7, 0.012335   );
  hPtMC_QCD -> SetBinContent( 8, 0.014894   );
  hPtMC_QCD -> SetBinContent( 9, 0.0280336  );
  hPtMC_QCD -> SetBinContent(10, 0.0207846  );
  hPtMC_QCD -> SetBinContent(11, 0.00949111 );
  hPtMC_QCD -> SetBinContent(12, 0.00197284 );
  hPtMC_QCD -> SetBinContent(13, 0.000230449);
  hPtMC_QCD -> SetBinContent(14, 0);
  hPtMC_QCD -> SetBinContent(15, 0);
  hPtMC_QCD -> SetBinContent(16, 0);
  hPtMC_QCD -> SetBinContent(17, 0);
  hPtMC_QCD -> SetBinContent(18, 0);

      
  // ------------------------------------------------------------------------------
  TH1F* hPtMC_allBkg    = new TH1F("hPtMC_allBkg",   "", nptBins, xbins_pt);
  

  hPtMC_ttbarjets   ->Sumw2();
  hPtMC_ztautau     ->Sumw2();
            
  hPtMC_WW          ->Sumw2();
  hPtMC_WZ          ->Sumw2();
  hPtMC_ZZ          ->Sumw2();
  hPtMC_WJetsToLNu  ->Sumw2();
  hPtMC_QCD         ->Sumw2();
               
  hPtMC_allBkg      ->Sumw2();
 
  // ------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------
  // Rescaling the MC
  hPtMC_ttbarjets -> Scale(intLumi_invpb*xsec_ttbarjets/nMCEvTTbarJets);
  hPtMC_ztautau   -> Scale(intLumi_invpb*xsec_ztautau/nMCEvZtautau);
  hPtMC_WW        -> Scale(intLumi_invpb*xsec_WW/nMCEvWW);
  hPtMC_WZ        -> Scale(intLumi_invpb*xsec_WZ/nMCEvWZ);
  hPtMC_ZZ        -> Scale(intLumi_invpb*xsec_ZZ/nMCEvZZ);
  hPtMC_WJetsToLNu -> Scale(intLumi_invpb*xsec_WJetsToLNu/nMCEvWJetsToLNu);
  hPtMC_QCD       -> Scale(intLumi_invpb*xsec_QCD/nMCEvQCD);


  // sum all bkgs
  hPtMC_allBkg  -> Add(hPtMC_ttbarjets);  
  hPtMC_allBkg  -> Add(hPtMC_ztautau);
  hPtMC_allBkg  -> Add(hPtMC_WJetsToLNu);
  hPtMC_allBkg  -> Add(hPtMC_QCD);
  hPtMC_allBkg  -> Add(hPtMC_ZZ);
  hPtMC_allBkg  -> Add(hPtMC_WZ); 
  hPtMC_allBkg  -> Add(hPtMC_WW);
           
 

  // ------------------------------------------------------------------------------
  // 1. Subtract the bkg
  hPtData -> Sumw2();
  for(int i(0); i<18; i++)
  {
    cout << "hPtData : " << hPtData->GetBinContent(i+1) << endl;
  }
    
  hPtData -> Add(hPtMC_allBkg, -1);
  
  
  // 2. Correct for  Eff and scale factors by using accEff_Lumi.C
  // For Nadeesha: probably you will need to change the location
  // TFile* fEffvsPt = new TFile("histo_EffvsPt.root");
  TFile* fEffvsPt = new TFile("histo_EffvsPt_withnewSF.root"); // this is the one with new SF from GP. 
  
  //  TH1F* hEffvsPt = (TH1F*) fEffvsPt -> GetListOfPrimitives()->At(1);
  TH1F *hEffvsPt = (TH1F*)fEffvsPt->Get("hEff");

  hEffvsPt->Sumw2();
  hPtData -> Divide(hEffvsPt);
    
  // --- test ---
  //   std::cout << "DATA\n";
  //   for(int bin=1; bin<=hPtData->GetNbinsX(); bin++){

  //     printf("%2d)  %8.8f +/- %8.8f +/- %8.8f | \n", 
  //            bin, 
  //            hPtData->GetBinContent(bin), hPtData->GetBinError(bin),
  //            sqrt(hPtData->GetBinContent(bin)) );
  //   }
  
  // 3. Unfolding Resolution
  hXSec -> Sumw2();
  applyMatrix(hXSec, hPtData);

  // FSR correction
  hXSecFSR -> Sumw2();
  applyMatrixFSR(hXSecFSR, hXSec);
 
  //=========================================================================
  // without FSR correction=============================  
  //4. Correct by Bin Width and Normalize
  /*
    TGraphAsymmErrors* gXSec= new TGraphAsymmErrors(hPtData->GetNbinsX());
  
  std::cout << "Stats Only\n";
   for(int bin=1; bin<=hPtData->GetNbinsX(); bin++){

    double value = hPtData->GetBinContent(bin)/(hPtData->GetBinWidth(bin) * hPtData -> Integral()); 
    double sigma = hPtData->GetBinError(bin)  /(hPtData->GetBinWidth(bin) * hPtData -> Integral()); 
    
    printf("%2d) %5.2f | %8.8f +/- %8.8f \n", 
          bin, hPtData->GetXaxis()->GetBinCenter(bin), value, sigma);
 
    gXSec->SetPoint(bin-1,hPtData->GetXaxis()->GetBinCenter(bin),value);
  //  // For Nadeesha: here I put the error on the x-axis 
    gXSec->SetPointError(bin-1,error_xbins[bin-1],error_xbins[bin-1],sigma,sigma);
   }
  */
  //============================================================================
  
  //4. Correct by Bin Width and Normalize
  TGraphAsymmErrors* gXSec= new TGraphAsymmErrors(hXSecFSR->GetNbinsX());
  
  std::cout << "Stats Only\n";
  for(int bin=1; bin<=hXSecFSR->GetNbinsX(); bin++){

    double value = hXSecFSR->GetBinContent(bin)/(hXSecFSR->GetBinWidth(bin) * hXSecFSR -> Integral()); 
    double sigma = hXSecFSR->GetBinError(bin)  /(hXSecFSR->GetBinWidth(bin) * hXSecFSR -> Integral()); 
     
    //for systematics fractional error on diffXSec
    //   double FerrorXSec = sigma/value;
    //    printf("%8.8f \n",
    //    FerrorXSec);  
    //end systematics fractional error on diffXSec
   
      printf("%2d) %5.2f | %8.8f +/- %8.8f | %6.6f | %6.6f\n", 
    //bin, hXSecFSR->GetXaxis()->GetBinCenter(bin), value, sigma);
      bin, point[bin-1], value, sigma, error_xbinsl[bin-1], error_xbinsh[bin-1] );
    
    //gXSec->SetPoint(bin-1,hXSecFSR->GetXaxis()->GetBinCenter(bin),value);
    //For Nadeesha: here I put the error on the x-axis 
    //gXSec->SetPointError(bin-1,error_xbins[bin-1],error_xbins[bin-1],sigma,sigma);
      gXSec->SetPoint(bin-1,point[bin-1],value);
  
      gXSec->SetPointError(bin-1,fabs(error_xbinsl[bin-1]),fabs(error_xbinsh[bin-1]),sigma,sigma);
  }


  

  // some style on the graph
  //  gXSec -> SetTitle("");
  // gXSec -> GetXaxis() -> SetTitle("q_{T} [GeV]");
  //gXSec -> GetYaxis() -> SetTitle("1/#sigma d#sigma/dp_{T} [GeV]^{-1}");
  //gXSec -> GetYaxis() -> SetTitleOffset(1.2);
  //gXSec -> SetMarkerSize(0.7);
  //gXSec -> SetMarkerStyle(20);


  //-----------------------------------------------------------------------------------------------
  TGraphAsymmErrors* grPowheg = Powheg(); 
  TGraphAsymmErrors* grResbos = Resbos();
  TGraphAsymmErrors* grMadgraph = Madgraph();
  TGraphAsymmErrors* grFEWZ = FEWZ();

 

  //-----------------------------------------------------------------------------------------------
  
  TGraphAsymmErrors* gXSecWithSyst= new TGraphAsymmErrors(hXSecFSR->GetNbinsX());
  fillSystErrorEff(gXSec,gXSecWithSyst);

  cout << "========= Print gXSecWithSyst values=========" << endl;
  gXSecWithSyst->Print();
  //for(int i(0); i<18; i++)
  //{
  //  cout << "gXSec point : " << gXSecWithSyst->GetBinContent(i) << endl;
  //}

  gXSecWithSyst -> SetTitle("");
  
  //gXSecWithSyst -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  gXSecWithSyst -> GetXaxis() -> SetTitle("p_{T}^{Z} [GeV]");
  gXSecWithSyst -> GetXaxis() -> SetTitleSize(0.04);
  gXSecWithSyst -> GetXaxis() -> SetLabelSize(0.03);

  //gXSecWithSyst -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} [GeV/c]^{-1}");
  gXSecWithSyst -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{Z}} [GeV]^{-1}");
  gXSecWithSyst -> GetYaxis() -> SetTitleSize(0.04);
  gXSecWithSyst -> GetYaxis() -> SetLabelSize(0.03);
  gXSecWithSyst -> GetYaxis() -> SetTitleOffset(1.02);

  //TLegend* tl11 = new TLegend(0.5,0.2,0.4,0.48);
  //TLegend* tl11 = new TLegend(0.5,0.25,0.4,0.53);
  //TLegend* tl11 = new TLegend(0.6,0.3,0.4,0.58);
  TLegend* tl11 = new TLegend(0.6,0.35,0.4,0.56);
  tl11 -> SetTextSize(0.045);
  tl11->SetTextFont(42);
  tl11->SetFillColor(0);
  tl11->SetBorderSize(0);
  tl11->AddEntry(grPowheg,"Pythia/Powheg", "f");
  tl11->AddEntry(gXSecWithSyst, "data",    "lp");
  
  //TLegend* tl1 = new TLegend(0.6,0.2,0.4,0.48);
  TLegend* tl1 = new TLegend(0.6,0.3,0.4,0.47);
  tl1 -> SetTextSize(0.045);
  tl1->SetTextFont(42);
  tl1->SetFillColor(0);
  tl1->SetBorderSize(0);
  tl1->AddEntry(grResbos,"Resbos", "f");
  tl1->AddEntry(gXSecWithSyst, "data",    "lp");

  //TLegend* tl2 = new TLegend(0.5,0.2,0.4,0.48);
  //TLegend* tl2 = new TLegend(0.6,0.11,0.4,0.39);
  TLegend* tl2 = new TLegend(0.6,0.23,0.4,0.39);
  tl2 -> SetTextSize(0.045);
  tl2->SetTextFont(42);
  tl2->SetFillColor(0);
  tl2->SetBorderSize(0);
  tl2->AddEntry(grMadgraph,"Madgraph", "f");
  tl2->AddEntry(gXSecWithSyst, "data", "lp");

  //TLatex *pave = new TLatex(0.25,0.1, "|#eta|<2.1, p_{T}>20 GeV");
  //TLatex *pave = new TLatex(0.26,0.08, "|#eta|<2.1, p_{T}>20 GeV");
  TLatex *pave = new TLatex(0.26,0.08, "");
  pave -> SetNDC();
  pave -> SetTextColor(kBlue);
  pave -> SetTextFont(42);
  pave->SetTextSize(0.05);
  pave -> SetTextAlign(22);

  gXSecWithSyst -> SetMarkerSize(0.8);
  gXSecWithSyst -> SetMarkerStyle(20);
  gXSecWithSyst -> SetMinimum(0.000001);
  gXSecWithSyst -> GetXaxis()->SetRangeUser(0,800);
  // gXSecWithSyst -> SetLineColor(1);
  gXSecWithSyst-> SetMarkerColor(kBlack);
  gXSecWithSyst -> SetLineColor (kBlack);
  //gXSecWithSyst -> Draw("AP");
 
  TH1F* hDummyResbos = new TH1F("hDummyResbos","", nptBins, xbins_pt);

  for(int j=1; j<=hDummyResbos->GetNbinsX(); j++){
    double x,y;
    grResbos->GetPoint(j-1,x,y);
    hDummyResbos->SetBinContent(j,y);
    hDummyResbos->SetBinError  (j-1,grResbos->GetErrorY(j));     
   }

  // draw powheg
  TCanvas* c1 = new TCanvas("c1","",0,0,700,500);
  c1->SetFillColor(0);
  c1->SetLogx();
  c1->SetLogy();

  hDummyResbos   -> SetLineColor (kBlue);
  //  hDummyPowheg   -> SetFillColor (kGreen+2);
  hDummyResbos   -> SetLineWidth (2);
  //hDummyResbos->Draw("histo same");

  // draw the white band and the green horizontal line
  const double x11[4]={2.46,2.55,2.55,2.46};
  const double y11[4]={0,0,0.02235400,0.02235400};
  TGraph* band111=new TGraph(4,x11,y11);
  band111->SetFillColor(kWhite);
  //band111->Draw("F SAME");



  const double x1[4]={2.45,2.55,2.55,2.45};
  const double y1[4]={0,0,0.02957739,0.02957739};
  TGraph* band11=new TGraph(4,x1,y1);
  band11->SetFillColor(kWhite);
  //band11->Draw("F SAME");
 
  const double x2[4]={2.46,2.55,2.55,2.46};
  const double y2[4]={0,0,0.02838694,0.02838694};
  TGraph* band12=new TGraph(4,x2,y2);
  band12->SetFillColor(kWhite);
  //band12->Draw("F SAME");
/*
  //   double Lxmax = canvas->GetUxmax();
  TLine* line111 = new TLine(0,0.02235400,2.52,0.02235400);
  line111 -> SetLineColor( kGreen+2 );
  line111 -> SetLineWidth( 2 );
  line111 -> Draw("L same");

  TLine* line11 = new TLine(0,0.02957739,2.52,0.02957739);
  line11 -> SetLineColor( kBlue );
  line11 -> SetLineWidth( 2 );
  line11 -> Draw("L same");

  TLine* line12 = new TLine(0,0.02838694,2.52,0.02838694);
  line12 -> SetLineColor( kViolet );
  line12 -> SetLineWidth( 2 );
  line12 -> Draw("L same");
*/
  PrintIt(c1, cTitleExtd);
  //tl1->Draw("same");
  //pave->Draw();



 
  // set powheg graph style
  //  grPowheg -> SetMarkerSize(0.8);
  //  grPowheg -> SetMarkerStyle(21);
  //grResbos -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  grResbos -> GetXaxis() -> SetTitle("p_{T}^{Z} [GeV]");
  //grResbos -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} [GeV/c]^{-1}");
  grResbos -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{Z}} [GeV]^{-1}");
  grResbos -> GetYaxis() -> SetTitleOffset(0.69);
  grResbos   -> SetLineColor (kBlue);
  grResbos   -> SetFillColor (kBlue);
  grResbos   -> SetLineWidth (2);


  // set madgraph style
  //grMadgraph -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} [GeV/c]^{-1}");
  grMadgraph -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{Z}} [GeV]^{-1}");
  //grMadgraph -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  grMadgraph -> GetXaxis() -> SetTitle("p_{T}^{Z} [GeV]");
  grMadgraph -> GetYaxis() -> SetTitleOffset(0.69);
  grMadgraph -> SetFillColor (kViolet);
  grMadgraph -> SetLineColor (kViolet);
  grMadgraph -> SetLineWidth(2);

  //grPowheg -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  grPowheg -> GetXaxis() -> SetTitle("p_{T}^{Z} [GeV]");
  //grPowheg -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} [GeV/c]^{-1}");
  grPowheg -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{Z}} [GeV]^{-1}");
  grPowheg -> GetYaxis() -> SetTitleOffset(0.69);
  grPowheg   -> SetLineColor (kGreen+2);
  grPowheg   -> SetFillColor (kGreen+2);
  grPowheg   -> SetLineWidth (2);


    //grResbos->Draw(" C same");
  
    //gXSec -> Draw("p same");

 

  TCanvas* c11 = new TCanvas("c11","",800,0,800,700);
  c11->SetFillColor(0);
  //DrawWithRatio(c11, cTitle, gXSecWithSyst, grPowheg, band111, line111, tl11, pave, pave1);
  
  TCanvas* c2 = new TCanvas("c2","",800,0,800,700);
  c2->SetFillColor(0);
  //DrawWithRatio(c2, cTitle, gXSecWithSyst, grResbos, band11, line11, tl1, pave, pave1);
  //DrawWithRatio(c11, cTitle, gXSecWithSyst, grResbos, band11, line11, tl1, pave, pave1);

  TCanvas* c3 = new TCanvas("c3","",0,0,800,700);
  c3->SetFillColor(0);
  //DrawWithRatio(c3, cTitle, gXSecWithSyst, grMadgraph, band12, line12, tl2, pave, pave1);
  //DrawWithRatio(c11, cTitle, gXSecWithSyst, grMadgraph, band12, line12, tl2, pave, pave1);
  
  
  //TCanvas* ctot = new TCanvas("ctot","",800,0,800,800);

  TCanvas* ctot = new TCanvas("ctot","",50,50,W,H);
  ctot->SetFillColor(0);
  DrawWithRatioTotal(ctot, cTitle, gXSecWithSyst, 
      				   grPowheg, 
				   grResbos, 
				   grMadgraph, 
				   grFEWZ, 
				   
				   band111,  band11,   band12, 
				   tl11,     tl1,      tl2,
				   
				   pave);
 
     ctot->SaveAs("diffXSec_PowhegResBosFEWZ."+format);
}



//---------------------------------------------------------------------
TGraphAsymmErrors* Powheg(){

  const int n = 18;

   double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
   double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
   double x[n];

   ///Powheg Pre fsr numbers
   double y[18] = {0.0232823816,
                   0.0534414648, 
                   0.0556310835, 
                   0.0454828550, 
                   0.0356700536, 
                   0.0280575476, 
                   0.0226186451, 
                   0.0182287311, 
                   0.0118469402, 
                   0.0064591971, 
                   0.0037984679, 
                   0.0019372698, 
                   0.0008350090, 
                   0.0003902566, 
                   0.0001503083, 
                   0.0000509369, 
                   0.0000184618, 
                   0.0000014680};

   //// stat error is normalized with Toy RMS method.
   double errorStat_y[18] = {0.0001395273, 
                             0.0002032366, 
                             0.0002069700, 
                             0.0001884388, 
                             0.0001705318, 
                             0.0001524742, 
                             0.0001370749,
                             0.0001241799, 
                             0.0000484067, 
                             0.0000366260, 
                             0.0000284210, 
                             0.0000144427, 
                             0.0000095287, 
                             0.0000065614, 
                             0.0000028981, 
                             0.0000016920, 
                             0.0000008289, 
                             0.0000000970};

   //// powheg norm PDF error in %, convert it to number
   double errorPDF_y[18] = {0.412089 , 
                            0.31525  , 
                            0.209474 , 
                            0.188146 , 
                            0.250735 , 
                            0.21012  , 
                            0.351397 ,
                            0.284557 , 
                            0.304525 , 
                            0.315051 , 
                            0.356816 , 
                            0.619421 , 
                            1.16534  , 
                            1.56789  , 
                            2.29144  , 
                            2.92386  , 
                            4.04627  , 
                            6.01494  };
  
  double error_y[n];
  for (int i=0;i<n;i++)
  { 
    errorPDF_y[i]= 0.01*errorPDF_y[i]*y[i];
    error_y[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_y[i]*errorPDF_y[i]); 
     //std::cout<<"Powhed:  errorPDF_y[i]\t"<<errorPDF_y[i]<<"\terrorStat_y[i]\t"<<errorStat_y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
     //std::cout<<"Powhed:  errorPDF_y[i]\t"<<100*errorPDF_y[i]/y[i]<<"\terrorStat_y[i]\t"<<100*errorStat_y[i]/y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
     //printf("Powheg PDF% : %.2f \t Stat % : %.2f \n ",100*errorPDF_y[i]/y[i] , 100*errorStat_y[i]/y[i]);
    std::cout << "Powheg PDF : " << 100*errorPDF_y[i]/y[i] << "\t Stat : " << 100*errorStat_y[i]/y[i] << endl;
  }

   //double error_x[18];

   //double error_x[18];
   double error_xlow[n];
   double error_xhigh[n];
  
   for (int i=0;i<n;i++) {
     error_xlow[i] = (point[i] - bin[i]);
     error_xhigh[i] = (bin[i+1]-point[i]);
   }

   //for (int i=0;i<n;i++) error_x[i]=0;

   //TGraphAsymmErrors* grPowheg = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
   TGraphAsymmErrors* grPowheg = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_y,error_y);


   return grPowheg;
}

TGraphAsymmErrors* Resbos(){
  
  const int n = 18;
  
   double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
   double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
   double x[n];
   

   double y[18] = {
 // Nadeesha PAS
        0.02957739, 
	0.05620224, 
	0.04961592, 
	0.04050349, 
	0.03325381, 
	0.02725463, 
	0.02253389, 
	0.01862596, 
	0.01237984, 
	0.00664593, 
	0.00393008, 
	0.00200304, 
	0.00085508, 
	0.00041032, 
	0.00017656, 
	0.00006332, 
	0.00002423, 
	0.00000115};

     // Hammid resbos cp
             //0.03223331, 
             //0.06286445,
             //0.05278694,
             //0.04188681,
             //0.03321407,
             //0.02650710,
             //0.02134491,
             //0.01716620,
             //0.01144344,
             //0.00580740,
             //0.00346440,
             //0.00182813,
             //0.00082492,
             //0.00040398,
             //0.00017512,
             //0.00006440,
             //0.00002296,
             //0.00000211}; 

	// Hammid resbos p
		  // 0.03801042,	
                  // 0.05823305,	
                  // 0.05171227,	
                  // 0.04158588,	
                  // 0.03318528,	
                  // 0.02647431,	
                  // 0.02130071,	
                  // 0.01728264,	
                  // 0.01103833,	
                  // 0.00592830,	
                  // 0.00353564,	
                  // 0.00188739,	
                  // 0.00085679,	
                  // 0.00042180,	
                  // 0.00018299,	
                  // 0.00006782,	
                  // 0.00002363,	
                  // 0.00000213};	

   //double error_y[18] = {0.00514001, 
   //                      0.00707116,
   //                      0.00664392, 
   //                      0.00600288,
   //                      0.00543919, 
   //                      0.00492418, 
   //                      0.00447746, 
   //                      0.00407074, 
   //                      0.00165936, 
   //                      0.00121580, 
   //                      0.00093494, 
   //                      0.00047197, 
   //                      0.00030837, 
   //                      0.00021361, 
   //                      0.00009908, 
   //                      0.00005934, 
   //                      0.00002997, 
   //                      0.00000270};

   double error_y[18] = {0, 
                         0,
                         0, 
                         0,
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0, 
                         0};

   //double error_x[18];
   double error_xlow[n];
   double error_xhigh[n];
  
   for (int i=0;i<n;i++) {
     error_xlow[i] = (point[i] - bin[i]);
     error_xhigh[i] = (bin[i+1]-point[i]);
   }

   //for (int i=0;i<n;i++) error_x[i]=0;

   //TGraphAsymmErrors* grResbos = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
   TGraphAsymmErrors* grResbos = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_y,error_y);
   
   
   return grResbos;
}

 


TGraphAsymmErrors* Madgraph(){
  
  const int n = 18;
  
  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
  double x[n];
  //  for (int i=0;i<n;i++) {
  //  x[i] = bin[i] + (bin[i+1]-bin[i])/2;
  // }

/*
  // PostFSR and StatErr
  double y[18] = {0.02838694, 0.05640728, 0.05274386, 0.04305959, 0.03456624, 0.02771716, 0.02242516, 0.01798233, 0.01161362, 0.00624185, 0.00375409, 0.00194832, 0.00084676, 0.00041085, 0.00018461, 0.00006234, 0.00002194, 0.00000187};
  double error_y[18] = {0.00022293, 0.00031362, 0.00030326, 0.00027401, 0.00024550, 0.00021984, 0.00019774, 0.00017707, 0.00007115, 0.00005216, 0.00004045, 0.00002061, 0.00001359, 0.00000946, 0.00000449, 0.00000261, 0.00000126, 0.00000015};
*/
  // PreFSR and Total Uncer(sqrt(StatErr^2 + PDFErr^2))
  double y[18] = { 0.0296683  , 
                   0.058359   , 
                   0.0532516  , 
                   0.042838   , 
                   0.03381    , 
                   0.0268681  , 
                   0.0216179  , 
                   0.0171872  , 
                   0.0113335  , 
                   0.00623927 , 
                   0.00378249 , 
                   0.0019837  , 
                   0.000862846,
                   0.000418884, 
                   0.000188615, 
                   6.45009e-05, 
                   2.23858e-05, 
                   1.91878e-06};
  double errorStat_y[18] = {0.0002159060, 
                            0.0002914961, 
                            0.0002808667, 
                            0.0002537420, 
                            0.0002301006, 
                            0.0002065859, 
                            0.0001854926, 
                            0.0001669025, 
                            0.0000656459, 
                            0.0000498163, 
                            0.0000392100, 
                            0.0000201929, 
                            0.0000133864, 
                            0.0000093946, 
                            0.0000044844, 
                            0.0000026313, 
                            0.0000012616, 
                            0.0000001533};
       /////normalized PDF syst in %, convert it to number
  double errorPDF_y[18] = {1.87987 , 
                           1.74473 , 
                           1.31112 , 
                           1.15229 , 
                           0.797034, 
                           0.369078, 
                           0.275328, 
                           0.691266, 
                           1.29233 , 
                           1.74525 , 
                           1.93177 , 
                           2.09596 , 
                           2.31831 , 
                           2.62671 , 
                           2.7677  , 
                           3.15532 , 
                           3.78032 , 
                           4.67155 };

 
  double error_y[n];
  for (int i=0;i<n;i++)
  { 
    errorPDF_y[i]= 0.01*errorPDF_y[i]*y[i];
    error_y[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_y[i]*errorPDF_y[i]); 
     //std::cout<<"Madgraph:   errorPDF_y[i]\t"<<errorPDF_y[i]<<"\terrorStat_y[i]\t"<<errorStat_y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
     //std::cout<<"Madgraph:  errorPDF_y[i]\t"<<100*errorPDF_y[i]/y[i]<<"\terrorStat_y[i]\t"<<100*errorStat_y[i]/y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
     //printf("Madgraph PDF % : %.2f \t Stat % : %.2f \n ",100*errorPDF_y[i]/y[i] , 100*errorStat_y[i]/y[i]);
  }


  //double error_x[18];
  double error_xlow[n];
  double error_xhigh[n];
  
  for (int i=0;i<n;i++) {
    error_xlow[i] = (point[i] - bin[i]);
    error_xhigh[i] = (bin[i+1]-point[i]);
  }

  //for (int i=0;i<n;i++) error_x[i]=0;
  
  //TGraphAsymmErrors* grMadgraph = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
  TGraphAsymmErrors* grMadgraph = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_y,error_y);
  
  
  return grMadgraph;
}

TGraphAsymmErrors* FEWZ(){

  const int n = 18;

   double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
   double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
   double x[n];


   //double y[18] = {0,0,0,0,0,0,0,0,0.010946908,0.005861624,0.00351576,0.00186525,0.000836974,0.000406356,0.000176403,6.37816E-05,2.25044E-05,2.08287E-06};
   //double error_ylow[18] = {0,0,0,0,0,0,0,0,0.000260501,0.00013212,7.60136E-05,3.86252E-05,1.65916E-05,7.78994E-06,3.33908E-06,1.20556E-06,3.89959E-07,4.13802E-08};
   //double error_yhigh[18] = {0,0,0,0,0,0,0,0,0.00031118,0.00015371,8.64209E-05,4.26788E-05,1.75198E-05,8.01257E-06,3.36195E-06,1.27319E-06,4.68827E-07,5.67557E-08};
double y[18];
y[0] = 100;//-0.08005877; 
y[1] = 100;//0.12748182 ; 
y[2] = 0.07877656 ; 
y[3] = 0.05332600 ; 
y[4] = 0.03845661 ; 
y[5] = 0.02904814 ; 
y[6] = 0.02265475 ; 
y[7] = 0.01805065 ; 
y[8] = 0.01135670 ; 
y[9] = 0.00594296 ; 
y[10] = 0.00353301; 
y[11] = 0.00185758; 
y[12] = 0.00082569; 
y[13] = 0.00039315; 
y[14] = 0.00016611; 
y[15] = 0.00005898; 
y[16] = 0.00001959; 
y[17] = 0.00000179; 
 
  double errorStat_y[18] = {0.0,//0.00094095, 
                            0.0,//0.00039520, 
                            0.00024179,
                            0.00016486,
                            0.00012112,
                            0.00009292,
                            0.00007464,
                            0.00006809,
                            0.00003028,
                            0.00001669,
                            0.00001118,
                            0.00000535,
                            0.00000286,
                            0.00000170,
                            0.00000067,
                            0.00000082,
                            0.00000031,
                            0.00000004};
       /////normalized PDF syst in %, convert it to number
  double errorPDF_ylow[18] = { 0.0,//-6.7678 , 
                               0.0,//1.0464, 
                                 0.7600 ,
                                 1.0219 ,
                                 1.2630 ,
                                 1.4623 ,
                                 1.6461 ,
                                 1.8262 ,
                                 2.1995 ,
                                 2.7651 ,
                                 3.2279 ,
                                 3.8587 ,
                                 4.6700 ,
                                 5.4480 ,
                                 6.4558 ,
                                 7.8312 ,
                                 9.3764 ,
                                 11.4845}; 
  double errorPDF_yhigh[18] = {0.0,// -6.4288, 
                               0.0,// 1.8289 , 
                               1.3007  ,
                               1.2656  ,
                               1.3277  ,
                               1.4084  ,
                               1.5036  ,
                               1.6060  ,
                               1.8347  ,
                               2.2177  ,
                                2.5401 ,
                                2.9723 ,
                                3.4925 ,
                                3.9540 ,
                                4.4901 ,
                                5.1299 ,
                                5.7526 ,
                                6.6129 }; 
  double errorScale[18] = {
	     0.0,//0.00842947        , 
	     0.0,//0.0106515         ,
	     0.00266521 ,
	     0.000428617,
	     0.000114661,
	     0.000205532,
	     0.000345275,
	     0.000297941,
	     0.000323374,
	     0.000210461,
	     0.000147318,
	     9.51447e-05,
	     4.73258e-05,
	     2.47657e-05,
	     1.03102e-05,
	     4.30353e-06,
	     1.45553e-06,
	     1.5117e-07 };
  
  double error_yhigh[n];
  double error_ylow[n];
  for (int i=0;i<n;i++)
  { 
    errorPDF_ylow[i]= 0.01*errorPDF_ylow[i]*y[i];
    errorPDF_yhigh[i]= 0.01*errorPDF_yhigh[i]*y[i];
    error_ylow[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_ylow[i]*errorPDF_ylow[i] + errorScale[i]*errorScale[i]); 
    error_yhigh[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_yhigh[i]*errorPDF_yhigh[i] + errorScale[i]*errorScale[i]); 
     //std::cout<<"FEWZ:   errorPDF_ylow[i]\t"<<errorPDF_ylow[i]<<"\terrorStat_y[i]\t"<<errorStat_y[i]<<"\t error_ylow[i]\t"<<error_ylow[i]<< std::endl;
     //std::cout<<"FEWZ:  errorPDF_y[i]\t"<<100*TMath::Max(errorPDF_ylow[i],errorPDF_yhigh[i])/y[i]<<"\terrorStat_y[i]\t"<<100*errorStat_y[i]/y[i]<<"\terrorScale_y[i]\t"<<100*errorScale[i]/y[i]<< std::endl;
     //printf("FEWZ PDF % : %.2f \t Stat % : %.2f Scale % : %.2f \n ",100*TMath::Max(errorPDF_ylow[i],errorPDF_yhigh[i])/y[i] , 100*errorStat_y[i]/y[i], 100*errorScale[i]/y[i]);
  }
 
   //double error_x[18];
   double error_xlow[n];
   double error_xhigh[n];
  
   for (int i=0;i<n;i++) {
     error_xlow[i] = (point[i] - bin[i]);
     error_xhigh[i] = (bin[i+1]-point[i]);
   }

   //double error_y[18];
   //for (int i=0;i<n;i++) error_x[i]=0;
   //for (int i=0;i<n;i++) error_y[i] = (error_yhigh[i]>error_ylow[i]) ? error_yhigh[i] : error_ylow[i];

   //TGraphAsymmErrors* grFEWZ = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
   TGraphAsymmErrors* grFEWZ = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_ylow,error_yhigh);


   return grFEWZ;
}
// ------------------------------------------------------------------------------


void fillHistos(TTree* tree,
                TH1F *hPt,
                double mLow,
                double mHigh,
                int mode) {


  typedef std::pair<int,int> pairOfInt;
  set<pairOfInt> uniqueSelectedEvents;

  typedef std::pair<int,int> pairOfInt;
  set<pairOfInt> uniqueEventsReco;


  _EventInfo eventInfo;
  tree->SetBranchAddress("eventInfo", &eventInfo);

  _MuonInfo reco1, reco2;

  tree->SetBranchAddress("reco1", &reco1);
  tree->SetBranchAddress("reco2", &reco2);

  float recoCandMass, recoCandPt, recoCandY;

  tree->SetBranchAddress("recoCandMass", &recoCandMass);
  tree->SetBranchAddress("recoCandPt"  , &recoCandPt );
  tree->SetBranchAddress("recoCandY"  ,  &recoCandY );

  float vertexNormChiSquare, angleDiMuons;

  tree->SetBranchAddress("vertexNormChiSquare", &vertexNormChiSquare);
  tree->SetBranchAddress("angleDiMuons"       , &angleDiMuons );

  _VertexInfo vertexInfo;
  tree->SetBranchAddress("vertexInfo", &vertexInfo);

  cout<<"Loop over the " << tree->GetEntries() << " entries ...\n";
  for (int iEvt=0; iEvt < tree->GetEntries(); iEvt++) {
   
    if ( (iEvt % 500000)==0 ) cout << "event " << iEvt << endl;
    tree -> GetEntry(iEvt);
/*
    // additional selection cuts
    if (reco1.charge == reco2.charge) continue;
    

    if (recoCandMass <  mLow) continue;
    if (recoCandMass > mHigh) continue;

    
    if (reco1.pt < 20 || reco2.pt < 20) continue;
    if (fabs(reco1.eta) > 2.1 || fabs(reco2.eta) > 2.1) continue;
     
    if (angleDiMuons > TMath::Pi()-0.02) continue;

    // both muons tight
    if( !isKinTight_2012(reco1) || !isKinTight_2012(reco2) ) continue;
  
    if( !reco1.isHltMatched[0] && !reco2.isHltMatched[0]) continue; 
    // remove duplicates
    pairOfInt runEvent(eventInfo.run,eventInfo.event);
    if( !uniqueEventsReco.insert( runEvent ).second ) continue;
*/
    
    // additional selection cuts
    if (reco1.charge == reco2.charge) false;
    

    if (recoCandMass <  mLow) false;
    if (recoCandMass > mHigh) false;

    
    if (reco1.pt < 20 || reco2.pt < 20) false;
    if (fabs(reco1.eta) > 2.1 || fabs(reco2.eta) > 2.1) false;
     
    if (angleDiMuons > TMath::Pi()-0.02) false;

    // both muons tight
    if( !isKinTight_2012(reco1) || !isKinTight_2012(reco2) ) false;
  
    if( !reco1.isHltMatched[0] && !reco2.isHltMatched[0]) false; 
    // remove duplicates
    pairOfInt runEvent(eventInfo.run,eventInfo.event);
    if( !uniqueEventsReco.insert( runEvent ).second ) false;
 
    hPt       -> Fill( recoCandPt );

  }
  std::cout << " DONE!" << std::endl;
  return;
}

/***
//------------------------------------------------------------------------------
// Draw projections and residuals
//------------------------------------------------------------------------------
void DrawWithRatio(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, TGraphAsymmErrors* gDen, TGraph* band, TLine* line,
		   TLegend* leg, TLatex* tex, TLatex* tex1) {
  
  // sanity check
  if (gNum->GetN() != gDen->GetN()){
    std::cout<< " *** Error: number of points not consistent between the"
             << " two graphs -> Exit!\n";
    return;
  }

  TH1F *hNum = ConvertGraphToHisto(gNum);
  TH1F *hDen = ConvertGraphToHisto(gDen);

  std::cout << "hDen->GetNbinsX() = " << hDen->GetNbinsX() << std::endl;
  for (int i=1; i<hDen->GetNbinsX()+1; i++)
    std::cout << i << ") content=\t " << hDen->GetBinContent(i) << std::endl;

  //hNum -> GetXaxis() -> SetLimits(-100,10000);
  //hNum -> Sumw2();
  //hDen -> Sumw2();

  TH1F *hPull = (TH1F*) hNum -> Clone("hPull");
  hPull->Sumw2();
  hPull->Reset();

  // Print the ratio
  std::cout << std::endl << std::endl;
  for (int iBin=0;iBin<gNum->GetN();iBin++) {
    Double_t x,yNum,yDen, error_high, error_low, error;
    gNum->GetPoint(iBin,x,yNum); //data
    gDen->GetPoint(iBin,x,yDen);

    
    // double sigmaN, sigmaD;
    //error_highN = gNum->GetErrorYhigh(iBin);
    //error_lowN  = gNum->GetErrorYlow(iBin);
    //if (error_highN > error_lowN) sigmaN = error_highN;
    //else sigmaN = error_lowN;

    //error_highD = gDen->GetErrorYhigh(iBin);
    //error_lowD  = gDen->GetErrorYlow(iBin);
    //if (error_highD > error_lowD) sigmaD = error_highD;
    //else sigmaD = error_lowD;
    
   //--this is for the final PAS--code start for below plot

    error_high = gNum->GetErrorYhigh(iBin);
    error_low  = gNum->GetErrorYlow(iBin);
    if (error_high > error_low) error = error_high;
    else error = error_low;
    double resid = (error) ? ((yNum-yDen)/error) : 0.0;
    if (fabs(resid) > 100) resid = 10;

    hPull->SetBinContent(iBin+1,resid);
    cout << iBin << ", " << resid << endl;
    
    //code end
    
    // double ratio = yNum/yDen;

    //double sigmaR = sqrt(((1/yNum)*(1/yNum))*(sigmaN*sigmaN)+ ((yNum/(yDen*yDen))*(yNum/(yDen*yDen)))*(sigmaD*sigmaD));
    //  std::cout << iBin << ") " << yNum << "/" << yDen << ", Ratio= " << ratio << std::endl;
    // hPull->SetBinContent(iBin+1,ratio);
    //hPull->SetBinError(iBin+1,sigmaR);
    
  }
  std::cout << std::endl << std::endl;

  //----------------------------------------------------------------------------
  // Create the pads
  //----------------------------------------------------------------------------
  TPad* pad1;
  TPad* pad2;

  pad1 = new TPad("pad1","This is pad1",0.02,0.30,0.98,0.98,0);
  pad2 = new TPad("pad2","This is pad2",0.02,0.01,0.98,0.29,0);
  
  pad1->SetLogx();
  pad1->SetLogy();
  pad1->SetBottomMargin(0.01);
  pad2->SetLogx();
  //pad2->SetBottomMargin(0.13);
  pad2->SetTopMargin   (0.01);


  pad1->Draw(); // Projections pad
  pad2->Draw(); // Residuals   pad

  pad1->cd();

  hNum -> SetMaximum(10E-2);
  hNum -> SetMinimum(10E-7);

  char yAxisName[200];
  sprintf(yAxisName,"%s",gDen->GetYaxis()->GetTitle());
  hNum -> GetYaxis() -> SetTitle( yAxisName );
  hNum -> GetYaxis() -> SetTitleSize  (0.06);
  hNum -> GetYaxis() -> SetTitleOffset(0.71);
  hNum -> GetYaxis() -> SetLabelSize  (0.045);
  
  //int lineColor = gDen -> GetLineColor();
  //hDen->SetLineColor(lineColor);
  //hDen->SetLineWidth (2);
  //hDen->SetTitle("");
  //hDen->Draw("histo ");
  
  //   Double_t ptScale[19] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  // TH1F* hDummy = new TH1F("hDummy","",18,ptScale);
  // hDummy -> Draw();
  //hDen -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} [GeV/c]^{-1}");
  hDen -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma^{Z}}{dp_{T}^{Z}} [GeV]^{-1}");
  hDen -> GetYaxis() -> SetLabelSize(0.045);
  gNum->Draw("P ");
  //gDen->Draw("C same");

//   const int nptBins = 18;
//   Double_t xbins_pt[nptBins+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  
//   TH1F* hDummyDen = new TH1F("hDummyDen","", nptBins, xbins_pt);
   
//   for(int k=1; k<=hDummyDen->GetNbinsX(); k++){
//     double x,y;
//     gDen->GetPoint(k-1,x,y);
//     hDummyDen->SetBinContent(k,y);
//     hDummyDen->SetBinError  (k-1,gDen->GetErrorY(k));     
//    }
  
//   int lineColor = gDen -> GetLineColor();
//   hDummyDen  -> SetLineColor (lineColor);
//   //  hDummyDen  -> SetFillColor (kGreen+2);
//   hDummyDen  -> SetLineWidth (2);
//   hDummyDen  ->Draw("histo same");

  

  PrintIt(pad1,cTitle);

  if (leg) leg ->Draw("same");
  if (tex) tex ->Draw("same");
  //if (band) band ->Draw("F same");
  //if (line) line ->Draw("L same");
  if (tex1) tex1 ->Draw("same");
  //----------------------------------------------------------------------------
  // Residuals pad
  //----------------------------------------------------------------------------
  pad2->cd();

  TAxis *xPull = NULL;
  TAxis *yPull = NULL;
  char xAxisName[200];
  sprintf(xAxisName,"%s",gDen->GetXaxis()->GetTitle());
  axis1F(hPull,xPull,yPull,xAxisName,"(data-MC)/#sigma_{data}");
  //  axis1F(hPull,xPull,yPull,xAxisName,"data/MC");
   

  //if (hPull->GetMaximum() > 100) {
  //  hPull->SetMinimum(-100);
  //  hPull->SetMaximum( 100);
  //}

  hPull->GetXaxis()->SetLabelOffset(0.005);
  hPull->GetXaxis()->SetLabelSize  (0.11);
  //hPull->GetXaxis()->CenterTitle(1);
  hPull->GetXaxis()->SetTitleOffset(1.10);
  hPull->GetXaxis()->SetTitleSize  (0.12);
  hPull->GetXaxis()->SetNdivisions(7);
  
  hPull->GetYaxis()->SetLabelSize  (0.11);
  hPull->GetYaxis()->CenterTitle(1);
  hPull->GetYaxis()->SetTitleOffset(0.38);
  hPull->GetYaxis()->SetTitleSize  (0.10);

  //hPull->SetMaximum(6);
  //hPull->SetMinimum(-6);

  hPull->SetMarkerSize(0.8);
  //   TH1F* hDummy2 = new TH1F("hDummy2","",18,ptScale);

  //  hDummy2->SetMaximum(6);
  // hDummy2->SetMinimum(-6);
  // hDummy -> Draw();
 
  hPull->GetXaxis()->SetRange(0.,600.);
  //hPull->Draw("PL");
  //hPull->Draw("axis");

  // draw green band
  // draw yellow band

  const double x[4]={0,600,600,0};
  const double y1[4]={1,1,-1,-1};
  const double y2[4]={2,2,-2,-2};
  TGraph* band1s=new TGraph(4,x,y1);
  TGraph* band2s=new TGraph(4,x,y2);
  band1s->SetFillColor(kGreen-7);
  band2s->SetFillColor(kYellow);
  //band2s->Draw("F SAME");
  //band1s->Draw("F SAME");

  //   double Lxmax = canvas->GetUxmax();
  TLine* lineDotted = new TLine(0,0,440,0);
  lineDotted -> SetLineColor( kRed );
  lineDotted -> SetLineWidth( 2 );
  lineDotted -> SetLineStyle( 2 );
  //lineDotted -> Draw("L same");
  

  pad2->Update();
  //  pad2->GetFrame()->Draw("same");
     
  //  pad2->GetFrame()->DrawClone();
  //hPull->Draw("PL same");
  //  line -> Draw("L same");
}
*/


////===============================================================================================================

//------------------------------------------------------------------------------
// Draw projections and residuals
//------------------------------------------------------------------------------
/*
void DrawWithRatioTotal(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, 
		   TGraphAsymmErrors* gDen, 
		   TGraphAsymmErrors* gDen2, 
		   TGraphAsymmErrors* gDen3, 
		   TGraph* band, 
		   TGraph* band2, 
		   TGraph* band3, 
		   TLine* line,
		   TLine* line2,
		   TLine* line3,
		   TLegend* leg, 
		   TLegend* leg2, 
		   TLegend* leg3, 
		   TLatex* tex, TLatex* tex1, TLatex* tex2) {
  */
void DrawWithRatioTotal(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, 
		   TGraphAsymmErrors* gDen, 
		   TGraphAsymmErrors* gDen2, 
		   TGraphAsymmErrors* gDen3, 
		   TGraphAsymmErrors* gDen4, 
		   TGraph* band, 
		   TGraph* band2, 
		   TGraph* band3, 
		   TLegend* leg, 
		   TLegend* leg2, 
		   TLegend* leg3, 
		   TLatex* tex) {
  
  // sanity check
  if (gNum->GetN() != gDen->GetN()){
    std::cout<< " *** Error: number of points not consistent between the"
             << " two graphs -> Exit!\n";
    return;
  }

  TH1F *hNum = ConvertGraphToHisto(gNum);
  
  TH1F *hDen = ConvertGraphToHisto(gDen);
  TH1F *hDen2 = ConvertGraphToHisto(gDen2);
  TH1F *hDen3 = ConvertGraphToHisto(gDen3);
  TH1F *hDen4 = ConvertGraphToHisto(gDen4);

    std::cout << "hDen->GetNbinsX() = " << hDen->GetNbinsX() << std::endl;
  for (int i=1; i<hDen->GetNbinsX()+1; i++)
    std::cout << i << ") content=\t Powheg\t" << hDen->GetBinContent(i) <<"\t ResBos \t "<< hDen2->GetBinContent(i)<< "\t Madgraph \t "<< hDen3->GetBinContent(i)<<std::endl;

  //hNum -> GetXaxis() -> SetLimits(-100,10000);
  //hNum -> Sumw2();
  //hDen -> Sumw2();

  TH1F *hPull_data = (TH1F*) hNum -> Clone("hPull_data");
  TH1F *hPull = (TH1F*) hNum -> Clone("hPull");
  TH1F *hPull2 = (TH1F*) hNum -> Clone("hPull2");
  TH1F *hPull3 = (TH1F*) hNum -> Clone("hPull3");
  TH1F *hPull4 = (TH1F*) hNum -> Clone("hPull4");
  
  hPull->Sumw2();
  hPull2->Sumw2();
  hPull3->Sumw2();
  hPull4->Sumw2();
  
  hPull->Reset();
  hPull2->Reset();
  hPull3->Reset();
  hPull4->Reset();

  // Print the ratio
  std::cout << std::endl << std::endl;
  double bin[19] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
  double error_xlow[18];
  double error_xhigh[18];
  
  double data_errhigh[18] = {0,};
  double data_errlow[18] = {0,};
  double powheg_err[18] = {0,};
  double resbos_err[18] = {0,};
  double madgraph_err[18] = {0,};
  double fewz_errhigh[18] = {0,};
  double fewz_errlow [18] = {0,};
    
  double ratioData[18] = {1,};
  double ratioPowheg[18] = {0,};
  double ratioResbos[18] = {0,};
  double ratioMadgraph[18] = {0,};
  double ratioFEWZ[18] = {0,};
   
  double data_errbandhigh[18] = {0,};
  double data_errbandlow[18] = {0,};
  double powheg_errband[18] = {0,};
  double resbos_errband[18] = {0,};
  double madgraph_errband[18] = {0,};
  double fewz_errbandhigh[18] = {0,};
  double fewz_errbandlow [18] = {0,};
    
  for (int i=0;i<18;i++) {
    error_xlow[i] = (point[i] - bin[i]);
    error_xhigh[i] = (bin[i+1]-point[i]);
  }

  for (int iBin=0;iBin<gNum->GetN();iBin++) {
    
    Double_t x,yNum,yDen,yDen2,yDen3,yDen4, error_high, error_low, error;
    
    gNum->GetPoint(iBin,x,yNum); //data
    
    gDen->GetPoint(iBin,x,yDen); // powheg
    gDen2->GetPoint(iBin,x,yDen2); // resbos
    gDen3->GetPoint(iBin,x,yDen3); // madgraph
    gDen4->GetPoint(iBin,x,yDen4); // fewz
    
    data_errhigh[iBin] = gNum->GetErrorYhigh(iBin);
    data_errlow[iBin] = gNum->GetErrorYlow(iBin);
    powheg_err[iBin] = gDen->GetErrorYhigh(iBin);
    resbos_err[iBin] = gDen2->GetErrorYhigh(iBin);
    madgraph_err[iBin] = gDen3->GetErrorYhigh(iBin);
    fewz_errhigh[iBin] = gDen4->GetErrorYhigh(iBin);
    fewz_errlow[iBin] = gDen4->GetErrorYlow(iBin);

    data_errbandhigh[iBin] = data_errhigh[iBin] / yNum ;
    data_errbandlow[iBin] = data_errlow[iBin] / yNum;
    powheg_errband[iBin] = powheg_err[iBin] / yNum; 
    resbos_errband[iBin] = resbos_err[iBin] / yNum; 
    madgraph_errband[iBin] = madgraph_err[iBin] / yNum;
    fewz_errbandhigh[iBin] = fewz_errhigh[iBin] / yNum;
    fewz_errbandlow[iBin] = fewz_errlow[iBin] / yNum;

    cout << "resbos error : " << resbos_err[iBin] << endl;
    //double data_err = gNum->GetErrorYhigh(iBin); // data error
    //double powheg_err = gDen->GetErrorYhigh(iBin); // powheg error
    //double resbos_err = gDen2->GetErrorYhigh(iBin); // resbos error
    //double madgraph_err = gDen3->GetErrorYhigh(iBin); // madgraph error
    //double fewz_errhigh = gDen4->GetErrorYhigh(iBin); // fewz error high
    //double fewz_errlow = gDen4->GetErrorYlow(iBin); // fewz error low

    // double sigmaN, sigmaD;
    //error_highN = gNum->GetErrorYhigh(iBin);
    //error_lowN  = gNum->GetErrorYlow(iBin);
    //if (error_highN > error_lowN) sigmaN = error_highN;
    //else sigmaN = error_lowN;

    //error_highD = gDen->GetErrorYhigh(iBin);
    //error_lowD  = gDen->GetErrorYlow(iBin);
    //if (error_highD > error_lowD) sigmaD = error_highD;
    //else sigmaD = error_lowD;
    
   //--this is for the final PAS--code start for below plot

    error_high = gNum->GetErrorYhigh(iBin);
    error_low  = gNum->GetErrorYlow(iBin);
    if (error_high > error_low) error = error_high;
    else error = error_low;
   
    // (Theory - data) / sigmaData style
    //double resid = (error) ? ((yNum-yDen)/error) : 0.0;
    //double resid2 = (error) ? ((yNum-yDen2)/error) : 0.0;
    //double resid3 = (error) ? ((yNum-yDen3)/error) : 0.0;
    //
    //if (fabs(resid) > 100) resid = 10;
    //if (fabs(resid2) > 100) resid2 = 10;
    //if (fabs(resid3) > 100) resid3 = 10;

    //hPull->SetBinContent(iBin+1,resid);
    //hPull2->SetBinContent(iBin+1,resid2);
    //hPull3->SetBinContent(iBin+1,resid3);
    //cout << iBin << "\t resid\t " << resid <<"\t resid2\t " << resid2 <<"\t resid3\t " << resid3 << endl;
   
    // Theory / Data style
    //double ratioData = 1;
    //double ratioPowheg = yDen / yNum;
    //double ratioResbos = yDen2 / yNum;
    //double ratioMadgraph = yDen3 / yNum;
    //double ratioFEWZ = yDen4 / yNum;

    ratioData[iBin] = 1;
    ratioPowheg[iBin] = yDen / yNum;
    ratioResbos[iBin] = yDen2 / yNum;
    ratioMadgraph[iBin] = yDen3 / yNum;
    ratioFEWZ[iBin] = yDen4 / yNum;

    hPull_data->SetBinContent(iBin+1, ratioData[iBin]);
    hPull_data->SetBinError(iBin+1,data_errhigh[iBin]/yNum);
   // 
   // hPull->SetBinContent(iBin+1,ratioPowheg);
   // hPull->SetBinError(iBin+1,powheg_err/yNum);

   // hPull2->SetBinContent(iBin+1,ratioResbos);
   // hPull2->SetBinError(iBin+1,resbos_err/yNum);
   // 
   // hPull3->SetBinContent(iBin+1,ratioMadgraph);
   // hPull3->SetBinError(iBin+1,madgraph_err/yNum);
   // 
   // hPull4->SetBinContent(iBin+1,ratioFEWZ);
    //hPull4->SetBinError(iBin+1,fewz_errhigh/yNum);
    //cout << iBin << "\t ratioPowheg : " << ratioPowheg <<"\t Resbos : " << ratioResbos <<"\t Madgraph : " << ratioMadgraph << endl;
    //code end
    
    hPull->SetBinContent(iBin+1,1);

    hPull2->SetBinContent(iBin+1,1);
    
    hPull3->SetBinContent(iBin+1,1);
    
    hPull4->SetBinContent(iBin+1,1);
    //hPull4->SetBinError(iBin+1,fewz_errhigh/yNum);
    //cout << iBin << "\t ratioPowheg : " << ratioPowheg <<"\t Resbos : " << ratioResbos <<"\t Madgraph : " << ratioMadgraph << endl;
    //code end
  
    // double ratio = yNum/yDen;

    //double sigmaR = sqrt(((1/yNum)*(1/yNum))*(sigmaN*sigmaN)+ ((yNum/(yDen*yDen))*(yNum/(yDen*yDen)))*(sigmaD*sigmaD));
    //  std::cout << iBin << ") " << yNum << "/" << yDen << ", Ratio= " << ratio << std::endl;
    // hPull->SetBinContent(iBin+1,ratio);
    //hPull->SetBinError(iBin+1,sigmaR);
   
    cout << "hDen Y : " << hDen->GetBinContent(iBin+1) << "\t error : " << hDen->GetBinError(iBin+1) << endl;
  }
  std::cout << std::endl << std::endl;

  //----------------------------------------------------------------------------
  // Create the pads
  //----------------------------------------------------------------------------
  TPad* pad1;
  TPad* pad2;

  pad1 = new TPad("pad1","This is pad1",0.00,0.42,1.0,0.98);
  pad2 = new TPad("pad2","This is pad2",0.00,0.01,1.0,0.49);
  
  pad1->SetLogx();
  pad1->SetLogy();
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetTopMargin(0.1);
  
  pad2->SetLogx();
  pad2->SetTickx(1);
  pad2->SetTicky(1);

  pad1->Draw(); // Projections pad
  pad2->Draw(); // Residuals   pad

  pad1->cd();
  
  hNum -> SetMaximum(10E-2);
  hNum -> SetMinimum(10E-7);

  char yAxisName[200];
  sprintf(yAxisName,"%s",gDen->GetYaxis()->GetTitle());
  hNum -> GetYaxis() -> SetTitle( yAxisName );
  hNum -> GetYaxis() -> SetTitleSize  (0.06);
  hNum -> GetYaxis() -> SetTitleOffset(0.71);
  hNum -> GetYaxis() -> SetLabelSize  (0.045);
  
/* 
  int lineColor = gDen -> GetLineColor();
  hDen->SetLineColor(lineColor);
  hDen->SetLineWidth (2);
  hDen->SetTitle("");
  //hDen->Draw("histo ");
  
  int lineColor2 = gDen2 -> GetLineColor();
  hDen2->SetLineColor(lineColor2);
  hDen2->SetLineWidth (2);
  hDen2->SetTitle("");
  //hDen2->Draw("histo same");
  
  int lineColor3= gDen3-> GetLineColor();
  hDen3->SetLineColor(lineColor3);
  hDen3->SetLineWidth (2);
  hDen3->SetTitle("");
  //hDen3->Draw("histo same");
*/

  //TH1F* hDen_dummy = (TH1F*)hDen->Clone();
  TH1F* hDen_dummy = new TH1F("","",18,0.5,600);
  hDen_dummy->Reset();
  hDen_dummy->SetMaximum(10E-2);
  hDen_dummy->SetMinimum(5*10E-8);
  hDen_dummy->SetTitle("");
  //hDen_dummy -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{Z}} [GeV]^{-1}");
  hDen_dummy -> GetYaxis() -> SetTitleSize(0.08);
  hDen_dummy -> GetYaxis() -> SetTitleOffset(0.7);
  hDen_dummy -> GetYaxis() -> SetTitle("1/#sigma d#sigma/dp_{T}^{Z} [GeV]^{-1}");
  //hDen_dummy -> GetYaxis() -> SetLabelSize(0.045);
  hDen_dummy -> GetYaxis() -> SetLabelSize(0.05);
  hDen_dummy -> GetXaxis() -> SetLabelSize(0);
  hDen_dummy -> GetXaxis() -> SetRangeUser(0.5,600); // here control Xrange
  hDen_dummy->Draw();
  
  //grMadgraphDist->Draw("P");

  //grMadgraphDist->Draw("5");
  //grResbosDist->Draw("5");
  //grPowhegDist->Draw("5");
  //grFEWZDist->Draw("5");
  
  //   Double_t ptScale[19] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  // TH1F* hDummy = new TH1F("hDummy","",18,ptScale);
  // hDummy -> Draw();
  //hDen -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} [GeV/c]^{-1}");
  hDen -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{Z}} [GeV]^{-1}");
  hDen -> GetYaxis() -> SetTitleSize(0.30);
  hDen -> GetYaxis() -> SetLabelSize(0.055);
  hDen -> GetYaxis() -> SetTitleOffset(1.2);
  //gDen->Draw("C same");

  //Color Transparent
  TColor *colRed = gROOT->GetColor(kRed);
  TColor *colBlue = gROOT->GetColor(kBlue);
  TColor *colViolet = gROOT->GetColor(kViolet);
  TColor *colGreen = gROOT->GetColor(kGreen);
  colRed->SetAlpha(0.2);
  colBlue->SetAlpha(0.2);
  colViolet->SetAlpha(0.2);
  colGreen->SetAlpha(0.2);

  //// Now design and Draw
  gStyle->SetLineWidth(2.);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesSpacing(0.75);
  gStyle->SetHatchesLineWidth(2);

// Powheg Draw   
  //gDen->SetFillStyle(3004);
  //gDen->SetFillColor(kGreen+2);
  //gDen->SetLineColor(kGreen+2);
  gDen->SetFillColor(kRed);
  gDen->SetLineColor(kRed+2);

// Resbos Draw   
  //gDen2->SetFillStyle(3005);
  //gDen2->SetFillColor(kBlue);
  //gDen2->SetLineColor(kBlue);
  gDen2->SetFillColor(kBlue);
  gDen2->SetLineColor(kBlue+2);
  
// Madgraph Draw   
  //gDen3->SetFillStyle(3013);
  //gDen3->SetFillColor(kViolet);
  //gDen3->SetLineColor(kViolet);
  gDen3->SetFillColor(kViolet);
  gDen3->SetLineColor(kViolet+2);
  
// FEWZ Draw   
  //gDen4->SetFillStyle(3002);
  //gDen4->SetFillColor(kGreen);
  //gDen4->SetLineColor(kGreen);
  gDen4->SetFillColor(kGreen);
  gDen4->SetLineColor(kGreen+2);
  
  gDen->Draw("5"); // Powheg Distribution
  gDen2->Draw("5"); // Resbos dist
  //gDen3->Draw("5"); // Madgraph dist
  gDen4->Draw("5"); // FEWZ dist
  gNum->Draw("P ");  // data dist

  cout << "gNum print" << endl;
  gNum->Print();

  cout << "gDen print" << endl;
  gDen->Print();

  CMS_lumi(pad1,iPeriod,iPos);
  pad1->RedrawAxis();







//   const int nptBins = 18;
//   Double_t xbins_pt[nptBins+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  
//   TH1F* hDummyDen = new TH1F("hDummyDen","", nptBins, xbins_pt);
   
//   for(int k=1; k<=hDummyDen->GetNbinsX(); k++){
//     double x,y;
//     gDen->GetPoint(k-1,x,y);
//     hDummyDen->SetBinContent(k,y);
//     hDummyDen->SetBinError  (k-1,gDen->GetErrorY(k));     
//    }
  
//   int lineColor = gDen -> GetLineColor();
//   hDummyDen  -> SetLineColor (lineColor);
//   //  hDummyDen  -> SetFillColor (kGreen+2);
//   hDummyDen  -> SetLineWidth (2);
//   hDummyDen  ->Draw("histo same");

  

  PrintIt(pad1,cTitle);

  TLegend* tlTotal = new TLegend(0.2,0.2,0.4,0.5);
  tlTotal -> SetTextFont(42);
  tlTotal -> SetFillColor(0);
  tlTotal -> SetBorderSize(0);
  tlTotal -> AddEntry(gDen,"Pythia/Powheg","f");
  tlTotal -> AddEntry(gDen2,"Resbos","f");
  //tlTotal -> AddEntry(gDen3,"Madgraph","f");
  tlTotal -> AddEntry(gDen4,"FEWZ","f");
  tlTotal -> AddEntry(gNum,"data","lp");
  tlTotal->Draw("same");
  
  //if (leg) leg ->Draw("same");
  //if (leg2) leg2 ->Draw("same");
  //if (leg3) leg3 ->Draw("same");
  
  if (tex) tex ->Draw("same");
  
//  if (band) band ->Draw("F same");
//  if (band2) band2 ->Draw("F same");
//  if (band3) band3 ->Draw("F same");
//  
  //if (line) line ->Draw("L same");
  //if (line2) line2 ->Draw("L same");
  //if (line3) line3 ->Draw("L same");
  
  
  //----------------------------------------------------------------------------
  // Residuals pad
  //----------------------------------------------------------------------------
  pad2->cd();

  TGraphAsymmErrors* grRatioData = new TGraphAsymmErrors(18,point,ratioData,error_xlow,error_xhigh,data_errbandlow,data_errbandhigh);
  TGraphAsymmErrors* grRatioPowheg = new TGraphAsymmErrors(18,point,ratioPowheg,error_xlow,error_xhigh,powheg_errband,powheg_errband);
  TGraphAsymmErrors* grRatioResbos = new TGraphAsymmErrors(18,point,ratioResbos,error_xlow,error_xhigh,resbos_errband,resbos_errband);
  TGraphAsymmErrors* grRatioMadgraph = new TGraphAsymmErrors(18,point,ratioMadgraph,error_xlow,error_xhigh,madgraph_errband,madgraph_errband);
  TGraphAsymmErrors* grRatioFEWZ = new TGraphAsymmErrors(18,point,ratioFEWZ,error_xlow,error_xhigh,fewz_errbandlow,fewz_errbandhigh);

  TAxis *xPull = NULL;
  TAxis *yPull = NULL;
  char xAxisName[200];
  sprintf(xAxisName,"%s",gDen->GetXaxis()->GetTitle());
  //axis1F(hPull,xPull,yPull,xAxisName,"(data-MC)/#sigma_{data}");
  axis1F(hPull,xPull,yPull,xAxisName,"Theory / Data");
  //  axis1F(hPull,xPull,yPull,xAxisName,"data/MC");
   

  //if (hPull->GetMaximum() > 100) {
  //  hPull->SetMinimum(-100);
  //  hPull->SetMaximum( 100);
  //}
  hPull->GetXaxis()->SetLabelOffset(0.005);
  //hPull->GetXaxis()->SetLabelSize  (0.11);
  hPull->GetXaxis()->SetLabelSize  (0.04);
  hPull->GetXaxis()->SetTitleOffset(1.10);
  //hPull->GetXaxis()->SetTitleSize  (0.12);
  hPull->GetXaxis()->SetTitleSize  (0.05);
  hPull->GetXaxis()->SetNdivisions(7);
  
  //hPull->GetYaxis()->SetLabelSize  (0.11);
  hPull->GetYaxis()->SetLabelSize  (0.04);
  hPull->GetYaxis()->CenterTitle(1);
  hPull->GetYaxis()->SetTitleOffset(0.38);
  //hPull->GetYaxis()->SetTitleSize  (0.10);
  hPull->GetYaxis()->SetTitleSize  (0.52);

  hPull->SetMaximum(2.7);
  hPull->SetMinimum(0.2);
  //hPull->SetMaximum(6.5);
  //hPull->SetMinimum(-1.5);

  hPull->SetMarkerSize(0.8);
  hPull->SetMarkerColor(kGreen+2);
  hPull->SetLineColor(kGreen+2);
  
  hPull2->SetMarkerSize(0.8);
  hPull2->SetMarkerColor(kBlue);
  hPull2->SetLineColor(kBlue);
  
  hPull3->SetMarkerSize(0.8);
  hPull3->SetMarkerColor(kViolet);
  hPull3->SetLineColor(kViolet);
 
  //TH1F *hPull_dummy = (TH1F*)hPull->Clone();
  TH1F *hPull_dummy = new TH1F("","",18,0.5,600);
  hPull_dummy->Reset();
  hPull_dummy->GetXaxis()->SetLabelOffset(0.005);
  hPull_dummy->GetXaxis()->SetLabelSize  (0.06);
  hPull_dummy->GetXaxis()->SetTitleOffset(0.50);
  hPull_dummy->GetXaxis()->SetTitleSize  (0.07);
  hPull_dummy->GetXaxis()->SetNdivisions(7);
  hPull_dummy->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  hPull_dummy->GetYaxis()->SetLabelSize  (0.06);
  hPull_dummy->GetYaxis()->CenterTitle(1);
  hPull_dummy->GetYaxis()->SetTitle("Theory / Data");
  hPull_dummy->GetYaxis()->SetTitleOffset(0.78);
  hPull_dummy->GetYaxis()->SetTitleSize  (0.08);
  hPull_dummy->SetMaximum(2.9);
  hPull_dummy->SetMinimum(0.2);
  hPull_dummy->GetXaxis()->SetRangeUser(0.5,600);
  hPull_dummy->Draw();

 // Data Ratio draw 
  grRatioData->SetFillStyle(3354);
  grRatioData->SetFillColor(kGray+2);
  grRatioData->SetLineColor(kBlack);

// Powheg Ratio draw 
  grRatioPowheg->SetMarkerStyle(22); // triangle
  grRatioPowheg->SetMarkerColor(kRed+2);
  grRatioPowheg->SetFillColor(kRed);
  grRatioPowheg->SetLineColor(kRed+2);
 
  // Resbos Ratio draw
  grRatioResbos->SetMarkerStyle(20); // circle
  grRatioResbos->SetMarkerColor(kBlue+2);
  grRatioResbos->SetFillColor(kBlue);
  grRatioResbos->SetLineColor(kBlue+2);
 
// Madgraph Ratio draw 
  grRatioMadgraph->SetMarkerStyle(34); // cross
  grRatioMadgraph->SetMarkerColor(kViolet+2);
  grRatioMadgraph->SetFillColor(kViolet);
  grRatioMadgraph->SetLineColor(kViolet+2);
 
// FEWZ Ratio draw 
  grRatioFEWZ->SetMarkerStyle(21); // squre
  grRatioFEWZ->SetMarkerColor(kGreen+2);
  grRatioFEWZ->SetFillColor(kGreen);
  grRatioFEWZ->SetLineColor(kGreen+2);
  
  grRatioData->Draw("2");
  grRatioResbos->Draw("5 P");
  grRatioPowheg->Draw("5 P");
  //grRatioMadgraph->Draw("5 P");
  grRatioFEWZ->Draw("5 P");

/*
  // draw green band
  // draw yellow band

  const double x[4]={0,600,600,0};
  const double y1[4]={1,1,-1,-1};
  const double y2[4]={2,2,-2,-2};
  TGraph* band1s=new TGraph(4,x,y1);
  TGraph* band2s=new TGraph(4,x,y2);
  band1s->SetFillColor(kGreen-7);
  band2s->SetFillColor(kYellow);
  band2s->Draw("F SAME");
  band1s->Draw("F SAME");

  //   double Lxmax = canvas->GetUxmax();
  TLine* lineDotted = new TLine(0,0,440,0);
  lineDotted -> SetLineColor( kRed );
  lineDotted -> SetLineWidth( 2 );
  lineDotted -> SetLineStyle( 2 );
  lineDotted -> Draw("L same");
  

  pad2->Update();
  //  pad2->GetFrame()->Draw("same");
     
  //  pad2->GetFrame()->DrawClone();
  hPull2->Draw("PL same");
  hPull3->Draw("PL same");
  //  line -> Draw("L same");

*/
}
////===============================================================================================================








void axis1F(TH1F  *histo,
            TAxis *xaxis,
            TAxis *yaxis,
            char  *xtitle,
            char  *ytitle)
{
  histo->SetMarkerSize(0.5);
  histo->SetMarkerStyle(kFullCircle);

  histo->SetTitle("");

  xaxis = histo->GetXaxis();
  yaxis = histo->GetYaxis();

  yaxis->SetTitle(ytitle);
  yaxis->SetLabelFont(42);
  yaxis->SetLabelOffset(0.01);
  yaxis->SetLabelSize(0.05);
  yaxis->SetNdivisions(505);
  yaxis->SetTitleColor(kBlack);
  yaxis->SetTitleFont(42);
  yaxis->SetTitleOffset(1.5);
  yaxis->SetTitleSize(0.17);
  yaxis->CenterTitle(kTRUE);

  xaxis->SetNdivisions(505);
  xaxis->SetLabelFont(42);
  xaxis->SetLabelOffset(0.005);
  xaxis->SetTitle(xtitle);
  xaxis->SetTitleColor(kBlack);
  xaxis->SetTitleFont(42);
  xaxis->SetTitleOffset(1.0);
  xaxis->SetTitleSize(0.045);
}



TH1F* ConvertGraphToHisto(TGraph *pGraph){
  // takes data from a graph, determines binning and fills data into histogram
  Int_t NPoints = pGraph->GetN();
  std::cout << "NPoints=" << NPoints << std::endl;
  Double_t BinLimits[NPoints+1];
  // sort graph
  pGraph->Sort();
  // determine lower limit of histogram: half the distance to next point
  Double_t x0,x1,y;
  pGraph->GetPoint(0,x0,y);
  pGraph->GetPoint(1,x1,y);
  Double_t Distance = TMath::Abs(x0-x1);
  BinLimits[0] = x0 - Distance/2.;
  
  //std::cout << "x0=" << x0 << ", x1=" << x1 
  //          << ", Distance=" << Distance
  //          << ", BinLimits[0]=" << BinLimits[0] << std::endl; 
  
  // now set upper limits for all the other points
  for (Int_t k = 0 ; k<NPoints;k++){
    pGraph->GetPoint(k  ,x0,y);
    pGraph->GetPoint(k+1,x1,y);
    
    Distance = TMath::Abs(x0-BinLimits[k]);
    BinLimits[k+1] = x0 + Distance;
    
    //std::cout << "x0=" << x0 << ", x1=" << x1
    //          << ", Distance=" << Distance
    //          << ", BinLimits[" << k+1 <<"]=" << BinLimits[k+1] << std::endl; 
  }

  //for (int i=0;i<NPoints+1;i++)
  //std::cout << "BinLimits[" << i <<"]=" << BinLimits[i] << std::endl;

  // now we know the binning and can create the histogram:
  TString Name = "ConvertedHisto"; 
  // make name unique 
  Name+= rand();
  TH1F *ThisHist = new TH1F(Name,"Converted Histogram",NPoints,BinLimits);
  // now fill the histogram
  for (Int_t i = 0; i<pGraph->GetN();i++){
    Double_t x,y,err;
    pGraph->GetPoint(i,x,y);
    err = pGraph->GetErrorYhigh(i);
    ThisHist->SetBinContent(i+1,y);
    ThisHist->SetBinError(i+1,err);
    
    std::cout << i << ")  x=" << x
	      << ", y=" << y << std::endl; 
  }

  return ThisHist;
}

void PrintIt(TPad *pad, TString title)//char *title)
{
  TLatex *latex = new TLatex();
  latex->SetTextFont(  42);
  latex->SetTextSize(0.05);
  
  // Get the most recent changes
  pad->Update();
  

  double xmin = pad->GetUxmin();
  double xmax = pad->GetUxmax();
  double ymin = pad->GetUymin();
  double ymax = pad->GetUymax();
 
  double xpos = xmin + 0.50*(xmax - xmin);
  double ypos = ymax + 0.05*(ymax - ymin);

  if (pad->GetLogy())  ypos = pow(10,ypos);
  if (pad->GetLogx())  xpos = pow(10,xpos);

  latex->SetTextAlign(22);
  latex->DrawLatex(xpos,ypos,title);
}


