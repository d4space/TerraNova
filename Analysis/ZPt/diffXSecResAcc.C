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

#include "selections.h"
#include "_unfoldData.h"
#include "systematics.h"

TH1F* ConvertGraphToHisto(TGraph *pGraph);

void PrintIt(TPad *pad, TString title);

void fillHistos(TTree* tree,
                TH1F *hPt,
                double mLow,
                double mHigh,
                int mode);

TGraphAsymmErrors* MuonChannel2010();
TGraphAsymmErrors* Combination2010();
TGraphAsymmErrors* PowhegRatio();
TGraphAsymmErrors* FEWZRatio();


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

 
 
void DrawWithRatio(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, 
                   TGraphAsymmErrors* gDen,
                   TLegend* leg=NULL);

void axis1F(TH1F  *histo,
            TAxis *xaxis,
            TAxis *yaxis,
            char  *xtitle,
            char  *ytitle);


void diffXSecResAcc(double intLumi_invpb  =   18.424,
                    int    mBins          =   30,  // Z mass
                    double mLow           =   60,  // Z mass
                    double mHigh          =  120,  // Z mass
                    bool   isSave         = !true) {
  
  gROOT->Clear();
  gStyle->SetOptStat(0);  
 
  // ---------------------------------------------------------------------------
  // general variables
  char cTitle[200];
  sprintf(cTitle, "CMS Preliminary");

  char cTitleExtd[200];
  sprintf(cTitleExtd, "CMS Preliminary, 18.4 pb^{-1} at #sqrt{s} = 8 TeV");

  //TString png      = "png/diffXSec/";
  //TString rootPlot = "rootPlot/diffXSec/";
  // ---------------------------------------------------------------------------


  const int nptBins = 18;
  Double_t xbins_pt[nptBins+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};

  // this is useful only if we stick the data in center of the bins
  Double_t error_xbins[nptBins];
  for (int i=0;i<nptBins;i++) error_xbins[i] = (xbins_pt[i+1] - xbins_pt[i]) / 2;
  

  
  TH1F* hXSec          = new TH1F("hXSec"      , "", nptBins, xbins_pt);
  TH1F* hXSecFSR       = new TH1F("hXSecFSR"      , "", nptBins, xbins_pt);

  // Fill the data
  //  TFile* fileData = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_5_ecalpatch1/Data/ntuplesWITHdz/SingleMuRun2012A_190949_191090_191367_193112_193116.root");
  //TFile* fileData = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_5_ecalpatch1/Data/ntuplesWITHdz/fullselection/data2012A_LowPU.root");
  TFile* fileData = new TFile("data2012A_LowPU.root");
  TTree* treeData = (TTree*) fileData -> Get("tree");

  // Fill the TTbarJets MC  
  // TFile* fileMC_ttbarjets = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/TTJets/mergeFile/TTJets_merged.root");
  //TFile* fileMC_ttbarjets = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/TTJets.root");
  TFile* fileMC_ttbarjets = new TFile("TTJets.root");
  TTree* treeMC_ttbarjets = (TTree*) fileMC_ttbarjets -> Get("tree");
 
  // Fill the ZTauTau MC  
  //  TFile* fileMC_ztautau = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/DYToTauTau/mergeFile/DYToTauTau_merged.root");
  //TFile* fileMC_ztautau = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/DYToTauTau.root");
  TFile* fileMC_ztautau = new TFile("DYToTauTau.root");
  TTree* treeMC_ztautau = (TTree*) fileMC_ztautau -> Get("tree");
 
  // Fill the WZ MC  
  //  TFile* fileMC_WZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WZ/mergeFile/WZ_merged.root");
  //TFile* fileMC_WZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/WZ.root");
  TFile* fileMC_WZ = new TFile("WZ.root");
  TTree* treeMC_WZ = (TTree*) fileMC_WZ -> Get("tree");

  // Fill the ZZ MC  
  //  TFile* fileMC_ZZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/ZZ/mergeFile/ZZ_merged.root");
  //TFile* fileMC_ZZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/ZZ.root");
  TFile* fileMC_ZZ = new TFile("ZZ.root");
  TTree* treeMC_ZZ = (TTree*) fileMC_ZZ -> Get("tree");

 
  // Fill the WW MC  
  // TFile* fileMC_WW = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WW/mergeFile/WW_merged.root");
  //TFile* fileMC_WW = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/WW.root");
  TFile* fileMC_WW = new TFile("WW.root");
  TTree* treeMC_WW = (TTree*) fileMC_WW -> Get("tree");

  // Fill the WJetsToLNu MC  
  //  TFile* fileMC_WJetsToLNu = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WJetsToLNu/mergeFile/WJetsToLNu_merged.root");
  //TFile* fileMC_WJetsToLNu = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/WJetsToLNu.root");
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
  // Data 
  std::cout << " Data " <<std::endl;
  for(int bin=1; bin<=hPtData->GetNbinsX(); bin++)
  {
    std::cout << bin<<"\t"<<hPtData->GetBinContent(bin)<<"\t\t"<<hPtData->GetBinError(bin)<<std::endl;
  }
  std::cout << " Total Data \t" <<hPtData->Integral()<<std::endl;
  

  
  
  // 1. Subtract the bkg
  hPtData -> Sumw2();
  hPtData -> Add(hPtMC_allBkg, -1);
  
   std::cout << " Subtract the bkg " <<std::endl;
  for(int bin=1; bin<=hPtData->GetNbinsX(); bin++)
  {
    std::cout << hPtData->GetBinContent(bin)<<"\t\t"<<hPtData->GetBinError(bin)<<std::endl;
  }
  std::cout << " Total \t" <<hPtData->Integral()<<std::endl;
  
  
 
  
  
  // 2. Correct for  Eff and scale factors by using accEff_Lumi.C
  // For Nadeesha: probably you will need to change the location
  // TFile* fEffvsPt = new TFile("histo_EffvsPt.root");
     TFile* fEffvsPt = new TFile("histo_EffvsPt_withnewSF.root"); // this is the one with new SF from GP. 
  
  //  TH1F* hEffvsPt = (TH1F*) fEffvsPt -> GetListOfPrimitives()->At(1);
  TH1F *hEffvsPt = (TH1F*)fEffvsPt->Get("hEff");

  hEffvsPt->Sumw2();
  hPtData -> Divide(hEffvsPt);
  
  for(int bin=1; bin<=hEffvsPt->GetNbinsX(); bin++)
  {
    std::cout << "hEffvsPt \t\t\t"<<bin<<"\t"<<hEffvsPt->GetBinContent(bin)<<std::endl;
  }
  std::cout << " After Correct for  Eff hPtData " <<std::endl;
  for(int bin=1; bin<=hEffvsPt->GetNbinsX(); bin++)
  {
    std::cout << hPtData->GetBinContent(bin)<<"\t\t"<<hPtData->GetBinError(bin)<<std::endl;
  }
  std::cout << " Total \t" <<hPtData->Integral()<<std::endl;
    
  
  
  
  
  
  // 3. Unfolding Resolution
  hXSec -> Sumw2();
  applyMatrix(hXSec, hPtData);
  
  std::cout << " After Correct for  Unfolding Resolution " <<std::endl;
  for(int bin=1; bin<=hXSec->GetNbinsX(); bin++)
  {
    std::cout << hXSec->GetBinContent(bin)<<"\t\t"<<hXSec->GetBinError(bin)<<std::endl;
  }
  std::cout << " Total \t" <<hXSec->Integral()<<std::endl;

  
  
  
  // FSR correction
  hXSecFSR -> Sumw2();
  applyMatrixFSR(hXSecFSR, hXSec);
 
  std::cout << " After FSR correction " <<std::endl;
  for(int bin=1; bin<=hXSecFSR->GetNbinsX(); bin++)
  {
    std::cout << hXSecFSR->GetBinContent(bin)<<"\t\t"<<hXSecFSR->GetBinError(bin)<<std::endl;
  }
  std::cout << " Total \t" <<hXSecFSR->Integral()<<std::endl;
  
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
    //double value = hXSecFSR->GetBinContent(bin)/(hXSecFSR->GetBinWidth(bin));// * hXSecFSR -> Integral()); 
    //double sigma = hXSecFSR->GetBinError(bin)  /(hXSecFSR->GetBinWidth(bin));// * hXSecFSR -> Integral()); 
     
    //for systematics fractional error on diffXSec
    //   double FerrorXSec = sigma/value;
    //   printf("%8.8f \n",
    //      FerrorXSec);  
    //end systematics fractional error on diffXSec
   
    	   printf("%2d) %5.2f | %8.8f +/- %8.8f \n", 
          bin, hXSecFSR->GetXaxis()->GetBinCenter(bin), value, sigma);
 
    
    gXSec->SetPoint(bin-1,hXSecFSR->GetXaxis()->GetBinCenter(bin),value);
    // For Nadeesha: here I put the error on the x-axis 
    gXSec->SetPointError(bin-1,error_xbins[bin-1],error_xbins[bin-1],sigma,sigma);
  }


  

  // some style on the graph
  //  gXSec -> SetTitle("");
  // gXSec -> GetXaxis() -> SetTitle("q_{T} [GeV]");
  //gXSec -> GetYaxis() -> SetTitle("1/#sigma d#sigma/dp_{T} [GeV]^{-1}");
  //gXSec -> GetYaxis() -> SetTitleOffset(1.2);
  //gXSec -> SetMarkerSize(0.7);
  //gXSec -> SetMarkerStyle(20);


  //-----------------------------------------------------------------------------------------------
  TGraphAsymmErrors* grMuon = MuonChannel2010();
  TGraphAsymmErrors* grComb = Combination2010();
  TGraphAsymmErrors* grPowhegRatio = PowhegRatio();
  TGraphAsymmErrors* grFEWZRatio = FEWZRatio();

  // some style on the graphs
  grMuon -> SetTitle("");
  grMuon -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  grMuon -> GetYaxis() -> SetTitle("1/#sigma d#sigma/dq_{T} [GeV/c]^{-1}");
  grMuon -> GetYaxis() -> SetTitleOffset(1.02);
  grMuon -> SetMarkerSize(0.7);
  grMuon -> SetMarkerStyle(20);
  grMuon -> SetMarkerColor(kBlue);
  grMuon -> SetLineColor  (kBlue);

  grComb -> SetTitle("");
  grComb -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  grComb -> GetYaxis() -> SetTitle("1/#sigma d#sigma/dq_{T} [GeV/c]^{-1}");
  grComb -> GetYaxis() -> SetTitleOffset(1.02);
  grComb -> SetMarkerSize(0.7);
  grComb -> SetMarkerStyle(20);
  grComb -> SetMarkerColor(kRed);
  grComb -> SetLineColor  (kRed);

  //-----------------------------------------------------------------------------------------------
  TCanvas* c1 = new TCanvas("c1","",0,0,700,500);
  c1->SetFillColor(0);
  c1->SetLogx();
  c1->SetLogy();


  //  gXSec->Draw("AP");
  // PrintIt(c1, cTitleExtd);

  
  TGraphAsymmErrors* gXSecWithSyst= new TGraphAsymmErrors(hXSecFSR->GetNbinsX());
  fillSystErrorEff(gXSec,gXSecWithSyst);

  gXSecWithSyst -> SetTitle("");
  
  gXSecWithSyst -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  gXSecWithSyst -> GetXaxis() -> SetTitleSize(0.04);
  gXSecWithSyst -> GetXaxis() -> SetLabelSize(0.03);

  gXSecWithSyst -> GetYaxis() -> SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} [GeV/c]^{-1}");
  gXSecWithSyst -> GetYaxis() -> SetTitleSize(0.04);
  gXSecWithSyst -> GetYaxis() -> SetLabelSize(0.03);
  gXSecWithSyst -> GetYaxis() -> SetTitleOffset(1.02);

  gXSecWithSyst -> SetMarkerSize(0.8);
  gXSecWithSyst -> SetMarkerStyle(20);
  gXSecWithSyst -> SetMinimum(0.000001);
  gXSecWithSyst-> SetMarkerColor(kBlack);
  gXSecWithSyst -> SetLineColor (kBlack);
  gXSecWithSyst -> Draw("AP");

  gXSec -> SetMarkerSize(0.7);
  gXSec -> SetMarkerStyle(21);
  // gXSec -> Draw("p same");
  PrintIt(c1, cTitleExtd);


  TCanvas* c2 = new TCanvas("c2","",800,0,700,500);
  c2->SetFillColor(0);
  // defining the legends
  TLegend* tl1 = new TLegend(0.15,0.06,0.48,0.48);
  tl1 -> SetTextSize(0.05);
  tl1->SetFillColor(0);
  tl1->SetBorderSize(0);
  tl1->AddEntry(grMuon,        "2010 #mu channel, #int L dt = 36 pb^{-1} at #sqrt{7} TeV",   "lep");
  tl1->AddEntry(gXSecWithSyst, "2012 #mu channel, #int L dt = 18.4 pb^{-1} at #sqrt{8} TeV", "lep");


  DrawWithRatio(c2, cTitle, gXSecWithSyst, grMuon, tl1); 
 

  TCanvas* c3 = new TCanvas("c3","",0,800,700,500);
  c3->SetFillColor(0);

  TLegend* tl2 = new TLegend(0.15,0.06,0.48,0.48);
  tl2 -> SetTextSize(0.05);
  tl2->SetFillColor(0);
  tl2->SetBorderSize(0);
  tl2->AddEntry(grComb,       "2010 e+#mu comb, #int L dt = 36 pb^{-1} at #sqrt{7} TeV",   "lep");
  tl2->AddEntry(gXSecWithSyst,"2012 #mu channel, #int L dt = 18.4 pb^{-1} at #sqrt{8} TeV","lep");

  DrawWithRatio(c3, cTitle, gXSecWithSyst, grComb, tl2); 

     c1->SaveAs("diffXSec_muons_2012.png");
     c2->SaveAs("diffXSec_muons_2012_vs_muons_2010.png");
     c3->SaveAs("diffXSec_muons_2012_vs_muonsEles_2010.png");

 //==================================================================================
 //to plot ratio only
   TGraphAsymmErrors* gXSecRatio= new TGraphAsymmErrors(hXSecFSR->GetNbinsX());
  
  for(int bin=1; bin<=hXSecFSR->GetNbinsX(); bin++){

 Double_t x,yNum,yDen,error_highN, error_lowN, error_highD, error_lowD;
 gXSecWithSyst->GetPoint(bin,x,yNum);
 grMuon->GetPoint(bin,x,yDen);


 double sigmaN, sigmaD;
    error_highN = gXSecWithSyst->GetErrorYhigh(bin);
    error_lowN  = gXSecWithSyst->GetErrorYlow(bin);
    if (error_highN > error_lowN) sigmaN = error_highN;
    else sigmaN = error_lowN;

    error_highD = grMuon->GetErrorYhigh(bin);
    error_lowD  = grMuon->GetErrorYlow(bin);
    if (error_highD > error_lowD) sigmaD = error_highD;
    else sigmaD = error_lowD;
    double ratio = yNum/yDen;

double sigmaR = sqrt(((1/yNum)*(1/yNum))*(sigmaN*sigmaN)+ ((yNum/(yDen*yDen))*(yNum/(yDen*yDen)))*(sigmaD*sigmaD));
    std::cout << bin << ") " << yNum << "/" << yDen << ", Ratio= " << ratio << std::endl;

    // gXSecRatio->SetBinContent(bin+1,ratio);
    //gXSecRatio->SetBinError(bin+1,sigmaR);

  gXSecRatio->SetPoint(bin-1,hXSecFSR->GetXaxis()->GetBinCenter(bin),ratio);
  
  gXSecRatio->SetPointError(bin-1,error_xbins[bin-1],error_xbins[bin-1],sigmaR, sigmaR);
  }
  std::cout << std::endl << std::endl;

  TCanvas* cRatio = new TCanvas("cRatio","",0,0,700,500);
  cRatio->SetFillColor(0);
  cRatio->SetLogx();


  gXSecRatio -> SetTitle("CMS Preliminary");
  
  gXSecRatio -> GetXaxis() -> SetTitle("q_{T} [GeV/c]");
  gXSecRatio -> GetXaxis() -> SetTitleSize(0.04);
  gXSecRatio -> GetXaxis() -> SetLabelSize(0.03);

  gXSecRatio -> GetYaxis() -> SetTitle("(#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} at 8 TeV) / (#frac{1}{#sigma} #frac{d#sigma}{dq_{T}} at 7 TeV) ");
  gXSecRatio -> GetYaxis() -> SetTitleSize(0.04);
  gXSecRatio -> GetYaxis() -> SetLabelSize(0.03);
  gXSecRatio -> GetYaxis() -> SetTitleOffset(1.02);
  gXSecRatio ->GetXaxis()->SetRangeUser(0, 20);
  gXSecRatio ->GetYaxis()->SetRangeUser(0.5, 1.5);
  gXSecRatio -> SetMarkerSize(0.8);
  gXSecRatio -> SetMarkerStyle(20);
  gXSecRatio->Draw("AP");
  grPowhegRatio->SetLineWidth(4);
  grPowhegRatio->SetLineColor (kGreen);
  grPowhegRatio->Draw("CP same");

  grFEWZRatio->SetLineWidth(4);
  grFEWZRatio->SetLineStyle(8);
  grFEWZRatio->SetLineColor (kBlue);
  //  grFEWZRatio->Draw("CP same");

  TLine* line = new TLine(0,1,20,1);
  line -> SetLineColor( kRed );
  line -> SetLineWidth( 4 );
  line -> SetLineStyle( 2 );
  line -> Draw("same");

  TLegend* tl3 = new TLegend(0.15,0.7,0.48,0.85);
  tl3 -> SetTextSize(0.035);
  tl3->SetFillColor(0);
  tl3->SetBorderSize(0);
  tl3->AddEntry(gXSecRatio,       "data",   "l");
  tl3->AddEntry(grPowhegRatio,"Powheg","l");
  //  tl3->AddEntry(grFEWZRatio,"FEWZ","l");
  tl3->Draw("same");

  TPaveText *pave = new TPaveText(200,2.0,300,2.8);
  pave->SetFillColor(0);
  TText *t1=pave->AddText("#int L dt = 36 pb^{-1} at #sqrt{7} TeV");
  t1->SetTextSize(0.035);
  pave->AddText("#int L dt = 18.4 pb^{-1} at #sqrt{8} TeV");
  pave->SetTextSize(0.035);
  pave->Draw("same");
}


// ------------------------------------------------------------------------------
TGraphAsymmErrors* MuonChannel2010(){
  
  const int n = 18;
  
  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  
  double x[n];
  for (int i=0;i<n;i++) {
    x[i] = bin[i] + (bin[i+1]-bin[i])/2;
  }
  
  double y[18] = {0.0321, 0.0589, 0.0551, 0.039, 0.0349, 0.0274, 0.0223, 0.0168, 0.0114, 0.00632, 0.00353, 0.00174, 0.000776, 0.000487, 0.000179, 0.000071, 0.0000117, 0.00000224};
  double error_y[18] = {0.0014, 0.0021, 0.0020, 0.0018, 0.0016, 0.0015, 0.0014, 0.0012, 0.0004, 0.00028, 0.00021,0.0001, 0.000071, 0.000055, 0.000022, 0.000014, 0.0000051, 0.00000078};   

  double error_x[18];
  for (int i=0;i<n;i++) error_x[i]=0;

  TGraphAsymmErrors* grMuon = new TGraphAsymmErrors(n,x,y,error_x,error_x,error_y,error_y);
  
  
  return grMuon;
}

TGraphAsymmErrors* Combination2010(){
  
  const int n = 18;
  
  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  
  double x[n];
  for (int i=0;i<n;i++) {
    x[i] = bin[i] + (bin[i+1]-bin[i])/2;
  }
  
  double y[18] = {0.0322, 0.0592, 0.055, 0.0396, 0.0353, 0.0272, 0.0216, 0.0165, 0.0116, 0.00598, 0.00338, 0.00181, 0.000779, 0.000475, 0.000193, 0.00006, 0.0000151, 0.00000129};
  double error_y[18] = {0.0013, 0.0017, 0.0016, 0.0014, 0.0012, 0.0012, 0.0010, 0.0009, 0.0004, 0.00027, 0.00018, 0.00009, 0.000054, 0.000042, 0.000017, 0.0000099, 0.0000043, 0.00000044};   

  double error_x[18];
  for (int i=0;i<n;i++) error_x[i]=0;

  TGraphAsymmErrors* grComb = new TGraphAsymmErrors(n,x,y,error_x,error_x,error_y,error_y);

  return grComb;
}

//=======Powheg ratio 2012 and 2010==============================================


TGraphAsymmErrors* PowhegRatio(){
  
  const int n = 18;
  
  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  
  double x[n];
  for (int i=0;i<n;i++) {
    x[i] = bin[i] + (bin[i+1]-bin[i])/2;
  }

 double y[18] = {
                  0.92267014,  
		  0.97829653,  
		  0.97398096,  
		  1.00130156,  
		  1.01119743,  
		  1.03806587,  
		  1.05950312,  
		  1.04055990,  
		  1.03731212,  
		  1.00456587,  
		  0.98004959,  
		  0.98354966,  
		  0.99338198,  
		  0.98937692,  
		  0.92094711,  
		  1.01964791,  
		  0.96149786,  
		  0.88050312};

 double error_x[18];
  for (int i=0;i<n;i++) error_x[i]=0;

  TGraphAsymmErrors* grPowhegRatio = new TGraphAsymmErrors(n,x,y,0,0,0,0);

  return grPowhegRatio;
}

//===========================Ratio FEWZ==================================
TGraphAsymmErrors* FEWZRatio(){
  
  const int n = 10;
  
  double bin[n+1] = {20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
  
  double x[n];
  for (int i=0;i<n;i++) {
    x[i] = bin[i] + (bin[i+1]-bin[i])/2;
  }

  double y[10] = {1.02915868, 
		  1.03787334, 
		  1.04711568, 
		  1.05834165, 
		  1.06827474, 
		  1.08422607, 
		  1.09767529, 
		  1.12226142, 
		  1.15839652, 
		  1.24345988 
                           };

 double error_x[10];
  for (int i=0;i<n;i++) error_x[i]=0;

  TGraphAsymmErrors* grFEWZRatio = new TGraphAsymmErrors(n,x,y,0,0,0,0);

  return grFEWZRatio;
}

//=====================filling histograms===============================================
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

    hPt       -> Fill( recoCandPt );

  }
  std::cout << " DONE!" << std::endl;
  return;
}

//------------------------------------------------------------------------------
// Draw projections and residuals
//------------------------------------------------------------------------------
void DrawWithRatio(TCanvas *canvas, char *cTitle,
                   TGraphAsymmErrors* gNum, TGraphAsymmErrors* gDen,
                   TLegend* leg) {
  
  // sanity check
  if (gNum->GetN() != gDen->GetN()){
    std::cout<< " *** Error: number of points not consistent between the"
             << " two graphs -> Exit!\n";
    return;
  }



  TH1F *hNum = ConvertGraphToHisto(gNum);
  TH1F *hDen = ConvertGraphToHisto(gDen);

  //hNum -> Sumw2();
  //hDen -> Sumw2();

  TH1F *hPull = (TH1F*) hNum -> Clone("hPull");
  hPull->Sumw2();
  hPull->Reset();

  // Print the ratio
  std::cout << std::endl << std::endl;
  for (int iBin=0;iBin<gNum->GetN();iBin++) {
    Double_t x,yNum,yDen, error_highN, error_lowN, error_highD, error_lowD;
    gNum->GetPoint(iBin,x,yNum);
    gDen->GetPoint(iBin,x,yDen);

    double sigmaN, sigmaD;
    error_highN = gNum->GetErrorYhigh(iBin);
    error_lowN  = gNum->GetErrorYlow(iBin);
    if (error_highN > error_lowN) sigmaN = error_highN;
    else sigmaN = error_lowN;

    error_highD = gDen->GetErrorYhigh(iBin);
    error_lowD  = gDen->GetErrorYlow(iBin);
    if (error_highD > error_lowD) sigmaD = error_highD;
    else sigmaD = error_lowD;

    double ratio = yNum/yDen;
 
    double sigmaR = sqrt(((1/yNum)*(1/yNum))*(sigmaN*sigmaN)+ ((yNum/(yDen*yDen))*(yNum/(yDen*yDen)))*(sigmaD*sigmaD));
       std::cout << iBin << ") " << yNum << "/" << yDen << ", Ratio= " << ratio << std::endl;
    hPull->SetBinContent(iBin+1,ratio);
    hPull->SetBinError(iBin+1,sigmaR);

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
  pad2->SetBottomMargin(0.33);
  pad2->SetTopMargin   (0.10);

  pad1->Draw(); // Projections pad
  pad2->Draw(); // Residuals   pad

  pad1->cd();

  hNum -> SetMaximum(10E-2);
  hNum -> SetMinimum(10E-7);

  char yAxisName[200];
  sprintf(yAxisName,"%s",gDen->GetYaxis()->GetTitle());
  hNum -> GetYaxis() -> SetTitle( yAxisName );
  hNum -> GetYaxis() -> SetTitleSize  (0.06);
  hNum -> GetYaxis() -> SetTitleOffset(0.85);
  
  hNum->Draw("axis");
  gNum->Draw("P");
  gDen->Draw("P");

  PrintIt(pad1,cTitle);

  if (leg) leg ->Draw("same");

  //----------------------------------------------------------------------------
  // Residuals pad
  //----------------------------------------------------------------------------
  pad2->cd();

  TAxis *xPull = NULL;
  TAxis *yPull = NULL;
  char xAxisName[200];
  sprintf(xAxisName,"%s",gDen->GetXaxis()->GetTitle());
  axis1F(hPull,xPull,yPull,xAxisName,"Ratio");

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
  
  hPull->GetYaxis()->SetLabelSize  (0.09);
  hPull->GetYaxis()->CenterTitle(1);
  hPull->GetYaxis()->SetTitleOffset(0.45);
  hPull->GetYaxis()->SetTitleSize  (0.12);

  hPull->SetMaximum(4);
  hPull->SetMinimum(0);
  hPull->SetMarkerSize(1.0);
  hPull->Draw("PEL");

  TLine* line = new TLine(0,1,600,1);
  line -> SetLineColor( kRed );
  line -> SetLineWidth( 4 );
  line -> SetLineStyle( 2 );
  line -> Draw("same");
  
  pad2->Update();
  //pad2->GetFrame()->DrawClone();

}

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

  xaxis->SetLabelFont(42);
  yaxis->SetLabelFont(42);
  xaxis->SetLabelOffset(0.005);
  yaxis->SetLabelOffset(0.005);
  //xaxis->SetLabelSize(0.04);
  yaxis->SetLabelSize(0.04);

  xaxis->SetNdivisions(505);
  yaxis->SetNdivisions(505);

  xaxis->SetTitle(xtitle);
  yaxis->SetTitle(ytitle);
  xaxis->SetTitleColor(kBlack);
  yaxis->SetTitleColor(kBlack);
  xaxis->SetTitleFont(42);
  yaxis->SetTitleFont(42);
  xaxis->SetTitleOffset(1.0);
  yaxis->SetTitleOffset(1.3);
  //xaxis->SetTitleSize(0.045);
  //yaxis->SetTitleSize(0.045);
  yaxis->CenterTitle(kTRUE);
}



TH1F* ConvertGraphToHisto(TGraph *pGraph){
  // takes data from a graph, determines binning and fills data into histogram
  Int_t NPoints = pGraph->GetN();
  //std::cout << "NPoints=" << NPoints << std::endl;
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
    Double_t x,y;
    pGraph->GetPoint(i,x,y);
    ThisHist->SetBinContent(i+1,y);
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


