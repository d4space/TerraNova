#include <iostream>
#include <set>

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
#include "TEfficiency.h"

#include "general_ratio.h"
#include "selections.h"
#include <time.h>
#include "zrapidityStandard.C"

 const int   nptBins  =   9;

void fillHistos(TTree* tree, 
                TH1F *hMass, 
                TH1F *hPt,
                TH1F *hPt20B, 
                TH1F *hRapidity,
                TH1F *hVertex,
                double mLow,
                double mHigh,
                double ptLow,
                double ptHigh,
                int mode);

// ---- CROSS SECTIONS ----

// in pb
double xsec_dymumu     = 5745.25/3.; //(*)
double xsec_ttbarjets  = 225.197;
double xsec_ztautau    = 5745.25/3.;
double xsec_WW         = 57.1097;
double xsec_WZ         = 32.3161;
double xsec_ZZ         = 8.25561;
double xsec_WJetsToLNu = 30400;  //(*) https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
double xsec_QCD        = 134680;
// - END OF CROSS SECTIONS -


// ------------------------------------------------------------------------------
// number of original MC events
// this you will get it from crab
int nMCEvDYToMuMu    = 1980204; // DYJetsToLL 
int nMCEvTTbarJets   = 6736135;
int nMCEvZtautau     = 1737776;
int nMCEvWW          = 2519385;
int nMCEvWZ          = 2182479;
int nMCEvZZ          = 1890152;
int nMCEvWJetsToLNu  = 2053500;
int nMCEvQCD         = 7529312;

//double intLumi_invpb = 19.378; //add your number
double intLumi_invpb = 18.424; //withecal reprocessing(190949) + adding 193112
//double intLumi_invpb = 18.729; // cross check with inclusive WZ

void dataMC_ratio(int    mBins  =   30, // Z mass
                   double mLow   =   60, // Z mass
                   double mHigh  =  120, // Z mass
                   int   ptBins  =   120, //  195, 
                   double ptLow  =   0, // 50, //   50,
                   double ptHigh =  600, // 950, // 2000,
                   int    rBins  =   48, // Z rapidity
                   double rLow   = -2.4, // Z rapidity
                   double rHigh  =  2.4, // Z rapidity 
                   bool isSave   = true) {

  //double xbins_pt[nptBins+1] = {0.001, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};

     //double xbins_pt[nptBins+1] = {0.001, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30}; // From bin1
     double xbins_pt[nptBins+1] = {30, 40, 50, 70, 90, 110, 150, 190, 250, 600}; // After bin9
  
  
 TFile *Hist_out = new TFile("histo_withecalreprocessing.root","RECREATE");

  gROOT->Clear();
  gStyle->SetOptStat(0);  
  //gStyle->SetOptStat(111111);  

  loadGeneral();

  // ---------------------------------------------------------------------------
  // local variables
  TString png = "png/";
  TString root= "root/";
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Fill the data
    
  //  TFile* fileData = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_5_ecalpatch1/Data/ntuplesWITHdz/SingleMuRun2012A_190949_191090_191367_193112_193116.root"); 
  //TFile* fileData = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_5_ecalpatch1/Data/ntuplesWITHdz/fullselection/data2012A_LowPU.root");   
  TFile* fileData = new TFile("./fullselection/data2012A_LowPU.root");   

     TTree* treeData = (TTree*) fileData -> Get("tree");

   
  TH1F* hMassData = new TH1F("hMassData", "",  mBins,  mLow,  mHigh);
  TH1F* hPtData   = new TH1F("hPtData"  , "", nptBins, xbins_pt);
  TH1F* hPt20BData   = new TH1F("hPt20BData"  , "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityData = new TH1F("hRapidityData", "", rBins,  rLow,  rHigh);
  TH1F* hVertexData = new TH1F("hVertexData", "", 25,  0,  25);
 

  std::cout << "Filling Collisions Data Events\n";
    fillHistos(treeData, hMassData, hPtData, hPt20BData, hRapidityData, hVertexData, mLow, mHigh, ptLow, ptHigh, 0);

   hPtData->SetBinContent(1,hPtData->GetBinContent(1)/2.5);
   hPtData->SetBinContent(2,hPtData->GetBinContent(2)/2.5);
   hPtData->SetBinContent(3,hPtData->GetBinContent(3)/2.5);
   hPtData->SetBinContent(4,hPtData->GetBinContent(4)/2.5);
   hPtData->SetBinContent(5,hPtData->GetBinContent(5)/2.5);
   hPtData->SetBinContent(6,hPtData->GetBinContent(6)/2.5);
   hPtData->SetBinContent(7,hPtData->GetBinContent(7)/2.5);
   hPtData->SetBinContent(8,hPtData->GetBinContent(8)/2.5);
   hPtData->SetBinContent(9,hPtData->GetBinContent(9)/10.0);
   hPtData->SetBinContent(10,hPtData->GetBinContent(10)/10.0);
   hPtData->SetBinContent(11,hPtData->GetBinContent(11)/10.0);
   hPtData->SetBinContent(12,hPtData->GetBinContent(12)/20.0);
   hPtData->SetBinContent(13,hPtData->GetBinContent(13)/20.0);
   hPtData->SetBinContent(14,hPtData->GetBinContent(14)/20.0);
   hPtData->SetBinContent(15,hPtData->GetBinContent(15)/40.0);
   hPtData->SetBinContent(16,hPtData->GetBinContent(16)/40.0);
   hPtData->SetBinContent(17,hPtData->GetBinContent(17)/60.0);
   hPtData->SetBinContent(18,hPtData->GetBinContent(18)/350.0);
   
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Fill the DYToMuMu MC
   //  TFile* fileMC_dymumu = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/DYToMuMu/mergeFile/DYToMuMu_merged.root");
  //TFile* fileMC_dymumu = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/DYToMuMu.root");
  TFile* fileMC_dymumu = new TFile("./fullselection/DYToMuMu.root");
  //TFile* fileMC_dymumu = new TFile("./DYToMuMu_merged.root");
  TTree* treeMC_dymumu = (TTree*) fileMC_dymumu -> Get("tree");
  
  TH1F* hMassMC_dymumu = new TH1F("hMassMC_dymumu", "", mBins, mLow, mHigh);
  TH1F* hPtMC_dymumu   = new TH1F("hPtMC_dymumu",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_dymumu   = new TH1F("hPt20BMC_dymumu",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_dymumu    = new TH1F("hRapidityMC_dymumu", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_dymumu      = new TH1F("hVertexMC_dymumu", "", 25, 0, 25);
  
  std::cout << "filling DYToMuMu events\n";
     fillHistos(treeMC_dymumu, hMassMC_dymumu, hPtMC_dymumu, hPt20BMC_dymumu, hRapidityMC_dymumu, hVertexMC_dymumu, mLow, mHigh, ptLow, ptHigh, 0);


   hPtMC_dymumu->SetBinContent(1,hPtMC_dymumu->GetBinContent(1)/2.5);
   hPtMC_dymumu->SetBinContent(2,hPtMC_dymumu->GetBinContent(2)/2.5);
   hPtMC_dymumu->SetBinContent(3,hPtMC_dymumu->GetBinContent(3)/2.5);
   hPtMC_dymumu->SetBinContent(4,hPtMC_dymumu->GetBinContent(4)/2.5);
   hPtMC_dymumu->SetBinContent(5,hPtMC_dymumu->GetBinContent(5)/2.5);
   hPtMC_dymumu->SetBinContent(6,hPtMC_dymumu->GetBinContent(6)/2.5);
   hPtMC_dymumu->SetBinContent(7,hPtMC_dymumu->GetBinContent(7)/2.5);
   hPtMC_dymumu->SetBinContent(8,hPtMC_dymumu->GetBinContent(8)/2.5);
   hPtMC_dymumu->SetBinContent(9,hPtMC_dymumu->GetBinContent(9)/10.0);
   hPtMC_dymumu->SetBinContent(10,hPtMC_dymumu->GetBinContent(10)/10.0);
   hPtMC_dymumu->SetBinContent(11,hPtMC_dymumu->GetBinContent(11)/10.0);
   hPtMC_dymumu->SetBinContent(12,hPtMC_dymumu->GetBinContent(12)/20.0);
   hPtMC_dymumu->SetBinContent(13,hPtMC_dymumu->GetBinContent(13)/20.0);
   hPtMC_dymumu->SetBinContent(14,hPtMC_dymumu->GetBinContent(14)/20.0);
   hPtMC_dymumu->SetBinContent(15,hPtMC_dymumu->GetBinContent(15)/40.0);
   hPtMC_dymumu->SetBinContent(16,hPtMC_dymumu->GetBinContent(16)/40.0);
   hPtMC_dymumu->SetBinContent(17,hPtMC_dymumu->GetBinContent(17)/60.0);
   hPtMC_dymumu->SetBinContent(18,hPtMC_dymumu->GetBinContent(18)/350.0);

  // Fill the TTbarJets MC  
   //  TFile* fileMC_ttbarjets = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/TTJets/mergeFile/TTJets_merged.root");
  //TFile* fileMC_ttbarjets = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/TTJets.root");
  TFile* fileMC_ttbarjets = new TFile("./fullselection/TTJets.root");
  TTree* treeMC_ttbarjets = (TTree*) fileMC_ttbarjets -> Get("tree");

  TH1F* hMassMC_ttbarjets = new TH1F("hMassMC_ttbarjets", "", mBins, mLow, mHigh);
  TH1F* hPtMC_ttbarjets   = new TH1F("hPtMC_ttbarjets",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_ttbarjets   = new TH1F("hPt20BMC_ttbarjets",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_ttbarjets    = new TH1F("hRapidityMC_ttbarjets", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_ttbarjets      = new TH1F("hVertexMC_ttbarjets", "", 25, 0, 25);

  std::cout << "filling TTbar+Jets events\n";
      fillHistos(treeMC_ttbarjets, hMassMC_ttbarjets, hPtMC_ttbarjets, hPt20BMC_ttbarjets, hRapidityMC_ttbarjets, hVertexMC_ttbarjets, mLow, mHigh, ptLow, ptHigh, 0);

   hPtMC_ttbarjets->SetBinContent(1,hPtMC_ttbarjets->GetBinContent(1)/2.5);
   hPtMC_ttbarjets->SetBinContent(2,hPtMC_ttbarjets->GetBinContent(2)/2.5);
   hPtMC_ttbarjets->SetBinContent(3,hPtMC_ttbarjets->GetBinContent(3)/2.5);
   hPtMC_ttbarjets->SetBinContent(4,hPtMC_ttbarjets->GetBinContent(4)/2.5);
   hPtMC_ttbarjets->SetBinContent(5,hPtMC_ttbarjets->GetBinContent(5)/2.5);
   hPtMC_ttbarjets->SetBinContent(6,hPtMC_ttbarjets->GetBinContent(6)/2.5);
   hPtMC_ttbarjets->SetBinContent(7,hPtMC_ttbarjets->GetBinContent(7)/2.5);
   hPtMC_ttbarjets->SetBinContent(8,hPtMC_ttbarjets->GetBinContent(8)/2.5);
   hPtMC_ttbarjets->SetBinContent(9,hPtMC_ttbarjets->GetBinContent(9)/10.0);
   hPtMC_ttbarjets->SetBinContent(10,hPtMC_ttbarjets->GetBinContent(10)/10.0);
   hPtMC_ttbarjets->SetBinContent(11,hPtMC_ttbarjets->GetBinContent(11)/10.0);
   hPtMC_ttbarjets->SetBinContent(12,hPtMC_ttbarjets->GetBinContent(12)/20.0);
   hPtMC_ttbarjets->SetBinContent(13,hPtMC_ttbarjets->GetBinContent(13)/20.0);
   hPtMC_ttbarjets->SetBinContent(14,hPtMC_ttbarjets->GetBinContent(14)/20.0);
   hPtMC_ttbarjets->SetBinContent(15,hPtMC_ttbarjets->GetBinContent(15)/40.0);
   hPtMC_ttbarjets->SetBinContent(16,hPtMC_ttbarjets->GetBinContent(16)/40.0);
   hPtMC_ttbarjets->SetBinContent(17,hPtMC_ttbarjets->GetBinContent(17)/60.0);
   hPtMC_ttbarjets->SetBinContent(18,hPtMC_ttbarjets->GetBinContent(18)/350.0);

  //// Fill the ZTauTau MC  
   // TFile* fileMC_ztautau = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/DYToTauTau/mergeFile/DYToTauTau_merged.root");
   //TFile* fileMC_ztautau = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/DYToTauTau.root");
   TFile* fileMC_ztautau = new TFile("./fullselection/DYToTauTau.root");
  TTree* treeMC_ztautau = (TTree*) fileMC_ztautau -> Get("tree");

  TH1F* hMassMC_ztautau = new TH1F("hMassMC_ztautau", "", mBins, mLow, mHigh);
  TH1F* hPtMC_ztautau   = new TH1F("hPtMC_ztautau",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_ztautau   = new TH1F("hPt20BMC_ztautau",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_ztautau    = new TH1F("hRapidityMC_ztautau", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_ztautau      = new TH1F("hVertexMC_ztautau", "", 25, 0, 25);

 std::cout << "filling ZToTauTau events\n";
  fillHistos(treeMC_ztautau, hMassMC_ztautau, hPtMC_ztautau, hPt20BMC_ztautau, hRapidityMC_ztautau, hVertexMC_ztautau, mLow, mHigh, ptLow, ptHigh, 0);

   hPtMC_ztautau->SetBinContent(1,hPtMC_ztautau->GetBinContent(1)/2.5);
   hPtMC_ztautau->SetBinContent(2,hPtMC_ztautau->GetBinContent(2)/2.5);
   hPtMC_ztautau->SetBinContent(3,hPtMC_ztautau->GetBinContent(3)/2.5);
   hPtMC_ztautau->SetBinContent(4,hPtMC_ztautau->GetBinContent(4)/2.5);
   hPtMC_ztautau->SetBinContent(5,hPtMC_ztautau->GetBinContent(5)/2.5);
   hPtMC_ztautau->SetBinContent(6,hPtMC_ztautau->GetBinContent(6)/2.5);
   hPtMC_ztautau->SetBinContent(7,hPtMC_ztautau->GetBinContent(7)/2.5);
   hPtMC_ztautau->SetBinContent(8,hPtMC_ztautau->GetBinContent(8)/2.5);
   hPtMC_ztautau->SetBinContent(9,hPtMC_ztautau->GetBinContent(9)/10.0);
   hPtMC_ztautau->SetBinContent(10,hPtMC_ztautau->GetBinContent(10)/10.0);
   hPtMC_ztautau->SetBinContent(11,hPtMC_ztautau->GetBinContent(11)/10.0);
   hPtMC_ztautau->SetBinContent(12,hPtMC_ztautau->GetBinContent(12)/20.0);
   hPtMC_ztautau->SetBinContent(13,hPtMC_ztautau->GetBinContent(13)/20.0);
   hPtMC_ztautau->SetBinContent(14,hPtMC_ztautau->GetBinContent(14)/20.0);
   hPtMC_ztautau->SetBinContent(15,hPtMC_ztautau->GetBinContent(15)/40.0);
   hPtMC_ztautau->SetBinContent(16,hPtMC_ztautau->GetBinContent(16)/40.0);
   hPtMC_ztautau->SetBinContent(17,hPtMC_ztautau->GetBinContent(17)/60.0);
   hPtMC_ztautau->SetBinContent(18,hPtMC_ztautau->GetBinContent(18)/350.0);

 

 // Fill the WZ MC  
   //  TFile* fileMC_WZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WZ/mergeFile/WZ_merged.root");
  //TFile* fileMC_WZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/WZ.root");
  TFile* fileMC_WZ = new TFile("./fullselection/WZ.root");
  TTree* treeMC_WZ = (TTree*) fileMC_WZ -> Get("tree");

  TH1F* hMassMC_WZ = new TH1F("hMassMC_WZ", "", mBins, mLow, mHigh);
  TH1F* hPtMC_WZ   = new TH1F("hPtMC_WZ",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_WZ   = new TH1F("hPt20BMC_WZ",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_WZ    = new TH1F("hRapidityMC_WZ", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_WZ      = new TH1F("hVertexMC_WZ", "", 25, 0, 25);
 
  std::cout << "filling WZ events\n";
    fillHistos(treeMC_WZ, hMassMC_WZ, hPtMC_WZ, hPt20BMC_WZ, hRapidityMC_WZ, hVertexMC_WZ, mLow, mHigh, ptLow, ptHigh, 0);

   hPtMC_WZ->SetBinContent(1,hPtMC_WZ->GetBinContent(1)/2.5);
   hPtMC_WZ->SetBinContent(2,hPtMC_WZ->GetBinContent(2)/2.5);
   hPtMC_WZ->SetBinContent(3,hPtMC_WZ->GetBinContent(3)/2.5);
   hPtMC_WZ->SetBinContent(4,hPtMC_WZ->GetBinContent(4)/2.5);
   hPtMC_WZ->SetBinContent(5,hPtMC_WZ->GetBinContent(5)/2.5);
   hPtMC_WZ->SetBinContent(6,hPtMC_WZ->GetBinContent(6)/2.5);
   hPtMC_WZ->SetBinContent(7,hPtMC_WZ->GetBinContent(7)/2.5);
   hPtMC_WZ->SetBinContent(8,hPtMC_WZ->GetBinContent(8)/2.5);
   hPtMC_WZ->SetBinContent(9,hPtMC_WZ->GetBinContent(9)/10.0);
   hPtMC_WZ->SetBinContent(10,hPtMC_WZ->GetBinContent(10)/10.0);
   hPtMC_WZ->SetBinContent(11,hPtMC_WZ->GetBinContent(11)/10.0);
   hPtMC_WZ->SetBinContent(12,hPtMC_WZ->GetBinContent(12)/20.0);
   hPtMC_WZ->SetBinContent(13,hPtMC_WZ->GetBinContent(13)/20.0);
   hPtMC_WZ->SetBinContent(14,hPtMC_WZ->GetBinContent(14)/20.0);
   hPtMC_WZ->SetBinContent(15,hPtMC_WZ->GetBinContent(15)/40.0);
   hPtMC_WZ->SetBinContent(16,hPtMC_WZ->GetBinContent(16)/40.0);
   hPtMC_WZ->SetBinContent(17,hPtMC_WZ->GetBinContent(17)/60.0);
   hPtMC_WZ->SetBinContent(18,hPtMC_WZ->GetBinContent(18)/350.0);

// Fill the ZZ MC  
//  TFile* fileMC_ZZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/ZZ/mergeFile/ZZ_merged.root");
  //TFile* fileMC_ZZ = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/ZZ.root");
  TFile* fileMC_ZZ = new TFile("./fullselection/ZZ.root");
  TTree* treeMC_ZZ = (TTree*) fileMC_ZZ -> Get("tree");

  TH1F* hMassMC_ZZ = new TH1F("hMassMC_ZZ", "", mBins, mLow, mHigh);
  TH1F* hPtMC_ZZ   = new TH1F("hPtMC_ZZ",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_ZZ   = new TH1F("hPt20BMC_ZZ",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_ZZ    = new TH1F("hRapidityMC_ZZ", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_ZZ      = new TH1F("hVertexMC_ZZ", "", 25, 0, 25);

  std::cout << "filling ZZ events\n";
    fillHistos(treeMC_ZZ, hMassMC_ZZ, hPtMC_ZZ, hPt20BMC_ZZ, hRapidityMC_ZZ, hVertexMC_ZZ, mLow, mHigh, ptLow, ptHigh, 0);

   hPtMC_ZZ->SetBinContent(1,hPtMC_ZZ->GetBinContent(1)/2.5);
   hPtMC_ZZ->SetBinContent(2,hPtMC_ZZ->GetBinContent(2)/2.5);
   hPtMC_ZZ->SetBinContent(3,hPtMC_ZZ->GetBinContent(3)/2.5);
   hPtMC_ZZ->SetBinContent(4,hPtMC_ZZ->GetBinContent(4)/2.5);
   hPtMC_ZZ->SetBinContent(5,hPtMC_ZZ->GetBinContent(5)/2.5);
   hPtMC_ZZ->SetBinContent(6,hPtMC_ZZ->GetBinContent(6)/2.5);
   hPtMC_ZZ->SetBinContent(7,hPtMC_ZZ->GetBinContent(7)/2.5);
   hPtMC_ZZ->SetBinContent(8,hPtMC_ZZ->GetBinContent(8)/2.5);
   hPtMC_ZZ->SetBinContent(9,hPtMC_ZZ->GetBinContent(9)/10.0);
   hPtMC_ZZ->SetBinContent(10,hPtMC_ZZ->GetBinContent(10)/10.0);
   hPtMC_ZZ->SetBinContent(11,hPtMC_ZZ->GetBinContent(11)/10.0);
   hPtMC_ZZ->SetBinContent(12,hPtMC_ZZ->GetBinContent(12)/20.0);
   hPtMC_ZZ->SetBinContent(13,hPtMC_ZZ->GetBinContent(13)/20.0);
   hPtMC_ZZ->SetBinContent(14,hPtMC_ZZ->GetBinContent(14)/20.0);
   hPtMC_ZZ->SetBinContent(15,hPtMC_ZZ->GetBinContent(15)/40.0);
   hPtMC_ZZ->SetBinContent(16,hPtMC_ZZ->GetBinContent(16)/40.0);
   hPtMC_ZZ->SetBinContent(17,hPtMC_ZZ->GetBinContent(17)/60.0);
   hPtMC_ZZ->SetBinContent(18,hPtMC_ZZ->GetBinContent(18)/350.0);

// Fill the WW MC  
//  TFile* fileMC_WW = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WW/mergeFile/WW_merged.root");
  //TFile* fileMC_WW = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/WW.root");
  TFile* fileMC_WW = new TFile("./fullselection/WW.root");
  TTree* treeMC_WW = (TTree*) fileMC_WW -> Get("tree");

  TH1F* hMassMC_WW = new TH1F("hMassMC_WW", "", mBins, mLow, mHigh);
  TH1F* hPtMC_WW   = new TH1F("hPtMC_WW",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_WW   = new TH1F("hPt20BMC_WW",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_WW    = new TH1F("hRapidityMC_WW", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_WW      = new TH1F("hVertexMC_WW", "", 25, 0, 25);

  std::cout << "filling WW events\n";
    fillHistos(treeMC_WW, hMassMC_WW, hPtMC_WW, hPt20BMC_WW, hRapidityMC_WW, hVertexMC_WW, mLow, mHigh, ptLow, ptHigh, 0);
 
   hPtMC_WW->SetBinContent(1,hPtMC_WW->GetBinContent(1)/2.5);
   hPtMC_WW->SetBinContent(2,hPtMC_WW->GetBinContent(2)/2.5);
   hPtMC_WW->SetBinContent(3,hPtMC_WW->GetBinContent(3)/2.5);
   hPtMC_WW->SetBinContent(4,hPtMC_WW->GetBinContent(4)/2.5);
   hPtMC_WW->SetBinContent(5,hPtMC_WW->GetBinContent(5)/2.5);
   hPtMC_WW->SetBinContent(6,hPtMC_WW->GetBinContent(6)/2.5);
   hPtMC_WW->SetBinContent(7,hPtMC_WW->GetBinContent(7)/2.5);
   hPtMC_WW->SetBinContent(8,hPtMC_WW->GetBinContent(8)/2.5);
   hPtMC_WW->SetBinContent(9,hPtMC_WW->GetBinContent(9)/10.0);
   hPtMC_WW->SetBinContent(10,hPtMC_WW->GetBinContent(10)/10.0);
   hPtMC_WW->SetBinContent(11,hPtMC_WW->GetBinContent(11)/10.0);
   hPtMC_WW->SetBinContent(12,hPtMC_WW->GetBinContent(12)/20.0);
   hPtMC_WW->SetBinContent(13,hPtMC_WW->GetBinContent(13)/20.0);
   hPtMC_WW->SetBinContent(14,hPtMC_WW->GetBinContent(14)/20.0);
   hPtMC_WW->SetBinContent(15,hPtMC_WW->GetBinContent(15)/40.0);
   hPtMC_WW->SetBinContent(16,hPtMC_WW->GetBinContent(16)/40.0);
   hPtMC_WW->SetBinContent(17,hPtMC_WW->GetBinContent(17)/60.0);
   hPtMC_WW->SetBinContent(18,hPtMC_WW->GetBinContent(18)/350.0);

// Fill the WJetsToLNu MC  
   // TFile* fileMC_WJetsToLNu = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/WJetsToLNu/mergeFile/WJetsToLNu_merged.root");
  //TFile* fileMC_WJetsToLNu = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/WJetsToLNu.root");
  TFile* fileMC_WJetsToLNu = new TFile("./fullselection/WJetsToLNu.root");
  TTree* treeMC_WJetsToLNu = (TTree*) fileMC_WJetsToLNu -> Get("tree");

  TH1F* hMassMC_WJetsToLNu = new TH1F("hMassMC_WJetsToLNu", "", mBins, mLow, mHigh);
  TH1F* hPtMC_WJetsToLNu   = new TH1F("hPtMC_WJetsToLNu",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_WJetsToLNu   = new TH1F("hPt20BMC_WJetsToLNu",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_WJetsToLNu   = new TH1F("hRapidityMC_WJetsToLNu", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_WJetsToLNu     = new TH1F("hVertexMC_WJetsToLNu", "", 25, 0, 25);

  std::cout << "filling WJetsToLNu events\n";
    fillHistos(treeMC_WJetsToLNu, hMassMC_WJetsToLNu, hPtMC_WJetsToLNu, hPt20BMC_WJetsToLNu, hRapidityMC_WJetsToLNu, hVertexMC_WJetsToLNu, mLow, mHigh, ptLow, ptHigh, 0);
 
   hPtMC_WJetsToLNu->SetBinContent(1,hPtMC_WJetsToLNu->GetBinContent(1)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(2,hPtMC_WJetsToLNu->GetBinContent(2)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(3,hPtMC_WJetsToLNu->GetBinContent(3)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(4,hPtMC_WJetsToLNu->GetBinContent(4)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(5,hPtMC_WJetsToLNu->GetBinContent(5)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(6,hPtMC_WJetsToLNu->GetBinContent(6)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(7,hPtMC_WJetsToLNu->GetBinContent(7)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(8,hPtMC_WJetsToLNu->GetBinContent(8)/2.5);
   hPtMC_WJetsToLNu->SetBinContent(9,hPtMC_WJetsToLNu->GetBinContent(9)/10.0);
   hPtMC_WJetsToLNu->SetBinContent(10,hPtMC_WJetsToLNu->GetBinContent(10)/10.0);
   hPtMC_WJetsToLNu->SetBinContent(11,hPtMC_WJetsToLNu->GetBinContent(11)/10.0);
   hPtMC_WJetsToLNu->SetBinContent(12,hPtMC_WJetsToLNu->GetBinContent(12)/20.0);
   hPtMC_WJetsToLNu->SetBinContent(13,hPtMC_WJetsToLNu->GetBinContent(13)/20.0);
   hPtMC_WJetsToLNu->SetBinContent(14,hPtMC_WJetsToLNu->GetBinContent(14)/20.0);
   hPtMC_WJetsToLNu->SetBinContent(15,hPtMC_WJetsToLNu->GetBinContent(15)/40.0);
   hPtMC_WJetsToLNu->SetBinContent(16,hPtMC_WJetsToLNu->GetBinContent(16)/40.0);
   hPtMC_WJetsToLNu->SetBinContent(17,hPtMC_WJetsToLNu->GetBinContent(17)/60.0);
   hPtMC_WJetsToLNu->SetBinContent(18,hPtMC_WJetsToLNu->GetBinContent(18)/350.0);


//Fill QCD MC
//  TFile* fileMC_QCD = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/QCD/mergeFile/QCD_Pt_20_MuEnrichedPt_15_merged.root");
  //TFile* fileMC_QCD = new TFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/fullselection/QCD.root");
  TFile* fileMC_QCD = new TFile("./fullselection/QCD.root");
  TTree* treeMC_QCD = (TTree*) fileMC_QCD -> Get("tree");

  TH1F* hMassMC_QCD = new TH1F("hMassMC_QCD", "", mBins, mLow, mHigh);
  TH1F* hPtMC_QCD   = new TH1F("hPtMC_QCD",   "", nptBins, xbins_pt);
  TH1F* hPt20BMC_QCD   = new TH1F("hPt20BMC_QCD",   "", ptBins, ptLow, ptHigh);
  TH1F* hRapidityMC_QCD    = new TH1F("hRapidityMC_QCD", "", rBins, rLow, rHigh);
  TH1F* hVertexMC_QCD      = new TH1F("hVertexMC_QCD", "", 25, 0, 25);

  std::cout << "filling QCD events\n";
     fillHistos(treeMC_QCD, hMassMC_QCD, hPtMC_QCD, hPt20BMC_QCD, hRapidityMC_QCD, hVertexMC_QCD, mLow, mHigh, ptLow, ptHigh, 0);

   hPtMC_QCD->SetBinContent(1,hPtMC_QCD->GetBinContent(1)/2.5);
   hPtMC_QCD->SetBinContent(2,hPtMC_QCD->GetBinContent(2)/2.5);
   hPtMC_QCD->SetBinContent(3,hPtMC_QCD->GetBinContent(3)/2.5);
   hPtMC_QCD->SetBinContent(4,hPtMC_QCD->GetBinContent(4)/2.5);
   hPtMC_QCD->SetBinContent(5,hPtMC_QCD->GetBinContent(5)/2.5);
   hPtMC_QCD->SetBinContent(6,hPtMC_QCD->GetBinContent(6)/2.5);
   hPtMC_QCD->SetBinContent(7,hPtMC_QCD->GetBinContent(7)/2.5);
   hPtMC_QCD->SetBinContent(8,hPtMC_QCD->GetBinContent(8)/2.5);
   hPtMC_QCD->SetBinContent(9,hPtMC_QCD->GetBinContent(9)/10.0);
   hPtMC_QCD->SetBinContent(10,hPtMC_QCD->GetBinContent(10)/10.0);
   hPtMC_QCD->SetBinContent(11,hPtMC_QCD->GetBinContent(11)/10.0);
   hPtMC_QCD->SetBinContent(12,hPtMC_QCD->GetBinContent(12)/20.0);
   hPtMC_QCD->SetBinContent(13,hPtMC_QCD->GetBinContent(13)/20.0);
   hPtMC_QCD->SetBinContent(14,hPtMC_QCD->GetBinContent(14)/20.0);
   hPtMC_QCD->SetBinContent(15,hPtMC_QCD->GetBinContent(15)/40.0); 
   hPtMC_QCD->SetBinContent(16,hPtMC_QCD->GetBinContent(16)/40.0);
   hPtMC_QCD->SetBinContent(17,hPtMC_QCD->GetBinContent(17)/60.0);
   hPtMC_QCD->SetBinContent(18,hPtMC_QCD->GetBinContent(18)/350.0);

  TH1F* hMassMC_allBkg    = new TH1F("hMassMC_allBkg",   "", mBins, mLow, mHigh);
  THStack* sMassMC_allBkg = new THStack("sMassMC_allBkg","");

  TH1F* hPt20BMC_allBkg    = new TH1F("hPt20BMC_allBkg",   "", ptBins, ptLow, ptHigh);
  THStack* sPt20BMC_allBkg = new THStack("sPt20BMC_allBkg","");

  TH1F* hPtMC_allBkg    = new TH1F("hPtMC_allBkg",   "", nptBins, xbins_pt);
  TH1F* hPtMC_EWKBkg    = new TH1F("hPtMC_EWKBkg",   "", nptBins, xbins_pt);
  THStack* sPtMC_allBkg = new THStack("sPtMC_allBkg","");
  THStack* sPtMC_EWKBkg = new THStack("sPtMC_EWKBkg","");

  TH1F* hRapidityMC_allBkg   = new TH1F("hRapidityMC_allBkg", "", rBins, rLow, rHigh);
  THStack* sRapidityMC_allBkg = new THStack("sRapidityMC_allBkg","");

  TH1F* hVertexMC_allBkg     = new TH1F("hVertexMC_allBkg", "", 25, 0, 25);
  THStack* sVertexMC_allBkg     = new THStack("sVertexMC_allBkg","");

 
//-------------------------------------------------------
//event counting before scaling
  std::cout << "Event Counting before luminosity rescaling:" << std::endl;
  outLatex(hMassData, "data");
  outLatex(hMassMC_dymumu, "DYToMuMu");
  outLatex(hMassMC_ttbarjets, "TTJets");
  outLatex(hMassMC_ztautau, "Ztautau");
  outLatex(hMassMC_WW, "WW");
  outLatex(hMassMC_WZ, "WZ");
  outLatex(hMassMC_ZZ, "ZZ");
  outLatex(hMassMC_WJetsToLNu, "WJetsToLNu");
  outLatex(hMassMC_QCD, "QCD");

  
 

// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
  // Rescaling the MC to match the luminosity
  std::cout << "xsec (Z/gamma + Jets, M_ll>50 GeV): "   << xsec_dymumu   << " pb" << std::endl;
 
  hMassMC_dymumu    -> Scale(intLumi_invpb*xsec_dymumu/nMCEvDYToMuMu);
  hMassMC_ttbarjets -> Scale(intLumi_invpb*xsec_ttbarjets/nMCEvTTbarJets);
  hMassMC_ztautau   -> Scale(intLumi_invpb*xsec_ztautau/nMCEvZtautau);
  hMassMC_WW        -> Scale(intLumi_invpb*xsec_WW/nMCEvWW);
  hMassMC_WZ        -> Scale(intLumi_invpb*xsec_WZ/nMCEvWZ);
  hMassMC_ZZ        -> Scale(intLumi_invpb*xsec_ZZ/nMCEvZZ);
  hMassMC_WJetsToLNu  -> Scale(intLumi_invpb*xsec_WJetsToLNu/nMCEvWJetsToLNu);
  hMassMC_QCD  -> Scale(intLumi_invpb*xsec_QCD/nMCEvQCD);

  hPtMC_dymumu    -> Scale(intLumi_invpb*xsec_dymumu/nMCEvDYToMuMu);
  hPtMC_ttbarjets -> Scale(intLumi_invpb*xsec_ttbarjets/nMCEvTTbarJets);
  hPtMC_ztautau   -> Scale(intLumi_invpb*xsec_ztautau/nMCEvZtautau);
  hPtMC_WW        -> Scale(intLumi_invpb*xsec_WW/nMCEvWW);
  hPtMC_WZ        -> Scale(intLumi_invpb*xsec_WZ/nMCEvWZ);
  hPtMC_ZZ        -> Scale(intLumi_invpb*xsec_ZZ/nMCEvZZ);
  hPtMC_WJetsToLNu  -> Scale(intLumi_invpb*xsec_WJetsToLNu/nMCEvWJetsToLNu);
  hPtMC_QCD  -> Scale(intLumi_invpb*xsec_QCD/nMCEvQCD);

  hPt20BMC_dymumu    -> Scale(intLumi_invpb*xsec_dymumu/nMCEvDYToMuMu);
  hPt20BMC_ttbarjets -> Scale(intLumi_invpb*xsec_ttbarjets/nMCEvTTbarJets);
  hPt20BMC_ztautau   -> Scale(intLumi_invpb*xsec_ztautau/nMCEvZtautau);
  hPt20BMC_WW        -> Scale(intLumi_invpb*xsec_WW/nMCEvWW);
  hPt20BMC_WZ        -> Scale(intLumi_invpb*xsec_WZ/nMCEvWZ);
  hPt20BMC_ZZ        -> Scale(intLumi_invpb*xsec_ZZ/nMCEvZZ);
  hPt20BMC_WJetsToLNu  -> Scale(intLumi_invpb*xsec_WJetsToLNu/nMCEvWJetsToLNu);
  hPt20BMC_QCD  -> Scale(intLumi_invpb*xsec_QCD/nMCEvQCD);

  hRapidityMC_dymumu    -> Scale(intLumi_invpb*xsec_dymumu/nMCEvDYToMuMu);
  hRapidityMC_ttbarjets -> Scale(intLumi_invpb*xsec_ttbarjets/nMCEvTTbarJets);
  hRapidityMC_ztautau   -> Scale(intLumi_invpb*xsec_ztautau/nMCEvZtautau);
  hRapidityMC_WW        -> Scale(intLumi_invpb*xsec_WW/nMCEvWW);
  hRapidityMC_WZ        -> Scale(intLumi_invpb*xsec_WZ/nMCEvWZ);
  hRapidityMC_ZZ        -> Scale(intLumi_invpb*xsec_ZZ/nMCEvZZ);
  hRapidityMC_WJetsToLNu  -> Scale(intLumi_invpb*xsec_WJetsToLNu/nMCEvWJetsToLNu);
  hRapidityMC_QCD  -> Scale(intLumi_invpb*xsec_QCD/nMCEvQCD);

  hVertexMC_dymumu    -> Scale(intLumi_invpb*xsec_dymumu/nMCEvDYToMuMu);
  hVertexMC_ttbarjets -> Scale(intLumi_invpb*xsec_ttbarjets/nMCEvTTbarJets);
  hVertexMC_ztautau   -> Scale(intLumi_invpb*xsec_ztautau/nMCEvZtautau);
  hVertexMC_WW        -> Scale(intLumi_invpb*xsec_WW/nMCEvWW);
  hVertexMC_WZ        -> Scale(intLumi_invpb*xsec_WZ/nMCEvWZ);
  hVertexMC_ZZ        -> Scale(intLumi_invpb*xsec_ZZ/nMCEvZZ);
  hVertexMC_WJetsToLNu  -> Scale(intLumi_invpb*xsec_WJetsToLNu/nMCEvWJetsToLNu);
  hVertexMC_QCD  -> Scale(intLumi_invpb*xsec_QCD/nMCEvQCD);

/*  
    double bkgSum =  hMassMC_dyjets -> Integral()+
                  hMassMC_ttbarjets -> Integral()+
                  hMassMC_ztautau -> Integral()+
                  hMassMC_WW -> Integral()+
                  hMassMC_WZ -> Integral()+
                  hMassMC_ZZ -> Integral();

  //rescale to the number of events found in data
  double rescale = hMassData -> Integral()/bkgSum;
  
  hMassMC_dyjets    -> Scale(rescale);
  hMassMC_ttbarjets -> Scale(rescale);
  hMassMC_ztautau   -> Scale(rescale);
  hMassMC_WW        -> Scale(rescale);
  hMassMC_WZ        -> Scale(rescale);
  hMassMC_ZZ        -> Scale(rescale);
  // end rescaling
  
*/
  std::cout << "Event Counting after luminosity rescaling:" << std::endl;
  outLatex(hMassData, "data");
  outLatex(hMassMC_dymumu, "DYToMuMu");
  outLatex(hMassMC_ttbarjets, "TTJets");
  outLatex(hMassMC_ztautau, "Ztautau");
  outLatex(hMassMC_WW, "WW");
  outLatex(hMassMC_WZ, "WZ");
  outLatex(hMassMC_ZZ, "ZZ");
  outLatex(hMassMC_WJetsToLNu, "WJetsToLNu");
  outLatex(hMassMC_QCD, "QCD");

  

  // ------------------------------------------------------------------------------
  
  double entriesPerBinMass = (mHigh-mLow)/mBins;
  char yAxisNameMass[128];
  sprintf(yAxisNameMass,"Entries/%3.2f GeV/c^{2}",entriesPerBinMass);

  double entriesPerBinPt = (ptHigh-ptLow)/ptBins;
  char yAxisNamePt[128];
  sprintf(yAxisNamePt,"Entries/%3.2f GeV",entriesPerBinPt);

  double entriesPerBinRapidity = (rHigh-rLow)/rBins;
  char yAxisNameRapidity[128];
  sprintf(yAxisNameRapidity,"Entries/%3.2f",entriesPerBinRapidity);

  double entriesPerBinVertex = 1;
  char yAxisNameVertex[128];
  sprintf(yAxisNameVertex,"Entries/%3.2f",entriesPerBinVertex);

  hPtMC_EWKBkg->Add(hPtMC_ttbarjets); 	sPtMC_EWKBkg->Add(hPtMC_ttbarjets);
  hPtMC_EWKBkg->Add(hPtMC_ztautau); 	sPtMC_EWKBkg->Add(hPtMC_ztautau);
  hPtMC_EWKBkg->Add(hPtMC_WJetsToLNu);  sPtMC_EWKBkg->Add(hPtMC_WJetsToLNu);
  hPtMC_EWKBkg->Add(hPtMC_ZZ);          sPtMC_EWKBkg->Add(hPtMC_ZZ);
  hPtMC_EWKBkg->Add(hPtMC_WZ);          sPtMC_EWKBkg->Add(hPtMC_WZ);
  hPtMC_EWKBkg->Add(hPtMC_WW);          sPtMC_EWKBkg->Add(hPtMC_WW);

  // define the tlegend
  //TLegend* tl = new TLegend(0.71, 0.62, 0.90, 0.90);
  //TLegend* tl = new TLegend(0.60, 0.57, 0.86, 0.86);
  TLegend* tl = new TLegend(0.60, 0.47, 0.86, 0.76);
  tl->SetFillColor(0);
  tl->SetLineColor(0);
  tl->AddEntry(hPtData ,"Data","lpe");
  //tl->AddEntry(hPtMC_dymumu,"DYToMuMu","f");
  tl->AddEntry(hPtMC_dymumu,"DY#rightarrow #mu^{+}#mu^{-}","f");
  //tl->AddEntry(hPtMC_QCD," QCD","f");
  tl->AddEntry(hPtMC_EWKBkg,"EW+t#bar{t}","f");
  //tl->AddEntry(hPtMC_WW," WW","f");
  //tl->AddEntry(hPtMC_WZ," WZ","f");
  //tl->AddEntry(hPtMC_ZZ," ZZ","f");
  //tl->AddEntry(hPtMC_WJetsToLNu," WJets","f");
  //tl->AddEntry(hPtMC_ztautau  ,"DYToTauTau","f");
  //tl->AddEntry(hPtMC_ttbarjets," t #bar{t}","f");

  // define MC contribution color and style
  //mass
  hMassMC_dymumu -> SetFillColor(798);//yellow
  hMassMC_dymumu -> SetFillStyle(1001);

  hMassMC_ttbarjets -> SetFillColor(30);//green
  hMassMC_ttbarjets -> SetFillStyle(1001);

  hMassMC_ztautau -> SetFillColor(51);//violet
  hMassMC_ztautau -> SetFillStyle(1001);

  hMassMC_WW   -> SetFillColor(41);//brown
  hMassMC_WW   -> SetFillStyle(1001);

  hMassMC_WZ   -> SetFillColor(44);//
  hMassMC_WZ   -> SetFillStyle(1001);

  hMassMC_ZZ   -> SetFillColor(47);//
  hMassMC_ZZ   -> SetFillStyle(1001);

  hMassMC_WJetsToLNu ->SetFillColor(2); //red
  hMassMC_WJetsToLNu ->SetFillStyle(1001);

  hMassMC_QCD   -> SetFillColor(kMagenta+3);//
  hMassMC_QCD   -> SetFillStyle(1001);
  
    
//pt
  hPtMC_dymumu -> SetFillColor(798);//yellow
  hPtMC_dymumu -> SetFillStyle(1001);

  hPtMC_ttbarjets -> SetFillColor(30);//green
  hPtMC_ttbarjets -> SetFillStyle(1001);

  hPtMC_ztautau -> SetFillColor(51);//violet
  hPtMC_ztautau -> SetFillStyle(1001);

  hPtMC_WW   -> SetFillColor(41);//brown
  hPtMC_WW   -> SetFillStyle(1001);

  hPtMC_WZ   -> SetFillColor(44);//
  hPtMC_WZ   -> SetFillStyle(1001);

  hPtMC_ZZ   -> SetFillColor(47);//
  hPtMC_ZZ   -> SetFillStyle(1001);

  hPtMC_WJetsToLNu ->SetFillColor(2); //red
  hPtMC_WJetsToLNu ->SetFillStyle(1001);

  hPtMC_QCD   -> SetFillColor(kMagenta+3);//
  hPtMC_QCD   -> SetFillStyle(1001);

  hPtMC_EWKBkg   -> SetFillColor(kOrange+7);//
  hPtMC_EWKBkg   -> SetFillStyle(1001);

//pt 20 bin
  hPt20BMC_dymumu -> SetFillColor(798);//yellow
  hPt20BMC_dymumu -> SetFillStyle(1001);

  hPt20BMC_ttbarjets -> SetFillColor(30);//green
  hPt20BMC_ttbarjets -> SetFillStyle(1001);

  hPt20BMC_ztautau -> SetFillColor(51);//violet
  hPt20BMC_ztautau -> SetFillStyle(1001);

  hPt20BMC_WW   -> SetFillColor(41);//brown
  hPt20BMC_WW   -> SetFillStyle(1001);

  hPt20BMC_WZ   -> SetFillColor(44);//
  hPt20BMC_WZ   -> SetFillStyle(1001);

  hPt20BMC_ZZ   -> SetFillColor(47);//
  hPt20BMC_ZZ   -> SetFillStyle(1001);

  hPt20BMC_WJetsToLNu ->SetFillColor(2); //red
  hPt20BMC_WJetsToLNu ->SetFillStyle(1001);

  hPt20BMC_QCD   -> SetFillColor(kMagenta+3);//
  hPt20BMC_QCD   -> SetFillStyle(1001);

 //Rapidity
  hRapidityMC_dymumu -> SetFillColor(798);//yellow
  hRapidityMC_dymumu -> SetFillStyle(1001);

  hRapidityMC_ttbarjets -> SetFillColor(30);//green
  hRapidityMC_ttbarjets -> SetFillStyle(1001);

  hRapidityMC_ztautau -> SetFillColor(51);//violet
  hRapidityMC_ztautau -> SetFillStyle(1001);

  hRapidityMC_WW   -> SetFillColor(41);//brown
  hRapidityMC_WW   -> SetFillStyle(1001);

  hRapidityMC_WZ   -> SetFillColor(44);//
  hRapidityMC_WZ   -> SetFillStyle(1001);

  hRapidityMC_ZZ   -> SetFillColor(47);//
  hRapidityMC_ZZ   -> SetFillStyle(1001);

  hRapidityMC_WJetsToLNu ->SetFillColor(2); //red
  hRapidityMC_WJetsToLNu ->SetFillStyle(1001);

  hRapidityMC_QCD   -> SetFillColor(kMagenta+3);//
  hRapidityMC_QCD   -> SetFillStyle(1001);
 
 //Vertex
  hVertexMC_dymumu -> SetFillColor(798);//yellow
  hVertexMC_dymumu -> SetFillStyle(1001);
  
  hVertexMC_ttbarjets -> SetFillColor(30);//green
  hVertexMC_ttbarjets -> SetFillStyle(1001);

  hVertexMC_ztautau -> SetFillColor(51);//violet
  hVertexMC_ztautau -> SetFillStyle(1001);

  hVertexMC_WW   -> SetFillColor(41);//brown
  hVertexMC_WW   -> SetFillStyle(1001);

  hVertexMC_WZ   -> SetFillColor(44);//
  hVertexMC_WZ   -> SetFillStyle(1001);

  hVertexMC_ZZ   -> SetFillColor(47);//
  hVertexMC_ZZ   -> SetFillStyle(1001);

  hVertexMC_WJetsToLNu ->SetFillColor(2); //red
  hVertexMC_WJetsToLNu ->SetFillStyle(1001);

  hVertexMC_QCD   -> SetFillColor(kMagenta+3);//
  hVertexMC_QCD   -> SetFillStyle(1001);

    
  // ------------------------------------------------------------------------------
  // Mass Plot Linear Scale
  //TCanvas* cmass1 = new TCanvas("cmass1","",0,0,750,700);
  TCanvas* cmass1 = new TCanvas("cmass1","",800,800);
  cmass1->cd();
    cmass1 -> SetLogy();
  
//  gPad->SetLogy();         // Log plots on both axes.
//  gPad->SetLogx();

  hMassData -> SetTitle("CMS Preliminary, 18.4 pb^{-1} at #sqrt{s}=8 TeV");
  hMassData -> SetTitle("#int 18.4 pb^{-1} at #sqrt{s}=8 TeV");
  hMassData -> GetXaxis() -> SetTitle("Mass (#mu#mu) [GeV/c^{2}]");
  hMassData -> GetYaxis() -> SetTitle(yAxisNameMass);
  hMassData -> GetYaxis() -> SetTitleOffset(1.5);

    hMassData -> SetMarkerSize(1);
  hMassData -> SetMarkerStyle(20);
  hMassData -> SetMinimum(0.01);
  //hMassData -> Draw("pe");
  //sMassMC_allBkg -> Draw("same");
  //hMassData -> Draw("pe same");
  //hMassData -> Draw("AXIS same");

  zrap_Prelim(0.85,0.9,0.4,0.17);
  zrap_Lumi(0.35,0.89,36);

  //tl->Draw("same");
 
  hMassMC_allBkg->Add(hMassMC_ttbarjets);   sMassMC_allBkg->Add(hMassMC_ttbarjets);
  hMassMC_allBkg->Add(hMassMC_ztautau);     sMassMC_allBkg->Add(hMassMC_ztautau);
  hMassMC_allBkg->Add(hMassMC_WJetsToLNu);  sMassMC_allBkg->Add(hMassMC_WJetsToLNu);
  hMassMC_allBkg->Add(hMassMC_QCD);          sMassMC_allBkg->Add(hMassMC_QCD);
  hMassMC_allBkg->Add(hMassMC_ZZ);          sMassMC_allBkg->Add(hMassMC_ZZ);
  hMassMC_allBkg->Add(hMassMC_WZ);          sMassMC_allBkg->Add(hMassMC_WZ);
  hMassMC_allBkg->Add(hMassMC_WW);          sMassMC_allBkg->Add(hMassMC_WW); 
  hMassMC_allBkg->Add(hMassMC_dymumu);      sMassMC_allBkg->Add(hMassMC_dymumu);
  
  // Mass Plot With Pulls
  TCanvas* cmass2 = new TCanvas("cmass2","",700,0,750,700);
  //DrawWithRes(cmass2, cmspreliminary,"", hMassData, hMassMC_allBkg, sMassMC_allBkg, tl, true);

 // ------------------------------------------------------------------------------
// Pt with equal bin with 5GeV
 
  TCanvas* cpt20B1 = new TCanvas("cpt20B1","",0,0,600,500);
  cpt20B1->cd();
  cpt20B1 -> SetLogy();
  cpt20B1 -> SetTicky(1);
  cpt20B1 -> SetTickx(1);
  
  //hPt20BData -> GetXaxis() -> SetTitle("q_{T}  [GeV/c]");
  hPt20BData -> GetXaxis() -> SetTitle("p_{T}^{Z}  [GeV]");
  hPt20BData -> GetYaxis() -> SetTitle(yAxisNamePt);
  hPt20BData -> GetYaxis() -> SetTitleOffset(1.2);

  hPt20BData -> SetMarkerSize(0.7);
  hPt20BData -> SetMarkerStyle(20);
  hPt20BData -> SetMinimum(0.01);
  //hPt20BData -> Draw("pe");
  //sPt20BMC_allBkg -> Draw("same");
  //hPt20BData -> Draw("pe same");
  //hPt20BData -> Draw("AXIS same");

  //tl->Draw("same");

  hPt20BMC_allBkg->Add(hPt20BMC_ttbarjets);   sPt20BMC_allBkg->Add(hPt20BMC_ttbarjets);
  hPt20BMC_allBkg->Add(hPt20BMC_ztautau);     sPt20BMC_allBkg->Add(hPt20BMC_ztautau);
  hPt20BMC_allBkg->Add(hPt20BMC_WJetsToLNu);  sPt20BMC_allBkg->Add(hPt20BMC_WJetsToLNu);
  hPt20BMC_allBkg->Add(hPt20BMC_QCD);         sPt20BMC_allBkg->Add(hPt20BMC_QCD);
  hPt20BMC_allBkg->Add(hPt20BMC_ZZ);          sPt20BMC_allBkg->Add(hPt20BMC_ZZ);
  hPt20BMC_allBkg->Add(hPt20BMC_WZ);          sPt20BMC_allBkg->Add(hPt20BMC_WZ);
  hPt20BMC_allBkg->Add(hPt20BMC_WW);          sPt20BMC_allBkg->Add(hPt20BMC_WW);
  hPt20BMC_allBkg->Add(hPt20BMC_dymumu);      sPt20BMC_allBkg->Add(hPt20BMC_dymumu);
 
  
 // pt Plot With Pulls
//  TCanvas* cpt20B2 = new TCanvas("cpt20B2","",700,0,600,500);
//  DrawWithRes(cpt20B2, cmspreliminaryextd, hPt20BData, hPt20BMC_allBkg, sPt20BMC_allBkg, tl, true, true);

 // ------------------------------------------------------------------------------
// Pt with dynamic binning

  //TCanvas* cpt1 = new TCanvas("cpt1","",0,0,750,600);
  TCanvas* cpt1 = new TCanvas("cpt1","",0,0,800,800);
  cpt1->cd();
  cpt1 -> SetLogy();
  cpt1 -> SetTicky(1);
  cpt1 -> SetTickx(1);
  //hPtData -> SetTitle("#int 18.4 pb^{-1} at #sqrt{s}=8 TeV"); 
  //hPtData -> SetTitle("L = 18.4 pb^{-1}, #sqrt{s}=8 TeV"); 
  //hPtData -> SetTitle("18.4 pb^{-1} (8 TeV)"); 
  //hPtData -> SetTitleSize(0.0); 
  //hPtData -> GetXaxis() -> SetTitle("q_{T}  [GeV/c]");
  hPtData -> GetXaxis() -> SetLabelSize(0);
  hPtData -> GetXaxis() -> SetTitle("p_{T}^{Z}  [GeV]");
  hPtData -> GetYaxis() -> SetTitle("Events");
  //hPtData -> GetYaxis() -> SetTitleSize(0.05);
  hPtData -> GetYaxis() -> SetTitleSize(0.07);
  //hPtData -> GetYaxis() -> SetTitleOffset(1.0);
  hPtData -> GetYaxis() -> SetTitleOffset(0.8);
  hPtData -> GetYaxis() -> SetLabelSize(0.05);
 
  hPtData -> SetMarkerSize(1);
  hPtData -> SetMarkerStyle(20);
  hPtData -> SetMinimum(5*0.001);

  hPtData -> Draw("pe");
  sPtMC_allBkg -> Draw("same");
  hPtData -> Draw("pe same");   
  hPtData -> Draw("AXIS same");

  //tl->Draw("same");
  
  hPtMC_allBkg->Add(hPtMC_QCD);         sPtMC_allBkg->Add(hPtMC_QCD);
  hPtMC_allBkg->Add(hPtMC_EWKBkg);      sPtMC_allBkg->Add(hPtMC_EWKBkg);
  hPtMC_allBkg->Add(hPtMC_dymumu);      sPtMC_allBkg->Add(hPtMC_dymumu);

  /*
  hPtMC_allBkg->Add(hPtMC_ttbarjets);   sPtMC_allBkg->Add(hPtMC_ttbarjets);
  hPtMC_allBkg->Add(hPtMC_ztautau);     sPtMC_allBkg->Add(hPtMC_ztautau);
  hPtMC_allBkg->Add(hPtMC_WJetsToLNu);  sPtMC_allBkg->Add(hPtMC_WJetsToLNu);
  hPtMC_allBkg->Add(hPtMC_QCD);         sPtMC_allBkg->Add(hPtMC_QCD);
  hPtMC_allBkg->Add(hPtMC_ZZ);          sPtMC_allBkg->Add(hPtMC_ZZ);
  hPtMC_allBkg->Add(hPtMC_WZ);          sPtMC_allBkg->Add(hPtMC_WZ);
  hPtMC_allBkg->Add(hPtMC_WW);          sPtMC_allBkg->Add(hPtMC_WW);
  hPtMC_allBkg->Add(hPtMC_dymumu);      sPtMC_allBkg->Add(hPtMC_dymumu);
   */
 // pt Plot With Pulls
  TCanvas* cpt2 = new TCanvas("cpt2","",700,0,750,700);
    cpt2 -> SetTicky(1);
    cpt2 -> SetTickx(1);
  //DrawWithRes(cpt2, cmspreliminary, hPtData, hPtMC_allBkg, sPtMC_allBkg, tl);
  //DrawWithRes(cpt2, "#font[61]{CMS}", "", hPtData, hPtMC_allBkg, sPtMC_allBkg, tl);
  DrawWithRes(cpt2, "#font[61]{CMS}", "18.4 pb^{-1} (8 TeV)", hPtData, hPtMC_allBkg, sPtMC_allBkg, tl);

 // ------------------------------------------------------------------------------
// Rapidity plots

 TCanvas* cra1 = new TCanvas("cra1","",0,0,600,500);
  cra1->cd();
  cra1 -> SetLogy();
  cra1 -> SetTicky(1);
  cra1 -> SetTickx(1);

  hRapidityData -> GetXaxis() -> SetTitle("Y (#mu#mu) ");
  hRapidityData -> GetYaxis() -> SetTitle(yAxisNameRapidity);
  hRapidityData -> GetYaxis() -> SetTitleOffset(1.2);

  hRapidityData -> SetMarkerSize(0.7);
  hRapidityData -> SetMarkerStyle(20);
  hRapidityData -> SetMinimum(0.01);

  //hRapidityData -> Draw("pe");
  //sRapidityMC_allBkg -> Draw("same");
  //hRapidityData -> Draw("pe same");
  //hRapidityData -> Draw("AXIS same");

  //tl->Draw("same");

  hRapidityMC_allBkg->Add(hRapidityMC_ttbarjets);   sRapidityMC_allBkg->Add(hRapidityMC_ttbarjets);
  hRapidityMC_allBkg->Add(hRapidityMC_ztautau);     sRapidityMC_allBkg->Add(hRapidityMC_ztautau);
  hRapidityMC_allBkg->Add(hRapidityMC_WJetsToLNu);  sRapidityMC_allBkg->Add(hRapidityMC_WJetsToLNu);
  hRapidityMC_allBkg->Add(hRapidityMC_QCD);          sRapidityMC_allBkg->Add(hRapidityMC_QCD);
  hRapidityMC_allBkg->Add(hRapidityMC_ZZ);          sRapidityMC_allBkg->Add(hRapidityMC_ZZ);
  hRapidityMC_allBkg->Add(hRapidityMC_WZ);          sRapidityMC_allBkg->Add(hRapidityMC_WZ);
  hRapidityMC_allBkg->Add(hRapidityMC_WW);          sRapidityMC_allBkg->Add(hRapidityMC_WW);
  hRapidityMC_allBkg->Add(hRapidityMC_dymumu);      sRapidityMC_allBkg->Add(hRapidityMC_dymumu);

 // Rapidity Plot With Pulls
//  TCanvas* cra2 = new TCanvas("cra2","",700,0,600,500);
//  DrawWithRes(cra2, cmspreliminaryextd, hRapidityData, hRapidityMC_allBkg, sRapidityMC_allBkg, tl, true);

 // ------------------------------------------------------------------------------

// # of Vertices plots

 TCanvas* cv1 = new TCanvas("cv1","",0,0,600,500);
  cv1->cd();
  cv1 -> SetLogy();

  hVertexData -> SetTitle("#int 18.4 pb^{-1} at #sqrt{s}=8 TeV");
  hVertexData -> GetXaxis() -> SetTitle("number of Vertices ");
  hVertexData -> GetYaxis() -> SetTitle(yAxisNameVertex);
  hVertexData -> GetYaxis() -> SetTitleOffset(1.2);

  hVertexData -> SetMarkerSize(0.7);
  hVertexData -> SetMarkerStyle(20);
  hVertexData -> SetMinimum(0.01);

  //hVertexData -> Draw("pe");
  //sVertexMC_allBkg -> Draw("same");
  //hVertexData -> Draw("pe same");
  //hVertexData -> Draw("AXIS same");

  //tl->Draw("same");

  hVertexMC_allBkg->Add(hVertexMC_ttbarjets);   sVertexMC_allBkg->Add(hVertexMC_ttbarjets);
  hVertexMC_allBkg->Add(hVertexMC_ztautau);     sVertexMC_allBkg->Add(hVertexMC_ztautau);
  hVertexMC_allBkg->Add(hVertexMC_WJetsToLNu);  sVertexMC_allBkg->Add(hVertexMC_WJetsToLNu);
  hVertexMC_allBkg->Add(hVertexMC_QCD);          sVertexMC_allBkg->Add(hVertexMC_QCD);
  hVertexMC_allBkg->Add(hVertexMC_ZZ);          sVertexMC_allBkg->Add(hVertexMC_ZZ);
  hVertexMC_allBkg->Add(hVertexMC_WZ);          sVertexMC_allBkg->Add(hVertexMC_WZ);
  hVertexMC_allBkg->Add(hVertexMC_WW);          sVertexMC_allBkg->Add(hVertexMC_WW);
  hVertexMC_allBkg->Add(hVertexMC_dymumu);      sVertexMC_allBkg->Add(hVertexMC_dymumu);

 // Vertex Plot With Pulls
//  TCanvas* cv2 = new TCanvas("cv2","",700,0,750,700);
//    DrawWithRes(cv2, cmspreliminary, hVertexData, hVertexMC_allBkg, sVertexMC_allBkg, tl, true);
/*
  
  //lumi->SetNDC();
  //energy->SetNDC();
  //lumi  ->DrawLatex(0.35,0.70, lumiString);
  //energy->DrawLatex(0.35,0.60, energyString);
    
//  time (&end);
//  dif = difftime (end,start);
  
//  std::cout << "the macro took " << dif << " secs to terminate\n";

//  if (!isSave) return;

//  char mSave[128];
//  c1   ->SaveAs(png+"dimuonMass.png");
*/

  Hist_out->Write();
  Hist_out->Close();  
}
void fillHistos(TTree* tree, 
                TH1F *hMass, 
                TH1F *hPt,
                TH1F *hPt20B,
                TH1F *hRapidity,
                TH1F *hVertex,
                double mLow,
                double mHigh,
                double ptLow,
                double ptHigh,
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

       //if (vertexNormChiSquare > 10) continue; 
     if (angleDiMuons > TMath::Pi()-0.02) continue;

    // both muons tight
    //if( !isKinTight_2012WZ(reco1) || !isKinTight_2012WZ(reco2) ) continue; 
     if( !isKinTight_2012(reco1) || !isKinTight_2012(reco2) ) continue; 

     if( !reco1.isHltMatched[0] && !reco2.isHltMatched[0]) continue; 
   
    // remove duplicates
    pairOfInt runEvent(eventInfo.run,eventInfo.event);
    if( !uniqueEventsReco.insert( runEvent ).second ) continue;

    hMass     -> Fill( recoCandMass); 
    hPt       -> Fill( recoCandPt );
    hPt20B     -> Fill( recoCandPt );
    hRapidity -> Fill( recoCandY );
    hVertex   -> Fill( vertexInfo.nVertices );
    /*
    if ( recoCandPt > 200 ){
    
     cout  << "run = "   << eventInfo.run      << "  " 
           << "event = " << eventInfo.event    << "  "
           << "lumi = "  << eventInfo.lumi     << " "
           << "Mass(dimuon) = " << recoCandMass << " "
           << "Pt(dimuon) = " << recoCandPt     << std::endl; }
    */
      
  
  }


  std::cout << " DONE!" << std::endl;
  return;
}


