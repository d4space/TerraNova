/*
 *  unfolding.C
 *  Created by Joseph Gartner on 12/2/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 */

#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
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
#include<TMatrixD.h>

#include "Math/VectorUtil_Cint.h"
#include "Math/GenVector/LorentzVector.h"
#include <TLorentzVector.h>

#include "TF1.h"
#include <TGraphErrors.h>
#include<modifiedStyle.C>
#include<userStyle.C>

/*
// ROOFIT
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooPlot.h"

#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

using namespace RooFit ;
using namespace std;
*/

//const float PTbin[] = {20., 30., 40., 50., 60., 70., 100., 150., 300.};
//const float PTbin[] = {20., 30., 40., 45., 50., 60., 70., 100.};
const float PTbin[] = {20., 30., 35., 40., 45., 50., 60., 70., 100.}; //default
const float ETAbin[] = {0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
const float ZmassEtaBins[] = {-2.1,-1.4,-0.7,0.,0.7,1.4,2.1};
const float ZmassPtBins[] = {20,30,40,45};
const int NPThist = (sizeof(PTbin)/sizeof(float)-1);
const int NETAhist = (sizeof(ETAbin)/sizeof(float)-1);
const int Nhist = NPThist*NETAhist;
const int nEtaBins = 6;
const int nPtBins = 4;

const float PTbin_inv[] = {20., 30., 50., 100.};
const float ETAbin_inv[] = {0., 0.9, 2.1};
const int NPThist_inv = (sizeof(PTbin_inv)/sizeof(float)-1);
const int NETAhist_inv = (sizeof(ETAbin_inv)/sizeof(float)-1);
const int Nhist_inv = NPThist_inv*NETAhist_inv;

//MC weight to Data
Double_t evtWeight = 0.01781811133263710876;

//const float ScaleFactor = 1.;
const float ScaleFactor = 1.1;
// variation Asig2 by 3 sigma;
const float ScaleAsig2 = 0;// by default + 0 sigma of Asig2;
//const float ScaleAsig2 = 3.;// + 3 sigma of Asig2;
//const float ScaleAsig2 = -3.;// + 3 sigma of Asig2;
TString Extra = "";
//TString Extra = "_AsigPlus3sigma";
//TString Extra = "_AsigMinus3sigma";
//TString Extra = "_MATGRAPH";

Double_t MASS_MUON = 0.105658367;    //GeV/c2
Double_t mean[Nhist];
Double_t sig1[Nhist];
Double_t sig2[Nhist];
Double_t Asig2[Nhist]; 
// ---- CP error ----

TH1F* hmuonRes[Nhist];
TH1F* hmuonResSim[Nhist];
TH1F* hDimuonRes[Nhist];

TH1F* hDimuonMass[Nhist_inv];
TH1F* hDimuonMassSim[Nhist_inv];
TH1F* hDimuonMassDATA[Nhist_inv];

TH1F* h1_ZmassData    = new TH1F("h1_ZmassData", "", 20,80.,100.);
TH1F* h1_ZmassData_BB = new TH1F("h1_ZmassData_BB", "", 20,80.,100.);
TH1F* h1_ZmassData_BE = new TH1F("h1_ZmassData_BE", "", 20,80.,100.);
TH1F* h1_ZmassData_EE = new TH1F("h1_ZmassData_EE", "", 20,80.,100.);

TH1F* h1_ZmassMC    = new TH1F("h1_ZmassMC", "", 20,80.,100.);
TH1F* h1_ZmassMC_BB = new TH1F("h1_ZmassMC_BB", "", 20,80.,100.);
TH1F* h1_ZmassMC_BE = new TH1F("h1_ZmassMC_BE", "", 20,80.,100.);
TH1F* h1_ZmassMC_EE = new TH1F("h1_ZmassMC_EE", "", 20,80.,100.);

TH1F* h1_ZmassMCcorr    = new TH1F("h1_ZmassMCcorr", "", 20,80.,100.);
TH1F* h1_ZmassMCcorr_BB = new TH1F("h1_ZmassMCcorr_BB", "", 20,80.,100.);
TH1F* h1_ZmassMCcorr_BE = new TH1F("h1_ZmassMCcorr_BE", "", 20,80.,100.);
TH1F* h1_ZmassMCcorr_EE = new TH1F("h1_ZmassMCcorr_EE", "", 20,80.,100.);

TH1F* h1_ZmassData_muEtaP[nEtaBins];
TH1F* h1_ZmassData_muEtaM[nEtaBins];
TH1F* h1_ZmassMC_muEtaP[nEtaBins];
TH1F* h1_ZmassMC_muEtaM[nEtaBins];
TH1F* h1_ZmassMCcorr_muEtaP[nEtaBins];
TH1F* h1_ZmassMCcorr_muEtaM[nEtaBins];

TH1F* h1_ZmassData_Leading_muEtaP[nEtaBins];
TH1F* h1_ZmassData_Leading_muEtaM[nEtaBins];
TH1F* h1_ZmassMC_Leading_muEtaP[nEtaBins];
TH1F* h1_ZmassMC_Leading_muEtaM[nEtaBins];

TH1F* h1_ZmassData_Trailing_muEtaP[nEtaBins];
TH1F* h1_ZmassData_Trailing_muEtaM[nEtaBins];
TH1F* h1_ZmassMC_Trailing_muEtaP[nEtaBins];
TH1F* h1_ZmassMC_Trailing_muEtaM[nEtaBins];

TH1F* h1_ZmassData_muPtP[nPtBins];
TH1F* h1_ZmassData_muPtM[nPtBins];
TH1F* h1_ZmassMC_muPtP[nPtBins];
TH1F* h1_ZmassMC_muPtM[nPtBins];
TH1F* h1_ZmassMCcorr_muPtP[nPtBins];
TH1F* h1_ZmassMCcorr_muPtM[nPtBins];

// unfolding matrix

TH1F *raw1;
TH1F *test1;
TH1F *key1;
TH1F *key1_MC;
TH1F *ratTest;
TH1F *ratRaw;


TH2F *rMatrixVisualization;
TH2F *rInvMatrixVisualization;
TH2F *rPriMatrixVisualization;

TCanvas *_c1;
TCanvas *_c2;
TCanvas *_c3;
TCanvas *_c4;
TCanvas *_c5;
TLegend *_leg1;
TLegend *_leg2;
float binning[19] = {0.,2.5,5,7.5,10,12.5,15,17.5,20,30,40,50,70,90,110,150,190,250,600};
// to plot from 1. but fill from 0.:
float binning_histo[19] = {1.,2.5,5,7.5,10,12.5,15,17.5,20,30,40,50,70,90,110,150,190,250,600};
//float binning_histo[19];
//binning_histo[0] = 1.;
//for(int j=1; j<19; j++)
//{
//    binning_histo[j] = binning[j];
//}

float mapping[18][18];
float primeMapping[18][18];
float recoCounter[18];
float realCounter[18];

using namespace std;

// muon info
typedef struct {

  int isTracker;
  int isStandAlone;
  int isGlobal;

  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;

  float trkPt;
  float trkPtErr;
  float trketa;
  float trkPhi;

  float normChiSquare;
  float d0_BM;
  float dz_BM;

  float d0_PV;
  float dz_PV;

  int numPixelLayers;   //number of pixel layers with valid hits
  int numTrackerLayers; //number of tracker layers with valid hits 
  int numStripLayers;   //number of strip layers with valid hits

  float validFracTracker; //valid fraction of tracker hits

  int numValidMuonHits;
  int numValidPixelHits;
  int numValidTrackerHits;
  int numValidStripHits;
  int numSegmentMatches;
  int numOfMatchedStations;

  float trackIsoSumPt;
  float hcalIso;
  float ecalIso;
  float relCombIso;

  // PF information
  int isPFMuon;

  float pfPt;
  float pfEta;
  float pfPhi;

  float sumChargedHadronPtR03; // sum-pt of charged Hadron 
  float sumChargedParticlePtR03; // sum-pt of charged Particles(inludes e/mu) 
  float sumNeutralHadronEtR03;  // sum pt of neutral hadrons
  float sumPhotonEtR03;  // sum pt of PF photons
  float sumPUPtR03;  // sum pt of charged Particles not from PV  (for Pu corrections)

  float sumChargedHadronPtR04; 
  float sumChargedParticlePtR04;
  float sumNeutralHadronEtR04;
  float sumPhotonEtR04;
  float sumPUPtR04;

  int isHltMatched[2];

} _MuonInfo;

// true info
typedef struct {
  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;

} _TrueInfo;

void bookhistos();     // to book histograms
void printhistos();
Double_t BWnonrel(Double_t*, Double_t* );
Double_t Gauss(Double_t*, Double_t* );
Double_t FuncVoigtian(Double_t*, Double_t* );
Double_t FuncVoigtianBG(Double_t*, Double_t* );
Double_t DoubleGauss(Double_t*, Double_t* );
bool isLoose(_MuonInfo& muon);
bool isKinTight(_MuonInfo& muon);

void DrawWithRatio(TCanvas *canvas, char *cTitle,
    TH1F *hNum, TH1F *hDen);

void axis1F(TH1F  *histo,
    TAxis *xaxis,
    TAxis *yaxis,
    char  *xtitle,
    char  *ytitle);


void ZPtZmassMakeHisto(){


  gROOT->Clear();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  userStyle();
  //modifiedStyle();
  // ---- open the MC files ----
  TChain* treeMC = new TChain("tree");
  //POHEG+PYTHIA
  //treeMC -> AddFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/DYToMuMu/mergeFile/DYToMuMu_merged.root");
  treeMC -> AddFile("../../DYToMuMu_merged.root");
  //MATGRAPH+PYTHIA 
  //treeMC -> AddFile("/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/MC/ntuplesWITHdz/DYJetsToLL/mergeFile/DYJetsToLL_merged.root");
  // output file
  //TFile *theFile    =new TFile(Form("Zmass_unfoldingScaleFactor%1.2f.root", ScaleFactor), "RECREATE");

  bookhistos();// histo init 

  float tMass;
  float rMass;
  float rPt;
  float tPt;
  float rEta;
  float tEta;
  _MuonInfo reco1, reco2;
  _TrueInfo true1mu, true2mu;
  treeMC->SetBranchAddress("reco1",&reco1);
  treeMC->SetBranchAddress("reco2",&reco2);
  treeMC->SetBranchAddress("true1",&true1mu);
  treeMC->SetBranchAddress("true2",&true2mu);
  treeMC->SetBranchAddress("truePt",&tPt);
  treeMC->SetBranchAddress("trueEta",&tEta);
  treeMC->SetBranchAddress("trueMass",&tMass);
  treeMC->SetBranchAddress("recoCandMass",&rMass);
  treeMC->SetBranchAddress("recoCandPt",&rPt);
  treeMC->SetBranchAddress("recoCandEta",&rEta);
  cout << "Nevent to process = " << treeMC->GetEntries() << endl;

  int nbad_tpt = 0;
  int nbad_tMass = 0;
  for( int k=0; k<treeMC->GetEntries(); k++)
  //for( int k=0; k<1000; k++)
  {
    //process progress
    if(k!=0 && (k%10000)==0)
      std::cout << "- processing event " << k << "\r" << std::flush;

    treeMC->GetEntry(k);
    if (reco1.charge == reco2.charge) continue;
    if (rMass <  60) continue;
    if (rMass > 120) continue;

    //check that generator part is filled 
    if(nbad_tpt < 10 && (true1mu.pt < -900. || true2mu.pt < -900)){//rejection for MATGRAPH
      cout << "CHECK ntuple, possible problem with fill gen level while RECO level is fine: true1mu.pt = " << true1mu.pt << endl; 
      nbad_tpt++;
    }
    if(true1mu.pt < -900. || true2mu.pt < -900) continue;//rejection for MATGRAPH
    if(nbad_tMass < 10 && tMass < -900.){//rejection for MATGRAPH
      cout << "CHECK ntuple, possible EXTRA problem with fill gen level while RECO level is fine: tMass = " << tMass << endl; 
      nbad_tMass++;
    }
    if(tMass < -900.) continue;//rejection for MATGRAPH

    if( !isKinTight(reco1) || !isKinTight(reco2) ) continue;
    if (reco1.pt < 20 )      continue; // pt cut
    if (reco2.pt < 20 )      continue; // pt cut
    if ((reco1.trackIsoSumPt)/reco1.pt >=0.1) continue; // isolation
    if ((reco2.trackIsoSumPt)/reco2.pt >=0.1) continue; // isolation


    ROOT::Math::PtEtaPhiEVector MuTrue1(true1mu.pt, true1mu.eta, true1mu.phi, true1mu.pt); 
    ROOT::Math::PtEtaPhiEVector MuTrue2(true2mu.pt, true2mu.eta, true2mu.phi, true2mu.pt); 
    ROOT::Math::PtEtaPhiEVector MuReco1(reco1.pt, reco1.eta, reco1.phi, reco1.pt); 
    ROOT::Math::PtEtaPhiEVector MuReco2(reco2.pt, reco2.eta, reco2.phi, reco2.pt); 
    Float_t muonRes_t1r1 = -10.; 
    Float_t muonRes_t1r2 = -10.;
    Float_t muonRes_t2r1 = -10.; 
    Float_t muonRes_t2r2 = -10.;
    if(true1mu.pt > 0.) muonRes_t1r1 = (reco1.pt-true1mu.pt)/true1mu.pt;  
    if(true1mu.pt > 0.) muonRes_t1r2 = (reco2.pt-true1mu.pt)/true1mu.pt;  
    if(true2mu.pt > 0.) muonRes_t2r2 = (reco2.pt-true2mu.pt)/true2mu.pt;
    if(true2mu.pt > 0.) muonRes_t2r1 = (reco1.pt-true2mu.pt)/true2mu.pt;
    float deltaR_t1r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco1);  
    float deltaR_t1r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco2);  
    float deltaR_t2r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco1);  
    float deltaR_t2r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco2);  
    Float_t muonResCorr_t1r1 = muonRes_t1r1;
    Float_t muonResCorr_t2r2 = muonRes_t2r2;
    if(deltaR_t1r1 > deltaR_t1r2) muonResCorr_t1r1 = muonRes_t1r2;
    if(deltaR_t2r2 > deltaR_t2r1) muonResCorr_t2r2 = muonRes_t2r1;
    if(muonResCorr_t1r1 < -0.5 && true1mu.pt >= 70.) cout << "muonResCorr_t1r1 = " << muonResCorr_t1r1 << " true pt = " << true1mu.pt << " reco pt = " << reco1.pt 
      << " true pt2 = " << true2mu.pt << " reco pt2 = " << reco2.pt << " deltaR_t1r1 = " << deltaR_t1r1 << " deltaR_t1r2 = " << deltaR_t1r2 << endl;  

    for(int iPT = 0; iPT < NPThist; iPT++){
      for(int iETA = 0; iETA < NETAhist; iETA++){
	int iK = iPT + iETA*NPThist;
	if(true1mu.pt >= PTbin[iPT] && true1mu.pt < PTbin[iPT+1] 
	    && fabs(true1mu.eta) >= ETAbin[iETA] && fabs(true1mu.eta) < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t1r1);  
	if(true2mu.pt >= PTbin[iPT] && true2mu.pt < PTbin[iPT+1] 
	    && fabs(true2mu.eta) >= ETAbin[iETA] && fabs(true2mu.eta) < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t2r2);  
	if(tPt >= PTbin[iPT] && tPt < PTbin[iPT+1] 
	    && fabs(tEta) >= ETAbin[iETA] && fabs(tEta) < ETAbin[iETA+1]) hDimuonRes[iK] -> Fill( (rPt-tPt)/tPt );  
      }} 

    for(int iPT = 0; iPT < NPThist_inv; iPT++){
      for(int iETA = 0; iETA < NETAhist_inv; iETA++){
	int iK = iPT + iETA*NPThist_inv;
	if(rPt >= PTbin_inv[iPT] && rPt < PTbin_inv[iPT+1]
	    && fabs(rEta) >= ETAbin_inv[iETA] && fabs(rEta) < ETAbin_inv[iETA+1]) hDimuonMass[iK] -> Fill(rMass);
      }}

    h1_ZmassMC -> Fill(rMass,evtWeight); 
    if((fabs(reco1.eta) >= 0. && fabs(reco1.eta) < 1.) && (fabs(reco2.eta) >= 0. && fabs(reco2.eta) < 1.))
      h1_ZmassMC_BB -> Fill(rMass,evtWeight); 
    if((fabs(reco1.eta) >= 0. && fabs(reco1.eta) < 1.) && (fabs(reco2.eta) >= 1. && fabs(reco2.eta) < 2.4))
      h1_ZmassMC_BE -> Fill(rMass,evtWeight); 
    if((fabs(reco1.eta) >= 1. && fabs(reco1.eta) < 2.4) && (fabs(reco2.eta) >= 0. && fabs(reco2.eta) < 1.))
      h1_ZmassMC_BE -> Fill(rMass,evtWeight); 
    if((fabs(reco1.eta) >= 1. && fabs(reco1.eta) < 2.4) && (fabs(reco2.eta) >= 1. && fabs(reco2.eta) < 2.4))
      h1_ZmassMC_EE -> Fill(rMass,evtWeight); 

    for(int iBin=0;iBin<nEtaBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(reco1.charge>0)
      {
	if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_muEtaP[iBin]->Fill(rMass,evtWeight);
      }
      if(reco2.charge>0)
      {
	if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_muEtaP[iBin]->Fill(rMass,evtWeight);
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(reco1.charge<0)
      {
	if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_muEtaM[iBin]->Fill(rMass,evtWeight);
      }
      if(reco2.charge<0)
      {
	if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_muEtaM[iBin]->Fill(rMass,evtWeight);
      }

      //========================================
      //Fill leading-trailing leptons
      //========================================
      if(reco1.pt>reco2.pt)
      {
	if(reco1.charge>0) if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Leading_muEtaP[iBin]->Fill(rMass,evtWeight);
	if(reco1.charge<0) if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Leading_muEtaM[iBin]->Fill(rMass,evtWeight);
	
	if(reco2.charge>0) if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Trailing_muEtaP[iBin]->Fill(rMass,evtWeight);
	if(reco2.charge<0) if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Trailing_muEtaM[iBin]->Fill(rMass,evtWeight);
      }else if(reco1.pt<=reco2.pt){
	if(reco1.charge>0) if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Trailing_muEtaP[iBin]->Fill(rMass,evtWeight);
	if(reco1.charge<0) if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Trailing_muEtaM[iBin]->Fill(rMass,evtWeight);
	
	if(reco2.charge>0) if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Leading_muEtaP[iBin]->Fill(rMass,evtWeight);
	if(reco2.charge<0) if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMC_Leading_muEtaM[iBin]->Fill(rMass,evtWeight);
      }
    }

    for(int iBin=0;iBin<nPtBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(reco1.charge>0)
      {
	if(iBin<3) if(reco1.pt >= ZmassPtBins[iBin] && reco1.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtP[iBin]->Fill(rMass,evtWeight);
	if(iBin==3) if(reco1.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtP[iBin]->Fill(rMass,evtWeight);
      }
      if(reco2.charge>0)
      {
	if(iBin<3) if(reco2.pt >= ZmassPtBins[iBin] && reco2.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtP[iBin]->Fill(rMass,evtWeight);
	if(iBin==3) if(reco2.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtP[iBin]->Fill(rMass,evtWeight);
      }
      //========================================
      //Check Charge Minus
      //========================================
      if(reco1.charge<0)
      {
	if(iBin<3) if(reco1.pt >= ZmassPtBins[iBin] && reco1.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtM[iBin]->Fill(rMass,evtWeight);
	if(iBin==3) if(reco1.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtM[iBin]->Fill(rMass,evtWeight);
      }
      if(reco2.charge<0)
      {
	if(iBin<3) if(reco2.pt >= ZmassPtBins[iBin] && reco2.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtM[iBin]->Fill(rMass,evtWeight);
	if(iBin==3) if(reco2.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMC_muPtM[iBin]->Fill(rMass,evtWeight);
      }
    }
  }

  TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 5);
  TF1* fitDoubleGauss2 = new TF1("fitDoubleGauss2", DoubleGauss, -0.1, 0.1, 5);
  for(int iPT = 0; iPT < NPThist; iPT++){
    for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      // 
      fitDoubleGauss->SetParameters(0., 0.011, 0.041, 0.25, 2300.);
      if(PTbin[iPT] > 45 && ETAbin[iETA] > 0.9) fitDoubleGauss->SetParameters(0., 0.019, 0.16, 0.03, 320.); 
      fitDoubleGauss->FixParameter(0,0.);
      fitDoubleGauss->SetParLimits(1, 0.005, 0.05);//restrict sigma1
      fitDoubleGauss->SetParName(0,"mean");
      fitDoubleGauss->SetParName(1,"sig1");
      fitDoubleGauss->SetParName(2,"Asig2");
      fitDoubleGauss->SetParName(3,"sig2");
      fitDoubleGauss->SetParName(4,"Norm");

      //hmuonRes[iK] -> Fit(fitDoubleGauss,"RLE");
      hmuonRes[iK] -> Fit(fitDoubleGauss,"RLE");

      fitDoubleGauss2->SetParameters(fitDoubleGauss->GetParameter(0),fitDoubleGauss->GetParameter(1),fitDoubleGauss->GetParameter(2),fitDoubleGauss->GetParameter(3),fitDoubleGauss->GetParameter(4));
      if(PTbin[iPT] > 45 && ETAbin[iETA] > 0.9) fitDoubleGauss->SetParameters(0., 0.019, 0.16, 0.03, 320.);
      fitDoubleGauss2->FixParameter(0,0.);
      fitDoubleGauss2->SetParLimits(1, 0.005, 0.05);//restrict sigma1
      fitDoubleGauss2->SetParName(0,"mean");
      fitDoubleGauss2->SetParName(1,"sig1");
      fitDoubleGauss2->SetParName(2,"Asig2");
      fitDoubleGauss2->SetParName(3,"sig2");
      fitDoubleGauss2->SetParName(4,"Norm");


      hmuonRes[iK] -> Fit(fitDoubleGauss2,"RLE");
      mean[iK] = fitDoubleGauss2->GetParameter(0);  
      sig1[iK] = fitDoubleGauss2->GetParameter(1);  
      sig2[iK] = fitDoubleGauss2->GetParameter(3);  
      Asig2[iK] = fitDoubleGauss2->GetParameter(2) + ScaleAsig2*fitDoubleGauss2->GetParError(2);  
      if(Asig2[iK] < 0.) Asig2[iK] = 0.;  
      //
    }}        

  //////////////////////////////////////////////     
  // unfolding matrix

  for(int j=0; j<18; j++)
    for(int k=0; k<18; k++)
    {
      mapping[j][k]=0;
      primeMapping[j][k]=0;
      recoCounter[k]=0;
      realCounter[k]=0;
    }
  int gPInd=-1;
  int rPind=-1;


  //simulate muon resolution

  for( int k=0; k<treeMC->GetEntries(); k++)
  //for( int k=0; k<1000; k++)
  {
    //process progress
    if(k!=0 && (k%10000)==0)
      std::cout << "- processing Simulate muon resolution event " << k << "\r" << std::flush;

    treeMC->GetEntry(k);
    if(true1mu.pt < -900. || true2mu.pt < -900) continue;//rejection for MATGRAPH
    if (reco1.charge == reco2.charge) continue;
    // mass cut before resimulation pt
    if (rMass < 40) continue;

    if( !isKinTight(reco1) || !isKinTight(reco2) ) continue;
    // find correct bin and make resolution
    int iK_cand1 = -10;
    int iK_cand2 = -10;
    Double_t pt1true = true1mu.pt;
    Double_t pt2true = true2mu.pt;
    Double_t eta1true = fabs(true1mu.eta);
    Double_t eta2true = fabs(true2mu.eta);
    // to resimulate all gen muon
    if(pt1true >= PTbin[NPThist]) pt1true = PTbin[NPThist]-0.5;
    if(pt2true >= PTbin[NPThist]) pt2true = PTbin[NPThist]-0.5;
    if(eta1true >= ETAbin[NETAhist]) eta1true = ETAbin[NETAhist]-0.001;
    if(eta2true >= ETAbin[NETAhist]) eta2true = ETAbin[NETAhist]-0.001;
    //
    Double_t pt1sim = 0.;
    Double_t pt2sim = 0.;
    Double_t pt1simgen = 0.;
    Double_t pt2simgen = 0.;
    for(int iPT = 0; iPT < NPThist; iPT++){
      for(int iETA = 0; iETA < NETAhist; iETA++){
	int iK = iPT + iETA*NPThist;
	if(pt1true >= PTbin[iPT] && pt1true < PTbin[iPT+1]
	    && fabs(eta1true) >= ETAbin[iETA] && fabs(eta1true) < ETAbin[iETA+1]) iK_cand1 = iK;
	if(pt2true >= PTbin[iPT] && pt2true < PTbin[iPT+1]
	    && fabs(eta2true) >= ETAbin[iETA] && fabs(eta2true) < ETAbin[iETA+1]) iK_cand2 = iK;
      }}
    //fitDoubleGauss->SetParameters(0,0.);
    //fitDoubleGauss->SetParameter(3,0.);
    fitDoubleGauss->SetParameter(4,1.);
    if(iK_cand1 > -1){
      fitDoubleGauss->SetParameter(0,mean[iK_cand1]);
      fitDoubleGauss->SetParameter(1,ScaleFactor*sig1[iK_cand1]);
      fitDoubleGauss->SetParameter(2,Asig2[iK_cand1]);
      fitDoubleGauss->SetParameter(3,ScaleFactor*sig2[iK_cand1]);
      //Double_t resSim = fitDoubleGauss->GetRandom(-0.15, 0.15);
      Double_t resSim = fitDoubleGauss->GetRandom();
      pt1simgen = true1mu.pt*(1+resSim); 
      //hmuonResSim[iK_cand1] -> Fill(resSim);
    }
    if(iK_cand2 > -1){
      fitDoubleGauss->SetParameter(0,mean[iK_cand1]);
      fitDoubleGauss->SetParameter(1,ScaleFactor*sig1[iK_cand2]);
      fitDoubleGauss->SetParameter(2,Asig2[iK_cand2]);
      fitDoubleGauss->SetParameter(3,ScaleFactor*sig2[iK_cand2]);
      //Double_t resSim = fitDoubleGauss->GetRandom(-0.15, 0.15); 
      Double_t resSim = fitDoubleGauss->GetRandom(); 
      pt2simgen = true2mu.pt*(1+resSim); 
      //hmuonResSim[iK_cand2] -> Fill(resSim);
    }
    // end: find correct bin and make resolution
    //create matching true and reco and sim reco 
    TLorentzVector MuTrue1, MuTrue2, MuReco1, MuReco2;
    MuTrue1.SetPtEtaPhiM(true1mu.pt, true1mu.eta, true1mu.phi, MASS_MUON);
    MuTrue2.SetPtEtaPhiM(true2mu.pt, true2mu.eta, true2mu.phi, MASS_MUON);
    MuReco1.SetPtEtaPhiM(reco1.pt, reco1.eta, reco1.phi, MASS_MUON);
    MuReco2.SetPtEtaPhiM(reco2.pt, reco2.eta, reco2.phi, MASS_MUON);
    //cout << "true1mu.pt =" << true1mu.pt << " True1:"  << endl;
    //MuTrue1.Print();
    //cout << "true2mu.pt =" << true2mu.pt <<" True2:" << endl;
    //MuTrue2.Print();
    //cout << "Reco1:" << endl;
    //MuReco1.Print();
    //cout << "Reco2:" << endl;
    //MuReco2.Print();

    Float_t muonRes_t1r1 = -10.;
    Float_t muonRes_t1r2 = -10.;
    Float_t muonRes_t2r1 = -10.;
    Float_t muonRes_t2r2 = -10.;
    if(true1mu.pt > 0.) muonRes_t1r1 = (reco1.pt-true1mu.pt)/true1mu.pt;
    if(true1mu.pt > 0.) muonRes_t1r2 = (reco2.pt-true1mu.pt)/true1mu.pt;
    if(true2mu.pt > 0.) muonRes_t2r2 = (reco2.pt-true2mu.pt)/true2mu.pt;
    if(true2mu.pt > 0.) muonRes_t2r1 = (reco1.pt-true2mu.pt)/true2mu.pt;
    float deltaR_t1r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco1);
    float deltaR_t1r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco2);
    float deltaR_t2r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco1);
    float deltaR_t2r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco2);
    Float_t muonResCorr_t1r1 = muonRes_t1r1;
    Float_t muonResCorr_t2r2 = muonRes_t2r2;
    pt1sim = pt1simgen;
    pt2sim = pt2simgen;
    if(deltaR_t1r1 > deltaR_t1r2) {
      muonResCorr_t1r1 = muonRes_t1r2;
      pt1sim = pt2simgen;
    }   
    if(deltaR_t2r2 > deltaR_t2r1) {
      muonResCorr_t2r2 = muonRes_t2r1;
      pt2sim = pt1simgen;
    }
    //end: create matching true and reco 
    // check extra tight selection
    if (pt1sim < 20 )      continue; // pt cut
    if (pt2sim < 20 )      continue; // pt cut
    if ((reco1.trackIsoSumPt)/pt1sim >=0.1) continue; // isolation
    if ((reco2.trackIsoSumPt)/pt2sim >=0.1) continue; // isolation
    //end: check extra tight selection
    // creat sim reco muons:
    TLorentzVector MuReco1sim, MuReco2sim;
    MuReco1sim.SetPtEtaPhiM(pt1sim, reco1.eta, reco1.phi, MASS_MUON);
    MuReco2sim.SetPtEtaPhiM(pt2sim, reco2.eta, reco2.phi, MASS_MUON);

    TLorentzVector MuRecoCandSim = MuReco1sim + MuReco2sim;
    Double_t rMassSim = MuRecoCandSim.M(); 
    //Double_t rMassSim = rMass; // old reco mass, shouldn't use if sim pt for cross-check only   
    Double_t rPtSim = MuRecoCandSim.Pt();
    Double_t rEtaSim = MuRecoCandSim.Eta();
    //cout << "rPtSim = " << rPtSim << " tPt = " << tPt <<endl;
    //recalculate new resolution
    if(true1mu.pt > 0.) muonRes_t1r1 = (pt1sim-true1mu.pt)/true1mu.pt;
    if(true1mu.pt > 0.) muonRes_t1r2 = (pt2sim-true1mu.pt)/true1mu.pt;
    if(true2mu.pt > 0.) muonRes_t2r2 = (pt2sim-true2mu.pt)/true2mu.pt;
    if(true2mu.pt > 0.) muonRes_t2r1 = (pt1sim-true2mu.pt)/true2mu.pt;
    muonResCorr_t1r1 = muonRes_t1r1;
    muonResCorr_t2r2 = muonRes_t2r2;
    if(deltaR_t1r1 > deltaR_t1r2) {
      muonResCorr_t1r1 = muonRes_t1r2;
    }
    if(deltaR_t2r2 > deltaR_t2r1) {
      muonResCorr_t2r2 = muonRes_t2r1;
    }
    //end: recalculate new resolution
    // mass cut:
    if (rMassSim <  60) continue;
    if (rMassSim > 120) continue;
    key1_MC->Fill(tPt);


    for(int iPT = 0; iPT < NPThist; iPT++){
      for(int iETA = 0; iETA < NETAhist; iETA++){
	int iK = iPT + iETA*NPThist;
	if(true1mu.pt >= PTbin[iPT] && true1mu.pt < PTbin[iPT+1]
	    && fabs(true1mu.eta) >= ETAbin[iETA] && fabs(true1mu.eta) < ETAbin[iETA+1]) hmuonResSim[iK] -> Fill(muonResCorr_t1r1);
	if(true2mu.pt >= PTbin[iPT] && true2mu.pt < PTbin[iPT+1]
	    && fabs(true2mu.eta) >= ETAbin[iETA] && fabs(true2mu.eta) < ETAbin[iETA+1]) hmuonResSim[iK] -> Fill(muonResCorr_t2r2);
      }}


    for(int iPT = 0; iPT < NPThist_inv; iPT++){
      for(int iETA = 0; iETA < NETAhist_inv; iETA++){
	int iK = iPT + iETA*NPThist_inv;
	if(rPtSim >= PTbin_inv[iPT] && rPtSim < PTbin_inv[iPT+1]
	    && fabs(rEtaSim) >= ETAbin_inv[iETA] && fabs(rEtaSim) < ETAbin_inv[iETA+1]) hDimuonMassSim[iK] -> Fill(rMassSim);
      }}

    h1_ZmassMCcorr -> Fill(rMassSim,evtWeight); 
    if((fabs(reco1.eta) >= 0. && fabs(reco1.eta) < 1.) && (fabs(reco2.eta) >= 0. && fabs(reco2.eta) < 1.))
      h1_ZmassMCcorr_BB -> Fill(rMassSim,evtWeight); 
    if((fabs(reco1.eta) >= 0. && fabs(reco1.eta) < 1.) && (fabs(reco2.eta) >= 1. && fabs(reco2.eta) < 2.4))
      h1_ZmassMCcorr_BE -> Fill(rMassSim,evtWeight); 
    if((fabs(reco1.eta) >= 1. && fabs(reco1.eta) < 2.4) && (fabs(reco2.eta) >= 0. && fabs(reco2.eta) < 1.))
      h1_ZmassMCcorr_BE -> Fill(rMassSim,evtWeight); 
    if((fabs(reco1.eta) >= 1. && fabs(reco1.eta) < 2.4) && (fabs(reco2.eta) >= 1. && fabs(reco2.eta) < 2.4))
      h1_ZmassMCcorr_EE -> Fill(rMassSim,evtWeight);

    //if((fabs(reco1.eta) < 0.9 || fabs(reco1.eta) > 1.2) && (fabs(reco2.eta) < 0.9 || fabs(reco2.eta) > 1.2))
    for(int iBin=0;iBin<nEtaBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(reco1.charge>0)
      {
	if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMCcorr_muEtaP[iBin]->Fill(rMassSim,evtWeight);
      }
      if(reco2.charge>0)
      {
	if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMCcorr_muEtaP[iBin]->Fill(rMassSim,evtWeight);
      }
      //========================================
      //Check Charge Minus
      //========================================
      if(reco1.charge<0)
      {
	if(reco1.eta >= ZmassEtaBins[iBin] && reco1.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMCcorr_muEtaM[iBin]->Fill(rMassSim,evtWeight);
      }
      if(reco2.charge<0)
      {
	if(reco2.eta >= ZmassEtaBins[iBin] && reco2.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassMCcorr_muEtaM[iBin]->Fill(rMassSim,evtWeight);
      }
    }

    for(int iBin=0;iBin<nPtBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(reco1.charge>0)
      {
	if(iBin<3) if(reco1.pt >= ZmassPtBins[iBin] && reco1.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtP[iBin]->Fill(rMassSim,evtWeight);
	if(iBin==3) if(reco1.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtP[iBin]->Fill(rMassSim,evtWeight);
      }
      if(reco2.charge>0)
      {
	if(iBin<3) if(reco2.pt >= ZmassPtBins[iBin] && reco2.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtP[iBin]->Fill(rMassSim,evtWeight);
	if(iBin==3) if(reco2.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtP[iBin]->Fill(rMassSim,evtWeight);
      }
      //========================================
      //Check Charge Minus
      //========================================
      if(reco1.charge<0)
      {
	if(iBin<3) if(reco1.pt >= ZmassPtBins[iBin] && reco1.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtM[iBin]->Fill(rMassSim,evtWeight);
	if(iBin==3) if(reco1.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtM[iBin]->Fill(rMassSim,evtWeight);
      }
      if(reco2.charge<0)
      {
	if(iBin<3) if(reco2.pt >= ZmassPtBins[iBin] && reco2.pt < ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtM[iBin]->Fill(rMassSim,evtWeight);
	if(iBin==3) if(reco2.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassMCcorr_muPtM[iBin]->Fill(rMassSim,evtWeight);
      }
    }
    //========================================

    // unfolding matrix:
    float rPt_corr = MuRecoCandSim.Pt();
    if(rPt_corr < binning_histo[0] && rPt_corr >= 0) rPt_corr = binning_histo[0];
    float tPt_corr = tPt;
    if(tPt_corr < binning_histo[0] && tPt_corr >= 0) tPt_corr = binning_histo[0];
    raw1->Fill(rPt_corr);
    key1->Fill(tPt_corr);
    for(int jjj=0; jjj<18; jjj++)
    {
      if( (rPt >= binning[jjj]) && (rPt < binning[jjj+1]) ) rPind = jjj;
      if( (tPt >= binning[jjj]) && (tPt < binning[jjj+1]) ) gPInd = jjj;
    }
    recoCounter[rPind]++;
    realCounter[gPInd]++;
    mapping[rPind][gPInd]++;
    primeMapping[gPInd][rPind]++;
    // end: unfolding matrix:

  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  TChain* treeData = new TChain("tree");
  //treeData -> AddFile("/data/uftrig01b/digiovan/root/CMSSW_5_2_5_ecalpatch1/NtuplesDataSingleMuRun2012A-23May2012-v2/mergedFiles/SingleMuRun2012A-23May2012-v2_190949.root");
  //treeData -> AddFile("/data/uftrig01b/digiovan/root/CMSSW_5_2_5_ecalpatch1/NtuplesDataSingleMuRun2012A-PromptReco-v1_LowPU/mergedFiles/SingleMuRun2012A-PromptReco-v1_191090_191367_193112_193116.root");
  treeData -> AddFile("../../fullselection/data2012A_LowPU.root");


  float rMass_data;
  float rPt_data;
  float rEta_data;
  _MuonInfo reco1_data, reco2_data;
  treeData->SetBranchAddress("reco1",&reco1_data);
  treeData->SetBranchAddress("reco2",&reco2_data);
  treeData->SetBranchAddress("recoCandMass",&rMass_data);
  treeData->SetBranchAddress("recoCandPt",&rPt_data);
  treeData->SetBranchAddress("recoCandEta",&rEta_data);
  cout << "Nevent DATA to process = " << treeData->GetEntries() << endl;

  for( int k=0; k<treeData->GetEntries(); k++)
  //for( int k=0; k<1000; k++)
  {
    //process progress
    if(k!=0 && (k%10000)==0)
      std::cout << "- processing DATA event " << k << "\r" << std::flush;

    treeData->GetEntry(k);
    if (reco1_data.charge == reco2_data.charge) continue;
    if (rMass_data <  60) continue;
    if (rMass_data > 120) continue;
    //cout << "Loose selection event = " << k << " reco1.pt = " << reco1.pt << " reco1.eta = " << reco1.eta << endl;
    //     << " reco1.numValidTrackerHits = " reco1.numValidTrackerHits << " reco1.trackIsoSumPt = " << reco1.trackIsoSumPt 
    //     << "reco1.d0 = " << reco1.d0 << endl;
    //if( !isLoose(reco1)    || !isLoose(reco2)    ) continue; // both loose

    //cout << "any selection event = " << k << " reco1_data.charge = " << reco1_data.charge << " rMass_data = " << rMass_data << endl;
    if( !isKinTight(reco1_data) || !isKinTight(reco2_data) ) continue;
    if (reco1_data.pt < 20 )      continue; // pt cut
    if (reco2_data.pt < 20 )      continue; // pt cut
    if ((reco1_data.trackIsoSumPt)/reco1_data.pt >=0.1) continue; // isolation
    if ((reco2_data.trackIsoSumPt)/reco2_data.pt >=0.1) continue; // isolation

    //cout << "pass selection event = " << k << " rEta_data = " << rEta_data << " rPt_data = " << rPt_data << " rMass_data = " << rMass_data << endl;

    for(int iPT = 0; iPT < NPThist_inv; iPT++){
      for(int iETA = 0; iETA < NETAhist_inv; iETA++){
	int iK = iPT + iETA*NPThist_inv;
	if(rPt_data >= PTbin_inv[iPT] && rPt_data < PTbin_inv[iPT+1]
	    && fabs(rEta_data) >= ETAbin_inv[iETA] && fabs(rEta_data) < ETAbin_inv[iETA+1]) hDimuonMassDATA[iK] -> Fill(rMass_data );
      }}

    h1_ZmassData -> Fill(rMass_data); 
    if((fabs(reco1_data.eta) >= 0. && fabs(reco1_data.eta) < 1.) && (fabs(reco2_data.eta) >= 0. && fabs(reco2_data.eta) < 1.))
      h1_ZmassData_BB -> Fill(rMass_data); 
    if((fabs(reco1_data.eta) >= 0. && fabs(reco1_data.eta) < 1.) && (fabs(reco2_data.eta) >= 1. && fabs(reco2_data.eta) < 2.4))
      h1_ZmassData_BE -> Fill(rMass_data); 
    if((fabs(reco1_data.eta) >= 1. && fabs(reco1_data.eta) < 2.4) && (fabs(reco2_data.eta) >= 0. && fabs(reco2_data.eta) < 1.))
      h1_ZmassData_BE -> Fill(rMass_data); 
    if((fabs(reco1_data.eta) >= 1. && fabs(reco1_data.eta) < 2.4) && (fabs(reco2_data.eta) >= 1. && fabs(reco2_data.eta) < 2.4))
      h1_ZmassData_EE -> Fill(rMass_data); 

    //if((fabs(reco1_data.eta) < 0.9 || fabs(reco1_data.eta) > 1.2) && (fabs(reco2_data.eta) < 0.9 || fabs(reco2_data.eta) > 1.2))
    
    //cout<<"Before, Lept1(Lept2): "<<reco1_data.pt<<"("<<reco2_data.pt<<")\t"<<reco1_data.eta<<"("<<reco2_data.eta<<")\t"<<reco1_data.charge<<"("<<reco2_data.charge<<")"<<endl;
    for(int iBin=0;iBin<nEtaBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(reco1_data.charge>0)
      {
	if(reco1_data.eta >= ZmassEtaBins[iBin] && reco1_data.eta < ZmassEtaBins[iBin+1])
	{
	  h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //if(reco1_data.pt > reco2_data.pt)
	  //{
	  //  if(reco1_data.eta < 0.35)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //  if(reco1_data.eta > 0.7 && reco1_data.eta < 1.4)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //  if(reco1_data.eta > 1.75)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //}

	  //if(reco1_data.pt <= reco2_data.pt)
	  //{
	  //  if(reco1_data.eta < 1.05)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //  if(reco1_data.eta > 1.4 && reco1_data.eta < 1.75)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //}
	}
      }
      if(reco2_data.charge>0)
      {
	if(reco2_data.eta >= ZmassEtaBins[iBin] && reco2_data.eta < ZmassEtaBins[iBin+1])
	{
	  h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //if(reco2_data.pt > reco1_data.pt)
	  //{
	  //  if(reco2_data.eta < 0.35)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //  if(reco2_data.eta > 0.7 && reco2_data.eta < 1.4)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //  if(reco2_data.eta > 1.75)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //}

	  //if(reco2_data.pt <= reco1_data.pt)
	  //{
	  //  if(reco2_data.eta < 1.05)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //  if(reco2_data.eta > 1.4 && reco2_data.eta < 1.75)
	  //    h1_ZmassData_muEtaP[iBin]->Fill(rMass_data);
	  //}
	}
      }
      //========================================
      //Check Charge Minus
      //========================================
      if(reco1_data.charge<0)
      {
	if(reco1_data.eta >= ZmassEtaBins[iBin] && reco1_data.eta < ZmassEtaBins[iBin+1])
	{
	  h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //if(reco1_data.pt > reco2_data.pt)
	  //{
	  //  if(reco1_data.eta < -1.4 && reco1_data.eta > -1.05)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //}
	  //    
	  //if(reco1_data.pt <= reco2_data.pt)
	  //{
	  //  if(reco1_data.eta < -1.05)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //  if(reco1_data.eta > -0.7 && reco1_data.eta < 1.05)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //  if(reco1_data.eta > 1.4)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //}
	}
      }

      if(reco2_data.charge<0)
      {
	if(reco2_data.eta >= ZmassEtaBins[iBin] && reco2_data.eta < ZmassEtaBins[iBin+1])
	{
	  h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //if(reco2_data.pt > reco1_data.pt)
	  //{
	  //  if(reco2_data.eta < -1.4 && reco2_data.eta > -1.05)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //}
	  //    
	  //if(reco2_data.pt <= reco1_data.pt)
	  //{
	  //  if(reco2_data.eta < -1.05)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //  if(reco2_data.eta > -0.7 && reco2_data.eta < 1.05)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //  if(reco2_data.eta > 1.4)
	  //    h1_ZmassData_muEtaM[iBin]->Fill(rMass_data);
	  //}
	}
      }

      //========================================
      //Fill leading-trailing leptons
      //========================================
      if(reco1_data.pt>reco2_data.pt)
      {
	if(reco1_data.charge>0) if(reco1_data.eta >= ZmassEtaBins[iBin] && reco1_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Leading_muEtaP[iBin]->Fill(rMass_data);
	if(reco1_data.charge<0) if(reco1_data.eta >= ZmassEtaBins[iBin] && reco1_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Leading_muEtaM[iBin]->Fill(rMass_data);
	
	if(reco2_data.charge>0) if(reco2_data.eta >= ZmassEtaBins[iBin] && reco2_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Trailing_muEtaP[iBin]->Fill(rMass_data);
	if(reco2_data.charge<0) if(reco2_data.eta >= ZmassEtaBins[iBin] && reco2_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Trailing_muEtaM[iBin]->Fill(rMass_data);
      }else if(reco1_data.pt<=reco2_data.pt){
	if(reco1_data.charge>0) if(reco1_data.eta >= ZmassEtaBins[iBin] && reco1_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Trailing_muEtaP[iBin]->Fill(rMass_data);
	if(reco1_data.charge<0) if(reco1_data.eta >= ZmassEtaBins[iBin] && reco1_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Trailing_muEtaM[iBin]->Fill(rMass_data);
	
	if(reco2_data.charge>0) if(reco2_data.eta >= ZmassEtaBins[iBin] && reco2_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Leading_muEtaP[iBin]->Fill(rMass_data);
	if(reco2_data.charge<0) if(reco2_data.eta >= ZmassEtaBins[iBin] && reco2_data.eta < ZmassEtaBins[iBin+1])
	  h1_ZmassData_Leading_muEtaM[iBin]->Fill(rMass_data);
      }
    }

    for(int iBin=0;iBin<nPtBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(reco1_data.charge>0)
      {
	if(iBin<3) if(reco1_data.pt >= ZmassPtBins[iBin] && reco1_data.pt < ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtP[iBin]->Fill(rMass_data);
	if(iBin==3) if(reco1_data.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtP[iBin]->Fill(rMass_data);
      }
      if(reco2_data.charge>0)
      {
	if(iBin<3) if(reco2_data.pt >= ZmassPtBins[iBin] && reco2_data.pt < ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtP[iBin]->Fill(rMass_data);
	if(iBin==3) if(reco2_data.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtP[iBin]->Fill(rMass_data);
      }
      //========================================
      //Check Charge Minus
      //========================================
      if(reco1_data.charge<0)
      {
	if(iBin<3) if(reco1_data.pt >= ZmassPtBins[iBin] && reco1_data.pt < ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtM[iBin]->Fill(rMass_data);
	if(iBin==3) if(reco1_data.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtM[iBin]->Fill(rMass_data);
      }
      if(reco2_data.charge<0)
      {
	if(iBin<3) if(reco2_data.pt >= ZmassPtBins[iBin] && reco2_data.pt < ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtM[iBin]->Fill(rMass_data);
	if(iBin==3) if(reco2_data.pt >= ZmassPtBins[iBin+1])
	  h1_ZmassData_muPtM[iBin]->Fill(rMass_data);
      }
    }
  }

  TFile *Hist_out = new TFile("histo_Zmass_Corr.root","RECREATE");

  h1_ZmassData    -> Write();
  h1_ZmassData_BB -> Write();
  h1_ZmassData_BE -> Write();
  h1_ZmassData_EE -> Write();

  h1_ZmassMC    -> Write();
  h1_ZmassMC_BB -> Write();
  h1_ZmassMC_BE -> Write();
  h1_ZmassMC_EE -> Write();

  h1_ZmassMCcorr    -> Write();
  h1_ZmassMCcorr_BB -> Write();
  h1_ZmassMCcorr_BE -> Write();
  h1_ZmassMCcorr_EE -> Write();

  for(int iBin=0;iBin<nEtaBins;iBin++){
    h1_ZmassData_muEtaP[iBin]   -> Write();
    h1_ZmassData_muEtaM[iBin]   -> Write();
    h1_ZmassMC_muEtaP[iBin]     -> Write();
    h1_ZmassMC_muEtaM[iBin]     -> Write();
    h1_ZmassMCcorr_muEtaP[iBin] -> Write();
    h1_ZmassMCcorr_muEtaM[iBin] -> Write();
    
    h1_ZmassData_Leading_muEtaP[iBin] -> Write();
    h1_ZmassData_Leading_muEtaM[iBin] -> Write();
    h1_ZmassMC_Leading_muEtaP[iBin] -> Write();
    h1_ZmassMC_Leading_muEtaM[iBin] -> Write();

    h1_ZmassData_Trailing_muEtaP[iBin] -> Write();
    h1_ZmassData_Trailing_muEtaM[iBin] -> Write();
    h1_ZmassMC_Trailing_muEtaP[iBin] -> Write();
    h1_ZmassMC_Trailing_muEtaM[iBin] -> Write();
  }

  for(int iBin=0;iBin<nPtBins;iBin++){
    h1_ZmassData_muPtP[iBin]   -> Write();
    h1_ZmassData_muPtM[iBin]   -> Write();
    h1_ZmassMC_muPtP[iBin]     -> Write();
    h1_ZmassMC_muPtM[iBin]     -> Write();
    h1_ZmassMCcorr_muPtP[iBin] -> Write();
    h1_ZmassMCcorr_muPtM[iBin] -> Write();
  }
  Hist_out->Close();  


  //  printhistos();
  //
  //   theFile->cd();
  //   theFile->Write();
  //   theFile->Close();
}


///////////////////////////////////////////
void bookhistos(){
  //	key1->GetXaxis()->SetTitle("P_{T}(#mu#mu) [GeV]");
  //	key1->GetYaxis()->SetTitle("# Events per GeV");
  //	test1->SetLineColor(2);
  //	test1->SetMarkerStyle(2);
  //	test1->SetMarkerColor(2);
  //	raw1->SetLineColor(3);

  for(int iPT = 0; iPT < NPThist; iPT++){
    for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      //hmuonRes[iK] = new TH1F(Form("hmuonRes%d",iK), Form("#mu resolution, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 200, -0.2, 0.2);
      //hmuonResSim[iK] = new TH1F(Form("hmuonResSim%d",iK), Form("#mu res., SF = %4.1f, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", ScaleFactor, PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 200, -0.2, 0.2);
      hmuonRes[iK] = new TH1F(Form("hmuonRes%d",iK), Form("MC Drell-Yan, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 80, -0.1, 0.1);
      hmuonResSim[iK] = new TH1F(Form("hmuonResSim%d",iK), Form("MC Drell-Yan, SF = %4.1f, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", ScaleFactor, PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 80, -0.1, 0.1);
      hDimuonRes[iK] = new TH1F(Form("hDimuonRes%d",iK), Form("#mu#mu resolution, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 100, -0.2, 0.2);
    }}
  for(int iPT = 0; iPT < NPThist_inv; iPT++){
    for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NPThist_inv;

      hDimuonMass[iK] = new TH1F(Form("hDimuonMass%d",iK), Form("MC Drell-Yan, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]), 60, 60., 120.);
      hDimuonMassSim[iK] = new TH1F(Form("hDimuonMassSim%d",iK), Form("MC Drell-Yan, SF = %4.1f, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", ScaleFactor, PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]), 60, 60., 120.);
      hDimuonMassDATA[iK] = new TH1F(Form("hDimuonMassDATA%d",iK), Form("DATA, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]), 60, 60., 120.);
    }}


  // unfolding matrix:
  key1_MC = new TH1F("key1_MC","",18,binning_histo);
  test1 = new TH1F("test1","",18,binning_histo);
  key1 = new TH1F("key1","",18,binning_histo);
  key1->GetXaxis()->SetTitle("P_{T}(#mu#mu) [GeV]");
  key1->GetYaxis()->SetTitle("# Events per GeV");
  raw1 = new TH1F("raw1","",18,binning_histo);
  test1->SetLineColor(2);
  test1->SetMarkerStyle(2);
  test1->SetMarkerColor(2);
  raw1->SetLineColor(3);

  for(int iBin=0;iBin<nEtaBins;iBin++){
    h1_ZmassData_muEtaP[iBin]   = new TH1F(Form("h1_ZmassData_muEtaP%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassData_muEtaM[iBin]   = new TH1F(Form("h1_ZmassData_muEtaM%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassMC_muEtaP[iBin]     = new TH1F(Form("h1_ZmassMC_muEtaP%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMC_muEtaM[iBin]     = new TH1F(Form("h1_ZmassMC_muEtaM%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMCcorr_muEtaP[iBin] = new TH1F(Form("h1_ZmassMCcorr_muEtaP%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMCcorr_muEtaM[iBin] = new TH1F(Form("h1_ZmassMCcorr_muEtaM%d",iBin), "Zmass MC", 40,80.,100);
    
    h1_ZmassData_Leading_muEtaP[iBin]   = new TH1F(Form("h1_ZmassData_Leading_muEtaP%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassData_Leading_muEtaM[iBin]   = new TH1F(Form("h1_ZmassData_Leading_muEtaM%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassMC_Leading_muEtaP[iBin]     = new TH1F(Form("h1_ZmassMC_Leading_muEtaP%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMC_Leading_muEtaM[iBin]     = new TH1F(Form("h1_ZmassMC_Leading_muEtaM%d",iBin), "Zmass MC", 40,80.,100);
    
    h1_ZmassData_Trailing_muEtaP[iBin]   = new TH1F(Form("h1_ZmassData_Trailing_muEtaP%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassData_Trailing_muEtaM[iBin]   = new TH1F(Form("h1_ZmassData_Trailing_muEtaM%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassMC_Trailing_muEtaP[iBin]     = new TH1F(Form("h1_ZmassMC_Trailing_muEtaP%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMC_Trailing_muEtaM[iBin]     = new TH1F(Form("h1_ZmassMC_Trailing_muEtaM%d",iBin), "Zmass MC", 40,80.,100);
  }

  for(int iBin=0;iBin<nPtBins;iBin++){
    h1_ZmassData_muPtP[iBin]   = new TH1F(Form("h1_ZmassData_muPtP%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassData_muPtM[iBin]   = new TH1F(Form("h1_ZmassData_muPtM%d",iBin), "Zmass Data", 40,80.,100);
    h1_ZmassMC_muPtP[iBin]     = new TH1F(Form("h1_ZmassMC_muPtP%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMC_muPtM[iBin]     = new TH1F(Form("h1_ZmassMC_muPtM%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMCcorr_muPtP[iBin] = new TH1F(Form("h1_ZmassMCcorr_muPtP%d",iBin), "Zmass MC", 40,80.,100);
    h1_ZmassMCcorr_muPtM[iBin] = new TH1F(Form("h1_ZmassMCcorr_muPtM%d",iBin), "Zmass MC", 40,80.,100);
  }
}
///////////////////////////////////////////
void printhistos(){
  /// print muon resolution data, MC and MC sim
  gStyle->SetOptFit(1111);
  TString gifname[NPThist];
  //TString gifname_inv[NPThist];
  TString gifmethod = "DoubleGauss";
  //TString gifmethodinv = "DoubleGauss";
  //TString gifmethod = "BWnonrel";
  //TString gifmethodinv = "Voigtian";
  gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
  gStyle->SetOptTitle(1); 

  TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 5);
  TF1* fitDoubleGauss2 = new TF1("fitDoubleGauss2", DoubleGauss, -0.1, 0.1, 5);
  for(int iPT = 0; iPT < NPThist; iPT++){
    gifname[iPT] = Form("zResolution/ResolutionPTScaleFactor%1.2f_"+gifmethod+"_%d", ScaleFactor, iPT);
    gifname[iPT] = gifname[iPT]+Extra;
    TCanvas *c20 = new TCanvas("c20","Resolution",3000,1700);
    TCanvas *c21 = new TCanvas("c21","Resolution mass",3000,1700);
    c20-> Divide(4,2);
    c21-> Divide(4,2);
    int Neta = NETAhist;
    if(Neta > 7) {Neta = 7; cout << "more then 8 bins, change canvas size" << endl;}
    for(int iETA = 0; iETA < Neta; iETA++){
      int iK = iPT + iETA*NPThist;
      c20 -> cd(iETA+1);
      // 
      //fitDoubleGauss->SetParameters(0., 0.01, 0.09, 0.24, 70.);
      fitDoubleGauss->SetParameters(0., 0.011, 0.041, 0.25, 2300.);
      if(PTbin[iPT] > 45 && ETAbin[iETA] > 0.9) fitDoubleGauss->SetParameters(0., 0.019, 0.16, 0.03, 320.); 
      fitDoubleGauss->FixParameter(0,0.);
      //fitDoubleGauss->SetParLimits(0, -0.01, 0.01);
      //fitDoubleGauss->SetParLimits(2, 0., 1.);
      fitDoubleGauss->SetParLimits(1, 0.005, 0.05);//restrict sigma1
      //fitDoubleGauss->SetParLimits(3, 0., 0.5);//restrict sigma2
      fitDoubleGauss->SetParName(0,"mean");
      fitDoubleGauss->SetParName(1,"sig1");
      fitDoubleGauss->SetParName(2,"Asig2");
      //fitDoubleGauss->SetParName(3,"mean2");
      fitDoubleGauss->SetParName(3,"sig2");
      fitDoubleGauss->SetParName(4,"Norm");

      //hmuonRes[iK] -> Fit(fitDoubleGauss,"RLE");
      hmuonRes[iK] -> Fit(fitDoubleGauss,"RLE");
      fitDoubleGauss2->SetParameters(fitDoubleGauss->GetParameter(0),fitDoubleGauss->GetParameter(1),fitDoubleGauss->GetParameter(2),fitDoubleGauss->GetParameter(3),fitDoubleGauss->GetParameter(4));
      if(PTbin[iPT] > 45 && ETAbin[iETA] > 0.9) fitDoubleGauss->SetParameters(0., 0.019, 0.16, 0.03, 320.);
      fitDoubleGauss2->FixParameter(0,0.);
      fitDoubleGauss2->SetParLimits(1, 0.005, 0.05);//restrict sigma1
      fitDoubleGauss2->SetParName(0,"mean");
      fitDoubleGauss2->SetParName(1,"sig1");
      fitDoubleGauss2->SetParName(2,"Asig2");
      fitDoubleGauss2->SetParName(3,"sig2");
      fitDoubleGauss2->SetParName(4,"Norm");

      hmuonRes[iK] -> Fit(fitDoubleGauss2,"RLE");

      TLegend* histinfo = SetLegend(.6,.57,1.,.73);
      histinfo->AddEntry(hmuonResSim[iK],Form("MC reco SF = %1.2f", ScaleFactor),"l");
      histinfo->AddEntry(hmuonRes[iK], "MC reco + Fit","lep");
      histinfo->AddEntry(fitDoubleGauss2, "Fit","l");


      hmuonRes[iK] ->  GetXaxis() ->SetNdivisions(505);// n = n1 + 100*n2 + 10000*n3
      hmuonRes[iK] -> GetXaxis()->SetTitle("p_{T}^{reco}(#mu)/p_{T}^{gen}(#mu)-1");
      hmuonRes[iK] -> GetYaxis()->SetTitle("Entries");

      hmuonRes[iK] -> Draw("e");
      hmuonResSim[iK] -> SetLineColor(kBlue);
      hmuonResSim[iK] -> Draw("samehist");
      histinfo -> Draw("same");
      //c21 -> cd(iETA+1);
      //hDimuonRes[iK] -> Draw("e");
    }
    c20 ->Print(gifname[iPT]+".png");
    c20 ->Print(gifname[iPT]+".root");
    //c21 ->Print(gifname_inv[iPT]); 
  }
  /// end: print muon resolution data, MC and MC sim

  /// print Mass reco in data, MC and MC sim     
  Double_t yScale[Nhist_inv];
  Double_t xScale[Nhist_inv];
  Double_t eyScale[Nhist_inv];
  Double_t exScale[Nhist_inv];

  TF1* fitBWnonrel = new TF1("fitBWnonrel", BWnonrel, -0.1, 0.1, 3);
  TF1* fitFuncVoigtian = new TF1("fitFuncVoigtian", FuncVoigtian, 80., 110., 4);
  //TF1* fitFuncVoigtianBG = new TF1("fitFuncVoigtianBG", FuncVoigtianBG, 70., 110., 7);
  TF1* fitFuncVoigtianBG = new TF1("fitFuncVoigtianBG", FuncVoigtianBG, 60., 120., 7);
  fitFuncVoigtianBG ->SetParName(0, "Norm");
  fitFuncVoigtianBG ->SetParName(1, "mass");
  fitFuncVoigtianBG ->SetParName(2, "width");
  fitFuncVoigtianBG ->SetParName(3, "resol");
  fitFuncVoigtianBG ->SetParName(4, "BG4");
  fitFuncVoigtianBG ->SetParName(5, "BG5");
  fitFuncVoigtianBG ->SetParName(6, "BG6");
  TF1* fitGauss = new TF1("fitGauss", Gauss, -0.1, 0.1, 3);
  //Double_t par[6];
  //TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 6);

  gStyle->SetOptFit(1111);
  TString gifname_inv[NPThist_inv];
  TString gifname_invDATA[NPThist_inv];
  TString gifmethodinv = "Voigtian";
  gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
  gStyle->SetOptTitle(1);
  for(int iPT = 0; iPT < NPThist_inv; iPT++){
    //gifname_inv[iPT] = "test.png";
    gifname_inv[iPT] = Form("zResolution/MassInv%deta%dptMCScaleFactor%1.2f_"+gifmethodinv+"_%d", ScaleFactor, NETAhist_inv, NPThist_inv, iPT);
    gifname_invDATA[iPT] = Form("zResolution/MassInv%deta%dptDATAScaleFactor%1.2f_"+gifmethodinv+"_%d", ScaleFactor, NETAhist_inv, NPThist_inv, iPT);
    gifname_inv[iPT] = gifname_inv[iPT]+Extra;
    gifname_invDATA[iPT] = gifname_invDATA[iPT]+Extra;
    //cout << "test" <<endl;
    //cout << gifname_inv[iPT] <<endl;
    TCanvas *c21 = new TCanvas("c21","Resolution mass",3000,1700);
    TCanvas *c21_data = new TCanvas("c21_data","Resolution mass",3000,1700);
    //c21-> Divide(4,2);
    //c21_data-> Divide(4,2);
    c21-> Divide(2,1);
    c21_data-> Divide(2,1);
    int Neta = NETAhist_inv;
    //if(Neta > 7) {Neta = 7; cout << "more then 8 bins, change canvas size" << endl;}
    if(Neta > 2) {Neta = 2; cout << "more then 2 bins, change canvas size" << endl;}
    for(int iETA = 0; iETA < Neta; iETA++){
      int iK = iPT + iETA*NPThist_inv;

      //fitFuncVoigtian ->SetParameters(200., 90., 2.495, 1.);
      //fitFuncVoigtian->SetParLimits(1, 80., 100.);
      //fitFuncVoigtian->FixParameter(2, 2.495);
      //fitFuncVoigtian->SetParLimits(3, 0.2, 5.);

      fitFuncVoigtianBG->SetParameters(200., 90., 2.495, 1., 2., -0.01, 0.);
      //fitFuncVoigtianBG->SetParameters(2000., 90., 2.495, 1., -20., -0.6, -0.0036);
      fitFuncVoigtianBG->SetParLimits(0, 0., 100000.);
      fitFuncVoigtianBG->SetParLimits(1, 80., 100.);
      fitFuncVoigtianBG->FixParameter(2, 2.495);
      fitFuncVoigtianBG->SetParLimits(3, 0.2, 5.);
      //fitFuncVoigtianBG->SetParLimits(4, -200., 200.);
      fitFuncVoigtianBG->SetParLimits(5, -20., 20.);
      fitFuncVoigtianBG->SetParLimits(6, -10., 10.);


      /////////////////
      c21 -> cd(iETA+1);
      hDimuonMassSim[iK] -> Fit(fitFuncVoigtianBG,"RLE");//Fit invmass sim with Scale Facotor
      //hDimuonMassSim[iK] -> Fit(fitFuncVoigtianBG,"RL");//Fit invmass sim with Scale Facotor
      Double_t par0 = fitFuncVoigtianBG ->GetParameter(0);
      Double_t par1 = fitFuncVoigtianBG ->GetParameter(1);
      Double_t par2 = fitFuncVoigtianBG ->GetParameter(2);
      Double_t par3 = fitFuncVoigtianBG ->GetParameter(3);
      Double_t par3err = fitFuncVoigtianBG ->GetParError(3);
      Double_t par4 = fitFuncVoigtianBG ->GetParameter(4);
      Double_t par5 = fitFuncVoigtianBG ->GetParameter(5);
      Double_t par6 = fitFuncVoigtianBG ->GetParameter(6);
      TF1* fitFuncVoigtianBGDraw = new TF1("fitFuncVoigtianBGDraw", FuncVoigtianBG, 60., 120., 7);
      fitFuncVoigtianBGDraw -> SetParameter(0, par0);
      fitFuncVoigtianBGDraw -> SetParameter(1, par1);
      fitFuncVoigtianBGDraw -> SetParameter(2, par2);
      fitFuncVoigtianBGDraw -> SetParameter(3, par3);
      fitFuncVoigtianBGDraw -> SetParameter(4, par4);
      fitFuncVoigtianBGDraw -> SetParameter(5, par5);
      fitFuncVoigtianBGDraw -> SetParameter(6, par6);

      //fitFuncVoigtianBG -> SetLineColor(2);//red
      hDimuonMassSim[iK] -> SetLineColor(kBlue);
      hDimuonMassSim[iK] -> GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
      hDimuonMassSim[iK] -> GetYaxis()->SetTitle("Entries");
      hDimuonMassSim[iK] -> Draw("hist");
      fitFuncVoigtianBGDraw -> Draw("same");
      hDimuonMass[iK] -> Draw("esame");
      //TF1 *fpol1 = new TF1 ("fpol1", "[3]*exp([0]-[1]*x+[2]*x*x)", 60., 120.);
      TF1 *fpol1 = new TF1 ("fpol1", "exp([0]-[1]*x+[2]*x*x)", 60., 120.);
      fpol1 -> SetParameter(0,par4);
      fpol1 -> SetParameter(1,par5);
      fpol1 -> SetParameter(2,par6);
      //fpol1 -> SetParameter(3,par0);
      fpol1 ->SetLineColor(4);
      fpol1 -> Draw("same");
      TLegend* histMCinfo = SetLegend(.56,.6,1.,.7);
      histMCinfo->AddEntry(hDimuonMassSim[iK],Form("MC reco SF = %1.2f + Fit", ScaleFactor),"l");
      histMCinfo->AddEntry(fitFuncVoigtianBGDraw, "Fit ","l");
      histMCinfo->AddEntry(hDimuonMass[iK], "MC reco","lep");
      histMCinfo-> Draw("same");

      /////////////////
      /////////////////
      //fitFuncVoigtianBG->FixParameter(4, par4); //fix from MC because not enough statistics
      fitFuncVoigtianBG->FixParameter(5, par5); //fix from MC because not enough statistics
      fitFuncVoigtianBG->FixParameter(6, par6); //fix from MC because not enough statistics
      c21_data -> cd(iETA+1);
      hDimuonMassDATA[iK] -> Fit(fitFuncVoigtianBG,"RLE");
      //hDimuonMassDATA[iK] -> Fit(fitFuncVoigtianBG,"RL");
      Double_t par0_data = fitFuncVoigtianBG ->GetParameter(0);
      Double_t par3_data = fitFuncVoigtianBG ->GetParameter(3);
      Double_t par3err_data = fitFuncVoigtianBG ->GetParError(3);
      Double_t par4_data = fitFuncVoigtianBG ->GetParameter(4);
      Double_t par5_data = fitFuncVoigtianBG ->GetParameter(5);
      Double_t par6_data = fitFuncVoigtianBG ->GetParameter(6);
      //Double_t err0 = fitFuncVoigtianBG ->GetParError(0);
      //TF1 *fpol1_data = new TF1 ("fpol1_data", "[0]+[1]*x", 70., 110.);
      //TF1 *fpol1_data = new TF1 ("fpol1_data", "[3]*exp([0]-[1]*x+[2]*x*x)", 60., 120.);
      TF1 *fpol1_data = new TF1 ("fpol1_data", "exp([0]-[1]*x+[2]*x*x)", 60., 120.);
      fpol1_data -> SetParameter(0,par4_data);
      fpol1_data -> SetParameter(1,par5_data);
      fpol1_data -> SetParameter(2,par6_data);
      //fpol1_data -> SetParameter(3,par0_data);
      fpol1_data ->SetLineColor(4);

      hDimuonMassDATA[iK] -> GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
      hDimuonMassDATA[iK] -> GetYaxis()->SetTitle("Entries");
      hDimuonMassDATA[iK] -> Draw("e");
      hDimuonMassSim[iK] -> SetLineColor(kBlue);
      hDimuonMassSim[iK] -> DrawNormalized("samehist", hDimuonMassDATA[iK]->Integral());
      fpol1_data -> Draw("same");
      TLegend* histDATAinfo = SetLegend(.6,.6, 1.,.7);
      histDATAinfo->AddEntry(hDimuonMassSim[iK],Form("MC reco SF = %1.2f", ScaleFactor),"l");
      histDATAinfo->AddEntry(hDimuonMassDATA[iK], "DATA reco + Fit","lep");
      histDATAinfo->AddEntry(fitFuncVoigtianBG, "Fit","l");
      histDATAinfo-> Draw("same");

      /////////////////
      //begin: calculate and fill scale factor
      xScale[iK] = iK+1;
      exScale[iK] = 0.5;
      yScale[iK] = 0;
      eyScale[iK] = 0;
      if(par3_data !=0 && par3 != 0){
	yScale[iK] = par3_data/par3;
	eyScale[iK] = yScale[iK]* pow(( par3err*par3err/par3/par3 + par3err_data*par3err_data/par3_data/par3_data  ),0.5);
      }
      //end: calculate and fill scale factor

    }
    c21 ->Print(gifname_inv[iPT]+".png");
    c21 ->Print(gifname_inv[iPT]+".root");
    c21_data ->Print(gifname_invDATA[iPT]+".png");
    c21_data ->Print(gifname_invDATA[iPT]+".root");
  }
  TF1 *fScale = new TF1 ("fScale", "[0]", 1., Nhist_inv );
  TGraphErrors* grScale = new TGraphErrors(Nhist_inv,xScale, yScale, exScale, eyScale);
  grScale->SetTitle(Form("Data/MC (SF = %1.2f)", ScaleFactor));
  grScale->SetMarkerColor(4);
  grScale->SetMarkerStyle(21);

  TString gifname_ScaleFaclor = Form("zResolution/ScaleFactor%deta%dptDATAScaleFactor%1.2f_"+gifmethodinv, ScaleFactor, NETAhist_inv, NPThist_inv);
  TCanvas *cScale = new TCanvas("cScale","Resolution mass",800,600);
  cScale -> cd();

  TH2F *h2 = new TH2F("h","Axes",Nhist_inv+1,0,Nhist_inv+1,100,0.,2.0);
  h2->GetXaxis()->SetNdivisions(Nhist_inv);
  h2->SetTitle(Form("Data/MC (SF = %1.2f)", ScaleFactor));
  for(int iPT = 0; iPT < NPThist_inv; iPT++){
    for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NPThist_inv+1;
      h2->GetXaxis()->SetBinLabel(iK, Form("p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f-%4.1f", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
    }}

  h2->Draw();
  //grScale->Draw("ALP");
  grScale->Draw("LP");
  grScale->Fit("fScale", "R");
  cScale ->Print(gifname_ScaleFaclor+Extra+".png");
  cScale ->Print(gifname_ScaleFaclor+Extra+".root");
  /// end: print Mass reco in data, MC and MC sim     


  ///////////////////
  // unfolding matrix
  //Make Matrix visualization plots
  rMatrixVisualization = new TH2F("rMatrixVisualization","Mapping True to Measured",18,binning_histo,18,binning_histo);
  rMatrixVisualization->GetYaxis()->SetTitle("p_{T}^{reco}(#mu#mu) [GeV]");
  rMatrixVisualization->GetXaxis()->SetTitle("p_{T}^{gen}(#mu#mu) [GeV]");
  rPriMatrixVisualization = new TH2F("rPriMatrixVisualization", "Direct Mapping Measured to True",18,binning_histo,18,binning_histo);
  rPriMatrixVisualization->GetYaxis()->SetTitle("p_{T}^{gen}(#mu#mu) [GeV]");
  rPriMatrixVisualization->GetXaxis()->SetTitle("p_{T}^{reco}(#mu#mu) [GeV]");
  rInvMatrixVisualization = new TH2F("rInvMatrixVisualization", "Inverted Matrix", 18, binning_histo, 18, binning_histo);
  rInvMatrixVisualization->GetYaxis()->SetTitle("p_{T}^{unfolded}(#mu#mu) [GeV]");
  rInvMatrixVisualization->GetXaxis()->SetTitle("p_{T}^{reco}(#mu#mu) [GeV]");
  for(int jj=0; jj<18; jj++)
    for(int kk=0; kk<18; kk++)
    {
      float temp = 0;
      if(realCounter[jj]!=0) temp = mapping[kk][jj]/realCounter[jj];
      mapping[kk][jj] = temp;
      rMatrixVisualization->SetBinContent(kk+1,jj+1,temp);
      temp = 0;
      if(realCounter[jj]!=0) temp = primeMapping[jj][kk]/realCounter[jj];
      primeMapping[jj][kk] = temp;
      //cout << "jj, kk = " << jj << ", " << kk << " value = " << temp 
      //<< " realCounter[jj] = " << realCounter[jj] << " primeMapping[jj][kk] = " << primeMapping[jj][kk]<< endl;
      rPriMatrixVisualization->SetBinContent(jj+1,kk+1,temp);
    }


  //Initialize Identy, Print matrix
  //TMatrixD inverse(kInverte,mapping);
  TMatrixD inverse(18,18);
  for(int c1=0;c1<18;c1++)
  {
    for(int c2=0; c2<18; c2++)
    {
      inverse(c1,c2)=mapping[c1][c2];
    }       
  }       
  inverse.Invert();

  float rowSum[18];
  std::cout << "**************************************" << endl;
  std::cout << "R Inverse with SF = " << ScaleFactor << endl;
  for(int pp1=0; pp1<18; pp1++)
  {
    rowSum[pp1]=0;
    for(int p3=0; p3<18; p3++)
    {
      std::cout << inverse(pp1,p3) << ", ";
      rInvMatrixVisualization->SetBinContent(pp1+1,p3+1,inverse(pp1,p3));
      rowSum[pp1]+=inverse(pp1,p3);
    }
    std::cout << std::endl;
  }
  std::cout << "**************************************" << endl;

  std::cout << "Row Sum: ";
  for(int psr=0; psr<18; psr++) std::cout << rowSum[psr] << ", " ;
  std::cout << std::endl;

  TString gifname_c1 = Form("zResolution/Matrix1ScaleFactor%1.2f_"+gifmethodinv+Extra, ScaleFactor);
  TString gifname_c2 = Form("zResolution/MatrixInvertScaleFactor%1.2f_"+gifmethodinv+Extra, ScaleFactor);
  TString gifname_c3 = Form("zResolution/Matrix2ScaleFactor%1.2f_"+gifmethodinv+Extra, ScaleFactor);

  gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
  gStyle->SetOptTitle(0);            // blue to red false color palette. Use 9 for b/w
  _c1 = new TCanvas("_c1","c1",100,0,600,600);
  _c2 = new TCanvas("_c2","c2",200,0,600,600);
  _c3 = new TCanvas("_c3","c3",300,0,600,600);
  _c1->cd();
  gPad -> SetLogx();
  gPad -> SetLogy();
  //gPad -> SetGrid();
  //rMatrixVisualization->SetTitle("Mapping Gen pT vs Reco pT"); 
  //rMatrixVisualization -> GetXaxis() -> SetRangeUser(1., binning[18]);
  //rMatrixVisualization -> GetYaxis() -> SetRangeUser(1., binning[18]);
  rMatrixVisualization->Draw("COLZ");
  // Title
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->SetTextColor(2);//red
  latex->SetNDC();
  latex->DrawLatex(0.3,0.94,Form("Matrix M: p_{T}^{reco} = M*p_{T}^{gen} (SF = %1.2f)", ScaleFactor));
  _c1 ->Print(gifname_c1+".png");
  _c1 ->Print(gifname_c1+".root");

  _c2->cd();
  gPad -> SetLogx();
  gPad -> SetLogy();
  rInvMatrixVisualization->Draw("COLZ");
  // Title
  latex->DrawLatex(0.15,0.94,Form("Matrix M_{invert}: p_{T}^{unfolded} = M_{invert}*p_{T}^{reco} (SF = %1.2f)", ScaleFactor));
  _c2 ->Print(gifname_c2+".png");
  _c2 ->Print(gifname_c2+".root");

  _c3->cd();
  gPad -> SetLogx();
  gPad -> SetLogy();
  rPriMatrixVisualization->Draw("COLZ");
  // Title
  latex->DrawLatex(0.3,0.94,Form("Matrix M: p_{T}^{reco} = M*p_{T}^{gen} (SF = %1.2f)", ScaleFactor));
  _c3 ->Print(gifname_c3+".png");
  _c3 ->Print(gifname_c3+".root");


  ////////////////////////////////////
  // Demonstrate ability of inversion to reproduce gen distribution.
  for(int uF1=0; uF1<18; uF1++){
    float binSum=0;
    for(int uF2=0; uF2<18; uF2++)
    {
      binSum+=inverse(uF1,uF2)*raw1->GetBinContent(uF2+1);
    }
    test1->SetBinContent(uF1+1,binSum);

  }
  _c4 = new TCanvas("_c4","c4",600,0,600,600);
  //_leg1 = new TLegend(.78,.13,.98,.33);
  _leg1 = new TLegend(.68,.77,.98,.93);
  _leg1->AddEntry(key1,"Gen ","l");
  _leg1->AddEntry(raw1,"Reco ","l");
  _leg1->AddEntry(test1,"Unfolded Reco","p");


  // make events per GeV
  float adjustedOcc = key1->GetBinContent(1)/(binning[1]-binning[0]);
  key1->SetBinContent(1,adjustedOcc);
  adjustedOcc = test1->GetBinContent(1)/(binning[1]-binning[0]);
  test1->SetBinContent(1, adjustedOcc);
  adjustedOcc = raw1->GetBinContent(1)/(binning[1]-binning[0]);
  raw1->SetBinContent(1, adjustedOcc);
  for(int binIt=2; binIt<19; binIt++)
  {
    float adjustedOcc = key1->GetBinContent(binIt)/key1->GetBinWidth(binIt);
    key1->SetBinContent(binIt,adjustedOcc);
    adjustedOcc = test1->GetBinContent(binIt)/test1->GetBinWidth(binIt);
    test1->SetBinContent(binIt, adjustedOcc);
    adjustedOcc = raw1->GetBinContent(binIt)/raw1->GetBinWidth(binIt);
    raw1->SetBinContent(binIt, adjustedOcc);
  }
  //end events per GeV

  _c4->cd();
  gPad -> SetLogx();
  //key1 -> GetXaxis() -> SetRangeUser(1., binning[18]);
  //key1 -> SetStats(0);
  //key1->GetXaxis()->TAxis::SetBinLabel(1,"0");
  key1->Draw();
  raw1->Draw("same");
  test1->Draw("same p");
  _leg1->Draw();

  _c5 = new TCanvas("_c5","c5",600,0,600,600);
  char* TitleX = "";

  DrawWithRatio(_c5, TitleX, test1, key1);

  ////////////////////////////////////
  std::cout << "Gen Occ by bin: \t";
  for(int g=0; g<16; g++) std::cout << key1->GetBinContent(g) << "\t";
  std::cout << "Integral: " << key1->Integral() << std::endl << "Reco Occ by bin:\t";
  for(int g=0; g<16; g++) std::cout << raw1->GetBinContent(g) << "\t";
  std::cout << "Integral: " << raw1->Integral() << std::endl << "Unfold Occ by bin:\t";
  for(int g=0; g<16; g++) std::cout << test1->GetBinContent(g) << "\t";
  std::cout << "Integral: " << test1->Integral() << std::endl;
  float rawChi=0;
  float testChi=0;
  for(int jg=1; jg<15; jg++)
  {
    rawChi +=  (raw1->GetBinContent(jg) - key1->GetBinContent(jg))*(raw1->GetBinContent(jg)  - key1->GetBinContent(jg));
    testChi+= (test1->GetBinContent(jg) - key1->GetBinContent(jg))*(test1->GetBinContent(jg) - key1->GetBinContent(jg));
    std::cout << "Step " << jg << ", raw: " << rawChi << ", test: " << testChi << std::endl;
  }
  float sqR = sqrt(rawChi);
  float sqT = sqrt(testChi);
  std::cout << "Chi2 Reco: " << sqR << ", unfolded: " << sqT << std::endl;

  // end: unfolding matrix
}
///////////////////////////////////////////

bool isLoose(_MuonInfo& muon) {
  bool isLoose=false;
  // kinematic cuts
  if (muon.pt < 20 )      return isLoose; // pt cut
  if (fabs(muon.eta) > 2.4) return isLoose; // eta cut
  if (muon.numValidTrackerHits < 10) return isLoose; // # hits in tracker
  // iso cut
  if (muon.trackIsoSumPt>10) return isLoose;
  // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return isLoose;
  if (fabs(muon.dz_PV) > 5) return isLoose;
  // if all the cuts are passed
  isLoose=true;
  return isLoose;
}


bool isKinTight(_MuonInfo& muon) {
  bool isKinTight=false;
  // minimum requirement to start investigating: to be loose
  if (!muon.isPFMuon) return isKinTight;
  if (!muon.isGlobal) return isKinTight;
  //if (muon.pt < 20 )      return isKinTight; // pt cut
  //if ((muon.trackIsoSumPt)/(muon.pt)>=0.1) return isKinTight;
  if (fabs(muon.eta) > 2.1) return isKinTight; // eta cut
  if ( muon.normChiSquare >= 10) return isKinTight;
  if (muon.numTrackerLayers < 6) return isKinTight; // # hits in tracker
  if ( muon.numValidPixelHits < 1 ) return isKinTight;
  if (fabs(muon.d0_PV) >= 0.2) return isKinTight;
  if (fabs(muon.dz_PV) > 5) return isKinTight;
  if ( muon.numOfMatchedStations < 2  ) return isKinTight;
  if ( muon.numValidMuonHits < 1 ) return isKinTight;
  isKinTight=true;
  return isKinTight;
}
Double_t BWnonrel(Double_t* x, Double_t* par){


  //Double_t eff = par[0]*(TMath::Erf(par[1]*x[0]-par[2]) - 1.) + par[3];
  //Non rel. Breit-Wigner                               mass     width 
  Double_t bw_nonrel = par[0]*(TMath::BreitWigner(x[0], par[1], par[2])); 
  return bw_nonrel;
}

Double_t Gauss(Double_t* x, Double_t* par){


  //Double_t eff = par[0]*(TMath::Erf(par[1]*x[0]-par[2]) - 1.) + par[3];
  //Non rel. Breit-Wigner                               mass     width 
  //Double_t bw_nonrel = par[0]*(TMath::BreitWigner(x[0], par[1], par[2])); 
  Double_t gauss = par[0]*(TMath::Gaus(x[0], par[1], par[2])); 
  return gauss;
}
Double_t FuncVoigtian(Double_t* x, Double_t* par){

  //--------------------------------------------------------------------//
  //Fit parameters:
  //par[0]=Width (scale) Breit-Wigner
  //par[1]=Most Probable (MP, location) Breit mean
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t twoPi = 6.2831853071795;//2Pi

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  //Double_t np = 60.0;      // number of convolution steps
  //Double_t sc =   3.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Breit and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    //fland = TMath::BreitWigner(xx,par[1],par[0]);
    fland = TMath::BreitWigner(xx,par[1],par[2]);
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    //fland = TMath::BreitWigner(xx,par[1],par[0]);
    fland = TMath::BreitWigner(xx,par[1],par[2]);
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[0] * step * sum * invsq2pi / par[3]);
}

Double_t FuncVoigtianBG(Double_t* x, Double_t* par){

  //--------------------------------------------------------------------//
  //Fit parameters:
  //par[0]=Width (scale) Breit-Wigner
  //par[1]=Most Probable (MP, location) Breit mean
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t twoPi = 6.2831853071795;//2Pi

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  //Double_t np = 60.0;      // number of convolution steps
  //Double_t sc =   3.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Breit and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    //fland = TMath::BreitWigner(xx,par[1],par[0]);
    fland = TMath::BreitWigner(xx,par[1],par[2]);
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    //fland = TMath::BreitWigner(xx,par[1],par[0]);
    fland = TMath::BreitWigner(xx,par[1],par[2]);
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  //return (par[0] * step * sum * invsq2pi / par[3] + par[4] + par[5]*x[0]);
  return (par[0] * step * sum * invsq2pi / par[3] + exp(par[4]-par[5]*x[0]+par[6]*x[0]*x[0]));
  //return (par[0] * step * sum * invsq2pi / par[3] + par[0] * exp(par[4]-par[5]*x[0]+par[6]*x[0]*x[0]));
}

// Sum of background and peak function
Double_t DoubleGauss(Double_t *x, Double_t *par) 
{ 

  //return BWnonrel(x,par) + Gauss(x,&par[3]);
  Double_t dgauss = 0.;
  //dgauss = Gauss(x,par) + Gauss(x,&par[3]);

  // sig1 < sig2 always
  //if(par[1] < par[3])dgauss = par[4]* ( TMath::Gaus(x[0],par[0],par[1],"kTRUE") + par[2] * TMath::Gaus(x[0],par[0],par[3],"kTRUE") );
  //if(par[1] < par[3]){//normalized gauss par[2] could correlated with par[3] here 
  //          dgauss =  1/par[1]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[1]/par[1]);
  //          dgauss = dgauss + par[2]/par[3]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[3]/par[3]);
  //          dgauss = par[4]*dgauss;
  //}
  if(par[1] < par[3]){//not normalized gauss both gauss are  = 1 at x[0]=par[0]  
    dgauss =  exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[1]/par[1]);
    dgauss = dgauss + par[2]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[3]/par[3]);
    dgauss = par[4]*dgauss;
  }
  //if (par[2]> 0.) dgauss = dgauss/(1+par[2]);  

  return dgauss;
}

//------------------------------------------------------------------------------
// Draw projections and residuals
//------------------------------------------------------------------------------
void DrawWithRatio(TCanvas *canvas, char *cTitle,
    TH1F *hNum, TH1F *hDen)
{

  // sanity check
  if (hNum->GetNbinsX() != hDen->GetNbinsX()){
    std::cout<< " *** Error: binning not consistent between data"
      << " and MC -> Exit!\n";
    return;
  }


  hNum->Sumw2();
  hDen->Sumw2();


  TH1F *hPull = (TH1F*)hNum ->Clone("hPull");
  hPull->Sumw2();
  hPull->Divide(hDen);

  //----------------------------------------------------------------------------
  // Create the pads
  //----------------------------------------------------------------------------
  TPad* pad1;
  TPad* pad2;

  pad1 = new TPad("pad1","This is pad1",0.02,0.30,0.98,0.98,0);
  pad2 = new TPad("pad2","This is pad2",0.02,0.01,0.98,0.29,0);

  pad1->SetLogx();
  pad2->SetLogx();
  pad1->SetBottomMargin(0.01);
  pad2->SetBottomMargin(0.33);
  pad2->SetTopMargin   (0.10);

  pad1->Draw(); // Projections pad
  pad2->Draw(); // Residuals   pad

  _leg2 = new TLegend(.68,.77,.98,.93);
  _leg2->AddEntry(hDen,"Gen ","l");
  _leg2->AddEntry(hNum,"Unfolded Reco","p");
  pad1->cd();
  hDen->Draw("histo");
  hNum->Draw("pe same");
  _leg2->Draw(); 
  PrintItLog(pad1,cTitle);

  //   TLegend* leg = SetLegend(0.73, 0.7, 0.92, 0.89);
  //   leg -> AddEntry(hDen," no mass cut","f");
  //   leg -> AddEntry(hNum," 60 < M^{#mu #mu} < 120","p");
  //   leg ->Draw("same");
  //----------------------------------------------------------------------------
  // Residuals pad
  //----------------------------------------------------------------------------
  pad2->cd();

  TAxis *xPull = NULL;
  TAxis *yPull = NULL;
  char xAxisName[200];
  sprintf(xAxisName,"%s",hDen->GetXaxis()->GetTitle());
  axis1F(hPull,xPull,yPull,xAxisName,"ratio");

  if (hPull->GetMaximum() > 100) {
    hPull->SetMinimum(-100);
    hPull->SetMaximum( 100);
  }

  hPull->GetXaxis()->SetLabelOffset(0.005);
  hPull->GetXaxis()->SetLabelSize  (0.11);
  hPull->GetXaxis()->CenterTitle(1);
  hPull->GetXaxis()->SetTitleOffset(1.10);
  hPull->GetXaxis()->SetTitleSize  (0.12);
  hPull->GetXaxis()->SetNdivisions(7);

  hPull->GetYaxis()->SetLabelSize  (0.09);
  hPull->GetYaxis()->CenterTitle(1);
  hPull->GetYaxis()->SetTitleOffset(0.5);
  hPull->GetYaxis()->SetTitleSize  (0.12);

  hPull->SetMaximum(1.5);
  hPull->SetMinimum(0.5);
  hPull->Draw("pe");

  pad2->Update();
  pad2->GetFrame()->DrawClone();

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

