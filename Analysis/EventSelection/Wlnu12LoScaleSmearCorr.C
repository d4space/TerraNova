// $Log: Wlnu12LoScaleSmearCorr.C,v $
// Revision 1.8  2013/09/13 00:09:32  salee
// *** empty log message ***
//
// Revision 1.7  2013/09/12 05:10:29  sangilpark
// *** empty log message ***
//
// Revision 1.6  2013/08/31 06:56:40  khakim
// *** empty log message ***
//
#define Wlnu12LoScaleSmearCorr_cxx
//#include <iostream>
//#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
//#include <TSystem.h>                      // interface to OS
#include <TBenchmark.h>                   // class to track macro running statistics
#include <algorithm>
#include "Wlnu12LoScaleSmearCorr.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TLorentzVector.h>
//#include "../Utils/MyTools.hh"	          // various helper functions
#include <TRandom3.h>
#include <TRandom.h>

#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMLorentzVectorD;
typedef PtEtaPhiMLorentzVectorD PtEtaPhiMLorentzVector;

void Wlnu12LoScaleSmearCorr::Loop()
{
  //gSystem->Load("libMathCore");
  //gSystem->Load("libPhysics");
  //using namespace ROOT::Math;

  Debug=false;
  cout<<"==================================================================="<<endl;
  cout<<"Wlnu12LoScaleSmearCorr Analysis with Mode: "<<Mode<<"  AnaChannel: "<<AnaChannel<<endl;
  cout<<"==================================================================="<<endl;
  gBenchmark->Start("Wlnu12LoScaleSmearCorr");

  if (fChain == 0) return;
   //int Ntries = fChain->GetEntriesFast(); this gives 1234567890 kkk
  Ntries = fChain->GetEntries();

  cout<<"Total: "<<Ntries<<endl;

  double leptMass;
  
  leptMass = 0.105658389;
  if(AnaChannel == "Electron2012LoPU")
    leptMass = 0.000510998902;
	
  //============================================
  // Looping for each Event 
  //============================================
  for (int i(0); i<Ntries;i++)
  {
    evtCnt = i;
    //===============================
    //W study
    //===============================
    if(i % 100000 == 0) cout<<i<<"th event"<<endl;
    if(Debug)cout<<"check point 1"<<endl;

    fChain->GetEntry(i);
    //===========================
    //Initialization of Variables
    //===========================
    InitVar4Evt();

    // Dump MET informations To put MET as TLorentz vector
    //DumpMETs();

    //===================
    // Check the channel : To check if the ntuple is for each lepton flavor
    //===================
    if(Wlnu12LoBase::CheckChannel()!=0) exit(-1);

    //============
    //Trigger Cut
    //============
    if(Wlnu12LoBase::TriggerCut() !=0) continue;

    //Vertex Study===========================
    if(VertexCut() !=0) continue;

    //===================
    // Calculate Event Weight
    //=====================
    mTTW = CalcEvtWeight();
    
    //cout<<"Muon size: "<<wMuons.pt->size()<<endl;
    //cout<<"W    size: "<<W_pt->size()<<endl;
  
    ZbestSelect();
    
    if(!Z.Pass) continue;
    
    mNZevt++;
    if(Mode == "ScaleMakeMC" || Mode == "ScaleMakeRD")
    {
      int etaRange1 = EtaRange(Z.Lep1etaSC);
      int etaRange2 = EtaRange(Z.Lep2etaSC);
      if(AnaChannel == "Electron2012LoPU") Fill_EleZmassDaughEta(etaRange1,etaRange2);
      if(AnaChannel == "Muon2012LoPU")     Fill_MuZmassDaughEta(etaRange1,etaRange2);
    }
    
    Fill_ZHisto();
    if(Mode == "ScaleMakeMC")
    {
      //===================
      //Wpt study apply shift correction to the Data
      //=====================
      //Scale_corrZlep1Pt = 1.0/GetScaleCorr(Z.Lep1etaSC)*Z.Lep1Pt;
      //Scale_corrZlep2Pt = 1.0/GetScaleCorr(Z.Lep2etaSC)*Z.Lep2Pt;

      //===================
      //Zpt study did not apply shift correction to the Data. We don't apply shift correction to the Data
      //=====================
      Scale_corrZlep1Pt = Z.Lep1Pt;
      Scale_corrZlep2Pt = Z.Lep2Pt;
      
      PtEtaPhiMLorentzVector Zlept1_4(Scale_corrZlep1Pt,Z.Lep1etaSC,Z.Lep1Phi,leptMass);
      PtEtaPhiMLorentzVector Zlept2_4(Scale_corrZlep2Pt,Z.Lep2etaSC,Z.Lep2Phi,leptMass);
      
      smearSFLep1 = gRandom->Gaus(Zlept1_4.E(), GetSmearCorr(Z.Lep1etaSC))/Zlept1_4.E();
      smearSFLep2 = gRandom->Gaus(Zlept2_4.E(), GetSmearCorr(Z.Lep2etaSC))/Zlept2_4.E();
      PtEtaPhiMLorentzVector Z_4 = smearSFLep1*Zlept1_4 + smearSFLep2*Zlept2_4;
      Z.mass=Z_4.M();
    }
    Fill_CorrectedZHisto();
  }//Ntries

  cout<<"Passed W evts: "<<mNWevt<<endl;
  Fout<<"Passed W evts: "<<mNWevt<<endl;
  //Results======================
  cout<<"selected converted: "<<evtSelected<<" +- "<<TMath::Sqrt(evtSelected)<<endl;
  Fout<<"selected converted: "<<evtSelected<<" +- "<<TMath::Sqrt(evtSelected)<<endl;
  Fout<<"selected events for each bin"<<endl;
  for( int i(0); i<NwPtBin; i++)
  {
    Fout<<i<<"   "<<mNselect4WptBin[i]<<endl;
  }

  // Notice: Use one of Write_Histo or myFile->Write
  // Write_Histo: to Save specific histograms
  // myFile->Write: to Save all Histograms
  Write_ZHisto();
  myFile->Write();
  myFile->Close();
  Fout.close();
  gBenchmark->Show("Wlnu12LoScaleSmearCorr");
}

void Wlnu12LoScaleSmearCorr::Nselected4Bin()
{
 // for(int i(0);i<NwPtBin;i++)
 // {
 //   if( W.pt >= WptBins[i] && W.pt <WptBins[i+1]) mNselect4WptBin[i]+=mTTW;
 // }
}
int Wlnu12LoScaleSmearCorr::InitVar()
{
//  cout<<"Initialize variable at Wlnu12LoScaleSmearCorr class ==========="<<endl;
//  evtCnt = 0;
//  mNWevt = 0;
//  mNZevt = 0;
//  TString FoutName = mResultDir+"/"+OutNameBase+".txt";
//  Fout.open(FoutName);
//  for(int i(0);i<NwPtBin;i++)
//  {
//    mNselect4WptBin[i]=0;
//  }
//  // Recoil CorrWptection initializaWpttion
//  // Recoil CorrWptection Parameter WptFiles
//  if( (  Mode == "AllCorrectionsMC"
//      || Mode == "RecoilCorrMC")
//      || Mode =="DumpUnfInfo" )
//  {
//    if(AnaChannel == "Muon2012LoPU" )
//    {
//      Rcl.ZRDfilename="../Recoil/ZmmData/fits.root";
//      Rcl.ZMCfilename="../Recoil/ZmmMC/fits.root";
//      Rcl.Wpfilename="../Recoil/WmpMC/fits.root";
//      Rcl.Wmfilename="../Recoil/WmmMC/fits.root";
//    }else if((AnaChannel == "Electron2012LoPU" ) || AnaChannel == "ElectronHighPU" )
//    {
//      Rcl.ZRDfilename="../Recoil/ZeeData/fits.root";
//      Rcl.ZMCfilename="../Recoil/ZeeMC/fits.root";
//      Rcl.Wpfilename="../Recoil/WepMC/fits.root";
//      Rcl.Wmfilename="../Recoil/WemMC/fits.root";
//    }
//    // RecoilCorrection Object.
//    RecoilCorr= new RecoilCorrector(
//      Rcl.ZRDfilename,
//      Rcl.Wpfilename,Rcl.Wmfilename,
//      Rcl.ZMCfilename,
//      0x1234);
//  //Int_t iSeed=0xDEADBEEF default seed for random number generator at constructor
//  }
//  return 0;
}
int Wlnu12LoScaleSmearCorr::InitVar4Evt()
{
  //cout<<"Wlnu12LoScaleSmearCorr::InitVar4Evt ==========================="<<endl;
  Wlnu12LoBase::InitVar4Evt();
  return 0;
}
int Wlnu12LoScaleSmearCorr::InitHistogram()
{
  myFile = new TFile(mResultDir+"/"+OutNameBase+".root","RECREATE");
  h1_Zmass   = new TH1D("h1_Zmass","Z Mass",20,80.,100);
  h1_Zmass_BB= new TH1D("h1_Zmass_BB","Inv Mass for dilepts BB",20,80.,100.);
  h1_Zmass_BE= new TH1D("h1_Zmass_BE","Inv Mass for dilepts BE",20,80.,100.);
  h1_Zmass_EE= new TH1D("h1_Zmass_EE","Inv Mass for dilepts EE",20,80.,100.);
  
  h1_ZmassCorr= new TH1D("h1_ZmassCorr","Inv Mass for dilepts after Scale&Smear",20,80.,100.);
  h1_ZmassCorr_BB= new TH1D("h1_ZmassCorr_BB","Inv Mass for dilepts after Scale&Smear BB",20,80.,100.);
  h1_ZmassCorr_BE= new TH1D("h1_ZmassCorr_BE","Inv Mass for dilepts after Scale&Smear BE",20,80.,100.);
  h1_ZmassCorr_EE= new TH1D("h1_ZmassCorr_EE","Inv Mass for dilepts after Scale&Smear EE",20,80.,100.);

  h1_ZLep1Pt = new TH1D("h1_ZLep1Pt","First lepton pt", 150,0,150);
  h1_ZLep2Pt = new TH1D("h1_ZLep2Pt","Second lepton pt",150,0,150);
  h1_ZLepPt_p = new TH1D("h1_ZLepPt_p","Lept^{+} pt", 150,0,150);
  h1_ZLepPt_m = new TH1D("h1_ZLepPt_m","Lept^{-} pt", 150,0,150);

  for(int iBin(0);iBin<ScaleBins;iBin++){
    sprintf(histName,"h1_Zmass_muEtaP_%d",iBin);
    h1_Zmass_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_Zmass_muEtaM_%d",iBin);
    h1_Zmass_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    
    sprintf(histName,"h1_ZmassCorr_muEtaP_%d",iBin);
    h1_ZmassCorr_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_ZmassCorr_muEtaM_%d",iBin);
    h1_ZmassCorr_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);

    sprintf(histName,"h1_Zmass_noOverLap_muEtaP_%d",iBin);
    h1_Zmass_noOverLap_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_Zmass_noOverLap_muEtaM_%d",iBin);
    h1_Zmass_noOverLap_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    
    sprintf(histName,"h1_ZmassCorr_noOverLap_muEtaP_%d",iBin);
    h1_ZmassCorr_noOverLap_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_ZmassCorr_noOverLap_muEtaM_%d",iBin);
    h1_ZmassCorr_noOverLap_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);

    sprintf(histName,"h1_Zmass_LeadingLept_noOverLap_muEtaP_%d",iBin);
    h1_Zmass_LeadingLept_noOverLap_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_Zmass_LeadingLept_noOverLap_muEtaM_%d",iBin);
    h1_Zmass_LeadingLept_noOverLap_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    
    sprintf(histName,"h1_ZmassCorr_LeadingLept_noOverLap_muEtaP_%d",iBin);
    h1_ZmassCorr_LeadingLept_noOverLap_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_ZmassCorr_LeadingLept_noOverLap_muEtaM_%d",iBin);
    h1_ZmassCorr_LeadingLept_noOverLap_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);

    sprintf(histName,"h1_Zmass_LeadingLept_muEtaP_%d",iBin);
    h1_Zmass_LeadingLept_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_Zmass_LeadingLept_muEtaM_%d",iBin);
    h1_Zmass_LeadingLept_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    
    sprintf(histName,"h1_ZmassCorr_LeadingLept_muEtaP_%d",iBin);
    h1_ZmassCorr_LeadingLept_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_ZmassCorr_LeadingLept_muEtaM_%d",iBin);
    h1_ZmassCorr_LeadingLept_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);

    sprintf(histName,"h1_Zmass_TrailingLept_noOverLap_muEtaP_%d",iBin);
    h1_Zmass_TrailingLept_noOverLap_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_Zmass_TrailingLept_noOverLap_muEtaM_%d",iBin);
    h1_Zmass_TrailingLept_noOverLap_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    
    sprintf(histName,"h1_ZmassCorr_TrailingLept_noOverLap_muEtaP_%d",iBin);
    h1_ZmassCorr_TrailingLept_noOverLap_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_ZmassCorr_TrailingLept_noOverLap_muEtaM_%d",iBin);
    h1_ZmassCorr_TrailingLept_noOverLap_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);

    sprintf(histName,"h1_Zmass_TrailingLept_muEtaP_%d",iBin);
    h1_Zmass_TrailingLept_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_Zmass_TrailingLept_muEtaM_%d",iBin);
    h1_Zmass_TrailingLept_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    
    sprintf(histName,"h1_ZmassCorr_TrailingLept_muEtaP_%d",iBin);
    h1_ZmassCorr_TrailingLept_muEtaP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_ZmassCorr_TrailingLept_muEtaM_%d",iBin);
    h1_ZmassCorr_TrailingLept_muEtaM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
  }

  for(int iBin(0);iBin<4;iBin++){
    sprintf(histName,"h1_Zmass_muPtP_%d",iBin);
    h1_Zmass_muPtP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_Zmass_muPtM_%d",iBin);
    h1_Zmass_muPtM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);

    sprintf(histName,"h1_ZmassCorr_muPtP_%d",iBin);
    h1_ZmassCorr_muPtP[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
    sprintf(histName,"h1_ZmassCorr_muPtM_%d",iBin);
    h1_ZmassCorr_muPtM[iBin]= new TH1D(histName,"Zmass",40,80.,100.);
  }

  if(Mode == "ScaleMakeMC" || Mode == "ScaleMakeRD")
  {
    if( AnaChannel=="Electron2012LoPU")
    {
      for(int i(0);i<ScElCombiBins;i++)
      {
	sprintf(histName,"h1_ZmassDaughEtaEle_%d",i);
	h1_ZmassDaughEtaEle[i]= new TH1D(histName,"ZmassDaughterEtaEle",20,80.,100.);
      }
      for(int i(0);i<ScElCombiBinsDiag;i++)
      {
	sprintf(histName,"h1_ZmassDaughEtaEleDiag_%d",i);
	h1_ZmassDaughEtaEleDiag[i]= new TH1D(histName,"ZmassDaughterEtaEleDiag",20,80.,100.);
      }
    }
    if(AnaChannel=="Muon2012LoPU")
    {
      for(int i(0);i<ScMuCombiBins;i++)
      {
	sprintf(histName,"h1_ZmassDaughEtaMu_%d",i);
	h1_ZmassDaughEtaMu[i]= new TH1D(histName,"ZmassDaughterEtaMu",20,80.,100.);
      }
    }
  }
  return 0;
}
int Wlnu12LoScaleSmearCorr::Fill_ZHisto()
{
  h1_Zmass->Fill(Z.mass,mTTW);
  h1_ZLep1Pt -> Fill(Z.Lep1Pt,mTTW);
  h1_ZLep2Pt -> Fill(Z.Lep2Pt,mTTW);
  if(Z.Lep1Charge>0)
    h1_ZLepPt_p -> Fill(Z.Lep1Pt,mTTW);
  if(Z.Lep2Charge>0)
    h1_ZLepPt_p -> Fill(Z.Lep2Pt,mTTW);
  if(Z.Lep1Charge<0)
    h1_ZLepPt_m -> Fill(Z.Lep1Pt,mTTW);
  if(Z.Lep2Charge<0)
    h1_ZLepPt_m -> Fill(Z.Lep2Pt,mTTW);
  
  if( AnaChannel=="Electron2012LoPU")
  {
    //Barrel Barrel
    if((fabs(Z.Lep1etaSC) >= 0.0 && fabs(Z.Lep1etaSC) < 1.4442) && (fabs(Z.Lep2etaSC) >= 0.0 && fabs(Z.Lep2etaSC) < 1.4442))
      h1_Zmass_BB->Fill(Z.mass,mTTW);
    //Barrel Endcap
    if((fabs(Z.Lep1etaSC) >= 0.0 && fabs(Z.Lep1etaSC) < 1.4442) && (fabs(Z.Lep2etaSC) >= 1.566 && fabs(Z.Lep2etaSC) < 2.5))
      h1_Zmass_BE->Fill(Z.mass,mTTW);
    if((fabs(Z.Lep1etaSC) >= 1.566 && fabs(Z.Lep1etaSC) < 2.5) && (fabs(Z.Lep2etaSC) >= 0.0 && fabs(Z.Lep2etaSC) < 1.4442))
      h1_Zmass_BE->Fill(Z.mass,mTTW);
    // Endcap Endcap
    if((fabs(Z.Lep1etaSC) >= 1.566 && fabs(Z.Lep1etaSC) < 2.5) && (fabs(Z.Lep2etaSC) >= 1.566 && fabs(Z.Lep2etaSC) < 2.5))
      h1_Zmass_EE->Fill(Z.mass,mTTW);
  }
  
  if( AnaChannel=="Muon2012LoPU")
  {
    //Barrel Barrel
    if((fabs(Z.Lep1etaSC) >= 0. && fabs(Z.Lep1etaSC) < 1.) && (fabs(Z.Lep2etaSC) >= 0. && fabs(Z.Lep2etaSC) < 1.))
      h1_Zmass_BB->Fill(Z.mass,mTTW);
    //Barrel Endcap
    if((fabs(Z.Lep1etaSC) >= 0. && fabs(Z.Lep1etaSC) < 1.) && (fabs(Z.Lep2etaSC) >= 1. && fabs(Z.Lep2etaSC) < 2.4))
      h1_Zmass_BE->Fill(Z.mass,mTTW);
    if((fabs(Z.Lep1etaSC) >= 1. && fabs(Z.Lep1etaSC) < 2.4) && (fabs(Z.Lep2etaSC) >= 0. && fabs(Z.Lep2etaSC) < 1.))
      h1_Zmass_BE->Fill(Z.mass,mTTW);
    // Endcap Endcap
    if((fabs(Z.Lep1etaSC) >= 1. && fabs(Z.Lep1etaSC) < 2.4) && (fabs(Z.Lep2etaSC) >= 1. && fabs(Z.Lep2etaSC) < 2.4))
      h1_Zmass_EE->Fill(Z.mass,mTTW);
    
    //========================================
    //Fill Zmass histograms at eta ranges
    //========================================
    for(int iBin(0);iBin<ScaleBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(Z.Lep1Charge>0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	{
          //========================================
          //Data: if low statistics
	  //========================================
	  //if(Mode == "ScaleMakeRD")
	  //{
	  //  if(Z.Lep1Pt > Z.Lep2Pt)
	  //  {
	  //    if(Z.Lep1etaSC < 0.35 || (Z.Lep1etaSC > 0.7 && Z.Lep1etaSC < 1.40))
	  //      h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //    if(Z.Lep1etaSC > 1.75)
	  //      h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //  }

	  //  if(Z.Lep1Pt <= Z.Lep2Pt)
	  //    if(Z.Lep1etaSC < 1.75)
	  //      h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //}

	  //if(Mode == "ScaleMakeMC")
	  //  h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  
	  h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}

        //========================================
        //Fill Charge Plus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Charge>0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	{
          //========================================
          //Data: if low statistics
	  //========================================
	  //if(Mode == "ScaleMakeRD")
	  //{
	  //  if(Z.Lep2Pt > Z.Lep1Pt)
	  //  {
	  //    if(Z.Lep2etaSC < 0.35 || (Z.Lep2etaSC > 0.7 && Z.Lep2etaSC < 1.40))
	  //      h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //    if(Z.Lep2etaSC > 1.75)
	  //      h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //  }

	  //  if(Z.Lep2Pt <= Z.Lep1Pt)
	  //    if(Z.Lep2etaSC < 1.75)
	  //      h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //}

	  //if(Mode == "ScaleMakeMC")
	  //  h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  
	  h1_Zmass_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
	
        //========================================
        //Fill Charge Plus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Charge<0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Charge<0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Fill Leading lepton histo
      //========================================
      if(Z.Lep1Pt>Z.Lep2Pt && Z.Lep1Charge>0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_LeadingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_LeadingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt>Z.Lep1Pt && Z.Lep2Charge>0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_LeadingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_LeadingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Pt>Z.Lep2Pt && Z.Lep1Charge<0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_LeadingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_LeadingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt>Z.Lep1Pt && Z.Lep2Charge<0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_LeadingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_LeadingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Fill Trailing lepton histo
      //========================================
      if(Z.Lep1Pt<Z.Lep2Pt && Z.Lep1Charge>0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_TrailingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_TrailingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt<Z.Lep1Pt && Z.Lep2Charge>0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_TrailingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_TrailingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Pt<Z.Lep2Pt && Z.Lep1Charge<0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_TrailingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_TrailingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt<Z.Lep1Pt && Z.Lep2Charge<0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_Zmass_TrailingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_Zmass_TrailingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
    }
    
    //========================================
    //Fill Zmass histograms at pT ranges
    //========================================
    for(int iBin(0);iBin<4;iBin++){
      //Check Charge Plus
      //========================================
      if(Z.Lep1Charge>0)
      {
	if(iBin<3) if(Z.Lep1Pt >= ZmassPtBins[iBin] && Z.Lep1Pt < ZmassPtBins[iBin+1])
	  h1_Zmass_muPtP[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep1Pt >= ZmassPtBins[iBin+1])
	  h1_Zmass_muPtP[iBin]->Fill(Z.mass,mTTW);
      }
      
      if(Z.Lep2Charge>0)
      {
	if(iBin<3) if(Z.Lep2Pt >= ZmassPtBins[iBin] && Z.Lep2Pt < ZmassPtBins[iBin+1])
	  h1_Zmass_muPtP[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep2Pt >= ZmassPtBins[iBin+1])
	  h1_Zmass_muPtP[iBin]->Fill(Z.mass,mTTW);
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Charge<0)
      {
	if(iBin<3) if(Z.Lep1Pt >= ZmassPtBins[iBin] && Z.Lep1Pt < ZmassPtBins[iBin+1])
	  h1_Zmass_muPtM[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep1Pt >= ZmassPtBins[iBin+1])
	  h1_Zmass_muPtM[iBin]->Fill(Z.mass,mTTW);
      }
      
      if(Z.Lep2Charge<0)
      {
	if(iBin<3) if(Z.Lep2Pt >= ZmassPtBins[iBin] && Z.Lep2Pt < ZmassPtBins[iBin+1])
	  h1_Zmass_muPtM[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep2Pt >= ZmassPtBins[iBin+1])
	  h1_Zmass_muPtM[iBin]->Fill(Z.mass,mTTW);
      }
    }
  }
  return 0;
}

int Wlnu12LoScaleSmearCorr::Fill_CorrectedZHisto()
{
  h1_ZmassCorr->Fill(Z.mass,mTTW);
  
  if( AnaChannel=="Electron2012LoPU")
  {
    //Barrel Barrel
    if((fabs(Z.Lep1etaSC) >= 0.0 && fabs(Z.Lep1etaSC) < 1.4442) && (fabs(Z.Lep2etaSC) >= 0.0 && fabs(Z.Lep2etaSC) < 1.4442))
      h1_ZmassCorr_BB->Fill(Z.mass,mTTW);
    //Barrel Endcap
    if((fabs(Z.Lep1etaSC) >= 0.0 && fabs(Z.Lep1etaSC) < 1.4442) && (fabs(Z.Lep2etaSC) >= 1.566 && fabs(Z.Lep2etaSC) < 2.5))
      h1_ZmassCorr_BE->Fill(Z.mass,mTTW);
    if((fabs(Z.Lep1etaSC) >= 1.566 && fabs(Z.Lep1etaSC) < 2.5) && (fabs(Z.Lep2etaSC) >= 0.0 && fabs(Z.Lep2etaSC) < 1.4442))
      h1_ZmassCorr_BE->Fill(Z.mass,mTTW);
    //Endcap Endcap
    if((fabs(Z.Lep1etaSC) >= 1.566 && fabs(Z.Lep1etaSC) < 2.5) && (fabs(Z.Lep2etaSC) >= 1.566 && fabs(Z.Lep2etaSC) < 2.5))
      h1_ZmassCorr_EE->Fill(Z.mass,mTTW);
  }
  
  if( AnaChannel=="Muon2012LoPU")
  {
    //Barrel Barrel
    if((fabs(Z.Lep1etaSC) >= 0.0 && fabs(Z.Lep1etaSC) < 1.) && (fabs(Z.Lep2etaSC) >= 0.0 && fabs(Z.Lep2etaSC) < 1.))
      h1_ZmassCorr_BB->Fill(Z.mass,mTTW);
    //Barrel Endcap
    if((fabs(Z.Lep1etaSC) >= 0.0 && fabs(Z.Lep1etaSC) < 1.) && (fabs(Z.Lep2etaSC) >= 1. && fabs(Z.Lep2etaSC) < 2.4))
      h1_ZmassCorr_BE->Fill(Z.mass,mTTW);
    if((fabs(Z.Lep1etaSC) >= 1. && fabs(Z.Lep1etaSC) < 2.4) && (fabs(Z.Lep2etaSC) >= 0.0 && fabs(Z.Lep2etaSC) < 1.))
      h1_ZmassCorr_BE->Fill(Z.mass,mTTW);
    //Endcap Endcap
    if((fabs(Z.Lep1etaSC) >= 1. && fabs(Z.Lep1etaSC) < 2.4) && (fabs(Z.Lep2etaSC) >= 1. && fabs(Z.Lep2etaSC) < 2.4))
      h1_ZmassCorr_EE->Fill(Z.mass,mTTW);
    
    //========================================
    //Fill Zmass histograms at eta ranges
    //========================================
    for(int iBin(0);iBin<ScaleBins;iBin++){
      //========================================
      //Check Charge Plus
      //========================================
      if(Z.Lep1Charge>0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	{
          //========================================
          //Data: if low statistics
	  //========================================
	  //if(Mode == "ScaleMakeRD")
	  //{
	  //  if(Z.Lep1Pt > Z.Lep2Pt)
	  //  {
	  //    if(Z.Lep1etaSC < 0.35 || (Z.Lep1etaSC > 0.7 && Z.Lep1etaSC < 1.40))
	  //      h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //    if(Z.Lep1etaSC > 1.75)
	  //      h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //  }
	  //  
	  //  if(Z.Lep1Pt <= Z.Lep2Pt)
	  //    if(Z.Lep1etaSC < 1.75)
	  //      h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //}
	  //
	  //if(Mode == "ScaleMakeMC")
	  //  h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //========================================
	  
	  h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}

        //========================================
        //Fill Charge Plus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Charge>0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	{
	  //========================================
          //Data: if low statistics
	  //========================================
	  //if(Mode == "ScaleMakeRD")
	  //{
	  //  if(Z.Lep2Pt > Z.Lep1Pt)
	  //  {
	  //    if(Z.Lep2etaSC < 0.35 || (Z.Lep2etaSC > 0.7 && Z.Lep2etaSC < 1.40))
	  //      h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //    if(Z.Lep2etaSC > 1.75)
	  //      h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //  }

	  //  if(Z.Lep2Pt <= Z.Lep1Pt)
	  //    if(Z.Lep2etaSC < 1.75)
	  //      h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //}

	  //if(Mode == "ScaleMakeMC")
	  //  h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	  //========================================

	  h1_ZmassCorr_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
        //========================================
        //Fill Charge Plus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Charge<0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Charge<0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Fill Leading lepton histo
      //========================================
      if(Z.Lep1Pt>Z.Lep2Pt && Z.Lep1Charge>0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_LeadingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_LeadingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt>Z.Lep1Pt && Z.Lep2Charge>0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_LeadingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_LeadingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Pt>Z.Lep2Pt && Z.Lep1Charge<0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_LeadingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_LeadingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt>Z.Lep1Pt && Z.Lep2Charge<0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_LeadingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_LeadingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Fill Trailing lepton histo
      //========================================
      if(Z.Lep1Pt<Z.Lep2Pt && Z.Lep1Charge>0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_TrailingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_TrailingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt<Z.Lep1Pt && Z.Lep2Charge>0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_TrailingLept_muEtaP[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Plus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_TrailingLept_noOverLap_muEtaP[iBin]->Fill(Z.mass,mTTW);
	}
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Pt<Z.Lep2Pt && Z.Lep1Charge<0)
      {
	if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_TrailingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept1 noOverLap histo
        //========================================
	if(fabs(Z.Lep1etaSC) < 0.9 || fabs(Z.Lep1etaSC) > 1.2)
	{
	  if(Z.Lep1etaSC >= ZmassEtaBins[iBin] && Z.Lep1etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_TrailingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
      
      if(Z.Lep2Pt<Z.Lep1Pt && Z.Lep2Charge<0)
      {
	if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	  h1_ZmassCorr_TrailingLept_muEtaM[iBin]->Fill(Z.mass,mTTW);
	
        //========================================
        //Fill Charge Minus Lept2 noOverLap histo
        //========================================
	if(fabs(Z.Lep2etaSC) < 0.9 || fabs(Z.Lep2etaSC) > 1.2)
	{
	  if(Z.Lep2etaSC >= ZmassEtaBins[iBin] && Z.Lep2etaSC < ZmassEtaBins[iBin+1])
	    h1_ZmassCorr_TrailingLept_noOverLap_muEtaM[iBin]->Fill(Z.mass,mTTW);
	}
      }
    }
    
    //========================================
    //Fill Zmass histograms at pT ranges
    //========================================
    for(int iBin(0);iBin<4;iBin++){
      //Check Charge Plus
      //========================================
      if(Z.Lep1Charge>0)
      {
	if(iBin<3) if(Z.Lep1Pt >= ZmassPtBins[iBin] && Z.Lep1Pt < ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtP[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep1Pt >= ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtP[iBin]->Fill(Z.mass,mTTW);
      }
      
      if(Z.Lep2Charge>0)
      {
	if(iBin<3) if(Z.Lep2Pt >= ZmassPtBins[iBin] && Z.Lep2Pt < ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtP[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep2Pt >= ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtP[iBin]->Fill(Z.mass,mTTW);
      }

      //========================================
      //Check Charge Minus
      //========================================
      if(Z.Lep1Charge<0)
      {
	if(iBin<3) if(Z.Lep1Pt >= ZmassPtBins[iBin] && Z.Lep1Pt < ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtM[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep1Pt >= ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtM[iBin]->Fill(Z.mass,mTTW);
      }
      
      if(Z.Lep2Charge<0)
      {
	if(iBin<3) if(Z.Lep2Pt >= ZmassPtBins[iBin] && Z.Lep2Pt < ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtM[iBin]->Fill(Z.mass,mTTW);
	if(iBin==3) if(Z.Lep2Pt >= ZmassPtBins[iBin+1])
	  h1_ZmassCorr_muPtM[iBin]->Fill(Z.mass,mTTW);
      }
    }
  }
  return 0;
}

int Wlnu12LoScaleSmearCorr::Write_ZHisto()
{
  h1_Zmass->Write();
  h1_Zmass_BB->Write();
  h1_Zmass_BE->Write();
  h1_Zmass_EE->Write();
  h1_ZLep1Pt -> Write();
  h1_ZLep2Pt -> Write();
  h1_ZLepPt_p-> Write();
  h1_ZLepPt_m-> Write();
  for(int ieta=0;ieta<ScaleBins;ieta++)
  {
    h1_Zmass_muEtaP[ieta]->Write();
    h1_Zmass_muEtaM[ieta]->Write();
    
    h1_Zmass_noOverLap_muEtaP[ieta]->Write();
    h1_Zmass_noOverLap_muEtaM[ieta]->Write();
    
    h1_Zmass_LeadingLept_noOverLap_muEtaP[ieta]->Write();
    h1_Zmass_LeadingLept_noOverLap_muEtaM[ieta]->Write();
    
    h1_Zmass_LeadingLept_muEtaP[ieta]->Write();
    h1_Zmass_LeadingLept_muEtaM[ieta]->Write();
    
    h1_Zmass_TrailingLept_noOverLap_muEtaP[ieta]->Write();
    h1_Zmass_TrailingLept_noOverLap_muEtaM[ieta]->Write();
    
    h1_Zmass_TrailingLept_muEtaP[ieta]->Write();
    h1_Zmass_TrailingLept_muEtaM[ieta]->Write();
  }
 
  for(int ipt=0;ipt<4;ipt++)
  {
    h1_Zmass_muPtP[ipt]->Write();
    h1_Zmass_muPtM[ipt]->Write();
  }
  if(Mode == "ScaleMakeMC")
  {
    h1_ZmassCorr->Write();
    h1_ZmassCorr_BB->Write();
    h1_ZmassCorr_BE->Write();
    h1_ZmassCorr_EE->Write(); 
    for(int ieta=0;ieta<ScaleBins;ieta++)
    {
      h1_ZmassCorr_muEtaP[ieta]->Write();
      h1_ZmassCorr_muEtaM[ieta]->Write();
    
      h1_ZmassCorr_noOverLap_muEtaP[ieta]->Write();
      h1_ZmassCorr_noOverLap_muEtaM[ieta]->Write();
      
      h1_ZmassCorr_LeadingLept_noOverLap_muEtaP[ieta]->Write();
      h1_ZmassCorr_LeadingLept_noOverLap_muEtaM[ieta]->Write();
      
      h1_ZmassCorr_LeadingLept_muEtaP[ieta]->Write();
      h1_ZmassCorr_LeadingLept_muEtaM[ieta]->Write();
      
      h1_ZmassCorr_TrailingLept_noOverLap_muEtaP[ieta]->Write();
      h1_ZmassCorr_TrailingLept_noOverLap_muEtaM[ieta]->Write();
      
      h1_ZmassCorr_TrailingLept_muEtaP[ieta]->Write();
      h1_ZmassCorr_TrailingLept_muEtaM[ieta]->Write();
    }
    for(int ipt=0;ipt<4;ipt++)
    {
      h1_ZmassCorr_muPtP[ipt]->Write();
      h1_ZmassCorr_muPtM[ipt]->Write();
    }
  }
  
  if(Mode =="ScaleMakeMC" || Mode =="ScaleMakeRD")
  {
    if(AnaChannel == "Electron2012LoPU")
    {
      for(int i(0);i<ScElCombiBins;i++)
      {
	h1_ZmassDaughEtaEle[i]->Write();
      }
      for(int i(0);i<ScElCombiBinsDiag;i++)
      {
	h1_ZmassDaughEtaEleDiag[i]->Write();
      }
    }
    
    if(AnaChannel == "Muon2012LoPU")
    {
      for(int i(0);i<ScMuCombiBins;i++)
      {
	h1_ZmassDaughEtaMu[i]->Write();
      }
    }
  }
  return 0;
}

int Wlnu12LoScaleSmearCorr::Fill_EleZmassDaughEta(int etaRange1, int etaRange2)
{
  //ScaleSmear Electron 6 only Diagobal category Fill 
  if( (etaRange1==0) && (etaRange2==0))
    h1_ZmassDaughEtaEleDiag[0]->Fill(Z.mass);
  if((etaRange1==1) && (etaRange2==1))
    h1_ZmassDaughEtaEleDiag[1]->Fill(Z.mass);
  if((etaRange1==2) && (etaRange2==2))
    h1_ZmassDaughEtaEleDiag[2]->Fill(Z.mass);
  if((etaRange1==3) && (etaRange2==3))
    h1_ZmassDaughEtaEleDiag[3]->Fill(Z.mass);
  if((etaRange1==4 && etaRange2==4))
    h1_ZmassDaughEtaEleDiag[4]->Fill(Z.mass);
  if((etaRange1==5 && etaRange2==5))
    h1_ZmassDaughEtaEleDiag[5]->Fill(Z.mass);

  //ScaleSmear Electron 21 category Fill 
  if( (etaRange1==0) && (etaRange2==0))
    h1_ZmassDaughEtaEle[0]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==1) || (etaRange1==1 && etaRange2==0))
    h1_ZmassDaughEtaEle[1]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==2) || (etaRange1==2 && etaRange2==0))
    h1_ZmassDaughEtaEle[2]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==3) || (etaRange1==3 && etaRange2==0))
    h1_ZmassDaughEtaEle[3]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==4) || (etaRange1==4 && etaRange2==0))
    h1_ZmassDaughEtaEle[4]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==5) || (etaRange1==5 && etaRange2==0))
    h1_ZmassDaughEtaEle[5]->Fill(Z.mass);
  if((etaRange1==1) && (etaRange2==1))
    h1_ZmassDaughEtaEle[6]->Fill(Z.mass);
  if((etaRange1==1 && etaRange2==2) || (etaRange1==2 && etaRange2==1))
    h1_ZmassDaughEtaEle[7]->Fill(Z.mass);
  if((etaRange1==1 && etaRange2==3) || (etaRange1==3 && etaRange2==1))
    h1_ZmassDaughEtaEle[8]->Fill(Z.mass);
  if((etaRange1==1 && etaRange2==4) || (etaRange1==4 && etaRange2==1))
    h1_ZmassDaughEtaEle[9]->Fill(Z.mass);
  if((etaRange1==1 && etaRange2==5) || (etaRange1==5 && etaRange2==1))
    h1_ZmassDaughEtaEle[10]->Fill(Z.mass);
  if((etaRange1==2) && (etaRange2==2))
    h1_ZmassDaughEtaEle[11]->Fill(Z.mass);
  if((etaRange1==2 && etaRange2==3) || (etaRange1==3 && etaRange2==2))
    h1_ZmassDaughEtaEle[12]->Fill(Z.mass);
  if((etaRange1==2 && etaRange2==4) || (etaRange1==4 && etaRange2==2))
    h1_ZmassDaughEtaEle[13]->Fill(Z.mass);
  if((etaRange1==2 && etaRange2==5) || (etaRange1==5 && etaRange2==2))
    h1_ZmassDaughEtaEle[14]->Fill(Z.mass);
  if((etaRange1==3) && (etaRange2==3))
    h1_ZmassDaughEtaEle[15]->Fill(Z.mass);
  if((etaRange1==3 && etaRange2==4) ||( etaRange1==4 && etaRange2==3))
    h1_ZmassDaughEtaEle[16]->Fill(Z.mass);
  if((etaRange1==3 && etaRange2==5) || ( etaRange1==5 && etaRange2==3))
    h1_ZmassDaughEtaEle[17]->Fill(Z.mass);
  if((etaRange1==4 && etaRange2==4))
    h1_ZmassDaughEtaEle[18]->Fill(Z.mass);
  if((etaRange1==4 && etaRange2==5) || (etaRange1==5 && etaRange2==4))
    h1_ZmassDaughEtaEle[19]->Fill(Z.mass);
  if((etaRange1==5 && etaRange2==5))
    h1_ZmassDaughEtaEle[20]->Fill(Z.mass);

  return 0;
}

int Wlnu12LoScaleSmearCorr::Fill_MuZmassDaughEta(int etaRange1, int etaRange2)
{
  //ScaleSmear Muon 15 category Fill
  if((etaRange1==0) && (etaRange2==0))
    h1_ZmassDaughEtaMu[0]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==1) || (etaRange1==1 && etaRange2==0))
    h1_ZmassDaughEtaMu[1]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==2) || (etaRange1==2 && etaRange2==0))
    h1_ZmassDaughEtaMu[2]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==3) || (etaRange1==3 && etaRange2==0))
    h1_ZmassDaughEtaMu[3]->Fill(Z.mass);
  if((etaRange1==0 && etaRange2==4) || (etaRange1==4 && etaRange2==0))
    h1_ZmassDaughEtaMu[4]->Fill(Z.mass);
  if((etaRange1==1) && (etaRange2==1))
    h1_ZmassDaughEtaMu[5]->Fill(Z.mass);
  if((etaRange1==1 && etaRange2==2) || (etaRange1==2 && etaRange2==1))
    h1_ZmassDaughEtaMu[6]->Fill(Z.mass);
  if((etaRange1==1 && etaRange2==3) || (etaRange1==3 && etaRange2==1))
    h1_ZmassDaughEtaMu[7]->Fill(Z.mass);
  if((etaRange1==1 && etaRange2==4) || (etaRange1==4 && etaRange2==1))
    h1_ZmassDaughEtaMu[8]->Fill(Z.mass);
  if((etaRange1==2) && (etaRange2==2))
    h1_ZmassDaughEtaMu[9]->Fill(Z.mass);
  if((etaRange1==2 && etaRange2==3) || (etaRange1==3 && etaRange2==2))
    h1_ZmassDaughEtaMu[10]->Fill(Z.mass);
  if((etaRange1==2 && etaRange2==4) || (etaRange1==4 && etaRange2==2))
    h1_ZmassDaughEtaMu[11]->Fill(Z.mass);
  if((etaRange1==3) && (etaRange2==3))
    h1_ZmassDaughEtaMu[12]->Fill(Z.mass);
  if((etaRange1==3 && etaRange2==4) || (etaRange1==4 && etaRange2==3))
    h1_ZmassDaughEtaMu[13]->Fill(Z.mass);
  if((etaRange1==4) && (etaRange2==4))
    h1_ZmassDaughEtaMu[14]->Fill(Z.mass);

  return 0;
}

Int_t Wlnu12LoScaleSmearCorr::EtaRange(double lep1Eta)
{
  int lep1Range(-1);
  if( AnaChannel=="Electron2012LoPU")
  {
    if( fabs(lep1Eta) >= 0.0   && fabs(lep1Eta) < 0.4)    lep1Range=0;
    if( fabs(lep1Eta) >= 0.4   && fabs(lep1Eta) < 0.8)    lep1Range=1;
    if( fabs(lep1Eta) >= 0.8   && fabs(lep1Eta) < 1.2)    lep1Range=2;
    if( fabs(lep1Eta) >= 1.2   && fabs(lep1Eta) < 1.4442) lep1Range=3;
    if( fabs(lep1Eta) >= 1.566 && fabs(lep1Eta) < 2.0)    lep1Range=4;
    if( fabs(lep1Eta) >= 2.0   && fabs(lep1Eta) < 2.5)    lep1Range=5;
  }
  if( AnaChannel=="Muon2012LoPU" )
  {
    if( fabs(lep1Eta) >= 0.0   && fabs(lep1Eta) < 0.4)    lep1Range=0;
    if( fabs(lep1Eta) >= 0.4   && fabs(lep1Eta) < 0.8)    lep1Range=1;
    if( fabs(lep1Eta) >= 0.8   && fabs(lep1Eta) < 1.2)    lep1Range=2;
    if( fabs(lep1Eta) >= 1.2   && fabs(lep1Eta) < 1.6)    lep1Range=3;
    if( fabs(lep1Eta) >= 1.6   && fabs(lep1Eta) < 2.1)    lep1Range=4;
  }
  return lep1Range;
}

Double_t Wlnu12LoScaleSmearCorr::GetScaleCorr(double LepEta)
{
  if(AnaChannel == "Electron2012LoPU")
  {
    ///Scale to be applied on MC to check on Z 
    if(fabs(LepEta) >= 0.0   && fabs(LepEta) < 0.4) 	 {return  0.999315;}
    if(fabs(LepEta) >= 0.4   && fabs(LepEta) < 0.8) 	 {return  1.00358;} 
    if(fabs(LepEta) >= 0.8   && fabs(LepEta) < 1.2) 	 {return  1.00325;}
    if(fabs(LepEta) >= 1.2   && fabs(LepEta) < 1.4442)   {return  1.00244;}
    if(fabs(LepEta) >= 1.566 && fabs(LepEta) < 2.0) 	 {return  1.0067;}
    if(fabs(LepEta) >= 2.0   && fabs(LepEta) < 2.5)	 {return  0.992737;}
  }

  if(AnaChannel == "Muon2012LoPU")
  {
    if(fabs(LepEta) >= 0.0 && fabs(LepEta) < 0.4) {return 1.0025;}
    if(fabs(LepEta) >= 0.4 && fabs(LepEta) < 0.8) {return 1.0022;}
    if(fabs(LepEta) >= 0.8 && fabs(LepEta) < 1.2) {return 0.9995;}
    if(fabs(LepEta) >= 1.2 && fabs(LepEta) < 1.6) {return 1.0018;}
    if(fabs(LepEta) >= 1.6 && fabs(LepEta) < 2.1) {return 1.0007;}
  }

}

Double_t Wlnu12LoScaleSmearCorr::GetSmearCorr(double LepEta)
{
  //double EndcapSystFactor =3.0;
  if(AnaChannel == "Electron2012LoPU")
  {
    ///21 category result: smear to be applied on MC to check on Z 
    if(fabs(LepEta) >= 0.0   && fabs(LepEta) < 0.4)	 {return  0.382443;}
    if(fabs(LepEta) >= 0.4   && fabs(LepEta) < 0.8)	 {return  0.356171;}
    if(fabs(LepEta) >= 0.8   && fabs(LepEta) < 1.2) 	 {return  0.559123;}
    if(fabs(LepEta) >= 1.2   && fabs(LepEta) < 1.4442)   {return  0.01;}
    //if(fabs(LepEta) >= 1.566 && fabs(LepEta) < 2.0)  	 {return  0.972944*EndcapSystFactor ;}
    //if(fabs(LepEta) >= 2.0   && fabs(LepEta) < 2.5)  	 {return  1.84788*EndcapSystFactor;}
    if(fabs(LepEta) >= 1.566 && fabs(LepEta) < 2.0)  	 {return  0.972944*2.5 ;}
    if(fabs(LepEta) >= 2.0   && fabs(LepEta) < 2.5)  	 {return  1.84788*2.5;}
    
    ///6 category result: smear to be applied on MC to check on Z 
    //if(fabs(LepEta) >= 0.0   && fabs(LepEta) < 0.4)	 {return  0.236791;}
    //if(fabs(LepEta) >= 0.4   && fabs(LepEta) < 0.8)	 {return  0.248012;}
    //if(fabs(LepEta) >= 0.8   && fabs(LepEta) < 1.2) 	 {return  0.334847;}
    //if(fabs(LepEta) >= 1.2   && fabs(LepEta) < 1.4442)   {return  1.39538;}
    //if(fabs(LepEta) >= 1.566 && fabs(LepEta) < 2.0)  	 {return  0.496947;}
    //if(fabs(LepEta) >= 2.0   && fabs(LepEta) < 2.5)  	 {return  1.88972;}
  }

  if(AnaChannel == "Muon2012LoPU")
  {
    if(fabs(LepEta) >= 0.0 && fabs(LepEta) < 0.4) {return  0.01; } 
    if(fabs(LepEta) >= 0.4 && fabs(LepEta) < 0.8) {return  0.381253;}
    if(fabs(LepEta) >= 0.8 && fabs(LepEta) < 1.2) {return  0.743902;}
    if(fabs(LepEta) >= 1.2 && fabs(LepEta) < 1.6) {return  0.637325;}
    if(fabs(LepEta) >= 1.6 && fabs(LepEta) < 2.1) {return  0.611946;}
  }
}



