#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include "TStyle.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "TMath.h"
#include <iomanip>

const TString format("pdf");
void FEWZ_PDFUncer(double *WmWpRatio, double *Err);
void Powheg_PDFUncer(double *WmWpRatio, double *Err);

double BinWidth[14] ={0.0, 7.5-0, 12.5-7.5,17.5-12.5, 24.0-17.5, 30.0-24.0, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0     , 190.0-150.0, 250.0-190.0, 600.0-250.0};
const int nBins = 13;
//double WptBins[nBins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double WptBins[nBins] = {1.,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

void MuEleCombineRatio()
{
  gROOT->LoadMacro("../../Utils/tdrstyle.C");
  setTDRStyle();
  gROOT->LoadMacro("../../Utils/CMS_lumi.C");
  //writeExtraText = "true";
  writeExtraText = false;
  extraText = "Preliminary";
  lumi_8TeV = "18.4 pb^{-1}";

  int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
  int iPos = 0;

  int W = 800;
  int H = 800;

  //
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
  // For instance:
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  //
  int H_ref = 600;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  // Read root file
  TFile* f_WpToMuNu = new TFile("./WpToMuNu_CommonFid.root");
  TFile* f_WmToMuNu = new TFile("./WmToMuNu_CommonFid.root");
  TFile* f_WpToEleNu = new TFile("./WpToEleNu_CommonFid.root");
  TFile* f_WmToEleNu = new TFile("./WmToEleNu_CommonFid.root");

  // Read histogram from input root file
  TH1D* h_WpToMuNu = (TH1D*)f_WpToMuNu->Get("h_NormDiffXsec_CommonFid");
  TH1D* h_WmToMuNu = (TH1D*)f_WmToMuNu->Get("h_NormDiffXsec_CommonFid");
  TH1D* h_WpToEleNu = (TH1D*)f_WpToEleNu->Get("h_NormDiffXsec_CommonFid");
  TH1D* h_WmToEleNu = (TH1D*)f_WmToEleNu->Get("h_NormDiffXsec_CommonFid");

  // Define combined histogram 
  TH1D* h_CombP = new TH1D("h_CombP","h_CombP",12,0,12);
  TH1D* h_CombM = new TH1D("h_CombM","h_CombM",12,0,12);

  // Define histogram to draw
  TH1D* hWmWpratio_RD = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  // Array for calculation
  double WpToMuNu[13] = {0};
  double WmToMuNu[13] = {0};
  double WpToEleNu[13] = {0};
  double WmToEleNu[13] = {0};
  double WpToMuNuErr[13] = {0};
  double WmToMuNuErr[13] = {0};
  double WpToEleNuErr[13] = {0};
  double WmToEleNuErr[13] = {0};

  double CombP[13] = {0};
  double CombM[13] = {0};
  double CombPErr[13] = {0};
  double CombMErr[13] = {0};
  double CombErr[13] = {0};
  double tmp_CombPErr[13] = {0};
  double tmp_CombMErr[13] = {0};
  double tmp_CombErr[13] = {0};
  for(int i(1); i<13; i++)
  {
    WpToMuNu[i] = h_WpToMuNu->GetBinContent(i);
    WmToMuNu[i] = h_WmToMuNu->GetBinContent(i);
    WpToEleNu[i] = h_WpToEleNu->GetBinContent(i);
    WmToEleNu[i] = h_WmToEleNu->GetBinContent(i);
    WpToMuNuErr[i] = h_WpToMuNu->GetBinError(i);
    WmToMuNuErr[i] = h_WmToMuNu->GetBinError(i);
    WpToEleNuErr[i] = h_WpToEleNu->GetBinError(i);
    WmToEleNuErr[i] = h_WmToEleNu->GetBinError(i);
    //cout << "WpToMuNu : " << WpToMuNu[i] << "\t WpToMuNuErr : " << WpToMuNuErr[i] << endl;

    // Mu,Ele combine
    tmp_CombPErr[i] = 1.0/(WpToMuNuErr[i]*WpToMuNuErr[i]) + 1.0/(WpToEleNuErr[i]*WpToEleNuErr[i]);
    tmp_CombMErr[i] = 1.0/(WmToMuNuErr[i]*WmToMuNuErr[i]) + 1.0/(WmToEleNuErr[i]*WmToEleNuErr[i]);
    CombP[i] = (WpToMuNu[i]/(WpToMuNuErr[i]*WpToMuNuErr[i])+WpToEleNu[i]/(WpToEleNuErr[i]*WpToEleNuErr[i]))/tmp_CombPErr[i]; 
    CombM[i] = (WmToMuNu[i]/(WmToMuNuErr[i]*WmToMuNuErr[i])+WmToEleNu[i]/(WmToEleNuErr[i]*WmToEleNuErr[i]))/tmp_CombMErr[i]; 
    CombPErr[i] = sqrt(1.0/tmp_CombPErr[i]); 
    CombMErr[i] = sqrt(1.0/tmp_CombMErr[i]); 

    cout << "CombP : " << CombP[i] << "\t CombPErr : " << CombPErr[i] << "\t CombM[i] : " << CombM[i] << "\t CombMErr : " << CombMErr[i] << endl;

  } // Data related thing finished

  ///Powheg
  TFile *f_WmToMu_Powheg = new TFile("../Wpt_PowhegPreFSR_12Bin/root/WmToMuNu_PowhegWpt.root");
  TFile *f_WpToMu_Powheg = new TFile("../Wpt_PowhegPreFSR_12Bin/root/WpToMuNu_PowhegWpt.root");
  TH1D* hWmToMu_Powheg = (TH1D*)f_WmToMu_Powheg->Get("hxsec_NormDiff")->Clone("hWmToMu_Powheg");
  TH1D* hWpToMu_Powheg = (TH1D*)f_WpToMu_Powheg->Get("hxsec_NormDiff")->Clone("hWpToMu_Powheg");
  TH1D* hWmWpratio_Powheg = new TH1D("W/Z ratio", "", nBins-1,WptBins);

  ///ResBos
  TFile *f_WmToMu_resbos = new TFile("../WpT_Resbos_12Bin/root/WmToMuNu_Resbos.root");
  TFile *f_WpToMu_resbos = new TFile("../WpT_Resbos_12Bin/root/WpToMuNu_Resbos.root");

  TH1D* hWmToMu_Resbos[7];
  TH1D* hWpToMu_Resbos[7];
  TH1D* hWmWpratio_Resbos[7];
  TH1D* hWmWpratio_Resbos_Central = new TH1D("W-/W+ ratio","",nBins-1,WptBins);
  char gridname[30];
  char histName[30];
  for(int i(0);i<7; i++)
  {
    sprintf(gridname,"hResbos%d",i+29);
    sprintf(histName,"hWmToMu_Resbos_%d",i);
    hWmToMu_Resbos[i] = (TH1D*)f_WmToMu_resbos->Get(gridname)->Clone(histName);
    sprintf(histName,"hWpToMu_Resbos_%d",i);
    hWpToMu_Resbos[i] = (TH1D*)f_WpToMu_resbos->Get(gridname)->Clone(histName);
    sprintf(histName,"hWmWpratio_Resbos%d",i);
    hWmWpratio_Resbos[i] = (TH1D*)hWmToMu_Resbos[i]->Clone(histName);
    hWmWpratio_Resbos[i] -> Divide(hWpToMu_Resbos[i]);
  }
  
  ///FEWZ
  TFile *f_WmToMu_fewz = new TFile("../WpT_FEWZ_12Bin/root/WmToMuNu_FEWZ.root");
  TFile *f_WpToMu_fewz = new TFile("../WpT_FEWZ_12Bin/root/WpToMuNu_FEWZ.root");
  TH1D* hWmToMu_FEWZ = (TH1D*)f_WmToMu_fewz->Get("hxsec_NormDiff")->Clone("hWmToMu_FEWZ");
  TH1D* hWpToMu_FEWZ = (TH1D*)f_WpToMu_fewz->Get("hxsec_NormDiff")->Clone("hWpToMu_FEWZ");
  TH1D* hWmToMu_Up_FEWZ = (TH1D*)f_WmToMu_fewz->Get("hxsec_NormDiff_up")->Clone("hWmToMu_Up_FEWZ");
  TH1D* hWmToMu_Down_FEWZ = (TH1D*)f_WmToMu_fewz->Get("hxsec_NormDiff_down")->Clone("hWmToMu_Down_FEWZ");
  TH1D* hWpToMu_Up_FEWZ = (TH1D*)f_WpToMu_fewz->Get("hxsec_NormDiff_up")->Clone("hWpToMu_Up_FEWZ");
  TH1D* hWpToMu_Down_FEWZ = (TH1D*)f_WpToMu_fewz->Get("hxsec_NormDiff_down")->Clone("hWpToMu_Down_FEWZ");
  TH1D* hWmWpratio_FEWZ = new TH1D("W+/W- ratio", "", nBins-1,WptBins);

  double WmToMu_Powheg[13]={0};
  double WpToMu_Powheg[13]={0};
  double WmToMuStatErr_Powheg[13]={0};
  double WpToMuStatErr_Powheg[13]={0};
 
  double WmToMu_Resbos[13]={0};
  double WpToMu_Resbos[13]={0};
  double WmToMuErr_Resbos[13]={0};
  double WpToMuErr_Resbos[13]={0};
 
  double WmToMu_FEWZ[13]={0};
  double WpToMu_FEWZ[13]={0};
  double WmToMuStatErr_FEWZ[13]={0};
  double WpToMuStatErr_FEWZ[13]={0};
 
  double WmWpratio_RD[13]={0};
  double WmWpratioErr_RD_NoCorr[13]={0};
  double WmWpratioErr_RD_FullCorr[13]={0};
  double WmWpratioErr_RD_Total[13]={0};

  double WmWpratio_Powheg[13]={0};
  double WmWpratioStatErr_Powheg[13]={0}; // Stat error by error propagation
  double WmWpratioPDFErr_Powheg[13]={0}; // W-/W+ PDF error
  double WmWpratioTotalErr_Powheg[13]={0}; // sqrt(StatErr^2 + PDFErr^2 + ScaleErr^2)
  
  double WmWpratio_FEWZ[13]={0};
  double WmWp_Up_ratio_FEWZ[13]={0};
  double WmWp_Down_ratio_FEWZ[13]={0};
  double WmWpratioStatErr_FEWZ[13]={0}; // Stat error by error propagation
  double WmWpratioPDFErr_FEWZ[13]={0}; // W-/W+ PDF error
  double WmWpratioScaleErr_FEWZ[13]={0}; // W-/W+ scale error
  double WmWpratioTotalErr_FEWZ[13]={0}; // sqrt(StatErr^2 + PDFErr^2 + ScaleErr^2)
 
  Double_t Resb_errMax[nBins];
  Double_t Resb_errMin[nBins];
  Double_t Resb_err[nBins];
 
  cout << fixed << setprecision(10) << endl;
  cout << " ===== Wpt minus and Wptplus Normalilzed Differntial cross-section and errors In Fiducial Volume ==="<< endl;
  cout << "Bin"<<
    "\t Wpt minus  " << " \t\tError " << 
    "\t\t\t Wpt plus  " << " \t\tError " << endl;
  for(int i=1; i<13; i++)
  {
    // RD Ratio and Error
    WmWpratio_RD[i] = CombM[i]/CombP[i];
    //tmp_CombErr[i] = 1.0/(CombMErr[i]*CombMErr[i]) + 1.0/(CombPErr[i]*CombPErr[i]);
    //CombErr[i] = sqrt(1.0/tmp_CombErr[i]);
    WmWpratioErr_RD_Total[i] = WmWpratio_RD[i] * sqrt(
	  CombMErr[i]*CombMErr[i]/CombM[i]/CombM[i]
	+ CombPErr[i]*CombPErr[i]/CombP[i]/CombP[i]
	);

    ///Powheg & Resbos ratio
    WmWpratio_Powheg[i] = hWmToMu_Powheg->GetBinContent(i) /hWpToMu_Powheg->GetBinContent(i) ; 
    //WmWpratio_Resbos[i] = hWmToMu_Resbos->GetBinContent(i) /hWpToMu_Resbos->GetBinContent(i) ; 
    WmWpratio_FEWZ[i] = hWmToMu_FEWZ->GetBinContent(i) /hWpToMu_FEWZ->GetBinContent(i) ; 
    WmWp_Up_ratio_FEWZ[i] = hWmToMu_Up_FEWZ->GetBinContent(i) /hWpToMu_Up_FEWZ->GetBinContent(i) ; 
    WmWp_Down_ratio_FEWZ[i] = hWmToMu_Down_FEWZ->GetBinContent(i) /hWpToMu_Down_FEWZ->GetBinContent(i) ; 
    WmWpratioScaleErr_FEWZ[i] = TMath::Max(fabs(WmWp_Up_ratio_FEWZ[i]-WmWpratio_FEWZ[i]),fabs(WmWp_Down_ratio_FEWZ[i]-WmWpratio_FEWZ[i]));

    // FEWZ ratio error propagation
    WmToMu_FEWZ[i] = hWmToMu_FEWZ->GetBinContent(i);
    WmToMuStatErr_FEWZ[i] = hWmToMu_FEWZ->GetBinError(i);
    
    WpToMu_FEWZ[i] = hWpToMu_FEWZ->GetBinContent(i);
    WpToMuStatErr_FEWZ[i] = hWpToMu_FEWZ->GetBinError(i);

    WmWpratioStatErr_FEWZ[i] = WmWpratio_FEWZ[i]* TMath::Sqrt(
	  (WmToMuStatErr_FEWZ[i]*WmToMuStatErr_FEWZ[i]/WmToMu_FEWZ[i]/WmToMu_FEWZ[i]) 
	+ (WpToMuStatErr_FEWZ[i]*WpToMuStatErr_FEWZ[i]/WpToMu_FEWZ[i]/WpToMu_FEWZ[i])
	); // FEWZ Stat Error propagation

    FEWZ_PDFUncer(WmWpratio_FEWZ,WmWpratioPDFErr_FEWZ);
    //cout << "WmWpratioPDFErr_FEWZ in main : " << WmWpratioPDFErr_FEWZ[i] << endl; // check PDF error number is correctly retured
 
    WmWpratioTotalErr_FEWZ[i] = sqrt(
	  WmWpratioStatErr_FEWZ[i]*WmWpratioStatErr_FEWZ[i] 
	+ WmWpratioPDFErr_FEWZ[i]*WmWpratioPDFErr_FEWZ[i] 
	+ WmWpratioScaleErr_FEWZ[i]*WmWpratioScaleErr_FEWZ[i]
	);

    // Powheg ratio error propagation
    WmToMu_Powheg[i] = hWmToMu_Powheg->GetBinContent(i);
    WmToMuStatErr_Powheg[i] = hWmToMu_Powheg->GetBinError(i);
    
    WpToMu_Powheg[i] = hWpToMu_Powheg->GetBinContent(i);
    WpToMuStatErr_Powheg[i] = hWpToMu_Powheg->GetBinError(i);

    WmWpratioStatErr_Powheg[i] =  WmWpratio_Powheg[i]* TMath::Sqrt(
	  (WmToMuStatErr_Powheg[i]*WmToMuStatErr_Powheg[i]/WmToMu_Powheg[i]/WmToMu_Powheg[i]) 
	+ (WpToMuStatErr_Powheg[i]*WpToMuStatErr_Powheg[i]/WpToMu_Powheg[i]/WpToMu_Powheg[i])
	);
 
    
    Powheg_PDFUncer(WmWpratio_Powheg,WmWpratioPDFErr_Powheg);
    //cout << "WmWpratioPDFErr_Powheg in main : " << WmWpratioPDFErr_Powheg[i] << "\t in %\t"<< 100*WmWpratioPDFErr_Powheg[i]/WmWpratio_Powheg[i]<<endl; // check PDF error number is correctly retured
    
    ///Total powheg error calculated here
    WmWpratioTotalErr_Powheg[i] = sqrt(WmWpratioStatErr_Powheg[i]*WmWpratioStatErr_Powheg[i] +WmWpratioPDFErr_Powheg[i]*WmWpratioPDFErr_Powheg[i]);

cout << "hahahahaha" << endl;
     // Resbos ratio error in normdiff stage
    double tmpVal,tmpDiff;
    double nomVal = hWmWpratio_Resbos[1]->GetBinContent(i);

    Resb_errMax[i] = -99999;
    Resb_errMin[i] = 990009;

    for (int j(0);j<7;j++)
    {
      tmpVal  = hWmWpratio_Resbos[j]->GetBinContent(i);
      tmpDiff = tmpVal - nomVal;
      if( tmpDiff > Resb_errMax[i]) Resb_errMax[i] = tmpDiff;
      if( tmpDiff < Resb_errMin[i]) Resb_errMin[i] = tmpDiff;
    }
//    cout << "Resbos center ratio : " << hWmWpratio_Resbos[1]->GetBinContent(i) << "\t error+ : " << Resb_errMax[i] << "\t error - : " << Resb_errMin[i] << endl;

    if (Resb_errMax[i] < 0) Resb_errMax[i] = 0.;
    if (Resb_errMin[i] > 0) Resb_errMin[i] = 0.;
    if (Resb_errMin[i] < 0) Resb_errMin[i] = -Resb_errMin[i];

    Resb_err[i] = TMath::Max(Resb_errMax[i],Resb_errMin[i]);
    //cout<<i<<" Bin Resbos ratio : " << hWmWpratio_Resbos[1]->GetBinContent(i) << " \t Error : " <<Resb_errMin[i]<<"\t"<<Resb_errMax[i]<<"\t"<<Resb_err[i] <<endl; 
    

    //printf("WmWpratioErr_FEWZ : %.8f\n",WmWpratioErr_FEWZ[i]);
    //printf("WmWpratioErr_Powheg : %.8f\n",WmWpratioErr_Powheg[i]);
    //printf("WmWpratioErr_Resbos : %.8f\n",WmWpratioErr_Resbos[i]);
  
  }

  cout << " \\"<<"\\" << endl;
  cout << endl;
 
  //// Print Powheg errors
  cout <<"=====  Print Powheg errors "<< endl;
  cout << fixed << setprecision(2) << endl;
  for(int i=1;i<13;i++)
  {
    cout << i<<"\tWmWpratioStatErr_Powheg: " << WmWpratioStatErr_Powheg[i] <<"\t in %\t"<< 100*WmWpratioStatErr_Powheg[i]/WmWpratio_Powheg[i]<<endl; 
  }
  for(int i=1;i<13;i++)
  {
    cout << i<<"\tWmWpratioPDFErr_Powheg : " << WmWpratioPDFErr_Powheg[i] << "\t in %\t"<< 100*WmWpratioPDFErr_Powheg[i]/WmWpratio_Powheg[i]<<endl;
  }
 //// Print FEWZ errors
  cout <<"=====  Print FEWZ errors "<< endl;

  for(int i=1;i<13;i++)
  {
    cout << i<<"\tWmWpratioStatErr_FEWZ : " << WmWpratioStatErr_FEWZ[i] << "\t in %\t"<< 100*WmWpratioStatErr_FEWZ[i]/WmWpratio_Powheg[i]<<endl;
  }

  for(int i=1;i<13;i++)
  {
    cout << i<<"\tWmWpratioPDFErr_FEWZ : " << WmWpratioPDFErr_FEWZ[i] << "\t in %\t"<< 100*WmWpratioPDFErr_FEWZ[i]/WmWpratio_Powheg[i]<<endl;
  }
  for(int i=1;i<13;i++)
  {
    cout << i<<"\tWmWpratioScaleErr_FEWZ : " << WmWpratioScaleErr_FEWZ[i] << "\t in %\t"<< 100*WmWpratioScaleErr_FEWZ[i]/WmWpratio_Powheg[i]<<endl;
  }

  cout << fixed << setprecision(5) << endl;
  cout<<"Wpt minus FEWZ\t"<<"Wpt plus FEWZ"<<endl; 
  for(int i=1;i<13;i++)
  {
   cout<<hWmToMu_FEWZ->GetBinContent(i) <<"\t"<<hWpToMu_FEWZ->GetBinContent(i)<<endl ; 
  }
  
  cout << fixed << setprecision(3) << endl;
  cout << " ============ W-/W+ ratio and errors(Numbers) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error " << endl;
  for(int i=1;i<13;i++)
  {
    cout <<i<< "\t" <<  WmWpratio_RD[i] <<  " \\pm " << WmWpratioErr_RD_Total[i] <<  endl;
  }
  cout << " ============ W-/W+ ratio and errors(in %) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error %" << endl;
  for(int i=1;i<13;i++)
  {

    cout <<i<< "\t" <<  WmWpratio_RD[i] <<  " \\pm " << WmWpratioErr_RD_Total[i]/WmWpratio_RD[i]*100 <<  endl;
  }

  for(int i=1;i<13;i++)
  {
    hWmWpratio_RD->SetBinContent(i,WmWpratio_RD[i]);
    hWmWpratio_RD->SetBinError(i,WmWpratioErr_RD_Total[i]);
    
    hWmWpratio_Powheg->SetBinContent(i,WmWpratio_Powheg[i]);
    hWmWpratio_Powheg->SetBinError(i,WmWpratioTotalErr_Powheg[i]);
    
    hWmWpratio_Resbos_Central->SetBinContent(i,hWmWpratio_Resbos[1]->GetBinContent(i));
    hWmWpratio_Resbos_Central->SetBinError(i,Resb_err[i]);
    //hWmWpratio_Resbos[1]->SetBinError(i,Resb_err[i]);
    
    hWmWpratio_FEWZ->SetBinContent(i,WmWpratio_FEWZ[i]);
    hWmWpratio_FEWZ->SetBinError(i,WmWpratioTotalErr_FEWZ[i]);
  }

  // Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataTotalErr = new TH1D("hRatioDataTotalErr","hRatioDataTotalErr",nBins-1,WptBins);
  TH1D *hRatioPowhegTotalErr = new TH1D("hRatioPowhegTotalErr","hRatioPowhegTotalErr",nBins-1,WptBins);
  TH1D *hRatioResbosTotalErr = new TH1D("hRatioResbosTotalErr","hRatioResbosTotalErr",nBins-1,WptBins);
  TH1D *hRatioFEWZTotalErr = new TH1D("hRatioFEWZTotalErr","hRatioFEWZTotalErr",nBins-1,WptBins);

  for(int i(1); i<nBins; i++)
  {
    hRatioDataTotalErr->SetBinContent(i,1.);
    hRatioDataTotalErr->SetBinError(i,WmWpratioErr_RD_Total[i] / WmWpratio_RD[i]);
    
    hRatioPowhegTotalErr->SetBinContent(i,WmWpratio_Powheg[i] / WmWpratio_RD[i]);
    hRatioPowhegTotalErr->SetBinError(i,WmWpratioTotalErr_Powheg[i] / WmWpratio_RD[i]);

    //hRatioResbosTotalErr->SetBinContent(i,hWmWpratio_Resbos[1]->GetBinContent(i) / WmWpratio_RD[i]);
    //hRatioResbosTotalErr->SetBinError(i,hWmWpratio_Resbos[1]->GetBinError(i) / WmWpratio_RD[i]);
    hRatioResbosTotalErr->SetBinContent(i,hWmWpratio_Resbos_Central->GetBinContent(i) / WmWpratio_RD[i]);
    hRatioResbosTotalErr->SetBinError(i,hWmWpratio_Resbos_Central->GetBinError(i) / WmWpratio_RD[i]);
    
    hRatioFEWZTotalErr->SetBinContent(i,WmWpratio_FEWZ[i] / WmWpratio_RD[i]);
    hRatioFEWZTotalErr->SetBinError(i,WmWpratioTotalErr_FEWZ[i] / WmWpratio_RD[i]);
  }

  //Color Transparent
  TColor *colRed = gROOT->GetColor(kRed);
  TColor *colBlue = gROOT->GetColor(kBlue);
  TColor *colGreen = gROOT->GetColor(kGreen);
  colRed->SetAlpha(0.2);
  colBlue->SetAlpha(0.2);
  colGreen->SetAlpha(0.2);
  //gStyle->SetOptStat(0); 

  // Draw 
  TGraphErrors *tgWmWpratio_Powheg = new TGraphErrors(hWmWpratio_Powheg);
  //TGraphErrors *tgWmWpratio_Resbos = new TGraphErrors(hWmWpratio_Resbos[1]);
  TGraphErrors *tgWmWpratio_Resbos = new TGraphErrors(hWmWpratio_Resbos_Central);
  TGraphErrors *tgWmWpratio_FEWZ = new TGraphErrors(hWmWpratio_FEWZ);
  
  TLegend *L1 = new TLegend(0.25,0.65,0.5,0.95);
  L1->SetFillColor(0);
  L1->SetBorderSize(0);
  L1->AddEntry(hWmWpratio_RD,"data","PL");
  L1->AddEntry(tgWmWpratio_Powheg,"Powheg","f");
  L1->AddEntry(tgWmWpratio_Resbos,"ResBos","f");
  L1->AddEntry(tgWmWpratio_FEWZ,"FEWZ","f");
  
  TCanvas *C1 = new TCanvas("can","can",50,50,W,H);
  CMS_lumi(C1,iPeriod,iPos);
  C1->Divide(1,2,0,0);
  C1->cd(1)->SetPad(0,0.52,0.98,0.95);
  C1->cd(1)->SetLogx();
  C1->cd(1)->SetTickx(1);
  C1->cd(1)->SetTicky(1);
 
  TPaveText *Wchannel = new TPaveText(0.2,0.35,0.4,0.35,"NDC");
  Wchannel->SetBorderSize(0);
  Wchannel->SetFillStyle(0);
  Wchannel->SetTextSize(0.04);
  Wchannel->SetTextColor(kBlue);
  Wchannel->AddText("W^{+} #rightarrow #mu^{+} #nu");
  
  TPaveText *Zchannel = new TPaveText(0.2,0.30,0.4,0.30,"NDC");
  Zchannel->SetBorderSize(0);
  Zchannel->SetFillStyle(0);
  Zchannel->SetTextSize(0.04);
  Zchannel->SetTextColor(kBlue);
  Zchannel->AddText("W^{-} #rightarrow #mu^{-} #bar{#nu}");
  
  TPaveText *LeptonCut = new TPaveText(0.45,0.32,0.65,0.32,"NDC");
  LeptonCut->SetBorderSize(0);
  LeptonCut->SetFillStyle(0);
  LeptonCut->SetTextSize(0.04);
  LeptonCut->SetTextColor(kBlue);
  LeptonCut->AddText("p_{T}>20 Gev, |#eta |<2.1");

  hWmWpratio_RD->GetYaxis()->SetRangeUser(0.5,1.5);
  hWmWpratio_RD->GetYaxis()->SetTitleOffset(0.8);
  hWmWpratio_RD->GetYaxis()->SetLabelSize(0.07);
  hWmWpratio_RD->GetYaxis()->SetTitleSize(0.07);
  //hWmWpratio_RD->GetYaxis()->SetTitle("(#frac{1}{#sigma^{W^{-}}} #frac{d#sigma^{W^{-}}}{p_{T}^{W^{-}}})/(#frac{1}{#sigma^{W^{+}}} #frac{d#sigma^{W^{+}}}{p_{T}^{W^{+}}})");
  hWmWpratio_RD->GetYaxis()->SetTitle("(#frac{1}{#sigma^{W^{#font[122]{\55}}}} #frac{d#sigma^{W^{#font[122]{\55}}}}{p_{T}^{W^{#font[122]{\55}}}})/(#frac{1}{#sigma^{W^{+}}} #frac{d#sigma^{W^{+}}}{p_{T}^{W^{+}}})");
  hWmWpratio_RD->GetYaxis()->SetNdivisions(405);
  hWmWpratio_RD->GetXaxis()->SetTitleOffset(1.);
  hWmWpratio_RD->GetXaxis()->SetLabelSize(0.04);
  hWmWpratio_RD->GetXaxis()->SetTitle("");
  hWmWpratio_RD->SetStats(0);
  hWmWpratio_RD->SetMarkerStyle(20);
  hWmWpratio_RD->SetMarkerColor(kBlack);
  hWmWpratio_RD->Draw("E1");

  tgWmWpratio_FEWZ->SetFillColor(kGreen);
  tgWmWpratio_FEWZ->SetLineColor(kGreen+2);
  tgWmWpratio_FEWZ->Draw("5");

  tgWmWpratio_Powheg->SetFillColor(kRed);
  tgWmWpratio_Powheg->SetLineColor(kRed+2);
  tgWmWpratio_Powheg->Draw("5");

  tgWmWpratio_Resbos->SetFillColor(kBlue);
  tgWmWpratio_Resbos->SetLineColor(kBlue+2);
  tgWmWpratio_Resbos->Draw("5");

  L1->Draw();
  //Wchannel->Draw();
  //Zchannel->Draw();
  //LeptonCut->Draw();
  hWmWpratio_RD->Draw("axis same");
  
  TLine *Line_XminYMax = new TLine(1.,0.5,1,1.5);
  Line_XminYMax->SetLineWidth(2);
  Line_XminYMax->SetLineColor(kBlack);
  Line_XminYMax->SetLineStyle(1);
  Line_XminYMax->Draw("same");
  
  TLine *Line_XminXMax = new TLine(1.,0.5,600,0.5);
  Line_XminXMax->SetLineWidth(2);
  Line_XminXMax->SetLineColor(kBlack);
  Line_XminXMax->SetLineStyle(1);
  Line_XminXMax->Draw("same");
  
  TLine *Line_XmaxYMin = new TLine(600.,0.5,600,1.5);
  Line_XmaxYMin->SetLineWidth(2);
  Line_XmaxYMin->SetLineColor(kBlack);
  Line_XmaxYMin->SetLineStyle(1);
  Line_XmaxYMin->Draw("same");
  
  TLine *Line_XmaxYMax = new TLine(1.,1.5,600,1.5);
  Line_XmaxYMax->SetLineWidth(2);
  Line_XmaxYMax->SetLineColor(kBlack);
  Line_XmaxYMax->SetLineStyle(1);
  Line_XmaxYMax->Draw("same");
  
  //gStyle->SetLineWidth(2); 
  gPad->RedrawAxis();

  C1->cd(2)->SetPad(0,0.1,0.98,0.5);
  C1->cd(2)->SetTickx(1);
  C1->cd(2)->SetTicky(1);
  C1->cd(2)->SetLogx(1);

  TH1D *hRatioDummy = new TH1D("hRatioDummy","",nBins-1,WptBins);
  
  TGraphErrors* tgRatioData = new TGraphErrors(hRatioDataTotalErr);
  TGraphErrors* tgRatioPowheg = new TGraphErrors(hRatioPowhegTotalErr);
  TGraphErrors* tgRatioResbos = new TGraphErrors(hRatioResbosTotalErr);
  TGraphErrors* tgRatioFEWZ = new TGraphErrors(hRatioFEWZTotalErr);
  
  // set canvas range and Y axis title
  hRatioDummy->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatioDummy->GetYaxis()->SetTitle("Theory / Data");
  hRatioDummy->GetYaxis()->CenterTitle();
  hRatioDummy->GetYaxis()->SetTitleSize(0.08);
  hRatioDummy->GetYaxis()->SetTitleOffset(0.78);
  hRatioDummy->GetYaxis()->SetLabelSize(0.09);
  hRatioDummy->GetYaxis()->SetNdivisions(605);
  
  hRatioDummy->GetXaxis()->SetTitle("p_{T}^{V} [GeV]");
  hRatioDummy->GetXaxis()->SetTitleOffset(0.6);
  hRatioDummy->GetXaxis()->SetTitleSize(0.08);
  hRatioDummy->GetXaxis()->SetLabelSize(0.07);

  // FEWZ Ratio plot setting
  tgRatioData->SetFillColor(kGray+2);
  tgRatioData->SetFillStyle(3354);

  tgRatioPowheg->SetFillColor(kRed);
  tgRatioPowheg->SetLineColor(kRed+2);
  tgRatioPowheg->SetMarkerStyle(22);
  tgRatioPowheg->SetMarkerColor(kRed+2);
  
  tgRatioResbos->SetFillColor(kBlue);
  tgRatioResbos->SetLineColor(kBlue+2);
  tgRatioResbos->SetMarkerStyle(20);
  tgRatioResbos->SetMarkerColor(kBlue+2);
  
  tgRatioFEWZ->SetFillColor(kGreen);
  tgRatioFEWZ->SetLineColor(kGreen+2);
  tgRatioFEWZ->SetMarkerStyle(21);
  tgRatioFEWZ->SetMarkerColor(kGreen+2);
  
  // Draw canvas
  hRatioDummy->Draw();
  tgRatioData->Draw("2");
  tgRatioFEWZ->Draw("5 P");
  tgRatioPowheg->Draw("5 P");
  tgRatioResbos->Draw("5 P");
  gPad->RedrawAxis();

  C1->SaveAs("WmMuWpMuNormFid12Bin."+format);
}
 
void FEWZ_PDFUncer(double *WmWpRatio, double *Err)
{
  // FEWZ Normalized PDF error in % unit
  double PDFErr[13] = {0.,};

  PDFErr[0] = 0.;  
  PDFErr[1] = 1.14;  
  PDFErr[2] = 0.36;
  PDFErr[3] = 0.59;
  PDFErr[4] = 0.58;
  PDFErr[5] = 0.80;
  PDFErr[6] = 1.06;
  PDFErr[7] = 1.23;
  PDFErr[8] = 1.58;
  PDFErr[9] = 1.80;
  PDFErr[10] = 2.52;
  PDFErr[11] = 2.24;
  PDFErr[12] = 3.87;

  for(int i(1); i<13; i++)
  {
    Err[i] = WmWpRatio[i] * PDFErr[i] * 0.01;
    //cout << "FEWZ PDFErr in function : " << Err[i] << endl;
  }
}
void Powheg_PDFUncer(double *WmWpRatio, double *Err)
{
  // Powheg Normalized PDF error in % unit
  double PDFErr[13] = {0.,};

  PDFErr[0] =0.0  ;  
  PDFErr[1] =0.34912  ;  
  PDFErr[2] =0.11815  ;
  PDFErr[3] =0.12652  ;
  PDFErr[4] =0.27171  ;
  PDFErr[5] =0.33522  ;
  PDFErr[6] =0.50769  ;
  PDFErr[7] =0.61623  ;
  PDFErr[8] =0.77008  ;
  PDFErr[9] =0.87044  ;
  PDFErr[10] =0.89737  ;
  PDFErr[11]=1.14706  ;
  PDFErr[12]=1.11763  ;

  for(int i(1); i<13; i++)
  {
    Err[i] = WmWpRatio[i] * PDFErr[i] * 0.01;
    //cout << "Powheg PDFErr in function : " << Err[i] << endl;
  }
}
