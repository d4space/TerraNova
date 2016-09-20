#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
using namespace std;

const TString format("pdf");

void FEWZ_PDFUncer(double *ZWRatio, double *Err);
void Powheg_PDFUncer(double *ZWRatio, double *Err);
char tmpName[30],tmpName_org[30];

int ZWratio_New()
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

  const int nBins = 13;
  double WptBins[nBins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
  //double WptBins[nBins] = {1.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

  ///Data
  TFile *f_WinclMu_RD = new TFile("../WptIncl_NormDiffXsec_InFid/Wpt_NormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  TFile *f_ZinclMu_RD = new TFile("../Zpt_RealData_12Bin/ZptXsecErrors/ZptXsecErrors_FidVolume.root");
  TH1D* hWpt_RD = (TH1D*)f_WinclMu_RD->Get("hData_Xsec_BornLogScaleNorm");
  TH1D* hZpt_RD = (TH1D*)f_ZinclMu_RD->Get("hZptDiffXsec12BinInFidNorm");
  TH1D* hWZratio_RD = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  ///Powheg
  TFile *f_WinclMu_Powheg = new TFile("../Wpt_PowhegPreFSR_12Bin/root/InclWToMuNu_PowhegWpt.root");
  TFile *f_ZinclMu_Powheg = new TFile("../Zpt_PowhegPreFSR_12Bin/root/ZToMuMu_Powheg.root");
  TH1D* hWpt_Powheg = (TH1D*)f_WinclMu_Powheg->Get("hxsec_NormDiff");
  TH1D* hZpt_Powheg = (TH1D*)f_ZinclMu_Powheg->Get("hxsec_NormDiff");
  TH1D* hWZratio_Powheg = new TH1D("W/Z ratio", "", nBins-1,WptBins);
 
  ///ResBos
  TFile *f_WinclMu_Resbos = new TFile("../WpT_Resbos_12Bin/root/WpTincl_Resbos.root");
  TFile *f_ZinclMu_Resbos = new TFile("../ZpT_Resbos_12Bin/root/ZToMuMu_Resbos.root");
  //TFile *f_ZinclMu_Resbos = new TFile("../ZpT_Resbos_12Bin/root_binMerge/ZToMuMu_Resbos.root");
  TH1D* hWpt_Resbos = (TH1D*)f_WinclMu_Resbos->Get("NormDiffXsec_Resbos");
  TH1D* hZpt_Resbos = (TH1D*)f_ZinclMu_Resbos->Get("NormDiffXsec");
  TH1D* hWZratio_Resbos = new TH1D("W/Z ratio", "", nBins-1,WptBins); 

  ///FEWZ
  TFile *f_ZinclMu_FEWZ = new TFile("../ZpT_FEWZ_12Bin/root/ZToMuMu_FEWZ.root");
  TFile *f_WinclMu_FEWZ = new TFile("../WpT_FEWZ_12Bin/root/WinclToMuNu_FEWZ.root");
  TH1D* hZpt_FEWZ = (TH1D*)f_ZinclMu_FEWZ->Get("hxsec_NormDiff");
  TH1D* hWpt_FEWZ = (TH1D*)f_WinclMu_FEWZ->Get("hxsec_NormDiff");
  TH1D* hZpt_Up_FEWZ = (TH1D*)f_ZinclMu_FEWZ->Get("hxsec_NormDiff_up");
  TH1D* hZpt_Down_FEWZ = (TH1D*)f_ZinclMu_FEWZ->Get("hxsec_NormDiff_down");
  TH1D* hWpt_Up_FEWZ = (TH1D*)f_WinclMu_FEWZ->Get("hxsec_NormDiff_up");
  TH1D* hWpt_Down_FEWZ = (TH1D*)f_WinclMu_FEWZ->Get("hxsec_NormDiff_down");
  TH1D* hWZratio_FEWZ = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  double Wpt_RD[12]={0};
  double Zpt_RD[12]={0};
  double WptErr_RD[12]={0};
  double ZptErr_RD[12]={0};
 
  double Wpt_Powheg[12]={0};
  double Zpt_Powheg[12]={0};
  double WptStatErr_Powheg[12]={0};
  double ZptStatErr_Powheg[12]={0};
 
  double Wpt_Resbos[12]={0};
  double Zpt_Resbos[12]={0};
  double WptErr_Resbos[12]={0};
  double ZptErr_Resbos[12]={0};
 
  double Wpt_FEWZ[12]={0};
  double Zpt_FEWZ[12]={0};
  double WptStatErr_FEWZ[12]={0};
  double ZptStatErr_FEWZ[12]={0};
 
  double WZratio_RD[12]={0};
  double WZratioErr_RD[12]={0};

  double WZratio_Powheg[12]={0};
  double WZratioStatErr_Powheg[12]={0};
  double WZratioPDFErr_Powheg[12]={0};
  double WZratioTotalErr_Powheg[12]={0}; // sqrt(Stat^2 + PDF^2 + Scale^2)
  
  double WZratio_Resbos[12]={0};
  double WZratioErr_Resbos[12]={0};
  
  double WZratio_FEWZ[12]={0};
  double WZratio_Up_FEWZ[12]={0};
  double WZratio_Down_FEWZ[12]={0};
  double WZratioStatErr_FEWZ[12]={0}; // Error propagation
  double WZratioPDFErr_FEWZ[12]={0}; // Norm PDF %
  double WZratioScaleErr_FEWZ[12]={0}; // Norm Scale subtraction
  double WZratioTotalErr_FEWZ[12]={0}; // sqrt(Stat^2 + PDF^2 + Scale^2)
  
  cout << fixed << setprecision(8) << endl;
  cout << " ===== Wpt and Zpt Normalilzed Differntial cross-section and errors In Fiducial Volume ==="<< endl;
  cout << "Bin"<<
    "\t Wpt  " << " \t\tError " << 
    "\t\t\t Zpt  " << " \t\tError " << endl;
  for(int i=0;i<12;i++)
  {
    Wpt_RD[i] = hWpt_RD->GetBinContent(i+1); 
    Zpt_RD[i] = hZpt_RD->GetBinContent(i+1); 
    WptErr_RD[i] = hWpt_RD->GetBinError(i+1); 
    ZptErr_RD[i] = hZpt_RD->GetBinError(i+1); 
    
//    cout <<i<< "\t" <<Wpt_RD[i] << "\t" <<  WptErr_RD[i] << "\t\t" << Zpt_RD[i] << "\t" << ZptErr_RD[i] << endl;
  
    // Data ratio and error propagation
    //WZratio_RD[i] = Wpt_RD[i] / Zpt_RD[i] ; 
    WZratio_RD[i] = Zpt_RD[i] / Wpt_RD[i] ; 
    WZratioErr_RD[i] = WZratio_RD[i] * sqrt(WptErr_RD[i]*WptErr_RD[i]/Wpt_RD[i]/Wpt_RD[i] + ZptErr_RD[i]*ZptErr_RD[i]/Zpt_RD[i]/Zpt_RD[i]); 
   
    ///Powheg,FEWZ,Resbos ratio
    WZratio_Powheg[i] = hZpt_Powheg->GetBinContent(i+1) / hWpt_Powheg->GetBinContent(i+1) ; 
    WZratio_Resbos[i] = hZpt_Resbos->GetBinContent(i+1) / hWpt_Resbos->GetBinContent(i+1) ; 
    WZratio_FEWZ[i]   = hZpt_FEWZ->GetBinContent(i+1)  / hWpt_FEWZ->GetBinContent(i+1) ; 
    WZratio_Up_FEWZ[i]   = hZpt_Up_FEWZ->GetBinContent(i+1)  / hWpt_Up_FEWZ->GetBinContent(i+1) ; 
    WZratio_Down_FEWZ[i]   = hZpt_Down_FEWZ->GetBinContent(i+1)  / hWpt_Down_FEWZ->GetBinContent(i+1) ; 

    WZratioScaleErr_FEWZ[i] = TMath::Max(fabs(WZratio_Up_FEWZ[i]-WZratio_FEWZ[i]),fabs(WZratio_Down_FEWZ[i]-WZratio_FEWZ[i]));

    // Powheg ratio error propagation
    Wpt_Powheg[i] = hWpt_Powheg->GetBinContent(i+1);
    WptStatErr_Powheg[i] = hWpt_Powheg->GetBinError(i+1);

    Zpt_Powheg[i] = hZpt_Powheg->GetBinContent(i+1);
    ZptStatErr_Powheg[i] = hZpt_Powheg->GetBinError(i+1);

    WZratioStatErr_Powheg[i] = WZratio_Powheg[i] * TMath::Sqrt((ZptStatErr_Powheg[i]*ZptStatErr_Powheg[i]/Zpt_Powheg[i]/Zpt_Powheg[i] + WptStatErr_Powheg[i]*WptStatErr_Powheg[i]/Wpt_Powheg[i]/Wpt_Powheg[i]));

    Powheg_PDFUncer(WZratio_Powheg,WZratioPDFErr_Powheg);
    
    WZratioTotalErr_Powheg[i] = sqrt(WZratioStatErr_Powheg[i]*WZratioStatErr_Powheg[i] + WZratioPDFErr_Powheg[i]*WZratioPDFErr_Powheg[i]); //  
    cout << Form("ZW Powheg Ratio : %.4f \t Stat : %.4f \t PDF : %.4f ",WZratio_Powheg[i],WZratioStatErr_Powheg[i],WZratioPDFErr_Powheg[i]) << endl;
    
    //*
    // Resbos ratio error propagation
    Wpt_Resbos[i] = hWpt_Resbos->GetBinContent(i+1);
    WptErr_Resbos[i] = hWpt_Resbos->GetBinError(i+1);

    Zpt_Resbos[i] = hZpt_Resbos->GetBinContent(i+1);
    ZptErr_Resbos[i] = hZpt_Resbos->GetBinError(i+1);

    WZratioErr_Resbos[i] = WZratio_Resbos[i] * TMath::Sqrt((ZptErr_Resbos[i]*ZptErr_Resbos[i]/Zpt_Resbos[i]/Zpt_Resbos[i] + WptErr_Resbos[i]*WptErr_Resbos[i]/Wpt_Resbos[i]/Wpt_Resbos[i]));
//*/
    // FEWZ ratio error propagation
    Wpt_FEWZ[i] = hWpt_FEWZ->GetBinContent(i+1);
    WptStatErr_FEWZ[i] = hWpt_FEWZ->GetBinError(i+1);

    Zpt_FEWZ[i] = hZpt_FEWZ->GetBinContent(i+1);
    ZptStatErr_FEWZ[i] = hZpt_FEWZ->GetBinError(i+1);

    WZratioStatErr_FEWZ[i] = WZratio_FEWZ[i] * TMath::Sqrt((ZptStatErr_FEWZ[i]*ZptStatErr_FEWZ[i]/Zpt_FEWZ[i]/Zpt_FEWZ[i] + WptStatErr_FEWZ[i]*WptStatErr_FEWZ[i]/Wpt_FEWZ[i]/Wpt_FEWZ[i]));

    FEWZ_PDFUncer(WZratio_FEWZ,WZratioPDFErr_FEWZ);
    //cout << "WZratioPDFErr_FEWZ in main : " << WZratioPDFErr_FEWZ[i] << endl; // check PDF error number is correctly returned

    WZratioTotalErr_FEWZ[i] = sqrt(WZratioStatErr_FEWZ[i]*WZratioStatErr_FEWZ[i] + WZratioPDFErr_FEWZ[i]*WZratioPDFErr_FEWZ[i] + WZratioScaleErr_FEWZ[i]*WZratioScaleErr_FEWZ[i]); //  
    
    //cout << Form("Wpt_FEWZ : %.8f +- %.8f \t Zpt_FEWZ : %.8f +- %.8f",Wpt_FEWZ[i],WptStatErr_FEWZ[i],Zpt_FEWZ[i],ZptStatErr_FEWZ[i]) << endl;
  }


  //// Print Powheg errors
  cout <<"=====  Print Powheg errors "<< endl;
  cout << fixed << setprecision(2) << endl;
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioStatErr_Powheg: " << WZratioStatErr_Powheg[i] <<"\t in %\t"<< 100*WZratioStatErr_Powheg[i]/WZratio_Powheg[i]<<endl; 
  }
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioPDFErr_Powheg : " << WZratioPDFErr_Powheg[i] << "\t in %\t"<< 100*WZratioPDFErr_Powheg[i]/WZratio_Powheg[i]<<endl;
  }
 //// Print FEWZ errors
  cout <<"=====  Print FEWZ errors "<< endl;

  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioStatErr_FEWZ : " << WZratioStatErr_FEWZ[i] << "\t in %\t"<< 100*WZratioStatErr_FEWZ[i]/WZratio_FEWZ[i]<<endl;
  }

  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioPDFErr_FEWZ : " << WZratioPDFErr_FEWZ[i] << "\t in %\t"<< 100*WZratioPDFErr_FEWZ[i]/WZratio_FEWZ[i]<<endl;
  }
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioScaleErr_FEWZ : " << WZratioScaleErr_FEWZ[i] << "\t in %\t"<< 100*WZratioScaleErr_FEWZ[i]/WZratio_FEWZ[i]<<endl;
  }

  cout<<"Wpt FEWZ\t"<<"Zpt FEWZ"<<endl; 
  for(int i=0;i<12;i++)
  {
   cout<<hWpt_FEWZ->GetBinContent(i+1) <<"\t"<<hZpt_FEWZ->GetBinContent(i+1)<<endl ; 
  }
  
  cout << fixed << setprecision(3) << endl;
  cout << " ============ Wpt/Zpt ratio and errors(Number) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error " << endl;
  for(int i=0;i<12;i++)
  {

    cout <<i+1<< "\t" <<  WZratio_RD[i] <<  " \\pm " << WZratioErr_RD[i] <<  endl;
  } 
  
  cout << " ============ Wpt/Zpt ratio and errors(in %) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error % " << endl;
  for(int i=0;i<12;i++)
  {

    cout <<i+1<< "\t" <<  WZratio_RD[i] <<  " \\pm " << WZratioErr_RD[i]/WZratio_RD[i]*100 <<  endl;
  } 
  
  for(int i=0;i<12;i++)
  {

    hWZratio_RD->SetBinContent(i+1,WZratio_RD[i]);
    hWZratio_RD->SetBinError(i+1,WZratioErr_RD[i]);
    
    hWZratio_Powheg->SetBinContent(i+1,WZratio_Powheg[i]);
    hWZratio_Powheg->SetBinError(i+1,WZratioTotalErr_Powheg[i]);
  
    hWZratio_Resbos->SetBinContent(i+1,WZratio_Resbos[i]);
    hWZratio_Resbos->SetBinError(i+1,WZratioErr_Resbos[i]);
    //hWZratio_Resbos->SetBinError(i+1,0);
    
    hWZratio_FEWZ->SetBinContent(i+1,WZratio_FEWZ[i]);
    hWZratio_FEWZ->SetBinError(i+1,WZratioTotalErr_FEWZ[i]);
  }
  
// Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataTotalErr = new TH1D("hRatioDataTotalErr","hRatioDataTotalErr",nBins-1,WptBins);
  
  TH1D *hRatioResbosTotalErr = new TH1D("hRatioResbosTotalErr","hRatioResbosTotalErr",nBins-1,WptBins);
  
  TH1D *hRatioPowhegStatErr = new TH1D("hRatioPowhegStatErr","hRatioPowhegStatErr",nBins-1,WptBins);
  TH1D *hRatioPowhegPDFErr = new TH1D("hRatioPowhegPDFErr","hRatioPowhegPDFErr",nBins-1,WptBins);
  TH1D *hRatioPowhegTotalErr = new TH1D("hRatioPowhegTotalErr","hRatioPowhegTotalErr",nBins-1,WptBins);
  
  TH1D *hRatioFEWZStatErr = new TH1D("hRatioFEWZStatErr","hRatioFEWZStatErr",nBins-1,WptBins);
  TH1D *hRatioFEWZPDFErr = new TH1D("hRatioFEWZPDFErr","hRatioFEWZPDFErr",nBins-1,WptBins);
  TH1D *hRatioFEWZScaleErr = new TH1D("hRatioFEWZScaleErr","hRatioFEWZScaleErr",nBins-1,WptBins);
  TH1D *hRatioFEWZTotalErr = new TH1D("hRatioFEWZTotalErr","hRatioFEWZTotalErr",nBins-1,WptBins);

  for(int i(0); i<nBins-1; i++)
  {
    hRatioDataTotalErr->SetBinContent(i+1,1.);
    hRatioDataTotalErr->SetBinError(i+1,WZratioErr_RD[i] / WZratio_RD[i]);

    hRatioResbosTotalErr->SetBinContent(i+1,WZratio_Resbos[i] / WZratio_RD[i]);
    hRatioResbosTotalErr->SetBinError(i+1,WZratioErr_Resbos[i] / WZratio_RD[i]);
    //hRatioResbosTotalErr->SetBinError(i+1,0);
    
    hRatioPowhegStatErr->SetBinContent(i+1,WZratio_Powheg[i] / WZratio_RD[i]);
    hRatioPowhegStatErr->SetBinError(i+1,WZratioStatErr_Powheg[i] /WZratio_RD[i]);
    hRatioPowhegPDFErr->SetBinContent(i+1,WZratio_Powheg[i] / WZratio_RD[i]);
    hRatioPowhegPDFErr->SetBinError(i+1,(WZratioStatErr_Powheg[i]+WZratioPDFErr_Powheg[i]) /WZratio_RD[i]);
    hRatioPowhegTotalErr->SetBinContent(i+1,WZratio_Powheg[i] / WZratio_RD[i]);
    hRatioPowhegTotalErr->SetBinError(i+1,WZratioTotalErr_Powheg[i] /WZratio_RD[i]);

    hRatioFEWZStatErr->SetBinContent(i+1,WZratio_FEWZ[i] / WZratio_RD[i]);
    hRatioFEWZStatErr->SetBinError(i+1,WZratioStatErr_FEWZ[i] / WZratio_RD[i]);
    hRatioFEWZPDFErr->SetBinContent(i+1,WZratio_FEWZ[i] / WZratio_RD[i]);
    hRatioFEWZPDFErr->SetBinError(i+1,(WZratioStatErr_FEWZ[i]+WZratioPDFErr_FEWZ[i]) / WZratio_RD[i]);
    hRatioFEWZScaleErr->SetBinContent(i+1,WZratio_FEWZ[i] / WZratio_RD[i]);
    hRatioFEWZScaleErr->SetBinError(i+1,(WZratioStatErr_FEWZ[i]+WZratioPDFErr_FEWZ[i]+WZratioScaleErr_FEWZ[i]) / WZratio_RD[i]);
    hRatioFEWZTotalErr->SetBinContent(i+1,WZratio_FEWZ[i] / WZratio_RD[i]);
    hRatioFEWZTotalErr->SetBinError(i+1,WZratioTotalErr_FEWZ[i] / WZratio_RD[i]);
  }

  //Define TGraph
  TGraphErrors *tgData = new TGraphErrors(hWZratio_RD);
  TGraphErrors *tgDataStatSystRatioBand = new TGraphErrors(hRatioDataTotalErr);
  
  TGraphErrors *tgResbos = new TGraphErrors(hWZratio_Resbos);
  TGraphErrors *tgResbosScaleRatioBand = new TGraphErrors(hRatioResbosTotalErr);
  
  TGraphErrors *tgPowheg = new TGraphErrors(hWZratio_Powheg);
  TGraphErrors *tgPowhegStatRatioBand = new TGraphErrors(hRatioPowhegStatErr);
  TGraphErrors *tgPowhegPDFRatioBand = new TGraphErrors(hRatioPowhegPDFErr);

  TGraphErrors *tgFEWZ = new TGraphErrors(hWZratio_FEWZ);
  TGraphErrors *tgFEWZStatRatioBand = new TGraphErrors(hRatioFEWZStatErr);
  TGraphErrors *tgFEWZPDFRatioBand = new TGraphErrors(hRatioFEWZPDFErr);
  TGraphErrors *tgFEWZScaleRatioBand = new TGraphErrors(hRatioFEWZScaleErr);

  //// Now design and Draw 
  gStyle->SetLineWidth(2.);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesSpacing(0.75);
  gStyle->SetHatchesLineWidth(2);

  // Marker Color
  Color_t MarkerColor_ResBos = kBlue+2;
  Color_t MarkerColor_POWHEG = kRed+3;
  Color_t MarkerColor_FEWZ = kGreen+4;

  // Marker Style (circle:20  square:21  triangle:22)
  Style_t MarkerStyle_ResBos = 21;
  Style_t MarkerStyle_POWHEG = 21;
  Style_t MarkerStyle_FEWZ = 21;

  Color_t BandColor_ResBos_scale = kBlue;
  Color_t BandColor_ResBos_PDF = kMagenta-9;
  Color_t BandColor_POWHEG_stat = kRed+2;
  Color_t BandColor_POWHEG_PDF = kYellow;
  Color_t BandColor_FEWZ_stat = kGreen+3;
  Color_t BandColor_FEWZ_PDF = kGreen+1;
  Color_t BandColor_FEWZ_scale = kCyan-9;
  
  // Band Color Transparent
  TColor *colResbos = gROOT->GetColor(kBlue);				// Used ResBos distribution and scale ratio band
  TColor *colResbosScale = gROOT->GetColor(BandColor_ResBos_scale);	// Used ResBos distribution and scale ratio band
  TColor *colResbosPDF = gROOT->GetColor(BandColor_ResBos_PDF);		// Used ResBos PDF ratio band
  TColor *colPowheg = gROOT->GetColor(kRed);				// Used Powheg distribution
  TColor *colPowhegStat = gROOT->GetColor(BandColor_POWHEG_stat);	// used Powheg ratio band
  TColor *colPowhegPDF = gROOT->GetColor(BandColor_POWHEG_PDF);		// used Powheg ratio band
  TColor *colFEWZ = gROOT->GetColor(kGreen);				// used FEWZ distribution
  TColor *colFEWZStat = gROOT->GetColor(BandColor_FEWZ_stat);		// used FEWZ stat ratio band
  TColor *colFEWZPDF = gROOT->GetColor(BandColor_FEWZ_PDF);		// used FEWZ PDF ratio band
  TColor *colFEWZScale = gROOT->GetColor(BandColor_FEWZ_scale);		// used FEWZ scale ratio band
  colResbosScale->SetAlpha(0.4);
  colResbosPDF->SetAlpha(0.4);
  colPowheg->SetAlpha(0.3);
  colPowhegStat->SetAlpha(0.6);
  colPowhegPDF->SetAlpha(0.6);
  colFEWZ->SetAlpha(0.5); 
  colFEWZStat->SetAlpha(0.6); 
  colFEWZPDF->SetAlpha(0.7); 
  colFEWZScale->SetAlpha(0.7); 


  TLegend *lL =new TLegend(0.25,0.37,0.65,0.75); lL->SetFillColor(0); lL->SetBorderSize(0);
  lL->AddEntry(tgData,"Data","PLE1");
  lL->AddEntry(tgResbos,"ResBos CT10 NNLL","f");
  lL->AddEntry(tgPowheg,"POWHEG CT10 NLO","f");
  lL->AddEntry(tgFEWZ,"FEWZ CT10 NNLO","f");

  TPaveText *channel = new TPaveText(0.27,0.78,0.67,0.84,"NDC");
  channel->SetBorderSize(0);
  channel->SetFillStyle(0);
  channel->AddText("Z #rightarrow #mu^{+}#mu^{-} / W #rightarrow #mu#nu_{#mu}");

  // Canvas for distribution
  TCanvas *lC1 = new TCanvas("Can","Can",50,50,W,H);
  lC1->SetFillColor(0);
  lC1->SetBorderMode(0);
  lC1->SetFrameFillStyle(0);
  lC1->SetFrameBorderMode(0);
  //lC1->SetLeftMargin( L/W );
  //lC1->SetLeftMargin( 0.15 );
  lC1->SetLeftMargin( 0.20 );
  lC1->SetRightMargin( R/W );
  lC1->SetTopMargin( T/H );
  lC1->SetBottomMargin( B/H );
  lC1->SetTickx(0);
  lC1->SetTicky(0);
  lC1->SetLogx();
  //lC1->SetLogy();

  // Frame setting 
  tgPowheg->GetYaxis()->SetRangeUser(0.,6.0);
  //tgPowheg->SetMinimum(1e-7);
  //tgPowheg->SetMaximum(5*TMath::MaxElement(n,tgPowheg->GetY()));
  tgPowheg->SetTitle("");
  tgPowheg->GetYaxis()->SetTitle("#left(#frac{1}{#lower[0.30]{#sigma}} #frac{d#sigma}{#lower[0.10]{p_{T}^{Z}}}#right) #lower[0.5]{#scale[2]{/}}#left(#frac{1}{#lower[0.30]{#sigma}} #frac{d#sigma}{#lower[0.07]{p_{T}^{W}}}#right)");
  tgPowheg->GetYaxis()->SetTitleOffset(2.2);
  tgPowheg->GetYaxis()->SetTitleSize(0.04);
  tgPowheg->GetYaxis()->SetLabelSize(0.04);

  tgPowheg->GetXaxis()->SetRangeUser(0.,600);
  tgPowheg->GetXaxis()->SetTitle("p_{T}^{V} [GeV]");
  tgPowheg->GetXaxis()->SetTitleSize(0.04);
  tgPowheg->GetXaxis()->SetTitleOffset(0.55);
  tgPowheg->GetXaxis()->SetLabelSize(0.04);

  //Data

  //Powheg
  tgPowheg->SetFillColor(kRed);
  tgFEWZ->SetFillColor(kGreen);
  tgResbos->SetFillColor(kBlue);

  //Draw Original Diff-Xsec Distribution
  tgPowheg->Draw("A2");
  tgResbos->Draw("2 same");
  tgFEWZ->Draw("2 same");
  tgData->Draw("P same");

  lL->Draw();
  channel->Draw();
  CMS_lumi(lC1,iPeriod,iPos);
  lC1->Update();
  lC1->RedrawAxis();
  lC1->GetFrame()->Draw();

  sprintf(tmpName,"RatioNormZW_Fid.");
  lC1->SaveAs(tmpName+format);

  // Ratio plot style 
  tgDataStatSystRatioBand->SetMarkerStyle(20);
  tgDataStatSystRatioBand->SetMarkerColor(kBlack);
  tgDataStatSystRatioBand->SetMarkerSize(0.7);
  tgDataStatSystRatioBand->SetLineWidth(2.0);
  tgDataStatSystRatioBand->SetLineColor(kBlack);

  tgDataStatSystRatioBand->SetFillStyle(3354);
  tgDataStatSystRatioBand->SetFillColor(kGray+1);

  // Resbos Ratio plot style
  //tgResbosPDFRatioBand->SetMarkerColor(MarkerColor_ResBos);
  //tgResbosPDFRatioBand->SetMarkerStyle(MarkerStyle_ResBos);
  //tgResbosPDFRatioBand->SetFillColor(BandColor_ResBos_PDF);

  tgResbosScaleRatioBand->SetMarkerColor(MarkerColor_ResBos);
  tgResbosScaleRatioBand->SetMarkerStyle(MarkerStyle_ResBos);
  tgResbosScaleRatioBand->SetFillColor(BandColor_ResBos_scale);

  // Powheg Ratio plot style
  tgPowhegStatRatioBand->SetMarkerColor(MarkerColor_POWHEG);
  tgPowhegStatRatioBand->SetMarkerStyle(MarkerStyle_POWHEG);
  tgPowhegStatRatioBand->SetFillColor(BandColor_POWHEG_stat);

  tgPowhegPDFRatioBand->SetMarkerColor(MarkerColor_POWHEG);
  tgPowhegPDFRatioBand->SetMarkerStyle(MarkerStyle_POWHEG);
  tgPowhegPDFRatioBand->SetFillColor(BandColor_POWHEG_PDF);

  // FEWZ Ratio plot style
  tgFEWZStatRatioBand->SetMarkerColor(MarkerColor_FEWZ);
  tgFEWZStatRatioBand->SetMarkerStyle(MarkerStyle_FEWZ);
  tgFEWZStatRatioBand->SetFillColor(BandColor_FEWZ_stat);

  tgFEWZPDFRatioBand->SetMarkerColor(MarkerColor_FEWZ);
  tgFEWZPDFRatioBand->SetMarkerStyle(MarkerStyle_FEWZ);
  tgFEWZPDFRatioBand->SetFillColor(BandColor_FEWZ_PDF);

  tgFEWZScaleRatioBand->SetMarkerColor(MarkerColor_FEWZ);
  tgFEWZScaleRatioBand->SetMarkerStyle(MarkerStyle_FEWZ);
  tgFEWZScaleRatioBand->SetFillColor(BandColor_FEWZ_scale);

  // Canvas for Theory/data Ratio
  TCanvas *lC2 = new TCanvas("Can","Can",50,50,W,H); 
  lC2->SetFillColor(0);
  lC2->SetBorderMode(0);
  lC2->SetFrameFillStyle(0);
  lC2->SetFrameBorderMode(0);
  lC2->SetLeftMargin( L/W );
  lC2->SetRightMargin( R/W );
  lC2->SetTopMargin( T/H );
  lC2->SetBottomMargin( B/H );
  CMS_lumi(lC2,iPeriod,iPos);

  lC2->Divide(1,3,0,0);
  lC2->cd(1)->SetPad(0,0.66,0.96,0.945);
  //lC2->cd(1)->SetFillColor(2);
  lC2->cd(1)->SetTickx(1);
  lC2->cd(1)->SetTicky(1);
  lC2->cd(1)->SetLogx(1);
  
  TLegend *lResbos =new TLegend(0.18,0.15,0.67,0.30); lResbos->SetFillColor(0); lResbos->SetBorderSize(0);
  lResbos-> SetNColumns(2);
  lResbos->AddEntry(tgResbosScaleRatioBand,"ResBos uncertainty","FP");
  lResbos->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  //lResbos->AddEntry(tgResbosPDFRatioBand,"ResBos PDF","F");
  lResbos->SetTextSize(0.08);

  TPaveText *tResBos = new TPaveText(0.18,0.82,0.34,0.92,"NDC");
  tResBos->SetBorderSize(0);
  tResBos->SetFillStyle(0);
  tResBos->SetTextSize(0.12);
  tResBos->AddText("#font[42]{ResBos}");
 
  TPaveText *tChannel = new TPaveText(0.40,0.82,0.72,0.93,"NDC");
  tChannel->SetBorderSize(0);
  tChannel->SetFillStyle(0);
  tChannel->SetTextSize(0.12);
  tChannel->AddText("Z #rightarrow #mu^{+}#mu^{-} / W #rightarrow #mu#nu_{#mu}");

  tgDataStatSystRatioBand->GetYaxis()->SetRangeUser(-0.2,2.6);
  tgDataStatSystRatioBand->GetYaxis()->SetTitle("Theory/Data");
  tgDataStatSystRatioBand->GetYaxis()->CenterTitle();
  tgDataStatSystRatioBand->GetYaxis()->SetTitleOffset(1.3);
  tgDataStatSystRatioBand->GetYaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetYaxis()->SetTitleSize(33);
  tgDataStatSystRatioBand->GetYaxis()->SetLabelFont(43);
  tgDataStatSystRatioBand->GetYaxis()->SetLabelSize(29);
  tgDataStatSystRatioBand->GetYaxis()->SetNdivisions(405);
  tgDataStatSystRatioBand->GetXaxis()->SetRangeUser(0,600);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleOffset(0.6);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleSize(20);
  tgDataStatSystRatioBand->GetXaxis()->SetLabelSize(0.11);
  tgDataStatSystRatioBand->Draw("2 A");
  //tgResbosPDFRatioBand->Draw("2 P");
  tgResbosScaleRatioBand->Draw("2 P");
  lResbos->Draw();
  tResBos->Draw();
  tChannel->Draw();

  //Powheg Ratio plot
  lC2->cd(2)->SetPad(0,0.39,0.96,0.64);
  //lC2->cd(2)->SetFillColor(4);
  lC2->cd(2)->SetTickx(1);
  lC2->cd(2)->SetTicky(1);
  lC2->cd(2)->SetLogx(1);

  TLegend *lPowheg =new TLegend(0.18,0.04,0.60,0.34); lPowheg->SetFillColor(0); lPowheg->SetBorderSize(0);
  lPowheg-> SetNColumns(2);
  lPowheg->AddEntry(tgPowhegStatRatioBand,"POWHEG stat","FP");
  lPowheg->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  lPowheg->AddEntry(tgPowhegPDFRatioBand,"POWHEG PDF    ","F");
  //lPowheg->AddEntry(hRatioDataStatErr,"","");
  lPowheg->SetTextSize(0.09);

  TPaveText *tPowheg = new TPaveText(0.18,0.82,0.34,0.92,"NDC");
  tPowheg->SetBorderSize(0);
  tPowheg->SetFillStyle(0);
  tPowheg->SetTextSize(0.12);
  tPowheg->AddText("#font[42]{POWHEG}");
  
  tgDataStatSystRatioBand->Draw("2 A");
  tgPowhegPDFRatioBand->Draw("2 P");
  tgPowhegStatRatioBand->Draw("2 P");
  lPowheg->Draw();
  tPowheg->Draw();

  // FEWZ Ratio Plot
  lC2->cd(3)->SetPad(0,0,0.96,0.37);
  lC2->cd(3)->SetBottomMargin(0.26);
  //lC2->cd(3)->SetFillColor(3);
  lC2->cd(3)->SetTickx(1);
  lC2->cd(3)->SetTicky(1);
  lC2->cd(3)->SetLogx(1);

  TLegend *lFewz =new TLegend(0.18,0.28,0.60,0.52); lFewz->SetFillColor(0); lFewz->SetBorderSize(0);
  lFewz-> SetNColumns(2);
  lFewz->AddEntry(tgFEWZStatRatioBand,"FEWZ stat","FP");
  lFewz->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  lFewz->AddEntry(tgFEWZPDFRatioBand,"FEWZ PDF","F");
  lFewz->AddEntry(tgFEWZScaleRatioBand,"FEWZ scales","F");
  lFewz->SetTextSize(0.06);

  TPaveText *tFewz = new TPaveText(0.18,0.82,0.34,0.92,"NDC");
  tFewz->SetBorderSize(0);
  tFewz->SetFillStyle(0);
  tFewz->SetTextSize(0.09);
  tFewz->AddText("#font[42]{FEWZ}");
 
  tgDataStatSystRatioBand->GetXaxis()->SetTitle("p_{T}^{V} [GeV]");
  tgDataStatSystRatioBand->GetXaxis()->SetTitleOffset(1.5);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleSize(33);
  tgDataStatSystRatioBand->Draw("2 A");
  tgFEWZScaleRatioBand->Draw("2");
  tgFEWZPDFRatioBand->Draw("2");
  tgFEWZStatRatioBand->Draw("2P");
  lFewz->Draw();
  tFewz->Draw();

  sprintf(tmpName,"ZWTheoryData.");
  lC2->SaveAs(tmpName+format);

  return 0;
}

void FEWZ_PDFUncer(double *ZWRatio, double *Err)
{
  // FEWZ Normalized PDF error in % unit
  double PDFErr[12] = {0.,};

  PDFErr[0] = 2.51;  
  PDFErr[1] = 0.61;
  PDFErr[2] = 0.96;
  PDFErr[3] = 1.13;
  PDFErr[4] = 1.44;
  PDFErr[5] = 1.61;
  PDFErr[6] = 1.82;
  PDFErr[7] = 2.07;
  PDFErr[8] = 2.42;
  PDFErr[9] = 2.34;
  PDFErr[10] =2.30;
  PDFErr[11] =2.47;

  for(int i(0); i<12; i++)
  {
    Err[i] = ZWRatio[i] * PDFErr[i] * 0.01;
    //cout << "PDFErr in function : " << Err[i] << endl;
  }
}
void Powheg_PDFUncer(double *ZWRatio, double *Err)
{
  // Powheg Normalized PDF error in % unit
  double PDFErr[12] = {0.,};

  PDFErr[0] = 0.70904 ;  
  PDFErr[1] = 0.08829 ;
  PDFErr[2] = 0.26802 ;
  PDFErr[3] = 0.48621 ;
  PDFErr[4] = 0.74545 ;
  PDFErr[5] = 0.78758 ;
  PDFErr[6] = 0.89615 ;
  PDFErr[7] = 1.05564 ;
  PDFErr[8] = 1.23073 ;
  PDFErr[9] = 1.09974 ;
  PDFErr[10]= 1.36105 ;
  PDFErr[11]= 2.64723 ;

  for(int i(0); i<12; i++)
  {
    Err[i] = ZWRatio[i] * PDFErr[i] * 0.01;
    //cout << "PDFErr in function : " << Err[i] << endl;
  }
}
