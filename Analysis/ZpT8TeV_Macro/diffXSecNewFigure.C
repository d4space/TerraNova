#include <iostream>
#include <iomanip>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "../Utils/const.h"
#include <TMath.h>

//const TString format("png");
const TString format("pdf");

const int n = 18;
double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};

int diffXSecNewFigure()
{
  gROOT->LoadMacro("../Utils/tdrstyle.C");
  setTDRStyle();
  gROOT->LoadMacro("../Utils/CMS_lumi.C");
  writeExtraText = "true";
  //extraText = "Preliminary";
  extraText = "";
  lumi_8TeV = "18.4 pb^{-1}";

  int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
  int iPos = 0;

  int H = 800;
  int W = 800;

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
  
  char tmpName[30],tmpName_org[30];

//  TFile *f_Data;
//  f_Data = new TFile("Zptout.root");
//  
//  TH1D* hData = (TH1D*)f_Data->Get("hdata")->Clone("hData");
//  double Data[18];
//  double DataStatSyst[18];
//  for(int i=0;i<18;i++)
//  {
//    Data[i] = hData->GetBinContent(i+1);
//    DataStatSyst[i] = hData->GetBinError(i+1);
//  }

  double Data[18];
  Data[0]	=	0.0333549   ;
  Data[1]	=       0.0553782   ; 	
  Data[2]	=       0.051878    ; 	
  Data[3]	=       0.0385965   ; 	
  Data[4]	=       0.0354796   ; 	
  Data[5]	=       0.0240797   ; 	
  Data[6]	=       0.0225014   ; 	
  Data[7]	=       0.0171383   ; 	
  Data[8]	=       0.0117612   ; 	
  Data[9]	=       0.00650165  ; 	
  Data[10]	=       0.00401844  ; 	
  Data[11]	=       0.00216843  ; 	
  Data[12]	=       0.000886498 ; 	
  Data[13]	=       0.000409235 ; 	
  Data[14]	=       0.000164556 ; 	
  Data[15]	=       8.04607e-05 ; 	
  Data[16]	=       8.91327e-06 ; 	
  Data[17]	=       4.43314e-06 ;

  double DataStatSyst[18];
  DataStatSyst[0]	=	0.00205861  ;
  DataStatSyst[1]	=       0.00265076  ; 	
  DataStatSyst[2]	=       0.00249936  ; 	
  DataStatSyst[3]	=       0.00229056  ; 	
  DataStatSyst[4]	=       0.00215988  ; 	
  DataStatSyst[5]	=       0.00189679  ; 	
  DataStatSyst[6]	=       0.00173219  ; 	
  DataStatSyst[7]	=       0.00154885  ; 	
  DataStatSyst[8]	=       0.000482165 ; 	
  DataStatSyst[9]	=       0.000359578 ; 	
  DataStatSyst[10]	=       0.000289945 ; 	
  DataStatSyst[11]	=       0.000144729 ; 	
  DataStatSyst[12]	=       9.5768e-05  ; 	
  DataStatSyst[13]	=       6.58497e-05 ; 	
  DataStatSyst[14]	=       2.83107e-05 ; 	
  DataStatSyst[15]	=       1.93958e-05 ; 	
  DataStatSyst[16]	=       6.26627e-06 ; 	
  DataStatSyst[17]	=       1.95776e-06 ;


  double Resbos[18];
  Resbos[0]	=	0.03335570;
  Resbos[1]	=        0.05909690; 	
  Resbos[2]	=        0.05080043; 	
  Resbos[3]	=        0.04145975; 	
  Resbos[4]	=        0.03394089; 	
  Resbos[5]	=        0.02773926; 	
  Resbos[6]	=        0.02220695; 	
  Resbos[7]	=        0.01779636; 	
  Resbos[8]	=        0.01131978; 	
  Resbos[9]	=        0.00601003; 	
  Resbos[10]	=        0.00357902; 	
  Resbos[11]	=        0.00188844; 	
  Resbos[12]	=        0.00085224; 	
  Resbos[13]	=        0.00041736; 	
  Resbos[14]	=        0.00018092; 	
  Resbos[15]	=        0.00006654; 	
  Resbos[16]	=        0.00002372; 	
  Resbos[17]	=        0.00000218;

  double ResbosSyst[18];
  ResbosSyst[0]		=	0.001694 ;
  ResbosSyst[1]		=       0.002232 ; 	
  ResbosSyst[2]		=       0.001961 ; 	
  ResbosSyst[3]		=       0.001568 ; 	
  ResbosSyst[4]		=       0.001186 ; 	
  ResbosSyst[5]		=       0.0008501; 	
  ResbosSyst[6]		=       0.0003468; 	
  ResbosSyst[7]		=       0.0006431; 	
  ResbosSyst[8]		=       0.0008231; 	
  ResbosSyst[9]		=       0.0005566; 	
  ResbosSyst[10]	=       0.0002914; 	
  ResbosSyst[11]	=       0.0001034; 	
  ResbosSyst[12]	=       5.643e-05; 	
  ResbosSyst[13]	=       3.159e-05; 	
  ResbosSyst[14]	=       1.899e-05; 	
  ResbosSyst[15]	=       7.658e-06; 	
  ResbosSyst[16]	=       2.876e-06; 	
  ResbosSyst[17]	=       2.919e-07;

  double Powheg[18];
  Powheg[0]	=	 0.0232823816;
  Powheg[1]	=        0.0534414648; 	
  Powheg[2]	=        0.0556310835; 	
  Powheg[3]	=        0.0454828550; 	
  Powheg[4]	=        0.0356700536; 	
  Powheg[5]	=        0.0280575476; 	
  Powheg[6]	=        0.0226186451; 	
  Powheg[7]	=        0.0182287311; 	
  Powheg[8]	=        0.0118469402; 	
  Powheg[9]	=        0.0064591971; 	
  Powheg[10]	=        0.0037984679; 	
  Powheg[11]	=        0.0019372698; 	
  Powheg[12]	=        0.0008350090; 	
  Powheg[13]	=        0.0003902566; 	
  Powheg[14]	=        0.0001503083; 	
  Powheg[15]	=        0.0000509369; 	
  Powheg[16]	=        0.0000184618; 	
  Powheg[17]	=        0.0000014680;

  double PowhegStat[18];
  PowhegStat[0]		=	 0.0001395273;
  PowhegStat[1]		=        0.0002032366; 	
  PowhegStat[2]		=        0.0002069700; 	
  PowhegStat[3]		=        0.0001884388; 	
  PowhegStat[4]		=        0.0001705318; 	
  PowhegStat[5]		=        0.0001524742; 	
  PowhegStat[6]		=        0.0001370749; 	
  PowhegStat[7]		=        0.0001241799; 	
  PowhegStat[8]		=        0.0000484067; 	
  PowhegStat[9]		=        0.0000366260; 	
  PowhegStat[10]	=        0.0000284210; 	
  PowhegStat[11]	=        0.0000144427; 	
  PowhegStat[12]	=        0.0000095287; 	
  PowhegStat[13]	=        0.0000065614; 	
  PowhegStat[14]	=        0.0000028981; 	
  PowhegStat[15]	=        0.0000016920; 	
  PowhegStat[16]	=        0.0000008289; 	
  PowhegStat[17]	=        0.0000000970;

  double PowhegPDF[18];
  PowhegPDF[0]		=	 0.412089;
  PowhegPDF[1]		=        0.31525 ; 	
  PowhegPDF[2]		=        0.209474; 	
  PowhegPDF[3]		=        0.188146; 	
  PowhegPDF[4]		=        0.250735; 	
  PowhegPDF[5]		=        0.21012 ; 	
  PowhegPDF[6]		=        0.351397; 	
  PowhegPDF[7]		=        0.284557; 	
  PowhegPDF[8]		=        0.304525; 	
  PowhegPDF[9]		=        0.315051; 	
  PowhegPDF[10]		=        0.356816; 	
  PowhegPDF[11]		=        0.619421; 	
  PowhegPDF[12]		=        1.16534 ; 	
  PowhegPDF[13]		=        1.56789 ; 	
  PowhegPDF[14]		=        2.29144 ; 	
  PowhegPDF[15]		=        2.92386 ; 	
  PowhegPDF[16]		=        4.04627 ; 	
  PowhegPDF[17]		=        6.01494 ;

  // Convert PDF% to Num, calculate total uncer.
  double PowhegPDFNum[18];
  double PowhegTotalUnc[18];
  for(int i=0;i<18;i++)
  {
    PowhegPDFNum[i] = 0.01*PowhegPDF[i]*Powheg[i];
    PowhegTotalUnc[i] = sqrt(PowhegStat[i]*PowhegStat[i] + PowhegPDFNum[i]*PowhegPDFNum[i]);
  }

  double FEWZ[18];
  FEWZ[0]	=	 100;//-0.08005877;
  FEWZ[1]	=        0.12748182 ; 	
  FEWZ[2]	=        0.07877656 ; 	
  FEWZ[3]	=        0.05332600 ; 	
  FEWZ[4]	=        0.03845661 ; 	
  FEWZ[5]	=        0.02904814 ; 	
  FEWZ[6]	=        0.02265475 ; 	
  FEWZ[7]	=        0.01805065 ; 	
  FEWZ[8]	=        0.01135670 ; 	
  FEWZ[9]	=        0.00594296 ; 	
  FEWZ[10]	=        0.00353301 ; 	
  FEWZ[11]	=        0.00185758 ; 	
  FEWZ[12]	=        0.00082569 ; 	
  FEWZ[13]	=        0.00039315 ; 	
  FEWZ[14]	=        0.00016611 ; 	
  FEWZ[15]	=        0.00005898 ; 	
  FEWZ[16]	=        0.00001959 ; 	
  FEWZ[17]	=        0.00000179 ;

  double FEWZStat[18];
  FEWZStat[0]	=	 0.0 	    ; //0.00094095;
  FEWZStat[1]	=        0.00039520 ; 	
  FEWZStat[2]	=        0.00024179 ; 	
  FEWZStat[3]	=        0.00016486 ; 	
  FEWZStat[4]	=        0.00012112 ; 	
  FEWZStat[5]	=        0.00009292 ; 	
  FEWZStat[6]	=        0.00007464 ; 	
  FEWZStat[7]	=        0.00006809 ; 	
  FEWZStat[8]	=        0.00003028 ; 	
  FEWZStat[9]	=        0.00001669 ; 	
  FEWZStat[10]	=        0.00001118 ; 	
  FEWZStat[11]	=        0.00000535 ; 	
  FEWZStat[12]	=        0.00000286 ; 	
  FEWZStat[13]	=        0.00000170 ; 	
  FEWZStat[14]	=        0.00000067 ; 	
  FEWZStat[15]	=        0.00000082 ; 	
  FEWZStat[16]	=        0.00000031 ; 	
  FEWZStat[17]	=        0.00000004 ;

  double FEWZPDFlow[18];
  FEWZPDFlow[0]		=	 0.0     ; //-6.7678;
  FEWZPDFlow[1]		=        1.0464  ; 	
  FEWZPDFlow[2]		=        0.7600  ; 	
  FEWZPDFlow[3]		=        1.0219  ; 	
  FEWZPDFlow[4]		=        1.2630  ; 	
  FEWZPDFlow[5]		=        1.4623  ; 	
  FEWZPDFlow[6]		=        1.6461  ; 	
  FEWZPDFlow[7]		=        1.8262  ; 	
  FEWZPDFlow[8]		=        2.1995  ; 	
  FEWZPDFlow[9]		=        2.7651  ; 	
  FEWZPDFlow[10]	=        3.2279  ; 	
  FEWZPDFlow[11]	=        3.8587  ; 	
  FEWZPDFlow[12]	=        4.6700  ; 	
  FEWZPDFlow[13]	=        5.4480  ; 	
  FEWZPDFlow[14]	=        6.4558  ; 	
  FEWZPDFlow[15]	=        7.8312  ; 	
  FEWZPDFlow[16]	=        9.3764  ; 	
  FEWZPDFlow[17]	=        11.4845 ;

  double FEWZPDFhigh[18];
  FEWZPDFhigh[0]	=	 0.0    ; //-6.4288;
  FEWZPDFhigh[1]	=        1.8289 ; 	
  FEWZPDFhigh[2]	=        1.3007 ; 	
  FEWZPDFhigh[3]	=        1.2656 ; 	
  FEWZPDFhigh[4]	=        1.3277 ; 	
  FEWZPDFhigh[5]	=        1.4084 ; 	
  FEWZPDFhigh[6]	=        1.5036 ; 	
  FEWZPDFhigh[7]	=        1.6060 ; 	
  FEWZPDFhigh[8]	=        1.8347 ; 	
  FEWZPDFhigh[9]	=        2.2177 ; 	
  FEWZPDFhigh[10]	=        2.5401 ; 	
  FEWZPDFhigh[11]	=        2.9723 ; 	
  FEWZPDFhigh[12]	=        3.4925 ; 	
  FEWZPDFhigh[13]	=        3.9540 ; 	
  FEWZPDFhigh[14]	=        4.4901 ; 	
  FEWZPDFhigh[15]	=        5.1299 ; 	
  FEWZPDFhigh[16]	=        5.7526 ; 	
  FEWZPDFhigh[17]	=        6.6129 ;

  double FEWZScale[18];
  FEWZScale[0]	=	 0.0	     ; //0.00842947;
  FEWZScale[1]	=        0.0106515   ; 	
  FEWZScale[2]	=        0.00266521  ; 	
  FEWZScale[3]	=        0.000428617 ; 	
  FEWZScale[4]	=        0.000114661 ; 	
  FEWZScale[5]	=        0.000205532 ; 	
  FEWZScale[6]	=        0.000345275 ; 	
  FEWZScale[7]	=        0.000297941 ; 	
  FEWZScale[8]	=        0.000323374 ; 	
  FEWZScale[9]	=        0.000210461 ; 	
  FEWZScale[10]	=        0.000147318 ; 	
  FEWZScale[11]	=        9.51447e-05 ; 	
  FEWZScale[12]	=        4.73258e-05 ; 	
  FEWZScale[13]	=        2.47657e-05 ; 	
  FEWZScale[14]	=        1.03102e-05 ; 	
  FEWZScale[15]	=        4.30353e-06 ; 	
  FEWZScale[16]	=        1.45553e-06 ; 	
  FEWZScale[17]	=        1.5117e-07  ;

  // Convert % to Num
  double FEWZPDFlowNum[18];
  double FEWZPDFhighNum[18];
  double FEWZTotalUnclow[18];
  double FEWZTotalUnchigh[18];
  for(int i=0;i<18;i++)
  {
    FEWZPDFlowNum[i] = 0.01*FEWZPDFlow[i]*FEWZ[i];
    FEWZPDFhighNum[i] = 0.01*FEWZPDFhigh[i]*FEWZ[i];
    FEWZTotalUnclow[i] = sqrt(FEWZStat[i]*FEWZStat[i] + FEWZPDFlowNum[i]*FEWZPDFlowNum[i] + FEWZScale[i]*FEWZScale[i]);
    FEWZTotalUnchigh[i] = sqrt(FEWZStat[i]*FEWZStat[i] + FEWZPDFhighNum[i]*FEWZPDFhighNum[i] + FEWZScale[i]*FEWZScale[i]);
  }

  for(int i=0;i<18;i++)
  {
    //cout << "hdata : " << hData->GetBinContent(i+1) << " +- " << hData->GetBinError(i+1) << endl;
    cout << "hdata : " << Data[i] << " +- " << DataStatSyst[i] << endl;
  }

  for(int i=0;i<18;i++)
  {
    cout << "Resbos : " << Resbos[i] << " +- " << ResbosSyst[i] << endl;
  }

  for(int i=0;i<18;i++)
  {
    cout << "Powheg : " << Powheg[i] << "\t stat : " << PowhegStat[i] << "\t pdf : " << PowhegPDFNum[i] << endl;
  }

  for(int i=0;i<18;i++)
  {
    cout << "FEWZ : " << FEWZ[i] << "\t stat : " << FEWZStat[i] << "\t pdflow : " << FEWZPDFlowNum[i] << "\t pdfhigh : " << FEWZPDFhighNum << "\t scale : " << FEWZScale[i] << endl;
  }

  // Binning
  double error_xlow[18];
  double error_xhigh[18];
  for (int i=0;i<n;i++) 
  {
    error_xlow[i] = (point[i] - bin[i]);
    error_xhigh[i] = (bin[i+1]-point[i]);
    cout << "xlow : " << error_xlow[i] << "\t xhigh : " << error_xhigh[i] << endl;
  }

  TH1D* hPowheg = new TH1D("","",18,1,18);
  for(int i=0;i<n;i++)
  {
    hPowheg->SetBinContent(i+1,Powheg[i]);
    hPowheg->SetBinError(i+1,PowhegTotalUnc[i]);
  }

  double DataRatioBand[18];
  double DataStatSystRatioBand[18];
  double ResbosRatioBand[18];
  double ResbosSystRatioBand[18];
  double PowhegRatioBand[18];
  double PowhegStatRatioBand[18];
  double PowhegPDFRatioBand[18];
  double FEWZRatioBand[18];
  double FEWZStatRatioBand[18];
  double FEWZPDFlowRatioBand[18];
  double FEWZPDFhighRatioBand[18];
  double FEWZScalelowRatioBand[18];
  double FEWZScalehighRatioBand[18];
  for(int i=0; i<18; i++)
  {
    DataRatioBand[i] = Data[i] / Data[i] ;
    DataStatSystRatioBand[i] = DataStatSyst[i] / Data[i] ;
    ResbosRatioBand[i] = Resbos[i] / Data[i] ;
    ResbosSystRatioBand[i] = ResbosSyst[i] / Data[i] ;
    PowhegRatioBand[i] = Powheg[i] / Data[i];
    PowhegStatRatioBand[i] = PowhegStat[i] / Data[i] ;
    PowhegPDFRatioBand[i] = (PowhegStat[i] + PowhegPDFNum[i]) / Data[i] ;
    FEWZRatioBand[i] = FEWZ[i] / Data[i] ;
    FEWZStatRatioBand[i] = FEWZStat[i] / Data[i] ;
    FEWZPDFlowRatioBand[i] = (FEWZStat[i] + FEWZPDFlowNum[i]) / Data[i];
    FEWZPDFhighRatioBand[i] = (FEWZStat[i] + FEWZPDFhighNum[i]) / Data[i];
    FEWZScalelowRatioBand[i] = (FEWZStat[i] + FEWZPDFlowNum[i] + FEWZScale[i]) / Data[i];
    FEWZScalehighRatioBand[i] = (FEWZStat[i] + FEWZPDFhighNum[i] + FEWZScale[i]) / Data[i];
  }

  //Define TGraph
  TGraphAsymmErrors *tgData = new TGraphAsymmErrors(n,point,Data,error_xlow,error_xhigh,DataStatSyst,DataStatSyst);
  TGraphAsymmErrors *tgDataStatSystRatioBand = new TGraphAsymmErrors(n,point,DataRatioBand,error_xlow,error_xhigh,DataStatSystRatioBand,DataStatSystRatioBand);
  
  TGraphAsymmErrors *tgResbos = new TGraphAsymmErrors(n,point,Resbos,error_xlow,error_xhigh,ResbosSyst,ResbosSyst);
  TGraphAsymmErrors *tgResbosSystRatioBand = new TGraphAsymmErrors(n,point,ResbosRatioBand,error_xlow,error_xhigh,ResbosSystRatioBand,ResbosSystRatioBand);
  
  TGraphAsymmErrors *tgPowheg = new TGraphAsymmErrors(n,point,Powheg,error_xlow,error_xhigh,PowhegTotalUnc,PowhegTotalUnc);
  TGraphAsymmErrors *tgPowhegStatRatioBand = new TGraphAsymmErrors(n,point,PowhegRatioBand,error_xlow,error_xhigh,PowhegStatRatioBand,PowhegStatRatioBand);
  TGraphAsymmErrors *tgPowhegPDFRatioBand = new TGraphAsymmErrors(n,point,PowhegRatioBand,error_xlow,error_xhigh,PowhegPDFRatioBand,PowhegPDFRatioBand);

  TGraphAsymmErrors *tgFEWZ = new TGraphAsymmErrors(n,point,FEWZ,error_xlow,error_xhigh,FEWZTotalUnclow,FEWZTotalUnchigh);
  TGraphAsymmErrors *tgFEWZStatRatioBand = new TGraphAsymmErrors(n,point,FEWZRatioBand,error_xlow,error_xhigh,FEWZStatRatioBand,FEWZStatRatioBand);
  TGraphAsymmErrors *tgFEWZPDFRatioBand = new TGraphAsymmErrors(n,point,FEWZRatioBand,error_xlow,error_xhigh,FEWZPDFlowRatioBand,FEWZPDFhighRatioBand);
  TGraphAsymmErrors *tgFEWZScaleRatioBand = new TGraphAsymmErrors(n,point,FEWZRatioBand,error_xlow,error_xhigh,FEWZScalelowRatioBand,FEWZScalehighRatioBand);

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

  // Band Color Transparent
  TColor *colBlue = gROOT->GetColor(kBlue);		// Used ResBos distribution and ratio band
  TColor *colRed = gROOT->GetColor(kRed);		// Used Powheg distribution
  TColor *colRedp2 = gROOT->GetColor(kRed+2);		// used Powheg ratio band
  TColor *colYellow = gROOT->GetColor(kYellow);		// used Powheg ratio band
  TColor *colGreen = gROOT->GetColor(kGreen);		// used FEWZ distribution
  TColor *colGreenp1 = gROOT->GetColor(kGreen+1);	// used FEWZ stat ratio band
  TColor *colGreenp3 = gROOT->GetColor(kGreen+3);	// used FEWZ PDF ratio band
  TColor *colCyanm9 = gROOT->GetColor(kCyan-9);		// used FEWZ scale ratio band
  colRedp2->SetAlpha(0.6);
  colRed->SetAlpha(0.2);
  colBlue->SetAlpha(0.3);
  colYellow->SetAlpha(0.6);
  colGreen->SetAlpha(0.5); 
  colGreenp1->SetAlpha(0.6); 
  colGreenp3->SetAlpha(0.7); 
  colCyanm9->SetAlpha(0.7); 

  Color_t BandColor_ResBos = kBlue;
  Color_t BandColor_POWHEG_stat = kRed+2;
  Color_t BandColor_POWHEG_PDF = kYellow;
  Color_t BandColor_FEWZ_stat = kGreen+3;
  Color_t BandColor_FEWZ_PDF = kGreen+1;
  Color_t BandColor_FEWZ_scale = kCyan-9;

  TLegend *lL =new TLegend(0.25,0.17,0.65,0.45); lL->SetFillColor(0); lL->SetBorderSize(0);
  lL->AddEntry(tgData,"Data","PLE1");
  lL->AddEntry(tgResbos,"ResBos CT10 NNLL","f");
  lL->AddEntry(tgPowheg,"POWHEG CT10 NLO","f");
  lL->AddEntry(tgFEWZ,"FEWZ CT10 NNLO","f");

  TPaveText *channel = new TPaveText(0.30,0.48,0.47,0.54,"NDC");
  channel->SetBorderSize(0);
  channel->SetFillStyle(0);
  channel->AddText("Z #rightarrow #mu^{+}#mu^{#font[122]{-}}");

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
  lC1->SetLogy();

  // Frame setting 
  //tgPowheg->GetYaxis()->SetRangeUser(1e-7,5*tgResbos[1]);
  tgPowheg->SetMinimum(1e-7);
  tgPowheg->SetMaximum(5*TMath::MaxElement(n,tgPowheg->GetY()));
  tgPowheg->SetTitle("");
  tgPowheg->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T}^{Z} [GeV]^{-1}");
  tgPowheg->GetYaxis()->SetTitleOffset(1.5);
  tgPowheg->GetYaxis()->SetTitleSize(0.05);
  tgPowheg->GetYaxis()->SetLabelSize(0.04);

  tgPowheg->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
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
  tgFEWZ->Draw("2 same");
  tgResbos->Draw("2 same");
  tgData->Draw("P same");

  lL->Draw();
  channel->Draw();
  CMS_lumi(lC1,iPeriod,iPos);
  lC1->Update();
  lC1->RedrawAxis();
  lC1->GetFrame()->Draw();

  sprintf(tmpName,"ZptNormDiffXsec.");
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
  tgResbosSystRatioBand->SetMarkerColor(MarkerColor_ResBos);
  tgResbosSystRatioBand->SetMarkerStyle(MarkerStyle_ResBos);
  tgResbosSystRatioBand->SetFillColor(BandColor_ResBos);

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
  
  TLegend *rL1 =new TLegend(0.18,0.15,0.60,0.30); rL1->SetFillColor(0); rL1->SetBorderSize(0);
  rL1-> SetNColumns(2);
  rL1->AddEntry(tgResbosSystRatioBand,"ResBos syst","FP");
  rL1->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  //rL1->AddEntry(tgDataStatSystRatioBand,"","");
  rL1->SetTextSize(0.08);

  TLegend *tL1 =new TLegend(0.18,0.82,0.34,0.92); tL1->SetFillColor(0); tL1->SetBorderSize(0);
  tL1->AddEntry(tgResbosSystRatioBand,"ResBos","");
  tL1->SetTextSize(0.12);
  tL1->SetTextFont(2);

  TPaveText *tb4 = new TPaveText(0.35,0.82,0.67,0.93,"NDC");
  tb4->SetBorderSize(0);
  tb4->SetFillStyle(0);
  tb4->SetTextSize(0.12);
  tb4->AddText("Z #rightarrow #mu^{+}#mu^{#font[122]{-}}");

  tgDataStatSystRatioBand->GetYaxis()->SetRangeUser(-0.8,3.2);
  tgDataStatSystRatioBand->GetYaxis()->SetTitle("Theory/Data");
  tgDataStatSystRatioBand->GetYaxis()->CenterTitle();
  tgDataStatSystRatioBand->GetYaxis()->SetTitleOffset(1.1);
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
  tgResbosSystRatioBand->Draw("2 P");
  rL1->Draw();
  tL1->Draw();
  tb4->Draw();

  //Powheg Ratio plot
  lC2->cd(2)->SetPad(0,0.39,0.96,0.64);
  //lC2->cd(2)->SetFillColor(4);
  lC2->cd(2)->SetTickx(1);
  lC2->cd(2)->SetTicky(1);
  lC2->cd(2)->SetLogx(1);

  TLegend *rL2 =new TLegend(0.18,0.04,0.60,0.34); rL2->SetFillColor(0); rL2->SetBorderSize(0);
  rL2-> SetNColumns(2);
  rL2->AddEntry(tgPowhegStatRatioBand,"POWHEG stat","FP");
  rL2->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  rL2->AddEntry(tgPowhegPDFRatioBand,"POWHEG PDF    ","FP");
  //rL2->AddEntry(hRatioDataStatErr,"","");
  rL2->SetTextSize(0.09);

  TLegend *tL2 =new TLegend(0.18,0.82,0.34,0.92); tL2->SetFillColor(0); tL2->SetBorderSize(0);
  tL2->AddEntry(tgPowhegStatRatioBand,"POWHEG","");
  tL2->SetTextSize(0.12);
  tL2->SetTextFont(2);

  tgDataStatSystRatioBand->Draw("2 A");
  tgPowhegPDFRatioBand->Draw("2 P");
  tgPowhegStatRatioBand->Draw("2 P");
  rL2->Draw();
  tL2->Draw();

  // FEWZ Ratio Plot
  lC2->cd(3)->SetPad(0,0,0.96,0.37);
  lC2->cd(3)->SetBottomMargin(0.26);
  //lC2->cd(3)->SetFillColor(3);
  lC2->cd(3)->SetTickx(1);
  lC2->cd(3)->SetTicky(1);
  lC2->cd(3)->SetLogx(1);

  TLegend *rL3 =new TLegend(0.18,0.28,0.60,0.52); rL3->SetFillColor(0); rL3->SetBorderSize(0);
  rL3-> SetNColumns(2);
  rL3->AddEntry(tgFEWZStatRatioBand,"FEWZ stat","FP");
  rL3->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  rL3->AddEntry(tgFEWZPDFRatioBand,"FEWZ PDF","FP");
  rL3->AddEntry(tgFEWZScaleRatioBand,"FEWZ scales","FP");
  rL3->SetTextSize(0.06);

  TLegend *tL3 =new TLegend(0.18,0.82,0.34,0.92); tL3->SetFillColor(0); tL3->SetBorderSize(0);
  tL3->AddEntry(tgFEWZStatRatioBand,"FEWZ","");
  tL3->SetTextSize(0.09);
  tL3->SetTextFont(2);

  tgDataStatSystRatioBand->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  tgDataStatSystRatioBand->GetXaxis()->SetTitleOffset(1.5);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleSize(33);
  tgDataStatSystRatioBand->Draw("2 A");
  tgFEWZScaleRatioBand->Draw("2 P");
  tgFEWZPDFRatioBand->Draw("2");
  tgFEWZStatRatioBand->Draw("2");
  rL3->Draw();
  tL3->Draw();

  sprintf(tmpName,"ZptNormDiffXsecRatio.");
  lC2->SaveAs(tmpName+format);
  
  return 0;
}

//TGraphAsymmErrors* Powheg(){
//
//  const int n = 18;
//
//  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
//  double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
//  double x[n];
//
//  ///Powheg Pre fsr numbers
//  double y[18] = {
//    0.0232823816,
//    0.0534414648, 
//    0.0556310835, 
//    0.0454828550, 
//    0.0356700536, 
//    0.0280575476, 
//    0.0226186451, 
//    0.0182287311, 
//    0.0118469402, 
//    0.0064591971, 
//    0.0037984679, 
//    0.0019372698, 
//    0.0008350090, 
//    0.0003902566, 
//    0.0001503083, 
//    0.0000509369, 
//    0.0000184618, 
//    0.0000014680};
//
//  //// stat error is normalized with Toy RMS method.
//  double errorStat_y[18] = {
//    0.0001395273, 
//    0.0002032366, 
//    0.0002069700, 
//    0.0001884388, 
//    0.0001705318, 
//    0.0001524742, 
//    0.0001370749,
//    0.0001241799, 
//    0.0000484067, 
//    0.0000366260, 
//    0.0000284210, 
//    0.0000144427, 
//    0.0000095287, 
//    0.0000065614, 
//    0.0000028981, 
//    0.0000016920, 
//    0.0000008289, 
//    0.0000000970};
//
//  //// powheg norm PDF error in %, convert it to number
//  double errorPDF_y[18] = {
//    0.412089 , 
//    0.31525  , 
//    0.209474 , 
//    0.188146 , 
//    0.250735 , 
//    0.21012  , 
//    0.351397 ,
//    0.284557 , 
//    0.304525 , 
//    0.315051 , 
//    0.356816 , 
//    0.619421 , 
//    1.16534  , 
//    1.56789  , 
//    2.29144  , 
//    2.92386  , 
//    4.04627  , 
//    6.01494  };
//
//  double error_y[n];
//  for (int i=0;i<n;i++)
//  { 
//    errorPDF_y[i]= 0.01*errorPDF_y[i]*y[i];
//    error_y[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_y[i]*errorPDF_y[i]); 
//    //std::cout<<"Powhed:  errorPDF_y[i]\t"<<errorPDF_y[i]<<"\terrorStat_y[i]\t"<<errorStat_y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
//    //std::cout<<"Powhed:  errorPDF_y[i]\t"<<100*errorPDF_y[i]/y[i]<<"\terrorStat_y[i]\t"<<100*errorStat_y[i]/y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
//    //printf("Powheg PDF% : %.2f \t Stat % : %.2f \n ",100*errorPDF_y[i]/y[i] , 100*errorStat_y[i]/y[i]);
//    std::cout << "Powheg PDF : " << 100*errorPDF_y[i]/y[i] << "\t Stat : " << 100*errorStat_y[i]/y[i] << endl;
//  }
//
//  //double error_x[18];
//
//  //double error_x[18];
//  double error_xlow[n];
//  double error_xhigh[n];
//
//  for (int i=0;i<n;i++) {
//    error_xlow[i] = (point[i] - bin[i]);
//    error_xhigh[i] = (bin[i+1]-point[i]);
//  }
//
//  //for (int i=0;i<n;i++) error_x[i]=0;
//
//  //TGraphAsymmErrors* grPowheg = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
//  TGraphAsymmErrors* grPowheg = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_y,error_y);
//
//
//  return grPowheg;
//}
//
//TGraphAsymmErrors* Resbos(){
//
//  const int n = 18;
//
//  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
//  double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
//  double x[n];
//
//
//  double y[18] = {
//    // Nadeesha PAS
//    //        0.02957739, 
//    //	0.05620224, 
//    //	0.04961592, 
//    //	0.04050349, 
//    //	0.03325381, 
//    //	0.02725463, 
//    //	0.02253389, 
//    //	0.01862596, 
//    //	0.01237984, 
//    //	0.00664593, 
//    //	0.00393008, 
//    //	0.00200304, 
//    //	0.00085508, 
//    //	0.00041032, 
//    //	0.00017656, 
//    //	0.00006332, 
//    //	0.00002423, 
//    //	0.00000115};
//
//    // Hammid resbos cp
//    //0.03223331, 
//    //0.06286445,
//    //0.05278694,
//    //0.04188681,
//    //0.03321407,
//    //0.02650710,
//    //0.02134491,
//    //0.01716620,
//    //0.01144344,
//    //0.00580740,
//    //0.00346440,
//    //0.00182813,
//    //0.00082492,
//    //0.00040398,
//    //0.00017512,
//    //0.00006440,
//    //0.00002296,
//    //0.00000211}; 
//
//    // Hammid resbos p
//    // 0.03801042,	
//    // 0.05823305,	
//    // 0.05171227,	
//    // 0.04158588,	
//    // 0.03318528,	
//    // 0.02647431,	
//    // 0.02130071,	
//    // 0.01728264,	
//    // 0.01103833,	
//    // 0.00592830,	
//    // 0.00353564,	
//    // 0.00188739,	
//    // 0.00085679,	
//    // 0.00042180,	
//    // 0.00018299,	
//    // 0.00006782,	
//    // 0.00002363,	
//    // 0.00000213};	
//
//    // Hammid resbos latest
//  0.03335570, 	
//  0.05909690, 	
//  0.05080043, 	
//  0.04145975, 	
//  0.03394089, 	
//  0.02773926, 	
//  0.02220695, 	
//  0.01779636, 	
//  0.01131978, 	
//  0.00601003, 	
//  0.00357902, 	
//  0.00188844, 	
//  0.00085224, 	
//  0.00041736, 	
//  0.00018092, 	
//  0.00006654, 	
//  0.00002372, 	
//  0.00000218
//  };
//// Nadeesha error	
////double error_y[18] = {0.00514001, 
////                      0.00707116,
////                      0.00664392, 
////                      0.00600288,
////                      0.00543919, 
////                      0.00492418, 
////                      0.00447746, 
////                      0.00407074, 
////                      0.00165936, 
////                      0.00121580, 
////                      0.00093494, 
////                      0.00047197, 
////                      0.00030837, 
////                      0.00021361, 
////                      0.00009908, 
////                      0.00005934, 
////                      0.00002997, 
////                      0.00000270};
//double error_y[18] = {
//  0.001694   , 
//  0.002232   ,
//  0.001961   ,
//  0.001568   , 
//  0.001186   ,
//  0.0008501  ,
//  0.0003468  ,
//  0.0006431  ,
//  0.0008231  ,
//  0.0005566  ,
//  0.0002914  ,
//  0.0001034  ,
//  5.643e-05  ,
//  3.159e-05  ,
//  1.899e-05  ,
//  7.658e-06  ,
//  2.876e-06  ,
//  2.919e-07  
//};
//
//
//
////double error_x[18];
//double error_xlow[n];
//double error_xhigh[n];
//
//for (int i=0;i<n;i++) {
//  error_xlow[i] = (point[i] - bin[i]);
//  error_xhigh[i] = (bin[i+1]-point[i]);
//}
//
////for (int i=0;i<n;i++) error_x[i]=0;
//
////TGraphAsymmErrors* grResbos = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
//TGraphAsymmErrors* grResbos = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_y,error_y);
//
//
//return grResbos;
//}
//
//
//
//
//TGraphAsymmErrors* Madgraph(){
//
//  const int n = 18;
//
//  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
//  double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
//  double x[n];
//  //  for (int i=0;i<n;i++) {
//  //  x[i] = bin[i] + (bin[i+1]-bin[i])/2;
//  // }
//
//  /*
//  // PostFSR and StatErr
//  double y[18] = {0.02838694, 0.05640728, 0.05274386, 0.04305959, 0.03456624, 0.02771716, 0.02242516, 0.01798233, 0.01161362, 0.00624185, 0.00375409, 0.00194832, 0.00084676, 0.00041085, 0.00018461, 0.00006234, 0.00002194, 0.00000187};
//  double error_y[18] = {0.00022293, 0.00031362, 0.00030326, 0.00027401, 0.00024550, 0.00021984, 0.00019774, 0.00017707, 0.00007115, 0.00005216, 0.00004045, 0.00002061, 0.00001359, 0.00000946, 0.00000449, 0.00000261, 0.00000126, 0.00000015};
//  */
//  // PreFSR and Total Uncer(sqrt(StatErr^2 + PDFErr^2))
//  double y[18] = { 0.0296683  , 
//    0.058359   , 
//    0.0532516  , 
//    0.042838   , 
//    0.03381    , 
//    0.0268681  , 
//    0.0216179  , 
//    0.0171872  , 
//    0.0113335  , 
//    0.00623927 , 
//    0.00378249 , 
//    0.0019837  , 
//    0.000862846,
//    0.000418884, 
//    0.000188615, 
//    6.45009e-05, 
//    2.23858e-05, 
//    1.91878e-06};
//  double errorStat_y[18] = {0.0002159060, 
//    0.0002914961, 
//    0.0002808667, 
//    0.0002537420, 
//    0.0002301006, 
//    0.0002065859, 
//    0.0001854926, 
//    0.0001669025, 
//    0.0000656459, 
//    0.0000498163, 
//    0.0000392100, 
//    0.0000201929, 
//    0.0000133864, 
//    0.0000093946, 
//    0.0000044844, 
//    0.0000026313, 
//    0.0000012616, 
//    0.0000001533};
//  /////normalized PDF syst in %, convert it to number
//  double errorPDF_y[18] = {1.87987 , 
//    1.74473 , 
//    1.31112 , 
//    1.15229 , 
//    0.797034, 
//    0.369078, 
//    0.275328, 
//    0.691266, 
//    1.29233 , 
//    1.74525 , 
//    1.93177 , 
//    2.09596 , 
//    2.31831 , 
//    2.62671 , 
//    2.7677  , 
//    3.15532 , 
//    3.78032 , 
//    4.67155 };
//
//
//  double error_y[n];
//  for (int i=0;i<n;i++)
//  { 
//    errorPDF_y[i]= 0.01*errorPDF_y[i]*y[i];
//    error_y[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_y[i]*errorPDF_y[i]); 
//    //std::cout<<"Madgraph:   errorPDF_y[i]\t"<<errorPDF_y[i]<<"\terrorStat_y[i]\t"<<errorStat_y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
//    //std::cout<<"Madgraph:  errorPDF_y[i]\t"<<100*errorPDF_y[i]/y[i]<<"\terrorStat_y[i]\t"<<100*errorStat_y[i]/y[i]<<"\t error_y[i]\t"<<error_y[i]<< std::endl;
//    //printf("Madgraph PDF % : %.2f \t Stat % : %.2f \n ",100*errorPDF_y[i]/y[i] , 100*errorStat_y[i]/y[i]);
//  }
//
//
//  //double error_x[18];
//  double error_xlow[n];
//  double error_xhigh[n];
//
//  for (int i=0;i<n;i++) {
//    error_xlow[i] = (point[i] - bin[i]);
//    error_xhigh[i] = (bin[i+1]-point[i]);
//  }
//
//  //for (int i=0;i<n;i++) error_x[i]=0;
//
//  //TGraphAsymmErrors* grMadgraph = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
//  TGraphAsymmErrors* grMadgraph = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_y,error_y);
//
//
//  return grMadgraph;
//}
//
//TGraphAsymmErrors* FEWZ(){
//
//  const int n = 18;
//
//  double bin[n+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
//  double point[18] = {1.25, 3.75, 6.25, 8.75, 11.75, 13.75, 16.25, 18.75, 24.25, 33.75, 44.75, 59.00, 79.00, 99.00, 126.00, 168.00, 214.00, 337.50};
//  double x[n];
//
//
//  //double y[18] = {0,0,0,0,0,0,0,0,0.010946908,0.005861624,0.00351576,0.00186525,0.000836974,0.000406356,0.000176403,6.37816E-05,2.25044E-05,2.08287E-06};
//  //double error_ylow[18] = {0,0,0,0,0,0,0,0,0.000260501,0.00013212,7.60136E-05,3.86252E-05,1.65916E-05,7.78994E-06,3.33908E-06,1.20556E-06,3.89959E-07,4.13802E-08};
//  //double error_yhigh[18] = {0,0,0,0,0,0,0,0,0.00031118,0.00015371,8.64209E-05,4.26788E-05,1.75198E-05,8.01257E-06,3.36195E-06,1.27319E-06,4.68827E-07,5.67557E-08};
//  double y[18];
//  y[0] = 100;//-0.08005877; 
//  y[1] = 0.12748182 ; 
//  y[2] = 0.07877656 ; 
//  y[3] = 0.05332600 ; 
//  y[4] = 0.03845661 ; 
//  y[5] = 0.02904814 ; 
//  y[6] = 0.02265475 ; 
//  y[7] = 0.01805065 ; 
//  y[8] = 0.01135670 ; 
//  y[9] = 0.00594296 ; 
//  y[10] = 0.00353301; 
//  y[11] = 0.00185758; 
//  y[12] = 0.00082569; 
//  y[13] = 0.00039315; 
//  y[14] = 0.00016611; 
//  y[15] = 0.00005898; 
//  y[16] = 0.00001959; 
//  y[17] = 0.00000179; 
//
//  double errorStat_y[18] = {
//    0.0,//0.00094095, 
//    0.00039520, 
//    0.00024179,
//    0.00016486,
//    0.00012112,
//    0.00009292,
//    0.00007464,
//    0.00006809,
//    0.00003028,
//    0.00001669,
//    0.00001118,
//    0.00000535,
//    0.00000286,
//    0.00000170,
//    0.00000067,
//    0.00000082,
//    0.00000031,
//    0.00000004};
//  /////normalized PDF syst in %, convert it to number
//  double errorPDF_ylow[18] = { 
//    0.0,//-6.7678 , 
//    1.0464, 
//    0.7600 ,
//    1.0219 ,
//    1.2630 ,
//    1.4623 ,
//    1.6461 ,
//    1.8262 ,
//    2.1995 ,
//    2.7651 ,
//    3.2279 ,
//    3.8587 ,
//    4.6700 ,
//    5.4480 ,
//    6.4558 ,
//    7.8312 ,
//    9.3764 ,
//    11.4845}; 
//  double errorPDF_yhigh[18] = {
//    0.0,     // -6.4288, 
//    1.8289 , 
//    1.3007  ,
//    1.2656  ,
//    1.3277  ,
//    1.4084  ,
//    1.5036  ,
//    1.6060  ,
//    1.8347  ,
//    2.2177  ,
//    2.5401 ,
//    2.9723 ,
//    3.4925 ,
//    3.9540 ,
//    4.4901 ,
//    5.1299 ,
//    5.7526 ,
//    6.6129 }; 
//  double errorScale[18] = {
//    0.0,//0.00842947        , 
//    0.0106515         ,
//    0.00266521 ,
//    0.000428617,
//    0.000114661,
//    0.000205532,
//    0.000345275,
//    0.000297941,
//    0.000323374,
//    0.000210461,
//    0.000147318,
//    9.51447e-05,
//    4.73258e-05,
//    2.47657e-05,
//    1.03102e-05,
//    4.30353e-06,
//    1.45553e-06,
//    1.5117e-07 };
//
//  double error_yhigh[n];
//  double error_ylow[n];
//  for (int i=0;i<n;i++)
//  { 
//    errorPDF_ylow[i]= 0.01*errorPDF_ylow[i]*y[i];
//    errorPDF_yhigh[i]= 0.01*errorPDF_yhigh[i]*y[i];
//    error_ylow[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_ylow[i]*errorPDF_ylow[i] + errorScale[i]*errorScale[i]); 
//    error_yhigh[i] = sqrt(errorStat_y[i]*errorStat_y[i] + errorPDF_yhigh[i]*errorPDF_yhigh[i] + errorScale[i]*errorScale[i]); 
    //std::cout<<"FEWZ:   errorPDF_ylow[i]\t"<<errorPDF_ylow[i]<<"\terrorStat_y[i]\t"<<errorStat_y[i]<<"\t error_ylow[i]\t"<<error_ylow[i]<< std::endl;
    //std::cout<<"FEWZ:  errorPDF_y[i]\t"<<100*TMath::Max(errorPDF_ylow[i],errorPDF_yhigh[i])/y[i]<<"\terrorStat_y[i]\t"<<100*errorStat_y[i]/y[i]<<"\terrorScale_y[i]\t"<<100*errorScale[i]/y[i]<< std::endl;
    //printf("FEWZ PDF % : %.2f \t Stat % : %.2f Scale % : %.2f \n ",100*TMath::Max(errorPDF_ylow[i],errorPDF_yhigh[i])/y[i] , 100*errorStat_y[i]/y[i], 100*errorScale[i]/y[i]);
//  }

//  //double error_x[18];
//  double error_xlow[n];
//  double error_xhigh[n];
//
//  for (int i=0;i<n;i++) {
//    error_xlow[i] = (point[i] - bin[i]);
//    error_xhigh[i] = (bin[i+1]-point[i]);
//  }
//
//  //double error_y[18];
//  //for (int i=0;i<n;i++) error_x[i]=0;
//  //for (int i=0;i<n;i++) error_y[i] = (error_yhigh[i]>error_ylow[i]) ? error_yhigh[i] : error_ylow[i];
//
//  //TGraphAsymmErrors* grFEWZ = new TGraphAsymmErrors(n,point,y,error_x,error_x,error_y,error_y);
//  TGraphAsymmErrors* grFEWZ = new TGraphAsymmErrors(n,point,y,error_xlow,error_xhigh,error_ylow,error_yhigh);
//
//
//  return grFEWZ;
//}
