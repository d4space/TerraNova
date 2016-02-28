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

  double ResbosPDF[18];
  ResbosPDF[0]		=	0.00110074 ;
  ResbosPDF[1]		=       0.0014006  ; 	
  ResbosPDF[2]		=       0.000624845; 	
  ResbosPDF[3]		=       0.000302656; 	
  ResbosPDF[4]		=       0.00024098 ; 	
  ResbosPDF[5]		=       0.000257975; 	
  ResbosPDF[6]		=       0.000322001; 	
  ResbosPDF[7]		=       0.000284742; 	
  ResbosPDF[8]		=       0.000196964; 	
  ResbosPDF[9]		=       0.000123807; 	
  ResbosPDF[10]		=       8.69702e-05; 	
  ResbosPDF[11]		=       5.13656e-05; 	
  ResbosPDF[12]		=       2.81239e-05; 	
  ResbosPDF[13]		=       1.49415e-05; 	
  ResbosPDF[14]		=       8.08712e-06; 	
  ResbosPDF[15]		=       3.5732e-06 ; 	
  ResbosPDF[16]		=       1.51808e-06; 	
  ResbosPDF[17]		=       1.83338e-07;

  double ResbosScale[18];
  ResbosScale[0]	=	0.00128753 ;
  ResbosScale[1]	=       0.00173745 ; 	
  ResbosScale[2]	=       0.0018593  ; 	
  ResbosScale[3]	=       0.00153816 ; 	
  ResbosScale[4]	=       0.00116078 ; 	
  ResbosScale[5]	=       0.000809986; 	
  ResbosScale[6]	=       0.0001288  ; 	
  ResbosScale[7]	=       0.000576602; 	
  ResbosScale[8]	=       0.000799176; 	
  ResbosScale[9]	=       0.000542706; 	
  ResbosScale[10]	=       0.00027809 ; 	
  ResbosScale[11]	=       8.97009e-05; 	
  ResbosScale[12]	=       4.89186e-05; 	
  ResbosScale[13]	=       2.78379e-05; 	
  ResbosScale[14]	=       1.71874e-05; 	
  ResbosScale[15]	=       6.77377e-06; 	
  ResbosScale[16]	=       2.44316e-06; 	
  ResbosScale[17]	=       2.27156e-07;

  double ResbosTotalUnc[18];
  for(int i=0;i<18;i++)
  {
    ResbosTotalUnc[i] = sqrt(ResbosPDF[i]*ResbosPDF[i] + ResbosScale[i]*ResbosScale[i]);
  }

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
    cout << "Resbos : " << Resbos[i] << "\t PDF : " << ResbosPDF[i] << "\t Scale : " << ResbosScale[i] <<  endl;
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
  double ResbosPDFRatioBand[18];
  double ResbosScaleRatioBand[18];
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
    ResbosScaleRatioBand[i] = ResbosScale[i] / Data[i] ;
    ResbosPDFRatioBand[i] = (ResbosPDF[i] + ResbosScale[i]) / Data[i] ;
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
  
  TGraphAsymmErrors *tgResbos = new TGraphAsymmErrors(n,point,Resbos,error_xlow,error_xhigh,ResbosTotalUnc,ResbosTotalUnc);
  TGraphAsymmErrors *tgResbosPDFRatioBand = new TGraphAsymmErrors(n,point,ResbosRatioBand,error_xlow,error_xhigh,ResbosPDFRatioBand,ResbosPDFRatioBand);
  TGraphAsymmErrors *tgResbosScaleRatioBand = new TGraphAsymmErrors(n,point,ResbosRatioBand,error_xlow,error_xhigh,ResbosScaleRatioBand,ResbosScaleRatioBand);
  
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


  TLegend *lL =new TLegend(0.25,0.17,0.65,0.45); lL->SetFillColor(0); lL->SetBorderSize(0);
  lL->AddEntry(tgData,"Data","PLE1");
  lL->AddEntry(tgResbos,"ResBos-CP CT10 NNLL","f");
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
  tgResbosPDFRatioBand->SetMarkerColor(MarkerColor_ResBos);
  tgResbosPDFRatioBand->SetMarkerStyle(MarkerStyle_ResBos);
  tgResbosPDFRatioBand->SetFillColor(BandColor_ResBos_PDF);

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
  
  TLegend *lResbos =new TLegend(0.18,0.04,0.60,0.30); lResbos->SetFillColor(0); lResbos->SetBorderSize(0);
  lResbos-> SetNColumns(2);
  lResbos->AddEntry(tgResbosScaleRatioBand,"ResBos scales","FP");
  lResbos->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  lResbos->AddEntry(tgResbosPDFRatioBand,"ResBos PDF","F");
  lResbos->SetTextSize(0.08);

  TPaveText *tResBos = new TPaveText(0.18,0.82,0.34,0.92,"NDC");
  tResBos->SetBorderSize(0);
  tResBos->SetFillStyle(0);
  tResBos->SetTextSize(0.12);
  tResBos->AddText("#font[42]{ResBos-CP}");
 
  TPaveText *tChannel = new TPaveText(0.35,0.82,0.67,0.93,"NDC");
  tChannel->SetBorderSize(0);
  tChannel->SetFillStyle(0);
  tChannel->SetTextSize(0.12);
  tChannel->AddText("Z #rightarrow #mu^{+}#mu^{#font[122]{-}}");

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
  tgResbosPDFRatioBand->Draw("2 P");
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
 
  tgDataStatSystRatioBand->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  tgDataStatSystRatioBand->GetXaxis()->SetTitleOffset(1.5);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleSize(33);
  tgDataStatSystRatioBand->Draw("2 A");
  tgFEWZScaleRatioBand->Draw("2");
  tgFEWZPDFRatioBand->Draw("2");
  tgFEWZStatRatioBand->Draw("2P");
  lFewz->Draw();
  tFewz->Draw();

  sprintf(tmpName,"ZptNormDiffXsecRatio.");
  lC2->SaveAs(tmpName+format);
  
  return 0;
}
