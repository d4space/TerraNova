#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>

const int nBins = 13;

int combineTheory(const TString BaseName)
{
  TString tmpTStr;
  char tmpName[30],tmpName_org[30];
  int Numb;

  TString resultDir = "Result"+BaseName;
  gSystem->mkdir(resultDir,kTRUE);


  TFile *f_Resbos_1;
  TFile *f_Resbos_2;

  if(BaseName == "WInclToMuNu")
  {
    f_Resbos_1 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wplus.root");
    f_Resbos_2 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wminus.root");
  }
  
  TH1D* lResbos1[7];
  TH1D* lResbos2[7];
  
  for( int i(0);i<7;i++)
  {
    Numb = 29+i;
    sprintf(tmpName_org,"hResbos%d",Numb);
    sprintf(tmpName,"lResbos1_%d",i);
    lResbos1[i] = (TH1D*)f_Resbos_1->Get(tmpName_org)->Clone(tmpName);
    sprintf(tmpName,"lResbos2_%d",i);
    lResbos2[i] = (TH1D*)f_Resbos_2->Get(tmpName_org)->Clone(tmpName);
  }

  double Resb1[7];
  double Resb2[7];
  double Resb[7];

  TH1D* hResbos29 = new TH1D("hResbos29","hResbos29",12,0,12);
  TH1D* hResbos30 = new TH1D("hResbos30","hResbos30",12,0,12);
  TH1D* hResbos31 = new TH1D("hResbos31","hResbos31",12,0,12);
  TH1D* hResbos32 = new TH1D("hResbos32","hResbos32",12,0,12);
  TH1D* hResbos33 = new TH1D("hResbos33","hResbos33",12,0,12);
  TH1D* hResbos34 = new TH1D("hResbos34","hResbos34",12,0,12);
  TH1D* hResbos35 = new TH1D("hResbos35","hResbos35",12,0,12);

  // make Histogram and save it to root file
  tmpTStr = resultDir+"/Result_"+BaseName+"_Theory.root";
  TFile *f_Out    = new TFile(tmpTStr,"recreate");

  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    for( int i(0);i<7;i++)
    {
      Resb1[i] = lResbos1[i]->GetBinContent(ipt+1);
      Resb2[i] = lResbos2[i]->GetBinContent(ipt+1);
      Resb[i]=Resb1[i]+Resb2[i];
    }
    
    hResbos29 -> SetBinContent(ipt+1, Resb[0]);
    hResbos30 -> SetBinContent(ipt+1, Resb[1]);
    hResbos31 -> SetBinContent(ipt+1, Resb[2]);
    hResbos32 -> SetBinContent(ipt+1, Resb[3]);
    hResbos33 -> SetBinContent(ipt+1, Resb[4]);
    hResbos34 -> SetBinContent(ipt+1, Resb[5]);
    hResbos35 -> SetBinContent(ipt+1, Resb[6]);
  }
  hResbos29 -> Write();
  hResbos30 -> Write();
  hResbos31 -> Write();
  hResbos32 -> Write();
  hResbos33 -> Write();
  hResbos34 -> Write();
  hResbos35 -> Write();

  return 0;
}
