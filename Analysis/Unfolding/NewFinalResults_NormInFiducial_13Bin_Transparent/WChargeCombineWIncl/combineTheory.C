#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>

const int nBins = 14;

int combineTheory(const TString BaseName)
{
  TString tmpTStr;
  char tmpName[30],tmpName_org[30];
  int Numb;

  TString resultDir = "Result"+BaseName;
  gSystem->mkdir(resultDir,kTRUE);


  TFile *f_Fewz_1;
  TFile *f_Fewz_2;
  TFile *f_Resbos_1;
  TFile *f_Resbos_2;

  if(BaseName == "WInclToMuNu")
  {
    f_Fewz_1 = new TFile("../../../RstFEWZ/WpToMuNu_13bin_dynamic_NNLO.root");
    f_Fewz_2 = new TFile("../../../RstFEWZ/WmToMuNu_13bin_dynamic_NNLO.root");
    f_Resbos_1 = new TFile("../../../RstResbos/Resbos_WpToMuNu.root");
    f_Resbos_2 = new TFile("../../../RstResbos/Resbos_WmToMuNu.root");
  }
  if(BaseName == "WInclToEleNu")
  {
    f_Fewz_1 = new TFile("../../../RstFEWZ/WpToEleNu_13bin_dynamic_NNLO.root");
    f_Fewz_2 = new TFile("../../../RstFEWZ/WpToEleNu_13bin_dynamic_NNLO.root");
    f_Resbos_1 = new TFile("../../../RstResbos/Resbos_WpToEleNu.root");
    f_Resbos_2 = new TFile("../../../RstResbos/Resbos_WmToEleNu.root");
  }
  
  TH1D* lFewz_Wp;
  TH1D* lFewz_Wp_up;
  TH1D* lFewz_Wp_down;
  TH1D* lFewz_Wm;
  TH1D* lFewz_Wm_up;
  TH1D* lFewz_Wm_down;
  
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

  lFewz_Wp = (TH1D*)f_Fewz_1->Get("hxsec")->Clone("lFewz_1");
  lFewz_Wp_up = (TH1D*)f_Fewz_1->Get("hxsec_up")->Clone("lFewz_1_up");
  lFewz_Wp_down = (TH1D*)f_Fewz_1->Get("hxsec_down")->Clone("lFewz_1_down");
  lFewz_Wm = (TH1D*)f_Fewz_2->Get("hxsec")->Clone("lFewz_2");
  lFewz_Wm_up = (TH1D*)f_Fewz_2->Get("hxsec_up")->Clone("lFewz_2_up");
  lFewz_Wm_down = (TH1D*)f_Fewz_2->Get("hxsec_down")->Clone("lFewz_2_down");

  double FEWZ_Wp[nBins-1];
  double FEWZ_Wp_StatErr[nBins-1];
  double FEWZ_Wp_up[nBins-1];
  double FEWZ_Wp_down[nBins-1];
  double FEWZ_Wm[nBins-1];
  double FEWZ_Wm_StatErr[nBins-1];
  double FEWZ_Wm_up[nBins-1];
  double FEWZ_Wm_down[nBins-1];
  double FEWZ_Wincl[nBins-1];
  double FEWZ_Wincl_StatErr[nBins-1];
  double FEWZ_Wincl_up[nBins-1];
  double FEWZ_Wincl_down[nBins-1];

  double Resb1[7];
  double Resb2[7];
  double Resb[7];

  TH1D* hxsec = new TH1D("hxsec","hxsec",13,0,13);
  TH1D* hxsec_up = new TH1D("hxsec_up","hxsec_up",13,0,13);
  TH1D* hxsec_down = new TH1D("hxsec_down","hxsec_down",13,0,13);

  TH1D* hResbos29 = new TH1D("hResbos29","hResbos29",13,0,13);
  TH1D* hResbos30 = new TH1D("hResbos30","hResbos30",13,0,13);
  TH1D* hResbos31 = new TH1D("hResbos31","hResbos31",13,0,13);
  TH1D* hResbos32 = new TH1D("hResbos32","hResbos32",13,0,13);
  TH1D* hResbos33 = new TH1D("hResbos33","hResbos33",13,0,13);
  TH1D* hResbos34 = new TH1D("hResbos34","hResbos34",13,0,13);
  TH1D* hResbos35 = new TH1D("hResbos35","hResbos35",13,0,13);

  // make Histogram and save it to root file
  tmpTStr = resultDir+"/Result_"+BaseName+"_Theory.root";
  TFile *f_Out    = new TFile(tmpTStr,"recreate");

  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    FEWZ_Wp[ipt] = lFewz_Wp->GetBinContent(ipt+1);
    FEWZ_Wp_StatErr[ipt] = lFewz_Wp->GetBinError(ipt+1);
    FEWZ_Wp_up[ipt] = lFewz_Wp_up->GetBinContent(ipt+1);
    FEWZ_Wp_down[ipt] = lFewz_Wp_down->GetBinContent(ipt+1);
    
    FEWZ_Wm[ipt] = lFewz_Wm->GetBinContent(ipt+1);
    FEWZ_Wm_StatErr[ipt] = lFewz_Wm->GetBinError(ipt+1);
    FEWZ_Wm_up[ipt] = lFewz_Wm_up->GetBinContent(ipt+1);
    FEWZ_Wm_down[ipt] = lFewz_Wm_down->GetBinContent(ipt+1);
    
    FEWZ_Wincl[ipt]=FEWZ_Wp[ipt]+FEWZ_Wm[ipt];
    FEWZ_Wincl_StatErr[ipt]=sqrt(FEWZ_Wp_StatErr[ipt]*FEWZ_Wp_StatErr[ipt] + FEWZ_Wm_StatErr[ipt]*FEWZ_Wm_StatErr[ipt]);
    FEWZ_Wincl_up[ipt]=FEWZ_Wp_up[ipt]+FEWZ_Wm_up[ipt];
    FEWZ_Wincl_down[ipt]=FEWZ_Wp_down[ipt]+FEWZ_Wm_down[ipt];
    
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
    
    hxsec->SetBinContent(ipt+1, FEWZ_Wincl[ipt]);
    hxsec->SetBinError(ipt+1, FEWZ_Wincl_StatErr[ipt]);
    hxsec_up->SetBinContent(ipt+1, FEWZ_Wincl_up[ipt]);
    hxsec_down->SetBinContent(ipt+1, FEWZ_Wincl_down[ipt]);
  }
  hxsec->Write();
  hxsec_up->Write();
  hxsec_down->Write();
  hResbos29 -> Write();
  hResbos30 -> Write();
  hResbos31 -> Write();
  hResbos32 -> Write();
  hResbos33 -> Write();
  hResbos34 -> Write();
  hResbos35 -> Write();

  return 0;
}
