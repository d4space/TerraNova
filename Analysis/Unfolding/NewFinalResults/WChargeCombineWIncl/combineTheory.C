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
    f_Fewz_2 = new TFile("../../../RstFEWZ/WmToEleNu_13bin_dynamic_NNLO.root");
    f_Resbos_1 = new TFile("../../../RstResbos/Resbos_WpToEleNu.root");
    f_Resbos_2 = new TFile("../../../RstResbos/Resbos_WmToEleNu.root");
  }
  
  TH1D* lFewz_1;
  TH1D* lFewz_2;
  
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

  lFewz_1 = (TH1D*)f_Fewz_1->Get("hxsec")->Clone("lFewz_1");
  lFewz_2 = (TH1D*)f_Fewz_2->Get("hxsec")->Clone("lFewz_2");
  Fewz1_PDFErr = (TH1D*)f_Fewz_1->Get("PDFErr")->Clone("Fewz_PDFErr");
  Fewz2_PDFErr = (TH1D*)f_Fewz_2->Get("PDFErr")->Clone("Fewz_PDFErr");
  Fewz1_ScaleErr = (TH1D*)f_Fewz_1->Get("ScaleErr")->Clone("Fewz_ScaleErr");
  Fewz2_ScaleErr = (TH1D*)f_Fewz_2->Get("ScaleErr")->Clone("Fewz_ScaleErr");

  Double_t MC1[nBins-1],MC1StatErr[nBins-1],MC1PDFErr[nBins-1],MC1ScaleErr;
  Double_t MC2[nBins-1],MC2StatErr[nBins-1],MC2PDFErr[nBins-1],MC2ScaleErr;
  Double_t MCmean[nBins-1],MCmeanStatErr[nBins-1],MCmeanPDFErr[nBins-1],MCmeanScaleErr[nBins-1],MCmeanTotalErr[nBins-1];

  Double_t Resb1[7], Resb1Err[7];
  Double_t Resb2[7], Resb2Err[7];
  Double_t Resb[7], ResbErr[7];

  TH1D* hxsec = new TH1D("hxsec","hxsec",13,0,13);
  TH1D* PDFError = new TH1D("PDFErr","PDFErr",13,0,13);
  TH1D* ScaleError = new TH1D("ScaleErr","ScaleErr",13,0,13);
  TH1D* hResbos29 = new TH1D("hResbos29","hResbos29",13,0,13);
  TH1D* hResbos30 = new TH1D("hResbos30","hResbos30",13,0,13);
  TH1D* hResbos31 = new TH1D("hResbos31","hResbos31",13,0,13);
  TH1D* hResbos32 = new TH1D("hResbos32","hResbos32",13,0,13);
  TH1D* hResbos33 = new TH1D("hResbos33","hResbos33",13,0,13);
  TH1D* hResbos34 = new TH1D("hResbos34","hResbos34",13,0,13);
  TH1D* hResbos35 = new TH1D("hResbos35","hResbos35",13,0,13);

  //if(BaseName=="WInclToMuNu")
  //  TFile f_out("Theory_Muon.root","recreate");
  //if(BaseName=="WInclToEleNu")
  //  //TFile f_out("Theory_Ele.root","recreate");
  //  TFile f_out("ResultWInclToEleNu/Result_WInclToEleNu_Theory.root","recreate");
  
  tmpTStr = resultDir+"/Result_"+BaseName+"_Theory.root";
  TFile *f_Out    = new TFile(tmpTStr,"recreate");

  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    MC1[ipt] = lFewz_1->GetBinContent(ipt+1);
    MC1StatErr[ipt] = lFewz_1->GetBinError(ipt+1);
    MC1PDFErr[ipt] = MC1[ipt]*(Fewz1_PDFErr->GetBinError(ipt+1))*0.01;
    MC1ScaleErr[ipt] = MC1[ipt]*(Fewz1_ScaleErr->GetBinError(ipt+1))*0.01;
    
    MC2[ipt] = lFewz_2->GetBinContent(ipt+1);
    MC2StatErr[ipt] = lFewz_2->GetBinError(ipt+1);
    MC2PDFErr[ipt] = MC2[ipt]*(Fewz2_PDFErr->GetBinError(ipt+1))*0.01;
    MC2ScaleErr[ipt] = MC2[ipt]*(Fewz2_ScaleErr->GetBinError(ipt+1))*0.01;
    
    MCmean[ipt]=MC1[ipt]+MC2[ipt];
    MCmeanStatErr[ipt]=sqrt(MC1StatErr[ipt]**2 + MC2StatErr[ipt]**2);
    MCmeanPDFErr[ipt]=MC1PDFErr[ipt] + MC2PDFErr[ipt];
    MCmeanScaleErr[ipt]=MC1ScaleErr[ipt] + MC2ScaleErr[ipt];
    MCmeanTotalErr[ipt] = sqrt(MCmeanPDFErr[ipt]*MCmeanPDFErr[ipt] + MCmeanScaleErr[ipt]*MCmeanScaleErr[ipt]);
    
    for( int i(0);i<7;i++)
    {
      Resb1[i] = lResbos1[i]->GetBinContent(ipt+1);
      Resb1Err[i] = lResbos1[i]->GetBinError(ipt+1);
      Resb2[i] = lResbos2[i]->GetBinContent(ipt+1);
      Resb2Err[i] = lResbos2[i]->GetBinError(ipt+1);
      ResbErr[i]=sqrt((Resb1Err[i]*Resb1Err[i])+(Resb2Err[i]*Resb2Err[i]));
      Resb[i]=Resb1[i]+Resb2[i];
    }
    
    hResbos29 -> SetBinContent(ipt+1, Resb[0]);
    hResbos30 -> SetBinContent(ipt+1, Resb[1]);
    hResbos31 -> SetBinContent(ipt+1, Resb[2]);
    hResbos32 -> SetBinContent(ipt+1, Resb[3]);
    hResbos33 -> SetBinContent(ipt+1, Resb[4]);
    hResbos34 -> SetBinContent(ipt+1, Resb[5]);
    hResbos35 -> SetBinContent(ipt+1, Resb[6]);
    
    hResbos29 -> SetBinError(ipt+1, ResbErr[0]);
    hResbos30 -> SetBinError(ipt+1, ResbErr[1]);
    hResbos31 -> SetBinError(ipt+1, ResbErr[2]);
    hResbos32 -> SetBinError(ipt+1, ResbErr[3]);
    hResbos33 -> SetBinError(ipt+1, ResbErr[4]);
    hResbos34 -> SetBinError(ipt+1, ResbErr[5]);
    hResbos35 -> SetBinError(ipt+1, ResbErr[6]);
    
    hxsec->SetBinContent(ipt+1, MCmean[ipt]);
    hxsec->SetBinError(ipt+1, MCmeanStatErr[ipt]);
    //hxsec->SetBinError(ipt+1, MCmeanTotalErr[ipt]);
    PDFError->SetBinError(ipt+1,MCmeanPDFErr[ipt]);
    ScaleError->SetBinError(ipt+1,MCmeanScaleErr[ipt]);
  }
  hxsec->Write();
  PDFError->Write();
  ScaleError->Write();
  
  hResbos29 -> Write();
  hResbos30 -> Write();
  hResbos31 -> Write();
  hResbos32 -> Write();
  hResbos33 -> Write();
  hResbos34 -> Write();
  hResbos35 -> Write();

  return 0;
}
