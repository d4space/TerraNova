#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>

const int nBins = 13;
double Wpt12Bins[nBins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

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
    f_Fewz_1 = new TFile("../../RstFEWZ_12Bin_Fiducial_dynamic/WpToMuNu_dynamic_NNLO.root");
    f_Fewz_2 = new TFile("../../RstFEWZ_12Bin_Fiducial_dynamic/WmToMuNu_dynamic_NNLO.root");

    f_Resbos_1 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wplus.root");
    f_Resbos_2 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wminus.root");

  }
  
  TH1D* lFewz_1;
  TH1D* lFewz_2;
  lFewz_1 = (TH1D*)f_Fewz_1->Get("hxsec")->Clone("lFewz_1");
  lFewz_2 = (TH1D*)f_Fewz_2->Get("hxsec")->Clone("lFewz_2");
   TH1D *lFewz_1_12Bin   = new TH1D("","",nBins-1,Wpt12Bins);lFewz_1_12Bin->Sumw2();
   TH1D *lFewz_2_12Bin   = new TH1D("","",nBins-1,Wpt12Bins);lFewz_2_12Bin->Sumw2();

  //FEWZ Errors
   //TH1D* Scale1 = (TH1D*)f_Fewz_1->Get("hScaleErr")->Clone("Scale1");
   //TH1D* PDFEr1   = (TH1D*)f_Fewz_1->Get("PDFErr")->Clone("PDF1");
   TH1D* Scale1 = (TH1D*)f_Fewz_1->Get("Norm_ScaleErr")->Clone("Scale1");
   TH1D* PDFEr1   = (TH1D*)f_Fewz_1->Get("Norm_PDFErr")->Clone("PDF1");
   
   TH1D* Scale2 = (TH1D*)f_Fewz_2->Get("Norm_ScaleErr")->Clone("Scale2");
   TH1D* PDF2   = (TH1D*)f_Fewz_2->Get("Norm_PDFErr")->Clone("PDF2");
 
   TH1D *ScaleError_12Bin   = new TH1D("ScaleErr","ScaleErr",nBins-1,Wpt12Bins);ScaleError_12Bin->Sumw2();
   TH1D *PDFError_12Bin   = new TH1D("PDFErr","PDFErr",nBins-1,Wpt12Bins);PDFError_12Bin->Sumw2();
   
   Double_t Scale11[nBins-1],PDF11[nBins-1];
   Double_t Scale22[nBins-1],PDF22[nBins-1];

   for( int ipt(1);ipt<nBins;ipt++)
    {
	lFewz_1_12Bin->SetBinContent(ipt,lFewz_1->GetBinContent(ipt));
	lFewz_2_12Bin->SetBinContent(ipt,lFewz_2->GetBinContent(ipt));
       

	Scale11[ipt]=0.01*Scale1->GetBinError(ipt)*lFewz_1->GetBinContent(ipt);
	Scale22[ipt]=0.01*Scale2->GetBinError(ipt)*lFewz_2->GetBinContent(ipt);
	PDF11[ipt]=0.01*PDF1->GetBinError(ipt)*lFewz_1->GetBinContent(ipt);
	PDF22[ipt]=0.01*PDF1->GetBinError(ipt)*lFewz_2->GetBinContent(ipt);

	//ScaleError_12Bin->SetBinError(ipt, sqrt(Scale11[ipt]*Scale11[ipt] + Scale22[ipt]*Scale22[ipt]));
	ScaleError_12Bin->SetBinError(ipt, Scale11[ipt] + Scale22[ipt]);
	//ScaleError_12Bin->SetBinError(ipt, 0);
	PDFError_12Bin->SetBinError(ipt, sqrt(PDF11[ipt]*PDF11[ipt] + PDF22[ipt]*PDF22[ipt]));

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


  Double_t MC1[nBins-1],MC1Err[nBins-1];
  Double_t MC2[nBins-1],MC2Err[nBins-1];
  Double_t MCmean[nBins-1],MCmeanErr[nBins-1];

  Double_t Resb1[7], Resb1Err[7];
  Double_t Resb2[7], Resb2Err[7];
  Double_t Resb[7], ResbErr[7];

  TH1D* hxsec = new TH1D("hxsec","hxsec",12,0,12);hxsec->Sumw2();
  TH1D* hResbos29 = new TH1D("hResbos29","hResbos29",12,0,12);hResbos29->Sumw2();
  TH1D* hResbos30 = new TH1D("hResbos30","hResbos30",12,0,12);hResbos30->Sumw2();
  TH1D* hResbos31 = new TH1D("hResbos31","hResbos31",12,0,12);hResbos31->Sumw2();
  TH1D* hResbos32 = new TH1D("hResbos32","hResbos32",12,0,12);hResbos32->Sumw2();
  TH1D* hResbos33 = new TH1D("hResbos33","hResbos33",12,0,12);hResbos33->Sumw2();
  TH1D* hResbos34 = new TH1D("hResbos34","hResbos34",12,0,12);hResbos34->Sumw2();
  TH1D* hResbos35 = new TH1D("hResbos35","hResbos35",12,0,12);hResbos35->Sumw2();

  //if(BaseName=="WInclToMuNu")
  //  TFile f_out("Theory_Muon.root","recreate");
  //if(BaseName=="WInclToEleNu")
  //  //TFile f_out("Theory_Ele.root","recreate");
  //  TFile f_out("ResultWInclToEleNu/Result_WInclToEleNu_Theory.root","recreate");
  
  tmpTStr = resultDir+"/Result_"+BaseName+"_Theory.root";
  TFile *f_Out    = new TFile(tmpTStr,"recreate");

  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    MC1[ipt] = lFewz_1_12Bin->GetBinContent(ipt+1);
    MC2[ipt] = lFewz_2_12Bin->GetBinContent(ipt+1);
    MCmean[ipt]=MC1[ipt]+MC2[ipt];
    
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
    //hxsec->SetBinError(ipt+1, sqrt(1.0/MCmeanErr[ipt]));
    hxsec->SetBinError(ipt+1, MCmeanErr[ipt]);
  }
  hxsec->Write();
  ScaleError_12Bin->Write();
  PDFError_12Bin->Write();
 
  hResbos29 -> Write();
  hResbos30 -> Write();
  hResbos31 -> Write();
  hResbos32 -> Write();
  hResbos33 -> Write();
  hResbos34 -> Write();
  hResbos35 -> Write();

  return 0;
}
