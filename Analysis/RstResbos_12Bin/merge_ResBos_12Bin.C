#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "../Utils/const.h"

const int nBins = 14;
double WptLogBins[nBins] = {1.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};
double WptBins[nBins] = {0.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};

const int n12Bins = 13;
double WptLog12Bins[n12Bins] = {1.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double Wpt12Bins[n12Bins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

int merge_ResBos_12Bin(const TString BaseName)
{
  TString tmpTStr;
  char tmpName[30],tmpName_org[30];
  int Numb;

  TFile *f_Resbos;

  f_Resbos = new TFile("../RstResbos/Resbos_"+BaseName+".root");

  TH1D* lResbos29;
  TH1D* lResbos30;
  TH1D* lResbos31;
  TH1D* lResbos32;
  TH1D* lResbos33;
  TH1D* lResbos34;
  TH1D* lResbos35;
 
  lResbos29 = (TH1D*)f_Resbos->Get("hResbos29")->Clone("lResbos29");
  lResbos30 = (TH1D*)f_Resbos->Get("hResbos30")->Clone("lResbos30");
  lResbos31 = (TH1D*)f_Resbos->Get("hResbos31")->Clone("lResbos31");
  lResbos32 = (TH1D*)f_Resbos->Get("hResbos32")->Clone("lResbos32");
  lResbos33 = (TH1D*)f_Resbos->Get("hResbos33")->Clone("lResbos33");
  lResbos34 = (TH1D*)f_Resbos->Get("hResbos34")->Clone("lResbos34");
  lResbos35 = (TH1D*)f_Resbos->Get("hResbos35")->Clone("lResbos35");


  for(int i(1);i<14;i++)
  {
    cout<<i<<"\t"<<lResbos29->GetBinContent(i)<<endl;
  }
  
  
  TH1D *hResbos29 = new TH1D("hResbos29","hResbos29",12,Wpt12Bins);hResbos29->Sumw2();
  TH1D *hResbos30 = new TH1D("hResbos30","hResbos30",12,Wpt12Bins);hResbos30->Sumw2();
  TH1D *hResbos31 = new TH1D("hResbos31","hResbos31",12,Wpt12Bins);hResbos31->Sumw2();
  TH1D *hResbos32 = new TH1D("hResbos32","hResbos32",12,Wpt12Bins);hResbos32->Sumw2();
  TH1D *hResbos33 = new TH1D("hResbos33","hResbos33",12,Wpt12Bins);hResbos33->Sumw2();
  TH1D *hResbos34 = new TH1D("hResbos34","hResbos34",12,Wpt12Bins);hResbos34->Sumw2();
  TH1D *hResbos35 = new TH1D("hResbos35","hResbos35",12,Wpt12Bins);hResbos35->Sumw2();
  
  for(int i(1);i<4;i++)
  {
    hResbos29->SetBinContent(i,lResbos29->GetBinContent(i));
    hResbos30->SetBinContent(i,lResbos30->GetBinContent(i));
    hResbos31->SetBinContent(i,lResbos31->GetBinContent(i));
    hResbos32->SetBinContent(i,lResbos32->GetBinContent(i));
    hResbos33->SetBinContent(i,lResbos33->GetBinContent(i));
    hResbos34->SetBinContent(i,lResbos34->GetBinContent(i));
    hResbos35->SetBinContent(i,lResbos35->GetBinContent(i));
  }

    hResbos29->SetBinContent(4,lResbos29->GetBinContent(4)+lResbos29->GetBinContent(5));
    hResbos30->SetBinContent(4,lResbos30->GetBinContent(4)+lResbos30->GetBinContent(5));
    hResbos31->SetBinContent(4,lResbos31->GetBinContent(4)+lResbos31->GetBinContent(5));
    hResbos32->SetBinContent(4,lResbos32->GetBinContent(4)+lResbos32->GetBinContent(5));
    hResbos33->SetBinContent(4,lResbos33->GetBinContent(4)+lResbos33->GetBinContent(5));
    hResbos34->SetBinContent(4,lResbos34->GetBinContent(4)+lResbos34->GetBinContent(5));
    hResbos35->SetBinContent(4,lResbos35->GetBinContent(4)+lResbos35->GetBinContent(5));

  for(int i(5);i<=12;i++)
  {
    hResbos29->SetBinContent(i,lResbos29->GetBinContent(i+1));
    hResbos30->SetBinContent(i,lResbos30->GetBinContent(i+1));
    hResbos31->SetBinContent(i,lResbos31->GetBinContent(i+1));
    hResbos32->SetBinContent(i,lResbos32->GetBinContent(i+1));
    hResbos33->SetBinContent(i,lResbos33->GetBinContent(i+1));
    hResbos34->SetBinContent(i,lResbos34->GetBinContent(i+1));
    hResbos35->SetBinContent(i,lResbos35->GetBinContent(i+1));
  }

  for(int i(1);i<13;i++)
  {
    cout<<i<<"\t"<<hResbos29->GetBinContent(i)<<endl;
  }


 // TFile* f_output = new TFile("./FEWZ_"+BaseName+"_12Bin"".root","RECREATE");
  TFile* f_output = new TFile("./Resbos_"+BaseName+"_12Bin.root","RECREATE");

  hResbos29->Write();
  hResbos30->Write();
  hResbos31->Write();
  hResbos32->Write();
  hResbos33->Write();
  hResbos34->Write();
  hResbos35->Write();



  return 0;
}
