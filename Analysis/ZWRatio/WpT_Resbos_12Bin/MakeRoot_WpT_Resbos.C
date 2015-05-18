#include <stdio>
#include <iostream>

static const int NB = 12;
double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

int MakeRoot_WpT_Resbos(TString BaseName)
{
  TFile* fResbos = new TFile("../../RstResbos_12Bin/Resbos_"+BaseName+"_12Bin.root");

  TH1D* hGrid[7];
  TH1D* hNormDiffGrid[7];
  char gridname[50];
  double TotalXsec[7] = {0.,};
 
  // NormDiff xsec calculation
  for(int i(0);i<7;i++)
  {
    sprintf(gridname,"hResbos%d",i+29);
    hGrid[i] = (TH1D*)fResbos->Get(gridname); // read xsec of each grid
    hNormDiffGrid[i] = (TH1D*)hGrid[i]->Clone(); // clone hGrid
    TotalXsec[i] = hGrid[i]->Integral(); // calculate total xsec of each grid
    hNormDiffGrid[i]->Scale(1.0/TotalXsec[i]); // calculate Normalized Xsec (divide total xsec)

    // divide bin width 
    for(int j(1);j<13;j++)
    {
      hNormDiffGrid[i]->SetBinContent(j,hNormDiffGrid[i]->GetBinContent(j)/BinWidth[j-1]);
    }
  }


  TFile* fOut = new TFile("./root/"+BaseName+"_Resbos.root","recreate");
  for(int i(0); i<7; i++)
  {
    hNormDiffGrid[i]->Write();
  }

  return 0;

}

