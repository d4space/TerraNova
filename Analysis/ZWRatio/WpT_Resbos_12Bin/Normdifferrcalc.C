#include <iostream>

double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};
const int nBins =  13;
void Normdifferrcalc()
{
  // ResBos X sec value and errors, making TGraph .

  TFile* f_Resbos = new TFile("./ResultWInclToMuNu/Result_WInclToMuNu_Theory.root"); 

  char tmpName[30], tmpName_org[30];
  int Numb;
  
  TH1D* hResBos30_Xsec;
  hResBos30_Xsec = (TH1D*)f_Resbos->Get("hResbos30")->Clone("hResBos30_Xsec");
  double ResBos_TotalXsec = hResBos30_Xsec->Integral(); //Resbos Total Xsec center grid
  cout << " ResBos_TotalXsec : " << ResBos_TotalXsec << endl;
  
  // Calculate NormDiff Xsec Resbos
  double ResBos_NormDiffXsec[12];
  for(int i(0); i<12; i++)
  {
    ResBos_NormDiffXsec[i] = hResBos30_Xsec->GetBinContent(i+1) / BinWidth[i] / ResBos_TotalXsec;
    cout << "ResBos NormDiffXsec : " << ResBos_NormDiffXsec[i] << endl;
  }

  // Calculate NormDiff Xsec Resbos other grids
  TH1D* AllResbos[7];
  double AllResbos_TotalXsec[7];
  double AllResbos_NormDiffXsec[7][13];
  for( int i(0);i<7;i++)
  {
    Numb = 29+i;
    sprintf(tmpName_org,"hResbos%d",Numb);
    sprintf(tmpName,"AllResbos_%d",i);
    AllResbos[i] = (TH1D*)f_Resbos->Get(tmpName_org)->Clone(tmpName);
    AllResbos_TotalXsec[i] = AllResbos[i]->Integral(); // Total Xsec for Resbos grid
    cout << "AllResbos TotalXsec : " << AllResbos_TotalXsec[i] << endl;
    for(int j(0); j<12; j++)
    {
      AllResbos_NormDiffXsec[i][j]=AllResbos[i]->GetBinContent(j+1) / BinWidth[j]/ AllResbos_TotalXsec[i];
      cout << "Resbos NormDiff : " << AllResbos[i]->GetBinContent(j+1) / BinWidth[j]/ AllResbos_TotalXsec[i]  << endl;
    }
  }

  //normdifferror calc
  Double_t Resb_errMax[nBins-1];
  Double_t Resb_errMin[nBins-1];
  Double_t NormResb_errMax[nBins-1];
  Double_t NormResb_errMin[nBins-1];
  double tmpVal,tmpDiff;

  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    double nomVal  = ResBos_NormDiffXsec[ipt];
    Resb_errMax[ipt] = -99999;
    Resb_errMin[ipt] = 990009;
    for (int i(0);i<7;i++)
    {
      tmpVal  = AllResbos_NormDiffXsec[i][ipt];
      tmpDiff = tmpVal - nomVal;
      if( tmpDiff > Resb_errMax[ipt] ) Resb_errMax[ipt] = tmpDiff;
      if( tmpDiff < Resb_errMin[ipt] ) Resb_errMin[ipt] = tmpDiff;
    }
    
    if (Resb_errMax[ipt] < 0) Resb_errMax[ipt] = 0.;
    if (Resb_errMin[ipt] > 0) Resb_errMin[ipt] = 0.;
    if (Resb_errMin[ipt] < 0) Resb_errMin[ipt] = -Resb_errMin[ipt];
    Resb_errMax[ipt] = Resb_errMax[ipt];
    Resb_errMin[ipt] = Resb_errMin[ipt];
  cout << "Resbos Xsec : " << ResBos_NormDiffXsec[ipt] << "\t error+ : " << Resb_errMax[ipt] << "\t error - : " << Resb_errMin[ipt] << endl; 
  }


}
