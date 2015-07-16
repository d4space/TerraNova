#include <stdio>
#include <iostream>

static const int NB = 12;
double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

int MakeRoot_WpTincl_Resbos()
{
  // cross-section and errors
  double Xsec[12] = {0.,};
    Xsec[0] = 2348.54  ;  
    Xsec[1] = 1196.82  ;
    Xsec[2] = 743.154  ;
    Xsec[3] = 948.342  ;
    Xsec[4] = 366.855  ;
    Xsec[5] = 214.112  ;
    Xsec[6] = 214.825  ;
    Xsec[7] = 141.191  ;
    Xsec[8] = 38.3433  ;
    Xsec[9] = 13.0015  ;
    Xsec[10] =6.58248  ;
    Xsec[11] =3.1293   ;
  
  double NormDiffError[12] = {0.,};
    NormDiffError[0] =  0.00231273 ; 
    NormDiffError[1] =  0.000777641;
    NormDiffError[2] =  0.000723605;
    NormDiffError[3] =  0.000354121;
    NormDiffError[4] =  0.000287119;
    NormDiffError[5] =  0.000202913;
    NormDiffError[6] =  0.000128583;
    NormDiffError[7] =  4.55758e-05;
    NormDiffError[8] =  1.72601e-05;
    NormDiffError[9] =  5.05825e-06;
    NormDiffError[10]=  1.64674e-06;
    NormDiffError[11]=  1.16567e-07;
  
  // Total Xsec 
  double TotalXsec = 0;
  for(int i(0); i<12; i++)
  {
    TotalXsec += Xsec[i];
    cout << "total Xsec : " << TotalXsec << "\t xsec : " << Xsec[i] << endl;
  }
  // Normalized differential cross-section
  double NormDiffXsec[12] = {0.,};
  for(int i(0); i<12; i++)
  {
    NormDiffXsec[i] = Xsec[i] / TotalXsec / BinWidth[i];
    printf("NormDiffXsec_12 %d : %.8f, Error_12 : %.8f \n",i, NormDiffXsec[i], NormDiffError[i]);
  }

  // Make histogram
  TH1D* NormDiffXsec_Resbos = new TH1D("NormDiffXsec_Resbos","NormDiffXsec_Resbos",12,0,600);

  for(int i(0); i<12; i++)
  {
    NormDiffXsec_Resbos->SetBinContent(i+1,NormDiffXsec[i]);
    NormDiffXsec_Resbos->SetBinError(i+1,NormDiffError[i]);
  }

  TFile* fResbos = new TFile("./root/WpTincl_Resbos.root","recreate");
  NormDiffXsec_Resbos->Write();
  
  return 0;

}

