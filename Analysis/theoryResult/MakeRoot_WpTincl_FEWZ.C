#include <stdio>
#include <iostream>

static const int NB = 12;
double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

void ErrPropaNormXsec(double *nn, double *dn, double ff[NB], double df[NB]){

  double ff[NB];
  //double df[NB];
  double dp[NB];

  ///Here is error propagation
  double NN=0;
  for(int i=0; i<NB; ++i) NN+=nn[i];
  for(int i=0; i<NB; ++i) {
    ff[i]=nn[i]/NN;
    dp[i]=dn[i]/NN;
    df[i]=0;
  }
  for(int i=0; i<NB; ++i){
    df[i] = pow( (1-ff[i])*dp[i], 2);
    for(int j=0; j<NB; ++j){
      if(j!=i) df[i]+=ff[i]*ff[i]*dp[j]*dp[j];
    }
    df[i]=sqrt(df[i]);
  }

  for(int i=0; i<NB; ++i){
    //cout << Form("%.8f  prop=%.8f", ff[i]/BinWidth[i], df[i]/BinWidth[i]) << endl;
    cout << Form("%.8f  prop=%.8f", ff[i], df[i]) << endl; // Print Normalized Xsec and TotalUncer.
  }
}

int MakeRoot_WpTincl_FEWZ()
{
  // cross-section and errors
  double Xsec[12] = {0.,};
  double Xsec_Wp[12] = {0.,};
  double Xsec_Wm[12] = {0.,};
  double Error[12] = {0.,};
  double Error_Wp[12] = {0.,};
  double Error_Wm[12] = {0.,};

  // W+ cross-section and total uncertainty
  Xsec_Wp[0] =  1454.93;  
  Xsec_Wp[1] =  811.191;
  Xsec_Wp[2] =  462.192;
  Xsec_Wp[3] =  542.741;
  Xsec_Wp[4] =  212.959;
  Xsec_Wp[5] =  121.442;
  Xsec_Wp[6] =  120.148;
  Xsec_Wp[7] =  75.4574;  
  Xsec_Wp[8] =  19.7655;
  Xsec_Wp[9] =  6.44498;
  Xsec_Wp[10] = 3.25329; 
  Xsec_Wp[11] = 1.58361;

  Error_Wp[0] = 135.55728530  ;
  Error_Wp[1] =  52.26503953  ;
  Error_Wp[2] =  29.58199419  ;
  Error_Wp[3] =  35.72014109  ;
  Error_Wp[4] =  14.88893591  ;
  Error_Wp[5] =   8.71445118  ;
  Error_Wp[6] =   8.92524355  ;
  Error_Wp[7] =   5.72036044  ;
  Error_Wp[8] =   1.57958237  ;
  Error_Wp[9] =   0.53131101  ;
  Error_Wp[10] =  0.50392133  ;
  Error_Wp[11] =  0.24155716  ;

  // W- cross-section and Total Uncertainty
  Xsec_Wm[0] = 909.233 ;  
  Xsec_Wm[1] = 537.571 ;
  Xsec_Wm[2] = 308.087 ;
  Xsec_Wm[3] = 366.409 ;
  Xsec_Wm[4] =  146.48 ;
  Xsec_Wm[5] =  84.5731;
  Xsec_Wm[6] =  85.4917;
  Xsec_Wm[7] =  54.6661;  
  Xsec_Wm[8] =  14.1036;
  Xsec_Wm[9] =  4.39228;
  Xsec_Wm[10] = 2.15351; 
  Xsec_Wm[11] = 0.940221;

  Error_Wm[0] =  96.51797048;
  Error_Wm[1] =  36.63957951;
  Error_Wm[2] =  20.00443721;
  Error_Wm[3] =  24.13938768;
  Error_Wm[4] =  10.16958671;
  Error_Wm[5] =   5.79590395;
  Error_Wm[6] =   6.11968783;
  Error_Wm[7] =   4.38428571;
  Error_Wm[8] =   1.33215735;
  Error_Wm[9] =   0.37390322;
  Error_Wm[10] =  0.21785016;
  Error_Wm[11] =  0.12688078;

  // Make Inclusive
  for(int i(0); i<12; i++)
  {
    Xsec[i] = Xsec_Wp[i] + Xsec_Wm[i];
    Error[i] = sqrt(Error_Wp[i]*Error_Wp[i] + Error_Wm[i]*Error_Wm[i]);
    printf("Xsec : %f \t +- %f, \t %f %\n", Xsec[i],Error[i],Error[i]/Xsec[i]*100);
  }
  
  // Calculate 12bin NormXsec and Errors
  double NormXsec[12];
  double NormError[12];
  ErrPropaNormXsec(Xsec, Error, NormXsec, NormError); 

  // Normalized 12Bin  print and check
  for(int i(0); i<12; i++)
  {
    printf("NormXsec %d : %.8f, Error : %.8f, \t %f % \n",i, NormXsec[i], NormError[i], NormError[i]/NormXsec[i]*100);
  }

  // Make Normalized Differential 12 bin
  double NormDiffXsec[12];
  double NormDiffError[12];
  for(int i(0); i<12; i++)
  {
    NormDiffXsec[i] = NormXsec[i] / BinWidth[i] ;
    NormDiffError[i] = NormError[i] / BinWidth[i] ;
    printf("NormDiffXsec %d : %.8f, Error : %.8f \n",i, NormDiffXsec[i], NormDiffError[i]);
  }

  // Make histogram
  TH1D* Xsec_FEWZ_13bin = new TH1D("Xsec_FEWZ_13bin","Xsec_FEWZ_13bin",13,0,600);
  TH1D* NormXsec_FEWZ = new TH1D("NormXsec_FEWZ","NormXsec_FEWZ",12,0,600);
  TH1D* NormDiffXsec_FEWZ = new TH1D("NormDiffXsec_FEWZ","NormDiffXsec_FEWZ",12,0,600);

  for(int i(0); i<13; i++)
  {
    Xsec_FEWZ_13bin->SetBinContent(i+1,Xsec[i]);
    Xsec_FEWZ_13bin->SetBinError(i+1,Error[i]);
  }

  for(int i(0); i<12; i++)
  {
    NormXsec_FEWZ->SetBinContent(i+1,NormXsec[i]);
    NormXsec_FEWZ->SetBinError(i+1,NormError[i]);

    NormDiffXsec_FEWZ->SetBinContent(i+1,NormDiffXsec[i]);
    NormDiffXsec_FEWZ->SetBinError(i+1,NormDiffError[i]);
  }

  TFile* fFEWZ = new TFile("./root/WinclToMuNu_FEWZ.root","recreate");
  Xsec_FEWZ_13bin->Write();
  NormXsec_FEWZ->Write();
  NormDiffXsec_FEWZ->Write();

  return 0;

}

