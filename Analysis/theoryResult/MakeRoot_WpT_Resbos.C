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

int MakeRoot_WpT_Resbos(TString BaseName)
{
  // cross-section and errors
  double Xsec[13] = {0.,};
  double Error[13] = {0.,};
  if(BaseName == "WpToMuNu")
  {
    Xsec[0] =   1440.91;  
    Xsec[1] =   725.499;
    Xsec[2] =   449.854;
    Xsec[3] =   358.457;
    Xsec[4] =   213.937;
    Xsec[5] =   219.716;
    Xsec[6] =   126.62 ;
    Xsec[7] =   126.39 ;
    Xsec[8] =   81.7365;
    Xsec[9] =   22.7662;
    Xsec[10] =  7.81976;
    Xsec[11] =  4.05342;
    Xsec[12] =  2.00304;

    Error[0] =   104.479 ; 
    Error[1] =   14.2278 ;
    Error[2] =   9.37335 ;
    Error[3] =   13.5594 ;
    Error[4] =   8.95398 ;
    Error[5] =   13.7904 ;
    Error[6] =   9.52116 ;
    Error[7] =   10.6327 ;
    Error[8] =   7.87614 ;
    Error[9] =   2.8959  ;
    Error[10]=   0.894395;
    Error[11]=   0.465657;
    Error[12]=   0.250487;
  }
  else if (BaseName == "WmToMuNu")
  {
    Xsec[0] =   907.199;  
    Xsec[1] =   471.034;
    Xsec[2] =   293.138;
    Xsec[3] =   234.647;
    Xsec[4] =   141.761;
    Xsec[5] =   146.899;
    Xsec[6] =   87.397 ;
    Xsec[7] =   88.5447;
    Xsec[8] =   59.2909;
    Xsec[9] =   15.5784;
    Xsec[10] =  5.14592;
    Xsec[11] =  2.52096;
    Xsec[12] =  1.12513;

    Error[0] =  64.3433   ; 
    Error[1] =  9.7948    ;
    Error[2] =  6.20416   ;
    Error[3] =  9.73717   ;
    Error[4] =  5.79605   ;
    Error[5] =  9.72729   ;
    Error[6] =  6.3833    ;
    Error[7] =  8.44543   ;
    Error[8] =  5.74321   ;
    Error[9] =  1.99148   ;
    Error[10]=  0.547323  ;
    Error[11]=  0.25781   ;
    Error[12]=  0.127705  ;
  }

  double Xsec_12bin[13] = {0,};
  double Error_12bin[13] = {0,};

  // Make 12 bin
  for(int i(0); i<3; i++)
  {
    Xsec_12bin[i] = Xsec[i];
    Error_12bin[i] = Error[i];
  }
  Xsec_12bin[3] = Xsec[3]+Xsec[4];
  Error_12bin[3] = sqrt(Error[3]**2 +Error[4]**2); // error : sqrt sum
  for(int i(4); i<13; i++)
  {
    Xsec_12bin[i] = Xsec[i+1];
    Error_12bin[i] = Error[i+1];
  }

  // 12Bin print and check
  for(int i(0); i<12; i++)
  {
    printf("Xsec_12 %d : %.8f, Error_12 : %.8f \n",i, Xsec_12bin[i], Error_12bin[i]);
  }

  // Calculate 12bin NormXsec and Errors
  double NormXsec_12bin[12];
  double NormError_12bin[12];
  ErrPropaNormXsec(Xsec_12bin, Error_12bin, NormXsec_12bin, NormError_12bin); 
  
  // Normalized 12Bin  print and check
  for(int i(0); i<12; i++)
  {
    printf("NormXsec_12 %d : %.8f, Error_12 : %.8f \n",i, NormXsec_12bin[i], NormError_12bin[i]);
  }

  // Make Normalized Differential 12 bin
  double NormDiffXsec_12bin[12];
  double NormDiffError_12bin[12];
  for(int i(0); i<12; i++)
  {
    NormDiffXsec_12bin[i] = NormXsec_12bin[i] / BinWidth[i] ;
    NormDiffError_12bin[i] = NormError_12bin[i] / BinWidth[i] ;
    printf("NormDiffXsec_12 %d : %.8f, Error_12 : %.8f \n",i, NormDiffXsec_12bin[i], NormDiffError_12bin[i]);
  }

  // Make histogram
  TH1D* Xsec_Resbos_13bin = new TH1D("Xsec_Resbos_13bin","Xsec_Resbos_13bin",13,0,600);
  TH1D* NormXsec_Resbos_12bin = new TH1D("NormXsec_Resbos_12bin","NormXsec_Resbos_12bin",12,0,600);
  TH1D* NormDiffXsec_Resbos_12bin = new TH1D("NormDiffXsec_Resbos_12bin","NormDiffXsec_Resbos_12bin",12,0,600);

  for(int i(0); i<13; i++)
  {
    Xsec_Resbos_13bin->SetBinContent(i+1,Xsec[i]);
    Xsec_Resbos_13bin->SetBinError(i+1,Error[i]);
  }

  for(int i(0); i<12; i++)
  {
    NormXsec_Resbos_12bin->SetBinContent(i+1,NormXsec_12bin[i]);
    NormXsec_Resbos_12bin->SetBinError(i+1,NormError_12bin[i]);
    
    NormDiffXsec_Resbos_12bin->SetBinContent(i+1,NormDiffXsec_12bin[i]);
    NormDiffXsec_Resbos_12bin->SetBinError(i+1,NormDiffError_12bin[i]);
  }

  TFile* fResbos = new TFile("./root/"+BaseName+"_Resbos.root","recreate");
  Xsec_Resbos_13bin->Write();
  NormXsec_Resbos_12bin->Write();
  NormDiffXsec_Resbos_12bin->Write();
  
  return 0;

}

