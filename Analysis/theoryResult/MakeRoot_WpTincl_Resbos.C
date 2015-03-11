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

int MakeRoot_WpTincl_Resbos()
{
  // cross-section and errors
  double Xsec[13] = {0.,};
  double Error[13] = {0.,};
    Xsec[0] =  2348.10833740 ;  
    Xsec[1] =  1196.53350830 ;
    Xsec[2] =  742.99261475  ;
    Xsec[3] =  593.10450745  ;
    Xsec[4] =  355.69859314  ;
    Xsec[5] =  366.61437988  ;
    Xsec[6] =  214.01686859  ;
    Xsec[7] =  214.93440247  ;
    Xsec[8] =  141.02732468  ;
    Xsec[9] =  38.34463024   ;
    Xsec[10] = 12.96568298   ;
    Xsec[11] = 6.57437801    ;
    Xsec[12] = 3.12816536    ;

    Error[0] =  167.50604248 ; 
    Error[1] =  23.35012817  ;
    Error[2] =  15.57751465 ;
    Error[3] =  23.29661560 ;
    Error[4] =  14.75003052 ;
    Error[5] =  23.51766968 ;
    Error[6] =  15.90446472 ;
    Error[7] =  19.07817841 ;
    Error[8] =  13.61935043 ;
    Error[9] =  4.88738346  ;
    Error[10]=  1.44171762  ;
    Error[11]=  0.72346711  ;
    Error[12]=  0.27354395  ;

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

  TFile* fResbos = new TFile("./root/WpTincl_Resbos.root","recreate");
  Xsec_Resbos_13bin->Write();
  NormXsec_Resbos_12bin->Write();
  NormDiffXsec_Resbos_12bin->Write();
  
  return 0;

}

