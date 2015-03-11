#include <stdio>
#include <iostream>

static const int NB = 12;
double ZptBins[19] = {0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,30.0,40.0,50.0,70.0,90.0,110.0,150.0,190.0,250.0,600.0};
double ZptBinWidth[18] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,10,10,10,20,20,20,40,40,60,350};
double Zpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
//double Zpt12BinWidth[12] = {7.5,5,5,12.5,10,10,20,20,40,40,60,350};
double Zpt12BinWidth[12] = {7.5-0,12.5-7.5,17.5-12.5,30-17.5,40-30,50-40,70-50,110-70,150-110,190-150,250-190,600-250};
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

int MakeRoot_ZpT_Resbos()
{
  // Normalized Differential cross-section and errors
  double NormDiffXsec[18] = {0.02957739, 0.05620224, 0.04961592, 0.04050349, 0.03325381, 0.02725463, 0.02253389, 0.01862596, 0.01237984, 0.00664593, 0.00393008, 0.00200304, 0.00085508, 0.00041032, 0.00017656, 0.00006332, 0.00002423, 0.00000115};
  double NormDiffError[18] = {0.00514001, 0.00707116, 0.00664392, 0.00600288, 0.00543919, 0.00492418, 0.00447746, 0.00407074, 0.00165936, 0.00121580, 0.00093494, 0.00047197, 0.00030837, 0.00021361, 0.00009908, 0.00005934, 0.00002997, 0.00000270};

  double NormXsec[18] ={0,};
  double NormError[18] ={0,};

  double NormXsec_12bin[12] = {0,};
  double NormError_12bin[12] = {0,};

  double NormDiffXsec_12bin[12] = {0,};
  double NormDiffError_12bin[12] = {0,};

  // Make NormDiff -> Norm
  for(int i(0); i<18; i++)
  {
    NormXsec[i] = NormDiffXsec[i] * ZptBinWidth[i];
    NormError[i] = NormDiffError[i] * ZptBinWidth[i];
  }
  
  // make 12 bin  
  NormXsec_12bin[0] = NormXsec[0] + NormXsec[1] + NormXsec[2];
  NormXsec_12bin[1] = NormXsec[3] + NormXsec[4];
  NormXsec_12bin[2] = NormXsec[5] + NormXsec[6];
  NormXsec_12bin[3] = NormXsec[7] + NormXsec[8];
  NormXsec_12bin[4] = NormXsec[9];
  NormXsec_12bin[5] = NormXsec[10];
  NormXsec_12bin[6] = NormXsec[11];
  NormXsec_12bin[7] = NormXsec[12] + NormXsec[13];
  NormXsec_12bin[8] = NormXsec[14];
  NormXsec_12bin[9] = NormXsec[15];
  NormXsec_12bin[10] = NormXsec[16];
  NormXsec_12bin[11] = NormXsec[17];

  NormError_12bin[0] =  sqrt(NormError[0]**2 + NormError[1]**2 + NormError[2]**2);
  NormError_12bin[1] =  sqrt(NormError[3]**2 + NormError[4]**2);
  NormError_12bin[2] =  sqrt(NormError[5]**2 + NormError[6]**2);
  NormError_12bin[3] =  sqrt(NormError[7]**2 + NormError[8]**2);
  NormError_12bin[4] =  NormError[9];
  NormError_12bin[5] =  NormError[10];
  NormError_12bin[6] =  NormError[11];
  NormError_12bin[7] =  sqrt(NormError[12]**2 + NormError[13]**2);
  NormError_12bin[8] =  NormError[14];
  NormError_12bin[9] = NormError[15];
  NormError_12bin[10] = NormError[16];
  NormError_12bin[11] = NormError[17];

  // NormXsec 12bin -> NormDiffXsec 12bin
  for(int i(0); i<12; i++)
  {
    NormDiffXsec_12bin[i] = NormXsec_12bin[i] / Zpt12BinWidth[i];
    NormDiffError_12bin[i] = NormError_12bin[i] / Zpt12BinWidth[i];
  }
  
  // Print Normalized Xsec and Error
  for(int i(0); i<12; i++)
  {
    printf("NormDiffXsec_12 %d : %.8f, Error : %.8f \n",i, NormDiffXsec_12bin[i], NormDiffError_12bin[i]);
  }

  //ErrPropaNormXsec(Xsec, TotalUncer, Norm_Xsec, Norm_TotalUncer); 

  // Make histogram
  TH1D* NormDiffXsec_Resbos = new TH1D("NormDiffXsec_Resbos","NormDiffXsec_Resbos",18,0,600);
  TH1D* NormDiffXsec_Resbos_12bin = new TH1D("NormDiffXsec_Resbos_12bin","NormDiffXsec_Resbos_12bin",12,0,600);

  for(int i(0); i<18; i++)
  {
    NormDiffXsec_Resbos->SetBinContent(i+1,NormDiffXsec[i]);
    NormDiffXsec_Resbos->SetBinError(i+1,NormDiffError[i]);
  }

  for(int i(0); i<12; i++)
  {
    NormDiffXsec_Resbos_12bin->SetBinContent(i+1,NormDiffXsec_12bin[i]);
    NormDiffXsec_Resbos_12bin->SetBinError(i+1,NormDiffError_12bin[i]);
  }

  TFile* fResbos = new TFile("./root/ZpT_Resbos.root","recreate");
  NormDiffXsec_Resbos->Write();
  NormDiffXsec_Resbos_12bin->Write();

  return 0;

}

