#include <stdio>
#include <iostream>

const double Zpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
const double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

static const int NB = 12;

void NormToyErr(double *Xsec, double *Err, double *Norm_Err);
void ErrPropaNormXsec(double *nn, double *dn, double ff[NB], double df[NB]);

int MakeRoot_ZpT_FEWZ()
{
  double Xsec[12];
  double Xsec_up[12];
  double Xsec_down[12];

  double Total_Xsec = 0;
  double Total_Xsec_up = 0;
  double Total_Xsec_down = 0;

  double StatErr[12];

  double TotalUncer[12];

  // Normalized cross-section and errors
  double Norm_Xsec[12];
  double Norm_Xsec_up[12];
  double Norm_Xsec_down[12];

  double Norm_StatErr[12];
  double Norm_PDFErrP[12];
  double Norm_PDFErrM[12];
  double Norm_PDFErr[12];

  double Norm_TotalUncer[12];

  // Normalized Differential Xsec and errors
  double NormDiff_Xsec[12];
  double NormDiff_Xsec_up[12];
  double NormDiff_Xsec_down[12];

  double NormDiff_StatErr[12];
  double NormDiff_PDFErr[12];
  double NormDiff_TotalUncer[12];

  Xsec[0] =  138.274;  
  Xsec[1] =  103.757;
  Xsec[2] =  59.4709;
  Xsec[3] =  68.4585;
  Xsec[4] =  26.5977;
  Xsec[5] =  15.9503;
  Xsec[6] =  16.6003;
  Xsec[7] =  10.8541;  
  Xsec[8] =  2.98891;
  Xsec[9] =  1.03651;
  Xsec[10] = 0.551997; 
  Xsec[11] = 0.281856;

  Xsec_up[0] =  142.959;  
  Xsec_up[1] =    104.9;
  Xsec_up[2] =  58.6561;
  Xsec_up[3] =  67.1108;
  Xsec_up[4] =  25.7478;
  Xsec_up[5] =  15.2494;
  Xsec_up[6] =  15.8542;
  Xsec_up[7] =  10.3382;  
  Xsec_up[8] =  2.79667;
  Xsec_up[9] =  0.962229;
  Xsec_up[10] = 0.495937; 
  Xsec_up[11] = 0.257324;

  Xsec_down[0] =  135.624;  
  Xsec_down[1] =  102.489;
  Xsec_down[2] =  59.4373;
  Xsec_down[3] =  69.8794;
  Xsec_down[4] =  27.6927;
  Xsec_down[5] =  16.6091;
  Xsec_down[6] =  17.4425;
  Xsec_down[7] =  11.5709;  
  Xsec_down[8] =  3.15182;
  Xsec_down[9] =  1.1378;
  Xsec_down[10] = 0.583274; 
  Xsec_down[11] = 0.314832;

  StatErr[0]	=  0.991013;  
  StatErr[1]	=  0.347511;
  StatErr[2]	=  0.197072;
  StatErr[3]	=  0.139355;
  StatErr[4]	= 0.0763465;
  StatErr[5]	= 0.0709948;
  StatErr[6]	= 0.0728004;
  StatErr[7]	= 0.0373067;
  StatErr[8]	= 0.0217453;
  StatErr[9]	= 0.0128423;
  StatErr[10]	= 0.0113465; 
  StatErr[11]	= 0.0180237;
 
  // Norm PDF Error in % unit
  Norm_PDFErrP[0]	=4.23  ;   
  Norm_PDFErrP[1]	=1.29  ;  
  Norm_PDFErrP[2]	=1.40  ;
  Norm_PDFErrP[3]	=1.76  ;
  Norm_PDFErrP[4]	=2.20  ;
  Norm_PDFErrP[5]	=2.53  ;
  Norm_PDFErrP[6]	=2.93  ;
  Norm_PDFErrP[7]	=3.66  ;
  Norm_PDFErrP[8]	=4.51  ;
  Norm_PDFErrP[9]	=5.06  ;
  Norm_PDFErrP[10]	= 5.65 ;
  Norm_PDFErrP[11]	= 6.50 ;

  Norm_PDFErrM[0]	=3.51 ; 
  Norm_PDFErrM[1]	=1.11 ; 
  Norm_PDFErrM[2]	=1.58 ; 
  Norm_PDFErrM[3]	=2.08 ; 
  Norm_PDFErrM[4]	=2.74 ; 
  Norm_PDFErrM[5]	=3.21 ;
  Norm_PDFErrM[6]	=3.81 ;
  Norm_PDFErrM[7]	=4.92 ;
  Norm_PDFErrM[8]	=6.51 ;
  Norm_PDFErrM[9]	=7.70 ;
  Norm_PDFErrM[10]	= 9.27;
  Norm_PDFErrM[11]	= 11.07;

  // Error Calculation 
  for(int i(0); i<12; i++)
  {
    Norm_PDFErr[i] = TMath::Max(Norm_PDFErrP[i],Norm_PDFErrM[i]); // choose bigger pdf error
  }

  // Calculate Nomianl, ScaleUp,Down case Total Xsec
  for(int i(0); i<12; i++)
  {
    Total_Xsec += Xsec[i];
    Total_Xsec_up += Xsec_up[i];
    Total_Xsec_down += Xsec_down[i];
  }
  
  // Print out Total Xsecs
  //printf("Total Xsec ScaleUp : %.4f ScaleDown : %.4f PDFp : %.4f PDFm : %.4f \n",Total_Xsec_up,Total_Xsec_down,Total_Xsec_PDFp,Total_Xsec_PDFm);

  for(int i(0); i<12; i++)
  {
    Norm_Xsec[i] = Xsec[i] / Total_Xsec;
    Norm_Xsec_up[i] = Xsec_up[i] / Total_Xsec_up; // Calculation of Scale up case  Normalized Xsec 
    Norm_Xsec_down[i] = Xsec_down[i] / Total_Xsec_down; // Calculation of Scale down case Normalized Xsec

    //Norm_StatErr[i] = StatErr[i] / Total_Xsec;
  }

  // Cross-section and Statistical error Normalization by Toy method
  NormToyErr(Xsec,StatErr,Norm_StatErr);

  // Print out All Normalized Xsec
  //for(int i(0); i<12; i++)
  //{
    //cout << "Norm_Xsec : " << Norm_Xsec[i] << "\t Norm_StatPDFErr : " << Norm_StatPDFErr[i] << endl; 
  //}

  // Calculate Normalized Xsec and Total uncertainty by using 1st order differential formula ??? 
  //ErrPropaNormXsec(Xsec, TotalUncer, Norm_Xsec, Norm_TotalUncer); 
  
  // Calculate Normalized Differential Xsec and Total Uncertainty
  for(int i(0); i<12; i++)
  {
    NormDiff_Xsec[i] = Norm_Xsec[i] / BinWidth[i] ; 
    NormDiff_Xsec_up[i] = Norm_Xsec_up[i] / BinWidth[i];
    NormDiff_Xsec_down[i] = Norm_Xsec_down[i] / BinWidth[i];
    
    NormDiff_StatErr[i] = Norm_StatErr[i] / BinWidth[i]; 
    NormDiff_PDFErr[i] = Norm_PDFErr[i]*NormDiff_Xsec[i] *0.01;
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatErr[i]**2 + NormDiff_PDFErr[i]**2) ; 
    NormDiff_TotalUncer[i] = NormDiff_StatErr[i] ; 
    
    printf("NormDiff_Xsec : %.8f \t +- %.8f \t  %.2f %\n ",NormDiff_Xsec[i],NormDiff_TotalUncer[i], NormDiff_TotalUncer[i]/NormDiff_Xsec[i]*100);
  }
  
  // Save in histogram
  TH1D* hxsec = new TH1D("hxsec","hxsec",12,Zpt12Bins);
  TH1D* hxsec_up = new TH1D("hxsec_up","hxsec Scale up",12,Zpt12Bins);
  TH1D* hxsec_down = new TH1D("hxsec_down","hxsec Scale down",12,Zpt12Bins);
  
  TH1D* hxsec_NormDiff = new TH1D("hxsec_NormDiff","hxsec_NormDiff",12,Zpt12Bins);
  TH1D* hxsec_NormDiff_up = new TH1D("hxsec_NormDiff_up","hxsec_NormDiff_up",12,Zpt12Bins);
  TH1D* hxsec_NormDiff_down = new TH1D("hxsec_NormDiff_down","hxsec_NormDiff_down",12,Zpt12Bins);

  TH1D* hStatErr = new TH1D("StatErr","Stat Error", 12,Zpt12Bins);
  TH1D* hPDFErr = new TH1D("PDFErr","PDF Error", 12,Zpt12Bins);

  for(int i(0); i<12; i++)
  {
    hxsec->SetBinContent(i+1,Xsec[i]);
    hxsec_up->SetBinContent(i+1, Xsec_up[i]);
    hxsec_down->SetBinContent(i+1, Xsec_down[i]);

    hxsec_NormDiff->SetBinContent(i+1,NormDiff_Xsec[i]);
    hxsec_NormDiff->SetBinError(i+1,NormDiff_TotalUncer[i]);
    
    hxsec_NormDiff_up->SetBinContent(i+1,NormDiff_Xsec_up[i]);
    hxsec_NormDiff_down->SetBinContent(i+1,NormDiff_Xsec_down[i]);
  }

  TFile* f_out = new TFile("./root/ZToMuMu_FEWZ.root","recreate");
  
  hxsec->Write();
  hxsec_up->Write();
  hxsec_down->Write();
  hxsec_NormDiff->Write();
  hxsec_NormDiff_up->Write();
  hxsec_NormDiff_down->Write();

  return 0;
}

void NormToyErr(double *_Xsec, double *_Err, double *_Norm_Err)
{
  int Ntoy = 100000;

  double mx2[12] = {0.,};
  double m2x[12] = {0.,};
  double temp[12] = {0.,};

  // Check _Xsec and _StatErr argument is correct 
//  for(int j=0; j<12; ++j)
//  {
    //cout << Form("_Xsec : %.4f \t _StatErr : %.4f",_Xsec[j],_StatErr[j]) << endl;
//  }

  // Here Calculate Normalized Stat Error by using toy
  for(int e=0; e<Ntoy; ++e)
  {
    for(int i=0; i<12; ++i)
    {
      	temp[i]=gRandom->Gaus(_Xsec[i], _Err[i]);
    }
    double sum=0;
    for(int i=0; i<12; ++i) sum+=temp[i];
    for(int i=0; i<12; ++i)
    {
      temp[i]/=sum;
      m2x[i]+=temp[i];
      mx2[i]+=temp[i]*temp[i];
    }
  }

  double rms[12];
  for(int i=0; i<12; ++i)
  {
    m2x[i]/=Ntoy;
    mx2[i]/=Ntoy;
    rms[i] = sqrt(mx2[i]-m2x[i]*m2x[i]);
  }
 
  // return value
  for(int i=0; i<12; ++i)
  {
    _Norm_Err[i] = rms[i]; // return rms to Norm_Err
  }
  for(int i=0; i<12; ++i)
  {
    //cout << Form("rms = %.8f",rms[i]) << endl;
  }
}

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

