#include <stdio>
#include <iostream>

#define Ntoy 50000

double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

static const int NB = 12;

void NormToyErr(double *Xsec, double *Err, double *Norm_Xsec, double *Norm_Err);
void ErrPropaNormXsec(double *nn, double *dn, double ff[NB], double df[NB]);

int MakeRoot_WpT_FEWZ(TString BaseName)
{
  double Xsec[12];
  double Xsec_up[12];
  double Xsec_down[12];
  double Xsec_PDFp[12];
  double Xsec_PDFm[12];
  
  double Total_Xsec = 0;
  double Total_Xsec_up = 0;
  double Total_Xsec_down = 0;
  double Total_Xsec_PDFp = 0;
  double Total_Xsec_PDFm = 0;
  
  double StatErr[12];
  
  double PDFErrP[12];
  double PDFErrM[12];
  double PDFErr[12];
  
  double StatPDFErr[12];

  double Norm_PDFErr[12];

  double ScaleErrP[12];
  double ScaleErrM[12];
  double ScaleErr[12];
  double ScaleErr_percent[12];

  double Norm_ScaleErrP[12];
  double Norm_ScaleErrM[12];
  double Norm_ScaleErr[12];
  
  double Norm_StatPDFErr[12];
  
  double TotalUncer[12];
 
  // Normalized cross-section and errors
  double Norm_Xsec[12];
  double Norm_Xsec_up[12];
  double Norm_Xsec_down[12];
  
  double Norm_StatErr[12];
  double Norm_TotalUncer[12];

  // Normalized Differential Xsec and errors
  double NormDiff_Xsec[12];
  double NormDiff_Xsec_up[12];
  double NormDiff_Xsec_down[12];

  double NormDiff_StatErr[12];
  double NormDiff_StatPDFErr[12];
  double NormDiff_ScaleErr[12];
  double NormDiff_TotalUncer[12];

  if(BaseName == "WpToMuNu")
  {
    Xsec[0] =  1454.93;  
    Xsec[1] =  811.191;
    Xsec[2] =  462.192;
    Xsec[3] =  542.741;
    Xsec[4] =  212.959;
    Xsec[5] =  121.442;
    Xsec[6] =  120.148;
    Xsec[7] =  75.4574;  
    Xsec[8] =  19.7655;
    Xsec[9] =  6.44498;
    Xsec[10] = 3.25329; 
    Xsec[11] = 1.58361;

    Xsec_up[0] =  1524.37;  
    Xsec_up[1] =  810.021;
    Xsec_up[2] =    454.4;
    Xsec_up[3] =  526.421;
    Xsec_up[4] =  204.797;
    Xsec_up[5] =  116.005;
    Xsec_up[6] =  114.095;
    Xsec_up[7] =  70.9633;  
    Xsec_up[8] =  18.3452;
    Xsec_up[9] =  6.01652;
    Xsec_up[10] = 2.95692; 
    Xsec_up[11] = 1.36274;

    Xsec_down[0] =  1388.61;  
    Xsec_down[1] =  810.397;
    Xsec_down[2] =  471.184;
    Xsec_down[3] =  560.601;
    Xsec_down[4] =  222.985;
    Xsec_down[5] =  128.017;
    Xsec_down[6] =  127.597;
    Xsec_down[7] =  80.7096;  
    Xsec_down[8] =  21.2982;
    Xsec_down[9] =  6.95912;
    Xsec_down[10] = 3.74546; 
    Xsec_down[11] = 1.73463;

    StatErr[0]	=   3.09186;  
    StatErr[1]	=    1.1629;
    StatErr[2]	=   0.65592;
    StatErr[3]	=  0.496853;
    StatErr[4]	=  0.340195;
    StatErr[5]	=  0.322259;
    StatErr[6]	=  0.243238;
    StatErr[7]	=  0.199159;
    StatErr[8]	=  0.116363;
    StatErr[9]	= 0.0786948;
    StatErr[10]	= 0.0698109; 
    StatErr[11]	= 0.0656297;

    PDFErrP[0]	=   116.38 ;   
    PDFErrP[1]	=   52.239;  
    PDFErrP[2]	=  28.1746;
    PDFErrP[3]	=  30.9306;
    PDFErrP[4]	=   11.002;
    PDFErrP[5]	=  5.71027;
    PDFErrP[6]	=  4.91052;
    PDFErrP[7]	=  2.25771;
    PDFErrP[8]	= 0.363828;
    PDFErrP[9]	= 0.093028;
    PDFErrP[10]	= 0.0627364;
    PDFErrP[11]	= 0.0451354;

    PDFErrM[0]	=  80.791; 
    PDFErrM[1]	= 29.0635; 
    PDFErrM[2]	=  14.896; 
    PDFErrM[3]	= 15.5722; 
    PDFErrM[4]	= 5.17971; 
    PDFErrM[5]	= 2.54266;
    PDFErrM[6]	= 2.08545;
    PDFErrM[7]	= 0.96614;
    PDFErrM[8]	= 0.249237;
    PDFErrM[9]	= 0.108437;
    PDFErrM[10]	= 0.0826549;
    PDFErrM[11]	= 0.0725193;
    
    // Norm_PDFErr made by toy - Hammid
    Norm_PDFErr[0]	= 5.2623; 
    Norm_PDFErr[1]	= 6.0062;
    Norm_PDFErr[2]	= 6.3814;
    Norm_PDFErr[3]	= 5.9908;
    Norm_PDFErr[4]	= 6.0217;
    Norm_PDFErr[5]	= 5.7650;
    Norm_PDFErr[6]	= 5.2885;
    Norm_PDFErr[7]	= 4.5724;
    Norm_PDFErr[8]	= 3.9804;
    Norm_PDFErr[9]	= 3.9192;
    Norm_PDFErr[10]	= 4.3573;
    Norm_PDFErr[11]	= 5.7849;

  }
  
  else if(BaseName == "WmToMuNu")
  {
    Xsec[0] = 909.233 ;  
    Xsec[1] = 537.571 ;
    Xsec[2] = 308.087 ;
    Xsec[3] = 366.409 ;
    Xsec[4] =  146.48 ;
    Xsec[5] =  84.5731;
    Xsec[6] =  85.4917;
    Xsec[7] =  54.6661;  
    Xsec[8] =  14.1036;
    Xsec[9] =  4.39228;
    Xsec[10] = 2.15351; 
    Xsec[11] = 0.940221;

    Xsec_up[0] = 953.265;  
    Xsec_up[1] = 538.044;
    Xsec_up[2] = 303.618;
    Xsec_up[3] = 356.334;
    Xsec_up[4] = 140.91;
    Xsec_up[5] = 81.6267;
    Xsec_up[6] = 80.8662;
    Xsec_up[7] = 51.0988;  
    Xsec_up[8] = 13.2132;
    Xsec_up[9] = 4.04412;
    Xsec_up[10] =1.95992; 
    Xsec_up[11] =0.835329;

    Xsec_down[0] = 862.83;  
    Xsec_down[1] = 535.113;
    Xsec_down[2] = 312.711;
    Xsec_down[3] = 378.113;
    Xsec_down[4] = 153.424; 
    Xsec_down[5] = 88.9128;
    Xsec_down[6] = 90.6541;
    Xsec_down[7] = 58.7801;  
    Xsec_down[8] = 15.4094;
    Xsec_down[9] = 4.72442;
    Xsec_down[10] =2.33955; 
    Xsec_down[11] =0.996848;

    StatErr[0]	=  1.95734;  
    StatErr[1]	= 0.707038;
    StatErr[2]	= 0.398771;
    StatErr[3]	= 0.300166;
    StatErr[4]	= 0.208696;
    StatErr[5]	= 0.233793;
    StatErr[6]	= 0.101461;
    StatErr[7]	= 0.123304;
    StatErr[8]	= 0.0665512;
    StatErr[9]	= 0.0516952;
    StatErr[10]	= 0.0438842; 
    StatErr[11]	= 0.0325658;

    PDFErrP[0]	= 84.6088;   
    PDFErrP[1]	= 36.5502;  
    PDFErrP[2]	= 19.4586;
    PDFErrP[3]	= 21.1101;
    PDFErrP[4]	= 7.42683;
    PDFErrP[5]	= 3.83469;
    PDFErrP[6]	= 3.2848;
    PDFErrP[7]	= 1.51055;
    PDFErrP[8]	= 0.255148;
    PDFErrP[9]	= 0.0865463;
    PDFErrP[10]	= 0.0520098;
    PDFErrP[11]	= 0.0323983;

    PDFErrM[0]	= 58.5766; 
    PDFErrM[1]	= 19.9108; 
    PDFErrM[2]	= 10.1198; 
    PDFErrM[3]	= 10.4207; 
    PDFErrM[4]	= 3.49216; 
    PDFErrM[5]	= 1.73801;
    PDFErrM[6]	= 1.48737;
    PDFErrM[7]	= 0.790125;
    PDFErrM[8]	= 0.244521;
    PDFErrM[9]	= 0.126158;
    PDFErrM[10]	= 0.089754;
    PDFErrM[11]	= 0.0635285;

    // Norm_PDFErr made by toy - Hammid
    Norm_PDFErr[0]	= 6.2259; 
    Norm_PDFErr[1]	= 6.4107;
    Norm_PDFErr[2]	= 6.7206;
    Norm_PDFErr[3]	= 6.2179;
    Norm_PDFErr[4]	= 6.1473;
    Norm_PDFErr[5]	= 5.8499;
    Norm_PDFErr[6]	= 5.3486;
    Norm_PDFErr[7]	= 4.7072;
    Norm_PDFErr[8]	= 4.2717;
    Norm_PDFErr[9]	= 4.8340;
    Norm_PDFErr[10]	= 5.6939;
    Norm_PDFErr[11]	= 7.7928;
  
  }

  // Error Calculation 
  for(int i(0); i<12; i++)
  {
    PDFErr[i] = TMath::Max(PDFErrP[i],PDFErrM[i]); // pdf error

    //ScaleErrP[i] = TMath::Abs(Xsec[i] - Xsec_up[i]);
    //ScaleErrM[i] = TMath::Abs(Xsec[i] - Xsec_down[i]);
    //ScaleErr[i] = TMath::Max(ScaleErrP[i],ScaleErrM[i]); // scale error
    
    StatPDFErr[i] = sqrt(StatErr[i]*StatErr[i] + PDFErr[i]*PDFErr[i]); // will be toy normalization
    
    //printf("Xsec : %f \t Stat : %.4f \t TotalUncer : %.8f, \t %.4f %\n",Xsec[i],StatErr[i],TotalUncer[i], TotalUncer[i]/Xsec[i]*100 );
  }
  
  // Calculate ScaleUp,Down case Total Xsec
  for(int i(0); i<12; i++)
  {
    Total_Xsec_up += Xsec_up[i];
    Total_Xsec_down += Xsec_down[i];
  }

  // Print out Total Xsecs
//  printf("Total Xsec ScaleUp : %.4f ScaleDown : %.4f PDFp : %.4f PDFm : %.4f \n",Total_Xsec_up,Total_Xsec_down,Total_Xsec_PDFp,Total_Xsec_PDFm);
  
  for(int i(0); i<12; i++)
  {
    Norm_Xsec_up[i] = Xsec_up[i] / Total_Xsec_up; // Calculation of Scale up case  Normalized Xsec 
    Norm_Xsec_down[i] = Xsec_down[i] / Total_Xsec_down; // Calculation of Scale down case Normalized Xsec
  }

  // Cross-section and Statistical error Normalization by Toy method
  //NormToyErr(Xsec,StatErr,Norm_Xsec,Norm_StatErr);
  NormToyErr(Xsec,StatPDFErr,Norm_Xsec,Norm_StatPDFErr);

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
    
    //NormDiff_StatErr[i] = Norm_StatErr[i] / BinWidth[i] ; 
    //NormDiff_ScaleErr[i] = TMath::Max(fabs(NormDiff_Xsec_up[i] - NormDiff_Xsec[i]),fabs(NormDiff_Xsec_down[i] - NormDiff_Xsec[i])); 
    
    NormDiff_StatPDFErr[i] = Norm_StatPDFErr[i] / BinWidth[i]; 
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatPDFErr[i]**2 + NormDiff_ScaleErr[i]**2) ; 
    NormDiff_TotalUncer[i] = NormDiff_StatPDFErr[i]; // Scale Error will be added at W-/W+ Ratio code.
    
    //printf("NormDiff_Xsec : %.8f \t ScaleUp %.8f \t ScaleDown %.8f \t PDFp %.8f \t PDFm %.8f %\n ",NormDiff_Xsec[i],NormDiff_Xsec_up[i],NormDiff_Xsec_down[i],NormDiff_Xsec_PDFp[i], NormDiff_Xsec_PDFm[i]);
    //printf("NormDiff_Xsec : %.8f \t +- %.8f \t +- %.8f \t +- %.8f \t Total Uncer : %.8f\t %.2f %\n ",NormDiff_Xsec[i],NormDiff_StatErr[i],NormDiff_ScaleErr[i],NormDiff_PDFErr[i],NormDiff_TotalUncer[i], NormDiff_TotalUncer[i]/NormDiff_Xsec[i]*100);
    //printf("NormDiff_PDFErr : %.8f\t %.2f %\n ",NormDiff_PDFErr[i], NormDiff_PDFErr[i]/NormDiff_Xsec[i]*100);
  }
  
  // Save in histogram
  TH1D* hxsec = new TH1D("hxsec","hxsec",12,Wpt12Bins);
  TH1D* hxsec_up = new TH1D("hxsec_up","hxsec Scale up",12,Wpt12Bins);
  TH1D* hxsec_down = new TH1D("hxsec_down","hxsec Scale down",12,Wpt12Bins);
  
  TH1D* hxsec_NormDiff = new TH1D("hxsec_NormDiff","hxsec_NormDiff",12,Wpt12Bins);
  TH1D* hxsec_NormDiff_up = new TH1D("hxsec_NormDiff_up","hxsec_NormDiff_up",12,Wpt12Bins);
  TH1D* hxsec_NormDiff_down = new TH1D("hxsec_NormDiff_down","hxsec_NormDiff_down",12,Wpt12Bins);

  TH1D* hStatErr = new TH1D("StatErr","Stat Error", 12,Wpt12Bins);
  TH1D* hPDFErr = new TH1D("PDFErr","PDF Error", 12,Wpt12Bins);

  for(int i(0); i<12; i++)
  {
    hxsec->SetBinContent(i+1,Xsec[i]);
    hxsec->SetBinError(i+1,StatPDFErr[i]);

    hxsec_up->SetBinContent(i+1, Xsec_up[i]);
    hxsec_down->SetBinContent(i+1, Xsec_down[i]);

    hxsec_NormDiff->SetBinContent(i+1,NormDiff_Xsec[i]);
    hxsec_NormDiff->SetBinError(i+1,NormDiff_TotalUncer[i]);
    
    hxsec_NormDiff_up->SetBinContent(i+1,NormDiff_Xsec_up[i]);
    hxsec_NormDiff_down->SetBinContent(i+1,NormDiff_Xsec_down[i]);

    hStatErr->SetBinContent(i+1, StatErr[i]);
    hPDFErr->SetBinContent(i+1, PDFErr[i]);
  }

  TFile* f_out = new TFile("./root/"+BaseName+"_FEWZ.root","recreate");
  
  hxsec->Write();
  hxsec_up->Write();
  hxsec_down->Write();
  hxsec_NormDiff->Write();
  hxsec_NormDiff_up->Write();
  hxsec_NormDiff_down->Write();

  hStatErr->Write();
  hPDFErr->Write();

  return 0;

}

void NormToyErr(double *_Xsec, double *_Err, double *_Norm_Xsec, double *_Norm_Err)
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

  ///Here calculate total Xsec
  double TotalXsec=0;
  for(int i=0; i<12; ++i) TotalXsec+=_Xsec[i];
  for(int i=0; i<12; ++i) 
  {
    _Norm_Xsec[i]=_Xsec[i]/TotalXsec;
    //cout << "Norm_Xsec : " << _Norm_Xsec[i] << endl;
  }

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

