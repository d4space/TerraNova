#include <stdio>
#include <iostream>

static const int NB = 12;
double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

void NormToyErr(double, double, double, double);
void ErrPropaNormXsec(double *nn, double *dn, double ff[NB], double df[NB]);

int MakeRoot_WpTincl_FEWZ()
{
  // cross-section and errors
  double Xsec[12] = {0.,};
  double Xsec_up[12];
  double Xsec_down[12];

  double Total_Xsec = 0;
  double Total_Xsec_up = 0;
  double Total_Xsec_down = 0;
  
  double StatErr[12];
  double PDFErr[12];
  
  double ScaleErrP[12];
  double ScaleErrM[12];
  double ScaleErr[12];

  double Norm_ScaleErrP[12];
  double Norm_ScaleErrM[12];
  double Norm_ScaleErr[12];
 
  double Norm_PDFErrP[12];
  double Norm_PDFErrM[12];
  double Norm_PDFErr[12];
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
  double NormDiff_ScaleErr[12];
  double NormDiff_PDFErr[12];
  double NormDiff_StatPDFErr[12];
  double NormDiff_TotalUncer[12];

  // Read W+ W- root file and histogram
  TFile* f_Wp = new TFile("./root/WpToMuNu_FEWZ.root");
  TFile* f_Wm = new TFile("./root/WmToMuNu_FEWZ.root");

  TH1D* Wp_hxsec = (TH1D*)f_Wp->Get("hxsec") ->Clone("Wp_hxsec");
  TH1D* Wm_hxsec = (TH1D*)f_Wm->Get("hxsec") ->Clone("Wm_hxsec");

  TH1D* Wp_hxsec_up = (TH1D*)f_Wp->Get("hxsec_up")->Clone("Wp_hxsec_up");
  TH1D* Wp_hxsec_down = (TH1D*)f_Wp->Get("hxsec_down")->Clone("Wp_hxsec_down");
  TH1D* Wm_hxsec_up = (TH1D*)f_Wm->Get("hxsec_up")->Clone("Wm_hxsec_up");
  TH1D* Wm_hxsec_down = (TH1D*)f_Wm->Get("hxsec_down")->Clone("Wm_hxsec_down");

  TH1D* Wp_StatErr = (TH1D*)f_Wp->Get("StatErr")->Clone("Wp_StatErr");
  TH1D* Wm_StatErr = (TH1D*)f_Wm->Get("StatErr")->Clone("Wm_StatErr");

    // Norm PDF uncer in % unit
    Norm_PDFErrP[0]	= 2.63 ; 
    Norm_PDFErrP[1]	= 1.03 ;
    Norm_PDFErrP[2]	= 1.22 ;
    Norm_PDFErrP[3]	= 1.44 ;
    Norm_PDFErrP[4]	= 1.85 ;
    Norm_PDFErrP[5]	= 2.20 ;
    Norm_PDFErrP[6]	= 2.57 ;
    Norm_PDFErrP[7]	= 3.25 ;
    Norm_PDFErrP[8]	= 3.99 ;
    Norm_PDFErrP[9]	= 4.61 ;
    Norm_PDFErrP[10]	=  5.56;
    Norm_PDFErrP[11]	=  5.91;
    
    Norm_PDFErrM[0]	= 2.16 ; 
    Norm_PDFErrM[1]	= 0.91 ;
    Norm_PDFErrM[2]	= 1.21 ;
    Norm_PDFErrM[3]	= 1.74 ;
    Norm_PDFErrM[4]	= 2.34 ;
    Norm_PDFErrM[5]	= 2.85 ;
    Norm_PDFErrM[6]	= 3.39 ;
    Norm_PDFErrM[7]	= 4.42 ;
    Norm_PDFErrM[8]	= 5.69 ;
    Norm_PDFErrM[9]	= 6.85 ;
    Norm_PDFErrM[10]	=  8.89;
    Norm_PDFErrM[11]	=  10.23;
  // Make Inclusive
  for(int i(0); i<12; i++)
  {
    Xsec[i] = Wp_hxsec->GetBinContent(i+1) + Wm_hxsec->GetBinContent(i+1); 
    Xsec_up[i] = Wp_hxsec_up->GetBinContent(i+1) + Wm_hxsec_up->GetBinContent(i+1); 
    Xsec_down[i] = Wp_hxsec_down->GetBinContent(i+1) + Wm_hxsec_down->GetBinContent(i+1); 

    StatErr[i] = sqrt(Wp_StatErr->GetBinContent(i+1)*Wp_StatErr->GetBinContent(i+1) + Wm_StatErr->GetBinContent(i+1) * Wm_StatErr->GetBinContent(i+1));
    Norm_PDFErr[i] = TMath::Max(Norm_PDFErrP[i],Norm_PDFErrM[i]);

    printf("Xsec : %.4f +- %.4f \t ScaleUp : %.4f \t ScaleDown : %.4f \n", Xsec[i],StatErr[i],Xsec_up[i],Xsec_down[i]);
  }
  
  // Calculate ScaleUp,Down
  for(int i(0); i<12; i++)
  {
    Total_Xsec += Xsec[i];
    Total_Xsec_up += Xsec_up[i];
    Total_Xsec_down += Xsec_down[i];
  }

  // Print out Total Xsecs
  printf("Total Xsec : %.4f ScaleUp : %.4f ScaleDown : %.4f \n",Total_Xsec, Total_Xsec_up,Total_Xsec_down);


  for(int i(0); i<12; i++)
  {
    Norm_Xsec[i] = Xsec[i] / Total_Xsec;
    Norm_Xsec_up[i] = Xsec_up[i] / Total_Xsec_up;
    Norm_Xsec_down[i] = Xsec_down[i] / Total_Xsec_down;
  }
  
  // Cross-section and Statistical error Normalization by Toy method
  NormToyErr(Xsec,StatErr,Norm_StatErr);
  
  // Print out All Normalized Xsec
  for(int i(0); i<12; i++)
  {
    //cout << "Norm_Xsec : " << Norm_Xsec[i] << "\t Norm_StatErr : " << Norm_StatErr[i] << endl; 
    printf("Norm_Xsec : %.8f \t +- %.8f \t ScaleUp : %.8f \t ScaleDown : %.8f \n",Norm_Xsec[i],Norm_StatErr[i],Norm_Xsec_up[i],Norm_Xsec_down[i]);
  }

  // Calculate Normalized Differential Xsec and Total Uncertainty
  for(int i(0); i<12; i++)
  {
    NormDiff_Xsec[i] = Norm_Xsec[i] / BinWidth[i] ; 
    NormDiff_Xsec_up[i] = Norm_Xsec_up[i] / BinWidth[i];
    NormDiff_Xsec_down[i] = Norm_Xsec_down[i] / BinWidth[i];
   
    NormDiff_ScaleErr[i] = TMath::Max(fabs(NormDiff_Xsec_up[i] - NormDiff_Xsec[i]),fabs(NormDiff_Xsec_down[i]-NormDiff_Xsec[i]));
    NormDiff_StatErr[i] = Norm_StatErr[i] / BinWidth[i] ; 
    NormDiff_PDFErr[i] = Norm_PDFErr[i] * NormDiff_Xsec[i] * 0.01;
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatPDFErr[i]**2 + NormDiff_PDFErr[i]**2 + NormDiff_ScaleErr[i]**2) ; 
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatPDFErr[i]**2 + NormDiff_ScaleErr[i]**2) ; 
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatPDFErr[i]**2 + NormDiff_PDFErr[i]**2) ; 
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatPDFErr[i]**2 + NormDiff_PDFErr[i]**2) ; 
    NormDiff_TotalUncer[i] = NormDiff_StatErr[i] ; 
    
    //printf("NormDiff_Xsec : %.8f \t ScaleUp %.8f \t ScaleDown %.8f \t PDFp %.8f \t PDFm %.8f %\n ",NormDiff_Xsec[i],NormDiff_Xsec_up[i],NormDiff_Xsec_down[i],NormDiff_Xsec_PDFp[i], NormDiff_Xsec_PDFm[i]);
    printf("NormDiff_Xsec : %.8f \t +-  %.8f \t , %.2f % \n ",NormDiff_Xsec[i],NormDiff_TotalUncer[i], NormDiff_TotalUncer[i]/NormDiff_Xsec[i]*100);
    //printf("NormDiff_PDFErr : %.8f\t %.2f %\n ",NormDiff_PDFErr[i], NormDiff_PDFErr[i]/NormDiff_Xsec[i]*100);
  }
  
  // Save in histogram
  TH1D* Wincl_hxsec = new TH1D("Wincl_hxsec","Wincl_hxsec",12,Wpt12Bins);
  TH1D* Wincl_hxsec_up = new TH1D("Wincl_hxsec_up","Wincl_hxsec Scale up",12,Wpt12Bins);
  TH1D* Wincl_hxsec_down = new TH1D("Wincl_hxsec_down","Wincl_hxsec Scale down",12,Wpt12Bins);
  
  TH1D* Wincl_hxsec_NormDiff = new TH1D("hxsec_NormDiff","hxsec_NormDiff",12,Wpt12Bins);
  TH1D* Wincl_hxsec_NormDiff_Up = new TH1D("hxsec_NormDiff_up","hxsec_NormDiff_up",12,Wpt12Bins);
  TH1D* Wincl_hxsec_NormDiff_Down = new TH1D("hxsec_NormDiff_down","hxsec_NormDiff_down",12,Wpt12Bins);
  
  for(int i(0); i<12; i++)
  {
    Wincl_hxsec->SetBinContent(i+1,Xsec[i]);

    Wincl_hxsec_up->SetBinContent(i+1, Xsec_up[i]);
    Wincl_hxsec_down->SetBinContent(i+1, Xsec_down[i]);

    Wincl_hxsec_NormDiff->SetBinContent(i+1,NormDiff_Xsec[i]);
    Wincl_hxsec_NormDiff->SetBinError(i+1,NormDiff_TotalUncer[i]);
    
    Wincl_hxsec_NormDiff_Up->SetBinContent(i+1,NormDiff_Xsec_up[i]);
    Wincl_hxsec_NormDiff_Down->SetBinContent(i+1,NormDiff_Xsec_down[i]);
  }
  
  TFile* fFEWZ = new TFile("./root/WinclToMuNu_FEWZ.root","recreate");
  Wincl_hxsec->Write();
  Wincl_hxsec_up->Write();
  Wincl_hxsec_down->Write();
  Wincl_hxsec_NormDiff->Write();
  Wincl_hxsec_NormDiff_Up->Write();
  Wincl_hxsec_NormDiff_Down->Write();

  return 0;

}

void NormToyErr(double *_Xsec, double *_Err, double *_Norm_Err)
{
  int Ntoy = 100000;

  double mx2[12] = {0.,};
  double m2x[12] = {0.,};
  double temp[12] = {0.,};

  // Check _Xsec and _StatErr argument is correct 
  for(int j=0; j<12; ++j)
  {
    //cout << Form("_Xsec : %.4f \t _StatErr : %.4f",_Xsec[j],_StatErr[j]) << endl;
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
    _Norm_Err[i] = rms[i]; // return rms to Norm_StatErr
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

