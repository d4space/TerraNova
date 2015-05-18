#include <stdio>
#include <iostream>

#define Ntoy 50000

double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

static const int NB = 12;

void NormToyErr(double *Xsec, double *Err, double *Norm_Err);

int MakeRoot_WpT_Powheg(TString BaseName)
{
  cout<<"Hello1"<<endl;
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
  
  double TotalXsec=0;
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
    Xsec[0] =1368.600  ;  
    Xsec[1] =808.316   ;
    Xsec[2] =484.039   ;
    Xsec[3] =616.798   ;
    Xsec[4] =232.447   ;
    Xsec[5] =131.952   ;
    Xsec[6] =130.105   ;
    Xsec[7] =77.975    ;  
    Xsec[8] =18.684    ;
    Xsec[9] =5.941     ;
    Xsec[10]= 3.092    ; 
    Xsec[11]= 1.534    ;

    StatErr[0]	=3.273 ;  
    StatErr[1]	=2.515 ;
    StatErr[2]	=1.946 ;
    StatErr[3]	=2.197 ;
    StatErr[4]	=1.349 ;
    StatErr[5]	=1.016 ;
    StatErr[6]	=1.009 ;
    StatErr[7]	=0.781 ;
    StatErr[8]	=0.382 ;
    StatErr[9]	=0.216 ;
    StatErr[10]	=0.156 ; 
    StatErr[11]	=0.110 ;


    //// normalized PDF sysy in %
    Norm_PDFErr[0]	=0.547639 ; 
    Norm_PDFErr[1]	=0.11084  ;
    Norm_PDFErr[2]	=0.423104 ;
    Norm_PDFErr[3]	=0.502498 ;
    Norm_PDFErr[4]	=0.523514 ;
    Norm_PDFErr[5]	=0.456424 ;
    Norm_PDFErr[6]	=0.390481 ;
    Norm_PDFErr[7]	=0.423705 ;
    Norm_PDFErr[8]	=1.07401  ;
    Norm_PDFErr[9]	=1.92522  ;
    Norm_PDFErr[10]	=2.57608  ;
    Norm_PDFErr[11]	=3.85306  ;

  }
  else if(BaseName == "WmToMuNu")
  {

    Xsec[0] =887.724 ;  
    Xsec[1] =541.313 ;
    Xsec[2] =331.008 ;
    Xsec[3] =426.589 ;
    Xsec[4] =165.197 ;
    Xsec[5] =94.502  ;
    Xsec[6] =92.338  ;
    Xsec[7] =54.933  ;  
    Xsec[8] =13.064  ;
    Xsec[9] =4.316   ;
    Xsec[10]=2.032   ; 
    Xsec[11]=0.916   ;

    StatErr[0]	=1.527 ;  
    StatErr[1]	=1.192 ;
    StatErr[2]	=0.932 ;
    StatErr[3]	=1.058 ;
    StatErr[4]	=0.659 ;
    StatErr[5]	=0.498 ;
    StatErr[6]	=0.492 ;
    StatErr[7]	=0.380 ;
    StatErr[8]	=0.185 ;
    StatErr[9]	=0.106 ;
    StatErr[10]	=0.073 ; 
    StatErr[11]	=0.049 ;

    // Norm_PDFErr in %
    Norm_PDFErr[0]	=0.805695 ; 
    Norm_PDFErr[1]	=0.193853 ;
    Norm_PDFErr[2]	=0.529425 ;
    Norm_PDFErr[3]	=0.705878 ;
    Norm_PDFErr[4]	=0.697846 ;
    Norm_PDFErr[5]	=0.707969 ;
    Norm_PDFErr[6]	=0.529504 ;
    Norm_PDFErr[7]	=0.780084 ;
    Norm_PDFErr[8]	=1.43238  ;
    Norm_PDFErr[9]	=2.21681  ;
    Norm_PDFErr[10]	=3.18584  ;
    Norm_PDFErr[11]	=3.50682  ;
  
  }

  cout<<"Hello2"<<endl;
  // Calculate ScaleUp,Down case Total Xsec
  for(int i(0); i<12; i++)
  {
    TotalXsec += Xsec[i];
  }

  
  cout<<"Hello3"<<endl;
  ///Normalized by total
  for(int i=0; i<12; ++i) 
  {
    Norm_Xsec[i]=Xsec[i]/TotalXsec;
    //cout << "Norm_Xsec : " << _Norm_Xsec[i] << endl;
  }
 

  NormToyErr(Xsec,StatErr,Norm_StatErr);


  // Calculate Normalized Differential Xsec and Total Uncertainty
  for(int i(0); i<12; i++)
  {
    NormDiff_Xsec[i] = Norm_Xsec[i] / BinWidth[i] ; 
    
    
    //NormDiff_StatErr[i] = StatErr[i] /TotalXsec/BinWidth[i]; 
    NormDiff_StatErr[i] = Norm_StatErr[i] /BinWidth[i]; 
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatErr[i]**2 + (0.01*Norm_PDFErr[i]*NormDiff_Xsec[i] )**2) ; 
    
    /// we calculated Norm PDF for Powheg directly for W-/W+. So here we save only stat error.
    NormDiff_TotalUncer[i] = NormDiff_StatErr[i] ; 
 
    cout<<"Xsec\t"<<Xsec[i]<<" \t TotXsec\t "<<TotalXsec <<"\tNorm_Xsec\t"<<Norm_Xsec[i]<<"\tNormDiff_Xsec \t"<<NormDiff_Xsec[i]<<endl;


    //printf("NormDiff_Xsec : %.8f \t ScaleUp %.8f \t ScaleDown %.8f \t PDFp %.8f \t PDFm %.8f %\n ",NormDiff_Xsec[i],NormDiff_Xsec_up[i],NormDiff_Xsec_down[i],NormDiff_Xsec_PDFp[i], NormDiff_Xsec_PDFm[i]);
    //printf("NormDiff_Xsec : %.8f \t +- %.8f \t +- %.8f \t +- %.8f \t Total Uncer : %.8f\t %.2f %\n ",NormDiff_Xsec[i],NormDiff_StatErr[i],NormDiff_ScaleErr[i],NormDiff_PDFErr[i],NormDiff_TotalUncer[i], NormDiff_TotalUncer[i]/NormDiff_Xsec[i]*100);
    //printf("NormDiff_PDFErr : %.8f\t %.2f %\n ",NormDiff_PDFErr[i], NormDiff_PDFErr[i]/NormDiff_Xsec[i]*100);
  }
 
  cout << fixed << setprecision(10);
  for(int i(0); i<12; i++)
  {
    cout<<"StatErr\t"<<StatErr[i]<<" \t TotXsec\t "<<TotalXsec <<
      "\tNor stat error\t"<<StatErr[i] /TotalXsec/BinWidth[i]<<
      "\t error in%\t"<<(StatErr[i] /TotalXsec/BinWidth[i])/NormDiff_Xsec[i]*100<<
      "\t NormDiff_StatErr[i] Toy \t"<<NormDiff_StatErr[i]<<
      "\t  Toy in % \t"<<NormDiff_StatErr[i]/NormDiff_Xsec[i]*100<<
      
      endl;
  } 
  for(int i(0); i<12; i++)
  {
    cout<<"Norm_PDFErr[i]\t"<<Norm_PDFErr[i]<<" \t NormDiff_Xsec[i]\t "<<NormDiff_Xsec[i] <<"\t Convert to number\t"<<0.01*Norm_PDFErr[i]*NormDiff_Xsec[i]<<endl;
  } 
    // Save in histogram
  TH1D* hxsec = new TH1D("hxsec","hxsec",12,Wpt12Bins);
  
  TH1D* hxsec_NormDiff = new TH1D("hxsec_NormDiff","hxsec_NormDiff",12,Wpt12Bins);

  TH1D* hStatErr = new TH1D("StatErr","Stat Error", 12,Wpt12Bins);
  TH1D* hPDFErr = new TH1D("PDFErr","PDF Error", 12,Wpt12Bins);

  for(int i(0); i<12; i++)
  {
    hxsec->SetBinContent(i+1,Xsec[i]);
    hxsec->SetBinError(i+1,StatErr[i]);


    hxsec_NormDiff->SetBinContent(i+1,NormDiff_Xsec[i]);
    hxsec_NormDiff->SetBinError(i+1,NormDiff_TotalUncer[i]);
    

    hStatErr->SetBinContent(i+1, StatErr[i]);
    hPDFErr->SetBinContent(i+1, PDFErr[i]);
  }

  TFile* f_out = new TFile("./root/"+BaseName+"_PowhegWpt.root","recreate");
  
  hxsec->Write();
  hxsec_NormDiff->Write();

  hStatErr->Write();
  hPDFErr->Write();

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



