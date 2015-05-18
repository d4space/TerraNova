#include <stdio>
#include <iostream>

#define Ntoy 50000

double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

static const int NB = 12;

void NormToyErr(double *Xsec, double *Err, double *Norm_Err);

int MakeRoot_InclWpT_Powheg(TString BaseName)
{
  cout<<"Hello1"<<endl;
  double Xsec[12];
  double XsecP[12];
  double XsecM[12];
  
  double Total_Xsec = 0;
  double Total_XsecP = 0;
  double Total_XsecM = 0;
  
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

  if(BaseName == "InclWToMuNu")
  {
    XsecP[0] =1368.600  ;  
    XsecP[1] =808.316   ;
    XsecP[2] =484.039   ;
    XsecP[3] =616.798   ;
    XsecP[4] =232.447   ;
    XsecP[5] =131.952   ;
    XsecP[6] =130.105   ;
    XsecP[7] =77.975    ;  
    XsecP[8] =18.684    ;
    XsecP[9] =5.941     ;
    XsecP[10]= 3.092    ; 
    XsecP[11]= 1.534    ;
    
    XsecM[0] =887.724    ;  
    XsecM[1] =541.313    ;
    XsecM[2] =331.008    ;
    XsecM[3] =426.589    ;
    XsecM[4] =165.197    ;
    XsecM[5] =94.502     ;
    XsecM[6] =92.338     ;
    XsecM[7] =54.933     ;  
    XsecM[8] =13.064     ;
    XsecM[9] =4.316      ;
    XsecM[10]=2.032      ; 
    XsecM[11]=0.916      ;

    StatErr[0]	=3.6115 ;  
    StatErr[1]	=2.7835 ;
    StatErr[2]	=2.1582 ;
    StatErr[3]	=1.9331 ;
    StatErr[4]	=1.5010 ;
    StatErr[5]	=1.1318 ;
    StatErr[6]	=1.1228 ;
    StatErr[7]	=0.8686 ;
    StatErr[8]	=0.4249 ;
    StatErr[9]	=0.2405 ;
    StatErr[10]	=0.1719 ; 
    StatErr[11]	=0.1201 ;


    //// normalized PDF sysy in %

    Norm_PDFErr[0]	= 0.71170  ; 
    Norm_PDFErr[1]	= 0.16135  ;
    Norm_PDFErr[2]	= 0.49645  ;
    Norm_PDFErr[3]	= 0.63787  ;
    Norm_PDFErr[4]	= 0.63341  ;
    Norm_PDFErr[5]	= 0.58773  ;
    Norm_PDFErr[6]	= 0.42823  ;
    Norm_PDFErr[7]	= 0.56209  ;
    Norm_PDFErr[8]	= 1.26718  ;
    Norm_PDFErr[9]	= 2.07182  ;
    Norm_PDFErr[10]	= 2.97034  ;
    Norm_PDFErr[11]	= 3.63010  ;

  }

  cout << fixed << setprecision(10);
  cout<<"Hello2"<<endl;
  // Calculate ScaleUp,Down case Total Xsec
  for(int i(0); i<12; i++)
  {
    Xsec[i] = ( XsecP[i] + XsecM[i] );
    TotalXsec += ( XsecP[i] + XsecM[i] );
  }

  
  cout<<"Hello3"<<endl;
  ///Normalized by total
  for(int i=0; i<12; ++i) 
  {
    Norm_Xsec[i]= ( XsecP[i] + XsecM[i] )/TotalXsec;
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
    
    /// we calculated Z/W Pdf error directly for Z/W ration. so take only stat error here.
    NormDiff_TotalUncer[i] = NormDiff_StatErr[i] ; 
 
    cout<<"Xsec\t"<<Xsec[i]<<" \t TotXsec\t "<<TotalXsec <<"\tNorm_Xsec\t"<<Norm_Xsec[i]<<"\tNormDiff_Xsec \t"<<NormDiff_Xsec[i]<<endl;


    //printf("NormDiff_Xsec : %.8f \t ScaleUp %.8f \t ScaleDown %.8f \t PDFp %.8f \t PDFm %.8f %\n ",NormDiff_Xsec[i],NormDiff_Xsec_up[i],NormDiff_Xsec_down[i],NormDiff_Xsec_PDFp[i], NormDiff_Xsec_PDFm[i]);
    //printf("NormDiff_Xsec : %.8f \t +- %.8f \t +- %.8f \t +- %.8f \t Total Uncer : %.8f\t %.2f %\n ",NormDiff_Xsec[i],NormDiff_StatErr[i],NormDiff_ScaleErr[i],NormDiff_PDFErr[i],NormDiff_TotalUncer[i], NormDiff_TotalUncer[i]/NormDiff_Xsec[i]*100);
    //printf("NormDiff_PDFErr : %.8f\t %.2f %\n ",NormDiff_PDFErr[i], NormDiff_PDFErr[i]/NormDiff_Xsec[i]*100);
  }
  
  for(int i(0); i<12; i++)
  {
    //cout<<"StatErr\t"<<StatErr[i]<<" \t TotXsec\t "<<TotalXsec <<"\tNorm_DiffStatErr[i]\t"<<NormDiff_StatErr[i]<<"\t NormDiff_TotalUncer[i] \t"<<NormDiff_TotalUncer[i]<<endl;
    cout<<"StatErr\t"<<StatErr[i]<<" \t TotXsec\t "<<TotalXsec <<"\tNor stat error\t"<<StatErr[i] /TotalXsec/BinWidth[i]<<"\t NormDiff_StatErr[i] Toy \t"<<NormDiff_StatErr[i]<<endl;
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


