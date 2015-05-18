#include <stdio>
#include <iostream>

const double Zpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
const double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

static const int NB = 12;

void NormToyErr(double *Xsec, double *Err, double *Norm_Xsec, double *Norm_Err);
void ErrPropaNormXsec(double *nn, double *dn, double ff[NB], double df[NB]);

//int MakeRoot_ZpT_FEWZ()
int MakeRoot_ZpT_Powheg()
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
  double normPDFErr[12];

  double StatPDFErr[12];

  double Norm_PDFErr[12];

  double ScaleErrP[12];
  double ScaleErrM[12];
  double ScaleErr[12];
  double ScaleErr_percent[12];

  double Norm_ScaleErrP[12];
  double Norm_ScaleErrM[12];
  double Norm_ScaleErr[12];

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


  Xsec[0] = 148759.0   ;  
  Xsec[1] = 91211.00   ;
  Xsec[2] = 56957.00   ;
  Xsec[3] = 73749.00   ;
  Xsec[4] = 29039.00   ;
  Xsec[5] = 17077.00   ;
  Xsec[6] = 17419.00   ;
  Xsec[7] = 11017.00   ;  
  Xsec[8] = 2703.000   ;
  Xsec[9] = 916.0000   ;
  Xsec[10] =498.0000   ; 
  Xsec[11] =231.0000   ;

  StatErr[0]	=385.69288300 ;  
  StatErr[1]	=302.01158918 ;
  StatErr[2]	=238.65665715 ;
  StatErr[3]	=271.56767112 ;
  StatErr[4]	=170.40833313 ;
  StatErr[5]	=130.67899602 ;
  StatErr[6]	=131.98105925 ;
  StatErr[7]	=104.96189785 ;
  StatErr[8]	=51.99038373  ;
  StatErr[9]	=30.26549190  ;
  StatErr[10]	=22.31591360  ; 
  StatErr[11]	=15.19868415  ;

  ////Zpt 12 powheg Norm PDF error in %
  normPDFErr[0]	       =0.280629    ;   
  normPDFErr[1]	       =0.205053    ;  
  normPDFErr[2]	       =0.272395    ;
  normPDFErr[3]	       =0.297853    ;
  normPDFErr[4]	       =0.315051    ;
  normPDFErr[5]	       =0.356816    ;
  normPDFErr[6]	       =0.619421    ;
  normPDFErr[7]	       =1.28964     ;
  normPDFErr[8]	       =2.29144     ;
  normPDFErr[9]	       =2.92386     ;
  normPDFErr[10]       =4.04627     ;
  normPDFErr[11]       =6.01494     ;


  // Error Calculation 
  for(int i(0); i<12; i++)
  {
    //PDFErr[i] = TMath::Max(PDFErrP[i],PDFErrM[i]); // pdf error

    //ScaleErrP[i] = TMath::Abs(Xsec[i] - Xsec_up[i]);
    //ScaleErrM[i] = TMath::Abs(Xsec[i] - Xsec_down[i]);
    //ScaleErr[i] = TMath::Max(ScaleErrP[i],ScaleErrM[i]); // scale error
    
    //StatPDFErr[i] = sqrt(StatErr[i]*StatErr[i] + PDFErr[i]*PDFErr[i]); // will be toy normalization
    
    //printf("Xsec : %f \t Stat : %.4f \t TotalUncer : %.8f, \t %.4f %\n",Xsec[i],StatErr[i],TotalUncer[i], TotalUncer[i]/Xsec[i]*100 );
  }

  // Calculate ScaleUp,Down case Total Xsec
  for(int i(0); i<12; i++)
  {
   // Total_Xsec_up += Xsec_up[i];
   // Total_Xsec_down += Xsec_down[i];
  }
  
  // Print out Total Xsecs
  //printf("Total Xsec ScaleUp : %.4f ScaleDown : %.4f PDFp : %.4f PDFm : %.4f \n",Total_Xsec_up,Total_Xsec_down,Total_Xsec_PDFp,Total_Xsec_PDFm);

  for(int i(0); i<12; i++)
  {
   // Norm_Xsec_up[i] = Xsec_up[i] / Total_Xsec_up; // Calculation of Scale up case  Normalized Xsec 
   // Norm_Xsec_down[i] = Xsec_down[i] / Total_Xsec_down; // Calculation of Scale down case Normalized Xsec
  }

  // Cross-section and Statistical error Normalization by Toy method
  //NormToyErr(Xsec,StatErr,Norm_Xsec,Norm_StatErr);
  NormToyErr(Xsec,StatErr,Norm_Xsec,Norm_StatErr);

  // Print out All Normalized Xsec
  //for(int i(0); i<12; i++)
  //{
    //cout << "Norm_Xsec : " << Norm_Xsec[i] << "\t Norm_StatPDFErr : " << Norm_StatPDFErr[i] << endl; 
  //}

  // Calculate Normalized Xsec and Total uncertainty by using 1st order differential formula ??? 
  //ErrPropaNormXsec(Xsec, TotalUncer, Norm_Xsec, Norm_TotalUncer); 
 
//   cout<<"60 Xsec\t"<< 9.42507+6.22349 <<"\tstat err\t"<<sqrt(0.00561647*0.00561647 + 0.00448416*0.00448416)  <<endl;
//   cout<<"80 Xsec\t"<< 4.18073+2.8357 <<"\tstat err\t"<<sqrt(0.00370456*0.00370456 + 0.00302071*0.00302071 )  <<endl;
//   cout<<"100 Xsec\t"<<1.9777 + 1.43425  <<
//     "\tstat err\t"<<sqrt(0.00249719*0.00249719+0.00215388*0.00215388 )  <<endl;
//   cout<<"120 Xsec\t"<<1.06091+0.800329+0.614703+0.480476 <<
//     "\tstat err\t"<<sqrt(0.00189177*0.00189177+0.00167684*0.00167684+0.00150472*0.00150472+0.00146742*0.00146742 )  <<endl;
//
//   cout<<"170 Xsec\t"<< 0.372581+0.293915+0.2241690+ 0.179449 <<
//     "\tstat err\t"<<sqrt(0.00134351*0.00134351+0.00116163*0.00116163+0.000947886*0.000947886+0.00102777*0.00102777 )  <<endl;
//   
//   cout<<"220 Xsec\t"<< 0.147707+0.119982+0.0980798+0.0810495+0.0667002+0.0558383 <<
//   "\tstat err\t"<<sqrt(0.000398149*0.000398149+0.000338095*0.000338095+0.000301994*0.000301994+0.000277197*0.000277197+0.000258883*0.000258883+0.000246317*0.000246317 )  <<endl;
//   cout<<"425 Xsec\t"<<0.167625+0.0732198+0.0343134+0.0175463+0.00856598+0.00447356+0.00265668  <<
//   "\tstat err\t"<<sqrt(0.000907514*0.000907514+0.000644972*0.000644972+0.00052556*0.00052556+0.000445827*0.000445827+0.000390905*0.000390905+0.000344287*0.000344287+0.000307968*0.000307968 )  <<endl;
   //cout<<"Xsec\t"<<  <<"\tstat err\t"<<sqrt( )  <<endl;
  // Calculate Normalized Differential Xsec and Total Uncertainty
  for(int i(0); i<12; i++)
  {
    NormDiff_Xsec[i] = Norm_Xsec[i] / BinWidth[i] ; 
   // NormDiff_Xsec_up[i] = Norm_Xsec_up[i] / BinWidth[i];
   // NormDiff_Xsec_down[i] = Norm_Xsec_down[i] / BinWidth[i];
    
    NormDiff_StatErr[i] = Norm_StatErr[i] / BinWidth[i]; 
    //NormDiff_ScaleErr[i] = TMath::Max(fabs(NormDiff_Xsec_up[i]-NormDiff_Xsec[i]),fabs(NormDiff_Xsec_down[i]-NormDiff_Xsec[i]));
    //NormDiff_TotalUncer[i] = sqrt(NormDiff_StatErr[i]**2 + (0.01*normPDFErr[i]*NormDiff_Xsec[i])**2) ; 
    
    // we calculated PDF error directly from Z/Winlc ratio. Here take only stat error
    NormDiff_TotalUncer[i] = NormDiff_StatErr[i] ; 
    
    printf("NormDiff_Xsec : %.8f \t +- %.10f \t  %.2f %\n ",NormDiff_Xsec[i],NormDiff_TotalUncer[i], NormDiff_TotalUncer[i]/NormDiff_Xsec[i]*100);
  }
   
  cout << fixed << setprecision(10);
  for(int i(0); i<12; i++)
   {
     cout<<"StatErr\t"<<StatErr[i]<<" \tNor stat error\t"<<StatErr[i]/449576/BinWidth[i]<<"\t NormDiff_StatErr[i] Toy \t"<<NormDiff_StatErr[i]<<endl;
   }
   for(int i(0); i<12; i++)
   {
     cout<<"normPDFErr[i]\t"<<normPDFErr[i]<<" \t NormDiff_Xsec[i]\t "<<NormDiff_Xsec[i] <<"\t Convert to number\t"<<0.01*normPDFErr[i]*NormDiff_Xsec[i]<<endl;
   }
 
 
  // Save in histogram
  TH1D* hxsec = new TH1D("hxsec","hxsec",12,Zpt12Bins);
  //TH1D* hxsec_up = new TH1D("hxsec_up","hxsec Scale up",12,Zpt12Bins);
  //TH1D* hxsec_down = new TH1D("hxsec_down","hxsec Scale down",12,Zpt12Bins);
  
  TH1D* hxsec_NormDiff = new TH1D("hxsec_NormDiff","hxsec_NormDiff",12,Zpt12Bins);
  //TH1D* hxsec_NormDiff_up = new TH1D("hxsec_NormDiff_up","hxsec_NormDiff_up",12,Zpt12Bins);
  //TH1D* hxsec_NormDiff_down = new TH1D("hxsec_NormDiff_down","hxsec_NormDiff_down",12,Zpt12Bins);

  TH1D* hStatErr = new TH1D("StatErr","Stat Error", 12,Zpt12Bins);
  TH1D* hPDFErr = new TH1D("PDFErr","PDF Error", 12,Zpt12Bins);

  for(int i(0); i<12; i++)
  {
    hxsec->SetBinContent(i+1,Xsec[i]);
    hxsec->SetBinError(i+1,StatPDFErr[i]);

    //hxsec_up->SetBinContent(i+1, Xsec_up[i]);
    //hxsec_down->SetBinContent(i+1, Xsec_down[i]);

    hxsec_NormDiff->SetBinContent(i+1,NormDiff_Xsec[i]);
    hxsec_NormDiff->SetBinError(i+1,NormDiff_TotalUncer[i]);
    
    //hxsec_NormDiff_up->SetBinContent(i+1,NormDiff_Xsec_up[i]);
    //hxsec_NormDiff_down->SetBinContent(i+1,NormDiff_Xsec_down[i]);
  }

  TFile* f_out = new TFile("./root/ZToMuMu_Powheg.root","recreate");
  
  hxsec->Write();
  //hxsec_up->Write();
  //hxsec_down->Write();
  hxsec_NormDiff->Write();
  //hxsec_NormDiff_up->Write();
  //hxsec_NormDiff_down->Write();

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
    cout << "ToatXsec: " << TotalXsec << endl;
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

