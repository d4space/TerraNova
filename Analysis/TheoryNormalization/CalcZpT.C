#include <iostream.h>

const double LUMI = 546.265815; // # of events / cross-section , by using ZpT twiki number.
const int nbin = 18;
double ZptBins[19] = {0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,30.0,40.0,50.0,70.0,90.0,110.0,150.0,190.0,250.0,600.0};
double ZptBinWidth[18] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,10,10,10,20,20,20,40,40,60,350};

void PDFUncer(double, double);
void NormToyErr(double *_Xsec, double *_Err, double *_Norm_Xsec, double *_Norm_Err);

void CalcZpT()
{
  TFile* f_madgraph = new TFile("./preFSR_Madgraph.root");
  //TFile* f_madgraph = new TFile("./postFSR_Madgraph.root");

  TH1D* hPreFSR_Madgraph = (TH1D*)f_madgraph->Get("hPtMCAcc")->Clone("hPreFSR_Madgraph");

  double Xsec[nbin] = {0.,};
  double StatErr[nbin] = {0.,};
  double PDFErr[nbin] = {0.,};
  double TotalXsec = 0;

  double Norm_Xsec[nbin] = {0.,};
  double Norm_StatErr[nbin] = {0.,};
  double Norm_PDFErr[nbin] = {0.,};
  double NormDiff_Xsec[nbin] = {0.,};
  double NormDiff_StatErr[nbin] = {0.,};
  double NormDiff_PDFErr[nbin] = {0.,};
  double NormDiff_TotalUncer[nbin] = {0.,};

  //Print out # of events in each bin
  for(int i(0); i<nbin; i++)
  {
    cout << i+1 <<" bin Events : " << hPreFSR_Madgraph->GetBinContent(i+1) << endl;
  }

  // Xsec = # of Events / LUMI
  for(int i(0); i<nbin; i++)
  {
    Xsec[i] = hPreFSR_Madgraph->GetBinContent(i+1) / LUMI ;
    StatErr[i] = sqrt(hPreFSR_Madgraph->GetBinContent(i+1)) / LUMI ;
    cout << i+1 <<" bin Xsec : " << Xsec[i] << "\t stat : " << StatErr[i] << " , \t " << StatErr[i]/Xsec[i]*100 << " %" <<endl; // print out Xsec and stat error
  }

  // PDF Error
  PDFUncer(Xsec,PDFErr);
  for(int i(0); i<nbin; i++)
  {
    cout << "PDFErr in main : " << PDFErr[i] << endl; // check PDF Error in main function
  }

/*  
  // Total Xsec 
  for(int i(0); i<nbin; i++)
  {
    TotalXsec += Xsec[i] ;
    //cout << i+1 <<" bin Xsec : " << Xsec[i]  << " Total Xsec : " << TotalXsec << endl;
  }
  // Norm Xsec
  for(int i(0); i<nbin; i++)
  {
    Norm_Xsec[i] = Xsec[i] / TotalXsec;
    //Norm_StatErr[i] = StatErr[i] / TotalXsec;
    cout << i+1 <<" bin Norm_Xsec : " << Norm_Xsec[i] << endl;
  }
*/

  // Norm Stat and PDF Error by using Toy
  NormToyErr(Xsec,StatErr,Norm_Xsec,Norm_StatErr);
  NormToyErr(Xsec,PDFErr,Norm_Xsec,Norm_PDFErr);

  for(int i(0); i<nbin; i++)
  {
    cout << "NormStat : " << Norm_StatErr[i] << "\t NormPDF : " << Norm_PDFErr[i] << endl;
  }
  
  // NormDiff Xsec
  for(int i(0); i<nbin; i++)
  {
    NormDiff_Xsec[i] = Norm_Xsec[i] / ZptBinWidth[i];
    NormDiff_StatErr[i] = Norm_StatErr[i] / ZptBinWidth[i];
    NormDiff_PDFErr[i] = Norm_PDFErr[i] / ZptBinWidth[i];
    NormDiff_TotalUncer[i] = sqrt(NormDiff_StatErr[i]**2 + NormDiff_PDFErr[i]**2);
    cout << i+1 <<" bin NormDiff_Xsec : " << NormDiff_Xsec[i] << "\t +- " << NormDiff_TotalUncer[i] <<  endl;
    //cout << i+1 <<" bin NormDiff_Xsec : " << NormDiff_Xsec[i] << endl;
  }


}

void NormToyErr(double *_Xsec, double *_Err, double *_Norm_Xsec, double *_Norm_Err)
{
  int Ntoy = 100000;

  double mx2[nbin] = {0.,};
  double m2x[nbin] = {0.,};
  double temp[nbin] = {0.,};

  // Check _Xsec and _StatErr argument is correct 
//  for(int j=0; j<nbin; ++j)
//  {
    //cout << Form("_Xsec : %.4f \t _StatErr : %.4f",_Xsec[j],_StatErr[j]) << endl;
//  }

  ///Here calculate total Xsec
  double TotalXsec=0;
  for(int i=0; i<nbin; ++i) TotalXsec+=_Xsec[i];
  for(int i=0; i<nbin; ++i) 
  {
    _Norm_Xsec[i]=_Xsec[i]/TotalXsec;
    //cout << "Norm_Xsec : " << _Norm_Xsec[i] << endl;
  }

  // Here Calculate Normalized Stat Error by using toy
  for(int e=0; e<Ntoy; ++e)
  {
    for(int i=0; i<nbin; ++i)
    {
      	temp[i]=gRandom->Gaus(_Xsec[i], _Err[i]);
    }
    double sum=0;
    for(int i=0; i<nbin; ++i) sum+=temp[i];
    for(int i=0; i<nbin; ++i)
    {
      temp[i]/=sum;
      m2x[i]+=temp[i];
      mx2[i]+=temp[i]*temp[i];
    }
  }

  double rms[nbin];
  for(int i=0; i<nbin; ++i)
  {
    m2x[i]/=Ntoy;
    mx2[i]/=Ntoy;
    rms[i] = sqrt(mx2[i]-m2x[i]*m2x[i]);
  }
 
  // return value
  for(int i=0; i<nbin; ++i)
  {
    _Norm_Err[i] = rms[i]; // return rms to Norm_Err
  }
  for(int i=0; i<nbin; ++i)
  {
    //cout << Form("rms = %.8f",rms[i]) << endl;
  }
}

void PDFUncer(double *Xsec, double *Err)
{
  // % errors
  double PDFErrP[nbin] = {0.,};
  double PDFErrM[nbin] = {0.,};
  double PDFErr[nbin] = {0.,};

  PDFErrP[0] =   5.01461; 
  PDFErrP[1] =   4.89219; 
  PDFErrP[2] =   4.66   ; 
  PDFErrP[3] =   4.55057; 
  PDFErrP[4] =   4.3323 ;
  PDFErrP[5] =   4.07014; 
  PDFErrP[6] =   3.91648; 
  PDFErrP[7] =   3.72249; 
  PDFErrP[8] =   3.54474; 
  PDFErrP[9] =   3.3922 ;
  PDFErrP[10] =  3.27882;
  PDFErrP[11] =  3.13866;
  PDFErrP[12] =  2.87089;
  PDFErrP[13] =  2.69051;
  PDFErrP[14] =  2.43785;
  PDFErrP[15] =  2.34029;
  PDFErrP[16] =  2.20104;
  PDFErrP[17] =  2.65944;

  PDFErrM[0] =   5.01402; 
  PDFErrM[1] =   4.853  ; 
  PDFErrM[2] =   4.56017; 
  PDFErrM[3] =   4.46248; 
  PDFErrM[4] =   4.25249;
  PDFErrM[5] =   4.01871; 
  PDFErrM[6] =   3.8786 ; 
  PDFErrM[7] =   3.74101; 
  PDFErrM[8] =   3.6636 ; 
  PDFErrM[9] =   3.60969;
  PDFErrM[10] =  3.50128;
  PDFErrM[11] =  3.35162;
  PDFErrM[12] =  3.08409;
  PDFErrM[13] =  2.99113;
  PDFErrM[14] =  2.75548;
  PDFErrM[15] =  2.59204;
  PDFErrM[16] =  2.36617;
  PDFErrM[17] =  2.68793;

  for(int i(0); i<nbin; i++)
  {
    PDFErr[i] = TMath::Max(PDFErrP[i],PDFErrM[i]);  // choosing bigger number
    Err[i] = Xsec[i] * PDFErr[i] * 0.01; //return PDF error as number
    cout << "PDFErr in function : " << Err[i] << endl; // print out PDF Error in PDFErr function
  }
}

