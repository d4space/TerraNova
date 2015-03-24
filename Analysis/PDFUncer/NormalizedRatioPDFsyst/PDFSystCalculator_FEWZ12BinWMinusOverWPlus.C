#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <TMath.h>

const int Nbin = 12; // Number of bin
const int Nev = 50; // Number of Eigen vector

void PDFSystCalculator_FEWZ12BinWMinusOverWPlus()
{
  //ofstream Fout;
  //Fout.open("PDF_Systematics.txt");

  ifstream FEWZ_Wp_Center("./FEWZ_Wp_Central/WpT_Center.dat");
  ifstream FEWZ_Wm_Center("./FEWZ_Wm_Central/WpT_Center.dat");

  double Wp_EigenVector[Nbin][Nev] = {0.,};
  double Wp_TotalEigenVector[Nev] = {0.,}; // total Xsec for EigenVector
  double Wp_Norm_EigenVector[Nbin][Nev] = {0.,}; // total Xsec for EigenVector
 
  double Wm_EigenVector[Nbin][Nev] = {0.,};
  double Wm_TotalEigenVector[Nev] = {0.,}; // total Xsec for EigenVector
  double Wm_Norm_EigenVector[Nbin][Nev] = {0.,}; // total Xsec for EigenVector
  
  double Wp_Center[Nbin] = {0.,};
  double Wp_TotalCenter = 0; // total Xsec for center value
  double Wp_Norm_Center[Nbin] = {0.,};
  
  double Wm_Center[Nbin] = {0.,};
  double Wm_TotalCenter = 0; // total Xsec for center value
  double Wm_Norm_Center[Nbin] = {0.,};
 
  double WmWpRatio_Norm_Center[Nbin] = {0.,};
  double WmWpRatio_Norm_EigenVector[Nbin][Nev] = {0.,};

  char tmpName[50];
  TString tmp;

  // read eigenvector dat file and save in array : EigenVector[Bin][eigenvector]
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    FEWZ_Wp_Center >> tmp >> Wp_Center[ibin] >> tmp >> tmp >> tmp; // Center Xsec read
    FEWZ_Wm_Center >> tmp >> Wm_Center[ibin] >> tmp >> tmp >> tmp; // Center Xsec read
    cout << "Wp_Center " << ibin << " : " << Wp_Center[ibin] << "\t Wm_Center : "<<Wm_Center[ibin] << endl;
    
    sprintf(tmpName,"./FEWZ_Wp_EigenVector/WpT_Bin%d.dat",ibin+1);
    cout << "file : " << tmpName << endl; // print out file Name
    ifstream FEWZ_Wp_EigenVector(tmpName);
    sprintf(tmpName,"./FEWZ_Wm_EigenVector/WpT_Bin%d.dat",ibin+1);
    ifstream FEWZ_Wm_EigenVector(tmpName);
    for(int j(0); j<Nev; j++)
    {
      FEWZ_Wp_EigenVector >> tmp >> Wp_EigenVector[ibin][j] >> tmp ; // Eigenvector Xsec read
      FEWZ_Wm_EigenVector >> tmp >> Wm_EigenVector[ibin][j] >> tmp ; // Eigenvector Xsec read
      //cout << "Wp_EigenVector["<<ibin<<"]["<<j<<"]" << Wp_EigenVector[ibin][j] << "\t Wp_EigenVector["<<ibin<<"]["<<j<<"]" << Wm_EigenVector[ibin][j] << endl; // print out EigenVector
      
      Wp_TotalEigenVector[j] += Wp_EigenVector[ibin][j]; // Calculate Total Eigenvector
      Wm_TotalEigenVector[j] += Wm_EigenVector[ibin][j]; // Calculate Total Eigenvector
    }
    Wp_TotalCenter += Wp_Center[ibin]; // Calculate Total Center (Total Xsec for center value)
    Wm_TotalCenter += Wm_Center[ibin]; // Calculate Total Center (Total Xsec for center value)
  }
  cout << "Wp_TotalCenter : " << Wp_TotalCenter << "\t Wp_TotalEgigen[1] : " << Wp_TotalEigenVector[1] << endl;

  // Normalized Eigenvector and center
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    for(int j(0); j<Nev; j++)
    {
      Wp_Norm_EigenVector[ibin][j] = Wp_EigenVector[ibin][j] / Wp_TotalEigenVector[j] ;
      Wm_Norm_EigenVector[ibin][j] = Wm_EigenVector[ibin][j] / Wm_TotalEigenVector[j] ;
      cout << Form("Wp_Norm_EigenVector[%d][%d] : %.8f \t Wm_Norm_EigenVector[%d][%d] : %.8f ",ibin,j,Wp_Norm_EigenVector[ibin][j], ibin,j,Wm_Norm_EigenVector[ibin][j])<< endl;
    }
    Wp_Norm_Center[ibin] = Wp_Center[ibin] / Wp_TotalCenter[ibin];
    Wm_Norm_Center[ibin] = Wm_Center[ibin] / Wm_TotalCenter[ibin];
  }

  // Normalized Wm/Wp ratio for Eigenvector and center
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    for(int j(0); j<Nev; j++)
    {
      WmWpRatio_Norm_EigenVector[ibin][j] = Wm_Norm_EigenVector[ibin][j] / Wp_Norm_EigenVector[ibin][j];
      cout << Form("Norm_WmWpRatio[%d][%d] : %.8f",ibin,j,WmWpRatio_Norm_EigenVector[ibin][j]) << endl;
    }
    WmWpRatio_Norm_Center[ibin] = Wm_Norm_Center[ibin] / Wp_Norm_Center[ibin];
  }


  // Xsec level PDF Error 
  double pluserr[Nbin]={0.,};
  double minuserr[Nbin]={0.,};
  // Calculate PDF uncertainty
  for(int ibin(0); ibin<Nbin; ibin++)
  {  
    for(int i=0; i<Nev/2;i++)
    {
      if (  (WmWpRatio_Norm_EigenVector[ibin][2*i] - WmWpRatio_Norm_Center[ibin]) <0 && (WmWpRatio_Norm_EigenVector[ibin][2*i+1] - WmWpRatio_Norm_Center[ibin]) <0 )
	pluserr[ibin] = pluserr[ibin] ;
      else
	pluserr[ibin] = pluserr[ibin] + TMath::Power( (TMath::Max(WmWpRatio_Norm_EigenVector[ibin][2*i] - WmWpRatio_Norm_Center[ibin], WmWpRatio_Norm_EigenVector[ibin][2*i+1] - WmWpRatio_Norm_Center[ibin])) ,2);
   
      if (  (WmWpRatio_Norm_Center[ibin] - WmWpRatio_Norm_EigenVector[ibin][2*i]) <0 && (WmWpRatio_Norm_Center[ibin] - WmWpRatio_Norm_EigenVector[ibin][2*i+1]) <0 )
	minuserr[ibin] = minuserr[ibin];
      else
	minuserr[ibin] = minuserr[ibin] + TMath::Power( (TMath::Max(WmWpRatio_Norm_Center[ibin] - WmWpRatio_Norm_EigenVector[ibin][2*i], WmWpRatio_Norm_Center[ibin] - WmWpRatio_Norm_EigenVector[ibin][2*i+1])) ,2);
    }
  }

  // make PDF uncer
  double PDFsystPlus[Nbin];
  double PDFsystMinus[Nbin];
  double PDFsyst[Nbin];
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    PDFsystPlus[ibin] = TMath::Sqrt(pluserr[ibin]);
    PDFsystMinus[ibin] = TMath::Sqrt(minuserr[ibin]);
    cout << Form("PDFsystPlus[%d] : %.8f, %.2f % \t PDFsystMinus[%d] : %.8f, %.2f % \t max : %.2f %",ibin,PDFsystPlus[ibin], PDFsystPlus[ibin] / WmWpRatio_Norm_Center[ibin] *100, ibin, PDFsystMinus[ibin], PDFsystMinus[ibin]/WmWpRatio_Norm_Center[ibin]*100, TMath::Max(PDFsystPlus[ibin] / WmWpRatio_Norm_Center[ibin] *100,PDFsystMinus[ibin]/WmWpRatio_Norm_Center[ibin]*100)) << endl;
  }
  
}
