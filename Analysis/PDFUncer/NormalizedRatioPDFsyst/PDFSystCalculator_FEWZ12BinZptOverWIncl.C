#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <TMath.h>

const int Nbin = 12; // Number of bin
const int Nev = 50; // Number of Eigen vector

void PDFSystCalculator_Powheg12BinZptOverWIncl()
{
  //ofstream Fout;
  //Fout.open("PDF_Systematics.txt");

  ifstream FEWZ_Z_Center("./Z_Central/ZpT_Center.dat");
  ifstream FEWZ_Wp_Center("./Wincl_Central/Wp_Center.dat");
  ifstream FEWZ_Wm_Center("./Wincl_Central/Wm_Center.dat");

  double Wp_EigenVector[Nbin][Nev] = {0.,};
  double Wm_EigenVector[Nbin][Nev] = {0.,};
  double Wincl_EigenVector[Nbin][Nev] = {0.,};
  double Wincl_TotalEigenVector[Nev] = {0.,}; // total Xsec for EigenVector
  double Wincl_Norm_EigenVector[Nbin][Nev] = {0.,}; 
 
  double Z_EigenVector[Nbin][Nev] = {0.,};
  double Z_TotalEigenVector[Nev] = {0.,}; // total Xsec for EigenVector
  double Z_Norm_EigenVector[Nbin][Nev] = {0.,}; 
  
  double Wp_Center[Nbin] = {0.,};
  double Wm_Center[Nbin] = {0.,};
  double Wincl_Center[Nbin] = {0.,};
  double Wincl_TotalCenter = 0;
  double Wincl_Norm_Center[Nbin] = {0.,};
  
  double Z_Center[Nbin] = {0.,};
  double Z_TotalCenter = 0;
  double Z_Norm_Center[Nbin] = {0.,};

  double ZWRatio_Norm_Center[Nbin] = {0.,};
  double ZWRatio_Norm_EigenVector[Nbin][Nev] = {0.,};

  char Wp_tmpName[100];
  char Wm_tmpName[100];
  char Z_tmpName[100];
  
  double tmp;

  // read eigenvector dat file and save in array : EigenVector[Bin][eigenvector]
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    FEWZ_Wp_Center >> Wp_Center[ibin] ; // Wp Xsec read
    FEWZ_Wm_Center >> Wm_Center[ibin] ; // Wm Xsec read
    FEWZ_Z_Center >> tmp >> Z_Center[ibin] >> tmp>> tmp >> tmp ; // Z Xsec read
    Wincl_Center[ibin] = Wp_Center[ibin] + Wm_Center[ibin];
    cout << "Wincl Center " << ibin << " : " << Wincl_Center[ibin] << "\t Z Center : " << Z_Center[ibin] <<  endl;
    
    sprintf(Wp_tmpName,"./Wincl_EigenVector/Wp/WpT_Bin%d.dat",ibin+1);
    cout << "Wp_tmpName : " << Wp_tmpName << endl;
    ifstream FEWZ_Wp_EigenVector(Wp_tmpName);
    sprintf(Wm_tmpName,"./Wincl_EigenVector/Wm/WpT_Bin%d.dat",ibin+1);
    cout << "Wm_tmpName : " << Wm_tmpName << endl;
    ifstream FEWZ_Wm_EigenVector(Wm_tmpName);
    sprintf(Z_tmpName,"./Z_EigenVector/ZpT_Bin%d.dat",ibin+1);
    cout << "Z_tmpName : " << Z_tmpName << endl;
    ifstream FEWZ_Z_EigenVector(Z_tmpName);
    for(int j(0); j<Nev; j++)
    {
      FEWZ_Wp_EigenVector >> tmp >> Wp_EigenVector[ibin][j] >> tmp ; // Eigenvector Xsec read
      FEWZ_Wm_EigenVector >> tmp >> Wm_EigenVector[ibin][j] >> tmp ; // Eigenvector Xsec read
      FEWZ_Z_EigenVector >> tmp >> Z_EigenVector[ibin][j] >> tmp ; // Eigenvector Xsec read
      
      Wincl_EigenVector[ibin][j] = Wp_EigenVector[ibin][j] + Wm_EigenVector[ibin][j];
      cout << "Wincl EigenVector["<<ibin<<"]["<<j<<"]" << Wincl_EigenVector[ibin][j] << endl; // print out EigenVector
      
      Wincl_TotalEigenVector[j] += Wincl_EigenVector[ibin][j]; // Calculate Total Eigenvector
      Z_TotalEigenVector[j] += Z_EigenVector[ibin][j]; // Calculate Total Eigenvector
      //cout << "TotalEigenVector : " << TotalEigenVector[j] <<endl;
    }
    FEWZ_Wp_EigenVector.close();
    FEWZ_Wm_EigenVector.close();
    Wincl_TotalCenter += Wincl_Center[ibin]; // Calculate Total Center (Total Xsec for center value)
    Z_TotalCenter += Z_Center[ibin]; // Calculate Total Center (Total Xsec for center value)
  }
  cout << "Wincl_TotalCenter : " << Wincl_TotalCenter << "\t Z_TotalCenter : " << Z_TotalCenter << endl;
  
  // Normalized Eigenvector and center
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    for(int j(0); j<Nev; j++)
    {
      Wincl_Norm_EigenVector[ibin][j] = Wincl_EigenVector[ibin][j] / Wincl_TotalEigenVector[j] ;
      Z_Norm_EigenVector[ibin][j] = Z_EigenVector[ibin][j] / Z_TotalEigenVector[j] ;
      cout << Form("Wincl_Norm_EigenVector[%d][%d] : %.8f \t Z_Norm_EigenVector[%d][%d] : %.8f ",ibin,j,Wincl_Norm_EigenVector[ibin][j],ibin,j,Z_Norm_EigenVector[ibin][j])<< endl;
    }
    Wincl_Norm_Center[ibin] = Wincl_Center[ibin] / Wincl_TotalCenter[ibin];
    Z_Norm_Center[ibin] = Z_Center[ibin] / Z_TotalCenter[ibin];
  }

  // Normalized Z/Wincl ratio for Eigenvector and center
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    for(int j(0); j<Nev; j++)
    {
      ZWRatio_Norm_EigenVector[ibin][j] = Z_Norm_EigenVector[ibin][j] / Wincl_Norm_EigenVector[ibin][j];
      cout << Form("Norm_ZWRatio[%d][%d] : %.8f",ibin,j,ZWRatio_Norm_EigenVector[ibin][j]) << endl;
    }
    ZWRatio_Norm_Center[ibin] = Z_Norm_Center[ibin] / Wincl_Norm_Center[ibin];
  }


/*
  // Xsec level PDF Error 
  double pluserr[Nbin]={0.,};
  double minuserr[Nbin]={0.,};
  // Calculate PDF uncertainty
  for(int ibin(0); ibin<Nbin; ibin++)
  {  
    for(int i=0; i<Nev/2;i++)
    {
      if (  (Wincl_EigenVector[ibin][2*i] - Wincl_Center[ibin]) <0 && (Wincl_EigenVector[ibin][2*i+1] - Wincl_Center[ibin]) <0 )
	pluserr[ibin] = pluserr[ibin] ;
      else
	pluserr[ibin] = pluserr[ibin] + TMath::Power( (TMath::Max(Wincl_EigenVector[ibin][2*i] - Wincl_Center[ibin], Wincl_EigenVector[ibin][2*i+1]-Wincl_Center[ibin])) ,2);
   
      if (  (Wincl_Center[ibin] - Wincl_EigenVector[ibin][2*i]) <0 && (Wincl_Center[ibin] - Wincl_EigenVector[ibin][2*i+1]) <0 )
	minuserr[ibin] = minuserr[ibin];
      else
	minuserr[ibin] = minuserr[ibin] + TMath::Power( (TMath::Max(Wincl_Center[ibin] - Wincl_EigenVector[ibin][2*i], Wincl_Center[ibin] - Wincl_EigenVector[ibin][2*i+1])) ,2);
    }
  }
*/

  // Normalized Xsec level calculation
  double Norm_pluserr[Nbin]={0.,};
  double Norm_minuserr[Nbin]={0.,};
  // Calculate PDF uncertainty
  for(int ibin(0); ibin<Nbin; ibin++)
  {  
    for(int i=0; i<Nev/2;i++)
    {
      if (  (ZWRatio_Norm_EigenVector[ibin][2*i] - ZWRatio_Norm_Center[ibin]) <0 && (ZWRatio_Norm_EigenVector[ibin][2*i+1] - ZWRatio_Norm_Center[ibin]) <0 )
	Norm_pluserr[ibin] = Norm_pluserr[ibin] ;
      else
	Norm_pluserr[ibin] = Norm_pluserr[ibin] + TMath::Power( (TMath::Max(ZWRatio_Norm_EigenVector[ibin][2*i] - ZWRatio_Norm_Center[ibin], ZWRatio_Norm_EigenVector[ibin][2*i+1]-ZWRatio_Norm_Center[ibin])) ,2);
   
      if (  (ZWRatio_Norm_Center[ibin] - ZWRatio_Norm_EigenVector[ibin][2*i]) <0 && (ZWRatio_Norm_Center[ibin] - ZWRatio_Norm_EigenVector[ibin][2*i+1]) <0 )
	Norm_minuserr[ibin] = Norm_minuserr[ibin];
      else
	Norm_minuserr[ibin] = Norm_minuserr[ibin] + TMath::Power( (TMath::Max(ZWRatio_Norm_Center[ibin] - ZWRatio_Norm_EigenVector[ibin][2*i], ZWRatio_Norm_Center[ibin] - ZWRatio_Norm_EigenVector[ibin][2*i+1])) ,2);
    }
  }

/*
  // make PDF uncer
  double PDFsystPlus[Nbin];
  double PDFsystMinus;
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    PDFsystPlus[ibin] = TMath::Sqrt(pluserr[ibin]);
    PDFsystMinus[ibin] = TMath::Sqrt(minuserr[ibin]);
    cout << Form("PDFsystPlus[%d] : %.8f, %.2f % \t PDFsystMinus[%d] : %.8f, %.2f % ",ibin,PDFsystPlus[ibin], PDFsystPlus[ibin] / Center[ibin] *100, ibin, PDFsystMinus[ibin], PDFsystMinus[ibin]/Center[ibin]*100) << endl;;
  }
  */
  // make Normalized PDF uncer
  double Norm_PDFsystPlus[Nbin];
  double Norm_PDFsystMinus;
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    Norm_PDFsystPlus[ibin] = TMath::Sqrt(Norm_pluserr[ibin]);
    Norm_PDFsystMinus[ibin] = TMath::Sqrt(Norm_minuserr[ibin]);
    cout << Form("Norm_PDFsystPlus[%d] : %.8f, %.2f % \t Norm_PDFsystMinus[%d] : %.8f, %.2f % \t max : %.2f %",ibin,Norm_PDFsystPlus[ibin], Norm_PDFsystPlus[ibin] / ZWRatio_Norm_Center[ibin] *100, ibin, Norm_PDFsystMinus[ibin], Norm_PDFsystMinus[ibin]/ZWRatio_Norm_Center[ibin]*100, TMath::Max(Norm_PDFsystPlus[ibin] / ZWRatio_Norm_Center[ibin] *100 , Norm_PDFsystMinus[ibin]/ZWRatio_Norm_Center[ibin]*100)) << endl;;
  }
  //Fout.close();
}
