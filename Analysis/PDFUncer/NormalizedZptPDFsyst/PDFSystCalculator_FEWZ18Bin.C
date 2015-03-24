#include <iostream.h>
#include <fstream>

void NormToyErr(double *_Center, double *_Err, double *_Norm_Err);

double ZptBinWidth[18] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,10,10,10,20,20,20,40,40,60,350};

const double Correction = 1.074;

const int Nev = 50; // number of CT10 eigenvector
const int Nbin = 38; // dimitri binning
const int Nbin_ZpT = 18; // ZpT binning

void PDFSystCalculator_FEWZ18Bin()
{
  ifstream FILE_EigenVector("./FEWZ_PDFInfo_fromDimitri.txt"); // 38 bin from Dimitri

  double EigenVector[Nbin][Nev]={0.,};
  double EigenVector_18Bin[Nbin_ZpT][Nev]={0.,};
  TString tmp;

  // read Dimitri 38 bin pdf eigenvectors and save it as array : EigenVector[bin][Nev]
  for(int ibin(0); ibin<Nbin; ibin++)
  {
    for(int iev(0); iev<Nev; iev++)
    {
      FILE_EigenVector >> tmp >> EigenVector[ibin][iev] >> tmp;
      //cout << "EigenVector["<<ibin<<"]["<<iev<<"] : "<<EigenVector[ibin][iev] << endl;
    }
  }

  //============== start merging : 38 Bin -> 18 bin ===============
  // bin 1 ~ 11 is same
  for(int ibin(0); ibin<11; ibin++)
  {
    for(int iev(0); iev<Nev; iev++)
    {
      EigenVector_18Bin[ibin][iev] = EigenVector[ibin][iev];
      cout << "EigenVector_18Bin["<<ibin<<"]["<<iev<<"] : "<<EigenVector_18Bin[ibin][iev] << endl;
    }
  }

  // bin 12~13 merge -> ZpTbin 12
  for(int iev(0); iev<Nev; iev++)
  {
    EigenVector_18Bin[11][iev] = EigenVector[11][iev]+ EigenVector[12][iev];
    cout << "EigenVector_18Bin[11]["<<iev<<"] : "<<EigenVector_18Bin[11][iev] << endl;
  }

  // bin 14 ~ 15 merge -> ZpTbin 13
  for(int iev(0); iev<Nev; iev++)
  {
    EigenVector_18Bin[12][iev] = EigenVector[13][iev]+ EigenVector[14][iev];
    cout << "EigenVector_18Bin[12]["<<iev<<"] : "<<EigenVector_18Bin[12][iev] << endl;
  }

  // bin 16 ~ 17 merge -> ZpTbin 14
  for(int iev(0); iev<Nev; iev++)
  {
    EigenVector_18Bin[13][iev] = EigenVector[15][iev]+ EigenVector[16][iev];
    cout << "EigenVector_18Bin[13]["<<iev<<"] : "<<EigenVector_18Bin[13][iev] << endl;
  }
  
  // bin 18 ~ 21 merge -> ZpTbin 15
  for(int iev(0); iev<Nev; iev++)
  {
    EigenVector_18Bin[14][iev] = EigenVector[17][iev]+ EigenVector[18][iev]+EigenVector[19][iev]+EigenVector[20][iev];
    cout << "EigenVector_18Bin[14]["<<iev<<"] : "<<EigenVector_18Bin[14][iev] << endl;
  }

  // bin 22 ~ 25 merge -> ZpTbin 16
  for(int iev(0); iev<Nev; iev++)
  {
    EigenVector_18Bin[15][iev] = EigenVector[21][iev]+ EigenVector[22][iev]+ EigenVector[23][iev] + EigenVector[24][iev] ;
    cout << "EigenVector_18Bin[15]["<<iev<<"] : "<<EigenVector_18Bin[15][iev] << endl;
  }

  // bin 26 ~ 31 merge -> ZpTbin 17
  for(int iev(0); iev<Nev; iev++)
  {
    EigenVector_18Bin[16][iev] = EigenVector[25][iev]+ EigenVector[26][iev]+ EigenVector[27][iev] + EigenVector[28][iev]+ EigenVector[29][iev]+EigenVector[30][iev] ;
    cout << "EigenVector_18Bin[16]["<<iev<<"] : "<<EigenVector_18Bin[16][iev] << endl;
  }

  // bin 32 ~ 38 merge -> ZpTbin 18
  for(int iev(0); iev<Nev; iev++)
  {
    EigenVector_18Bin[17][iev] = EigenVector[31][iev]+ EigenVector[32][iev]+ EigenVector[33][iev] + EigenVector[34][iev] + EigenVector[35][iev] + EigenVector[36][iev] + EigenVector[37][iev];
    cout << "EigenVector_18Bin[17]["<<iev<<"] : "<<EigenVector_18Bin[17][iev] << endl;
  }
// ============== merging finish here ========================


  // Read center value
  ifstream FILE_Central("./FEWZ_CentralValue_18bin.txt");

  double Center[Nbin_ZpT] = {0.};
  double StatErr[Nbin_ZpT] = {0.};

  double TotalCenter = 0;

  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {
    FILE_Central >> Center[ibin] >> StatErr[ibin];
    Center[ibin] = Center[ibin] * Correction; // correction factor by dimitri 1.074
    StatErr[ibin] = StatErr[ibin] * Correction; // correction factor by dimitri 1.074
    TotalCenter += Center[ibin]; // Total Center value
    cout << "Center["<<ibin<<"] : " << Center[ibin] <<"\t TotalCenter : " << TotalCenter << endl;
  }

  // calculate Total eigenvector
  double TotalEigenVector_18Bin[Nev] = 0.;
  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {
    for(int iev(0); iev<Nev; iev++)
    {
      EigenVector_18Bin[ibin][iev] = EigenVector_18Bin[ibin][iev] * Correction; // correction factor by dimitri 1.074
      TotalEigenVector_18Bin[iev] += EigenVector_18Bin[ibin][iev];
      //cout << "TotalEigenVector["<<iev<<"] : " << TotalEigenVector[iev] << endl;
    }
  }

  // Normalized EigenVector and Center value
  double NormDiff_EigenVector_18Bin[Nbin_ZpT][Nev];
  double NormDiff_Center[Nbin_ZpT] = {0.,};
  double NormDiff_StatErr[Nbin_ZpT] = {0.,};
  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {
    for(int iev(0); iev<Nev; iev++)
    {
      NormDiff_EigenVector_18Bin[ibin][iev] = EigenVector_18Bin[ibin][iev] / TotalEigenVector_18Bin[iev] / ZptBinWidth[ibin] ;
//      cout << "Norm_EigenVector["<<ibin <<"]["<<iev<<"] : " << NormDiff_EigenVector[ibin][iev] << endl;
    }
    NormDiff_Center[ibin] = Center[ibin] / TotalCenter / ZptBinWidth[ibin];
//    cout << "Norm_Center[ibin]" << NormDiff_Center[ibin]<< endl;
  }
 
  NormToyErr(Center,StatErr,NormDiff_StatErr);

  // Xsec level PDF Error 
  double pluserr[Nbin_ZpT]={0.,};
  double minuserr[Nbin_ZpT]={0.,};
  // Calculate PDF uncertainty
  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {   
    for(int i=0; i<Nev/2;i++)
    {   
      if (  (EigenVector_18Bin[ibin][2*i] - Center[ibin]) <0 && (EigenVector_18Bin[ibin][2*i+1] - Center[ibin]) <0 )
        pluserr[ibin] = pluserr[ibin] ;
      else
        pluserr[ibin] = pluserr[ibin] + TMath::Power( (TMath::Max(EigenVector_18Bin[ibin][2*i] - Center[ibin], EigenVector_18Bin[ibin][2*i+1]-Center[ibin])) ,2);
   
      if (  (Center[ibin] - EigenVector_18Bin[ibin][2*i]) <0 && (Center[ibin] - EigenVector_18Bin[ibin][2*i+1]) <0 )
        minuserr[ibin] = minuserr[ibin];
      else
        minuserr[ibin] = minuserr[ibin] + TMath::Power( (TMath::Max(Center[ibin] - EigenVector_18Bin[ibin][2*i], Center[ibin] - EigenVector_18Bin[ibin][2*i+1])) ,2);
    }   
  }


  // Normalized Xsec level calculation
  double NormDiff_pluserr[Nbin_ZpT];
  double NormDiff_minuserr[Nbin_ZpT];
  // Calculate PDF uncertainty
  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {   
    for(int i=0; i<Nev/2;i++)
    {   
      if (  (NormDiff_EigenVector_18Bin[ibin][2*i] - NormDiff_Center[ibin]) <0 && (NormDiff_EigenVector_18Bin[ibin][2*i+1] - NormDiff_Center[ibin]) <0 )
        NormDiff_pluserr[ibin] = NormDiff_pluserr[ibin] ;
      else
        NormDiff_pluserr[ibin] = NormDiff_pluserr[ibin] + TMath::Power( (TMath::Max(NormDiff_EigenVector_18Bin[ibin][2*i] - NormDiff_Center[ibin], NormDiff_EigenVector_18Bin[ibin][2*i+1]-NormDiff_Center[ibin])) ,2);
   
      if (  (NormDiff_Center[ibin] - NormDiff_EigenVector_18Bin[ibin][2*i]) <0 && (NormDiff_Center[ibin] - NormDiff_EigenVector_18Bin[ibin][2*i+1]) <0 )
        NormDiff_minuserr[ibin] = NormDiff_minuserr[ibin];
      else
        NormDiff_minuserr[ibin] = NormDiff_minuserr[ibin] + TMath::Power( (TMath::Max(NormDiff_Center[ibin] - NormDiff_EigenVector_18Bin[ibin][2*i], NormDiff_Center[ibin] - NormDiff_EigenVector_18Bin[ibin][2*i+1])) ,2);
   
    
    }   
  }

   // make PDF uncer
  double PDFsystPlus[Nbin_ZpT];
  double PDFsystMinus[Nbin_ZpT];
  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {
    PDFsystPlus[ibin] = TMath::Sqrt(pluserr[ibin]);
    PDFsystMinus[ibin] = TMath::Sqrt(minuserr[ibin]);
    cout << Form("PDFsystPlus[%d] : %.8f, %.2f % \t PDFsystMinus[%d] : %.8f, %.2f % ",ibin,PDFsystPlus[ibin], PDFsystPlus[ibin] / Center[ibin] *100, ibin, PDFsystMinus[ibin], PDFsystMinus[ibin]/Center[ibin]*100) << endl;;
  }

  // make Normalized PDF uncer
  double NormDiff_PDFsystPlus[Nbin_ZpT];
  double NormDiff_PDFsystMinus[Nbin_ZpT];
  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {
    NormDiff_PDFsystPlus[ibin] = TMath::Sqrt(NormDiff_pluserr[ibin]);
    NormDiff_PDFsystMinus[ibin] = TMath::Sqrt(NormDiff_minuserr[ibin]);
    
    cout << Form("NormDiff_PDFsystPlus[%d] : %.8f, %.4f % \t Minus[%d] : %.8f, %.4f % \t max : %.4f % ",ibin,NormDiff_PDFsystPlus[ibin], NormDiff_PDFsystPlus[ibin] / NormDiff_Center[ibin] *100, ibin, NormDiff_PDFsystMinus[ibin], NormDiff_PDFsystMinus[ibin]/NormDiff_Center[ibin]*100, TMath::Max(fabs(NormDiff_PDFsystPlus[ibin] / NormDiff_Center[ibin] *100),fabs(NormDiff_PDFsystMinus[ibin]/NormDiff_Center[ibin]*100))) << endl;;
  }


  for(int ibin(0); ibin<Nbin_ZpT; ibin++)
  {
    cout << Form("NormDiff Xsec[%d] : %.8f \t +- %.8f",ibin,NormDiff_Center[ibin], NormDiff_StatErr[ibin]) << endl;
  }



}

void NormToyErr(double *_Xsec, double *_Err, double *_Norm_Err)
{
  int Ntoy = 100000;

  double mx2[Nbin_ZpT] = {0.,};
  double m2x[Nbin_ZpT] = {0.,};
  double temp[Nbin_ZpT] = {0.,};

  // Check _Xsec and _StatErr argument is correct 
//  for(int j=0; j<12; ++j)
//  {
    //cout << Form("_Xsec : %.4f \t _StatErr : %.4f",_Xsec[j],_StatErr[j]) << endl;
//  }

  // Here Calculate Normalized Stat Error by using toy
  for(int e=0; e<Ntoy; ++e)
  {
    for(int i=0; i<Nbin_ZpT; ++i)
    {
      	temp[i]=gRandom->Gaus(_Xsec[i], _Err[i]);
    }
    double sum=0;
    for(int i=0; i<Nbin_ZpT; ++i) sum+=temp[i];
    for(int i=0; i<Nbin_ZpT; ++i)
    {
      temp[i]/=sum;
      m2x[i]+=temp[i];
      mx2[i]+=temp[i]*temp[i];
    }
  }

  double rms[Nbin];
  for(int i=0; i<Nbin_ZpT; ++i)
  {
    m2x[i]/=Ntoy;
    mx2[i]/=Ntoy;
    rms[i] = sqrt(mx2[i]-m2x[i]*m2x[i]);
  }
 
  // return value
  for(int i=0; i<Nbin_ZpT; ++i)
  {
    _Norm_Err[i] = rms[i]/ZptBinWidth[i]; // return rms to NormDiff_Err
  }
  for(int i=0; i<Nbin_ZpT; ++i)
  {
    //cout << Form("rms = %.8f",rms[i]) << endl;
  }
}

