#include <iostream.h>
#include <fstream>
double ZpTbin[19] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};

void ScaleSystCalculator_FEWZ18Bin()
{
  ifstream FILE_ScaleUp("./FEWZ_ScaleUp_18bin.txt");

  double XsecUp[18] = {0.,};
  double TotalXsecUp = 0;
  double NormDiff_XsecUp[18] = {0.,};
  double ScaleErr[18] = {0.,};
  TString tmp;

  //Nominal Norm Xsec from ZpT code
  double y[18];
  y[0] = -0.04927787; 
  y[1] = 0.11569336 ; 
  y[2] = 0.07234848 ;
  y[3] = 0.04951937 ;
  y[4] = 0.03593429 ;
  y[5] = 0.02730082 ;
  y[6] = 0.02129024 ;
  y[7] = 0.01716449 ;
  y[8] = 0.01089575 ;
  y[9] = 0.00582537 ;
  y[10] = 0.00348453;
  y[11] = 0.00184388;
  y[12] = 0.00082675;
  y[13] = 0.00040203;
  y[14] = 0.00017418;
  y[15] = 0.00006305;
  y[16] = 0.00002236;
  y[17] = 0.00000208;

  for(int ibin(0); ibin<18; ibin++)
  {
    FILE_ScaleUp >> tmp >> XsecUp[ibin];
    TotalXsecUp += XsecUp[ibin] ;
    cout << XsecUp[ibin] << endl;
  }
  cout << TotalXsecUp << endl;

  //Make NormDiff Xsec scale Up case
  for(int ibin(0); ibin<18; ibin++)
  {
    NormDiff_XsecUp[ibin] = XsecUp[ibin] / TotalXsecUp /(ZpTbin[ibin+1] - ZpTbin[ibin]) ;
    cout << "NormDiff_XsecUp : " << NormDiff_XsecUp[ibin] << endl;
  }

  //Calc Scale syst
  for(int ibin(0); ibin<18; ibin++)
  {
    ScaleErr[ibin] = fabs(y[ibin]-NormDiff_XsecUp[ibin]);
    cout << "ScaleErr : " << ScaleErr[ibin] <<endl;
  }

}
