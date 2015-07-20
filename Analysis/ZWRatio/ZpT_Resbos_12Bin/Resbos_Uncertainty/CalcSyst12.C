#include <stdio>

void CalcSyst12()
{
  double PDFsyst[12];
  PDFsyst[0] = 3.70; 
  PDFsyst[1] = 2.64; 
  PDFsyst[2] = 3.28;
  PDFsyst[3] = 1.84;
  PDFsyst[4] = 0.87;
  PDFsyst[5] = 1.18;
  PDFsyst[6] = 1.45;
  PDFsyst[7] = 2.16;
  PDFsyst[8] = 3.35;
  PDFsyst[9] = 4.35;
  PDFsyst[10] =5.36;
  PDFsyst[11] =7.26;

  double Scalesyst[12];
  Scalesyst[0] = 2.68  ; 
  Scalesyst[1] = 3.58  ; 
  Scalesyst[2] = 1.68  ;
  Scalesyst[3] = 5.98  ;
  Scalesyst[4] = 9.03  ;
  Scalesyst[5] = 7.77  ;
  Scalesyst[6] = 4.75  ;
  Scalesyst[7] = 5.93  ;
  Scalesyst[8] = 9.50  ;
  Scalesyst[9] = 10.18 ;
  Scalesyst[10] =10.30 ;
  Scalesyst[11] =10.42 ;

  double totalsyst[12] = {0.,};
  for(int i(0); i<12 ; i++)
  {
    totalsyst[i] = sqrt(PDFsyst[i]*PDFsyst[i] + Scalesyst[i]*Scalesyst[i]);
    cout << "totalsyst : " << totalsyst[i] << endl;
  }

  double normdiffxsec[12] = {
      0.03941847,
      0.04359641,
      0.02712281,
      0.01397478,
      0.00619133,
      0.00369245,
      0.00194821,
      0.00065490,
      0.00018665,
      0.00006864,
      0.00002447,
      0.00000225
  };

  double normdifferr[12] = {0.,};
  cout.precision(4);
  for(int j(0); j<12; j++)
  {
    normdifferr[j] = normdiffxsec[j] * totalsyst[j] / 100.0 ;
    cout << "normdifferr: " << normdifferr[j] << endl;
  }



}

