#include <stdio>

void CalcSyst()
{
  // PDF syst in % unit
  double PDFsyst[18];
  PDFsyst[0] = 2.85; 
  PDFsyst[1] = 8.22; 
  PDFsyst[2] = 3.31;
  PDFsyst[3] = 1.71;
  PDFsyst[4] = 1.19;
  PDFsyst[5] = 2.56;
  PDFsyst[6] = 3.20;
  PDFsyst[7] = 3.55;
  PDFsyst[8] = 2.91;
  PDFsyst[9] = 4.16;
  PDFsyst[10] = 4.59;
  PDFsyst[11] = 4.94;
  PDFsyst[12] = 5.56;
  PDFsyst[13] = 5.85;
  PDFsyst[14] = 6.75;
  PDFsyst[15] = 7.63;
  PDFsyst[16] = 8.55;
  PDFsyst[17] = 10.33;

  // Scale syst in % unit
  double Scalesyst[18];
  Scalesyst[0] = 3.86 ; 
  Scalesyst[1] = 2.94 ; 
  Scalesyst[2] = 3.66 ;
  Scalesyst[3] = 3.71 ;
  Scalesyst[4] = 3.42 ;
  Scalesyst[5] = 2.92 ;
  Scalesyst[6] = 0.58 ;
  Scalesyst[7] = 3.24 ;
  Scalesyst[8] = 7.06 ;
  Scalesyst[9] = 9.03 ;
  Scalesyst[10] =7.77 ;
  Scalesyst[11] =4.75 ;
  Scalesyst[12] =5.74 ;
  Scalesyst[13] =6.67 ;
  Scalesyst[14] =9.50 ;
  Scalesyst[15] =10.18;
  Scalesyst[16] =10.30;
  Scalesyst[17] =10.42;

  // Calc total syst in %
  double totalsyst[18] = {0.,};
  for(int i(0); i<18 ; i++)
  {
    totalsyst[i] = sqrt(PDFsyst[i]*PDFsyst[i] + Scalesyst[i]*Scalesyst[i]);
    cout << "totalsyst : " << totalsyst[i] << endl;
  }


  double normdiffxsec[18] = {
		 0.03335570, 	
                 0.05909690, 	
                 0.05080043, 	
                 0.04145975, 	
                 0.03394089, 	
                 0.02773926, 	
                 0.02220695, 	
                 0.01779636, 	
                 0.01131978, 	
                 0.00601003, 	
                 0.00357902, 	
                 0.00188844, 	
                 0.00085224, 	
                 0.00041736, 	
                 0.00018092, 	
                 0.00006654, 	
                 0.00002372, 	
                 0.00000218
  };

  double normdifferr[18] = {0.,};
  cout.precision(4);
  for(int j(0); j<18; j++)
  {
    normdifferr[j] = normdiffxsec[j] * totalsyst[j] / 100.0 ;
    cout << "normdifferr : " << normdifferr[j] << endl;
  }




}

