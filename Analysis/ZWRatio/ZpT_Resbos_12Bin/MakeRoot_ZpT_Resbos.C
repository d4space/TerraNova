#include <iostream>

int MakeRoot_ZpT_Resbos()
{
  // Resbos Normalized differential Xsec 
  double normdiffxsec[12] = {
 0.04777023,  
 0.03770032,
 0.02497310,
 0.01261509,
 0.00601003,
 0.00357902,
 0.00188844,
 0.00063480,
 0.00018092,
 0.00006654,
 0.00002372,
 0.00000218
  };

  double normdifferr[12] = {
	0.001801  , 
	0.001939  ,
	0.0009995 ,
	0.0008744 ,
	0.0005617 ,
	0.0002902 ,
	9.676e-05 ,
	4.133e-05 ,
	1.88e-05  ,
	7.599e-06 ,
	2.841e-06 ,
	2.857e-07 
  };
  
  TH1D* NormDiffXsec = new TH1D("NormDiffXsec","",12,0,12);
  for(int i(0); i<12;i++)
  {
    NormDiffXsec->SetBinContent(i+1,normdiffxsec[i]);
    NormDiffXsec->SetBinError(i+1,normdifferr[i]);
  }

  TFile* outfile = new TFile("./root/ZToMuMu_Resbos.root","recreate");
  NormDiffXsec->Write();


  return 0;
}

