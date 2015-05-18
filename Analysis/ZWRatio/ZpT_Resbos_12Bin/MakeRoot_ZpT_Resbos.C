#include <iostream>

int MakeRoot_ZpT_Resbos()
{
  double Zpt18Bins[19] = {0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,30.0,40.0,50.0,70.0,90.0,110.0,150.0,190.0,250.0,600.0};
  double Zpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

  // Resbos Normalized differential Xsec 
  double NormDiffXsec_18bin[18] = {0.02957739,  
                             0.05620224,
                             0.04961592,
                             0.04050349,
                             0.03325381,
                             0.02725463,
                             0.02253389,
                             0.01862596,
                             0.01237984,
                             0.00664593,
                             0.00393008,
                             0.00200304,
                             0.00085508,
                             0.00041032,
                             0.00017656,
                             0.00006332,
                             0.00002423,
                             0.00000115};

  //Make Resbos Normalized Xsec
  double NormXsec_18bin[18];
  for(int i(0);i<18;i++)
  {
    NormXsec_18bin[i] = NormDiffXsec_18bin[i] * (Zpt18Bins[i+1]-Zpt18Bins[i]) ;  // multiply bin width
    cout << "NormXsec_18bin : " << NormXsec_18bin[i] << endl;
  }

  //Merge bins 18 -> 12
  double NormXsec_12bin[12];
  NormXsec_12bin[0] = NormXsec_18bin[0]	+ NormXsec_18bin[1] + NormXsec_18bin[2]; // 0.0 ~ 7.5
  NormXsec_12bin[1] = NormXsec_18bin[3] + NormXsec_18bin[4]; 			// 7.5 ~ 12.5
  NormXsec_12bin[2] = NormXsec_18bin[5] + NormXsec_18bin[6]; 			// 12.5 ~ 17.5
  NormXsec_12bin[3] = NormXsec_18bin[7] + NormXsec_18bin[8]; 			// 17.5 ~ 30
  NormXsec_12bin[4] = NormXsec_18bin[9];
  NormXsec_12bin[5] = NormXsec_18bin[10];
  NormXsec_12bin[6] = NormXsec_18bin[11];
  NormXsec_12bin[7] = NormXsec_18bin[12] + NormXsec_18bin[13];
  NormXsec_12bin[8] = NormXsec_18bin[14];
  NormXsec_12bin[9] = NormXsec_18bin[15];
  NormXsec_12bin[10] = NormXsec_18bin[16];
  NormXsec_12bin[11] = NormXsec_18bin[17];

  // Make Norm to NormDiff
  double NormDiffXsec_12bin[12];
  for(int i(0);i<12;i++)
  {
    NormDiffXsec_12bin[i] = NormXsec_12bin[i] / (Zpt12Bins[i+1]-Zpt12Bins[i]); // divide bin width
    cout << "NormDiffXsec_12bin : " << NormDiffXsec_12bin[i] << endl;
  }


  TH1D* NormDiffXsec = new TH1D("NormDiffXsec","",12,0,12);
  for(int i(0); i<12;i++)
  {
    NormDiffXsec->SetBinContent(i+1,NormDiffXsec_12bin[i]);
  }

  TFile* outfile = new TFile("./root/ZToMuMu_Resbos.root","recreate");
  NormDiffXsec->Write();


  return 0;
}

