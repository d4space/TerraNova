#include <stdio>
#include <iostream>

const double Zpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
const double BinWidth[12] ={7.5-0, 12.5-7.5,17.5-12.5, 30.0-17.5, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

static const int NB = 12;

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

int MakeRoot_ZpT_FEWZ()
{
  // cross-section and errors
  double Xsec[12];
  double Xsec_up[12];
  double Xsec_down[12];
  
  double StatErr[12];
  double StatErr_percent[12];

  double PDFErrP[12];
  double PDFErrM[12];
  double PDFErr[12];
  double PDFErr_percent[12];

  double ScaleErrP[12];
  double ScaleErrM[12];
  double ScaleErr[12];
  double ScaleErr_percent[12];

  double TotalUncer[12];
  double TotalUncer_percent[12];
 
  // Normalized cross-section and errors
  double Norm_Xsec[12];
  double Norm_TotalUncer[12];

  // Normalized Differential Xsec and errors
  double NormDiff_Xsec[12];
  double NormDiff_TotalUncer[12];

    Xsec[0] =  138.274;  
    Xsec[1] =  103.757;
    Xsec[2] =  59.4709;
    Xsec[3] =  68.4585;
    Xsec[4] =  26.5977;
    Xsec[5] =  15.9503;
    Xsec[6] =  16.6003;
    Xsec[7] =  10.8541;  
    Xsec[8] =  2.98891;
    Xsec[9] =  1.03651;
    Xsec[10] = 0.551997; 
    Xsec[11] = 0.281856;

    Xsec_up[0] =  142.959;  
    Xsec_up[1] =    104.9;
    Xsec_up[2] =  58.6561;
    Xsec_up[3] =  67.1108;
    Xsec_up[4] =  25.7478;
    Xsec_up[5] =  15.2494;
    Xsec_up[6] =  15.8542;
    Xsec_up[7] =  10.3382;  
    Xsec_up[8] =  2.79667;
    Xsec_up[9] =  0.962229;
    Xsec_up[10] = 0.495937; 
    Xsec_up[11] = 0.257324;

    Xsec_down[0] =  135.624;  
    Xsec_down[1] =  102.489;
    Xsec_down[2] =  59.4373;
    Xsec_down[3] =  69.8794;
    Xsec_down[4] =  27.6927;
    Xsec_down[5] =  16.6091;
    Xsec_down[6] =  17.4425;
    Xsec_down[7] =  11.5709;  
    Xsec_down[8] =  3.15182;
    Xsec_down[9] =  1.1378;
    Xsec_down[10] = 0.583274; 
    Xsec_down[11] = 0.314832;

    StatErr[0]	=  0.991013;  
    StatErr[1]	=  0.347511;
    StatErr[2]	=  0.197072;
    StatErr[3]	=  0.139355;
    StatErr[4]	= 0.0763465;
    StatErr[5]	= 0.0709948;
    StatErr[6]	= 0.0728004;
    StatErr[7]	= 0.0373067;
    StatErr[8]	= 0.0217453;
    StatErr[9]	= 0.0128423;
    StatErr[10]	= 0.0113465; 
    StatErr[11]	= 0.0180237;

    PDFErrP[0]	= 14.0653;   
    PDFErrP[1]	= 7.32285;  
    PDFErrP[2]	=  3.8426;
    PDFErrP[3]	=  4.1238;
    PDFErrP[4]	= 1.42682;
    PDFErrP[5]	= 0.77231;
    PDFErrP[6]	= 0.686723;
    PDFErrP[7]	= 0.306618;
    PDFErrP[8]	= 0.0488802;
    PDFErrP[9]	= 0.0184541;
    PDFErrP[10]	= 0.0124737;
    PDFErrP[11]	= 0.00830217;

    PDFErrM[0]	=  10.0929; 
    PDFErrM[1]	=  3.83243; 
    PDFErrM[2]	=  1.92563; 
    PDFErrM[3]	=  1.90581; 
    PDFErrM[4]	= 0.605502; 
    PDFErrM[5]	= 0.310734;
    PDFErrM[6]	= 0.267817;
    PDFErrM[7]	= 0.140892;
    PDFErrM[8]	= 0.0480252;
    PDFErrM[9]	= 0.0236329;
    PDFErrM[10]	= 0.0199344;
    PDFErrM[11]	= 0.0149747;

  // Error Calculation 
  for(int i(0); i<12; i++)
  {
    StatErr_percent[i] = StatErr[i] / Xsec[i] *100 ; // stat error

    PDFErr[i] = TMath::Max(PDFErrP[i],PDFErrM[i]); // pdf error
    PDFErr_percent[i] = PDFErr[i] / Xsec[i] * 100;

    ScaleErrP[i] = TMath::Abs(Xsec[i] - Xsec_up[i]);
    ScaleErrM[i] = TMath::Abs(Xsec[i] - Xsec_down[i]);
    ScaleErr[i] = TMath::Max(ScaleErrP[i],ScaleErrM[i]); // scale error
    ScaleErr_percent[i] = ScaleErr[i] / Xsec[i] * 100;
    
    TotalUncer[i] = sqrt(StatErr[i]*StatErr[i] + PDFErr[i]*PDFErr[i] + ScaleErr[i]*ScaleErr[i]); // total uncertainty
    TotalUncer_percent[i] = TotalUncer[i] / Xsec[i] * 100; 
    
    //cout << Form("PDFErr + : %.8f PDFErr- : %.8f PDFErr : %.8f",PDFErrP[i],PDFErrM[i],PDFErr[i]) << endl;
    //cout << Form("ScaleErr + : %.8f ScaleErr- : %.8f ScaleErr : %.8f",ScaleErrP[i],ScaleErrM[i],ScaleErr[i]) << endl;
    cout << Form("Xsec : %f TotalUncer : %.8f",Xsec[i],TotalUncer[i] ) << endl;
  }

  // Calculate Normalized Xsec and Total Uncertainty
  for(int i(0); i<12; i++)
  {
    double TotalXsec;
    TotalXsec += Xsec[i];
  }
 // for(int i(0); i<12; i++)
 // {
 //   Norm_Xsec[i] = Xsec[i] / TotalXsec;
 // }
  
  // Calculate Normalized Xsec and Total uncertainty by using propagation formala ??? 
  ErrPropaNormXsec(Xsec, TotalUncer, Norm_Xsec, Norm_TotalUncer); 
  for(int i(0); i<12; i++)
  {
    cout << "Norm_Xsec : " << Norm_Xsec[i] << "  Norm_TotalUncer : " << Norm_TotalUncer[i] << endl;
  }
  
  // Calculate Normalized Differential Xsec and Total Uncertainty
  for(int i(0); i<12; i++)
  {
    NormDiff_Xsec[i] = Norm_Xsec[i] / BinWidth[i] ; 
    NormDiff_TotalUncer[i] = Norm_TotalUncer[i] / BinWidth[i] ; 
    cout << Form("NormDiff_Xsec : %.8f TotalUncer : %.8f", NormDiff_Xsec[i], NormDiff_TotalUncer[i]) << endl;
  }
  
  // Save in histogram
  TH1D* hxsec = new TH1D("hxsec","hxsec",12,Zpt12Bins);
  TH1D* hxsec_Norm = new TH1D("hxsec_Norm","hxsec_Norm",12,Zpt12Bins);
  TH1D* hxsec_NormDiff = new TH1D("hxsec_NormDiff","hxsec_NormDiff",12,Zpt12Bins);
  
  TH1D* hStatErr_percent = new TH1D("StatErr_percent","StatErr_percent",12,Zpt12Bins); 
  TH1D* hPDFErr_percent = new TH1D("PDFErr_percent","PDFErr_percent",12,Zpt12Bins); 
  TH1D* hScaleErr_percent = new TH1D("ScaleErr_percent","ScaleErr_percent",12,Zpt12Bins);

  for(int i(0); i<12; i++)
  {
    hxsec->SetBinContent(i+1,Xsec[i]);
    hxsec->SetBinError(i+1,TotalUncer[i]);

    hxsec_Norm->SetBinContent(i+1,Norm_Xsec[i]);
    hxsec_Norm->SetBinError(i+1,Norm_TotalUncer[i]);
    
    hxsec_NormDiff->SetBinContent(i+1,NormDiff_Xsec[i]);
    hxsec_NormDiff->SetBinError(i+1,NormDiff_TotalUncer[i]);
    
    hStatErr_percent->SetBinError(i+1,StatErr_percent[i]);
    hPDFErr_percent->SetBinError(i+1,PDFErr_percent[i]);
    hScaleErr_percent->SetBinError(i+1,ScaleErr_percent[i]);
  }

  TFile* f_out = new TFile("./root/ZpT_FEWZ.root","recreate");
  hxsec->Write();
  hxsec_Norm->Write();
  hxsec_NormDiff->Write();
  hStatErr_percent->Write();
  hPDFErr_percent->Write();
  hScaleErr_percent->Write();

  return 0;

}

