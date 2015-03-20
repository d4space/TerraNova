{
#include "./Utils/const.h"
  //*

  
    //cout << fixed << setprecision(10);
    cout << fixed << setprecision(8);
  
  const int nBinsZpt=19;
  double ZptBins[19]   = {0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,30.0,40.0,50.0,70.0,90.0,110.0,150.0,190.0,250.0,600.0}; 
  
  TH1D *hZptBins_LinScale   = new TH1D("hZptBins_LinScale","hZptBins_LinScale",nBinsZpt-1,ZptBins);hZptBins_LinScale->Sumw2();                     
  TH1D *hZptXsec_LinScale   = new TH1D("hZptXsec_LinScale","hZptXsec_LinScale",nBinsZpt-1,ZptBins);hZptXsec_LinScale->Sumw2();                     
  
  TH1D *hZptYield   = new TH1D("hZptYield","hZptYield",nBinsZpt-1,ZptBins);hZptYield->Sumw2();                     
 
  /////Zpt Yield after FSR corr
  double ZptYield[19]={0};
	ZptYield[1] =0.03338098     ;//682.757 ;   //1	   0.0-2.5    
	ZptYield[2] =0.05532895     ;//1131.67 ;   //2	   2.5-5.0    
	ZptYield[3] =0.05192021     ;//1061.95 ;   //3	   5.0-7.5    
	ZptYield[4] =0.03861802     ;//789.873 ;   //4	   7.5-10.0   
	ZptYield[5] =0.03550062     ;//726.111 ;   //5	   10.0-12.5  
	ZptYield[6] =0.02409441     ;//492.814 ;   //6	   12.5-15.0  
	ZptYield[7] =0.02251570     ;//460.524 ;   //7	   15.0-17.5  
	ZptYield[8] =0.01715371     ;//350.853 ;   //8	   17.5-20.0  
	ZptYield[9] =0.01175338     ;//961.589 ;   //9	   20.0-30.0  
	ZptYield[10]=0.00650825     ;//532.466 ;  //10    30.0-40.0  
	ZptYield[11]=0.00402494     ;//329.296 ;  //11    40.0-50.0  
	ZptYield[12]=0.00215579     ;//352.747 ;  //12    50.0-70.0  
	ZptYield[13]=0.00088889     ;//145.448 ;  //13    70.0-90.0  
	ZptYield[14]=0.00040997     ;//67.082  ;  //14    90.0-110.0 
	ZptYield[15]=0.00016520     ;//54.0626 ;  //15    110.0-150.0
	ZptYield[16]=0.00007652     ;//25.0404 ;  //16    150.0-190.0
	ZptYield[17]=0.00000898     ;//4.40732 ;  //17    190.0-250.0
	ZptYield[18]=0.00000444     ;//12.7016 ;  //18    250.0-600.0
	for(int i(1);i<nBinsZpt;i++)
	{
	  hZptYield->SetBinContent(i, ZptYield[i]);
	  cout<<""<<hZptYield->GetBinContent(i)<<endl;
	}


  /////Zpt Normalized stat error after FSR corr  =====================================================================
  double ZptStatErr[19]={0};
	ZptStatErr[1] =  0.00184508     ;  //  37.7383 ; 
	ZptStatErr[2] =  0.00253791     ;  //  51.9091 ; 
	ZptStatErr[3] =  0.00248749     ;  //  50.8778 ; 
	ZptStatErr[4] =  0.00223300     ;  //  45.6727 ; 
	ZptStatErr[5] =  0.00209782     ;  //  42.9076 ; 
	ZptStatErr[6] =  0.00181309     ;  //  37.0839 ; 
	ZptStatErr[7] =  0.00170830     ;  //  34.9406 ; 
	ZptStatErr[8] =  0.00152366     ;  //  31.1642 ; 
	ZptStatErr[9] =  0.00047966     ;  //  39.2429 ; 
	ZptStatErr[10]=  0.00035710     ;  //  29.216   ; 
	ZptStatErr[11]=  0.00028541     ;  //  23.3506  ; 
	ZptStatErr[12]=  0.00014035     ;  //  22.9646  ; 
	ZptStatErr[13]=  0.00009272     ;  //  15.171   ; 
	ZptStatErr[14]=  0.00006424     ;  //  10.5114  ; 
	ZptStatErr[15]=  0.00002765     ;  //  9.04859  ; 
	ZptStatErr[16]=  0.00001888     ;  //  6.1784   ; 
	ZptStatErr[17]=  0.00000618     ;  //  3.0343   ; 
	ZptStatErr[18]=  0.00000196     ;  //  5.60278  ; 
	///// Zpt stat error converting from NormalizedError to Error 
	 cout<<" 18 Bin Zpt converting from Number to % Error:\t"<<endl;
    
	for(int i(1);i<nBinsZpt;i++)
	{
	  ///Zpt 18 bin stat error in %  to check 18 bin stat in % from Note
	  //ZptStatErr_Percent[i]=(ZptStatErr[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	  
	  ZptStatErr[i]=(ZptStatErr[i]*100.)/hZptYield->GetBinContent(i);
	  cout<<ZptStatErr[i]<<endl;
	}
    cout << "=======ZptStatErr Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptStatErr[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptStatErr[i]  << endl;
    }
    cout << fixed << setprecision(8);

  /////Zpt Bkg syst ==============================================================================
  double ZptSystBkg_plus[19]={0};
	ZptSystBkg_plus[1] =0.00014512 ; 
	ZptSystBkg_plus[2] =0.00023114 ; 
	ZptSystBkg_plus[3] =0.00021415 ; 
	ZptSystBkg_plus[4] =0.00011173 ; 
	ZptSystBkg_plus[5] =0.00010321 ; 
	ZptSystBkg_plus[6] =0.00005561 ; 
	ZptSystBkg_plus[7] =0.00003265 ; 
	ZptSystBkg_plus[8] =0.00003759 ; 
	ZptSystBkg_plus[9] =0.00000169 ; 
	ZptSystBkg_plus[10]=0.00002375  ; 
	ZptSystBkg_plus[11]=0.00003118  ; 
	ZptSystBkg_plus[12]=0.00003287  ; 
	ZptSystBkg_plus[13]=0.00002372  ; 
	ZptSystBkg_plus[14]=0.00001423  ; 
	ZptSystBkg_plus[15]=0.00000578  ; 
	ZptSystBkg_plus[16]=0.00000152  ; 
	ZptSystBkg_plus[17]=0.00000055  ; 
	ZptSystBkg_plus[18]=0.00000009  ; 
	///// Zpt Bkg syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptSystBkg_plus[i]=(ZptSystBkg_plus[i])*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptSystBkg_plus[i]<<endl;
	}
    cout << "=======ZptSystBkg_plus Put in note horisontal=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystBkg_plus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystBkg_plus[i] << endl;
    }
    cout << fixed << setprecision(8);
    
  double ZptSystBkg_minus[19]={0};
	ZptSystBkg_minus[1] = 0.00014363; 
	ZptSystBkg_minus[2] = 0.00022876; 
	ZptSystBkg_minus[3] = 0.00021198; 
	ZptSystBkg_minus[4] = 0.00011058; 
	ZptSystBkg_minus[5] = 0.00010216; 
	ZptSystBkg_minus[6] = 0.00005504; 
	ZptSystBkg_minus[7] = 0.00003232; 
	ZptSystBkg_minus[8] = 0.00003720; 
	ZptSystBkg_minus[9] = 0.00000167; 
	ZptSystBkg_minus[10]= 0.00002399 ; 
	ZptSystBkg_minus[11]= 0.00003150 ; 
	ZptSystBkg_minus[12]= 0.00003321 ; 
	ZptSystBkg_minus[13]= 0.00002397 ; 
	ZptSystBkg_minus[14]= 0.00001438 ; 
	ZptSystBkg_minus[15]= 0.00000584 ; 
	ZptSystBkg_minus[16]= 0.00000153 ; 
	ZptSystBkg_minus[17]= 0.00000055 ; 
	ZptSystBkg_minus[18]= 0.00000009 ; 
	///// Zpt Bkg syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptSystBkg_minus[i]=(ZptSystBkg_minus[i])*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptSystBkg_minus[i]<<endl;
	}
    cout << "=======ZptSystBkg_minus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystBkg_minus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystBkg_minus[i]<< endl ;
    }
    cout << fixed << setprecision(8);
  
 
  /////Zpt Eff syst  =====================================================================
  double ZptSystEff_plus[19]={0};
	ZptSystEff_plus[1] = 0.000001941; 
	ZptSystEff_plus[2] = 0.000002089; 
	ZptSystEff_plus[3] = 0.000000319; 
	ZptSystEff_plus[4] = 0.000000908; 
	ZptSystEff_plus[5] = 0.000000698; 
	ZptSystEff_plus[6] = 0.000000612; 
	ZptSystEff_plus[7] = 0.000000661; 
	ZptSystEff_plus[8] = 0.000000691; 
	ZptSystEff_plus[9] = 0.000000108; 
	ZptSystEff_plus[10]= 0.000000248 ; 
	ZptSystEff_plus[11]= 0.000000021 ; 
	ZptSystEff_plus[12]= 0.000000075 ; 
	ZptSystEff_plus[13]= 0.000000014 ; 
	ZptSystEff_plus[14]= 0.000000017 ; 
	ZptSystEff_plus[15]= 0.000000001 ; 
	ZptSystEff_plus[16]= 0.000000006 ; 
	ZptSystEff_plus[17]= 0.000000001 ; 
	ZptSystEff_plus[18]= 0.000000000001 ; 
	///// Zpt Eff syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptSystEff_plus[i]=(ZptSystEff_plus[i])*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptSystEff_plus[i]<<endl;
	}
    cout << "=======ZptSystEff_plus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystEff_plus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystEff_plus[i]<< endl ;
    }
    cout << fixed << setprecision(8);

  double ZptSystEff_minus[19]={0};
	ZptSystEff_minus[1] = 0.000000593; 
	ZptSystEff_minus[2] = 0.000000889; 
	ZptSystEff_minus[3] = 0.000002373; 
	ZptSystEff_minus[4] = 0.000000590; 
	ZptSystEff_minus[5] = 0.000000118; 
	ZptSystEff_minus[6] = 0.000000229; 
	ZptSystEff_minus[7] = 0.000000033; 
	ZptSystEff_minus[8] = 0.000000018; 
	ZptSystEff_minus[9] = 0.000000237; 
	ZptSystEff_minus[10]= 0.000000002 ; 
	ZptSystEff_minus[11]= 0.000000043 ; 
	ZptSystEff_minus[12]= 0.000000021 ; 
	ZptSystEff_minus[13]= 0.000000011 ; 
	ZptSystEff_minus[14]= 0.000000003 ; 
	ZptSystEff_minus[15]= 0.000000004 ; 
	ZptSystEff_minus[16]= 0.000000004 ; 
	ZptSystEff_minus[17]= 0.000000001 ; 
	ZptSystEff_minus[18]= 0.000000000001 ; 
	///// Zpt Eff syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptSystEff_minus[i]=(ZptSystEff_minus[i])*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptSystEff_minus[i]<<endl;
	}
    cout << "=======ZptSystEff_minus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystEff_minus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystEff_minus[i] << endl;
    }
    cout << fixed << setprecision(8);
	
  /////Zpt FSR syst  =====================================================================
  double ZptSystFSRErr[19]={0};
	ZptSystFSRErr[1] = 0.00000851; 
	ZptSystFSRErr[2] = 0.00001015 ; 
	ZptSystFSRErr[3] = 0.00000583 ; 
	ZptSystFSRErr[4] = 0.00000278 ; 
	ZptSystFSRErr[5] = 0.00000458 ; 
	ZptSystFSRErr[6] = 0.00000755 ; 
	ZptSystFSRErr[7] = 0.00000487 ; 
	ZptSystFSRErr[8] = 0.00000601 ; 
	ZptSystFSRErr[9] = 0.00000191 ; 
	ZptSystFSRErr[10]= 0.00000009  ; 
	ZptSystFSRErr[11]= 0.00000053  ; 
	ZptSystFSRErr[12]= 0.00000051  ; 
	ZptSystFSRErr[13]= 0.00000039  ; 
	ZptSystFSRErr[14]= 0.00000004  ; 
	ZptSystFSRErr[15]= 0.00000022  ; 
	ZptSystFSRErr[16]= 0.00000014  ; 
	ZptSystFSRErr[17]= 0.00000003  ; 
	ZptSystFSRErr[18]= 0.00000001 ; 
	///// Zpt FSR syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptSystFSRErr[i]=(ZptSystFSRErr[i] )*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptSystFSRErr[i]<<endl;
	}
    cout << "=======ZptSystFSRErr Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystFSRErr[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptSystFSRErr[i] << endl;
    }
    cout << fixed << setprecision(8);

  /////Zpt ZptUnfoldSyst_doubleGaussin_plus  =====================================================================
  double ZptUnfoldSyst_doubleGaussin_plus[19]={0};
	ZptUnfoldSyst_doubleGaussin_plus[1] = 0.000006077  ; 
	ZptUnfoldSyst_doubleGaussin_plus[2] = 0.000008190  ; 
	ZptUnfoldSyst_doubleGaussin_plus[3] = 0.000004431  ; 
	ZptUnfoldSyst_doubleGaussin_plus[4] = 0.000003394  ; 
	ZptUnfoldSyst_doubleGaussin_plus[5] = 0.000001564  ; 
	ZptUnfoldSyst_doubleGaussin_plus[6] = 0.000001001  ; 
	ZptUnfoldSyst_doubleGaussin_plus[7] = 0.000000411  ; 
	ZptUnfoldSyst_doubleGaussin_plus[8] = 0.000001419  ; 
	ZptUnfoldSyst_doubleGaussin_plus[9] = 0.000000177  ; 
	ZptUnfoldSyst_doubleGaussin_plus[10]= 0.000000489   ; 
	ZptUnfoldSyst_doubleGaussin_plus[11]= 0.000000110   ; 
	ZptUnfoldSyst_doubleGaussin_plus[12]= 0.000000163   ; 
	ZptUnfoldSyst_doubleGaussin_plus[13]= 0.000000010   ; 
	ZptUnfoldSyst_doubleGaussin_plus[14]= 0.000000007   ; 
	ZptUnfoldSyst_doubleGaussin_plus[15]= 0.000000007   ; 
	ZptUnfoldSyst_doubleGaussin_plus[16]= 0.000000001   ; 
	ZptUnfoldSyst_doubleGaussin_plus[17]= 0.00000000004 ; 
	ZptUnfoldSyst_doubleGaussin_plus[18]= 0.00000000123  ; 
	///// Zpt UnfoldSyst_doubleGaussin syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptUnfoldSyst_doubleGaussin_plus[i]=(ZptUnfoldSyst_doubleGaussin_plus[i])*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptUnfoldSyst_doubleGaussin_plus[i]<<endl;
	}
    cout << "=======ZptUnfoldSyst_doubleGaussin_plus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldSyst_doubleGaussin_plus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldSyst_doubleGaussin_plus[i] << endl;
    }
    cout << fixed << setprecision(8);


  double ZptUnfoldSyst_doubleGaussin_minus[19]={0};
	ZptUnfoldSyst_doubleGaussin_minus[1] = 0.000000856 ; 
	ZptUnfoldSyst_doubleGaussin_minus[2] = 0.000001552 ; 
	ZptUnfoldSyst_doubleGaussin_minus[3] = 0.000003008 ; 
	ZptUnfoldSyst_doubleGaussin_minus[4] = 0.000002112 ; 
	ZptUnfoldSyst_doubleGaussin_minus[5] = 0.000000641 ; 
	ZptUnfoldSyst_doubleGaussin_minus[6] = 0.000001026 ; 
	ZptUnfoldSyst_doubleGaussin_minus[7] = 0.000000980 ; 
	ZptUnfoldSyst_doubleGaussin_minus[8] = 0.000000133 ; 
	ZptUnfoldSyst_doubleGaussin_minus[9] = 0.000000339 ; 
	ZptUnfoldSyst_doubleGaussin_minus[10]= 0.000000307  ; 
	ZptUnfoldSyst_doubleGaussin_minus[11]= 0.000000092  ; 
	ZptUnfoldSyst_doubleGaussin_minus[12]= 0.000000048  ; 
	ZptUnfoldSyst_doubleGaussin_minus[13]= 0.000000148  ; 
	ZptUnfoldSyst_doubleGaussin_minus[14]= 0.000000146  ; 
	ZptUnfoldSyst_doubleGaussin_minus[15]= 0.000000002  ; 
	ZptUnfoldSyst_doubleGaussin_minus[16]= 0.000000003  ; 
	ZptUnfoldSyst_doubleGaussin_minus[17]= 0.00000000693 ; 
	ZptUnfoldSyst_doubleGaussin_minus[18]= 0.00000000000001 ; 
	///// Zpt UnfoldSyst_doubleGaussin syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptUnfoldSyst_doubleGaussin_minus[i]=(ZptUnfoldSyst_doubleGaussin_minus[i])*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptUnfoldSyst_doubleGaussin_minus[i]<<endl;
	}
    cout << "=======ZptUnfoldSyst_doubleGaussin_minus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldSyst_doubleGaussin_minus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldSyst_doubleGaussin_minus[i] << endl;
    }
    cout << fixed << setprecision(8);
	
  /////Zpt Unfolding Scale variation syst  =====================================================================
  double ZptUnfoldingShapeSyst_plus[19]={0};
	ZptUnfoldingShapeSyst_plus[1] = 0.00000400   ; 
	ZptUnfoldingShapeSyst_plus[2] = 0.00000515   ; 
	ZptUnfoldingShapeSyst_plus[3] = 0.00000775   ; 
	ZptUnfoldingShapeSyst_plus[4] = 0.00000173   ; 
	ZptUnfoldingShapeSyst_plus[5] = 0.00000044   ; 
	ZptUnfoldingShapeSyst_plus[6] = 0.00000058   ; 
	ZptUnfoldingShapeSyst_plus[7] = 0.00000191   ; 
	ZptUnfoldingShapeSyst_plus[8] = 0.00000135   ; 
	ZptUnfoldingShapeSyst_plus[9] = 0.00000154   ; 
	ZptUnfoldingShapeSyst_plus[10]= 0.00000047    ; 
	ZptUnfoldingShapeSyst_plus[11]= 0.00000008    ; 
	ZptUnfoldingShapeSyst_plus[12]= 0.00000020    ; 
	ZptUnfoldingShapeSyst_plus[13]= 0.00000014    ; 
	ZptUnfoldingShapeSyst_plus[14]= 0.00000014    ; 
	ZptUnfoldingShapeSyst_plus[15]= 0.00000000    ; 
	ZptUnfoldingShapeSyst_plus[16]= 0.00000001    ; 
	ZptUnfoldingShapeSyst_plus[17]= 0.00000001100 ; 
	ZptUnfoldingShapeSyst_plus[18]= 0.00000000001  ; 
	///// Zpt Unfolding Scale syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptUnfoldingShapeSyst_plus[i]=(ZptUnfoldingShapeSyst_plus[i] )*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptUnfoldingShapeSyst_plus[i]<<endl;
	}
    cout << "=======ZptUnfoldingShapeSyst_plus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldingShapeSyst_plus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldingShapeSyst_plus[i] << endl;
    }
    cout << fixed << setprecision(8);
	

  double ZptUnfoldingShapeSyst_minus[19]={0};
	ZptUnfoldingShapeSyst_minus[1] = 0.00000281   ; 
	ZptUnfoldingShapeSyst_minus[2] = 0.00000930   ; 
	ZptUnfoldingShapeSyst_minus[3] = 0.00000364   ; 
	ZptUnfoldingShapeSyst_minus[4] = 0.00000196   ; 
	ZptUnfoldingShapeSyst_minus[5] = 0.00000383   ; 
	ZptUnfoldingShapeSyst_minus[6] = 0.00000031   ; 
	ZptUnfoldingShapeSyst_minus[7] = 0.00000454   ; 
	ZptUnfoldingShapeSyst_minus[8] = 0.00000171   ; 
	ZptUnfoldingShapeSyst_minus[9] = 0.00000088   ; 
	ZptUnfoldingShapeSyst_minus[10]= 0.00000062    ; 
	ZptUnfoldingShapeSyst_minus[11]= 0.00000045    ; 
	ZptUnfoldingShapeSyst_minus[12]= 0.00000016    ; 
	ZptUnfoldingShapeSyst_minus[13]= 0.00000002    ; 
	ZptUnfoldingShapeSyst_minus[14]= 0.00000000    ; 
	ZptUnfoldingShapeSyst_minus[15]= 0.00000008    ; 
	ZptUnfoldingShapeSyst_minus[16]= 0.00000001    ; 
	ZptUnfoldingShapeSyst_minus[17]= 0.00000000007 ; 
	ZptUnfoldingShapeSyst_minus[18]= 0.00000000136  ; 
	///// Zpt Unfolding Scale syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptUnfoldingShapeSyst_minus[i]=(ZptUnfoldingShapeSyst_minus[i] )*100./hZptYield->GetBinContent(i);
	  cout<<""<<ZptUnfoldingShapeSyst_minus[i]<<endl;
	}
    cout << "=======ZptUnfoldingShapeSyst_minus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldingShapeSyst_minus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptUnfoldingShapeSyst_minus[i] << endl;
    }
    cout << fixed << setprecision(8);
	
  //////// Zpt SystMadgraph 
  double ZptMadgraphSyst_plus[19]={0};
	ZptMadgraphSyst_plus[1] = 0.00090357; 
	ZptMadgraphSyst_plus[2] = 0.00073026; 
	ZptMadgraphSyst_plus[3] = 0.00014385; 
	ZptMadgraphSyst_plus[4] = 0.00050344; 
	ZptMadgraphSyst_plus[5] = 0.00050843; 
	ZptMadgraphSyst_plus[6] = 0.00055764; 
	ZptMadgraphSyst_plus[7] = 0.00029025; 
	ZptMadgraphSyst_plus[8] = 0.00028008; 
	ZptMadgraphSyst_plus[9] = 0.00004866; 
	ZptMadgraphSyst_plus[10]= 0.00003656 ; 
	ZptMadgraphSyst_plus[11]= 0.00004129 ; 
	ZptMadgraphSyst_plus[12]= 0.00000561 ; 
	ZptMadgraphSyst_plus[13]= 0.00000326 ; 
	ZptMadgraphSyst_plus[14]= 0.00000273 ; 
	ZptMadgraphSyst_plus[15]= 0.00000189 ; 
	ZptMadgraphSyst_plus[16]= 0.00000011 ; 
	ZptMadgraphSyst_plus[17]= 0.00000089 ; 
	ZptMadgraphSyst_plus[18]= 0.00000002  ; 
	///// Zpt Madgraph syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptMadgraphSyst_plus[i]=(ZptMadgraphSyst_plus[i])*100./hZptYield->GetBinContent(i);
	  
	  cout<<
             ZptMadgraphSyst_plus[i]<<"\t %"
	    <<endl;
	}
    cout << fixed << setprecision(8);
    
    cout << "=======ZptMadgraphSyst_plus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptMadgraphSyst_plus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<ZptMadgraphSyst_plus[i] << endl;
    }
    cout << fixed << setprecision(8);
  
    
  double ZptMadgraphSyst_minus[19]={0};
	ZptMadgraphSyst_minus[1] = 0.00090357; 
	ZptMadgraphSyst_minus[2] = 0.00073026; 
	ZptMadgraphSyst_minus[3] = 0.00014385; 
	ZptMadgraphSyst_minus[4] = 0.00050344; 
	ZptMadgraphSyst_minus[5] = 0.00050843; 
	ZptMadgraphSyst_minus[6] = 0.00055764; 
	ZptMadgraphSyst_minus[7] = 0.00029025; 
	ZptMadgraphSyst_minus[8] = 0.00028008; 
	ZptMadgraphSyst_minus[9] = 0.00004866; 
	ZptMadgraphSyst_minus[10]= 0.00003656 ; 
	ZptMadgraphSyst_minus[11]= 0.00004129 ; 
	ZptMadgraphSyst_minus[12]= 0.00000561 ; 
	ZptMadgraphSyst_minus[13]= 0.00000326 ; 
	ZptMadgraphSyst_minus[14]= 0.00000273 ; 
	ZptMadgraphSyst_minus[15]= 0.00000189 ; 
	ZptMadgraphSyst_minus[16]= 0.00000011 ; 
	ZptMadgraphSyst_minus[17]= 0.00000089 ; 
	ZptMadgraphSyst_minus[18]= 0.00000002  ; 
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptMadgraphSyst_minus[i]=(ZptMadgraphSyst_minus[i])*100./hZptYield->GetBinContent(i);
	}
	

  //////// Zpt Accept Syst 
  double ZptAcceptSyst_plus[19]={0};
	ZptAcceptSyst_plus[1] = 0.00013368; 
	ZptAcceptSyst_plus[2] = 0.00016758; 
	ZptAcceptSyst_plus[3] = 0.00006861; 
	ZptAcceptSyst_plus[4] = 0.00002677; 
	ZptAcceptSyst_plus[5] = 0.00001256; 
	ZptAcceptSyst_plus[6] = 0.00000900; 
	ZptAcceptSyst_plus[7] = 0.00002268; 
	ZptAcceptSyst_plus[8] = 0.00002440; 
	ZptAcceptSyst_plus[9] = 0.00001646; 
	ZptAcceptSyst_plus[10]= 0.00001719 ; 
	ZptAcceptSyst_plus[11]= 0.00001693 ; 
	ZptAcceptSyst_plus[12]= 0.00001264 ; 
	ZptAcceptSyst_plus[13]= 0.00000672 ; 
	ZptAcceptSyst_plus[14]= 0.00000352 ; 
	ZptAcceptSyst_plus[15]= 0.00000181 ; 
	ZptAcceptSyst_plus[16]= 0.00000085 ; 
	ZptAcceptSyst_plus[17]= 0.00000010 ; 
	ZptAcceptSyst_plus[18]= 0.00000005  ; 
	///// Zpt Accept syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptAcceptSyst_plus[i]=(ZptAcceptSyst_plus[i] )*100./hZptYield->GetBinContent(i);
	}
	
  //////// Zpt Accept Syst 
  double ZptAcceptSyst_minus[19]={0};
	ZptAcceptSyst_minus[1] = 0.00015581; 
	ZptAcceptSyst_minus[2] = 0.00020012; 
	ZptAcceptSyst_minus[3] = 0.00007948; 
	ZptAcceptSyst_minus[4] = 0.00003783; 
	ZptAcceptSyst_minus[5] = 0.00001998; 
	ZptAcceptSyst_minus[6] = 0.00000828; 
	ZptAcceptSyst_minus[7] = 0.00001917; 
	ZptAcceptSyst_minus[8] = 0.00002046; 
	ZptAcceptSyst_minus[9] = 0.00001741; 
	ZptAcceptSyst_minus[10]= 0.00001449 ; 
	ZptAcceptSyst_minus[11]= 0.00001380 ; 
	ZptAcceptSyst_minus[12]= 0.00000985 ; 
	ZptAcceptSyst_minus[13]= 0.00000515 ; 
	ZptAcceptSyst_minus[14]= 0.00000262 ; 
	ZptAcceptSyst_minus[15]= 0.00000135 ; 
	ZptAcceptSyst_minus[16]= 0.00000061 ; 
	ZptAcceptSyst_minus[17]= 0.00000008 ; 
	ZptAcceptSyst_minus[18]= 0.00000003  ; 
	///// Zpt Accept syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptAcceptSyst_minus[i]=(ZptAcceptSyst_minus[i] )*100./hZptYield->GetBinContent(i);
	}
	
    ///// Calculating Zpt total syst for 18 bins 
    double syst2plus[19]={0};  
    double syst2minus[19]={0};
    for(int i(1);i<nBinsZpt;i++)
    {
      syst2plus[i] = 
      TMath::Power(ZptSystBkg_plus[i],2)+
      TMath::Power(ZptUnfoldingShapeSyst_plus[i],2)+
      TMath::Power(ZptUnfoldSyst_doubleGaussin_plus[i],2)+
      TMath::Power(ZptSystEff_plus[i],2)+
      TMath::Power(ZptMadgraphSyst_plus[i],2)+
      TMath::Power(ZptSystFSRErr[i],2);
      //TMath::Power(ZptAcceptSyst_plus[i],2); //// Do not use this for Fid Volume Xsec 
    }
    for(int i(1);i<nBinsZpt;i++)
    {
      syst2minus[i] = 
      TMath::Power(ZptSystBkg_minus[i],2)+
      TMath::Power(ZptUnfoldingShapeSyst_minus[i],2)+
      TMath::Power(ZptUnfoldSyst_doubleGaussin_minus[i],2)+
      TMath::Power(ZptSystEff_minus[i],2)+
      TMath::Power(ZptMadgraphSyst_minus[i],2) +
      TMath::Power(ZptSystFSRErr[i],2);
      //TMath::Power(ZptAcceptSyst_minus[i],2);  //// Do not use this for Fid Volume Xsec 
    }

    double systplus[19]={0};
    double systminus[19]={0};
    cout<<"=================================================================="<<endl;
    cout<<"Zpt Total syst Errors 18 Bin as in Note in %:"<<endl;
    cout<<"Bin:\t"<<"systplus\t"<<"systminus"<<endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      systplus[i]  = TMath::Sqrt(syst2plus[i] );
      systminus[i] = TMath::Sqrt(syst2minus[i]);
      cout<<i<<"\t+"<< systplus[i]<<"\t-"<< systminus[i]<<endl;
    }
    cout << "=======systplus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<systplus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    cout << fixed << setprecision(8);
    
    cout << "=======systminus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<systminus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    cout << fixed << setprecision(8);
    
    
    double errorplus[19]={0};
    double errorminus[19]={0};
    cout<<"=================================================================="<<endl;
    cout<<"Zpt Total Uncertainty 18 Bin as in Note in %:"<<endl;
    cout<<"Bin:\t"<<"errorplus\t"<<"errorminus"<<endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      errorplus[i]= TMath::Sqrt( TMath::Power(ZptStatErr[i],2) + TMath::Power(systplus[i] ,2) );
      errorminus[i]= TMath::Sqrt( TMath::Power(ZptStatErr[i] ,2) + TMath::Power(systminus[i],2) );
      cout<<i<<"\t+"<<errorplus[i]<<"\t-"<<errorminus[i]<<endl;
    }
    
    
    
    cout << "=======Total uncertainty errorplus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<errorplus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    cout << fixed << setprecision(8);
    cout << "=======Total uncertainty errorminus Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<errorminus[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    cout << fixed << setprecision(8);
    
    cout << "=======Total uncertainty TMath::Max(errorplus,errorminus) Put in note=======" << endl;
    cout << fixed << setprecision(2);
    for(int i(1);i<nBinsZpt;i++)
    {
      cout<<
	" & "<<TMath::Max(errorplus[i],errorminus[i] );
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    cout << fixed << setprecision(8);


}
