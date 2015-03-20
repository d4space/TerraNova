{
#include "./Utils/const.h"
  //*

  
    //cout << fixed << setprecision(10);
    cout << fixed << setprecision(8);
  
  const int nBinsZpt=19;
  const int nBins12=13;
  double ZptBins[19]   = {0.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,30.0,40.0,50.0,70.0,90.0,110.0,150.0,190.0,250.0,600.0}; 
                   //Merge  Bin1    0.0,2.5,5.0,7.5 
		   //       Bin2    7.5,10.0,12.5 
		   //       Bin3    12.5,15.0,17.5 
		   //       Bin4    17.5,20.0,30.0 
		   //       Bin8    70.0,90.0,110.0   
  double Zpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

  double Wpt12Bins[13] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
  
  TH1D *hZptBins_LinScale   = new TH1D("hZptBins_LinScale","hZptBins_LinScale",nBinsZpt-1,ZptBins);hZptBins_LinScale->Sumw2();                     
  TH1D *hWptBins_LinScale   = new TH1D("hWptBins_LinScale","hWptBins_LinScale",nBins12-1,Wpt12Bins);hWptBins_LinScale->Sumw2(); 
 

  TH1D *hZptXsec_LinScale   = new TH1D("hZptXsec_LinScale","hZptXsec_LinScale",nBinsZpt-1,ZptBins);hZptXsec_LinScale->Sumw2();                     
  
  TH1D *hZptYield   = new TH1D("hZptYield","hZptYield",nBinsZpt-1,ZptBins);hZptYield->Sumw2();                     
  TH1D *hZptYield12Bin   = new TH1D("hZptYield12Bin","hZptYield12Bin",nBins12-1,Wpt12Bins);hZptYield12Bin->Sumw2();                     
  
  /////Zpt Yield after FSR corr
  double ZptYield[19]={0};
	ZptYield[1] =682.757 ;   //1	   0.0-2.5    
	ZptYield[2] =1131.67 ;   //2	   2.5-5.0    
	ZptYield[3] =1061.95 ;   //3	   5.0-7.5    
	ZptYield[4] =789.873 ;   //4	   7.5-10.0   
	ZptYield[5] =726.111 ;   //5	   10.0-12.5  
	ZptYield[6] =492.814 ;   //6	   12.5-15.0  
	ZptYield[7] =460.524 ;   //7	   15.0-17.5  
	ZptYield[8] =350.853 ;   //8	   17.5-20.0  
	ZptYield[9] =961.589 ;   //9	   20.0-30.0  
	ZptYield[10]=532.466 ;  //10    30.0-40.0  
	ZptYield[11]=329.296 ;  //11    40.0-50.0  
	ZptYield[12]=352.747 ;  //12    50.0-70.0  
	ZptYield[13]=145.448 ;  //13    70.0-90.0  
	ZptYield[14]=67.082  ;  //14    90.0-110.0 
	ZptYield[15]=54.0626 ;  //15    110.0-150.0
	ZptYield[16]=25.0404 ;  //16    150.0-190.0
	ZptYield[17]=4.40732 ;  //17    190.0-250.0
	ZptYield[18]=12.7016 ;  //18    250.0-600.0
	for(int i(1);i<nBinsZpt;i++)
	{
	  hZptYield->SetBinContent(i, ZptYield[i]);
	}

  ///// Making 12 Bin Zpt N of Yield
  double ZptYield12[13]={0};
     ZptYield12[1]=ZptYield[1]+ZptYield[2]+ZptYield[3];      // 0.0-7.5    
     ZptYield12[2]=ZptYield[4]+ZptYield[5];                  // 7.5-12.5   
     ZptYield12[3]=ZptYield[6]+ZptYield[7];                  // 12.5-17.5  
     ZptYield12[4]=ZptYield[8]+ZptYield[9];                  // 17.5-30.0  
     ZptYield12[5]=ZptYield[10];                             // 30.0-40.0  
     ZptYield12[6]=ZptYield[11];                             // 40.0-50.0  
     ZptYield12[7]=ZptYield[12];                             // 50.0-70.0  
     ZptYield12[8]=ZptYield[13]+ZptYield[14];                // 70.0-110.0 
     ZptYield12[9]=ZptYield[15];                             // 110.0-150.0
     ZptYield12[10]=ZptYield[16];                            // 150.0-190.0
     ZptYield12[11]=ZptYield[17];                            // 190.0-250.0
     ZptYield12[12]=ZptYield[18];                            // 250.0-600.0
	for(int i(1);i<nBins12;i++)
	{
	  hZptYield12Bin->SetBinContent(i, ZptYield12[i]);
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
	 cout<<" 18 Bin Zpt converting from NormalizedError to Error:\t"<<endl;
    
	double ZptStatErr_Percent[19]={0};
	for(int i(1);i<nBinsZpt;i++)
	{
	  ///Zpt 18 bin stat error in %  to check 18 bin stat in % from Note
	  ZptStatErr_Percent[i]=(ZptStatErr[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	  
	  
	  ZptStatErr[i]=(ZptStatErr[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) );
	  cout<<""<<ZptStatErr[i]<<endl;
	}
	for(int i(1);i<nBinsZpt;i++)
	{
	  //ZptStatErr[i]=ZptStatErr[i]*100./hZptYield->GetBinContent(i);
	  //cout<<"converting from to %:\t\t"<<ZptStatErr[i]<<endl;
	}

  double ZptStatErr12[13]={0};
     //cout<<" making 18 Bin stat error 12 bin stat error: sqrt(value1^2 + value2^2...)"<<endl;
     //ZptStatErr12[1] =TMath::Sqrt(TMath::Power(ZptStatErr[1],2)+ TMath::Power(ZptStatErr[2],2)+TMath::Power(ZptStatErr[3],2) );
     //ZptStatErr12[2] =TMath::Sqrt(TMath::Power(ZptStatErr[4],2)+ TMath::Power(ZptStatErr[5],2) );                	  // 7.5-12.5   
     //ZptStatErr12[3] =TMath::Sqrt(TMath::Power(ZptStatErr[6],2)+ TMath::Power(ZptStatErr[7],2) );                         // 12.5-17.5  
     //ZptStatErr12[4] =TMath::Sqrt(TMath::Power(ZptStatErr[8],2)+ TMath::Power(ZptStatErr[9],2) );                         // 17.5-30.0  
     //ZptStatErr12[5] =ZptStatErr[10];                           		     					  // 30.0-40.0  
     //ZptStatErr12[6] =ZptStatErr[11];                        				 		   	          // 40.0-50.0  
     //ZptStatErr12[7] =ZptStatErr[12];                        							          // 50.0-70.0  
     //ZptStatErr12[8] =TMath::Sqrt(TMath::Power(ZptStatErr[13],2)+ TMath::Power(ZptStatErr[14],2) );               	  // 70.0-110.0 
     //ZptStatErr12[9] =ZptStatErr[15];                            	        				          // 110.0-150.0
     //ZptStatErr12[10]=ZptStatErr[16];                                   	  					  // 150.0-190.0
     //ZptStatErr12[11]=ZptStatErr[17];             		            					          // 190.0-250.0
     //ZptStatErr12[12]=ZptStatErr[18];                       							          // 250.0-600.0

        /////   Zpt 12 Bin RMS numbers
	ZptStatErr12[1] =0.00105518;     
	ZptStatErr12[2] =0.00135422;
	ZptStatErr12[3] =0.00115714;
	ZptStatErr12[4] =0.00045697;
	ZptStatErr12[5] =0.00034962;
	ZptStatErr12[6] =0.00028205;
	ZptStatErr12[7] =0.00013820;
	ZptStatErr12[8] =0.00005546;
	ZptStatErr12[9] =0.00002776;
	ZptStatErr12[10]=0.00001897;
	ZptStatErr12[11]=0.00000615;
	ZptStatErr12[12]=0.00000196;	 
       
	for(int i(1);i<nBins12;i++)
	{
	  cout<<""<<ZptStatErr12[i]<<endl;
	}

	  cout<<"converting 12 bin stat error from number to % \t\t"<<endl;
	for(int i(1);i<nBins12;i++)
	{
	  //ZptStatErr12[i]=ZptStatErr12[i]*100./hZptYield12Bin->GetBinContent(i);
	  ZptStatErr12[i]=(ZptStatErr12[i]*hZptYield12Bin->Integral()*hWptBins_LinScale->GetXaxis()->GetBinWidth(i))*100./hZptYield12Bin->GetBinContent(i);
	  cout<<ZptStatErr12[i]<<endl;
	}

   
     
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
	  ZptSystBkg_plus[i]=(ZptSystBkg_plus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}


  double ZptSystBkg_plus12[13]={0};
    cout <<"ZptSystBkg_plus 1,2,3 \t"<<ZptSystBkg_plus[1]<<"\t"<<ZptSystBkg_plus[2]<<"\t"<<ZptSystBkg_plus[3]<<endl;
    if (ZptSystBkg_plus[1]>=ZptSystBkg_plus[2] && ZptSystBkg_plus[1]>=ZptSystBkg_plus[3])
       ZptSystBkg_plus12[1]=ZptSystBkg_plus[1];
     if (ZptSystBkg_plus[2]>=ZptSystBkg_plus[1] && ZptSystBkg_plus[2]>=ZptSystBkg_plus[3])
       ZptSystBkg_plus12[1]=ZptSystBkg_plus[2];
     if (ZptSystBkg_plus[3]>=ZptSystBkg_plus[1] && ZptSystBkg_plus[3]>=ZptSystBkg_plus[2])
       ZptSystBkg_plus12[1]=ZptSystBkg_plus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptSystBkg_plus 1,2,3 \t"<<ZptSystBkg_plus12[1]<<endl;
     ZptSystBkg_plus12[2] =TMath::Max(ZptSystBkg_plus[4],ZptSystBkg_plus[5]);               	      // 7.5-12.5   
     ZptSystBkg_plus12[3] =TMath::Max(ZptSystBkg_plus[6],ZptSystBkg_plus[7]);                        // 12.5-17.5  
     ZptSystBkg_plus12[4] =TMath::Max(ZptSystBkg_plus[8],ZptSystBkg_plus[9]);                        // 17.5-30.0  
     ZptSystBkg_plus12[5] =ZptSystBkg_plus[10];                           		              // 30.0-40.0  
     ZptSystBkg_plus12[6] =ZptSystBkg_plus[11];                         			              // 40.0-50.0  
     ZptSystBkg_plus12[7] =ZptSystBkg_plus[12];                        			      	      // 50.0-70.0  
     ZptSystBkg_plus12[8] =TMath::Max(ZptSystBkg_plus[13],ZptSystBkg_plus[14]);            	      // 70.0-110.0 
    cout <<"Biggest one of ZptSystBkg_plus 13,14 \t"<<ZptSystBkg_plus12[8]<<endl;
     ZptSystBkg_plus12[9] =ZptSystBkg_plus[15];                            	        	      // 110.0-150.0
     ZptSystBkg_plus12[10]=ZptSystBkg_plus[16];                                   	              // 150.0-190.0
     ZptSystBkg_plus12[11]=ZptSystBkg_plus[17];             		            		      // 190.0-250.0
     ZptSystBkg_plus12[12]=ZptSystBkg_plus[18];                       			    	      // 250.0-600.0

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
	  ZptSystBkg_minus[i]=(ZptSystBkg_minus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
  
  double ZptSystBkg_minus12[13]={0};
    cout <<"ZptSystBkg_minus 1,2,3 \t"<<ZptSystBkg_minus[1]<<"\t"<<ZptSystBkg_minus[2]<<"\t"<<ZptSystBkg_minus[3]<<endl;
    if (ZptSystBkg_minus[1]>=ZptSystBkg_minus[2] && ZptSystBkg_minus[1]>=ZptSystBkg_minus[3])
       ZptSystBkg_minus12[1]=ZptSystBkg_minus[1];
     if (ZptSystBkg_minus[2]>=ZptSystBkg_minus[1] && ZptSystBkg_minus[2]>=ZptSystBkg_minus[3])
       ZptSystBkg_minus12[1]=ZptSystBkg_minus[2];
     if (ZptSystBkg_minus[3]>=ZptSystBkg_minus[1] && ZptSystBkg_minus[3]>=ZptSystBkg_minus[2])
       ZptSystBkg_minus12[1]=ZptSystBkg_minus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptSystBkg_minus 1,2,3 \t"<<ZptSystBkg_minus12[1]<<endl;
     ZptSystBkg_minus12[2] =TMath::Max(ZptSystBkg_minus[4],ZptSystBkg_minus[5]);               	      // 7.5-12.5   
     ZptSystBkg_minus12[3] =TMath::Max(ZptSystBkg_minus[6],ZptSystBkg_minus[7]);                        // 12.5-17.5  
     ZptSystBkg_minus12[4] =TMath::Max(ZptSystBkg_minus[8],ZptSystBkg_minus[9]);                        // 17.5-30.0  
     ZptSystBkg_minus12[5] =ZptSystBkg_minus[10];                           		              // 30.0-40.0  
     ZptSystBkg_minus12[6] =ZptSystBkg_minus[11];                         			              // 40.0-50.0  
     ZptSystBkg_minus12[7] =ZptSystBkg_minus[12];                        			      	      // 50.0-70.0  
     ZptSystBkg_minus12[8] =TMath::Max(ZptSystBkg_minus[13],ZptSystBkg_minus[14]);            	      // 70.0-110.0 
    cout <<"Biggest one of ZptSystBkg_minus 13,14 \t"<<ZptSystBkg_minus12[8]<<endl;
     ZptSystBkg_minus12[9] =ZptSystBkg_minus[15];                            	        	      // 110.0-150.0
     ZptSystBkg_minus12[10]=ZptSystBkg_minus[16];                                   	              // 150.0-190.0
     ZptSystBkg_minus12[11]=ZptSystBkg_minus[17];             		            		      // 190.0-250.0
     ZptSystBkg_minus12[12]=ZptSystBkg_minus[18];                       			    	      // 250.0-600.0
 
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
	  ZptSystEff_plus[i]=(ZptSystEff_plus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}

  double ZptSystEff_plus12[13]={0};
    cout <<"ZptSystEff_plus 1,2,3 \t"<<ZptSystEff_plus[1]<<"\t"<<ZptSystEff_plus[2]<<"\t"<<ZptSystEff_plus[3]<<endl;
    if (ZptSystEff_plus[1]>=ZptSystEff_plus[2] && ZptSystEff_plus[1]>=ZptSystEff_plus[3])
       ZptSystEff_plus12[1]=ZptSystEff_plus[1];
     if (ZptSystEff_plus[2]>=ZptSystEff_plus[1] && ZptSystEff_plus[2]>=ZptSystEff_plus[3])
       ZptSystEff_plus12[1]=ZptSystEff_plus[2];
     if (ZptSystEff_plus[3]>=ZptSystEff_plus[1] && ZptSystEff_plus[3]>=ZptSystEff_plus[2])
       ZptSystEff_plus12[1]=ZptSystEff_plus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptSystEff_plus 1,2,3 \t"<<ZptSystEff_plus12[1]<<endl;
     ZptSystEff_plus12[2] =TMath::Max(ZptSystEff_plus[4],ZptSystEff_plus[5]);               	      // 7.5-12.5   
     ZptSystEff_plus12[3] =TMath::Max(ZptSystEff_plus[6],ZptSystEff_plus[7]);                        // 12.5-17.5  
     ZptSystEff_plus12[4] =TMath::Max(ZptSystEff_plus[8],ZptSystEff_plus[9]);                        // 17.5-30.0  
     ZptSystEff_plus12[5] =ZptSystEff_plus[10];                           		              // 30.0-40.0  
     ZptSystEff_plus12[6] =ZptSystEff_plus[11];                         			              // 40.0-50.0  
     ZptSystEff_plus12[7] =ZptSystEff_plus[12];                        			      	      // 50.0-70.0  
     ZptSystEff_plus12[8] =TMath::Max(ZptSystEff_plus[13],ZptSystEff_plus[14]);            	      // 70.0-110.0 
    cout <<"Biggest one of ZptSystEff_plus 13,14 \t"<<ZptSystEff_plus12[8]<<endl;
     ZptSystEff_plus12[9] =ZptSystEff_plus[15];                            	        	      // 110.0-150.0
     ZptSystEff_plus12[10]=ZptSystEff_plus[16];                                   	              // 150.0-190.0
     ZptSystEff_plus12[11]=ZptSystEff_plus[17];             		            		      // 190.0-250.0
     ZptSystEff_plus12[12]=ZptSystEff_plus[18];                       			    	      // 250.0-600.0
	
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
	  ZptSystEff_minus[i]=(ZptSystEff_minus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
	
  double ZptSystEff_minus12[13]={0};
    cout <<"ZptSystEff_minus 1,2,3 \t"<<ZptSystEff_minus[1]<<"\t"<<ZptSystEff_minus[2]<<"\t"<<ZptSystEff_minus[3]<<endl;
    if (ZptSystEff_minus[1]>=ZptSystEff_minus[2] && ZptSystEff_minus[1]>=ZptSystEff_minus[3])
       ZptSystEff_minus12[1]=ZptSystEff_minus[1];
     if (ZptSystEff_minus[2]>=ZptSystEff_minus[1] && ZptSystEff_minus[2]>=ZptSystEff_minus[3])
       ZptSystEff_minus12[1]=ZptSystEff_minus[2];
     if (ZptSystEff_minus[3]>=ZptSystEff_minus[1] && ZptSystEff_minus[3]>=ZptSystEff_minus[2])
       ZptSystEff_minus12[1]=ZptSystEff_minus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptSystEff_minus 1,2,3 \t"<<ZptSystEff_minus12[1]<<endl;
     ZptSystEff_minus12[2] =TMath::Max(ZptSystEff_minus[4],ZptSystEff_minus[5]);               	      // 7.5-12.5   
     ZptSystEff_minus12[3] =TMath::Max(ZptSystEff_minus[6],ZptSystEff_minus[7]);                        // 12.5-17.5  
     ZptSystEff_minus12[4] =TMath::Max(ZptSystEff_minus[8],ZptSystEff_minus[9]);                        // 17.5-30.0  
     ZptSystEff_minus12[5] =ZptSystEff_minus[10];                           		              // 30.0-40.0  
     ZptSystEff_minus12[6] =ZptSystEff_minus[11];                         			              // 40.0-50.0  
     ZptSystEff_minus12[7] =ZptSystEff_minus[12];                        			      	      // 50.0-70.0  
     ZptSystEff_minus12[8] =TMath::Max(ZptSystEff_minus[13],ZptSystEff_minus[14]);            	      // 70.0-110.0 
    cout <<"Biggest one of ZptSystEff_minus 13,14\t"<<ZptSystEff_minus12[8]<<endl;
     ZptSystEff_minus12[9] =ZptSystEff_minus[15];                            	        	      // 110.0-150.0
     ZptSystEff_minus12[10]=ZptSystEff_minus[16];                                   	              // 150.0-190.0
     ZptSystEff_minus12[11]=ZptSystEff_minus[17];             		            		      // 190.0-250.0
     ZptSystEff_minus12[12]=ZptSystEff_minus[18];                       			    	      // 250.0-600.0

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
	  ZptSystFSRErr[i]=(ZptSystFSRErr[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}

  double ZptSystFSRErr12[13]={0};
    cout <<"ZptSystFSRErr 1,2,3 \t"<<ZptSystFSRErr[1]<<"\t"<<ZptSystFSRErr[2]<<"\t"<<ZptSystFSRErr[3]<<endl;
    if (ZptSystFSRErr[1]>=ZptSystFSRErr[2] && ZptSystFSRErr[1]>=ZptSystFSRErr[3])
       ZptSystFSRErr12[1]=ZptSystFSRErr[1];
     if (ZptSystFSRErr[2]>=ZptSystFSRErr[1] && ZptSystFSRErr[2]>=ZptSystFSRErr[3])
       ZptSystFSRErr12[1]=ZptSystFSRErr[2];
     if (ZptSystFSRErr[3]>=ZptSystFSRErr[1] && ZptSystFSRErr[3]>=ZptSystFSRErr[2])
       ZptSystFSRErr12[1]=ZptSystFSRErr[3];						      // 0.0-7.5
    cout <<"Biggest one of ZptSystFSRErr 1,2,3 \t"<<ZptSystFSRErr12[1]<<endl;
     ZptSystFSRErr12[2] =TMath::Max(ZptSystFSRErr[4],ZptSystFSRErr[5]);               	      // 7.5-12.5   
     ZptSystFSRErr12[3] =TMath::Max(ZptSystFSRErr[6],ZptSystFSRErr[7]);                       // 12.5-17.5  
     ZptSystFSRErr12[4] =TMath::Max(ZptSystFSRErr[8],ZptSystFSRErr[9]);                       // 17.5-30.0  
     ZptSystFSRErr12[5] =ZptSystFSRErr[10];                           		              // 30.0-40.0  
     ZptSystFSRErr12[6] =ZptSystFSRErr[11];                         			      // 40.0-50.0  
     ZptSystFSRErr12[7] =ZptSystFSRErr[12];                        			      // 50.0-70.0  
     ZptSystFSRErr12[8] =TMath::Max(ZptSystFSRErr[13],ZptSystFSRErr[14]);            	      // 70.0-110.0 
    cout <<"Biggest one of ZptSystFSRErr 13,14 \t"<<ZptSystFSRErr12[8]<<endl;
     ZptSystFSRErr12[9] =ZptSystFSRErr[15];                            	        	      // 110.0-150.0
     ZptSystFSRErr12[10]=ZptSystFSRErr[16];                                   	              // 150.0-190.0
     ZptSystFSRErr12[11]=ZptSystFSRErr[17];             		            	      // 190.0-250.0
     ZptSystFSRErr12[12]=ZptSystFSRErr[18];                       			      // 250.0-600.0

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
	  ZptUnfoldSyst_doubleGaussin_plus[i]=(ZptUnfoldSyst_doubleGaussin_plus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}


  double ZptUnfoldSyst_doubleGaussin_plus12[13]={0};
    cout <<"ZptUnfoldSyst_doubleGaussin_plus 1,2,3 \t"<<ZptUnfoldSyst_doubleGaussin_plus[1]<<"\t"<<ZptUnfoldSyst_doubleGaussin_plus[2]<<"\t"<<ZptUnfoldSyst_doubleGaussin_plus[3]<<endl;
    if (ZptUnfoldSyst_doubleGaussin_plus[1]>=ZptUnfoldSyst_doubleGaussin_plus[2] && ZptUnfoldSyst_doubleGaussin_plus[1]>=ZptUnfoldSyst_doubleGaussin_plus[3])
       ZptUnfoldSyst_doubleGaussin_plus12[1]=ZptUnfoldSyst_doubleGaussin_plus[1];
     if (ZptUnfoldSyst_doubleGaussin_plus[2]>=ZptUnfoldSyst_doubleGaussin_plus[1] && ZptUnfoldSyst_doubleGaussin_plus[2]>=ZptUnfoldSyst_doubleGaussin_plus[3])
       ZptUnfoldSyst_doubleGaussin_plus12[1]=ZptUnfoldSyst_doubleGaussin_plus[2];
     if (ZptUnfoldSyst_doubleGaussin_plus[3]>=ZptUnfoldSyst_doubleGaussin_plus[1] && ZptUnfoldSyst_doubleGaussin_plus[3]>=ZptUnfoldSyst_doubleGaussin_plus[2])
       ZptUnfoldSyst_doubleGaussin_plus12[1]=ZptUnfoldSyst_doubleGaussin_plus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptUnfoldSyst_doubleGaussin_plus 1,2,3 \t"<<ZptUnfoldSyst_doubleGaussin_plus12[1]<<endl;
     ZptUnfoldSyst_doubleGaussin_plus12[2] =TMath::Max(ZptUnfoldSyst_doubleGaussin_plus[4],ZptUnfoldSyst_doubleGaussin_plus[5]);          	      // 7.5-12.5   
     ZptUnfoldSyst_doubleGaussin_plus12[3] =TMath::Max(ZptUnfoldSyst_doubleGaussin_plus[6],ZptUnfoldSyst_doubleGaussin_plus[7]);                 // 12.5-17.5  
     ZptUnfoldSyst_doubleGaussin_plus12[4] =TMath::Max(ZptUnfoldSyst_doubleGaussin_plus[8],ZptUnfoldSyst_doubleGaussin_plus[9]);                 // 17.5-30.0  
     ZptUnfoldSyst_doubleGaussin_plus12[5] =ZptUnfoldSyst_doubleGaussin_plus[10];                           		              // 30.0-40.0  
     ZptUnfoldSyst_doubleGaussin_plus12[6] =ZptUnfoldSyst_doubleGaussin_plus[11];                         		              // 40.0-50.0  
     ZptUnfoldSyst_doubleGaussin_plus12[7] =ZptUnfoldSyst_doubleGaussin_plus[12];                        			      	      // 50.0-70.0  
     ZptUnfoldSyst_doubleGaussin_plus12[8] =TMath::Max(ZptUnfoldSyst_doubleGaussin_plus[13],ZptUnfoldSyst_doubleGaussin_plus[14]);        	      // 70.0-110.0 
    cout <<"Biggest one of ZptUnfoldSyst_doubleGaussin_plus 13,14 \t"<<ZptUnfoldSyst_doubleGaussin_plus12[8]<<endl;
     ZptUnfoldSyst_doubleGaussin_plus12[9] =ZptUnfoldSyst_doubleGaussin_plus[15];                            	        	      // 110.0-150.0
     ZptUnfoldSyst_doubleGaussin_plus12[10]=ZptUnfoldSyst_doubleGaussin_plus[16];                                   	              // 150.0-190.0
     ZptUnfoldSyst_doubleGaussin_plus12[11]=ZptUnfoldSyst_doubleGaussin_plus[17];             		            		      // 190.0-250.0
     ZptUnfoldSyst_doubleGaussin_plus12[12]=ZptUnfoldSyst_doubleGaussin_plus[18];                       			    	      // 250.0-600.0


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
	  ZptUnfoldSyst_doubleGaussin_minus[i]=(ZptUnfoldSyst_doubleGaussin_minus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
	
  double ZptUnfoldSyst_doubleGaussin_minus12[13]={0};
    cout <<"ZptUnfoldSyst_doubleGaussin_minus 1,2,3 \t"<<ZptUnfoldSyst_doubleGaussin_minus[1]<<"\t"<<ZptUnfoldSyst_doubleGaussin_minus[2]<<"\t"<<ZptUnfoldSyst_doubleGaussin_minus[3]<<endl;
    if (ZptUnfoldSyst_doubleGaussin_minus[1]>=ZptUnfoldSyst_doubleGaussin_minus[2] && ZptUnfoldSyst_doubleGaussin_minus[1]>=ZptUnfoldSyst_doubleGaussin_minus[3])
       ZptUnfoldSyst_doubleGaussin_minus12[1]=ZptUnfoldSyst_doubleGaussin_minus[1];
     if (ZptUnfoldSyst_doubleGaussin_minus[2]>=ZptUnfoldSyst_doubleGaussin_minus[1] && ZptUnfoldSyst_doubleGaussin_minus[2]>=ZptUnfoldSyst_doubleGaussin_minus[3])
       ZptUnfoldSyst_doubleGaussin_minus12[1]=ZptUnfoldSyst_doubleGaussin_minus[2];
     if (ZptUnfoldSyst_doubleGaussin_minus[3]>=ZptUnfoldSyst_doubleGaussin_minus[1] && ZptUnfoldSyst_doubleGaussin_minus[3]>=ZptUnfoldSyst_doubleGaussin_minus[2])
       ZptUnfoldSyst_doubleGaussin_minus12[1]=ZptUnfoldSyst_doubleGaussin_minus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptUnfoldSyst_doubleGaussin_minus 1,2,3 \t"<<ZptUnfoldSyst_doubleGaussin_minus12[1]<<endl;
     ZptUnfoldSyst_doubleGaussin_minus12[2] =TMath::Max(ZptUnfoldSyst_doubleGaussin_minus[4],ZptUnfoldSyst_doubleGaussin_minus[5]);         	      // 7.5-12.5   
     ZptUnfoldSyst_doubleGaussin_minus12[3] =TMath::Max(ZptUnfoldSyst_doubleGaussin_minus[6],ZptUnfoldSyst_doubleGaussin_minus[7]);                 // 12.5-17.5  
     ZptUnfoldSyst_doubleGaussin_minus12[4] =TMath::Max(ZptUnfoldSyst_doubleGaussin_minus[8],ZptUnfoldSyst_doubleGaussin_minus[9]);                 // 17.5-30.0  
     ZptUnfoldSyst_doubleGaussin_minus12[5] =ZptUnfoldSyst_doubleGaussin_minus[10];                           		              // 30.0-40.0  
     ZptUnfoldSyst_doubleGaussin_minus12[6] =ZptUnfoldSyst_doubleGaussin_minus[11];                         		              // 40.0-50.0  
     ZptUnfoldSyst_doubleGaussin_minus12[7] =ZptUnfoldSyst_doubleGaussin_minus[12];                        			      	      // 50.0-70.0  
     ZptUnfoldSyst_doubleGaussin_minus12[8] =TMath::Max(ZptUnfoldSyst_doubleGaussin_minus[13],ZptUnfoldSyst_doubleGaussin_minus[14]);        	      // 70.0-110.0 
    cout <<"Biggest one of ZptUnfoldSyst_doubleGaussin_minus 13,14 \t"<<ZptUnfoldSyst_doubleGaussin_minus[8]<<endl;
     ZptUnfoldSyst_doubleGaussin_minus12[9] =ZptUnfoldSyst_doubleGaussin_minus[15];                            	        	      // 110.0-150.0
     ZptUnfoldSyst_doubleGaussin_minus12[10]=ZptUnfoldSyst_doubleGaussin_minus[16];                                   	              // 150.0-190.0
     ZptUnfoldSyst_doubleGaussin_minus12[11]=ZptUnfoldSyst_doubleGaussin_minus[17];             		            		      // 190.0-250.0
     ZptUnfoldSyst_doubleGaussin_minus12[12]=ZptUnfoldSyst_doubleGaussin_minus[18];                       			    	      // 250.0-600.0

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
	  ZptUnfoldingShapeSyst_plus[i]=(ZptUnfoldingShapeSyst_plus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
	
  double ZptUnfoldingShapeSyst_plus12[13]={0};
    cout <<"ZptUnfoldingShapeSyst_plus 1,2,3 \t"<<ZptUnfoldingShapeSyst_plus[1]<<"\t"<<ZptUnfoldingShapeSyst_plus[2]<<"\t"<<ZptUnfoldingShapeSyst_plus[3]<<endl;
    if (ZptUnfoldingShapeSyst_plus[1]>=ZptUnfoldingShapeSyst_plus[2] && ZptUnfoldingShapeSyst_plus[1]>=ZptUnfoldingShapeSyst_plus[3])
       ZptUnfoldingShapeSyst_plus12[1]=ZptUnfoldingShapeSyst_plus[1];
     if (ZptUnfoldingShapeSyst_plus[2]>=ZptUnfoldingShapeSyst_plus[1] && ZptUnfoldingShapeSyst_plus[2]>=ZptUnfoldingShapeSyst_plus[3])
       ZptUnfoldingShapeSyst_plus12[1]=ZptUnfoldingShapeSyst_plus[2];
     if (ZptUnfoldingShapeSyst_plus[3]>=ZptUnfoldingShapeSyst_plus[1] && ZptUnfoldingShapeSyst_plus[3]>=ZptUnfoldingShapeSyst_plus[2])
       ZptUnfoldingShapeSyst_plus12[1]=ZptUnfoldingShapeSyst_plus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptUnfoldingShapeSyst_plus 1,2,3 \t"<<ZptUnfoldingShapeSyst_plus12[1]<<endl;
     ZptUnfoldingShapeSyst_plus12[2] =TMath::Max(ZptUnfoldingShapeSyst_plus[4],ZptUnfoldingShapeSyst_plus[5]);                 // 7.5-12.5   
     ZptUnfoldingShapeSyst_plus12[3] =TMath::Max(ZptUnfoldingShapeSyst_plus[6],ZptUnfoldingShapeSyst_plus[7]);                 // 12.5-17.5  
     ZptUnfoldingShapeSyst_plus12[4] =TMath::Max(ZptUnfoldingShapeSyst_plus[8],ZptUnfoldingShapeSyst_plus[9]);                 // 17.5-30.0  
     ZptUnfoldingShapeSyst_plus12[5] =ZptUnfoldingShapeSyst_plus[10];                           		              // 30.0-40.0  
     ZptUnfoldingShapeSyst_plus12[6] =ZptUnfoldingShapeSyst_plus[11];                         			      // 40.0-50.0  
     ZptUnfoldingShapeSyst_plus12[7] =ZptUnfoldingShapeSyst_plus[12];                        			      	      // 50.0-70.0  
     ZptUnfoldingShapeSyst_plus12[8] =TMath::Max(ZptUnfoldingShapeSyst_plus[13],ZptUnfoldingShapeSyst_plus[14]);               // 70.0-110.0 
    cout <<"Biggest one of ZptUnfoldingShapeSyst_plus 13,14 \t"<<ZptUnfoldingShapeSyst_plus[8]<<endl;
     ZptUnfoldingShapeSyst_plus12[9] =ZptUnfoldingShapeSyst_plus[15];                            	        	      // 110.0-150.0
     ZptUnfoldingShapeSyst_plus12[10]=ZptUnfoldingShapeSyst_plus[16];                                   	              // 150.0-190.0
     ZptUnfoldingShapeSyst_plus12[11]=ZptUnfoldingShapeSyst_plus[17];             		            		      // 190.0-250.0
     ZptUnfoldingShapeSyst_plus12[12]=ZptUnfoldingShapeSyst_plus[18];                       			    	      // 250.0-600.0


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
	  ZptUnfoldingShapeSyst_minus[i]=(ZptUnfoldingShapeSyst_minus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
	
  double ZptUnfoldingShapeSyst_minus12[13]={0};
    cout <<"ZptUnfoldingShapeSyst_minus 1,2,3 \t"<<ZptUnfoldingShapeSyst_minus[1]<<"\t"<<ZptUnfoldingShapeSyst_minus[2]<<"\t"<<ZptUnfoldingShapeSyst_minus[3]<<endl;
    if (ZptUnfoldingShapeSyst_minus[1]>=ZptUnfoldingShapeSyst_minus[2] && ZptUnfoldingShapeSyst_minus[1]>=ZptUnfoldingShapeSyst_minus[3])
       ZptUnfoldingShapeSyst_minus12[1]=ZptUnfoldingShapeSyst_minus[1];
     if (ZptUnfoldingShapeSyst_minus[2]>=ZptUnfoldingShapeSyst_minus[1] && ZptUnfoldingShapeSyst_minus[2]>=ZptUnfoldingShapeSyst_minus[3])
       ZptUnfoldingShapeSyst_minus12[1]=ZptUnfoldingShapeSyst_minus[2];
     if (ZptUnfoldingShapeSyst_minus[3]>=ZptUnfoldingShapeSyst_minus[1] && ZptUnfoldingShapeSyst_minus[3]>=ZptUnfoldingShapeSyst_minus[2])
       ZptUnfoldingShapeSyst_minus12[1]=ZptUnfoldingShapeSyst_minus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptUnfoldingShapeSyst_minus 1,2,3 \t"<<ZptUnfoldingShapeSyst_minus12[1]<<endl;
     ZptUnfoldingShapeSyst_minus12[2] =TMath::Max(ZptUnfoldingShapeSyst_minus[4],ZptUnfoldingShapeSyst_minus[5]);                 // 7.5-12.5   
     ZptUnfoldingShapeSyst_minus12[3] =TMath::Max(ZptUnfoldingShapeSyst_minus[6],ZptUnfoldingShapeSyst_minus[7]);                 // 12.5-17.5  
     ZptUnfoldingShapeSyst_minus12[4] =TMath::Max(ZptUnfoldingShapeSyst_minus[8],ZptUnfoldingShapeSyst_minus[9]);                 // 17.5-30.0  
     ZptUnfoldingShapeSyst_minus12[5] =ZptUnfoldingShapeSyst_minus[10];                           		              // 30.0-40.0  
     ZptUnfoldingShapeSyst_minus12[6] =ZptUnfoldingShapeSyst_minus[11];                         			      // 40.0-50.0  
     ZptUnfoldingShapeSyst_minus12[7] =ZptUnfoldingShapeSyst_minus[12];                        			      	      // 50.0-70.0  
     ZptUnfoldingShapeSyst_minus12[8] =TMath::Max(ZptUnfoldingShapeSyst_minus[13],ZptUnfoldingShapeSyst_minus[14]);               // 70.0-110.0 
    cout <<"Biggest one of ZptUnfoldingShapeSyst_minus 13,14 \t"<<ZptUnfoldingShapeSyst_minus12[8]<<endl;
     ZptUnfoldingShapeSyst_minus12[9] =ZptUnfoldingShapeSyst_minus[15];                            	        	      // 110.0-150.0
     ZptUnfoldingShapeSyst_minus12[10]=ZptUnfoldingShapeSyst_minus[16];                                   	              // 150.0-190.0
     ZptUnfoldingShapeSyst_minus12[11]=ZptUnfoldingShapeSyst_minus[17];             		            		      // 190.0-250.0
     ZptUnfoldingShapeSyst_minus12[12]=ZptUnfoldingShapeSyst_minus[18];                       			    	      // 250.0-600.0
       
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
	  cout<<
             ZptMadgraphSyst_plus[i]<<"\t"<<
             hZptYield->GetBinContent(i)<<"\t"<<
             hZptYield->GetBinContent(i)/(hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i))<<"\t"<<
             100.*ZptMadgraphSyst_plus[i]/(hZptYield->GetBinContent(i)/(hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i)))<<"\t"<<
	    endl;
	}
	
    cout << fixed << setprecision(4);
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptMadgraphSyst_plus[i]=(ZptMadgraphSyst_plus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	  
	  cout<<
             ZptMadgraphSyst_plus[i]<<"\t %"
	    <<endl;
	}
    cout << fixed << setprecision(8);
  
    
    double ZptMadgraphSyst_plus12[13]={0};
    cout <<"ZptMadgraphSyst_plus 1,2,3 \t"<<ZptMadgraphSyst_plus[1]<<"\t"<<ZptMadgraphSyst_plus[2]<<"\t"<<ZptMadgraphSyst_plus[3]<<endl;
    //if (ZptMadgraphSyst_plus[1]>=ZptMadgraphSyst_plus[2] && ZptMadgraphSyst_plus[1]>=ZptMadgraphSyst_plus[3])
    //   ZptMadgraphSyst_plus12[1]=ZptMadgraphSyst_plus[1];
    // if (ZptMadgraphSyst_plus[2]>=ZptMadgraphSyst_plus[1] && ZptMadgraphSyst_plus[2]>=ZptMadgraphSyst_plus[3])
    //   ZptMadgraphSyst_plus12[1]=ZptMadgraphSyst_plus[2];
    // if (ZptMadgraphSyst_plus[3]>=ZptMadgraphSyst_plus[1] && ZptMadgraphSyst_plus[3]>=ZptMadgraphSyst_plus[2])
    //   ZptMadgraphSyst_plus12[1]=ZptMadgraphSyst_plus[3];							      // 0.0-7.5
    //cout <<"Biggest one of ZptMadgraphSyst_plus 1,2,3 \t"<<ZptMadgraphSyst_plus12[1]<<endl;
    // ZptMadgraphSyst_plus12[2] =TMath::Max(ZptMadgraphSyst_plus[4],ZptMadgraphSyst_plus[5]);                 // 7.5-12.5   
    // ZptMadgraphSyst_plus12[3] =TMath::Max(ZptMadgraphSyst_plus[6],ZptMadgraphSyst_plus[7]);                 // 12.5-17.5  
    // ZptMadgraphSyst_plus12[4] =TMath::Max(ZptMadgraphSyst_plus[8],ZptMadgraphSyst_plus[9]);                 // 17.5-30.0  
    // ZptMadgraphSyst_plus12[5] =ZptMadgraphSyst_plus[10];                           		              // 30.0-40.0  
    // ZptMadgraphSyst_plus12[6] =ZptMadgraphSyst_plus[11];                         			      // 40.0-50.0  
    // ZptMadgraphSyst_plus12[7] =ZptMadgraphSyst_plus[12];                        			      	      // 50.0-70.0  
    // ZptMadgraphSyst_plus12[8] =TMath::Max(ZptMadgraphSyst_plus[13],ZptMadgraphSyst_plus[14]);               // 70.0-110.0 
    //cout <<"Biggest one of ZptMadgraphSyst_plus 13,14 \t"<<ZptMadgraphSyst_plus12[8]<<endl;
    // ZptMadgraphSyst_plus12[9] =ZptMadgraphSyst_plus[15];                            	        	      // 110.0-150.0
    // ZptMadgraphSyst_plus12[10]=ZptMadgraphSyst_plus[16];                                   	              // 150.0-190.0
    // ZptMadgraphSyst_plus12[11]=ZptMadgraphSyst_plus[17];             		            		      // 190.0-250.0
    // ZptMadgraphSyst_plus12[12]=ZptMadgraphSyst_plus[18];                       			    	      // 250.0-600.0

    /////   12 bin syst from merged Powheg and Madgraph.
     ZptMadgraphSyst_plus12[1] =0.2255   ;	      // 0.0-7.5
     ZptMadgraphSyst_plus12[2] =0.0067   ;          // 7.5-12.5   
     ZptMadgraphSyst_plus12[3] =0.5737   ;          // 12.5-17.5  
     ZptMadgraphSyst_plus12[4] =0.1332   ;          // 17.5-30.0  
     ZptMadgraphSyst_plus12[5] =0.5616   ;           // 30.0-40.0  
     ZptMadgraphSyst_plus12[6] =1.0259   ;           // 40.0-50.0  
     ZptMadgraphSyst_plus12[7] =0.2607   ;   	      // 50.0-70.0  
     ZptMadgraphSyst_plus12[8] =0.4619   ;          // 70.0-110.0 
     ZptMadgraphSyst_plus12[9] =1.1441   ;           // 110.0-150.0
     ZptMadgraphSyst_plus12[10]=0.1438   ;           // 150.0-190.0
     ZptMadgraphSyst_plus12[11]=9.9109   ;           // 190.0-250.0
     ZptMadgraphSyst_plus12[12]=0.4505   ; 	      // 250.0-600.0


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
	///// Zpt Madgraph syst error converting from NormalizedError to Error and to %
	for(int i(1);i<nBinsZpt;i++)
	{
	  ZptMadgraphSyst_minus[i]=(ZptMadgraphSyst_minus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
	
  double ZptMadgraphSyst_minus12[13]={0};
    //cout <<"ZptMadgraphSyst_minus 1,2,3 \t"<<ZptMadgraphSyst_minus[1]<<"\t"<<ZptMadgraphSyst_minus[2]<<"\t"<<ZptMadgraphSyst_minus[3]<<endl;
    //if (ZptMadgraphSyst_minus[1]>=ZptMadgraphSyst_minus[2] && ZptMadgraphSyst_minus[1]>=ZptMadgraphSyst_minus[3])
    //   ZptMadgraphSyst_minus12[1]=ZptMadgraphSyst_minus[1];
    // if (ZptMadgraphSyst_minus[2]>=ZptMadgraphSyst_minus[1] && ZptMadgraphSyst_minus[2]>=ZptMadgraphSyst_minus[3])
    //   ZptMadgraphSyst_minus12[1]=ZptMadgraphSyst_minus[2];
    // if (ZptMadgraphSyst_minus[3]>=ZptMadgraphSyst_minus[1] && ZptMadgraphSyst_minus[3]>=ZptMadgraphSyst_minus[2])
    //   ZptMadgraphSyst_minus12[1]=ZptMadgraphSyst_minus[3];							      // 0.0-7.5
    //cout <<"Biggest one of ZptMadgraphSyst_minus 1,2,3 \t"<<ZptMadgraphSyst_minus12[1]<<endl;
    // ZptMadgraphSyst_minus12[2] =TMath::Max(ZptMadgraphSyst_minus[4],ZptMadgraphSyst_minus[5]);                 // 7.5-12.5   
    // ZptMadgraphSyst_minus12[3] =TMath::Max(ZptMadgraphSyst_minus[6],ZptMadgraphSyst_minus[7]);                 // 12.5-17.5  
    // ZptMadgraphSyst_minus12[4] =TMath::Max(ZptMadgraphSyst_minus[8],ZptMadgraphSyst_minus[9]);                 // 17.5-30.0  
    // ZptMadgraphSyst_minus12[5] =ZptMadgraphSyst_minus[10];                           		              // 30.0-40.0  
    // ZptMadgraphSyst_minus12[6] =ZptMadgraphSyst_minus[11];                         			      // 40.0-50.0  
    // ZptMadgraphSyst_minus12[7] =ZptMadgraphSyst_minus[12];                        			      	      // 50.0-70.0  
    // ZptMadgraphSyst_minus12[8] =TMath::Max(ZptMadgraphSyst_minus[13],ZptMadgraphSyst_minus[14]);               // 70.0-110.0 
    //cout <<"Biggest one of ZptMadgraphSyst_minus 13,14 \t"<<ZptMadgraphSyst_minus12[8]<<endl;
    // ZptMadgraphSyst_minus12[9] =ZptMadgraphSyst_minus[15];                            	        	      // 110.0-150.0
    // ZptMadgraphSyst_minus12[10]=ZptMadgraphSyst_minus[16];                                   	              // 150.0-190.0
    // ZptMadgraphSyst_minus12[11]=ZptMadgraphSyst_minus[17];             		            		      // 190.0-250.0
    // ZptMadgraphSyst_minus12[12]=ZptMadgraphSyst_minus[18];                       			    	      // 250.0-600.0
     
     ZptMadgraphSyst_minus12[1] =0.2255   ;	      // 0.0-7.5
     ZptMadgraphSyst_minus12[2] =0.0067   ;          // 7.5-12.5   
     ZptMadgraphSyst_minus12[3] =0.5737   ;          // 12.5-17.5  
     ZptMadgraphSyst_minus12[4] =0.1332   ;          // 17.5-30.0  
     ZptMadgraphSyst_minus12[5] =0.5616   ;           // 30.0-40.0  
     ZptMadgraphSyst_minus12[6] =1.0259   ;           // 40.0-50.0  
     ZptMadgraphSyst_minus12[7] =0.2607   ;   	      // 50.0-70.0  
     ZptMadgraphSyst_minus12[8] =0.4619   ;          // 70.0-110.0 
     ZptMadgraphSyst_minus12[9] =1.1441   ;           // 110.0-150.0
     ZptMadgraphSyst_minus12[10]=0.1438   ;           // 150.0-190.0
     ZptMadgraphSyst_minus12[11]=9.9109   ;           // 190.0-250.0
     ZptMadgraphSyst_minus12[12]=0.4505   ; 	      // 250.0-600.0

       

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
	  ZptAcceptSyst_plus[i]=(ZptAcceptSyst_plus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
	
  double ZptAcceptSyst_plus12[13]={0};
    cout <<"ZptAcceptSyst_plus 1,2,3 \t"<<ZptAcceptSyst_plus[1]<<"\t"<<ZptAcceptSyst_plus[2]<<"\t"<<ZptAcceptSyst_plus[3]<<endl;
    if (ZptAcceptSyst_plus[1]>=ZptAcceptSyst_plus[2] && ZptAcceptSyst_plus[1]>=ZptAcceptSyst_plus[3])
       ZptAcceptSyst_plus12[1]=ZptAcceptSyst_plus[1];
     if (ZptAcceptSyst_plus[2]>=ZptAcceptSyst_plus[1] && ZptAcceptSyst_plus[2]>=ZptAcceptSyst_plus[3])
       ZptAcceptSyst_plus12[1]=ZptAcceptSyst_plus[2];
     if (ZptAcceptSyst_plus[3]>=ZptAcceptSyst_plus[1] && ZptAcceptSyst_plus[3]>=ZptAcceptSyst_plus[2])
       ZptAcceptSyst_plus12[1]=ZptAcceptSyst_plus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptAcceptSyst_plus 1,2,3 \t"<<ZptAcceptSyst_plus12[1]<<endl;
     ZptAcceptSyst_plus12[2] =TMath::Max(ZptAcceptSyst_plus[4],ZptAcceptSyst_plus[5]);                 // 7.5-12.5   
     ZptAcceptSyst_plus12[3] =TMath::Max(ZptAcceptSyst_plus[6],ZptAcceptSyst_plus[7]);                 // 12.5-17.5  
     ZptAcceptSyst_plus12[4] =TMath::Max(ZptAcceptSyst_plus[8],ZptAcceptSyst_plus[9]);                 // 17.5-30.0  
     ZptAcceptSyst_plus12[5] =ZptAcceptSyst_plus[10];                           		              // 30.0-40.0  
     ZptAcceptSyst_plus12[6] =ZptAcceptSyst_plus[11];                         			      // 40.0-50.0  
     ZptAcceptSyst_plus12[7] =ZptAcceptSyst_plus[12];                        			      	      // 50.0-70.0  
     ZptAcceptSyst_plus12[8] =TMath::Max(ZptAcceptSyst_plus[13],ZptAcceptSyst_plus[14]);               // 70.0-110.0 
    cout <<"Biggest one of ZptAcceptSyst_plus 13,14 \t"<<ZptAcceptSyst_plus12[8]<<endl;
     ZptAcceptSyst_plus12[9] =ZptAcceptSyst_plus[15];                            	        	      // 110.0-150.0
     ZptAcceptSyst_plus12[10]=ZptAcceptSyst_plus[16];                                   	              // 150.0-190.0
     ZptAcceptSyst_plus12[11]=ZptAcceptSyst_plus[17];             		            		      // 190.0-250.0
     ZptAcceptSyst_plus12[12]=ZptAcceptSyst_plus[18];                       			    	      // 250.0-600.0
       
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
	  ZptAcceptSyst_minus[i]=(ZptAcceptSyst_minus[i]*hZptYield->Integral()*hZptBins_LinScale->GetXaxis()->GetBinWidth(i) )*100./hZptYield->GetBinContent(i);
	}
	
  double ZptAcceptSyst_minus12[13]={0};
    cout <<"ZptAcceptSyst_minus 1,2,3 \t"<<ZptAcceptSyst_minus[1]<<"\t"<<ZptAcceptSyst_minus[2]<<"\t"<<ZptAcceptSyst_minus[3]<<endl;
    if (ZptAcceptSyst_minus[1]>=ZptAcceptSyst_minus[2] && ZptAcceptSyst_minus[1]>=ZptAcceptSyst_minus[3])
       ZptAcceptSyst_minus12[1]=ZptAcceptSyst_minus[1];
     if (ZptAcceptSyst_minus[2]>=ZptAcceptSyst_minus[1] && ZptAcceptSyst_minus[2]>=ZptAcceptSyst_minus[3])
       ZptAcceptSyst_minus12[1]=ZptAcceptSyst_minus[2];
     if (ZptAcceptSyst_minus[3]>=ZptAcceptSyst_minus[1] && ZptAcceptSyst_minus[3]>=ZptAcceptSyst_minus[2])
       ZptAcceptSyst_minus12[1]=ZptAcceptSyst_minus[3];							      // 0.0-7.5
    cout <<"Biggest one of ZptAcceptSyst_minus 1,2,3 \t"<<ZptAcceptSyst_minus12[1]<<endl;
     ZptAcceptSyst_minus12[2] =TMath::Max(ZptAcceptSyst_minus[4],ZptAcceptSyst_minus[5]);                 // 7.5-12.5   
     ZptAcceptSyst_minus12[3] =TMath::Max(ZptAcceptSyst_minus[6],ZptAcceptSyst_minus[7]);                 // 12.5-17.5  
     ZptAcceptSyst_minus12[4] =TMath::Max(ZptAcceptSyst_minus[8],ZptAcceptSyst_minus[9]);                 // 17.5-30.0  
     ZptAcceptSyst_minus12[5] =ZptAcceptSyst_minus[10];                           		              // 30.0-40.0  
     ZptAcceptSyst_minus12[6] =ZptAcceptSyst_minus[11];                         			      // 40.0-50.0  
     ZptAcceptSyst_minus12[7] =ZptAcceptSyst_minus[12];                        			      	      // 50.0-70.0  
     ZptAcceptSyst_minus12[8] =TMath::Max(ZptAcceptSyst_minus[13],ZptAcceptSyst_minus[14]);               // 70.0-110.0 
    cout <<"Biggest one of ZptAcceptSyst_minus 13,14 \t"<<ZptAcceptSyst_minus12[8]<<endl;
     ZptAcceptSyst_minus12[9] =ZptAcceptSyst_minus[15];                            	        	      // 110.0-150.0
     ZptAcceptSyst_minus12[10]=ZptAcceptSyst_minus[16];                                   	              // 150.0-190.0
     ZptAcceptSyst_minus12[11]=ZptAcceptSyst_minus[17];             		            		      // 190.0-250.0
     ZptAcceptSyst_minus12[12]=ZptAcceptSyst_minus[18];                       			    	      // 250.0-600.0
      


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
    double errorplus[19]={0};
    double errorminus[19]={0};
    cout<<"=================================================================="<<endl;
    cout<<"Zpt Total Uncertainty 18 Bin as in Note in %:"<<endl;
    cout<<"Bin:\t"<<"errorplus\t"<<"errorminus"<<endl;
    for(int i(1);i<nBinsZpt;i++)
    {
      errorplus[i]= TMath::Sqrt( TMath::Power(ZptStatErr_Percent[i],2) + TMath::Power(systplus[i] ,2) );
      errorminus[i]= TMath::Sqrt( TMath::Power(ZptStatErr_Percent[i] ,2) + TMath::Power(systminus[i],2) );
      cout<<i<<"\t+"<<errorplus[i]<<"\t-"<<errorminus[i]<<endl;
    }
    
    
    ///// Calculating Zpt total syst for 12 bins 
    double syst2plus12[13]={0};  
    double syst2minus12[13]={0};
    for(int i(1);i<nBins12;i++)
    {
      syst2plus12[i] = 
      TMath::Power(ZptSystBkg_plus12[i],2)+
      TMath::Power(ZptUnfoldingShapeSyst_plus12[i],2)+
      TMath::Power(ZptUnfoldSyst_doubleGaussin_plus12[i],2)+
      TMath::Power(ZptSystEff_plus12[i],2)+
      TMath::Power(ZptMadgraphSyst_plus12[i],2) +
      TMath::Power(ZptSystFSRErr12[i],2);
      //TMath::Power(ZptAcceptSyst_plus12[i],2); //// Do not use this for Fid Volume Xsec 
    }
    for(int i(1);i<nBins12;i++)
    {
      syst2minus12[i] = 
      TMath::Power(ZptSystBkg_minus12[i],2)+
      TMath::Power(ZptUnfoldingShapeSyst_minus12[i],2)+
      TMath::Power(ZptUnfoldSyst_doubleGaussin_minus12[i],2)+
      TMath::Power(ZptSystEff_minus12[i],2)+
      TMath::Power(ZptMadgraphSyst_minus12[i],2) +
      TMath::Power(ZptSystFSRErr12[i],2);
      //TMath::Power(ZptAcceptSyst_minus12[i],2);  //// Do not use this for Fid Volume Xsec
    }
    
    cout<<"======================Printing 12 bin Errors and X-sections============================================"<<endl;
    cout<<" All systematics in %:"<<endl;

    cout<<"Bin:\t"<<
      "+BkgSyst\t"<<"-BkgSyst\t"<<
      "+UnfShape\t"<<"-UnfShape\t"<<
      "+UnfDoubleGaus\t"<<"-UnfDoubleGaus\t"<<
      "+SystEff\t"<<"-SystEff\t"<<
      endl;
    
    for(int i(1);i<nBins12;i++)
    {
      cout<<i<<
	"\t+ "<<ZptSystBkg_plus12[i] <<"\t- "<< ZptSystBkg_minus12[i]<<" | "<<
	"\t+ "<<ZptUnfoldingShapeSyst_plus12[i] <<"\t- "<< ZptUnfoldingShapeSyst_minus12[i]<<" | "<<
	"\t+ "<<ZptUnfoldSyst_doubleGaussin_plus12[i] <<"\t- "<< ZptUnfoldSyst_doubleGaussin_minus12[i]<<" | "<<
	"\t+ "<<ZptSystEff_plus12[i] <<"\t- "<< ZptSystEff_minus12[i]<<
	endl;
    }

    cout<<"=================================================================="<<endl;
    cout<<"Bin:\t"<<
      "SystFSR\t\t"<<
      "+SystMadgraph\t"<<"-SystMadgraph\t"<<
      endl;
    for(int i(1);i<nBins12;i++)
    {
      cout<<i<<
	"\t"<<ZptSystFSRErr12[i]<<" | "<<
	"\t+ "<<ZptMadgraphSyst_plus12[i] <<"\t- "<< ZptMadgraphSyst_minus12[i]<<" | "<<
	endl;
    }
    cout << "======= Put in note=======" << endl;
    cout << fixed << setprecision(4);
    for(int i(1);i<nBins12;i++)
    {
      cout<<
	" & "<<ZptMadgraphSyst_plus12[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    cout << fixed << setprecision(8);


    cout<<"=================================================================="<<endl;
    cout<<" 12 Bin Zpt Statistic Error in %:"<<endl;
    cout<<"Bin:\t"<<"Statistic error in %"<< endl;
    for(int i(1);i<nBins12;i++)
    {
      cout<<i<<	"\t"<<ZptStatErr12[i]<<endl;
    }

    cout<<"=================================================================="<<endl;
    cout<<" 12 Bin Zpt Total syst Errors in %:"<<endl;
    cout<<"Bin:\t"<<"systplus12\t"<<"systminus12"<<endl;
    double systplus12[13]={0};
    double systminus12[13]={0};
    for(int i(1);i<nBins12;i++)
    {
      systplus12[i]  = TMath::Sqrt(syst2plus12[i] );
      systminus12[i] = TMath::Sqrt(syst2minus12[i]);
      cout<<i<<"\t+"<< systplus12[i]<<"\t-"<< systminus12[i]<<endl;
    }
    cout << "======= Put in note=======" << endl;
    cout << fixed << setprecision(4);
    for(int i(1);i<nBins12;i++)
    {
      cout<<
	" & "<< systplus12[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    for(int i(1);i<nBins12;i++)
    {
      cout<<
	" & "<< systminus12[i] ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
   cout << fixed << setprecision(8);



    double errorplus12[13]={0};
    double errorminus12[13]={0};
    cout<<"=================================================================="<<endl;
    cout<<" 12 Bin Zpt Total Uncertainty in %:"<<endl;
    cout<<"Bin:\t"<<"errorplus12\t"<<"errorminus12"<<endl;
    for(int i(1);i<nBins12;i++)
    {
      errorplus12[i]= TMath::Sqrt( TMath::Power(ZptStatErr12[i],2) + TMath::Power(systplus12[i] ,2) );
      errorminus12[i]= TMath::Sqrt( TMath::Power(ZptStatErr12[i] ,2) + TMath::Power(systminus12[i],2) );
      cout<<i<<"\t+"<<errorplus12[i]<<"\t-"<<errorminus12[i]<<endl;
    }

    cout << "======= Put in note=======" << endl;
    cout << fixed << setprecision(4);
    for(int i(1);i<nBins12;i++)
    {
      cout<<
	" & "<<TMath::Max( errorplus12[i], errorminus12[i]) ;
    }
   cout << " \\"<<"\\" << endl;
   cout << endl;
    cout << fixed << setprecision(8);
   



    //// Make cross-section from N of Yield
    TH1D* hZptXsec12Bin = (TH1D*)hZptYield12Bin->Clone("hZptXsec12Bin");  
        //// Set bin errors to 12 Bin cross sectio: hZptXsec12Bin
	for(int i(1);i<nBins12;i++)
	{
	  hZptXsec12Bin->SetBinError(i, TMath::Max(errorplus12[i],errorminus12[i])*0.01*hZptYield12Bin->GetBinContent(i) );
	}
    //// Divide by Luminocity = Make it cross section now
    hZptXsec12Bin->Scale(1./18.429);

    //// make Differential cross-section from cross-section
    TH1D* hZptDiffXsec12Bin = (TH1D*)hZptXsec12Bin->Clone("hZptDiffXsec12Bin");  
	for(int i(1);i<nBins12;i++)
	{
	  hZptDiffXsec12Bin->SetBinContent(i, hZptXsec12Bin->GetBinContent(i)/hWptBins_LinScale->GetXaxis()->GetBinWidth(i) );
	  hZptDiffXsec12Bin->SetBinError(i, hZptXsec12Bin->GetBinError(i)/hWptBins_LinScale->GetXaxis()->GetBinWidth(i) );
	}
    
    //// Normalize Differential cross-section with Total cross-section  in Fid volume
    TH1D* hZptDiffXsec12BinInFidNorm = (TH1D*)hZptDiffXsec12Bin->Clone("hZptDiffXsec12BinInFidNorm");  
    hZptDiffXsec12BinInFidNorm->Scale(1./hZptXsec12Bin->Integral());
   



    cout<<"=================================================================="<<endl;
    cout<<" 12 Bin Zpt Differential cross-section in Fiducial volume"<<endl;
    cout<<"Bin:\t"<<"Value\t\t"<<" Total Uncertainty error\t"<<" Total Uncertainty error in %"<<endl;
    for(int i(1);i<nBins12;i++)
    {
      cout<<i<<"\t"<<hZptDiffXsec12Bin->GetBinContent(i)<<"\t+- "<<hZptDiffXsec12Bin->GetBinError(i)<<"\t\t\t+- "<<TMath::Max(errorplus12[i],errorminus12[i])<<endl;
    }
    
    cout<<"=================================================================="<<endl;
    cout<<" 12 Bin Zpt Normalized Differential cross-section in Fiducial volume"<<endl;
    cout<<"Bin:\t"<<"Value\t\t"<<" Total Uncertainty error\t"<<" Total Uncertainty error in %"<<endl;
    for(int i(1);i<nBins12;i++)
    {
      cout<<i<<"\t"<<hZptDiffXsec12BinInFidNorm->GetBinContent(i)<<"\t+- "<<hZptDiffXsec12BinInFidNorm->GetBinError(i)<<"\t\t\t+- "<<TMath::Max(errorplus12[i],errorminus12[i])<<endl;
    }


//===================Write to root file ========================================================================
      
  TH1D *hZptSystBkg_plus = new TH1D("hZptSystBkg_plus","hZptSystBkg_plus",nBins12-1,Wpt12Bins);hZptSystBkg_plus->Sumw2(); 
  TH1D *hZptSystBkg_minus = new TH1D("hZptSystBkg_minus","hZptSystBkg_minus",nBins12-1,Wpt12Bins);hZptSystBkg_minus->Sumw2(); 

  TH1D *hZptSystEff_plus = new TH1D("hZptSystEff_plus","hZptSystEff_plus",nBins12-1,Wpt12Bins);hZptSystEff_plus->Sumw2(); 
  TH1D *hZptSystEff_minus = new TH1D("hZptSystEff_minus","hZptSystEff_minus",nBins12-1,Wpt12Bins);hZptSystEff_minus->Sumw2(); 
  
  TH1D *hZptSystFSRErr = new TH1D("hZptSystFSRErr","hZptSystFSRErr",nBins12-1,Wpt12Bins);hZptSystFSRErr->Sumw2(); 

  TH1D *hZptUnfoldSyst_doubleGaussin_plus = new TH1D("hZptUnfoldSyst_doubleGaussin_plus","hZptUnfoldSyst_doubleGaussin_plus",nBins12-1,Wpt12Bins);hZptUnfoldSyst_doubleGaussin_plus->Sumw2(); 
  TH1D *hZptUnfoldSyst_doubleGaussin_minus = new TH1D("hZptUnfoldSyst_doubleGaussin_minus","hZptUnfoldSyst_doubleGaussin_minus",nBins12-1,Wpt12Bins);hZptUnfoldSyst_doubleGaussin_minus->Sumw2(); 

  TH1D *hZptUnfoldingShapeSyst_plus = new TH1D("hZptUnfoldingShapeSyst_plus","hZptUnfoldingShapeSyst_plus",nBins12-1,Wpt12Bins);hZptUnfoldingShapeSyst_plus->Sumw2(); 
  TH1D *hZptUnfoldingShapeSyst_minus = new TH1D("hZptUnfoldingShapeSyst_minus","hZptUnfoldingShapeSyst_minus",nBins12-1,Wpt12Bins);hZptUnfoldingShapeSyst_minus->Sumw2(); 
  
  TH1D *hZptMadgraphSyst_plus = new TH1D("hZptMadgraphSyst_plus","hZptMadgraphSyst_plus",nBins12-1,Wpt12Bins);hZptMadgraphSyst_plus->Sumw2(); 
  TH1D *hZptMadgraphSyst_minus = new TH1D("hZptMadgraphSyst_minus","hZptMadgraphSyst_minus",nBins12-1,Wpt12Bins);hZptMadgraphSyst_minus->Sumw2(); 

  TH1D *hZptStatErr = new TH1D("hZptStatErr","hZptStatErr",nBins12-1,Wpt12Bins);hZptStatErr->Sumw2(); 
  
  TH1D *hZptTotalSyst_plus = new TH1D("hZptTotalSyst_plus","hZptTotalSyst_plus",nBins12-1,Wpt12Bins);hZptTotalSyst_plus->Sumw2(); 
  TH1D *hZptTotalSyst_minus = new TH1D("hZptTotalSyst_minus","hZptTotalSyst_minus",nBins12-1,Wpt12Bins);hZptTotalSyst_minus->Sumw2();

  TH1D *hZptTotalUncertainty_plus = new TH1D("hZptTotalUncertainty_plus","hZptTotalUncertainty_plus",nBins12-1,Wpt12Bins);hZptTotalUncertainty_plus->Sumw2(); 
  TH1D *hZptTotalUncertainty_minus = new TH1D("hZptTotalUncertainty_minus","hZptTotalUncertainty_minus",nBins12-1,Wpt12Bins);hZptTotalUncertainty_minus->Sumw2();

	for(int i(1);i<nBins12;i++)
	{
	  hZptSystBkg_plus->SetBinContent(i,ZptSystBkg_plus12[i]);
	  hZptSystBkg_minus->SetBinContent(i,ZptSystBkg_minus12[i]);
	  
	  hZptSystEff_plus->SetBinContent(i,ZptSystEff_plus12[i]);
	  hZptSystEff_minus->SetBinContent(i,ZptSystEff_minus12[i]);
	  
	  hZptSystFSRErr->SetBinContent(i,ZptSystFSRErr12[i]);
	  
	  hZptUnfoldSyst_doubleGaussin_plus->SetBinContent(i,ZptUnfoldSyst_doubleGaussin_plus12[i]);
	  hZptUnfoldSyst_doubleGaussin_minus->SetBinContent(i,ZptUnfoldSyst_doubleGaussin_minus12[i]);
	  
	  hZptUnfoldingShapeSyst_plus->SetBinContent(i,ZptUnfoldingShapeSyst_plus12[i]);
	  hZptUnfoldingShapeSyst_minus->SetBinContent(i,ZptUnfoldingShapeSyst_minus12[i]);
	  
	  hZptMadgraphSyst_plus->SetBinContent(i,ZptMadgraphSyst_plus12[i]);
	  hZptMadgraphSyst_minus->SetBinContent(i,ZptMadgraphSyst_minus12[i]);
	  
	  hZptStatErr->SetBinContent(i,ZptStatErr12[i]);
	  
	  hZptTotalSyst_plus->SetBinContent(i,systplus12[i]);
	  hZptTotalSyst_minus->SetBinContent(i,systminus12[i]);
	  
	  hZptTotalUncertainty_plus->SetBinContent(i,errorplus12[i]);
	  hZptTotalUncertainty_minus->SetBinContent(i,errorminus12[i]);
	}

  TString resultDir = "ZptXsecErrors";
  gSystem->mkdir(resultDir,kTRUE);
  //TFile f_out(resultDir+"/ZptXsecErrors_FidVolume.root","recreate");
  TFile f_out(resultDir+"/ZptXsecErrors_FidVolume_MergeChooseBig.root","recreate");

  hZptSystBkg_plus->Write();
  hZptSystBkg_minus->Write();
  hZptSystEff_plus->Write();
  hZptSystEff_minus->Write();
  hZptSystFSRErr->Write();
  hZptUnfoldSyst_doubleGaussin_plus->Write();
  hZptUnfoldSyst_doubleGaussin_minus->Write();
  hZptUnfoldingShapeSyst_plus->Write();
  hZptUnfoldingShapeSyst_minus->Write();
  
  hZptMadgraphSyst_plus->Write();
  hZptMadgraphSyst_minus->Write();
  
  hZptStatErr->Write();
  
  hZptTotalSyst_plus->Write();
  hZptTotalSyst_minus->Write();
  hZptTotalUncertainty_plus->Write();
  hZptTotalUncertainty_minus->Write();

  //In Fiducial volume 
  hZptDiffXsec12Bin->Write(); 
  //Normalized In Fiducial volume 
  hZptDiffXsec12BinInFidNorm->Write();
 
}
