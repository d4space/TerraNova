{
#include "../../Utils/const.h"
#include "TMathBase.h"
double BinWidth[14] ={0.0, 7.5-0, 12.5-7.5,17.5-12.5, 24.0-17.5, 30.0-24.0, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};
double Diff_Xsec[14] = {0.,};
Diff_Xsec[1] = 306.833    ;
Diff_Xsec[2] = 262.008    ;
Diff_Xsec[3] = 149.336    ;
Diff_Xsec[4] = 89.8894    ;
Diff_Xsec[5] = 58.2311    ;
Diff_Xsec[6] = 37.2108    ;
Diff_Xsec[7] = 22.0426    ;
Diff_Xsec[8] = 11.1543    ;
Diff_Xsec[9] = 3.38976    ;
Diff_Xsec[10] =0.816256   ;
Diff_Xsec[11] =0.264671   ;
Diff_Xsec[12] =0.0881628  ;
Diff_Xsec[13] =0.00722373 ;

double Diff_Xsec_Wp[14] = {0.,};
Diff_Xsec_Wp[1] = 185.589   ;
Diff_Xsec_Wp[2] = 156.357   ;
Diff_Xsec_Wp[3] = 88.2991   ;
Diff_Xsec_Wp[4] = 53.2227   ;
Diff_Xsec_Wp[5] = 33.9276   ;
Diff_Xsec_Wp[6] = 21.6567   ;
Diff_Xsec_Wp[7] = 12.7669   ;
Diff_Xsec_Wp[8] = 6.46915   ;
Diff_Xsec_Wp[9] = 1.96717   ;
Diff_Xsec_Wp[10] =0.474284  ;
Diff_Xsec_Wp[11] =0.151101  ;
Diff_Xsec_Wp[12] =0.0524342 ;
Diff_Xsec_Wp[13] =0.00445951;

double Diff_Xsec_Wm[14] = {0.,};
Diff_Xsec_Wm[1] = 121.244    ;
Diff_Xsec_Wm[2] = 105.651    ;
Diff_Xsec_Wm[3] = 61.0366    ;
Diff_Xsec_Wm[4] = 36.6667    ;
Diff_Xsec_Wm[5] = 24.3035    ;
Diff_Xsec_Wm[6] = 15.5541    ;
Diff_Xsec_Wm[7] = 9.27561    ;
Diff_Xsec_Wm[8] = 4.68518    ;
Diff_Xsec_Wm[9] = 1.42259    ;
Diff_Xsec_Wm[10] =0.341972   ;
Diff_Xsec_Wm[11] =0.11357    ;
Diff_Xsec_Wm[12] =0.0357286  ;
Diff_Xsec_Wm[13] =0.00276423 ;

//*
  double mutracksigErri[14]={0};
  mutracksigErri[1]= 0.0347;
  mutracksigErri[2]= 0.0047;
  mutracksigErri[3]= 0.0315;
  mutracksigErri[4]= 0.0495;
  mutracksigErri[5]= 0.0465;
  mutracksigErri[6]= 0.0306;
  mutracksigErri[7]= 0.0157;
  mutracksigErri[8]= 0.0206;
  mutracksigErri[9]= 0.0377;
  mutracksigErri[10]=0.0536;
  mutracksigErri[11]=0.0646;
  mutracksigErri[12]=0.0768;
  mutracksigErri[13]=0.0851;

  double mutracksigErriMerge[14]={0};
  //mutracksigErriMerge[4]=TMath::Max(mutracksigErri[4],mutracksigErri[5]);
  //mutracksigErriMerge[4]=(BinWidth[4]*mutracksigErri[4] + BinWidth[5]*mutracksigErri[5]) / (BinWidth[4] + BinWidth[5]);
  mutracksigErriMerge[4]=(BinWidth[4]*mutracksigErri[4]*Diff_Xsec[4] + BinWidth[5]*mutracksigErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double mutrackbckErri[14] = {0};
  mutrackbckErri[1 ]= 0.0227;
  mutrackbckErri[2 ]= 0.0101;
  mutrackbckErri[3 ]= 0.0289;
  mutrackbckErri[4 ]= 0.0332;
  mutrackbckErri[5 ]= 0.0302;
  mutrackbckErri[6 ]= 0.0298;
  mutrackbckErri[7 ]= 0.0291;
  mutrackbckErri[8 ]= 0.0224;
  mutrackbckErri[9 ]= 0.0125;
  mutrackbckErri[10]= 0.0051;
  mutrackbckErri[11]= 0.0074;
  mutrackbckErri[12]= 0.0152;
  mutrackbckErri[13]= 0.0215;

  double mutrackbckErriMerge[14]={0};
  //mutrackbckErriMerge[4]=TMath::Max(mutrackbckErri[4],mutrackbckErri[5]);
  //mutrackbckErriMerge[4]=(BinWidth[4]*mutrackbckErri[4] + BinWidth[5]*mutrackbckErri[5]) / (BinWidth[4] + BinWidth[5]);
  mutrackbckErriMerge[4]=(BinWidth[4]*mutrackbckErri[4]*Diff_Xsec[4] + BinWidth[5]*mutrackbckErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);


  double muidisosigErri[14] = {0};
  muidisosigErri[1 ] = 0.0243;
  muidisosigErri[2 ] = 0.0066;
  muidisosigErri[3 ] = 0.0390;
  muidisosigErri[4 ] = 0.0478;
  muidisosigErri[5 ] = 0.0321;
  muidisosigErri[6 ] = 0.0099;
  muidisosigErri[7 ] = 0.0216;
  muidisosigErri[8 ] = 0.0358;
  muidisosigErri[9 ] = 0.0452;
  muidisosigErri[10] = 0.0517;
  muidisosigErri[11] = 0.0572;
  muidisosigErri[12] = 0.0586;
  muidisosigErri[13] = 0.0588;

  double muidisosigErriMerge[14]={0};
  //muidisosigErriMerge[4]=TMath::Max(muidisosigErri[4],muidisosigErri[5]);
  //muidisosigErriMerge[4]=(BinWidth[4]*muidisosigErri[4] + BinWidth[5]*muidisosigErri[5]) / (BinWidth[4] + BinWidth[5]);
  //muidisosigErriMerge[4]=(BinWidth[4]*muidisosigErri[4] + BinWidth[5]*muidisosigErri[5]) / (BinWidth[4] + BinWidth[5]);
  muidisosigErriMerge[4]=(BinWidth[4]*muidisosigErri[4]*Diff_Xsec[4] + BinWidth[5]*muidisosigErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double muidisobckErri[14]={0};
  muidisobckErri[1 ] = 0.0211;
  muidisobckErri[2 ] = 0.0135;
  muidisobckErri[3 ] = 0.0405;
  muidisobckErri[4 ] = 0.0455;
  muidisobckErri[5 ] = 0.0249;
  muidisobckErri[6 ] = 0.0210;
  muidisobckErri[7 ] = 0.0461;
  muidisobckErri[8 ] = 0.0664;
  muidisobckErri[9 ] = 0.0807;
  muidisobckErri[10] = 0.0910;
  muidisobckErri[11] = 0.0984;
  muidisobckErri[12] = 0.1037;
  muidisobckErri[13] = 0.1074;

  double muidisobckErriMerge[14]={0};
  //muidisobckErriMerge[4]=TMath::Max(muidisobckErri[4],muidisobckErri[5]);
  //muidisobckErriMerge[4]=(BinWidth[4]*muidisobckErri[4] + BinWidth[5]*muidisobckErri[5]) / (BinWidth[4] + BinWidth[5]);
  muidisobckErriMerge[4]=(BinWidth[4]*muidisobckErri[4]*Diff_Xsec[4] + BinWidth[5]*muidisobckErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double mutoyErri[14] = {0};
  mutoyErri[1 ] = 0.190374893;
  mutoyErri[2 ] = 0.095197952;
  mutoyErri[3 ] = 0.123205722;
  mutoyErri[4 ] = 0.250625877;
  mutoyErri[5 ] = 0.321657116;
  mutoyErri[6 ] = 0.33107028 ;
  mutoyErri[7 ] = 0.294405112;
  mutoyErri[8 ] = 0.2508887  ;
  mutoyErri[9 ] = 0.228813505;
  mutoyErri[10] = 0.231403652;
  mutoyErri[11] = 0.244196908;
  mutoyErri[12] = 0.262318661;
  mutoyErri[13] = 0.275364195;

  double mutoyErriMerge[14]={0};
  //mutoyErriMerge[4]=TMath::Max(mutoyErri[4],mutoyErri[5]);
  //mutoyErriMerge[4]=(BinWidth[4]*mutoyErri[4] + BinWidth[5]*mutoyErri[5])/ (BinWidth[4] + BinWidth[5]);
  mutoyErriMerge[4]=(BinWidth[4]*mutoyErri[4]*Diff_Xsec[4] + BinWidth[5]*mutoyErri[5]*Diff_Xsec[5])/ (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double muPOGErri[14] = {0};
  muPOGErri[1 ] = 0.08744747;
  muPOGErri[2 ] = 0.05771940;
  muPOGErri[3 ] = 0.10647995;
  muPOGErri[4 ] = 0.18673989;
  muPOGErri[5 ] = 0.23551753;
  muPOGErri[6 ] = 0.18493942;
  muPOGErri[7 ] = 0.09281697;
  muPOGErri[8 ] = 0.12348279;
  muPOGErri[9 ] = 0.19797320;
  muPOGErri[10] = 0.24525617;
  muPOGErri[11] = 0.27925962;
  muPOGErri[12] = 0.28776855;
  muPOGErri[13] = 0.31057025;

  double muPOGErriMerge[14]={0};
  //muPOGErriMerge[4]=TMath::Max(muPOGErri[4],muPOGErri[5]);
  //muPOGErriMerge[4]=(BinWidth[4]*muPOGErri[4] + BinWidth[5]*muPOGErri[5]) / (BinWidth[4] + BinWidth[5]);
  muPOGErriMerge[4]=(BinWidth[4]*muPOGErri[4]*Diff_Xsec[4] + BinWidth[5]*muPOGErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double muEffErri[14]={0};
  for(int i(1);i<14;i++)
  {
    muEffErri[i]  = sqrt(
	mutracksigErri[i] *mutracksigErri[i] 
	+mutrackbckErri[i] * mutrackbckErri[i] 
	+mutoyErri[i] *mutoyErri[i] 
	+muidisosigErri[i] *muidisosigErri[i] 
	+muidisobckErri[i]* muidisobckErri[i] 
	+muPOGErri[i] *muPOGErri[i]);
  } 
  
  double muEffErriMerge[14]={0};
  muEffErriMerge[4]  = sqrt(
      mutracksigErriMerge[4] *mutracksigErriMerge[4]
      +mutrackbckErriMerge[4]*mutrackbckErriMerge[4] 
      +mutoyErriMerge[4] *mutoyErriMerge[4] 
      +muidisosigErriMerge[4] *muidisosigErriMerge[4]
      +muidisobckErriMerge[4]*muidisobckErriMerge[4] 
      +muPOGErriMerge[4] *muPOGErriMerge[4]);

  double muStatErri[13] = {0};
  muStatErri[1 ]= 0.5112 ; 
  muStatErri[2 ]= 0.6499 ;
  muStatErri[3 ]= 0.7834 ;
  muStatErri[4 ]= 0.6807 ;
  muStatErri[5 ]= 1.1344 ;
  muStatErri[6 ]= 1.5846 ;
  muStatErri[7 ]= 1.5717 ;
  muStatErri[8 ]= 2.0272 ;
  muStatErri[9 ]= 4.0856 ;
  muStatErri[10]= 7.9529 ;
  muStatErri[11]= 12.7258;
  muStatErri[12]= 19.6641;

  double muMetErri[14] = {0};
  muMetErri[1 ] = 0.035468296;
  muMetErri[2 ] = 0.01959209 ;
  muMetErri[3 ] = 0.040611821;
  muMetErri[4 ] = 0.060981965;
  muMetErri[5 ] = 0.063530308;
  muMetErri[6 ] = 0.062677348;
  muMetErri[7 ] = 0.063817944;
  muMetErri[8 ] = 0.072547433;
  muMetErri[9 ] = 0.093860588;
  muMetErri[10] = 0.123264675;
  muMetErri[11] = 0.152265065;
  muMetErri[12] = 0.17241914 ;
  muMetErri[13] = 0.182738721;

  double muMetErriMerge[14] = {0};
  //muMetErriMerge[4]=TMath::Max(muMetErri[4 ],muMetErri[5 ]);
  //muMetErriMerge[4]=(BinWidth[4]*muMetErri[4 ] + BinWidth[5]*muMetErri[5]) / (BinWidth[4] + BinWidth[5]);
  muMetErriMerge[4]=(BinWidth[4]*muMetErri[4 ]*Diff_Xsec[4] + BinWidth[5]*muMetErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);


  double muscaleErri[14]={0};
  muscaleErri[1] = 0.1751;
  muscaleErri[2] = 0.1069;
  muscaleErri[3] = 0.0573;
  muscaleErri[4] = 0.1907;
  muscaleErri[5] = 0.2999;
  muscaleErri[6] = 0.3399;
  muscaleErri[7] = 0.3249;
  muscaleErri[8] = 0.2882;
  muscaleErri[9] = 0.2489;
  muscaleErri[10]= 0.2155;
  muscaleErri[11]= 0.2340;
  muscaleErri[12]= 0.2637;
  muscaleErri[13]= 0.2794;

  double muscaleErriMerge[14]={0};
  //muscaleErriMerge[4] = TMath::Max(muscaleErri[4],muscaleErri[5]);
  //muscaleErriMerge[4] = (BinWidth[4]*muscaleErri[4] + BinWidth[5]*muscaleErri[5]) / (BinWidth[4] + BinWidth[5]);
  muscaleErriMerge[4] = (BinWidth[4]*muscaleErri[4]*Diff_Xsec[4] + BinWidth[5]*muscaleErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double musmearErri[14]={0};
  musmearErri[1] = 0.1587;
  musmearErri[2] = 0.0816;
  musmearErri[3] = 0.0722;
  musmearErri[4] = 0.2135;
  musmearErri[5] = 0.2908;
  musmearErri[6] = 0.2855;
  musmearErri[7] = 0.2282;
  musmearErri[8] = 0.1728;
  musmearErri[9] = 0.1610;
  musmearErri[10]= 0.2093;
  musmearErri[11]= 0.2533;
  musmearErri[12]= 0.2864;
  musmearErri[13]= 0.3062;

  double musmearErriMerge[14]={0};
  //musmearErriMerge[4] = TMath::Max(musmearErri[4],musmearErri[5]);
  //musmearErriMerge[4] = (BinWidth[4]*musmearErri[4] + BinWidth[5]*musmearErri[5]) / (BinWidth[4] + BinWidth[5]);
  musmearErriMerge[4] = (BinWidth[4]*musmearErri[4]*Diff_Xsec[4] + BinWidth[5]*musmearErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double muMomResErri[14]={0};
  muMomResErri[1] = 0.236317;
  muMomResErri[2] = 0.134485;
  muMomResErri[3] = 0.092174;
  muMomResErri[4] = 0.286267;
  muMomResErri[5] = 0.417738;
  muMomResErri[6] = 0.443894;
  muMomResErri[7] = 0.397033;
  muMomResErri[8] = 0.336034;
  muMomResErri[9] = 0.296432;
  muMomResErri[10]= 0.300411;
  muMomResErri[11]= 0.344843;
  muMomResErri[12]= 0.38931;
  muMomResErri[13]= 0.414515;

  double muMomResErriMerge[14]={0};
  //muMomResErriMerge[4]  = TMath::Max(muMomResErri[4],muMomResErri[5]);
  //muMomResErriMerge[4]  = (BinWidth[4]*muMomResErri[4] + BinWidth[5]*muMomResErri[5]) / (BinWidth[4] + BinWidth[5]);
  muMomResErriMerge[4]  = (BinWidth[4]*muMomResErri[4]*Diff_Xsec[4] + BinWidth[5]*muMomResErri[5]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);
  
  double muQCDBckErri[14]={0};
  muQCDBckErri[1 ] = 0.6246;
  muQCDBckErri[2 ] = 0.9519;
  muQCDBckErri[3 ] = 0.8705;
  muQCDBckErri[4 ] = 0.9374;
  muQCDBckErri[5 ] = 0.9376;
  muQCDBckErri[6 ] = 1.5180;
  muQCDBckErri[7 ] = 0.8863;
  muQCDBckErri[8 ] = 1.4678;
  muQCDBckErri[9 ] = 0.6777; 
  muQCDBckErri[10] = 0.6800;
  muQCDBckErri[11] = 0.6951;
  muQCDBckErri[12] = 0.7067;
  muQCDBckErri[13] = 0.7198;

  double muQCDBckErriMerge[14]={0};
   //muQCDBckErriMerge[4]=(BinWidth[4]*muQCDBckErri[4 ] + BinWidth[5]*muQCDBckErri[5 ]) / (BinWidth[4] + BinWidth[5]);
   muQCDBckErriMerge[4]=(BinWidth[4]*muQCDBckErri[4 ]*Diff_Xsec[4] + BinWidth[5]*muQCDBckErri[5 ]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);


  double muQCDShapeErri[14] = {0};
  muQCDShapeErri[1 ] = 0.1654;
  muQCDShapeErri[2 ] = 0.2559;
  muQCDShapeErri[3 ] = 0.2184;
  muQCDShapeErri[4 ] = 0.2669;
  muQCDShapeErri[5 ] = 0.2800;
  muQCDShapeErri[6 ] = 0.2574;
  muQCDShapeErri[7 ] = 0.1527;
  muQCDShapeErri[8 ] = 0.3100;
  muQCDShapeErri[9 ] = 0.2504; 
  muQCDShapeErri[10] = 0.1431;
  muQCDShapeErri[11] = 0.6232;
  muQCDShapeErri[12] = 0.6740;
  muQCDShapeErri[13] = 0.6749;

  double muQCDShapeErriMerge[14] = {0};
  //muQCDShapeErriMerge[4]=TMath::Max(muQCDShapeErri[4 ],muQCDShapeErri[5 ]);
  //muQCDShapeErriMerge[4]=(BinWidth[4]*muQCDShapeErri[4 ] + BinWidth[5]*muQCDShapeErri[5 ]) / (BinWidth[4] + BinWidth[5]);
  muQCDShapeErriMerge[4]=(BinWidth[4]*muQCDShapeErri[4 ]*Diff_Xsec[4] + BinWidth[5]*muQCDShapeErri[5 ]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);
  
  double muEWKErri[14] ={0};
  muEWKErri[1 ] = 0.004762989;
  muEWKErri[2 ] = 0.015008568;
  muEWKErri[3 ] = 0.032965654;
  muEWKErri[4 ] = 0.036521741;
  muEWKErri[5 ] = 0.019697343;
  muEWKErri[6 ] = 0.027481851;
  muEWKErri[7 ] = 0.063564663;
  muEWKErri[8 ] = 0.096868134;
  muEWKErri[9 ] = 0.124563018; 
  muEWKErri[10] = 0.14618646 ;
  muEWKErri[11] = 0.162271334;
  muEWKErri[12] = 0.173472133;
  muEWKErri[13] = 0.179434836;

  double muEWKErriMerge[14] ={0};
  //muEWKErriMerge[4]=TMath::Max(muEWKErri[4 ] , muEWKErri[5 ] ); 
  //muEWKErriMerge[4]=(BinWidth[4]*muEWKErri[4 ] + BinWidth[5]*muEWKErri[5 ] ) / (BinWidth[4] + BinWidth[5]); 
  muEWKErriMerge[4]=(BinWidth[4]*muEWKErri[4 ]*Diff_Xsec[4] + BinWidth[5]*muEWKErri[5 ]*Diff_Xsec[5] ) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]); 
  
  double muFSRErri[14] ={0};
  muFSRErri[1 ] = 0.0020 ;
  muFSRErri[2 ] = 0.0030 ;
  muFSRErri[3 ] = 0.0032 ;
  muFSRErri[4 ] = 0.0013 ;
  muFSRErri[5 ] = 0.0021 ;
  muFSRErri[6 ] = 0.0060 ;
  muFSRErri[7 ] = 0.0099 ;
  muFSRErri[8 ] = 0.0133 ;
  muFSRErri[9 ] = 0.0159 ;
  muFSRErri[10] =  0.0180;
  muFSRErri[11] =  0.0194;
  muFSRErri[12] =  0.0203;
  muFSRErri[13] =  0.0207;

  double muFSRErriMerge[14] ={0};
  //muFSRErriMerge[4]=(BinWidth[4]*muFSRErri[4 ] + BinWidth[5]* muFSRErri[5 ]) / (BinWidth[4] + BinWidth[5]);
  muFSRErriMerge[4]=(BinWidth[4]*muFSRErri[4 ]*Diff_Xsec[4] + BinWidth[5]* muFSRErri[5 ]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  double musvdunfp[14] ={0};
  musvdunfp[1 ] = 0.2036;
  musvdunfp[2 ] = 0.1786;
  musvdunfp[3 ] = 0.2121;
  musvdunfp[4 ] = 0.2623;
  musvdunfp[5 ] = 0.2858;
  musvdunfp[6 ] = 0.2961;
  musvdunfp[7 ] = 0.3214;
  musvdunfp[8 ] = 0.3971;
  musvdunfp[9 ] = 0.5205;
  musvdunfp[10] = 0.6648;
  musvdunfp[11] = 0.8008;
  musvdunfp[12] = 0.9051;
  musvdunfp[13] = 0.9609;
  
  double musvdunfm[14] ={0};
  musvdunfm[1 ] = 0.1421;
  musvdunfm[2 ] = 0.1275;
  musvdunfm[3 ] = 0.1514;
  musvdunfm[4 ] = 0.1769;
  musvdunfm[5 ] = 0.1907;
  musvdunfm[6 ] = 0.2037;
  musvdunfm[7 ] = 0.2255;
  musvdunfm[8 ] = 0.2802;
  musvdunfm[9 ] = 0.3717;
  musvdunfm[10] = 0.4873;
  musvdunfm[11] = 0.6025;
  musvdunfm[12] = 0.6930;
  musvdunfm[13] = 0.7420;

  double musvdunfpMerge[14] ={0};
  double musvdunfmMerge[14] ={0};
  //musvdunfpMerge[4]=TMath::Max(musvdunfp[4 ], musvdunfp[5 ]);
  //musvdunfmMerge[4]=TMath::Max(musvdunfm[4 ], musvdunfm[5 ]);
  //musvdunfpMerge[4]=(BinWidth[4]*musvdunfp[4 ] + BinWidth[5]* musvdunfp[5 ])/ (BinWidth[4] + BinWidth[5]);
  //musvdunfmMerge[4]=(BinWidth[4]*musvdunfm[4 ] + BinWidth[5]* musvdunfm[5 ])/ (BinWidth[4] + BinWidth[5]);
  musvdunfpMerge[4]=(BinWidth[4]*musvdunfp[4 ]*Diff_Xsec_Wp[4] + BinWidth[5]* musvdunfp[5 ]*Diff_Xsec_Wp[5])/ (BinWidth[4]*Diff_Xsec_Wp[4] + BinWidth[5]*Diff_Xsec_Wp[5]);
  musvdunfmMerge[4]=(BinWidth[4]*musvdunfm[4 ]*Diff_Xsec_Wm[4] + BinWidth[5]* musvdunfm[5 ]*Diff_Xsec_Wm[5])/ (BinWidth[4]*Diff_Xsec_Wm[4] + BinWidth[5]*Diff_Xsec_Wm[5]);

  double muUnfBiasErri[14] = {0};
  muUnfBiasErri[1 ] =0.1790   ;//0.9022;
  muUnfBiasErri[2 ] =0.1820   ;//0.6571;
  muUnfBiasErri[3 ] =0.1319   ;//1.1301;
  muUnfBiasErri[4 ] =0.0370   ;//3.3043;
  muUnfBiasErri[5 ] =0.1919   ;//2.7469;
  muUnfBiasErri[6 ] =0.3983   ;//3.0764;
  muUnfBiasErri[7 ] =0.5934   ;//3.2467;
  muUnfBiasErri[8 ] =0.7881   ;//5.5057;
  muUnfBiasErri[9 ] =0.9423   ;//5.3121;
  muUnfBiasErri[10] = 1.0577 ;//1.8476;
  muUnfBiasErri[11] = 1.2444 ;//1.6692;
  muUnfBiasErri[12] = 1.2193 ;//4.4869;
  muUnfBiasErri[13] = 1.4388 ;//5.9061;

  double muUnfBiasErriMerge[14] = {0};
  //muUnfBiasErriMerge[4]=TMath::Max(muUnfBiasErri[4 ] , muUnfBiasErri[5 ]);
  //muUnfBiasErriMerge[4]=(BinWidth[4]*muUnfBiasErri[4 ] + BinWidth[5]* muUnfBiasErri[5 ]) / (BinWidth[4] + BinWidth[5]);
  muUnfBiasErriMerge[4]=(BinWidth[4]*muUnfBiasErri[4 ]*Diff_Xsec[4] + BinWidth[5]* muUnfBiasErri[5 ]*Diff_Xsec[5]) / (BinWidth[4]*Diff_Xsec[4] + BinWidth[5]*Diff_Xsec[5]);

  // Make Incl
  //*
  TFile *fp = new TFile("../ResultWpToMuNu/Result_WpToMuNu.root");
  TFile *fm = new TFile("../ResultWmToMuNu/Result_WmToMuNu.root");
  TH1D *h_data_p;
  TH1D *h_data_m;
  TH1D *h_MC_p;
  TH1D *h_MC_m;
  TH1D *h_MC_p_UnWeighted;
  TH1D *h_MC_m_UnWeighted;
  TH1D *h_dataRec_p;
  TH1D *h_dataRec_m;
  h_data_p = (TH1D*)fp->Get("BornEffCorr")->Clone("h_data_p");
  h_data_m = (TH1D*)fm->Get("BornEffCorr")->Clone("h_data_m");
  h_MC_p = (TH1D*)fp->Get("SVD_Born.Gen")->Clone("h_MC_p");
  h_MC_m = (TH1D*)fm->Get("SVD_Born.Gen")->Clone("h_MC_m");
  h_MC_p_UnWeighted = (TH1D*)fp->Get("SVD_Born.Gen")->Clone("h_MC_p_UnWeighted");
  h_MC_m_UnWeighted = (TH1D*)fm->Get("SVD_Born.Gen")->Clone("h_MC_m_UnWeighted");
  // To make SVD_Born.Gen to h1_Born_BornFid(Not LumiWeighted)
    h_MC_p_UnWeighted->Scale(1./LumiWeight_Muon_WpToMuNu_S8);
    h_MC_m_UnWeighted->Scale(1./LumiWeight_Muon_WmToMuNu_S8);
  
  h_dataRec_p = (TH1D*)fp->Get("data_Rec")->Clone("h_dataRec_p");
  h_dataRec_m = (TH1D*)fm->Get("data_Rec")->Clone("h_dataRec_m");
  cout << "Inclusive Cross-section" << endl;
  cout << "bin\tW+\t\tW-\t\t Wincl"<<endl;
  for( int ipt(1);ipt<14;ipt++)
  {
    cout<<ipt<<"\t"<<h_data_p->GetBinContent(ipt)<<"\t\t"<<h_data_m->GetBinContent(ipt)<< "\t\t" << h_data_p->GetBinContent(ipt)+h_data_m->GetBinContent(ipt) << endl;
  }
  
  
  
  ///Syst Errors
  //
/*
  //Mu track signal Error
  double mutracksigErrp[14]={0};
  double mutracksigErrm[14]={0};
  double mutracksigErri[14]={0};
  cout <<"W+ tracksigErr \t W- tracksigErr \t W tracksigErr" <<endl;
  for(int i(1);i<14;i++)
  {
    mutracksigErrp[i] = h_data_p->GetBinContent(i)*0.01*mutracksigp[i];
    mutracksigErrm[i] = h_data_m->GetBinContent(i)*0.01*mutracksigm[i];
    mutracksigErri[i] = sqrt(mutracksigErrp[i]*mutracksigErrp[i] + mutracksigErrm[i]*mutracksigErrm[i]);
    cout << mutracksigErrp[i] << "\t" << mutracksigErrm[i] << "\t" << mutracksigErri[i] << endl;
  }
  
  double mutracksigErrpMerge[14]={0};
  double mutracksigErrmMerge[14]={0};
  double mutracksigErriMerge[14]={0};
    mutracksigErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*mutracksigpMerge[4];
    mutracksigErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*mutracksigmMerge[4];
    mutracksigErriMerge[4] = sqrt(mutracksigErrpMerge[4]*mutracksigErrpMerge[4] + mutracksigErrmMerge[4]*mutracksigErrmMerge[4]);
  cout <<"Merged bin:  W+ tracksigErr \t W- tracksigErr \t W tracksigErr" <<endl;
    cout <<"\t\t" <<mutracksigErrpMerge[4] << "\t" << mutracksigErrmMerge[4] << "\t" << mutracksigErriMerge[4] << endl;
 
    // Mu track background Error
  double mutrackbckErrp[14]={0};
  double mutrackbckErrm[14]={0};
  double mutrackbckErri[14]={0};
  cout <<"W+ trackbckErr \t W- trackbckErr \t W trackbckErr" <<endl;
  for(int i(1);i<14;i++)
  {
    mutrackbckErrp[i] = h_data_p->GetBinContent(i)*0.01*mutrackbckp[i];
    mutrackbckErrm[i] = h_data_m->GetBinContent(i)*0.01*mutrackbckm[i];
    mutrackbckErri[i] = sqrt(mutrackbckErrp[i]*mutrackbckErrp[i] + mutrackbckErrm[i]*mutrackbckErrm[i]);
    cout << mutrackbckErrp[i] << "\t" << mutrackbckErrm[i] << "\t" << mutrackbckErri[i] << endl;
  }
  
  double mutrackbckErrpMerge[14]={0};
  double mutrackbckErrmMerge[14]={0};
  double mutrackbckErriMerge[14]={0};
    mutrackbckErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*mutrackbckpMerge[4];
    mutrackbckErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*mutrackbckmMerge[4];
    mutrackbckErriMerge[4] = sqrt(mutrackbckErrpMerge[4]*mutrackbckErrpMerge[4] + mutrackbckErrmMerge[4]*mutrackbckErrmMerge[4]);
  cout <<"Merged bin:  W+ trackbckErr \t W- trackbckErr \t W trackbckErr" <<endl;
    cout <<"\t\t" <<mutrackbckErrpMerge[4] << "\t" << mutrackbckErrmMerge[4] << "\t" << mutrackbckErriMerge[4] << endl;
 
  //Mu idiso signal Error
  double muidisosigErrp[14]={0};
  double muidisosigErrm[14]={0};
  double muidisosigErri[14]={0};
  cout <<"W+ idisosigErr \t W- idisosigErr \t W idisosigErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muidisosigErrp[i] = h_data_p->GetBinContent(i)*0.01*muidisosigp[i];
    muidisosigErrm[i] = h_data_m->GetBinContent(i)*0.01*muidisosigm[i];
    muidisosigErri[i] = sqrt(muidisosigErrp[i]*muidisosigErrp[i] + muidisosigErrm[i]*muidisosigErrm[i]);
    cout << muidisosigErrp[i] << "\t" << muidisosigErrm[i] << "\t" << muidisosigErri[i] << endl;
  }
  
  double muidisosigErrpMerge[14]={0};
  double muidisosigErrmMerge[14]={0};
  double muidisosigErriMerge[14]={0};
    muidisosigErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muidisosigpMerge[4];
    muidisosigErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muidisosigmMerge[4];
    muidisosigErriMerge[4] = sqrt(muidisosigErrpMerge[4]*muidisosigErrpMerge[4] + muidisosigErrmMerge[4]*muidisosigErrmMerge[4]);
  cout <<"Merged bin:  W+ idisosigErr \t W- idisosigErr \t W idisosigErr" <<endl;
    cout <<"\t\t" <<muidisosigErrpMerge[4] << "\t" << muidisosigErrmMerge[4] << "\t" << muidisosigErriMerge[4] << endl;
 
    // Mu idiso background Error
  double muidisobckErrp[14]={0};
  double muidisobckErrm[14]={0};
  double muidisobckErri[14]={0};
  cout <<"W+ idisobckErr \t W- idisobckErr \t W idisobckErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muidisobckErrp[i] = h_data_p->GetBinContent(i)*0.01*muidisobckp[i];
    muidisobckErrm[i] = h_data_m->GetBinContent(i)*0.01*muidisobckm[i];
    muidisobckErri[i] = sqrt(muidisobckErrp[i]*muidisobckErrp[i] + muidisobckErrm[i]*muidisobckErrm[i]);
    cout << muidisobckErrp[i] << "\t" << muidisobckErrm[i] << "\t" << muidisobckErri[i] << endl;
  }
  
  double muidisobckErrpMerge[14]={0};
  double muidisobckErrmMerge[14]={0};
  double muidisobckErriMerge[14]={0};
    muidisobckErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muidisobckpMerge[4];
    muidisobckErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muidisobckmMerge[4];
    muidisobckErriMerge[4] = sqrt(muidisobckErrpMerge[4]*muidisobckErrpMerge[4] + muidisobckErrmMerge[4]*muidisobckErrmMerge[4]);
  cout <<"Merged bin:  W+ idisobckErr \t W- idisobckErr \t W idisobckErr" <<endl;
    cout <<"\t\t" <<muidisobckErrpMerge[4] << "\t" << muidisobckErrmMerge[4] << "\t" << muidisobckErriMerge[4] << endl;
 
  //Mu tracking Error
  double mutrackErrp[14]={0};
  double mutrackErrm[14]={0};
  double mutrackErri[14]={0};
  cout <<"W+ trackErr \t W- trackErr \t W trackErr" <<endl;
  for(int i(1);i<14;i++)
  {
    mutrackErrp[i] = h_data_p->GetBinContent(i)*0.01*mutrackp[i];
    mutrackErrm[i] = h_data_m->GetBinContent(i)*0.01*mutrackm[i];
    mutrackErri[i] = sqrt(mutrackErrp[i]*mutrackErrp[i] + mutrackErrm[i]*mutrackErrm[i]);
    cout << mutrackErrp[i] << "\t" << mutrackErrm[i] << "\t" << mutrackErri[i] << endl;
  }
  
  double mutrackErrpMerge[14]={0};
  double mutrackErrmMerge[14]={0};
  double mutrackErriMerge[14]={0};
    mutrackErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*mutrackpMerge[4];
    mutrackErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*mutrackmMerge[4];
    mutrackErriMerge[4] = sqrt(mutrackErrpMerge[4]*mutrackErrpMerge[4] + mutrackErrmMerge[4]*mutrackErrmMerge[4]);
  cout <<"Merged bin:  W+ trackErr \t W- trackErr \t W trackErr" <<endl;
    cout <<"\t\t" <<mutrackErrpMerge[4] << "\t" << mutrackErrmMerge[4] << "\t" << mutrackErriMerge[4] << endl;
 
    // Mu idiso Error
  double muidisoErrp[14]={0};
  double muidisoErrm[14]={0};
  double muidisoErri[14]={0};
  cout <<"W+ idisoErr \t W- idisoErr \t W idisoErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muidisoErrp[i] = h_data_p->GetBinContent(i)*0.01*muidisop[i];
    muidisoErrm[i] = h_data_m->GetBinContent(i)*0.01*muidisom[i];
    muidisoErri[i] = sqrt(muidisoErrp[i]*muidisoErrp[i] + muidisoErrm[i]*muidisoErrm[i]);
    cout << muidisoErrp[i] << "\t" << muidisoErrm[i] << "\t" << muidisoErri[i] << endl;
  }
  
  double muidisoErrpMerge[14]={0};
  double muidisoErrmMerge[14]={0};
  double muidisoErriMerge[14]={0};
    muidisoErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muidisopMerge[4];
    muidisoErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muidisomMerge[4];
    muidisoErriMerge[4] = sqrt(muidisoErrpMerge[4]*muidisoErrpMerge[4] + muidisoErrmMerge[4]*muidisoErrmMerge[4]);
  cout <<"Merged bin:  W+ idisobckErr \t W- idisobckErr \t W idisobckErr" <<endl;
    cout <<"\t\t" <<muidisobckErrpMerge[4] << "\t" << muidisobckErrmMerge[4] << "\t" << muidisobckErriMerge[4] << endl;
 
    // Mu toy Error
  double mutoyErrp[14]={0};
  double mutoyErrm[14]={0};
  double mutoyErri[14]={0};
  cout <<"W+ toyErr \t W- toyErr \t W toyErr" <<endl;
  for(int i(1);i<14;i++)
  {
    mutoyErrp[i] = h_data_p->GetBinContent(i)*0.01*mutoyp[i];
    mutoyErrm[i] = h_data_m->GetBinContent(i)*0.01*mutoym[i];
    mutoyErri[i] = sqrt(mutoyErrp[i]*mutoyErrp[i] + mutoyErrm[i]*mutoyErrm[i]);
    cout << mutoyErrp[i] << "\t" << mutoyErrm[i] << "\t" << mutoyErri[i] << endl;
  }
  
  double mutoyErrpMerge[14]={0};
  double mutoyErrmMerge[14]={0};
  double mutoyErriMerge[14]={0};
    mutoyErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*mutoypMerge[4];
    mutoyErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*mutoymMerge[4];
    mutoyErriMerge[4] = sqrt(mutoyErrpMerge[4]*mutoyErrpMerge[4] + mutoyErrmMerge[4]*mutoyErrmMerge[4]);
  cout <<"Merged bin:  W+ toyErr \t W- toyErr \t W toyErr" <<endl;
    cout <<"\t\t" <<mutoyErrpMerge[4] << "\t" << mutoyErrmMerge[4] << "\t" << mutoyErriMerge[4] << endl;
 
    // Mu POG Error
  double muPOGErrp[14]={0};
  double muPOGErrm[14]={0};
  double muPOGErri[14]={0};
  cout <<"W+ toyErr \t W- toyErr \t W toyErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muPOGErrp[i] = h_data_p->GetBinContent(i)*0.01*muPOGp[i];
    muPOGErrm[i] = h_data_m->GetBinContent(i)*0.01*muPOGm[i];
    muPOGErri[i] = sqrt(muPOGErrp[i]*muPOGErrp[i] + muPOGErrm[i]*muPOGErrm[i]);
    cout << muPOGErrp[i] << "\t" << muPOGErrm[i] << "\t" << muPOGErri[i] << endl;
  }
  
  double muPOGErrpMerge[14]={0};
  double muPOGErrmMerge[14]={0};
  double muPOGErriMerge[14]={0};
    muPOGErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muPOGpMerge[4];
    muPOGErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muPOGmMerge[4];
    muPOGErriMerge[4] = sqrt(muPOGErrpMerge[4]*muPOGErrpMerge[4] + muPOGErrmMerge[4]*muPOGErrmMerge[4]);
  cout <<"Merged bin:  W+ POGErr \t W- POGErr \t W POGErr" <<endl;
    cout <<"\t\t" <<muPOGErrpMerge[4] << "\t" << muPOGErrmMerge[4] << "\t" << muPOGErriMerge[4] << endl;

    //Mu Total Effi Error
  double muEffErri[14]={0};
  cout <<"W+ EffErr \t W- EffErr \t W EffErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muEffErrp[i] = h_data_p->GetBinContent(i)*0.01*mutotaleffp[i];
    muEffErrm[i] = h_data_m->GetBinContent(i)*0.01*mutotaleffm[i];
    muEffErri[i] = sqrt(muEffErrp[i]*muEffErrp[i] + muEffErrm[i]*muEffErrm[i]);
    cout << muEffErrp[i] << "\t" << muEffErrm[i] << "\t" << muEffErri[i] << endl;
  }

  double muEffErrpMerge[14]={0};
  double muEffErrmMerge[14]={0};
  double muEffErriMerge[14]={0};
    muEffErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*mutotaleffpMerge[4];
    muEffErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*mutotaleffmMerge[4];
    muEffErriMerge[4] = sqrt(muEffErrpMerge[4]*muEffErrpMerge[4] + muEffErrmMerge[4]*muEffErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muEffErrpMerge[4] << "\t" << muEffErrmMerge[4] << "\t" << muEffErriMerge[4] << endl;

// Mu ScaleErr
  double muscaleErrp[14]={0};
  double muscaleErrm[14]={0};
  double muscaleErri[14]={0};
  cout <<"W+ scaleErr \t W- scaleErr \t W EffErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muscaleErrp[i] = h_data_p->GetBinContent(i)*0.01*muscalep[i];
    muscaleErrm[i] = h_data_m->GetBinContent(i)*0.01*muscalem[i];
    muscaleErri[i] = sqrt(muscaleErrp[i]*muscaleErrp[i] + muscaleErrm[i]*muscaleErrm[i]);
    cout << muscaleErrp[i] << "\t" << muscaleErrm[i] << "\t" << muscaleErri[i]<<endl;
  }

  double muscaleErrpMerge[14]={0};
  double muscaleErrmMerge[14]={0};
  double muscaleErriMerge[14]={0};
    muscaleErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muscalepMerge[4];
    muscaleErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muscalemMerge[4];
    muscaleErriMerge[4] = sqrt(muscaleErrpMerge[4]*muscaleErrpMerge[4] + muscaleErrmMerge[4]*muscaleErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muscaleErrpMerge[4] << "\t" << muscaleErrmMerge[4] << "\t" << muscaleErriMerge[4] << endl;

// Mu SmearErr
  double musmearErrp[14]={0};
  double musmearErrm[14]={0};
  double musmearErri[14]={0};
  cout <<"W+ smearErr \t W- smearErr \t W EffErr" <<endl;
  for(int i(1);i<14;i++)
  {
    musmearErrp[i] = h_data_p->GetBinContent(i)*0.01*musmearp[i];
    musmearErrm[i] = h_data_m->GetBinContent(i)*0.01*musmearm[i];
    musmearErri[i] = sqrt(musmearErrp[i]*musmearErrp[i] + musmearErrm[i]*musmearErrm[i]);
    cout << musmearErrp[i] << "\t" << musmearErrm[i] << "\t" << musmearErri[i]<<endl;
  }

  double musmearErrpMerge[14]={0};
  double musmearErrmMerge[14]={0};
  double musmearErriMerge[14]={0};
    musmearErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*musmearpMerge[4];
    musmearErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*musmearmMerge[4];
    musmearErriMerge[4] = sqrt(musmearErrpMerge[4]*musmearErrpMerge[4] + musmearErrmMerge[4]*musmearErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<musmearErrpMerge[4] << "\t" << musmearErrmMerge[4] << "\t" << musmearErriMerge[4] << endl;





// Mu MomResErr
  double muMomResErrp[14]={0};
  double muMomResErrm[14]={0};
  double muMomResErri[14]={0};
  cout <<"W+ MomResErr \t W- MomResErr \t W EffErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muMomResErrp[i] = h_data_p->GetBinContent(i)*0.01*muMomResp[i];
    muMomResErrm[i] = h_data_m->GetBinContent(i)*0.01*muMomResm[i];
    muMomResErri[i] = sqrt(muMomResErrp[i]*muMomResErrp[i] + muMomResErrm[i]*muMomResErrm[i]);
    cout << muMomResErrp[i] << "\t" << muMomResErrm[i] << "\t" << muMomResErri[i]<<endl;
  }

  double muMomResErrpMerge[14]={0};
  double muMomResErrmMerge[14]={0};
  double muMomResErriMerge[14]={0};
    muMomResErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muMomRespMerge[4];
    muMomResErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muMomResmMerge[4];
    muMomResErriMerge[4] = sqrt(muMomResErrpMerge[4]*muMomResErrpMerge[4] + muMomResErrmMerge[4]*muMomResErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muMomResErrpMerge[4] << "\t" << muMomResErrmMerge[4] << "\t" << muMomResErriMerge[4] << endl;

// Mu Met Error
  double muMetErrp[14]={0};
  double muMetErrm[14]={0};
  double muMetErri[14]={0};
  cout <<"W+ MetErr \t W- MetErr \t W MetErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muMetErrp[i] = h_data_p->GetBinContent(i)*0.01*mumetp[i];
    muMetErrm[i] = h_data_m->GetBinContent(i)*0.01*mumetm[i];
    muMetErri[i] = sqrt(muMetErrp[i]*muMetErrp[i] + muMetErrm[i]*muMetErrm[i]);
    cout << muMetErrp[i] << "\t" << muMetErrm[i] <<"\t" << muMetErri[i]<<endl;
  }

  double muMetErrpMerge[14]={0};
  double muMetErrmMerge[14]={0};
  double muMetErriMerge[14]={0};
    muMetErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*mumetpMerge[4];
    muMetErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*mumetmMerge[4];
    muMetErriMerge[4] = sqrt(muMetErrpMerge[4]*muMetErrpMerge[4] + muMetErrmMerge[4]*muMetErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muMetErrpMerge[4] << "\t" << muMetErrmMerge[4] << "\t" << muMetErriMerge[4] << endl;
 

// Mu QCD Background Error
  double muQCDBckErrp[14]={0};
  double muQCDBckErrm[14]={0};
  double muQCDBckErri[14]={0};
  cout <<"W+ QCDBckErr \t W- QCDBckErr \t W QCDBckErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muQCDBckErrp[i] = h_data_p->GetBinContent(i)*0.01*muqcdbckgrp[i];
    muQCDBckErrm[i] = h_data_m->GetBinContent(i)*0.01*muqcdbckgrm[i];
    muQCDBckErri[i] = sqrt(muQCDBckErrp[i]*muQCDBckErrp[i] + muQCDBckErrm[i]*muQCDBckErrm[i]);
    cout << muQCDBckErrp[i] << "\t" << muQCDBckErrm[i] <<"\t" << muQCDBckErri[i] <<endl;
  }
  
  double muQCDBckErrpMerge[14]={0};
  double muQCDBckErrmMerge[14]={0};
  double muQCDBckErriMerge[14]={0};
    muQCDBckErrpMerge[4] =( h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muqcdbckgrpMerge[4];
    muQCDBckErrmMerge[4] =( h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muqcdbckgrmMerge[4];
    muQCDBckErriMerge[4] = sqrt(muQCDBckErrpMerge[4]*muQCDBckErrpMerge[4] + muQCDBckErrmMerge[4]*muQCDBckErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muQCDBckErrpMerge[4] << "\t" << muQCDBckErrmMerge[4] << "\t" << muQCDBckErriMerge[4] << endl; 

// Mu QCD Shape Error 
  double muQCDShapeErrp[14]={0};
  double muQCDShapeErrm[14]={0};
  double muQCDShapeErri[14]={0};
  cout <<"W+ QCDShapeErr \t W- QCDShapeErr \t W QCDShapeErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muQCDShapeErrp[i] = h_data_p->GetBinContent(i)*0.01*muqcdshapep[i];
    muQCDShapeErrm[i] = h_data_m->GetBinContent(i)*0.01*muqcdshapem[i];
    muQCDShapeErri[i] = sqrt(muQCDShapeErrp[i]*muQCDShapeErrp[i] + muQCDShapeErrm[i]*muQCDShapeErrm[i]);
    cout << muQCDShapeErrp[i] << "\t" << muQCDShapeErrm[i] <<"\t" << muQCDShapeErri[i] <<endl;
  }
 
  double muQCDShapeErrpMerge[14]={0};
  double muQCDShapeErrmMerge[14]={0};
  double muQCDShapeErriMerge[14]={0};
    muQCDShapeErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muqcdshapepMerge[4];
    muQCDShapeErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muqcdshapemMerge[4];
    muQCDShapeErriMerge[4] = sqrt(muQCDShapeErrpMerge[4]*muQCDShapeErrpMerge[4] + muQCDShapeErrmMerge[4]*muQCDShapeErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muQCDShapeErrpMerge[4] << "\t" << muQCDShapeErrmMerge[4] << "\t" << muQCDShapeErriMerge[4] << endl;

// Mu EWK Error
  double muEWKErrp[14]={0};
  double muEWKErrm[14]={0};
  double muEWKErri[14]={0};
  cout <<"W+ EWKErr \t W- EWKErr \t W EWKKErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muEWKErrp[i] = h_data_p->GetBinContent(i)*0.01*muewkp[i];
    muEWKErrm[i] = h_data_m->GetBinContent(i)*0.01*muewkm[i];
    //muEWKErri[i] = sqrt(muEWKErrp[i]*muEWKErrp[i] + muEWKErrm[i]*muEWKErrm[i]);
    muEWKErri[i] = muEWKErrp[i] + muEWKErrm[i];
    cout << muEWKErrp[i] << "\t" << muEWKErrm[i] <<"\t" << muEWKErri[i] <<endl;
  }
 
  double muEWKErrpMerge[14]={0};
  double muEWKErrmMerge[14]={0};
  double muEWKErriMerge[14]={0};
    muEWKErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*muewkpMerge[4];
    muEWKErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*muewkmMerge[4];
    //muEWKErriMerge[4] = sqrt(muEWKErrpMerge[4]*muEWKErrpMerge[4] + muEWKErrmMerge[4]*muEWKErrmMerge[4]);
    muEWKErriMerge[4] = muEWKErrpMerge[4] + muEWKErrmMerge[4];
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muEWKErrpMerge[4] << "\t" << muEWKErrmMerge[4] << "\t" << muEWKErriMerge[4] << endl;
*/

    // Mu SVDUnf Error
  double muSVDUnfErrp[14]={0};
  double muSVDUnfErrm[14]={0};
  double muSVDUnfErri[14]={0};
  cout <<"W+ SVDUnfErr \t W- SVDUnfErr \t W SVDUnfErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muSVDUnfErrp[i] = h_data_p->GetBinContent(i)*0.01*musvdunfp[i];
    muSVDUnfErrm[i] = h_data_m->GetBinContent(i)*0.01*musvdunfm[i];
    //muSVDUnfErri[i] = sqrt(muSVDUnfErrp[i]*muSVDUnfErrp[i] + muSVDUnfErrm[i]*muSVDUnfErrm[i]);
    muSVDUnfErri[i] = (sqrt(muSVDUnfErrp[i]*muSVDUnfErrp[i] + muSVDUnfErrm[i]*muSVDUnfErrm[i])/(h_data_p->GetBinContent(i)+h_data_m->GetBinContent(i))) * 100;
    cout << muSVDUnfErrp[i] << "\t" << muSVDUnfErrm[i] <<"\t" << muSVDUnfErri[i] <<endl;
  }
 
  double muSVDUnfErrpMerge[14]={0};
  double muSVDUnfErrmMerge[14]={0};
  double muSVDUnfErriMerge[14]={0};
    muSVDUnfErrpMerge[4] = ( h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*musvdunfpMerge[4];
    muSVDUnfErrmMerge[4] = ( h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*musvdunfmMerge[4];
    //muSVDUnfErriMerge[4] = sqrt(muSVDUnfErrpMerge[4]*muSVDUnfErrpMerge[4] + muSVDUnfErrmMerge[4]*muSVDUnfErrmMerge[4]);
    muSVDUnfErriMerge[4] = (sqrt(muSVDUnfErrpMerge[4]*muSVDUnfErrpMerge[4] + muSVDUnfErrmMerge[4]*muSVDUnfErrmMerge[4])/(h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5) + h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5)))*100;
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muSVDUnfErrpMerge[4] << "\t" << muSVDUnfErrmMerge[4] << "\t" << muSVDUnfErriMerge[4] << endl;
/* 
    // Mu FSR ERROR
  double muFSRErrp[14]={0};
  double muFSRErrm[14]={0};
  double muFSRErri[14]={0};
  cout <<"W+ FSRErr \t W- FSRErr \t W FSRErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muFSRErrp[i] = h_data_p->GetBinContent(i)*0.01*mufsrp[i];
    muFSRErrm[i] = h_data_m->GetBinContent(i)*0.01*mufsrm[i];
    //muFSRErri[i] = sqrt(muFSRErrp[i]*muFSRErrp[i] + muFSRErrm[i]*muFSRErrm[i]);
    muFSRErri[i] = muFSRErrp[i] + muFSRErrm[i];
    cout << muFSRErrp[i] << "\t" << muFSRErrm[i] <<"\t" << muFSRErri[i] <<endl;
  }

  double muFSRErrpMerge[14]={0};
  double muFSRErrmMerge[14]={0};
  double muFSRErriMerge[14]={0};
    muFSRErrpMerge[4] = ( h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*mufsrpMerge[4];
    muFSRErrmMerge[4] = ( h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*mufsrmMerge[4];
    //muFSRErriMerge[4] = sqrt(muFSRErrpMerge[4]*muFSRErrpMerge[4] + muFSRErrmMerge[4]*muFSRErrmMerge[4]);
    muFSRErriMerge[4] = muFSRErrpMerge[4] + muFSRErrmMerge[4];
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muFSRErrpMerge[4] << "\t" << muFSRErrmMerge[4] << "\t" << muFSRErriMerge[4] << endl;

 // Mu LUMI Error 
  double muLumiErrp[14]={0};
  double muLumiErrm[14]={0};
  double muLumiErri[14]={0};
  cout <<"W+ LumiErr \t W- LumiErr \t W LumiErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muLumiErrp[i] = h_data_p->GetBinContent(i)*0.01*2.6;
    muLumiErrm[i] = h_data_m->GetBinContent(i)*0.01*2.6;
    //muLumiErri[i] = sqrt(muLumiErrp[i]*muLumiErrp[i] + muLumiErrm[i]*muLumiErrm[i]);
    muLumiErri[i] = muLumiErrp[i]+muLumiErrm[i];
    cout << muLumiErrp[i] << "\t" << muLumiErrm[i] <<"\t" << muLumiErri[i] <<endl;
  }
 
  double muLumiErrpMerge[14]={0};
  double muLumiErrmMerge[14]={0};
  double muLumiErriMerge[14]={0};
    muLumiErrpMerge[4] = (h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5))*0.01*2.6;
    muLumiErrmMerge[4] = (h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5))*0.01*2.6;
    //muLumiErriMerge[4] = sqrt(muLumiErrpMerge[4]*muLumiErrpMerge[4] + muLumiErrmMerge[4]*muLumiErrmMerge[4]);
    muLumiErriMerge[4] = muLumiErrpMerge[4]+muLumiErrmMerge[4];
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muLumiErrpMerge[4] << "\t" << muLumiErrmMerge[4] << "\t" << muLumiErriMerge[4] << endl;

    // Mu UnfBias Error
  double muUnfBiasErrp[14]={0};
  double muUnfBiasErrm[14]={0};
  double muUnfBiasErri[14]={0};
  cout <<"W+ UnfBiasErr \t W- UnfBiasErr \t W UnfBiasErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muUnfBiasErrp[i] = h_data_p->GetBinContent(i)*0.01*muBiasp[i];
    muUnfBiasErrm[i] = h_data_m->GetBinContent(i)*0.01*muBiasm[i];
    muUnfBiasErri[i] = sqrt(muUnfBiasErrp[i]*muUnfBiasErrp[i] + muUnfBiasErrm[i]*muUnfBiasErrm[i]);
    cout << muUnfBiasErrp[i] << "\t" << muUnfBiasErrm[i] <<"\t" << muUnfBiasErri[i] <<endl;
  }
 
 double muUnfBiasErrpMerge[14]={0};
  double muUnfBiasErrmMerge[14]={0};
  double muUnfBiasErriMerge[14]={0};
    muUnfBiasErrpMerge[4] = ( h_data_p->GetBinContent(4)+h_data_p->GetBinContent(5) )*0.01*muBiaspMerge[4];
    muUnfBiasErrmMerge[4] = ( h_data_m->GetBinContent(4)+h_data_m->GetBinContent(5) )*0.01*muBiasmMerge[4];
    muUnfBiasErriMerge[4] = sqrt(muUnfBiasErrpMerge[4]*muUnfBiasErrpMerge[4] + muUnfBiasErrmMerge[4]*muUnfBiasErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muUnfBiasErrpMerge[4] << "\t" << muUnfBiasErrmMerge[4] << "\t" << muUnfBiasErriMerge[4] << endl;
  

    // Mu Stat Error
  double muStatErrp[14]={0};
  double muStatErrm[14]={0};
  double muStatErri[14]={0};
  cout <<"W+ StatErr \t W- StatErr \t W StatErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muStatErrp[i] = h_data_p->GetBinContent(i)*0.01*mustatp[i];
    muStatErrm[i] = h_data_m->GetBinContent(i)*0.01*mustatm[i];
    muStatErri[i] = sqrt(muStatErrp[i]*muStatErrp[i] + muStatErrm[i]*muStatErrm[i]);
    cout << muStatErrp[i] << "\t" << muStatErrm[i] <<"\t" << muStatErri[i] <<endl;
  }
 
  double muStatErrpMerge[14]={0};
  double muStatErrmMerge[14]={0};
  double muStatErriMerge[14]={0};
    muStatErrpMerge[4] = h_data_p->GetBinContent(4)*0.01*mustatpMerge[4];
    muStatErrmMerge[4] = h_data_m->GetBinContent(4)*0.01*mustatmMerge[4];
    muStatErriMerge[4] = sqrt(muStatErrpMerge[4]*muStatErrpMerge[4] + muStatErrmMerge[4]*muStatErrmMerge[4]);
  cout <<"Merged bin:  W+ EffErr \t W- EffErr \t W EffErr" <<endl;
    cout <<"\t\t" <<muStatErrpMerge[4] << "\t" << muStatErrmMerge[4] << "\t" << muStatErriMerge[4] << endl;
 */ 

 /* 
// Mu Total Syst
  double muTotalSysti[14]={0};
  cout <<"W+ TotSyst \t W- TotSyst \t W TotSyst" <<endl;
  for(int i(1);i<14;i++)
  {
    //muTotalSysti[i] = sqrt(muEffErri[i]*muEffErri[i] + muMomResErri[i]*muMomResErri[i]+muMetErri[i]*muMetErri[i]+ muQCDBckErri[i]*muQCDBckErri[i] + muQCDShapeErri[i]*muQCDShapeErri[i]+ muEWKErri[i]*muEWKErri[i] + muSVDUnfErri[i]*muSVDUnfErri[i] + muFSRErri[i]*muFSRErri[i]+ muLumiErri[i]*muLumiErri[i] + muUnfBiasErri[i]*muUnfBiasErri[i]);
    muTotalSysti[i] = sqrt(muEffErri[i]*muEffErri[i] + muMomResErri[i]*muMomResErri[i]+muMetErri[i]*muMetErri[i]+ muQCDBckErri[i]*muQCDBckErri[i] + muQCDShapeErri[i]*muQCDShapeErri[i]+ muEWKErri[i]*muEWKErri[i] + muSVDUnfErri[i]*muSVDUnfErri[i] + muFSRErri[i]*muFSRErri[i]+ muUnfBiasErri[i]*muUnfBiasErri[i]);
    cout << muTotalSysti[i] <<endl;
  }
    
  double muTotalSystiMerge[14]={0};
  cout <<"Merged bin4:   W+ TotSyst \t W- TotSyst \t W TotSyst" <<endl;
    //muTotalSystiMerge[4] = sqrt(muEffErriMerge[4]*muEffErriMerge[4] + muMomResErriMerge[4]*muMomResErriMerge[4]+muMetErriMerge[4]*muMetErriMerge[4]+ muQCDBckErriMerge[4]*muQCDBckErriMerge[4] + muQCDShapeErriMerge[4]*muQCDShapeErriMerge[4]+ muEWKErriMerge[4]*muEWKErriMerge[4] + muSVDUnfErriMerge[4]*muSVDUnfErriMerge[4] + muFSRErriMerge[4]*muFSRErriMerge[4]+ muLumiErriMerge[4]*muLumiErriMerge[4] + muUnfBiasErriMerge[4]*muUnfBiasErriMerge[4]);
    muTotalSystiMerge[4] = sqrt(muEffErriMerge[4]*muEffErriMerge[4] + muMomResErriMerge[4]*muMomResErriMerge[4]+muMetErriMerge[4]*muMetErriMerge[4]+ muQCDBckErriMerge[4]*muQCDBckErriMerge[4] + muQCDShapeErriMerge[4]*muQCDShapeErriMerge[4]+ muEWKErriMerge[4]*muEWKErriMerge[4] + muSVDUnfErriMerge[4]*muSVDUnfErriMerge[4] + muFSRErriMerge[4]*muFSRErriMerge[4]+ muUnfBiasErriMerge[4]*muUnfBiasErriMerge[4]);

    cout  <<  muTotalSystiMerge[4] <<endl; 

     
 
 
  
  
 // Mu Total Uncertainty 
  double muTotalUnceri[14]={0};
  cout <<" Wincl TotalUncertainty" <<endl;
  for(int i(1);i<14;i++)
  {
    //muTotalUnceri[i] = sqrt(muEffErri[i]*muEffErri[i] + muMomResErri[i]*muMomResErri[i]+muMetErri[i]*muMetErri[i]+ muQCDBckErri[i]*muQCDBckErri[i] + muQCDShapeErri[i]*muQCDShapeErri[i]+ muEWKErri[i]*muEWKErri[i] + muSVDUnfErri[i]*muSVDUnfErri[i] + muFSRErri[i]*muFSRErri[i]+ muLumiErri[i]*muLumiErri[i] + muUnfBiasErri[i]*muUnfBiasErri[i]+ muStatErri[i]*muStatErri[i]);
    muTotalUnceri[i] = sqrt(muEffErri[i]*muEffErri[i] + muMomResErri[i]*muMomResErri[i]+muMetErri[i]*muMetErri[i]+ muQCDBckErri[i]*muQCDBckErri[i] + muQCDShapeErri[i]*muQCDShapeErri[i]+ muEWKErri[i]*muEWKErri[i] + muSVDUnfErri[i]*muSVDUnfErri[i] + muFSRErri[i]*muFSRErri[i]+muUnfBiasErri[i]*muUnfBiasErri[i]+ muStatErri[i]*muStatErri[i]);
    cout << muTotalUnceri[i] <<endl;
  }



  double muTotalUnceriMerge[14]={0};
  cout <<"Merged bin4:    Wincl TotalUncertainty" <<endl;
    //muTotalUnceriMerge[4] = sqrt(muEffErriMerge[4]*muEffErriMerge[4] + muMomResErriMerge[4]*muMomResErriMerge[4]+muMetErriMerge[4]*muMetErriMerge[4]+ muQCDBckErriMerge[4]*muQCDBckErriMerge[4] + muQCDShapeErriMerge[4]*muQCDShapeErriMerge[4]+ muEWKErriMerge[4]*muEWKErriMerge[4] + muSVDUnfErriMerge[4]*muSVDUnfErriMerge[4] + muFSRErriMerge[4]*muFSRErriMerge[4]+ muLumiErriMerge[4]*muLumiErriMerge[4] + muUnfBiasErriMerge[4]*muUnfBiasErriMerge[4]+ muStatErriMerge[4]*muStatErriMerge[4]);
    muTotalUnceriMerge[4] = sqrt(muEffErriMerge[4]*muEffErriMerge[4] + muMomResErriMerge[4]*muMomResErriMerge[4]+muMetErriMerge[4]*muMetErriMerge[4]+ muQCDBckErriMerge[4]*muQCDBckErriMerge[4] + muQCDShapeErriMerge[4]*muQCDShapeErriMerge[4]+ muEWKErriMerge[4]*muEWKErriMerge[4] + muSVDUnfErriMerge[4]*muSVDUnfErriMerge[4] + muFSRErriMerge[4]*muFSRErriMerge[4]+ muUnfBiasErriMerge[4]*muUnfBiasErriMerge[4]+ muStatErriMerge[4]*muStatErriMerge[4]);
    
    cout <<"\t\t" << muTotalUnceriMerge[4] <<endl;
*/

  TString resultDir = "ResultWInclToMuNu";
  gSystem->mkdir(resultDir,kTRUE);

  TFile f_out(resultDir+"/Result_WInclToMuNu_12Bin.root","recreate");
 
  TH1D* BornEffCorr13 = new TH1D("BornEffCorr13","BornEffCorr13",13,0,13);
  TH1D* SVD_BornGen13 = new TH1D("SVD_Born.Gen13","SVD_Born.Gen13",13,0,13);SVD_BornGen13->Sumw2();
  TH1D* SVD_BornGen_UnWeighted13 = new TH1D("SVD_Born.Gen_UnWeighted13","SVD_Born.Gen_UnWeighted13",13,0,13);SVD_BornGen_UnWeighted13->Sumw2();
  TH1D* PowhegErr13 = new TH1D("PowhegErr13","PowhegErr13",13,0,13);
  TH1D* data_Rec13 = new TH1D("data_Rec13","data_Rec13",13,0,13);data_Rec13->Sumw2();
  
  TH1D* BornEffCorr12 = new TH1D("BornEffCorr","BornEffCorr",12,0,12);
  TH1D* SVD_BornGen12 = new TH1D("SVD_Born.Gen","SVD_Born.Gen",12,0,12);SVD_BornGen12->Sumw2();
  TH1D* hPowheg_Yield_BornAfterFidCut12 = new TH1D("hPowheg_Yield_BornAfterFidCut","hPowheg_Yield_BornAfterFidCut",12,0,12);hPowheg_Yield_BornAfterFidCut12->Sumw2();
  TH1D* data_Rec12 = new TH1D("data_Rec","data_Rec",12,0,12);data_Rec12->Sumw2();


  TH1D* h_MC_p12 = new TH1D("h_MC_p12","h_MC_p12",12,0,12);h_MC_p12->Sumw2();
  TH1D* h_MC_m12 = new TH1D("h_MC_m12","h_MC_m12",12,0,12);h_MC_m12->Sumw2();

  h_MC_p->Scale(1./18.429);
  h_MC_m->Scale(1./18.429);
  
  double PowhegPDFp[14];
    PowhegPDFp[0] = 0; 
    PowhegPDFp[1] =3.0116; //4.268; 
    PowhegPDFp[2] =3.7087; //4.147; 
    PowhegPDFp[3] =4.0534; //4.122; 
    PowhegPDFp[4] =4.1509; //4.123; 
    PowhegPDFp[5] =4.3433; //4.132; 
    PowhegPDFp[6] =4.3076; //4.126; 
    PowhegPDFp[7] =4.4361; //4.143; 
    PowhegPDFp[8] =5.4356; //5.279; 
    PowhegPDFp[9] =4.5321; //4.222; 
    PowhegPDFp[10]=4.8190; //4.426;
    PowhegPDFp[11]=5.1662; //4.819; 
    PowhegPDFp[12]=5.3727; //5.075; 
    PowhegPDFp[13]=6.3736; //6.084; 

    double PowhegPDFm[14];
    PowhegPDFm[0] = 0; 
    PowhegPDFm[1] =3.1886;// 4.4  ; 
    PowhegPDFm[2] =3.9034;// 4.385; 
    PowhegPDFm[3] =4.3073;// 4.426;
    PowhegPDFm[4] =4.4302;// 4.444;
    PowhegPDFm[5] =4.6290;// 4.451;
    PowhegPDFm[6] =4.6619;// 4.532;
    PowhegPDFm[7] =7.7632;// 7.792;
    PowhegPDFm[8] =6.2674;// 6.181;
    PowhegPDFm[9] =4.9092;// 4.622;
    PowhegPDFm[10]=4.7558;// 4.33 ;
    PowhegPDFm[11]=4.9631;// 4.573;
    PowhegPDFm[12]=5.3300;// 5.005;
    PowhegPDFm[13]=5.6229;// 5.269;
    
  double PowhegPDFpMerge[14]={0};
  double PowhegPDFmMerge[14]={0};
  //PowhegPDFpMerge[4]=TMath::Max(PowhegPDFp[4],PowhegPDFp[5]);
  //PowhegPDFmMerge[4]=TMath::Max(PowhegPDFm[4],PowhegPDFm[5]);
  PowhegPDFpMerge[4]=(BinWidth[4]*PowhegPDFp[4]*(h_MC_p->GetBinContent(4)/BinWidth[4]) + BinWidth[5]*PowhegPDFp[5]*(h_MC_p->GetBinContent(5)/BinWidth[5]) )
    		/ (BinWidth[4]*(h_MC_p->GetBinContent(4)/BinWidth[4]) + BinWidth[5]*(h_MC_p->GetBinContent(5)/BinWidth[5]));
  PowhegPDFmMerge[4]=(BinWidth[4]*PowhegPDFm[4]*(h_MC_m->GetBinContent(4)/BinWidth[4]) + BinWidth[5]*PowhegPDFm[5]*(h_MC_m->GetBinContent(5)/BinWidth[5]) )
    		/ (BinWidth[4]*(h_MC_m->GetBinContent(4)/BinWidth[4]) + BinWidth[5]*(h_MC_m->GetBinContent(5)/BinWidth[5]) );
  
  //Incl W Powheg stat error
    
  
  double PowhegPDFErrp[14]={0};
  double PowhegPDFErrm[14]={0};
  double PowhegPDFErri[14]={0};
  cout <<"W+ PowhegPDFErr \t W- PowhegPDFErr \t W PowhegPDFErr" <<endl;
  for(int i(1);i<14;i++)
  {
    PowhegPDFErrp[i] = h_MC_p->GetBinContent(i)*0.01*PowhegPDFp[i];
    PowhegPDFErrm[i] = h_MC_m->GetBinContent(i)*0.01*PowhegPDFm[i];
    //PowhegPDFErri[i] = sqrt(PowhegPDFErrp[i]*PowhegPDFErrp[i] + PowhegPDFErrm[i]*PowhegPDFErrm[i]);
    PowhegPDFErri[i] = PowhegPDFErrp[i] + PowhegPDFErrm[i];
    cout << PowhegPDFErrp[i] << "\t" << PowhegPDFErrm[i] <<"\t" << PowhegPDFErri[i] <<endl;
  }
  double PowhegPDFErrpMerge[14]={0};
  double PowhegPDFErrmMerge[14]={0};
  double PowhegPDFErriMerge[14]={0};
    PowhegPDFErrpMerge[4] = (h_MC_p->GetBinContent(4)+h_MC_p->GetBinContent(5))*0.01*PowhegPDFpMerge[4];
    PowhegPDFErrmMerge[4] = (h_MC_m->GetBinContent(4)+h_MC_m->GetBinContent(5))*0.01*PowhegPDFpMerge[4];
    //PowhegPDFErriMerge[4] = sqrt(PowhegPDFErrpMerge[4]*PowhegPDFErrpMerge[4] + PowhegPDFErrmMerge[4]*PowhegPDFErrmMerge[4]);
    PowhegPDFErriMerge[4] = PowhegPDFErrpMerge[4] + PowhegPDFErrmMerge[4];
  cout <<"Merged bin:  W+ StatErr \t W- StatErr \t W StatErr" <<endl;
    cout <<"\t\t" <<PowhegPDFErrpMerge[4] << "\t" << PowhegPDFErrmMerge[4] << "\t" << PowhegPDFErriMerge[4] << endl;


    //Powheg stat error
  TH1D* PowhegStatErrRatioP = new TH1D("PowhegStatErrRatioP","PowhegStatErrRatioP",12,0,12);
  TH1D* PowhegStatErrRatioM = new TH1D("PowhegStatErrRatioM","PowhegStatErrRatioM",12,0,12);
   double PowhegStatErrP[13];
   double PowhegStatErrM[13];
   double PowhegStatErrI[13];
   for(int i(1);i<13;i++)
   {
     if (i<4)
      {
        PowhegStatErrRatioP->SetBinContent(i,sqrt(h_MC_p_UnWeighted->GetBinContent(i))/h_MC_p_UnWeighted->GetBinContent(i));
        PowhegStatErrRatioM->SetBinContent(i,sqrt(h_MC_m_UnWeighted->GetBinContent(i))/h_MC_m_UnWeighted->GetBinContent(i));
      
	PowhegStatErrP[i] = h_MC_p->GetBinContent(i) * PowhegStatErrRatioP->GetBinContent(i);
        PowhegStatErrM[i] = h_MC_m->GetBinContent(i) * PowhegStatErrRatioM->GetBinContent(i);
        
	PowhegStatErrI[i] = sqrt(PowhegStatErrP[i]*PowhegStatErrP[i] + PowhegStatErrM[i]*PowhegStatErrM[i]);
      }else if (i==4)
      {
        PowhegStatErrRatioP->SetBinContent(i,sqrt(h_MC_p_UnWeighted->GetBinContent(i)+h_MC_p_UnWeighted->GetBinContent(i+1))/(h_MC_p_UnWeighted->GetBinContent(i)+h_MC_p_UnWeighted->GetBinContent(i)));
        PowhegStatErrRatioM->SetBinContent(i,sqrt(h_MC_m_UnWeighted->GetBinContent(i)+h_MC_m_UnWeighted->GetBinContent(i+1))/(h_MC_m_UnWeighted->GetBinContent(i)+h_MC_m_UnWeighted->GetBinContent(i)));
	
	PowhegStatErrP[i] = (h_MC_p->GetBinContent(i)+h_MC_p->GetBinContent(i+1)) * PowhegStatErrRatioP->GetBinContent(i);
        PowhegStatErrM[i] = (h_MC_m->GetBinContent(i)+h_MC_m->GetBinContent(i+1)) * PowhegStatErrRatioM->GetBinContent(i);
	
	PowhegStatErrI[i] = sqrt(PowhegStatErrP[i]*PowhegStatErrP[i] + PowhegStatErrM[i]*PowhegStatErrM[i]);
      }else if (i>=5)
      {
        PowhegStatErrRatioP->SetBinContent(i,sqrt(h_MC_p_UnWeighted->GetBinContent(i+1))/h_MC_p_UnWeighted->GetBinContent(i+1));
        PowhegStatErrRatioM->SetBinContent(i,sqrt(h_MC_m_UnWeighted->GetBinContent(i+1))/h_MC_m_UnWeighted->GetBinContent(i+1));
       
	PowhegStatErrP[i] = h_MC_p->GetBinContent(i+1) * PowhegStatErrRatioP->GetBinContent(i);
        PowhegStatErrM[i] = h_MC_m->GetBinContent(i+1) * PowhegStatErrRatioM->GetBinContent(i);
	

       PowhegStatErrI[i] = sqrt(PowhegStatErrP[i]*PowhegStatErrP[i] + PowhegStatErrM[i]*PowhegStatErrM[i]);
      }
       
   }


  ///12 bin Inclusive 
  for(int i(1);i<13;i++)
    {
      if (i<4)
      {
	BornEffCorr12->SetBinContent(i,h_data_p->GetBinContent(i) + h_data_m->GetBinContent(i));
        SVD_BornGen12->SetBinContent(i,h_MC_p->GetBinContent(i) + h_MC_m->GetBinContent(i));
        hPowheg_Yield_BornAfterFidCut12->SetBinContent(i,h_MC_p_UnWeighted->GetBinContent(i) + h_MC_m_UnWeighted->GetBinContent(i));
        data_Rec12->SetBinContent(i,h_dataRec_p->GetBinContent(i) + h_dataRec_m->GetBinContent(i));
        
      }else if (i==4)
      {
	BornEffCorr12->SetBinContent(i, (h_data_p->GetBinContent(i)+ h_data_p->GetBinContent(i+1) ) + (h_data_m->GetBinContent(i)+h_data_m->GetBinContent(i+1)) );
        SVD_BornGen12->SetBinContent(i, (h_MC_p->GetBinContent(i)+h_MC_p->GetBinContent(i+1)) + (h_MC_m->GetBinContent(i)+h_MC_m->GetBinContent(i+1)) );
        hPowheg_Yield_BornAfterFidCut12->SetBinContent(i, (h_MC_p_UnWeighted->GetBinContent(i)+h_MC_p_UnWeighted->GetBinContent(i+1)) + ( h_MC_m_UnWeighted->GetBinContent(i)+h_MC_m_UnWeighted->GetBinContent(i+1))  );
        
        data_Rec12->SetBinContent(i, (h_dataRec_p->GetBinContent(i) +h_dataRec_p->GetBinContent(i+1)) + (h_dataRec_m->GetBinContent(i)+h_dataRec_m->GetBinContent(i+1) ) );
        
      }else if (i>=5)
      {
	BornEffCorr12->SetBinContent(i,h_data_p->GetBinContent(i+1) + h_data_m->GetBinContent(i+1));
        SVD_BornGen12->SetBinContent(i,h_MC_p->GetBinContent(i+1) + h_MC_m->GetBinContent(i+1));
        hPowheg_Yield_BornAfterFidCut12->SetBinContent(i,h_MC_p_UnWeighted->GetBinContent(i+1) + h_MC_m_UnWeighted->GetBinContent(i+1));
        data_Rec12->SetBinContent(i,h_dataRec_p->GetBinContent(i+1) + h_dataRec_m->GetBinContent(i+1));

      }    
    }
  
  // make Born level to diff Xsec 13bin
  TH1D* Diff_Xsec13 = new TH1D("Diff_Xsec13","Diff_Xsec13",13,0,13);
  TH1D* Diff_Xsec13_Wp = new TH1D("Diff_Xsec13_Wp","Diff_Xsec13_Wp",13,0,13);
  TH1D* Diff_Xsec13_Wm = new TH1D("Diff_Xsec13_Wm","Diff_Xsec13_Wm",13,0,13);

  for(int i(1); i<14; i++)
  {
    Diff_Xsec13->SetBinContent(i,h_data_p->GetBinContent(i) + h_data_m->GetBinContent(i));
    Diff_Xsec13_Wp->SetBinContent(i,h_data_p->GetBinContent(i));
    Diff_Xsec13_Wm->SetBinContent(i,h_data_m->GetBinContent(i));
  }
  Diff_Xsec13->Scale(1.0/18.429);
  Diff_Xsec13_Wp->Scale(1.0/18.429);
  Diff_Xsec13_Wm->Scale(1.0/18.429);
  for(int i(1); i<14; i++)
  {
    cout <<"Diff_Xsec : " <<  Diff_Xsec13->GetBinContent(i) / BinWidth[i] << endl; 
  }
  for(int i(1); i<14; i++)
  {
    cout <<"Diff_Xsec_Wp : " <<  Diff_Xsec13_Wp->GetBinContent(i) / BinWidth[i] << endl; 
  }
  for(int i(1); i<14; i++)
  {
    cout <<"Diff_Xsec_Wm : " <<  Diff_Xsec13_Wm->GetBinContent(i) / BinWidth[i] << endl; 
  }

    BornEffCorr12->Write();
    SVD_BornGen12->Write();
    hPowheg_Yield_BornAfterFidCut12->Write();
    data_Rec12->Write();
    ///Write all errors to root
    
    TH1D* h_tracksig = new TH1D("h_tracksig","h_tracksig",12,0,12);
    TH1D* h_trackbck = new TH1D("h_trackbck","h_trackbck",12,0,12);
    
    TH1D* h_idisosig = new TH1D("h_idisosig","h_idisosig",12,0,12);
    TH1D* h_idisobck = new TH1D("h_idisobck","h_idisobck",12,0,12);
  
    TH1D* h_track = new TH1D("h_track","h_track",12,0,12);
   
    TH1D* h_idiso = new TH1D("h_idiso","h_idiso",12,0,12);
    
    TH1D* h_toy = new TH1D("h_toy","h_toy",12,0,12);
    
    TH1D* h_POG = new TH1D("h_POG","h_POG",12,0,12);
    
    TH1D* h_TotalEff = new TH1D("h_Totaleff","h_TotalEff",12,0,12);
    
    TH1D* h_met = new TH1D("h_met","h_met",12,0,12);
    
    TH1D* h_scale = new TH1D("h_scale","h_scale",12,0,12);
    
    TH1D* h_smear = new TH1D("h_smear","h_smear",12,0,12);
    
    TH1D* h_EnMomRes = new TH1D("h_EnMomRes","h_EnMomRes",12,0,12);
    
    TH1D* h_qcdbckgr = new TH1D("h_qcdbckgr","h_qcdbckgr",12,0,12);
    
    TH1D* h_qcdshape = new TH1D("h_qcdshape","h_qcdshape",12,0,12);
    
    TH1D* h_ewk = new TH1D("h_ewk","h_ewk",12,0,12);
    
    TH1D* h_fsr = new TH1D("h_fsr","h_fsr",12,0,12);
    
    TH1D* h_SvdUnf = new TH1D("h_SvdUnf","h_SvdUnf",12,0,12);
    
    TH1D* h_UnfoldBias = new TH1D("h_UnfoldBias","h_UnfoldBias",12,0,12);
    
    TH1D* h_TotalSyst = new TH1D("h_TotalSyst","h_TotalSyst",12,0,12);
    
    TH1D* h_Stat = new TH1D("h_Stat","h_Stat",12,0,12);
    
    TH1D* h_TotalUncer = new TH1D("h_TotalUncer","h_TotalUncer",12,0,12);
    
    TH1D* h_PowhegPDF = new TH1D("h_PowhegPDF","h_PowhegPDF",12,0,12);
    
    TH1D* h_PowhegStat = new TH1D("h_PowhegStat","h_PowhegStat",12,0,12);

    TH1D* h_LumiSyst = new TH1D("h_LumiSyst","h_LumiSyst",12,0,12);


    for(int i(1);i<13;i++)
    {
	h_Stat->SetBinContent(i,muStatErri[i]);
       if(i<4){
	h_tracksig->SetBinContent(i,mutracksigErri[i]);
	h_trackbck->SetBinContent(i,mutrackbckErri[i]);
        h_idisosig->SetBinContent(i,muidisosigErri[i]);
        h_idisobck->SetBinContent(i,muidisobckErri[i]);
	//h_track->SetBinContent(i,mutrackErri[i]);
        //h_idiso->SetBinContent(i,muidisoErri[i]);
        h_toy->SetBinContent(i,mutoyErri[i]);
        h_POG->SetBinContent(i,muPOGErri[i]);
        //h_TotalEff->SetBinContent(i,muEffErri[i]);
        h_met->SetBinContent(i,muMetErri[i]);
        h_scale->SetBinContent(i,muscaleErri[i]);
        h_smear->SetBinContent(i,musmearErri[i]);
        //h_EnMomRes->SetBinContent(i,muMomResErri[i]);
        h_qcdbckgr->SetBinContent(i,muQCDBckErri[i]);
        h_qcdshape->SetBinContent(i,muQCDShapeErri[i]);
        h_ewk->SetBinContent(i,muEWKErri[i]);
        h_fsr->SetBinContent(i,muFSRErri[i]);
        h_SvdUnf->SetBinContent(i,muSVDUnfErri[i]);
        h_UnfoldBias->SetBinContent(i,muUnfBiasErri[i]);
        //h_TotalSyst->SetBinContent(i,muTotalSysti[i]);
        //h_LumiSyst->SetBinContent(i,muLumiErri[i]);
        
        h_PowhegPDF->SetBinContent(i,PowhegPDFErri[i]);
	h_PowhegStat->SetBinContent(i,PowhegStatErrI[i]);
      }else  if(i==4) {
	h_tracksig->SetBinContent(i,mutracksigErriMerge[i]);
	h_trackbck->SetBinContent(i,mutrackbckErriMerge[i]);
        h_idisosig->SetBinContent(i,muidisosigErriMerge[i]);
        h_idisobck->SetBinContent(i,muidisobckErriMerge[i]);
	//h_track->SetBinContent(i,mutrackErriMerge[i]);
        //h_idiso->SetBinContent(i,muidisoErriMerge[i]);
        h_toy->SetBinContent(i,mutoyErriMerge[i]);
        h_POG->SetBinContent(i,muPOGErriMerge[i]);
        //h_TotalEff->SetBinContent(i,muEffErriMerge[i]);
        h_met->SetBinContent(i,muMetErriMerge[i]);
        h_scale->SetBinContent(i,muscaleErriMerge[i]);
        h_smear->SetBinContent(i,musmearErriMerge[i]);
        //h_EnMomRes->SetBinContent(i,muMomResErriMerge[i]);
        h_qcdbckgr->SetBinContent(i,muQCDBckErriMerge[i]);
        h_qcdshape->SetBinContent(i,muQCDShapeErriMerge[i]);
        h_ewk->SetBinContent(i,muEWKErriMerge[i]);
        h_fsr->SetBinContent(i,muFSRErriMerge[i]);
        h_SvdUnf->SetBinContent(i,muSVDUnfErriMerge[i]);
        h_UnfoldBias->SetBinContent(i,muUnfBiasErriMerge[i]);
        //h_TotalSyst->SetBinContent(i,muTotalSystiMerge[i]);
        //h_LumiSyst->SetBinContent(i,muLumiErriMerge[i]);
        
        h_PowhegPDF->SetBinContent(i,PowhegPDFErriMerge[i]);
	h_PowhegStat->SetBinContent(i,PowhegStatErrI[i]);
      }else if(i>=5){
	h_tracksig->SetBinContent(i,mutracksigErri[i+1]);
	h_trackbck->SetBinContent(i,mutrackbckErri[i+1]);
        h_idisosig->SetBinContent(i,muidisosigErri[i+1]);
        h_idisobck->SetBinContent(i,muidisobckErri[i+1]);
	//h_track->SetBinContent(i,mutrackErri[i+1]);
        //h_idiso->SetBinContent(i,muidisoErri[i+1]);
        h_toy->SetBinContent(i,mutoyErri[i+1]);
        h_POG->SetBinContent(i,muPOGErri[i+1]);
        //h_TotalEff->SetBinContent(i,muEffErri[i+1]);
        h_met->SetBinContent(i,muMetErri[i+1]);
        h_scale->SetBinContent(i,muscaleErri[i+1]);
        h_smear->SetBinContent(i,musmearErri[i+1]);
        //h_EnMomRes->SetBinContent(i,muMomResErri[i+1]);
        h_qcdbckgr->SetBinContent(i,muQCDBckErri[i+1]);
        h_qcdshape->SetBinContent(i,muQCDShapeErri[i+1]);
        h_ewk->SetBinContent(i,muEWKErri[i+1]);
        h_fsr->SetBinContent(i,muFSRErri[i+1]);
        h_SvdUnf->SetBinContent(i,muSVDUnfErri[i+1]);
        h_UnfoldBias->SetBinContent(i,muUnfBiasErri[i+1]);
        //h_TotalSyst->SetBinContent(i,muTotalSysti[i+1]);
        //h_LumiSyst->SetBinContent(i,muLumiErri[i+1]);
        
        h_PowhegPDF->SetBinContent(i,PowhegPDFErri[i+1]);
	h_PowhegStat->SetBinContent(i,PowhegStatErrI[i]);
      }
    }

    h_Stat->Write();

    h_tracksig->Write();
    h_trackbck->Write();

    h_idisosig->Write();
    h_idisobck->Write();
    
   // h_track->Write();
   // h_idiso->Write();

    h_toy->Write();
    //
    h_POG->Write();
    
    //h_TotalEff->Write();

    h_met->Write();
                   
    h_scale->Write();
                   
    h_smear->Write();
                   
    //h_EnMomRes->Write();
                   
    h_qcdbckgr->Write();
                   
    h_qcdshape->Write();
                   
    h_ewk->Write();
                   
    h_fsr->Write();
                   
    h_SvdUnf->Write();
                   
    h_UnfoldBias->Write();
                   
    //h_TotalSyst->Write();
                   
    h_PowhegPDF->Write();
    h_PowhegStat->Write();
   
    //h_LumiSyst->Write();

}
