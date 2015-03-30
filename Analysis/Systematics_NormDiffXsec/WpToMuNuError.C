{
#include "Utils/const.h"
#include "TMathBase.h"
//*

  double mutracksigp[14]={0};
  mutracksigp[1] =0.02454  ;
  mutracksigp[2] =0.00450  ;
  mutracksigp[3] =0.02576  ;
  mutracksigp[4] =0.04923  ;
  mutracksigp[5] =0.05218  ;
  mutracksigp[6] =0.03528  ;
  mutracksigp[7] =0.00713  ;
  mutracksigp[8] =0.02326  ;
  mutracksigp[9] =0.05037  ;
  mutracksigp[10]=0.07215  ;
  mutracksigp[11]=0.08839  ;
  mutracksigp[12]=0.09914  ;
  mutracksigp[13]=0.10449  ;

  double mutracksigm[14] ={0};
  mutracksigm[1] =0.07757 ; 
  mutracksigm[2] =0.00913 ;
  mutracksigm[3] =0.06609 ;
  mutracksigm[4] =0.09579 ;
  mutracksigm[5] =0.07758 ;
  mutracksigm[6] =0.04768 ;
  mutracksigm[7] =0.03414 ;
  mutracksigm[8] =0.03959 ;
  mutracksigm[9] =0.05648 ;
  mutracksigm[10]=0.07550 ;
  mutracksigm[11]=0.09211 ;
  mutracksigm[12]=0.10406 ;
  mutracksigm[13]=0.11023 ;
 
  double mutracksigpMerge[14]={0};
  double mutracksigmMerge[14]={0};

  mutracksigpMerge[4]=TMath::Max(mutracksigp[4],mutracksigp[5]);
  mutracksigmMerge[4]=TMath::Max(mutracksigm[4],mutracksigm[5]);

  double mutrackbckp[14] = {0};
  mutrackbckp[1 ]=0.02409 ;
  mutrackbckp[2 ]=0.01113 ;
  mutrackbckp[3 ]=0.01077 ;
  mutrackbckp[4 ]=0.03225 ;
  mutrackbckp[5 ]=0.04367 ;
  mutrackbckp[6 ]=0.04293 ;
  mutrackbckp[7 ]=0.03359 ;
  mutrackbckp[8 ]=0.02077 ;
  mutrackbckp[9 ]=0.00819 ;
  mutrackbckp[10]=0.00231 ;
  mutrackbckp[11]=0.01027 ;
  mutrackbckp[12]=0.01558 ;
  mutrackbckp[13]=0.01823 ;
 
  double mutrackbckm[14] ={0};
  mutrackbckm[1 ]=0.03392 ;
  mutrackbckm[2 ]=0.01909 ;
  mutrackbckm[3 ]=0.06511 ;
  mutrackbckm[4 ]=0.06048 ;
  mutrackbckm[5 ]=0.01440 ;
  mutrackbckm[6 ]=0.03182 ;
  mutrackbckm[7 ]=0.05172 ;
  mutrackbckm[8 ]=0.04698 ;
  mutrackbckm[9 ]=0.02869 ; 
  mutrackbckm[10]=0.00800 ;
  mutrackbckm[11]=0.01004 ;
  mutrackbckm[12]=0.02301 ;
  mutrackbckm[13]=0.02970 ;
 
  double mutrackbckpMerge[14]={0};
  double mutrackbckmMerge[14]={0};

  mutrackbckpMerge[4]=TMath::Max(mutrackbckp[4],mutrackbckp[5]);
  mutrackbckmMerge[4]=TMath::Max(mutrackbckm[4],mutrackbckm[5]);



  double muidisosigp[14] = {0};
  muidisosigp[1 ] =0.00881 ;
  muidisosigp[2 ] =0.00245 ;
  muidisosigp[3 ] =0.00741 ;
  muidisosigp[4 ] =0.01564 ;
  muidisosigp[5 ] =0.01783 ;
  muidisosigp[6 ] =0.01380 ;
  muidisosigp[7 ] =0.00594 ;
  muidisosigp[8 ] =0.00294 ;
  muidisosigp[9 ] =0.01102 ;
  muidisosigp[10] =0.01754 ;
  muidisosigp[11] =0.02240 ;
  muidisosigp[12] =0.02561 ;
  muidisosigp[13] =0.02721 ;

  double muidisosigm[14] = {0};
  muidisosigm[1 ] =0.05884 ;
  muidisosigm[2 ] =0.01613 ; 
  muidisosigm[3 ] =0.09438 ;
  muidisosigm[4 ] =0.11377 ;
  muidisosigm[5 ] =0.06960 ;
  muidisosigm[6 ] =0.00328 ;
  muidisosigm[7 ] =0.05140 ;
  muidisosigm[8 ] =0.08792 ;
  muidisosigm[9 ] =0.11097 ;
  muidisosigm[10] =0.12694 ;
  muidisosigm[11] =0.13823 ;
  muidisosigm[12] =0.14560 ;
  muidisosigm[13] =0.14926 ;

  double muidisosigpMerge[14]={0};
  double muidisosigmMerge[14]={0};

  muidisosigpMerge[4]=TMath::Max(muidisosigp[4],muidisosigp[5]);
  muidisosigmMerge[4]=TMath::Max(muidisosigm[4],muidisosigm[5]);

  double muidisobckp[14]={0};
  muidisobckp[1 ] =0.01475 ;
  muidisobckp[2 ] =0.01739 ;
  muidisobckp[3 ] =0.01639 ;
  muidisobckp[4 ] =0.00618 ;
  muidisobckp[5 ] =0.01343 ;
  muidisobckp[6 ] =0.03767 ;
  muidisobckp[7 ] =0.06185 ;
  muidisobckp[8 ] =0.08327 ;
  muidisobckp[9 ] =0.10085 ;
  muidisobckp[10] =0.11446 ;
  muidisobckp[11] =0.12443 ;
  muidisobckp[12] =0.13096 ;
  muidisobckp[13] =0.13419 ;

  double muidisobckm[14] = {0};
  muidisobckm[1 ] =0.05311 ;
  muidisobckm[2 ] =0.02141 ;
  muidisobckm[3 ] =0.09772 ;
  muidisobckm[4 ] =0.11284 ;
  muidisobckm[5 ] =0.06258 ;
  muidisobckm[6 ] =0.00998 ;
  muidisobckm[7 ] =0.06961 ;
  muidisobckm[8 ] =0.10916 ;
  muidisobckm[9 ] =0.13341 ;
  muidisobckm[10] =0.14944 ;
  muidisobckm[11] =0.16032 ;
  muidisobckm[12] =0.16722 ;
  muidisobckm[13] =0.17060 ;

  double muidisobckpMerge[14]={0};
  double muidisobckmMerge[14]={0};

  muidisobckpMerge[4]=TMath::Max(muidisobckp[4],muidisobckp[5]);
  muidisobckmMerge[4]=TMath::Max(muidisobckm[4],muidisobckm[5]);

  double mutrackp[14]={0};
  double mutrackm[14]={0};
  double muidisop[14]={0};
  double muidisom[14]={0};
  for(int i(1);i<14;i++)
  {
    mutrackp[i]  = TMath::Sqrt(mutracksigp[i] *mutracksigp[i]  + mutrackbckp[i]*mutrackbckp[i]);
    mutrackm[i]  = TMath::Sqrt(mutracksigm[i] *mutracksigm[i]  + mutrackbckm[i]*mutrackbckm[i]);

    muidisop[i]  = TMath::Sqrt(muidisosigp[i] *muidisosigp[i]  + muidisobckp[i] *muidisobckp[i]);
    muidisom[i]  = TMath::Sqrt(muidisosigm[i] *muidisosigm[i]  + muidisobckm[i] *muidisobckm[i]);
  }
 
  double mutrackpMerge[14]={0};
  double mutrackmMerge[14]={0};

  //mutrackpMerge[4]=TMath::Max(mutrackp[4],mutrackp[5]);
  //mutrackmMerge[4]=TMath::Max(mutrackm[4],mutrackm[5]);
  mutrackpMerge[4]=TMath::Sqrt(mutracksigpMerge[4] *mutracksigpMerge[4]  + mutrackbckpMerge[4]*mutrackbckpMerge[4]);
  mutrackmMerge[4]=TMath::Sqrt(mutracksigmMerge[4] *mutracksigmMerge[4]  + mutrackbckmMerge[4]*mutrackbckmMerge[4]);

  double muidisopMerge[14]={0};
  double muidisomMerge[14]={0};

  //muidisopMerge[4]=TMath::Max(muidisop[4],muidisop[5]);
  //muidisomMerge[4]=TMath::Max(muidisom[4],muidisom[5]);
  muidisopMerge[4]=TMath::Sqrt(muidisosigpMerge[4] *muidisosigpMerge[4]  + muidisobckpMerge[4] *muidisobckpMerge[4]);
  muidisomMerge[4]=TMath::Sqrt(muidisosigmMerge[4] *muidisosigmMerge[4]  + muidisobckmMerge[4] *muidisobckmMerge[4]);

 // //// Variate stat error by Up and Down syst
 // double mutoyp[14] = {0};
 // mutoyp[1 ] =0.2257 ;
 // mutoyp[2 ] =0.1141 ;
 // mutoyp[3 ] =0.1515 ;
 // mutoyp[4 ] =0.3936 ;
 // mutoyp[5 ] =0.4923 ;
 // mutoyp[6 ] =0.4173 ;
 // mutoyp[7 ] =0.3111 ;
 // mutoyp[8 ] =0.1999 ;
 // mutoyp[9 ] =0.2369 ;
 // mutoyp[10] =0.4190 ;
 // mutoyp[11] =0.5568 ;
 // mutoyp[12] =0.6494 ;
 // mutoyp[13] =0.6959 ;

 // double mutoym[14] = {0};
 // mutoym[1 ] =0.2003 ;
 // mutoym[2 ] =0.1167 ;
 // mutoym[3 ] =0.0321 ;
 // mutoym[4 ] =0.1583 ;
 // mutoym[5 ] =0.2338 ;
 // mutoym[6 ] =0.2848 ;
 // mutoym[7 ] =0.3722 ;
 // mutoym[8 ] =0.4894 ;
 // mutoym[9 ] =0.6084 ;
 // mutoym[10] =0.7083 ;
 // mutoym[11] =0.7846 ;
 // mutoym[12] =0.8359 ;
 // mutoym[13] =0.8617 ;

  //// Toy, variate stat error by gauss, CovMat syst
  double mutoyp[14] = {0};
  mutoyp[1 ] =0.2893 ;
  mutoyp[2 ] =0.1382 ;
  mutoyp[3 ] =0.1891 ;
  mutoyp[4 ] =0.3915 ;
  mutoyp[5 ] =0.5051 ;
  mutoyp[6 ] =0.5148 ;
  mutoyp[7 ] =0.4482 ;
  mutoyp[8 ] =0.3696 ;
  mutoyp[9 ] =0.3287 ;
  mutoyp[10] =0.3298 ;
  mutoyp[11] =0.3503 ;
  mutoyp[12] =0.3716 ;
  mutoyp[13] =0.3842 ;

  double mutoym[14] = {0};
  mutoym[1 ] =0.1801 ;
  mutoym[2 ] =0.1167 ;
  mutoym[3 ] =0.1246 ;
  mutoym[4 ] =0.2328 ;
  mutoym[5 ] =0.3069 ;
  mutoym[6 ] =0.3302 ;
  mutoym[7 ] =0.3192 ;
  mutoym[8 ] =0.3006 ;
  mutoym[9 ] =0.2925 ;
  mutoym[10] =0.2990 ;
  mutoym[11] =0.3129 ;
  mutoym[12] =0.3265 ;
  mutoym[13] =0.3345 ;
  
  double mutoypMerge[14]={0};
  double mutoymMerge[14]={0};

  mutoypMerge[4]=TMath::Max(mutoyp[4],mutoyp[5]);
  mutoymMerge[4]=TMath::Max(mutoym[4],mutoym[5]);

  double muPOGp[14] = {0};
  muPOGp[1 ] =0.133661663 ;
  muPOGp[2 ] =0.064522038 ;
  muPOGp[3 ] =0.127689268 ;
  muPOGp[4 ] =0.33167295  ;
  muPOGp[5 ] =0.372987565 ;
  muPOGp[6 ] =0.231466785 ;
  muPOGp[7 ] =0.018961462 ;
  muPOGp[8 ] =0.200546027 ;
  muPOGp[9 ] =0.318957936 ;
  muPOGp[10] =0.372069411 ;
  muPOGp[11] =0.389084157 ;
  muPOGp[12] =0.39102901  ;
  muPOGp[13] =0.38978387  ;

  double muPOGm[14] = {0};
  muPOGm[1 ] =0.058208724 ;
  muPOGm[2 ] =0.102663987 ;
  muPOGm[3 ] =0.18257866  ;
  muPOGm[4 ] =0.054788539 ;
  muPOGm[5 ] =0.160966607 ;
  muPOGm[6 ] =0.263767445 ;
  muPOGm[7 ] =0.193651562 ;
  muPOGm[8 ] =0.034530617 ;
  muPOGm[9 ] =0.146917268 ;
  muPOGm[10] =0.286757072 ;
  muPOGm[11] =0.390624939 ;
  muPOGm[12] =0.459443251 ;
  muPOGm[13] =0.493780731 ;

  double muPOGpMerge[14]={0};
  double muPOGmMerge[14]={0};

  muPOGpMerge[4]=TMath::Max(muPOGp[4],muPOGp[5]);
  muPOGmMerge[4]=TMath::Max(muPOGm[4],muPOGm[5]);

  double mutotaleffp[14]={0};
  double mutotaleffm[14]={0};
  for(int i(1);i<14;i++)
  {
    mutotaleffp[i]  = sqrt(mutrackp[i] *mutrackp[i] +mutoyp[i] *mutoyp[i] +muidisop[i] *muidisop[i] +muPOGp[i] *muPOGp[i]);
    mutotaleffm[i]  = sqrt(mutrackm[i] *mutrackm[i] +mutoym[i] *mutoym[i] +muidisom[i] *muidisom[i] +muPOGm[i] *muPOGm[i]);
  } 
  

  double mutotaleffpMerge[14]={0};
  double mutotaleffmMerge[14]={0};

  //mutotaleffpMerge[4]=mutotaleffp[4]+mutotaleffp[5];
  //mutotaleffmMerge[4]=mutotaleffm[4]+mutotaleffm[5];
  //mutotaleffpMerge[4]=TMath::Max(mutotaleffp[4],mutotaleffp[5]);
  //mutotaleffmMerge[4]=TMath::Max(mutotaleffm[4],mutotaleffm[5]);
  mutotaleffpMerge[4]  = sqrt(mutrackpMerge[4] *mutrackpMerge[4] +mutoypMerge[4] *mutoypMerge[4] +muidisopMerge[4] *muidisopMerge[4] +muPOGpMerge[4] *muPOGpMerge[4]);
  mutotaleffmMerge[4]  = sqrt(mutrackmMerge[4] *mutrackmMerge[4] +mutoymMerge[4] *mutoymMerge[4] +muidisomMerge[4] *muidisomMerge[4] +muPOGmMerge[4] *muPOGmMerge[4]);

  //// Stat error from Toy RMS
  double mustatp[14] = {0};
  mustatp[1 ]=0.6551   ; 
  mustatp[2 ]=0.8448   ;
  mustatp[3 ]=1.0135   ;
  mustatp[4 ]=1.1148   ;
  mustatp[5 ]=1.4970   ;
  mustatp[6 ]=1.4750   ;
  mustatp[7 ]=2.0859   ;
  mustatp[8 ]=2.0329   ;
  mustatp[9 ]=2.6863   ;
  mustatp[10]=5.4655   ;
  mustatp[11]=10.4859  ;
  mustatp[12]=18.0951  ;
  mustatp[13]=24.2087  ;

  double mustatm[14] = {0};
  mustatm[1 ] =0.8302  ;
  mustatm[2 ] =1.0233  ;
  mustatm[3 ] =1.2421  ;
  mustatm[4 ] =1.3381  ;
  mustatm[5 ] =1.7568  ;
  mustatm[6 ] =1.7759  ;
  mustatm[7 ] =2.4416  ;
  mustatm[8 ] =2.4307  ;
  mustatm[9 ] =3.0754  ;
  mustatm[10] =6.2777  ;
  mustatm[11] =12.1005 ;
  mustatm[12] =16.2326 ;
  mustatm[13] =33.3904 ;

  double mustatpMerge[14] = {0};
  double mustatmMerge[14] = {0};
   mustatpMerge[4] = sqrt( mustatp[4]* mustatp[4] + mustatp[5]* mustatp[5]);
   mustatmMerge[4] = sqrt( mustatm[4]* mustatm[4] + mustatm[5]* mustatm[5]);

  double mumetp[14] = {0};
  mumetp[1 ] = 0.05468997; 
  mumetp[2 ] = 0.02971837;
  mumetp[3 ] = 0.06181037;
  mumetp[4 ] = 0.09433043;
  mumetp[5 ] = 0.10000429;
  mumetp[6 ] = 0.09811283;
  mumetp[7 ] = 0.09758622;
  mumetp[8 ] = 0.10308453;
  mumetp[9 ] = 0.12197750;
  mumetp[10] = 0.15185826;
  mumetp[11] = 0.18209261;
  mumetp[12] = 0.20525547;
  mumetp[13] = 0.21753188;

  double mumetm[14] ={0};
  mumetm[1 ] = 0.03276728 ;
  mumetm[2 ] = 0.02071402 ;
  mumetm[3 ] = 0.04285863 ;
  mumetm[4 ] = 0.05990856 ;
  mumetm[5 ] = 0.06130073 ;
  mumetm[6 ] = 0.06091243 ;
  mumetm[7 ] = 0.06802467 ;
  mumetm[8 ] = 0.09529188 ;
  mumetm[9 ] = 0.14614713 ;
  mumetm[10] = 0.20751724 ;
  mumetm[11] = 0.26424269 ;
  mumetm[12] = 0.30651657 ;
  mumetm[13] = 0.32874582 ;
  
  double mumetpMerge[14] = {0};
  double mumetmMerge[14] = {0};
  //mumetpMerge[4]=mumetp[4 ] + mumetp[5 ];
  //mumetmMerge[4]=mumetm[4 ] + mumetm[5 ];
  mumetpMerge[4]=TMath::Max(mumetp[4 ],mumetp[5 ]);
  mumetmMerge[4]=TMath::Max(mumetm[4 ],mumetm[5 ]);

  double muscalep[14] = {0};
  muscalep[1 ] =0.2619  ; 
  muscalep[2 ] =0.1636  ;
  muscalep[3 ] =0.0816  ;
  muscalep[4 ] =0.2992  ;
  muscalep[5 ] =0.4800  ;
  muscalep[6 ] =0.5406  ;
  muscalep[7 ] =0.5054  ;
  muscalep[8 ] =0.4262  ;
  muscalep[9 ] =0.3421  ;
  muscalep[10] = 0.2710 ;
  muscalep[11] = 0.3115 ;
  muscalep[12] = 0.3539 ;
  muscalep[13] = 0.3750 ;

  double muscalem[14] ={0};
  muscalem[1 ] =0.1943  ; 
  muscalem[2 ] =0.1093  ;
  muscalem[3 ] =0.0760  ;
  muscalem[4 ] =0.1752  ;
  muscalem[5 ] =0.2675  ;
  muscalem[6 ] =0.3160  ;
  muscalem[7 ] =0.3443  ;
  muscalem[8 ] =0.3613  ;
  muscalem[9 ] =0.3654  ;
  muscalem[10] = 0.3585 ;
  muscalem[11] = 0.3517 ;
  muscalem[12] = 0.3919 ;
  muscalem[13] = 0.4124 ;

  double muscalepMerge[14] = {0};
  double muscalemMerge[14] = {0};
  muscalepMerge[4]=TMath::Max(muscalep[4 ],muscalep[5 ]);
  muscalemMerge[4]=TMath::Max(muscalem[4 ],muscalem[5 ]);

  double musmearp[14] ={0};
  musmearp[1 ] =0.2345   ; 
  musmearp[2 ] =0.1202   ;
  musmearp[3 ] =0.1179   ;
  musmearp[4 ] =0.3430   ;
  musmearp[5 ] =0.4710   ;
  musmearp[6 ] =0.4525   ;
  musmearp[7 ] =0.3312   ;
  musmearp[8 ] =0.1726   ;
  musmearp[9 ] =0.0958 ;
  musmearp[10] =0.2179  ;
  musmearp[11] =0.3076  ;
  musmearp[12] =0.3668  ;
  musmearp[13] =0.3963  ;

  double musmearm[14] = {0};
  musmearm[1 ] =0.1812  ; 
  musmearm[2 ] =0.0966 ;
  musmearm[3 ] =0.0446  ;
  musmearm[4 ] =0.1715  ;
  musmearm[5 ] =0.2374  ;
  musmearm[6 ] =0.2642 ;
  musmearm[7 ] =0.2959  ;
  musmearm[8 ] =0.3387  ;
  musmearm[9 ] =0.3866  ;
  musmearm[10] =0.4280 ;
  musmearm[11] =0.4599 ;
  musmearm[12] =0.4815 ;
  musmearm[13] =0.4923 ;

  double musmearpMerge[14] = {0};
  double musmearmMerge[14] = {0};
  musmearpMerge[4]=TMath::Max(musmearp[4 ],musmearp[5 ]);
  musmearmMerge[4]=TMath::Max(musmearm[4 ],musmearm[5 ]);

  double muMomResp[14]={0};
  double muMomResm[14]={0};
  for(int i(1);i<14;i++)
  {
    muMomResp[i]  = sqrt(musmearp[i] *musmearp[i] +muscalep[i] *muscalep[i] ); 
    muMomResm[i]  = sqrt(musmearm[i] *musmearm[i] +muscalem[i] *muscalem[i] ); 
  }
  
  double muMomRespMerge[14]={0};
  double muMomResmMerge[14]={0};
  //muMomRespMerge[4]=muMomResp[4] +muMomResp[5]; 
  //muMomResmMerge[4]=muMomResm[4] +muMomResm[5]; 
  //muMomRespMerge[4]=TMath::Max(muMomResp[4],muMomResp[5]); 
  //muMomResmMerge[4]=TMath::Max(muMomResm[4],muMomResm[5]); 
  muMomRespMerge[4]  = sqrt(musmearpMerge[4] *musmearpMerge[4] +muscalepMerge[4] *muscalepMerge[4] ); 
  muMomResmMerge[4]  = sqrt(musmearmMerge[4] *musmearmMerge[4] +muscalemMerge[4] *muscalemMerge[4] ); 
  
  double muqcdbckgrp[14]={0};
  muqcdbckgrp[1 ] =0.8581 ;
  muqcdbckgrp[2 ] =1.3639 ;
  muqcdbckgrp[3 ] =1.1574 ;
  muqcdbckgrp[4 ] =1.3953 ;
  muqcdbckgrp[5 ] =1.3276 ;
  muqcdbckgrp[6 ] =2.2050 ;
  muqcdbckgrp[7 ] =1.2643 ;
  muqcdbckgrp[8 ] =1.8920 ;
  muqcdbckgrp[9 ] =0.9038 ; 
  muqcdbckgrp[10] =0.9099 ;
  muqcdbckgrp[11] =0.9221 ;
  muqcdbckgrp[12] =0.9389 ;
  muqcdbckgrp[13] =0.9572 ;
 
  double muqcdbckgrm[14]={0};
  muqcdbckgrm[1 ] =0.8909 ;
  muqcdbckgrm[2 ] =1.2335 ;
  muqcdbckgrm[3 ] =1.3089 ;
  muqcdbckgrm[4 ] =1.0940 ;
  muqcdbckgrm[5 ] =1.2754 ;
  muqcdbckgrm[6 ] =1.9296 ;
  muqcdbckgrm[7 ] =1.1988 ;
  muqcdbckgrm[8 ] =2.2915 ;
  muqcdbckgrm[9 ] =1.0062 ; 
  muqcdbckgrm[10] =1.0250 ;
  muqcdbckgrm[11] =1.0400 ;
  muqcdbckgrm[12] =1.0571 ;
  muqcdbckgrm[13] =1.0790 ;

  double muqcdbckgrpMerge[14]={0};
  double muqcdbckgrmMerge[14]={0};
   //muqcdbckgrpMerge[4]=muqcdbckgrp[4 ] + muqcdbckgrp[5 ];
   //muqcdbckgrmMerge[4]=muqcdbckgrm[4 ] + muqcdbckgrm[5 ];
   muqcdbckgrpMerge[4]=TMath::Max(muqcdbckgrp[4 ],muqcdbckgrp[5 ]);
   muqcdbckgrmMerge[4]=TMath::Max(muqcdbckgrm[4 ],muqcdbckgrm[5 ]);
  
  double muqcdshapep[14] = {0};
  muqcdshapep[1 ] =0.2599 ;
  muqcdshapep[2 ] =0.3915 ;
  muqcdshapep[3 ] =0.2549 ;
  muqcdshapep[4 ] =0.3903 ;
  muqcdshapep[5 ] =0.4437 ;
  muqcdshapep[6 ] =0.4100 ;
  muqcdshapep[7 ] =0.2388 ;
  muqcdshapep[8 ] =0.3322 ;
  muqcdshapep[9 ] =0.3805 ; 
  muqcdshapep[10] =0.2195 ;
  muqcdshapep[11] =0.8909 ;
  muqcdshapep[12] =0.9472 ;
  muqcdshapep[13] =0.9389 ;
 
  double muqcdshapem[14] ={0};
  muqcdshapem[1 ] =0.1270 ;
  muqcdshapem[2 ] =0.2613 ;
  muqcdshapem[3 ] =0.3828 ;
  muqcdshapem[4 ] =0.3288 ;
  muqcdshapem[5 ] =0.2562 ;
  muqcdshapem[6 ] =0.2243 ;
  muqcdshapem[7 ] =0.1486 ;
  muqcdshapem[8 ] =0.5703 ;
  muqcdshapem[9 ] =0.2721 ; 
  muqcdshapem[10] =0.1488 ;
  muqcdshapem[11] =0.8452 ;
  muqcdshapem[12] =0.8953 ;
  muqcdshapem[13] =0.8980 ;

  double muqcdshapepMerge[14] = {0};
  double muqcdshapemMerge[14] = {0};
  //muqcdshapepMerge[4]=muqcdshapep[4 ] +muqcdshapep[5 ];
  //muqcdshapemMerge[4]=muqcdshapem[4 ] +muqcdshapem[5 ];
  muqcdshapepMerge[4]=TMath::Max(muqcdshapep[4 ],muqcdshapep[5 ]);
  muqcdshapemMerge[4]=TMath::Max(muqcdshapem[4 ],muqcdshapem[5 ]);
  
  double muewkp[14] ={0};
  muewkp[1 ] = 0.007944514;
  muewkp[2 ] = 0.015402191;
  muewkp[3 ] = 0.029667537;
  muewkp[4 ] = 0.035121542;
  muewkp[5 ] = 0.023169638;
  muewkp[6 ] = 0.030295625;
  muewkp[7 ] = 0.063217132;
  muewkp[8 ] = 0.099424685;
  muewkp[9 ] = 0.132093526; 
  muewkp[10] = 0.158809094;
  muewkp[11] = 0.179019613;
  muewkp[12] = 0.192534523;
  muewkp[13] = 0.199288826;
 
  double muewkm[14] ={0};
  muewkm[1 ] = 0.008630046;
  muewkm[2 ] = 0.018503291;
  muewkm[3 ] = 0.040513953;
  muewkm[4 ] = 0.04565623;
  muewkm[5 ] = 0.025722909;
  muewkm[6 ] = 0.033626406;
  muewkm[7 ] = 0.075325799;
  muewkm[8 ] = 0.108321835;
  muewkm[9 ] = 0.129368051; 
  muewkm[10] = 0.1428973;
  muewkm[11] = 0.152708486;
  muewkm[12] = 0.159829944;
  muewkm[13] = 0.16350779;

  double muewkpMerge[14] ={0};
  double muewkmMerge[14] ={0};
  //muewkpMerge[4]=muewkp[4 ] + muewkp[5 ] ; 
  //muewkmMerge[4]=muewkm[4 ] + muewkm[5 ] ; 
  muewkpMerge[4]=TMath::Max(muewkp[4 ], muewkp[5 ] ); 
  muewkmMerge[4]=TMath::Max(muewkm[4 ] , muewkm[5 ] ); 
  
  double mufsrp[14] ={0};
  mufsrp[1 ] =0.0034  ;
  mufsrp[2 ] =0.0038  ;
  mufsrp[3 ] =0.0027  ;
  mufsrp[4 ] =0.0006  ;
  mufsrp[5 ] =0.0052  ;
  mufsrp[6 ] =0.0093  ;
  mufsrp[7 ] =0.0125  ;
  mufsrp[8 ] =0.0147  ;
  mufsrp[9 ] =0.0163  ;
  mufsrp[10] = 0.0175 ;
  mufsrp[11] = 0.0183 ;
  mufsrp[12] = 0.0188 ;
  mufsrp[13] = 0.0190 ;
 
  double mufsrm[14] ={0};
  mufsrm[1 ] =0.0002  ;
  mufsrm[2 ] =0.0017  ;
  mufsrm[3 ] =0.0038  ;
  mufsrm[4 ] =0.0041  ;
  mufsrm[5 ] =0.0022  ;
  mufsrm[6 ] =0.0015  ;
  mufsrm[7 ] =0.0063  ;
  mufsrm[8 ] =0.0111  ;
  mufsrm[9 ] =0.0153  ;
  mufsrm[10] = 0.0187 ;
  mufsrm[11] = 0.0212 ;
  mufsrm[12] = 0.0228 ;
  mufsrm[13] = 0.0236 ;


  double mufsrpMerge[14] ={0};
  double mufsrmMerge[14] ={0};
  //mufsrpMerge[4]=mufsrp[4 ] + mufsrp[5 ];
  //mufsrmMerge[4]=mufsrm[4 ] + mufsrm[5 ];
  mufsrpMerge[4]=TMath::Max(mufsrp[4 ] , mufsrp[5 ]);
  mufsrmMerge[4]=TMath::Max(mufsrm[4 ] , mufsrm[5 ]);

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
  //musvdunfpMerge[4]=musvdunfp[4 ]+ musvdunfp[5 ];
  //musvdunfmMerge[4]=musvdunfm[4 ]+ musvdunfm[5 ];
  musvdunfpMerge[4]=TMath::Max(musvdunfp[4 ], musvdunfp[5 ]);
  musvdunfmMerge[4]=TMath::Max(musvdunfm[4 ], musvdunfm[5 ]);
  
  double muBiasp[14] = {0};
  muBiasp[1 ] =0.1529   ;
  muBiasp[2 ] =0.1290   ;
  muBiasp[3 ] =0.1524   ;
  muBiasp[4 ] =0.1167   ;
  muBiasp[5 ] =0.1466   ;
  muBiasp[6 ] =0.2724   ;
  muBiasp[7 ] =0.4469   ;
  muBiasp[8 ] =0.6328   ;
  muBiasp[9 ] =0.7921   ;
  muBiasp[10] = 0.9076  ;
  muBiasp[11] = 0.9938  ;
  muBiasp[12] = 1.0502  ;
  muBiasp[13] = 1.0794  ;

  double muBiasm[14]={0};
  muBiasm[1 ] =0.3856   ;
  muBiasm[2 ] =0.2950   ;
  muBiasm[3 ] =0.1037   ;
  muBiasm[4 ] =0.1694   ;
  muBiasm[5 ] =0.4341   ;
  muBiasm[6 ] =0.7133   ;
  muBiasm[7 ] =0.9562   ;
  muBiasm[8 ] =1.1874   ;
  muBiasm[9 ] =1.3588   ;
  muBiasm[10] = 1.4642  ;
  muBiasm[11] = 1.5321  ;
  muBiasm[12] = 1.5702  ;
  muBiasm[13] = 1.5877  ;

  double muBiaspMerge[14] = {0};
  double muBiasmMerge[14] = {0};
  muBiaspMerge[4]=TMath::Max(muBiasp[4 ], muBiasp[5 ]);
  muBiasmMerge[4]=TMath::Max(muBiasm[4 ] , muBiasm[5 ]);

  double systtotalp[14] = {0};
  double systtotalm[14] = {0};
  for(int i(1);i<14;i++)
  {
    systtotalp[i] =sqrt(
	mutracksigp[i] *mutracksigp[i]
	+ mutrackbckp[i]*mutrackbckp[i]
	+ muidisosigp[i]*muidisosigp[i]
	+ muidisobckp[i]*muidisobckp[i]
	+ mutoyp[i]*mutoyp[i]
	+ muPOGp[i]*muPOGp[i]
	+ mumetp[i] *mumetp[i] 
	+ musmearp[i] *musmearp[i] 
	+ muscalep[i] *muscalep[i] 
	+ muqcdshapep[i] *muqcdshapep[i] 
	+ muqcdbckgrp[i] *muqcdbckgrp[i] 
	+ muewkp[i] *muewkp[i]  
	+ muBiasp[i] *muBiasp[i] 
	+ musvdunfp[i] *musvdunfp[i]  
	+ mufsrp[i] *mufsrp[i] 
	);

    systtotalm[i] =sqrt(
	mutracksigm[i] *mutracksigm[i]
	+ mutrackbckm[i]*mutrackbckm[i]
	+ muidisosigm[i]*muidisosigm[i]
	+ muidisobckm[i]*muidisobckm[i]
	+ mutoym[i]*mutoym[i]
	+ muPOGm[i]*muPOGm[i] 
	+ mumetm[i] *mumetm[i] 
	+ musmearm[i] *musmearm[i] 
	+ muscalem[i] *muscalem[i] 
	+ muqcdshapem[i] *muqcdshapem[i] 
	+ muqcdbckgrm[i] *muqcdbckgrm[i] 
	+ muewkm[i] *muewkm[i]  
	+ muBiasm[i] *muBiasm[i] 
	+ musvdunfm[i] *musvdunfm[i]  
	+ mufsrm[i] *mufsrm[i] 
	);
  }

  
  double systtotalpMerge[14] = {0};
  double systtotalmMerge[14] = {0};
    systtotalpMerge[4] =sqrt(
	mutracksigpMerge[4] *mutracksigpMerge[4]
	+ mutrackbckpMerge[4]*mutrackbckpMerge[4]
	+ muidisosigpMerge[4]*muidisosigpMerge[4]
	+ muidisobckpMerge[4]*muidisobckpMerge[4]
	+ mutoypMerge[4]*mutoypMerge[4]
	+ muPOGpMerge[4]*muPOGpMerge[4] 
	+ mumetpMerge[4] *mumetpMerge[4] 
	+ musmearpMerge[4] *musmearpMerge[4] 
	+ muscalepMerge[4] *muscalepMerge[4] 
	+ muqcdshapepMerge[4] *muqcdshapepMerge[4] 
	+ muqcdbckgrpMerge[4] *muqcdbckgrpMerge[4] 
	+ muewkpMerge[4] *muewkpMerge[4]  
	+ muBiaspMerge[4] *muBiaspMerge[4] 
	+ musvdunfpMerge[4] *musvdunfpMerge[4]  
	+ mufsrpMerge[4] *mufsrpMerge[4] 
	);
    
    systtotalmMerge[4] =sqrt(
	mutracksigmMerge[4] *mutracksigmMerge[4]
	+ mutrackbckmMerge[4]*mutrackbckmMerge[4]
	+ muidisosigmMerge[4]*muidisosigmMerge[4]
	+ muidisobckmMerge[4]*muidisobckmMerge[4]
	+ mutoymMerge[4]*mutoymMerge[4]
	+ muPOGmMerge[4]*muPOGmMerge[4] 
	+ mumetmMerge[4] *mumetmMerge[4] 
	+ musmearmMerge[4] *musmearmMerge[4] 
	+ muscalemMerge[4] *muscalemMerge[4] 
	+ muqcdshapemMerge[4] *muqcdshapemMerge[4] 
	+ muqcdbckgrmMerge[4] *muqcdbckgrmMerge[4] 
	+ muewkmMerge[4] *muewkmMerge[4]  
	+ muBiasmMerge[4] *muBiasmMerge[4] 
	+ musvdunfmMerge[4] *musvdunfmMerge[4]  
	+ mufsrmMerge[4] *mufsrmMerge[4] 
	);
    

//  Print out result
  cout<<fixed<<setprecision(3);

  cout<<"TrackSigShape:\t\t"<< mutracksigpMerge[4]    << "\t\t" << mutracksigmMerge[4] <<endl;  
  cout<<"TrackBckShape:\t\t"<< mutrackbckpMerge[4]    << "\t\t" << mutrackbckmMerge[4] <<endl;  
  cout<<"IdisoSigShape:\t\t"<< muidisosigpMerge[4]    << "\t\t" << muidisosigmMerge[4] <<endl;  
  cout<<"IdisoBckShape:\t\t"<< muidisobckpMerge[4]    << "\t\t" << muidisobckmMerge[4] <<endl;  
  cout<<"toy:\t\t\t"<< mutoypMerge[4]    << "\t\t" << mutoymMerge[4] <<endl;  
  cout<<"POG:\t\t\t"<< muPOGpMerge[4]    << "\t\t" << muPOGmMerge[4] <<endl;  
  cout<<"TotalEffMerge:\t\t"<< mutotaleffpMerge[4] 	<< "\t\t" << mutotaleffmMerge[4] <<endl;

  cout<<"" <<endl;  
  cout<<"scale:\t\t\t"<< muscalepMerge[4]      << "\t\t" << muscalemMerge[4] <<endl;
  cout<<"smear:\t\t\t"<< musmearpMerge[4]      << "\t\t" << musmearmMerge[4] <<endl;
  cout<<"EnMomRes:\t\t"<< muMomRespMerge[4]	<< "\t\t" << muMomResmMerge[4] <<endl;
  
  cout<<"" <<endl;  
  cout<<"MET:\t\t\t"<< mumetpMerge[4] 	<< "\t\t" << mumetmMerge[4] <<endl;
  cout<<"QcdBckgr:\t\t"<< muqcdbckgrpMerge[4] 	<< "\t\t" << muqcdbckgrmMerge[4] <<endl;
  cout<<"QcdShape:\t\t"<< muqcdshapepMerge[4]	<< "\t\t" << muqcdshapemMerge[4] <<endl;
  cout<<"EWK:\t\t\t"<< muewkpMerge[4] 	<< "\t\t" << muewkmMerge[4] <<endl;
  cout<<"SvdUnf:\t\t\t"<< musvdunfpMerge[4] 	<< "\t\t" << musvdunfmMerge[4] <<endl;
  cout<<"FSR:\t\t\t"<< mufsrpMerge[4] 	<< "\t\t" << mufsrmMerge[4] <<endl;
  cout<<"UnfBias:\t\t"<< muBiaspMerge[4] 	<< "\t\t" << muBiasmMerge[4] <<endl;
  
  cout<<"" <<endl;  
  cout<<"TotalSyst:\t\t"<< systtotalpMerge[4] 	<< "\t\t" << systtotalmMerge[4] <<endl;
  
  cout<<"" <<endl;  
  cout<<"Stat:\t\t\t"<< mustatpMerge[4] 	<< "\t\t" << mustatmMerge[4] <<endl;
  
  double totaluncerpMerge[14]={0};
  double totaluncermMerge[14]={0};
  totaluncerpMerge[4]= sqrt( mustatpMerge[4]*mustatpMerge[4] +systtotalpMerge[4]*systtotalpMerge[4] );
  totaluncermMerge[4]= sqrt( mustatmMerge[4]*mustatmMerge[4] +systtotalmMerge[4]*systtotalmMerge[4] );
  cout<< "Merged Tota Uncer:\t"<< totaluncerpMerge[4] 	<< "\t\t" << totaluncermMerge[4] <<endl;
  
  
  double totaluncerp[14]={0};
  double totaluncerm[14]={0};
  for(int i(1);i<14;i++)
  {
    totaluncerp[i] = sqrt(mustatp[i] *mustatp[i] +systtotalp[i] *systtotalp[i]);
    totaluncerm[i] = sqrt(mustatm[i] *mustatm[i] +systtotalm[i] *systtotalm[i]);
  }
  
  
  // Recon. 
  cout<<" Lepton Track Sig Shape Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< mutracksigp[i] << "\t" << mutracksigm[i] <<endl;
  }
  cout<<" Lepton Track Bkgr Shape Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< mutrackbckp[i] << "\t" << mutrackbckm[i] <<endl;
  }
  cout<<" Lepton ID-ISO Sig Shape Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< muidisosigp[i] << "\t" << muidisosigm[i] <<endl;
  }
  cout<<" Lepton ID-ISO Bkgr Shape Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< muidisobckp[i] << "\t" << muidisobckm[i] <<endl;
  }
  cout<<" Lepton Toy CovMat Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< mutoyp[i] << "\t" << mutoym[i] <<endl;
  }
  cout<<" Lepton MuPOG Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< muPOGp[i] << "\t" << muPOGm[i] <<endl;
  }
  cout<<"Total Lepton Recon. Effi Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< mutotaleffp[i] << "\t" << mutotaleffm[i] <<endl;
  }
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mutotaleffp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mutotaleffm[i] ;
  }
    cout << " \\"<<"\\"<<endl;
    cout << endl;
  
  //Momentum Resolution
  cout<<"MomRes.Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<muMomResp[i] <<"\t"<< muMomResm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muMomResp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muMomResm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //MET Resolution
  cout<<"METRes.Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<mumetp[i] <<"\t"<< mumetm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mumetp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mumetm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //QCD Background
  cout<<"QCD Bckgr.Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<muqcdbckgrp[i] <<"\t"<< muqcdbckgrm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muqcdbckgrp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muqcdbckgrm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //QCD Shape
  cout<<"QCD Shape. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<muqcdshapep[i] <<"\t"<< muqcdshapem[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muqcdshapep[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muqcdshapem[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //EWK
  cout<<"EWK. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<muewkp[i] <<"\t"<< muewkm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muewkp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muewkm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //SVD.Unf
  cout<<"SVD.Unf. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<musvdunfp[i] <<"\t"<< musvdunfm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << musvdunfp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << musvdunfm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //FSR
  cout<<"FSR. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<mufsrp[i] <<"\t"<< mufsrm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mufsrp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mufsrm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //Unf.Bias
  cout<<"Unf.Bias. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<muBiasp[i] <<"\t"<< muBiasm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muBiasp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << muBiasm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //Total syst
  cout<<"Total. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<systtotalp[i] <<"\t"<< systtotalm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << systtotalp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << systtotalm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //Statistical error
  cout<<"Statistical Error"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<mustatp[i] <<"\t"<< mustatm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mustatp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << mustatm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //Total Uncertainty
  cout<<"Total Uncertainty"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<totaluncerp[i] <<"\t"<< totaluncerm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << totaluncerp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << totaluncerm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  TString resultDir = "WptXsecErrors";
  gSystem->mkdir(resultDir,kTRUE);

  TFile f_out(resultDir+"/WpToMuNuErrors.root","recreate");

  ///Write all errors to root
    
    TH1D* h_tracksig = new TH1D("h_tracksig","h_tracksig",13,0,13);
    //TH1D* h_mutracksigm = new TH1D("h_mutracksigm","h_mutracksigm",13,0,13);
    
    TH1D* h_trackbck = new TH1D("h_trackbck","h_trackbck",13,0,13);
    //TH1D* h_mutrackbckm = new TH1D("h_mutrackbckm","h_mutrackbckm",13,0,13);

    TH1D* h_idisosig = new TH1D("h_idisosig","h_idisosig",13,0,13);
    //TH1D* h_muidisosigm = new TH1D("h_muidisosigm","h_muidisosigm",13,0,13);
  
    TH1D* h_idisobck = new TH1D("h_idisobck","h_idisobck",13,0,13);
    //TH1D* h_muidisobckm = new TH1D("h_muidisobckm","h_muidisobckm",13,0,13);
  
    TH1D* h_track = new TH1D("h_track","h_track",13,0,13);
    //TH1D* h_mutrackm = new TH1D("h_mutrackm","h_mutrackm",13,0,13);

    TH1D* h_idiso = new TH1D("h_idiso","h_idiso",13,0,13);
    //TH1D* h_muidisom = new TH1D("h_muidisom","h_muidisom",13,0,13);
    
    TH1D* h_toy = new TH1D("h_toy","h_toy",13,0,13);
    //TH1D* h_mutoym = new TH1D("h_mutoym","h_mutoym",13,0,13);
    
    TH1D* h_POG = new TH1D("h_POG","h_POG",13,0,13);
    //TH1D* h_muPOGm = new TH1D("h_muPOGm","h_muPOGm",13,0,13);
    
    TH1D* h_TotalEff = new TH1D("h_TotalEff","h_TotalEff",13,0,13);
    //TH1D* h_TotalEffm = new TH1D("h_TotalEffm","h_TotalEffm",13,0,13);
    
    TH1D* h_Stat = new TH1D("h_Stat","h_Stat",13,0,13);
    //TH1D* h_Statm = new TH1D("h_Statm","h_Statm",13,0,13);
    
    TH1D* h_met = new TH1D("h_met","h_met",13,0,13);
    //TH1D* h_metm = new TH1D("h_metm","h_metm",13,0,13);
    
    TH1D* h_scale = new TH1D("h_scale","h_scale",13,0,13);
    //TH1D* h_scalem = new TH1D("h_scalem","h_scalem",13,0,13);
    
    TH1D* h_smear = new TH1D("h_smear","h_smear",13,0,13);
    //TH1D* h_smearm = new TH1D("h_smearm","h_smearm",13,0,13);
    
    TH1D* h_EnMomRes = new TH1D("h_EnMomRes","h_EnMomRes",13,0,13);
    //TH1D* h_EnMomResm = new TH1D("h_EnMomResm","h_EnMomResm",13,0,13);
    
    TH1D* h_qcdbckgr = new TH1D("h_qcdbckgr","h_qcdbckgr",13,0,13);
    //TH1D* h_qcdbckgrm = new TH1D("h_qcdbckgrm","h_qcdbckgrm",13,0,13);
    
    TH1D* h_qcdshape = new TH1D("h_qcdshape","h_qcdshape",13,0,13);
    //TH1D* h_qcdshapem = new TH1D("h_qcdshapem","h_qcdshapem",13,0,13);
    
    TH1D* h_ewk = new TH1D("h_ewk","h_ewk",13,0,13);
    //TH1D* h_ewkm = new TH1D("h_ewkm","h_ewkm",13,0,13);
    
    TH1D* h_fsr = new TH1D("h_fsr","h_fsr",13,0,13);
    //TH1D* h_fsrm = new TH1D("h_fsrm","h_fsrm",13,0,13);
    
    TH1D* h_SvdUnf = new TH1D("h_SvdUnf","h_SvdUnf",13,0,13);
    //TH1D* h_SvdUnfm = new TH1D("h_SvdUnfm","h_SvdUnfm",13,0,13);
    
    TH1D* h_UnfoldBias = new TH1D("h_UnfoldBias","h_UnfoldBias",13,0,13);
    //TH1D* h_UnfoldBiasm = new TH1D("h_UnfoldBiasm","h_UnfoldBiasm",13,0,13);
    
    TH1D* h_TotalSyst = new TH1D("h_TotalSyst","h_TotalSyst",13,0,13);
    //TH1D* h_TotalSystm = new TH1D("h_TotalSystm","h_TotalSystm",13,0,13);
    
    TH1D* h_TotalUncer = new TH1D("h_TotalUncer","h_TotalUncer",13,0,13);
    //TH1D* h_TotalUncerm = new TH1D("h_muTotalUncerm","h_muTotalUncerm",13,0,13);
    

    TH1D* h_LumiSyst = new TH1D("h_LumiSyst","h_LumiSyst",13,0,13);
    //TH1D* h_LumiSystm = new TH1D("h_LumiSystm","h_LumiSystm",13,0,13);
   
    
    TH1D* h_PowhegPDF = new TH1D("h_PowhegPDF","h_PowhegPDF",12,0,12);
    //TH1D* h_PowhegPDFm = new TH1D("h_muPowhegPDFm","h_muPowhegPDFm",13,0,13);
   
    //normalized PowhegPDF 12 Bin error
    double PowhegPDFp[13];
    PowhegPDFp[0] = 0; 
    PowhegPDFp[1] =0.547639;   //3.0116; //4.268; 
    PowhegPDFp[2] =0.11084 ;   //3.7087; //4.147; 
    PowhegPDFp[3] =0.423104;   //4.0534; //4.122; 
    PowhegPDFp[4] =0.502498;   //4.1509; //4.123; 
    PowhegPDFp[5] =0.523514;   //4.3433; //4.132; 
    PowhegPDFp[6] =0.456424;   //4.3076; //4.126; 
    PowhegPDFp[7] =0.390481;   //4.4361; //4.143; 
    PowhegPDFp[8] =0.423705;   //5.4356; //5.279; 
    PowhegPDFp[9] =1.07401 ;  //4.5321; //4.222; 
    PowhegPDFp[10]=1.92522 ;   //4.8190; //4.426;
    PowhegPDFp[11]=2.57608 ;   //5.1662; //4.819; 
    PowhegPDFp[12]=3.85306 ;   //5.3727; //5.075; 
    //PowhegPDFp[13]=   //6.3736; //6.084; 

    for(int i(1);i<13;i++)
    {
      h_PowhegPDF->SetBinContent(i,PowhegPDFp[i]);
    }

    for(int i(1);i<14;i++)
    {
      h_tracksig->SetBinContent(i,mutracksigp[i]);
      //h_mutracksigm->SetBinContent(i,mutracksigm[i]);
   
      h_trackbck->SetBinContent(i,mutrackbckp[i]);
      //h_mutrackbckm->SetBinContent(i,mutrackbckm[i]);
      
      h_idisosig->SetBinContent(i,muidisosigp[i]);
      //h_muidisosigm->SetBinContent(i,muidisosigm[i]);
   
      h_idisobck->SetBinContent(i,muidisobckp[i]);
      //h_muidisobckm->SetBinContent(i,muidisobckm[i]);
   
      h_track->SetBinContent(i,mutrackp[i]);
      //h_mutrackm->SetBinContent(i,mutrackm[i]);
   
      h_idiso->SetBinContent(i,muidisop[i]);
      //h_muidisom->SetBinContent(i,muidisom[i]);
   
      h_toy->SetBinContent(i,mutoyp[i]);
      //h_mutoym->SetBinContent(i,mutoym[i]);
   
      h_POG->SetBinContent(i,muPOGp[i]);
      //h_muPOGm->SetBinContent(i,muPOGm[i]);
   
      h_TotalEff->SetBinContent(i,mutotaleffp[i]);
      //h_TotalEffm->SetBinContent(i,mutotaleffm[i]);
   
      h_Stat->SetBinContent(i,mustatp[i]);
      //h_Statm->SetBinContent(i,mustatm[i]);
   
      h_met->SetBinContent(i,mumetp[i]);
      //h_metm->SetBinContent(i,mumetm[i]);
   
      h_scale->SetBinContent(i,muscalep[i]);
      //h_scalem->SetBinContent(i,muscalem[i]);
   
      h_smear->SetBinContent(i,musmearp[i]);
      //h_smearm->SetBinContent(i,musmearm[i]);
   
      h_EnMomRes->SetBinContent(i,muMomResp[i]);
      //h_EnMomResm->SetBinContent(i,muMomResm[i]);
   
      h_qcdbckgr->SetBinContent(i,muqcdbckgrp[i]);
      //h_qcdbckgrm->SetBinContent(i,muqcdbckgrm[i]);
   
      h_qcdshape->SetBinContent(i,muqcdshapep[i]);
      //h_qcdshapem->SetBinContent(i,muqcdshapem[i]);
   
      h_ewk->SetBinContent(i,muewkp[i]);
      //h_ewkm->SetBinContent(i,muewkm[i]);
   
      h_fsr->SetBinContent(i,mufsrp[i]);
      //h_fsrm->SetBinContent(i,mufsrm[i]);
   
      h_SvdUnf->SetBinContent(i,musvdunfp[i]);
      //h_SvdUnfm->SetBinContent(i,musvdunfm[i]);
   
      h_UnfoldBias->SetBinContent(i,muBiasp[i]);
      //h_UnfoldBiasm->SetBinContent(i,muBiasm[i]);
   
      h_TotalSyst->SetBinContent(i,systtotalp[i]);
      //h_TotalSystm->SetBinContent(i,systtotalm[i]);
   
      h_TotalUncer->SetBinContent(i,totaluncerp[i]);
      //h_TotalUncerm->SetBinContent(i,totaluncerm[i]);
   
     
      h_LumiSyst->SetBinContent(i,2.6);
      //h_LumiSystm->SetBinContent(i,2.6);

    }

    h_tracksig->Write();
    //h_mutracksigm->Write();
                     // */
    h_trackbck->Write();
    //h_mutrackbckm->Write();
    
    h_idisosig->Write();
    //h_muidisosigm->Write();
                   
    h_idisobck->Write();
    //h_muidisobckm->Write();
                   
    h_track->Write();
    //h_mutrackm->Write();
                   
    h_idiso->Write();
    //h_muidisom->Write();
                   
    h_toy->Write();
    //h_mutoym->Write();
                   
    h_POG->Write();
    //h_muPOGm->Write();
                   
    h_TotalEff->Write();
    //h_TotalEffm->Write();
                   
    h_Stat->Write();
    //h_Statm->Write();
                   
    h_met->Write();
    //h_metm->Write();
                   
    h_scale->Write();
    //h_scalem->Write();
                   
    h_smear->Write();
    //h_smearm->Write();
                   
    h_EnMomRes->Write();
    //h_EnMomResm->Write();
                   
    h_qcdbckgr->Write();
    //h_qcdbckgrm->Write();
                   
    h_qcdshape->Write();
    //h_qcdshapem->Write();
                   
    h_ewk->Write();
    //h_ewkm->Write();
                   
    h_fsr->Write();
    //h_fsrm->Write();
                   
    h_SvdUnf->Write();
    //h_SvdUnfm->Write();
                   
    h_UnfoldBias->Write();
    //h_UnfoldBiasm->Write();
                   
    h_TotalSyst->Write();
    //h_TotalSystm->Write();
                   
    h_TotalUncer->Write();
    //h_TotalUncerm->Write();

    h_PowhegPDF->Write();
    //h_PowhegPDFm->Write();
    
    h_LumiSyst->Write();
    //h_LumiSystm->Write();


}
