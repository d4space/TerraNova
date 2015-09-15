{
#include "../Utils/const.h"
 //*
  double elestatp[14] = {0};
  elestatp[1 ]= 0.7750 ; 
  elestatp[2 ]= 0.9651 ;
  elestatp[3 ]= 1.1508 ;
  elestatp[4 ]= 1.2460 ;
  elestatp[5 ]= 1.6704 ;
  elestatp[6 ]= 1.6506 ;
  elestatp[7 ]= 2.2474 ;
  elestatp[8 ]= 2.3459 ;
  elestatp[9 ]= 2.8381 ;
  elestatp[10]= 5.9416 ;
  elestatp[11]= 10.0630;
  elestatp[12]= 14.1571;
  elestatp[13]= 24.4459;
  
  double elestatm[14] = {0};
  elestatm[1 ] = 0.9535 ;
  elestatm[2 ] = 1.1501 ;
  elestatm[3 ] = 1.3974 ;
  elestatm[4 ] = 1.4909 ;
  elestatm[5 ] = 2.0097 ;
  elestatm[6 ] = 2.0297 ;
  elestatm[7 ] = 2.6665 ;
  elestatm[8 ] = 2.6099 ;
  elestatm[9 ] = 3.3108 ;
  elestatm[10] = 6.8110 ;
  elestatm[11] = 12.2189;
  elestatm[12] = 17.6516;
  elestatm[13] = 25.4084;

  double elebinp[14]={0};
  elebinp[1 ]= 0.03877;
  elebinp[2 ]= 0.01183;
  elebinp[3 ]= 0.02302;
  elebinp[4 ]= 0.04825;
  elebinp[5 ]= 0.05766;
  elebinp[6 ]= 0.05428;
  elebinp[7 ]= 0.04226;
  elebinp[8 ]= 0.02474;
  elebinp[9 ]= 0.00508;
  elebinp[10]= 0.01294;
  elebinp[11]= 0.02730;
  elebinp[12]= 0.03714;
  elebinp[13]= 0.04210;

  double elebinm[14] ={0};
  elebinm[1 ]= 0.25280; 
  elebinm[2 ]= 0.19785;
  elebinm[3 ]= 0.06880;
  elebinm[4 ]= 0.11966;
  elebinm[5 ]= 0.31634;
  elebinm[6 ]= 0.48077;
  elebinm[7 ]= 0.60386;
  elebinm[8 ]= 0.69643;
  elebinm[9 ]= 0.76928;
  elebinm[10]= 0.82633;
  elebinm[11]= 0.86829;
  elebinm[12]= 0.89576;
  elebinm[13]= 0.90930;
 
  double elesigp[14] = {0};
  elesigp[1 ]= 0.39865;
  elesigp[2 ]= 0.33049;
  elesigp[3 ]= 0.15210;
  elesigp[4 ]= 0.13138;
  elesigp[5 ]= 0.46286;
  elesigp[6 ]= 0.78814;
  elesigp[7 ]= 1.08308;
  elesigp[8 ]= 1.34320;
  elesigp[9 ]= 1.56515;
  elesigp[10]= 1.74286;
  elesigp[11]= 1.87481;
  elesigp[12]= 1.96179;
  elesigp[13]= 2.00490;

  double elesigm[14] ={0};
  elesigm[1 ]= 0.37072;
  elesigm[2 ]= 0.36742;
  elesigm[3 ]= 0.28057;
  elesigm[4 ]= 0.04663;
  elesigm[5 ]= 0.32093;
  elesigm[6 ]= 0.74984;
  elesigm[7 ]= 1.17463;
  elesigm[8 ]= 1.55596;
  elesigm[9 ]= 1.87348; 
  elesigm[10]= 2.11925;
  elesigm[11]= 2.29682;
  elesigm[12]= 2.41168;
  elesigm[13]= 2.46798;
  
  double elebckgrp[14] = {0};
  elebckgrp[1 ] = 0.04517;
  elebckgrp[2 ] = 0.05155;
  elebckgrp[3 ] = 0.04341;
  elebckgrp[4 ] = 0.00934;
  elebckgrp[5 ] = 0.04491;
  elebckgrp[6 ] = 0.10735;
  elebckgrp[7 ] = 0.17016;
  elebckgrp[8 ] = 0.22973;
  elebckgrp[9 ] = 0.28301;
  elebckgrp[10] = 0.32680;
  elebckgrp[11] = 0.35980;
  elebckgrp[12] = 0.38175;
  elebckgrp[13] = 0.39267;

  double elebckgrm[14] = {0};
  elebckgrm[1 ] = 0.14033;
  elebckgrm[2 ] = 0.08632; 
  elebckgrm[3 ] = 0.00934;
  elebckgrm[4 ] = 0.11442;
  elebckgrm[5 ] = 0.19382;
  elebckgrm[6 ] = 0.23734;
  elebckgrm[7 ] = 0.25406;
  elebckgrm[8 ] = 0.25822;
  elebckgrm[9 ] = 0.25960;
  elebckgrm[10] = 0.26125;
  elebckgrm[11] = 0.26291;
  elebckgrm[12] = 0.26418;
  elebckgrm[13] = 0.26485;

  double eletoyp[14] = {0};
  eletoyp[1 ] = 0.09991;
  eletoyp[2 ] = 0.05228;
  eletoyp[3 ] = 0.08039;
  eletoyp[4 ] = 0.13596;
  eletoyp[5 ] = 0.15650;
  eletoyp[6 ] = 0.15925;
  eletoyp[7 ] = 0.16929;
  eletoyp[8 ] = 0.19459;
  eletoyp[9 ] = 0.22651;
  eletoyp[10] = 0.25662;
  eletoyp[11] = 0.28079;
  eletoyp[12] = 0.29735;
  eletoyp[13] = 0.30573;

  double eletoym[14] = {0};
  eletoym[1 ] = 0.11073;
  eletoym[2 ] = 0.05699;
  eletoym[3 ] = 0.08784;
  eletoym[4 ] = 0.15023;
  eletoym[5 ] = 0.17076;
  eletoym[6 ] = 0.16855;
  eletoym[7 ] = 0.16933;
  eletoym[8 ] = 0.18348;
  eletoym[9 ] = 0.20621;
  eletoym[10] = 0.22998;
  eletoym[11] = 0.24993;
  eletoym[12] = 0.26388;
  eletoym[13] = 0.27100;

  double eletotaleffp[14]={0};
  double eletotaleffm[14]={0};
  for(int i(1);i<14;i++)
  {
    eletotaleffp[i]  = sqrt(elebinp[i] *elebinp[i] + elesigp[i]*elesigp[i] + elebckgrp[i]*elebckgrp[i] + eletoyp[i] *eletoyp[i] );
    eletotaleffm[i]  = sqrt(elebinm[i] *elebinm[i] + elesigm[i]*elesigm[i] + elebckgrm[i]*elebckgrm[i] + eletoym[i] *eletoym[i] );
  } 
  
  double elemetp[14] = {0};
  elemetp[1 ] = 0.03088258; 
  elemetp[2 ] = 0.01641456;
  elemetp[3 ] = 0.02595842;
  elemetp[4 ] = 0.04142589;
  elemetp[5 ] = 0.04663606;
  elemetp[6 ] = 0.04826601;
  elemetp[7 ] = 0.05480307;
  elemetp[8 ] = 0.06866801;
  elemetp[9 ] = 0.08643269;
  elemetp[10] = 0.10408737;
  elemetp[11] = 0.11890443;
  elemetp[12] = 0.12938180;
  elemetp[13] = 0.13476167;

  double elemetm[14] ={0};
  elemetm[1 ] = 0.03826224;
  elemetm[2 ] = 0.02026341;
  elemetm[3 ] = 0.03113793;
  elemetm[4 ] = 0.05174631;
  elemetm[5 ] = 0.05822277;
  elemetm[6 ] = 0.05796734;
  elemetm[7 ] = 0.06263454;
  elemetm[8 ] = 0.07688163;
  elemetm[9 ] = 0.09644810;
  elemetm[10] = 0.11586899;
  elemetm[11] = 0.13192741;
  elemetm[12] = 0.14312494;
  elemetm[13] = 0.14881747;

  double elescalep[14] = {0};
  elescalep[1 ] = 0.22; 
  elescalep[2 ] = 0.11;
  elescalep[3 ] = 0.24;
  elescalep[4 ] = 0.32;
  elescalep[5 ] = 0.25;
  elescalep[6 ] = 0.23;
  elescalep[7 ] = 0.25;
  elescalep[8 ] = 0.42;
  elescalep[9 ] = 0.38;
  elescalep[10] = 0.22;
  elescalep[11] = 0.22;
  elescalep[12] = 0.30;
  elescalep[13] = 0.23;

  double elescalem[14] ={0};
  elescalem[1 ] = 0.13; 
  elescalem[2 ] = 0.07;
  elescalem[3 ] = 0.09;
  elescalem[4 ] = 0.21;
  elescalem[5 ] = 0.17;
  elescalem[6 ] = 0.17;
  elescalem[7 ] = 0.23;
  elescalem[8 ] = 0.46;
  elescalem[9 ] = 0.55;
  elescalem[10] = 0.60;
  elescalem[11] = 0.57;
  elescalem[12] = 1.17;
  elescalem[13] = 0.46;

  double elesmearp[14] ={0};
  elesmearp[1 ] = 0.1732; 
  elesmearp[2 ] = 0.0974;
  elesmearp[3 ] = 0.2148;
  elesmearp[4 ] = 0.2825;
  elesmearp[5 ] = 0.3039;
  elesmearp[6 ] = 0.2855;
  elesmearp[7 ] = 0.3677;
  elesmearp[8 ] = 0.5378;
  elesmearp[9 ] = 1.0118;
  elesmearp[10] = 0.9671;
  elesmearp[11] = 1.3643;
  elesmearp[12] = 1.5346;
  elesmearp[13] = 1.2907;

  double elesmearm[14] = {0};
  elesmearm[1 ] = 0.1992; 
  elesmearm[2 ] = 0.0922;
  elesmearm[3 ] = 0.2937;
  elesmearm[4 ] = 0.4304;
  elesmearm[5 ] = 0.4309;
  elesmearm[6 ] = 0.4920;
  elesmearm[7 ] = 0.3853;
  elesmearm[8 ] = 0.5217;
  elesmearm[9 ] = 1.3134;
  elesmearm[10] = 0.8498;
  elesmearm[11] = 1.2776;
  elesmearm[12] = 2.0242;
  elesmearm[13] = 1.5236;

  double eleEnResp[14]={0};
  double eleEnResm[14]={0};
  for(int i(1);i<14;i++)
  {
    eleEnResp[i]  = sqrt(elesmearp[i] *elesmearp[i] +elescalep[i] *elescalep[i] ); 
    eleEnResm[i]  = sqrt(elesmearm[i] *elesmearm[i] +elescalem[i] *elescalem[i] ); 
  }
  
  double eleqcdbckgrp[14]={0};
  eleqcdbckgrp[1 ] = 0.7566;
  eleqcdbckgrp[2 ] = 0.9815;
  eleqcdbckgrp[3 ] = 0.7080;
  eleqcdbckgrp[4 ] = 0.9895;
  eleqcdbckgrp[5 ] = 1.1917;
  eleqcdbckgrp[6 ] = 1.8995;
  eleqcdbckgrp[7 ] = 1.2961;
  eleqcdbckgrp[8 ] = 2.6422;
  eleqcdbckgrp[9 ] = 1.2377; 
  eleqcdbckgrp[10] = 1.7518;
  eleqcdbckgrp[11] = 1.8791;
  eleqcdbckgrp[12] = 1.7617;
  eleqcdbckgrp[13] = 1.9678;
 
  double eleqcdbckgrm[14]={0};
  eleqcdbckgrm[1 ] = 0.5666;
  eleqcdbckgrm[2 ] = 0.6639;
  eleqcdbckgrm[3 ] = 0.5736;
  eleqcdbckgrm[4 ] = 0.7592;
  eleqcdbckgrm[5 ] = 0.9046;
  eleqcdbckgrm[6 ] = 1.4389;
  eleqcdbckgrm[7 ] = 1.0041;
  eleqcdbckgrm[8 ] = 1.9908;
  eleqcdbckgrm[9 ] = 0.8273; 
  eleqcdbckgrm[10] = 2.1828;
  eleqcdbckgrm[11] = 1.4805;
  eleqcdbckgrm[12] = 1.4136;
  eleqcdbckgrm[13] = 3.4543;

  double eleqcdshapep[14] = {0};
  eleqcdshapep[1 ] = 0.2751;
  eleqcdshapep[2 ] = 0.4092;
  eleqcdshapep[3 ] = 0.5718;
  eleqcdshapep[4 ] = 0.4207;
  eleqcdshapep[5 ] = 0.6810;
  eleqcdshapep[6 ] = 0.5946;
  eleqcdshapep[7 ] = 0.6205;
  eleqcdshapep[8 ] = 0.8766;
  eleqcdshapep[9 ] = 0.8689; 
  eleqcdshapep[10] = 0.9719;
  eleqcdshapep[11] = 0.4631;
  eleqcdshapep[12] = 0.8889;
  eleqcdshapep[13] = 0.9202;
 
  double eleqcdshapem[14] ={0};
  eleqcdshapem[1 ] = 0.2662;
  eleqcdshapem[2 ] = 0.2565;
  eleqcdshapem[3 ] = 0.3511;
  eleqcdshapem[4 ] = 0.8790;
  eleqcdshapem[5 ] = 0.7882;
  eleqcdshapem[6 ] = 0.4798;
  eleqcdshapem[7 ] = 0.6443;
  eleqcdshapem[8 ] = 0.6762;
  eleqcdshapem[9 ] = 0.8812; 
  eleqcdshapem[10] = 0.6523;
  eleqcdshapem[11] = 0.9204;
  eleqcdshapem[12] = 0.7787;
  eleqcdshapem[13] = 0.8676;

  double eleewkp[14] ={0};
  eleewkp[1 ] = 0.03630;
  eleewkp[2 ] = 0.03531;
  eleewkp[3 ] = 0.02500;
  eleewkp[4 ] = 0.02548;
  eleewkp[5 ] = 0.03906;
  eleewkp[6 ] = 0.07970;
  eleewkp[7 ] = 0.12190;
  eleewkp[8 ] = 0.15987;
  eleewkp[9 ] = 0.19155; 
  eleewkp[10] = 0.21635;
  eleewkp[11] = 0.23452;
  eleewkp[12] = 0.24640;
  eleewkp[13] = 0.25226;
 
  double eleewkm[14] ={0};
  eleewkm[1 ] = 0.06766;
  eleewkm[2 ] = 0.03786;
  eleewkm[3 ] = 0.02815;
  eleewkm[4 ] = 0.06349;
  eleewkm[5 ] = 0.07465;
  eleewkm[6 ] = 0.09850;
  eleewkm[7 ] = 0.12432;
  eleewkm[8 ] = 0.14779;
  eleewkm[9 ] = 0.16639; 
  eleewkm[10] = 0.17964;
  eleewkm[11] = 0.18861;
  eleewkm[12] = 0.19418;
  eleewkm[13] = 0.19686;

  double elefsrp[14] ={0};
  elefsrp[1 ] = 0.13016;
  elefsrp[2 ] = 0.12063;
  elefsrp[3 ] = 0.08499;
  elefsrp[4 ] = 0.00678;
  elefsrp[5 ] = 0.12112;
  elefsrp[6 ] = 0.27081;
  elefsrp[7 ] = 0.43071;
  elefsrp[8 ] = 0.57278;
  elefsrp[9 ] = 0.68939;
  elefsrp[10] = 0.78659;
  elefsrp[11] = 0.87403;
  elefsrp[12] = 0.88221;
  elefsrp[13] = 0.92890;
 
  double elefsrm[14] ={0};
  elefsrm[1 ] = 0.07837;
  elefsrm[2 ] = 0.04943;
  elefsrm[3 ] = 0.02817;
  elefsrm[4 ] = 0.00768;
  elefsrm[5 ] = 0.03028;
  elefsrm[6 ] = 0.10332;
  elefsrm[7 ] = 0.20393;
  elefsrm[8 ] = 0.30426;
  elefsrm[9 ] = 0.38436;
  elefsrm[10] = 0.45585;
  elefsrm[11] = 0.50665;
  elefsrm[12] = 0.55795;
  elefsrm[13] = 0.65364;

  double elesvdunfp[14] ={0};
  elesvdunfp[1 ] = 0.1153;
  elesvdunfp[2 ] = 0.1043;
  elesvdunfp[3 ] = 0.1063;
  elesvdunfp[4 ] = 0.1227;
  elesvdunfp[5 ] = 0.1418;
  elesvdunfp[6 ] = 0.1621;
  elesvdunfp[7 ] = 0.1875;
  elesvdunfp[8 ] = 0.2223;
  elesvdunfp[9 ] = 0.2622;
  elesvdunfp[10] = 0.3000;
  elesvdunfp[11] = 0.3312;
  elesvdunfp[12] = 0.3531;
  elesvdunfp[13] = 0.3644;
  
  double elesvdunfm[14] ={0};
  elesvdunfm[1 ] = 0.1196;
  elesvdunfm[2 ] = 0.1094;
  elesvdunfm[3 ] = 0.1138;
  elesvdunfm[4 ] = 0.1303;
  elesvdunfm[5 ] = 0.1481;
  elesvdunfm[6 ] = 0.1684;
  elesvdunfm[7 ] = 0.1979;
  elesvdunfm[8 ] = 0.2364;
  elesvdunfm[9 ] = 0.2783;
  elesvdunfm[10] = 0.3174;
  elesvdunfm[11] = 0.3489;
  elesvdunfm[12] = 0.3706;
  elesvdunfm[13] = 0.3816;

  double eleBiasp[14] = {0};
  eleBiasp[1 ] = 0.73817 ;
  eleBiasp[2 ] = 1.49013 ;
  eleBiasp[3 ] = 1.34147 ;
  eleBiasp[4 ] = 0.34694 ;
  eleBiasp[5 ] = 0.40953 ;
  eleBiasp[6 ] = 0.17184 ;
  eleBiasp[7 ] = 0.30586 ;
  eleBiasp[8 ] = 0.26019 ;
  eleBiasp[9 ] = 2.07064 ;
  eleBiasp[10] = 1.64592 ;
  eleBiasp[11] = 5.16036 ;
  eleBiasp[12] = 2.27049 ;
  eleBiasp[13] = 4.99595 ;
 
  double eleBiasm[14]={0};
  eleBiasm[1 ] = 0.77536 ;
  eleBiasm[2 ] = 1.34162 ;
  eleBiasm[3 ] = 0.77528 ;
  eleBiasm[4 ] = 0.39036 ;
  eleBiasm[5 ] = 0.82576 ;
  eleBiasm[6 ] = 0.46733 ;
  eleBiasm[7 ] = 0.61969 ;
  eleBiasm[8 ] = 0.74028 ;
  eleBiasm[9 ] = 2.59330 ;
  eleBiasm[10] = 3.17404 ;
  eleBiasm[11] = 3.77210 ;
  eleBiasm[12] = 3.91817 ;
  eleBiasm[13] = 2.53789 ;
/*
  double eleWptCorrp[14] = {0};
  eleWptCorrp[1 ] = 0.915318; 
  eleWptCorrp[2 ] = 0.0294394;
  eleWptCorrp[3 ] = 0.793438; 
  eleWptCorrp[4 ] = 1.00581;  
  eleWptCorrp[5 ] = 0.659224; 
  eleWptCorrp[6 ] = 0.121046; 
  eleWptCorrp[7 ] = 0.299264; 
  eleWptCorrp[8 ] = 0.525082; 
  eleWptCorrp[9 ] = 0.592397; 
  eleWptCorrp[10] = 0.57529;  
  eleWptCorrp[11] = 0.529585; 
  eleWptCorrp[12] = 0.486141; 
  eleWptCorrp[13] = 0.461337; 

  double eleWptCorrm[14] ={0};
  eleWptCorrm[1 ] = 1.21547; 
  eleWptCorrm[2 ] = 0.262792;
  eleWptCorrm[3 ] = 0.743959;
  eleWptCorrm[4 ] = 1.17545; 
  eleWptCorrm[5 ] = 0.937213;
  eleWptCorrm[6 ] = 0.336006;
  eleWptCorrm[7 ] = 0.252154;
  eleWptCorrm[8 ] = 0.603749;
  eleWptCorrm[9 ] = 0.671394;
  eleWptCorrm[10] = 0.576984;
  eleWptCorrm[11] = 0.437997;
  eleWptCorrm[12] = 0.319664;
  eleWptCorrm[13] = 0.255065;

  for(int i(1);i<14;i++)
  {
    eleBiasp[i] = sqrt(eleBiasp[i] * eleBiasp[i] + eleWptCorrp[i] * eleWptCorrp[i]) ; 
    eleBiasm[i] = sqrt(eleBiasm[i] * eleBiasm[i] + eleWptCorrm[i] * eleWptCorrm[i]) ; 
  }
 */ 

  double systtotalp[14] = {0};
  double systtotalm[14] = {0};
  for(int i(1);i<14;i++)
  {
    systtotalp[i] =sqrt(eletotaleffp[i] *eletotaleffp[i] +elemetp[i] *elemetp[i] +eleEnResp[i] *eleEnResp[i] +eleqcdshapep[i] *eleqcdshapep[i] +eleqcdbckgrp[i] *eleqcdbckgrp[i] + eleewkp[i] *eleewkp[i]  +eleBiasp[i] *eleBiasp[i] +elesvdunfp[i] *elesvdunfp[i]  +elefsrp[i] *elefsrp[i] + 2.6*2.6);
    systtotalm[i] =sqrt(eletotaleffm[i] *eletotaleffm[i] +elemetm[i] *elemetm[i] +eleEnResm[i] *eleEnResm[i] +eleqcdshapem[i] *eleqcdshapem[i] +eleqcdbckgrm[i] *eleqcdbckgrm[i] + eleewkm[i] *eleewkm[i]  +eleBiasm[i] *eleBiasm[i] +elesvdunfm[i] *elesvdunfm[i]  +elefsrm[i] *elefsrm[i] + 2.6*2.6);
  }

  double totaluncerp[14]={0};
  double totaluncerm[14]={0};
  for(int i(1);i<14;i++)
  {
    totaluncerp[i] = sqrt(elestatp[i] *elestatp[i] +systtotalp[i] *systtotalp[i]);
    totaluncerm[i] = sqrt(elestatm[i] *elestatm[i] +systtotalm[i] *systtotalm[i]);
  }
 
//  Print out result
  cout<<fixed<<setprecision(2);
  
  // Recon. 
  cout<<"Recon. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<< eletotaleffp[i] << "\t" << eletotaleffm[i] <<endl;
  }
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eletotaleffp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eletotaleffm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;
  
  //Energy Resolution
  cout<<"EnRes.Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<eleEnResp[i] <<"\t"<< eleEnResm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleEnResp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleEnResm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //MET Resolution
  cout<<"METRes.Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<elemetp[i] <<"\t"<< elemetm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elemetp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elemetm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //QCD Background
  cout<<"QCD Bckgr.Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<eleqcdbckgrp[i] <<"\t"<< eleqcdbckgrm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleqcdbckgrp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleqcdbckgrm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //QCD Shape
  cout<<"QCD Shape. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<eleqcdshapep[i] <<"\t"<< eleqcdshapem[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleqcdshapep[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleqcdshapem[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //EWK
  cout<<"EWK. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<eleewkp[i] <<"\t"<< eleewkm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleewkp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleewkm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //SVD.Unf
  cout<<"SVD.Unf. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<elesvdunfp[i] <<"\t"<< elesvdunfm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elesvdunfp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elesvdunfm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //FSR
  cout<<"FSR. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<elefsrp[i] <<"\t"<< elefsrm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elefsrp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elefsrm[i] ;
  }
    cout << " \\"<<"\\" << endl;
    cout << endl;

  //Unf.Bias
  cout<<"Unf.Bias. Syst"<<endl;
  cout<<"W+ \t W-"<<endl;
  for(int i(1);i<14;i++)
  {
    cout<<eleBiasp[i] <<"\t"<< eleBiasm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleBiasp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << eleBiasm[i] ;
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
    cout<<elestatp[i] <<"\t"<< elestatm[i]<<endl; 
  } 
  cout << "======= Put in note=======" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elestatp[i] ;
  }
    cout << " \\"<<"\\" << endl;
  for(int i(1);i<14;i++)
  {
    cout<< " & " << elestatm[i] ;
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

  TFile f_out(resultDir+"/WpToEleNuErrors.root","recreate");
    
    ///Write all errors to root
    
    TH1D* h_bin = new TH1D("h_bin","h_bin",13,0,13);
    //TH1D* h_elebinm = new TH1D("h_elebinm","h_elebinm",13,0,13);
    
    TH1D* h_idisosig = new TH1D("h_idisosig","h_idisosig",13,0,13);
    //TH1D* h_sig = new TH1D("h_sig","h_sig",13,0,13);
    //TH1D* h_elesigm = new TH1D("h_elesigm","h_elesigm",13,0,13);

    TH1D* h_idisobck = new TH1D("h_idisobck","h_idisobck",13,0,13);
    //TH1D* h_bckgr = new TH1D("h_bckgr","h_bckgr",13,0,13);
    //TH1D* h_elebckgrm = new TH1D("h_elebckgrm","h_elebckgrm",13,0,13);
  
    TH1D* h_toy = new TH1D("h_toy","h_toy",13,0,13);
    //TH1D* h_eletoym = new TH1D("h_eletoym","h_eletoym",13,0,13);
    
    TH1D* h_TotalEff = new TH1D("h_Totaleff","h_TotalEff",13,0,13);
    //TH1D* h_eleTotalEffm = new TH1D("h_eleTotaleffm","h_eleTotalEffm",13,0,13);
    
    TH1D* h_met = new TH1D("h_met","h_met",13,0,13);
    //TH1D* h_elemetm = new TH1D("h_elemetm","h_elemetm",13,0,13);
    
    TH1D* h_scale = new TH1D("h_scale","h_scale",13,0,13);
    //TH1D* h_elescalem = new TH1D("h_elescalem","h_elescalem",13,0,13);
    
    TH1D* h_smear = new TH1D("h_smear","h_smear",13,0,13);
    //TH1D* h_elesmearm = new TH1D("h_elesmearm","h_elesmearm",13,0,13);
    
    TH1D* h_EnMomRes = new TH1D("h_EnMomRes","h_EnMomRes",13,0,13);
    //TH1D* h_eleEnResm = new TH1D("h_eleEnResm","h_eleEnResm",13,0,13);
    
    TH1D* h_qcdbckgr = new TH1D("h_qcdbckgr","h_qcdbckgr",13,0,13);
    //TH1D* h_eleqcdbckgrm = new TH1D("h_eleqcdbckgrm","h_eleqcdbckgrm",13,0,13);
    
    TH1D* h_qcdshape = new TH1D("h_qcdshape","h_qcdshape",13,0,13);
    //TH1D* h_eleqcdshapem = new TH1D("h_eleqcdshapem","h_eleqcdshapem",13,0,13);
    
    TH1D* h_ewk = new TH1D("h_ewk","h_ewk",13,0,13);
    //TH1D* h_eleewkm = new TH1D("h_eleewkm","h_eleewkm",13,0,13);
    
    TH1D* h_fsr = new TH1D("h_fsr","h_fsr",13,0,13);
    //TH1D* h_elefsrm = new TH1D("h_elefsrm","h_elefsrm",13,0,13);
    
    TH1D* h_SvdUnf = new TH1D("h_SvdUnf","h_SvdUnf",13,0,13);
    //TH1D* h_eleSvdUnfm = new TH1D("h_eleSvdUnfm","h_eleSvdUnfm",13,0,13);
    
    TH1D* h_UnfoldBias = new TH1D("h_UnfoldBias","h_UnfoldBias",13,0,13);
    //TH1D* h_eleUnfoldBiasm = new TH1D("h_eleUnfoldBiasm","h_eleUnfoldBiasm",13,0,13);
    
    TH1D* h_TotalSyst = new TH1D("h_TotalSyst","h_TotalSyst",13,0,13);
    //TH1D* h_eleTotalSystm = new TH1D("h_eleTotalSystm","h_eleTotalSystm",13,0,13);
    
    TH1D* h_Stat = new TH1D("h_Stat","h_Stat",13,0,13);
    //TH1D* h_eleStatm = new TH1D("h_eleStatm","h_eleStatm",13,0,13);
    
    TH1D* h_TotalUncer = new TH1D("h_TotalUncer","h_TotalUncer",13,0,13);
    //TH1D* h_eleTotalUncerm = new TH1D("h_eleTotalUncerm","h_eleTotalUncerm",13,0,13);
    
    TH1D* h_PowhegPDF = new TH1D("h_PowhegPDF","h_PowhegPDF",13,0,13);
    //TH1D* h_elePowhegPDFm = new TH1D("h_elePowhegPDFm","h_elePowhegPDFm",13,0,13);

    TH1D* h_LumiSyst = new TH1D("h_LumiSyst","h_LumiSyst",13,0,13);
    //TH1D* h_eleLumiSystm = new TH1D("h_eleLumiSystm","h_eleLumiSystm",13,0,13);

    double PowhegPDFp[14];
    PowhegPDFp[0] = 0;
    PowhegPDFp[1] = 4.285;
    PowhegPDFp[2] = 4.156;
    PowhegPDFp[3] = 4.132 ;
    PowhegPDFp[4] = 5.839;
    PowhegPDFp[5] = 4.127;
    PowhegPDFp[6] = 5.549;
    PowhegPDFp[7] = 3.104;
    PowhegPDFp[8] = 4.887;
    PowhegPDFp[9] = 2.932;
    PowhegPDFp[10]= 2.299;
    PowhegPDFp[11]= 2.881;
    PowhegPDFp[12]= 2.793;
    PowhegPDFp[13]= 5.133;


    double PowhegPDFm[14];
    PowhegPDFm[0] = 0;
    PowhegPDFm[1] = 4.369;
    PowhegPDFm[2] = 4.373;
    PowhegPDFm[3] = 4.392;
    PowhegPDFm[4] = 4.416;
    PowhegPDFm[5] = 4.475;
    PowhegPDFm[6] = 4.811;
    PowhegPDFm[7] = 4.394;
    PowhegPDFm[8] = 4.216;
    PowhegPDFm[9] = 4.728 ;
    PowhegPDFm[10]= 4.818;
    PowhegPDFm[11]= 4.679;
    PowhegPDFm[12]= 4.765;
    PowhegPDFm[13]= 5.371;
    

    for(int i(1);i<14;i++)
    {
      h_Stat->SetBinContent(i,elestatp[i]);
      //h_eleStatm->SetBinContent(i,elestatm[i]);
      
      h_bin->SetBinContent(i,elebinp[i]);
      //h_elebinm->SetBinContent(i,elebinm[i]);
      
      h_idisosig->SetBinContent(i,elesigp[i]);
      //h_elesigm->SetBinContent(i,elesigm[i]);
      
      h_idisobck->SetBinContent(i,elebckgrp[i]);
      //h_bckgr->SetBinContent(i,elebckgrp[i]);
      //h_elebckgrm->SetBinContent(i,elebckgrm[i]);
   
      h_toy->SetBinContent(i,eletoyp[i]);
      //h_eletoym->SetBinContent(i,eletoym[i]);
   
      h_TotalEff->SetBinContent(i,eletotaleffp[i]);
      //h_eleTotalEffm->SetBinContent(i,eletotaleffm[i]);
   
      h_met->SetBinContent(i,elemetp[i]);
      //h_elemetm->SetBinContent(i,elemetm[i]);
   
      h_scale->SetBinContent(i,elescalep[i]);
      //h_elescalem->SetBinContent(i,elescalem[i]);
   
      h_smear->SetBinContent(i,elesmearp[i]);
      //h_elesmearm->SetBinContent(i,elesmearm[i]);
   
      h_EnMomRes->SetBinContent(i,eleEnResp[i]);
      //h_eleEnResm->SetBinContent(i,eleEnResm[i]);
   
      h_qcdbckgr->SetBinContent(i,eleqcdbckgrp[i]);
      //h_eleqcdbckgrm->SetBinContent(i,eleqcdbckgrm[i]);
   
      h_qcdshape->SetBinContent(i,eleqcdshapep[i]);
      //h_eleqcdshapem->SetBinContent(i,eleqcdshapem[i]);
   
      h_ewk->SetBinContent(i,eleewkp[i]);
      //h_eleewkm->SetBinContent(i,eleewkm[i]);
   
      h_fsr->SetBinContent(i,elefsrp[i]);
      //h_elefsrm->SetBinContent(i,elefsrm[i]);
   
      h_SvdUnf->SetBinContent(i,elesvdunfp[i]);
      //h_eleSvdUnfm->SetBinContent(i,elesvdunfm[i]);
   
      h_UnfoldBias->SetBinContent(i,eleBiasp[i]);
      //h_eleUnfoldBiasm->SetBinContent(i,eleBiasm[i]);
   
      h_TotalSyst->SetBinContent(i,systtotalp[i]);
      //h_eleTotalSystm->SetBinContent(i,systtotalm[i]);
   
      h_TotalUncer->SetBinContent(i,totaluncerp[i]);
      //h_eleTotalUncerm->SetBinContent(i,totaluncerm[i]);
   
      h_PowhegPDF->SetBinContent(i,PowhegPDFp[i]);
      //h_elePowhegPDFm->SetBinContent(i,PowhegPDFm[i]);
  
      h_LumiSyst->SetBinContent(i,2.6);
      //h_eleLumiSystm->SetBinContent(i,2.6);

    }
    h_Stat->Write();
    //h_eleStatm->Write();

    h_bin->Write();
    //h_elebinm->Write();

    h_idisosig->Write();
    //h_elesigm->Write();
    
    h_idisobck->Write();
    //h_bckgr->Write();
    //h_elebckgrm->Write();
                   
    h_toy->Write();
    //h_eletoym->Write();
                   
    h_TotalEff->Write();
    //h_eleTotalEffm->Write();
                   
    h_met->Write();
    //h_elemetm->Write();
                   
    h_scale->Write();
    //h_elescalem->Write();
                   
    h_smear->Write();
    //h_elesmearm->Write();
                   
    h_EnMomRes->Write();
    //h_eleEnResm->Write();
                   
    h_qcdbckgr->Write();
    //h_eleqcdbckgrm->Write();
                   
    h_qcdshape->Write();
    //h_eleqcdshapem->Write();
                   
    h_ewk->Write();
    //h_eleewkm->Write();
                   
    h_fsr->Write();
    //h_elefsrm->Write();
                   
    h_SvdUnf->Write();
    //h_eleSvdUnfm->Write();
                   
    h_UnfoldBias->Write();
    //h_eleUnfoldBiasm->Write();
                   
    h_TotalSyst->Write();
    //h_eleTotalSystm->Write();
                   
    h_TotalUncer->Write();
    //h_eleTotalUncerm->Write();

    h_PowhegPDF->Write();
    //h_elePowhegPDFm->Write();
   
    h_LumiSyst->Write();
    //h_eleLumiSystm->Write();


}
