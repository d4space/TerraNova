#include "../../../Utils/const.h"

double BinWidth[14] ={0.0, 7.5-0, 12.5-7.5,17.5-12.5, 24.0-17.5, 30.0-24.0, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

void CombineChargeToIncl_Mu()
{
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
 
  double muStatErri[14] = {0};
  muStatErri[1 ]= 0.51230950  ; 
  muStatErri[2 ]= 0.64972669  ;
  muStatErri[3 ]= 0.78613389  ;
  muStatErri[4 ]= 0.85383044  ;
  muStatErri[5 ]= 1.13627363  ;
  muStatErri[6 ]= 1.13678485  ;
  muStatErri[7 ]= 1.58058107  ;
  muStatErri[8 ]= 1.56710703  ;
  muStatErri[9 ]= 2.03175297  ;
  muStatErri[10]= 4.10703670  ;
  muStatErri[11]= 7.88607215  ;
  muStatErri[12]= 12.69383560 ;
  muStatErri[13]= 19.66763898 ;

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

  double muscaleErri[14]={0};
  muscaleErri[1] = 0.0833  ;
  muscaleErri[2] = 0.0539  ;
  muscaleErri[3] =  0.0791 ;
  muscaleErri[4] =  0.1508 ;
  muscaleErri[5] =   0.1703;
  muscaleErri[6] =   0.1849;
  muscaleErri[7] =   0.1506;
  muscaleErri[8] =  0.0981 ;
  muscaleErri[9] =   0.1885;
  muscaleErri[10]=   0.2679;
  muscaleErri[11]=   0.3299;
  muscaleErri[12]=   0.3693;
  muscaleErri[13]=   0.3883;

  double musmearErri[14]={0};
  musmearErri[1] = 0.0720 ;
  musmearErri[2] = 0.0331 ;
  musmearErri[3] = 0.0523 ;
  musmearErri[4] = 0.1379 ;
  musmearErri[5] = 0.1804 ;
  musmearErri[6] = 0.1553 ;
  musmearErri[7] = 0.0701 ;
  musmearErri[8] = 0.0959 ;
  musmearErri[9] = 0.2050 ;
  musmearErri[10]=  0.2930;
  musmearErri[11]=  0.3581;
  musmearErri[12]=  0.4037;
  musmearErri[13]=  0.4283;
  
  double muMomResErri[14]={0};
  muMomResErri[1] = 0.1101;
  muMomResErri[2] = 0.0633;
  muMomResErri[3] = 0.0948;
  muMomResErri[4] = 0.2043;
  muMomResErri[5] = 0.2481;
  muMomResErri[6] = 0.2415;
  muMomResErri[7] = 0.1661;
  muMomResErri[8] = 0.1372;
  muMomResErri[9] = 0.2785;
  muMomResErri[10]= 0.3970;
  muMomResErri[11]= 0.4869;
  muMomResErri[12]= 0.5471;
  muMomResErri[13]= 0.5781;

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

  double muUnfBiasErri[14] = {0};
  muUnfBiasErri[1 ] = 0.9277 ;  //0.31   ;
  muUnfBiasErri[2 ] = 1.7211 ;  //0.33   ;
  muUnfBiasErri[3 ] = 1.1504 ;  //0.40   ;
  muUnfBiasErri[4 ] = 0.3039 ;  //0.38   ;
  muUnfBiasErri[5 ] = 0.6472 ;  //0.36   ;
  muUnfBiasErri[6 ] = 0.2676 ;  //0.67   ;
  muUnfBiasErri[7 ] = 0.4361 ;  //1.20   ;
  muUnfBiasErri[8 ] = 0.7751 ;  //1.67   ;
  muUnfBiasErri[9 ] = 1.9656 ;  //1.72   ;
  muUnfBiasErri[10] = 1.9656 ;  //1.43   ;
  muUnfBiasErri[11] = 1.9656 ;  //1.72   ;
  muUnfBiasErri[12] = 5.4563 ;  //1.81   ;
  muUnfBiasErri[13] = 5.4563 ;  //2.04   ;

  // Make Incl
  //*
  TFile *fp = new TFile("../../ResultWpToMuNu/Result_WpToMuNu.root");
  TFile *fm = new TFile("../../ResultWmToMuNu/Result_WmToMuNu.root");

  TH1D *h_data_p;
  TH1D *h_data_m;
  TH1D *h_powheg_p;
  TH1D *h_powheg_m;
  TH1D *h_powheg_p_UnWeighted;
  TH1D *h_powheg_m_UnWeighted;
  h_data_p = (TH1D*)fp->Get("BornEffCorr")->Clone("h_data_p");
  h_data_m = (TH1D*)fm->Get("BornEffCorr")->Clone("h_data_m");
  h_powheg_p = (TH1D*)fp->Get("SVD_Born.Gen")->Clone("h_MC_p");
  h_powheg_m = (TH1D*)fm->Get("SVD_Born.Gen")->Clone("h_MC_m");
  h_powheg_p_UnWeighted = (TH1D*)fp->Get("SVD_Born.Gen")->Clone("h_MC_p_UnWeighted");
  h_powheg_m_UnWeighted = (TH1D*)fm->Get("SVD_Born.Gen")->Clone("h_MC_m_UnWeighted");

  // To make SVD_Born.Gen to h1_Born_BornFid(Not LumiWeighted)
  h_powheg_p_UnWeighted->Scale(1./LumiWeight_Muon_WpToMuNu_S8);
  h_powheg_m_UnWeighted->Scale(1./LumiWeight_Muon_WmToMuNu_S8);

  h_powheg_p->Scale(1./18.429);
  h_powheg_m->Scale(1./18.429);

  cout << "Inclusive Cross-section" << endl;
  cout << "bin\tW+\t\tW-\t\t Wincl"<<endl;
  for( int ipt(1);ipt<14;ipt++)
  {
    cout<<ipt<<"\t"<<h_data_p->GetBinContent(ipt)<<"\t\t"<<h_data_m->GetBinContent(ipt)<< "\t\t" << h_data_p->GetBinContent(ipt)+h_data_m->GetBinContent(ipt) << endl;
  }

    // Mu SVDUnf Error
  double muSVDUnfErrp[14]={0};
  double muSVDUnfErrm[14]={0};
  double muSVDUnfErri[14]={0};
  cout <<"W+ SVDUnfErr \t W- SVDUnfErr \t W SVDUnfErr" <<endl;
  for(int i(1);i<14;i++)
  {
    muSVDUnfErrp[i] = h_data_p->GetBinContent(i)*0.01*musvdunfp[i];
    muSVDUnfErrm[i] = h_data_m->GetBinContent(i)*0.01*musvdunfm[i];
    muSVDUnfErri[i] = (sqrt(muSVDUnfErrp[i]*muSVDUnfErrp[i] + muSVDUnfErrm[i]*muSVDUnfErrm[i])/(h_data_p->GetBinContent(i)+h_data_m->GetBinContent(i))) * 100; // % unit
    cout << muSVDUnfErrp[i] << "\t" << muSVDUnfErrm[i] <<"\t" << muSVDUnfErri[i] <<endl;
  }

  // Save files
  TString resultDir = "ResultWInclToMuNu";
  gSystem->mkdir(resultDir,kTRUE);

  TFile f_out(resultDir+"/Result_WInclToMuNu_DataPowheg.root","recreate");
 
  TH1D* BornEffCorr = new TH1D("BornEffCorr","BornEffCorr13",13,0,13);
  TH1D* SVD_BornGen = new TH1D("SVD_Born.Gen","SVD_Born.Gen13",13,0,13);
  TH1D* SVD_BornGen_UnWeighted = new TH1D("SVD_Born.Gen_UnWeighted","SVD_Born.Gen13_UnWeighted",13,0,13);

  double PowhegStatErrP[14] ={0.,};
  double PowhegStatErrM[14] ={0.,};
  double PowhegStatErrI[14] ={0.,};
  for(int i(1); i<14; i++)
  {
    BornEffCorr->SetBinContent(i,h_data_p->GetBinContent(i) + h_data_m->GetBinContent(i));
    SVD_BornGen->SetBinContent(i,h_powheg_p->GetBinContent(i) + h_powheg_m->GetBinContent(i));
    SVD_BornGen_UnWeighted->SetBinContent(i,h_powheg_p_UnWeighted->GetBinContent(i) + h_powheg_m_UnWeighted->GetBinContent(i));
    PowhegStatErrP[i] = h_powheg_p->GetBinContent(i) * sqrt(h_powheg_p_UnWeighted->GetBinContent(i))/h_powheg_p_UnWeighted->GetBinContent(i);
    PowhegStatErrM[i] = h_powheg_m->GetBinContent(i) * sqrt(h_powheg_m_UnWeighted->GetBinContent(i))/h_powheg_m_UnWeighted->GetBinContent(i);
    PowhegStatErrI[i] = sqrt(PowhegStatErrP[i]*PowhegStatErrP[i] + PowhegStatErrM[i]*PowhegStatErrM[i]);

    cout << Form("SVD_BornGen : %.2f, StatErr : %.4f",SVD_BornGen->GetBinContent(i),PowhegStatErrI[i]) << endl;
  }
  BornEffCorr->Write();
  SVD_BornGen->Write();


    ///Write all errors to root

  
    TH1D* h_Stat = new TH1D("h_Stat","h_Stat",13,0,13);
    
    TH1D* h_tracksig = new TH1D("h_tracksig","h_tracksig",13,0,13);
    TH1D* h_trackbck = new TH1D("h_trackbck","h_trackbck",13,0,13);
    
    TH1D* h_idisosig = new TH1D("h_idisosig","h_idisosig",13,0,13);
    TH1D* h_idisobck = new TH1D("h_idisobck","h_idisobck",13,0,13);
  
    TH1D* h_track = new TH1D("h_track","h_track",13,0,13);
   
    TH1D* h_idiso = new TH1D("h_idiso","h_idiso",13,0,13);
    
    TH1D* h_toy = new TH1D("h_toy","h_toy",13,0,13);
    
    TH1D* h_POG = new TH1D("h_POG","h_POG",13,0,13);
    
    TH1D* h_TotalEff = new TH1D("h_Totaleff","h_TotalEff",13,0,13);
    
    TH1D* h_met = new TH1D("h_met","h_met",13,0,13);
    
    TH1D* h_scale = new TH1D("h_scale","h_scale",13,0,13);
    
    TH1D* h_smear = new TH1D("h_smear","h_smear",13,0,13);
    
    TH1D* h_EnMomRes = new TH1D("h_EnMomRes","h_EnMomRes",13,0,13);
    
    TH1D* h_qcdbckgr = new TH1D("h_qcdbckgr","h_qcdbckgr",13,0,13);
    
    TH1D* h_qcdshape = new TH1D("h_qcdshape","h_qcdshape",13,0,13);
    
    TH1D* h_ewk = new TH1D("h_ewk","h_ewk",13,0,13);
    
    TH1D* h_fsr = new TH1D("h_fsr","h_fsr",13,0,13);
    
    TH1D* h_SvdUnf = new TH1D("h_SvdUnf","h_SvdUnf",13,0,13);
    
    TH1D* h_UnfoldBias = new TH1D("h_UnfoldBias","h_UnfoldBias",13,0,13);
    
    TH1D* h_TotalSyst = new TH1D("h_TotalSyst","h_TotalSyst",13,0,13);
    
    TH1D* h_TotalUncer = new TH1D("h_TotalUncer","h_TotalUncer",13,0,13);
    
    TH1D* h_PowhegStat = new TH1D("h_PowhegStat","h_PowhegStat",13,0,13);
    
    for(int i(1);i<14;i++)
    {
      h_Stat->SetBinContent(i,muStatErri[i]);
      h_tracksig->SetBinContent(i,mutracksigErri[i]);
      h_trackbck->SetBinContent(i,mutrackbckErri[i]);
      h_idisosig->SetBinContent(i,muidisosigErri[i]);
      h_idisobck->SetBinContent(i,muidisobckErri[i]);
      h_toy->SetBinContent(i,mutoyErri[i]);
      h_POG->SetBinContent(i,muPOGErri[i]);
      h_met->SetBinContent(i,muMetErri[i]);
      h_scale->SetBinContent(i,muscaleErri[i]);
      h_smear->SetBinContent(i,musmearErri[i]);
      h_qcdbckgr->SetBinContent(i,muQCDBckErri[i]);
      h_qcdshape->SetBinContent(i,muQCDShapeErri[i]);
      h_ewk->SetBinContent(i,muEWKErri[i]);
      h_fsr->SetBinContent(i,muFSRErri[i]);
      h_SvdUnf->SetBinContent(i,muSVDUnfErri[i]);
      h_UnfoldBias->SetBinContent(i,muUnfBiasErri[i]);
      h_PowhegStat->SetBinContent(i,PowhegStatErrI[i]);
    }

    h_Stat->Write();

    h_tracksig->Write();
    h_trackbck->Write();

    h_idisosig->Write();
    h_idisobck->Write();
    
    h_toy->Write();
    
    h_POG->Write();
    
    h_met->Write();
                   
    h_scale->Write();
                   
    h_smear->Write();
                   
    h_qcdbckgr->Write();
                   
    h_qcdshape->Write();
                   
    h_ewk->Write();
                   
    h_fsr->Write();
                   
    h_SvdUnf->Write();
                   
    h_UnfoldBias->Write();

    h_PowhegStat->Write();
                   
}
