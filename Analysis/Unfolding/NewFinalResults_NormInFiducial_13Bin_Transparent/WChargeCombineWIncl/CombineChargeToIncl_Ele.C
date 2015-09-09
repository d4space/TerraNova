#include "../../../Utils/const.h"

double BinWidth[14] ={0.0, 7.5-0, 12.5-7.5,17.5-12.5, 24.0-17.5, 30.0-24.0, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

void CombineChargeToIncl_Ele()
{
//*

 
  double eleidisosigErri[14] = {0};
  eleidisosigErri[1 ] =0.2802 ;
  eleidisosigErri[2 ] =0.2430 ;
  eleidisosigErri[3 ] =0.1419 ;
  eleidisosigErri[4 ] =0.0796 ;
  eleidisosigErri[5 ] =0.3059 ;
  eleidisosigErri[6 ] =0.5594 ;
  eleidisosigErri[7 ] =0.8011 ;
  eleidisosigErri[8 ] =1.0190 ;
  eleidisosigErri[9 ] =1.2036 ;
  eleidisosigErri[10] =1.3475 ;
  eleidisosigErri[11] =1.4520 ;
  eleidisosigErri[12] =1.5189 ;
  eleidisosigErri[13] =1.5429 ;


  double eleidisobckErri[14]={0};
  eleidisobckErri[1 ] = 0.0590;
  eleidisobckErri[2 ] = 0.0460;
  eleidisobckErri[3 ] = 0.0262;
  eleidisobckErri[4 ] = 0.0474;
  eleidisobckErri[5 ] = 0.0827;
  eleidisobckErri[6 ] = 0.1130;
  eleidisobckErri[7 ] = 0.1391;
  eleidisobckErri[8 ] = 0.1601;
  eleidisobckErri[9 ] = 0.1822;
  eleidisobckErri[10] = 0.2065;
  eleidisobckErri[11] = 0.2278;
  eleidisobckErri[12] = 0.2456;
  eleidisobckErri[13] = 0.2815;

  double eletoyErri[14] = {0};
  eletoyErri[1 ] =0.075997697 ;
  eletoyErri[2 ] =0.040651814 ;
  eletoyErri[3 ] =0.07560172  ;
  eletoyErri[4 ] =0.122198036 ;
  eletoyErri[5 ] =0.129290564 ;
  eletoyErri[6 ] =0.117781238 ;
  eletoyErri[7 ] =0.13797652  ;
  eletoyErri[8 ] =0.1936125   ;
  eletoyErri[9 ] =0.252288248 ;
  eletoyErri[10] =0.300506656 ;
  eletoyErri[11] =0.336240688 ;
  eletoyErri[12] =0.359969013 ;
  eletoyErri[13] =0.376577495 ;

  double eleBinningErri[14] = {0};
  eleBinningErri[1 ] =0.1027 ;
  eleBinningErri[2 ] =0.0801 ;
  eleBinningErri[3 ] =0.0312 ;
  eleBinningErri[4 ] =0.0563 ;
  eleBinningErri[5 ] =0.1331 ;
  eleBinningErri[6 ] =0.1994 ;
  eleBinningErri[7 ] =0.2502 ;
  eleBinningErri[8 ] =0.2922 ;
  eleBinningErri[9 ] =0.3249 ;
  eleBinningErri[10] =0.3471 ;
  eleBinningErri[11] =0.3617 ;
  eleBinningErri[12] =0.3688 ;
  eleBinningErri[13] =0.3527 ;

  double eleEffErri[14]={0};
  for(int i(1);i<14;i++)
  {
    eleEffErri[i]  = sqrt(
	eletoyErri[i] *eletoyErri[i] 
	+eleidisosigErri[i] *eleidisosigErri[i] 
	+eleidisobckErri[i]* eleidisobckErri[i] 
	+eleBinningErri[i] *eleBinningErri[i]);
  } 


  double eleStatErri[14] = {0};
  eleStatErri[1 ]=0.5989   ; 
  eleStatErri[2 ]=0.7374   ;
  eleStatErri[3 ]=0.8890   ;
  eleStatErri[4 ]=0.9533   ;
  eleStatErri[5 ]=1.2811   ;
  eleStatErri[6 ]=1.2831   ;
  eleStatErri[7 ]=1.7124   ;
  eleStatErri[8 ]=1.7542   ;
  eleStatErri[9 ]=2.1637   ;
  eleStatErri[10]=4.4597   ;
  eleStatErri[11]=7.7354   ;
  eleStatErri[12]=11.1409  ;
  eleStatErri[13]=18.0656  ;

  double eleMetErri[14] = {0};
  eleMetErri[1 ] = 0.215627016;
  eleMetErri[2 ] = 0.096979276;
  eleMetErri[3 ] = 0.098676542;
  eleMetErri[4 ] = 0.266666177;
  eleMetErri[5 ] = 0.350696821;
  eleMetErri[6 ] = 0.339040956;
  eleMetErri[7 ] = 0.263083428;
  eleMetErri[8 ] = 0.174917809;
  eleMetErri[9 ] = 0.126110626;
  eleMetErri[10] = 0.139374783;
  eleMetErri[11] = 0.172770513;
  eleMetErri[12] = 0.199154312;
  eleMetErri[13] = 0.203652474;


  double elescaleErri[14]={0};
  elescaleErri[1] = 0.1558;
  elescaleErri[2] = 0.0541;
  elescaleErri[3] = 0.1587;
  elescaleErri[4] = 0.1878;
  elescaleErri[5] = 0.1655;
  elescaleErri[6] = 0.1431;
  elescaleErri[7] = 0.1636;
  elescaleErri[8] = 0.3178;
  elescaleErri[9] = 0.2118;
  elescaleErri[10]= 0.2514;
  elescaleErri[11]= 0.2668;
  elescaleErri[12]= 0.5656;
  elescaleErri[13]= 0.2333;


  double elesmearErri[14]={0};
  elesmearErri[1] =0.1376 ;
  elesmearErri[2] =0.0656 ;
  elesmearErri[3] =0.1787 ;
  elesmearErri[4] =0.2375 ;
  elesmearErri[5] =0.1972 ;
  elesmearErri[6] =0.1775 ;
  elesmearErri[7] =0.2816 ;
  elesmearErri[8] =0.3302 ;
  elesmearErri[9] =0.2826 ;
  elesmearErri[10]=0.4392 ;
  elesmearErri[11]=1.2074 ;
  elesmearErri[12]=0.8745 ;
  elesmearErri[13]=0.5695 ;

  double eleMomResErri[14]={0};
  eleMomResErri[1] =0.2079 ;
  eleMomResErri[2] =0.0850 ;
  eleMomResErri[3] =0.2390 ;
  eleMomResErri[4] =0.3028 ;
  eleMomResErri[5] =0.2574 ;
  eleMomResErri[6] =0.2280 ;
  eleMomResErri[7] =0.3257 ;
  eleMomResErri[8] =0.4583 ;
  eleMomResErri[9] =0.3532 ;
  eleMomResErri[10]=0.5061 ;
  eleMomResErri[11]=1.2365 ;
  eleMomResErri[12]=1.0415 ; 
  eleMomResErri[13]=0.6154 ;


  double eleQCDBckErri[14]={0};
  eleQCDBckErri[1 ] =0.5064  ;
  eleQCDBckErri[2 ] =0.6421  ;
  eleQCDBckErri[3 ] =0.4825  ;
  eleQCDBckErri[4 ] =0.6647  ;
  eleQCDBckErri[5 ] =0.7967  ;
  eleQCDBckErri[6 ] =1.2676  ;
  eleQCDBckErri[7 ] =0.8609  ;
  eleQCDBckErri[8 ] =1.7436  ;
  eleQCDBckErri[9 ] =0.7949  ; 
  eleQCDBckErri[10] =1.3658  ;
  eleQCDBckErri[11] =1.2523  ;
  eleQCDBckErri[12] =1.1945  ;
  eleQCDBckErri[13] =1.7797  ;


  double eleQCDShapeErri[14] = {0};
  eleQCDShapeErri[1 ] =0.1955 ;
  eleQCDShapeErri[2 ] =0.2646 ;
  eleQCDShapeErri[3 ] =0.3677 ;
  eleQCDShapeErri[4 ] =0.4328 ;
  eleQCDShapeErri[5 ] =0.5147 ;
  eleQCDShapeErri[6 ] =0.4029 ;
  eleQCDShapeErri[7 ] =0.4488 ;
  eleQCDShapeErri[8 ] =0.5818 ;
  eleQCDShapeErri[9 ] =0.6259 ; 
  eleQCDShapeErri[10] =0.6202 ;
  eleQCDShapeErri[11] =0.4710 ;
  eleQCDShapeErri[12] =0.6166 ;
  eleQCDShapeErri[13] =0.6632 ;

  double eleEWKErri[14] ={0};
  eleEWKErri[1 ] =0.0490   ;
  eleEWKErri[2 ] =0.0364   ;
  eleEWKErri[3 ] =0.0153   ;
  eleEWKErri[4 ] =0.0407   ;
  eleEWKErri[5 ] =0.0536   ;
  eleEWKErri[6 ] =0.0877   ;
  eleEWKErri[7 ] =0.1233   ;
  eleEWKErri[8 ] =0.1555   ;
  eleEWKErri[9 ] =0.1816   ; 
  eleEWKErri[10] =0.2014   ;
  eleEWKErri[11] =0.2157   ;
  eleEWKErri[12] =0.2251   ;
  eleEWKErri[13] =0.2305   ;

  double eleFSRErri[14] ={0};
  eleFSRErri[1 ] =0.0458 ;
  eleFSRErri[2 ] =0.0512 ;
  eleFSRErri[3 ] =0.0387 ;
  eleFSRErri[4 ] =0.0009 ;
  eleFSRErri[5 ] =0.0592 ;
  eleFSRErri[6 ] =0.1160 ;
  eleFSRErri[7 ] =0.1651 ;
  eleFSRErri[8 ] =0.1968 ;
  eleFSRErri[9 ] =0.2241 ;
  eleFSRErri[10] =0.2530 ;
  eleFSRErri[11] =0.2929 ;
  eleFSRErri[12] =0.2860 ;
  eleFSRErri[13] =0.3427 ;

  double elesvdunfp[14] ={0};
  elesvdunfp[1 ] =0.1153 ;
  elesvdunfp[2 ] =0.1043 ;
  elesvdunfp[3 ] =0.1063 ;
  elesvdunfp[4 ] =0.1227 ;
  elesvdunfp[5 ] =0.1418 ;
  elesvdunfp[6 ] =0.1621 ;
  elesvdunfp[7 ] =0.1875 ;
  elesvdunfp[8 ] =0.2223 ;
  elesvdunfp[9 ] =0.2622 ;
  elesvdunfp[10] =0.3000 ;
  elesvdunfp[11] =0.3312 ;
  elesvdunfp[12] =0.3531 ;
  elesvdunfp[13] =0.3644 ;
  
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

  double eleUnfBiasErri[14] = {0};
  eleUnfBiasErri[1 ] = 0.7531 ; // 0.39 ;
  eleUnfBiasErri[2 ] = 1.4299 ; // 0.39 ;
  eleUnfBiasErri[3 ] = 1.1139 ; // 0.46 ;
  eleUnfBiasErri[4 ] = 0.3641 ; // 0.39 ;
  eleUnfBiasErri[5 ] = 0.5756 ; // 0.61 ;
  eleUnfBiasErri[6 ] = 0.2923 ; // 0.83 ;
  eleUnfBiasErri[7 ] = 0.3374 ; // 1.14 ;
  eleUnfBiasErri[8 ] = 0.4686 ; // 1.55 ;
  eleUnfBiasErri[9 ] = 2.3000 ; // 1.92 ;
  eleUnfBiasErri[10] = 2.3104 ; // 2.15 ;
  eleUnfBiasErri[11] = 4.5665 ; // 2.29 ;
  eleUnfBiasErri[12] = 2.9609 ; // 2.43 ;
  eleUnfBiasErri[13] = 4.0689 ; // 2.42 ;

  // Make Incl
  //*
  TFile *fp = new TFile("../../ResultWpToEleNu/Result_WpToEleNu.root");
  TFile *fm = new TFile("../../ResultWmToEleNu/Result_WmToEleNu.root");

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
  double eleSVDUnfErrp[14]={0};
  double eleSVDUnfErrm[14]={0};
  double eleSVDUnfErri[14]={0};
  cout <<"W+ SVDUnfErr \t W- SVDUnfErr \t W SVDUnfErr" <<endl;
  for(int i(1);i<14;i++)
  {
    eleSVDUnfErrp[i] = h_data_p->GetBinContent(i)*0.01*elesvdunfp[i];
    eleSVDUnfErrm[i] = h_data_m->GetBinContent(i)*0.01*elesvdunfm[i];
    eleSVDUnfErri[i] = (sqrt(eleSVDUnfErrp[i]*eleSVDUnfErrp[i] + eleSVDUnfErrm[i]*eleSVDUnfErrm[i])/(h_data_p->GetBinContent(i)+h_data_m->GetBinContent(i))) * 100; // % unit
    cout << eleSVDUnfErrp[i] << "\t" << eleSVDUnfErrm[i] <<"\t" << eleSVDUnfErri[i] <<endl;
  }

  // Save files
  TString resultDir = "ResultWInclToEleNu";
  gSystem->mkdir(resultDir,kTRUE);

  TFile f_out(resultDir+"/Result_WInclToEleNu_DataPowheg.root","recreate");
 
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
    
    
    TH1D* h_idisosig = new TH1D("h_idisosig","h_idisosig",13,0,13);
    TH1D* h_idisobck = new TH1D("h_idisobck","h_idisobck",13,0,13);
  
    
    TH1D* h_toy = new TH1D("h_toy","h_toy",13,0,13);
    
    TH1D* h_bin = new TH1D("h_bin","h_bin",13,0,13);
    
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
      h_Stat->SetBinContent(i,eleStatErri[i]);
      
      h_idisosig->SetBinContent(i,eleidisosigErri[i]);
      h_idisobck->SetBinContent(i,eleidisobckErri[i]);
      h_toy->SetBinContent(i,eletoyErri[i]);
      h_bin->SetBinContent(i,eleBinningErri[i]);
      
      h_met->SetBinContent(i,eleMetErri[i]);
      
      h_scale->SetBinContent(i,elescaleErri[i]);
      h_smear->SetBinContent(i,elesmearErri[i]);
      
      h_qcdbckgr->SetBinContent(i,eleQCDBckErri[i]);
      h_qcdshape->SetBinContent(i,eleQCDShapeErri[i]);
      
      h_ewk->SetBinContent(i,eleEWKErri[i]);
      h_fsr->SetBinContent(i,eleFSRErri[i]);
      
      h_SvdUnf->SetBinContent(i,eleSVDUnfErri[i]);
      h_UnfoldBias->SetBinContent(i,eleUnfBiasErri[i]);
      
      h_PowhegStat->SetBinContent(i,PowhegStatErrI[i]);
    }

    h_Stat->Write();


    h_idisosig->Write();
    h_idisobck->Write();
    
    h_toy->Write();
    
    h_bin->Write();
    
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
