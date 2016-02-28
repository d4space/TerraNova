#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include "TStyle.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "TMath.h"
#include <iomanip>

using namespace std;

const TString format("pdf");
void FEWZ_PDFUncer(double *WmWpRatio, double *Err);
void Powheg_PDFUncer(double *WmWpRatio, double *Err);
char tmpName[30];

int WMinusWplusRatio_New()
{
  gROOT->LoadMacro("../../Utils/tdrstyle.C");
  setTDRStyle();
  gROOT->LoadMacro("../../Utils/CMS_lumi.C");
  //writeExtraText = "true";
  writeExtraText = false;
  extraText = "Preliminary";
  lumi_8TeV = "18.4 pb^{-1}";

  int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
  int iPos = 0;

  int W = 800;
  int H = 800;

  //
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
  // For instance:
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  //
  int H_ref = 600;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  double BinWidth[14] ={0.0, 7.5-0, 12.5-7.5,17.5-12.5, 24.0-17.5, 30.0-24.0, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0     , 190.0-150.0, 250.0-190.0, 600.0-250.0};
  const int nBins = 13;
  //double WptBins[nBins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
  double WptBins[nBins] = {1.,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

  //Data
  TFile *f_WmToMu = new TFile("../WptCharge_NormDiffXsec_InFid/Wpt_WmToMuNuNormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  TFile *f_WpToMu = new TFile("../WptCharge_NormDiffXsec_InFid/Wpt_WpToMuNuNormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  TH1D* hWmToMu_RD = (TH1D*)f_WmToMu->Get("hData_Xsec_BornLogScaleNorm12")->Clone("hWmToMu_RD");
  TH1D* hWpToMu_RD = (TH1D*)f_WpToMu->Get("hData_Xsec_BornLogScaleNorm12")->Clone("hWpToMu_RD");
  TH1D* hWmWpratio_RD = new TH1D("W/Z ratio", "", nBins-1,WptBins);
 
  ///Powheg
  TFile *f_WmToMu_Powheg = new TFile("../Wpt_PowhegPreFSR_12Bin/root/WmToMuNu_PowhegWpt.root");
  TFile *f_WpToMu_Powheg = new TFile("../Wpt_PowhegPreFSR_12Bin/root/WpToMuNu_PowhegWpt.root");
  TH1D* hWmToMu_Powheg = (TH1D*)f_WmToMu_Powheg->Get("hxsec_NormDiff")->Clone("hWmToMu_Powheg");
  TH1D* hWpToMu_Powheg = (TH1D*)f_WpToMu_Powheg->Get("hxsec_NormDiff")->Clone("hWpToMu_Powheg");
  TH1D* hWmWpratio_Powheg = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  TH1D* hWmWpratioRatioBand_Powheg = new TH1D("W/Z ratio", "", nBins-1,WptBins);

  ///ResBos
  TFile *f_WmToMu_resbos = new TFile("../WpT_Resbos_12Bin/root/WmToMuNu_Resbos.root");
  TFile *f_WpToMu_resbos = new TFile("../WpT_Resbos_12Bin/root/WpToMuNu_Resbos.root");

  TH1D* hWmToMu_Resbos[7];
  TH1D* hWpToMu_Resbos[7];
  TH1D* hWmWpratio_Resbos[7];
  TH1D* hWmWpratio_Resbos_Central = new TH1D("W-/W+ ratio","",nBins-1,WptBins);
  char gridname[30];
  char histName[30];
  for(int i(0);i<7; i++)
  {
    sprintf(gridname,"hResbos%d",i+29);
    sprintf(histName,"hWmToMu_Resbos_%d",i);
    hWmToMu_Resbos[i] = (TH1D*)f_WmToMu_resbos->Get(gridname)->Clone(histName);
    sprintf(histName,"hWpToMu_Resbos_%d",i);
    hWpToMu_Resbos[i] = (TH1D*)f_WpToMu_resbos->Get(gridname)->Clone(histName);
    sprintf(histName,"hWmWpratio_Resbos%d",i);
    hWmWpratio_Resbos[i] = (TH1D*)hWmToMu_Resbos[i]->Clone(histName);
    hWmWpratio_Resbos[i] -> Divide(hWpToMu_Resbos[i]);
  }
  
  ///FEWZ
  TFile *f_WmToMu_fewz = new TFile("../WpT_FEWZ_12Bin/root/WmToMuNu_FEWZ.root");
  TFile *f_WpToMu_fewz = new TFile("../WpT_FEWZ_12Bin/root/WpToMuNu_FEWZ.root");
  TH1D* hWmToMu_FEWZ = (TH1D*)f_WmToMu_fewz->Get("hxsec_NormDiff")->Clone("hWmToMu_FEWZ");
  TH1D* hWpToMu_FEWZ = (TH1D*)f_WpToMu_fewz->Get("hxsec_NormDiff")->Clone("hWpToMu_FEWZ");
  TH1D* hWmToMu_Up_FEWZ = (TH1D*)f_WmToMu_fewz->Get("hxsec_NormDiff_up")->Clone("hWmToMu_Up_FEWZ");
  TH1D* hWmToMu_Down_FEWZ = (TH1D*)f_WmToMu_fewz->Get("hxsec_NormDiff_down")->Clone("hWmToMu_Down_FEWZ");
  TH1D* hWpToMu_Up_FEWZ = (TH1D*)f_WpToMu_fewz->Get("hxsec_NormDiff_up")->Clone("hWpToMu_Up_FEWZ");
  TH1D* hWpToMu_Down_FEWZ = (TH1D*)f_WpToMu_fewz->Get("hxsec_NormDiff_down")->Clone("hWpToMu_Down_FEWZ");
  TH1D* hWmWpratio_FEWZ = new TH1D("W/Z ratio", "", nBins-1,WptBins);


  //Systematics
  TH1D* mTrackSigShapeSyst12 = (TH1D*)f_WmToMu->Get("hData_TrackSigShapeSystErr12")->Clone("mTrackSigShapeSyst12");
  TH1D* mTrackBkgrShapeSyst12 = (TH1D*)f_WmToMu->Get("hData_TrackBkgrShapeSystErr12")->Clone("mTrackBkgrShapeSyst12");
  TH1D* mIDIsoSigShapeSyst12 = (TH1D*)f_WmToMu->Get("hData_IDIsoSigShapeSystErr12")->Clone("mIDIsoSigShapeSyst12");
  TH1D* mIDIsoBkgrShapeSyst12 = (TH1D*)f_WmToMu->Get("hData_IDIsoBkgrShapeSystErr12")->Clone("mIDIsoBkgrShapeSyst12");
  TH1D* mMuonPOGSyst12 = (TH1D*)f_WmToMu->Get("hData_MuonPOGSystErr12")->Clone("mMuonPOGSyst12");
  TH1D* mEffiToySyst12 = (TH1D*)f_WmToMu->Get("hData_EffiToySystErr12")->Clone("mEffiToySyst12");
  TH1D* mStat12 = (TH1D*)f_WmToMu->Get("hData_StatErr12")->Clone("mStat12");
  TH1D* mMetResolSyst12 = (TH1D*)f_WmToMu->Get("hData_MetResolSystErr12")->Clone("mMetResolSyst12");
  TH1D* mEnMomScaleSyst12 = (TH1D*)f_WmToMu->Get("hData_EnMomScaleSystErr12")->Clone("mEnMomScaleSyst12");
  TH1D* mEnMomSmearSyst12 = (TH1D*)f_WmToMu->Get("hData_EnMomSmearSystErr12")->Clone("mEnMomSmearSyst12");
  TH1D* mQcdBckgrSyst12 = (TH1D*)f_WmToMu->Get("hData_QcdBckgrSystErr12")->Clone("mQcdBckgrSyst12");
  TH1D* mQcdShapeSyst12 = (TH1D*)f_WmToMu->Get("hData_QcdShapeSystErr12")->Clone("mQcdShapeSyst12");
  TH1D* mEwkSyst12 = (TH1D*)f_WmToMu->Get("hData_EwkSystErr12")->Clone("mEwkSyst12");
  TH1D* mFsrSyst12 = (TH1D*)f_WmToMu->Get("hData_FsrSystErr12")->Clone("mFsrSyst12");
  TH1D* mSvdUnfSyst12 = (TH1D*)f_WmToMu->Get("hData_SvdUnfSystErr12")->Clone("mSvdUnfSyst12");
  TH1D* mUnfoldBiasSyst12 = (TH1D*)f_WmToMu->Get("hData_UnfoldBiasSystErr12")->Clone("mUnfoldBiasSyst12");
  TH1D* mTotalSyst12 = new TH1D("TotalSyst","TotalSyst",nBins-1,0,nBins-1);
  TH1D* mTotalUncer12 = new TH1D("TotalUncer","TotalUncer",nBins-1,0,nBins-1);

  TH1D* pTrackSigShapeSyst12 = (TH1D*)f_WpToMu->Get("hData_TrackSigShapeSystErr12")->Clone("pTrackSigShapeSyst12");
  TH1D* pTrackBkgrShapeSyst12 = (TH1D*)f_WpToMu->Get("hData_TrackBkgrShapeSystErr12")->Clone("pTrackBkgrShapeSyst12");
  TH1D* pIDIsoSigShapeSyst12 = (TH1D*)f_WpToMu->Get("hData_IDIsoSigShapeSystErr12")->Clone("pIDIsoSigShapeSyst12");
  TH1D* pIDIsoBkgrShapeSyst12 = (TH1D*)f_WpToMu->Get("hData_IDIsoBkgrShapeSystErr12")->Clone("pIDIsoBkgrShapeSyst12");
  TH1D* pMuonPOGSyst12 = (TH1D*)f_WpToMu->Get("hData_MuonPOGSystErr12")->Clone("pMuonPOGSyst12");
  TH1D* pEffiToySyst12 = (TH1D*)f_WpToMu->Get("hData_EffiToySystErr12")->Clone("pEffiToySyst12");
  TH1D* pStat12 = (TH1D*)f_WpToMu->Get("hData_StatErr12")->Clone("pStat12");
  TH1D* pMetResolSyst12 = (TH1D*)f_WpToMu->Get("hData_MetResolSystErr12")->Clone("pMetResolSyst12");
  TH1D* pEnMomScaleSyst12 = (TH1D*)f_WpToMu->Get("hData_EnMomScaleSystErr12")->Clone("pEnMomScaleSyst12");
  TH1D* pEnMomSmearSyst12 = (TH1D*)f_WpToMu->Get("hData_EnMomSmearSystErr12")->Clone("pEnMomSmearSyst12");
  TH1D* pQcdBckgrSyst12 = (TH1D*)f_WpToMu->Get("hData_QcdBckgrSystErr12")->Clone("pQcdBckgrSyst12");
  TH1D* pQcdShapeSyst12 = (TH1D*)f_WpToMu->Get("hData_QcdShapeSystErr12")->Clone("pQcdShapeSyst12");
  TH1D* pEwkSyst12 = (TH1D*)f_WpToMu->Get("hData_EwkSystErr12")->Clone("pEwkSyst12");
  TH1D* pFsrSyst12 = (TH1D*)f_WpToMu->Get("hData_FsrSystErr12")->Clone("pFsrSyst12");
  TH1D* pSvdUnfSyst12 = (TH1D*)f_WpToMu->Get("hData_SvdUnfSystErr12")->Clone("pSvdUnfSyst12");
  TH1D* pUnfoldBiasSyst12 = (TH1D*)f_WpToMu->Get("hData_UnfoldBiasSystErr12")->Clone("pUnfoldBiasSyst12");
  TH1D* pTotalSyst12 = new TH1D("TotalSyst","TotalSyst",nBins-1,0,nBins-1);
  TH1D* pTotalUncer12 = new TH1D("TotalUncer","TotalUncer",nBins-1,0,nBins-1);

  cout<<"Bin#: "<<pTrackSigShapeSyst12->GetNbinsX()<<endl;
  double TrackSigSyst[14] = {0.,};
  TrackSigSyst[1] =  0.1021;
  TrackSigSyst[2] =  0.0136;
  TrackSigSyst[3] =  0.0919;
  TrackSigSyst[4] =  0.1451;
  TrackSigSyst[5] =  0.1298;
  TrackSigSyst[6] =  0.0830;
  TrackSigSyst[7] =  0.0413;
  TrackSigSyst[8] =  0.0163;
  TrackSigSyst[9] =  0.0061;
  TrackSigSyst[10] = 0.0033;
  TrackSigSyst[11] = 0.0037;
  TrackSigSyst[12] = 0.0049;
  TrackSigSyst[13] = 0.0057;

  double TrackSigSystMerge[14] = {0};
  //TrackSigSystMerge[4]=(BinWidth[4]*TrackSigSyst[4] + BinWidth[5]*TrackSigSyst[5]) / (BinWidth[4] + BinWidth[5]);
  TrackSigSystMerge[4]=(
      BinWidth[4]*TrackSigSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4))+
      BinWidth[5]*TrackSigSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double TrackSigSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) TrackSigSyst12[i] = TrackSigSyst[i];
    if(i==4) TrackSigSyst12[4] = TrackSigSystMerge[4];
    if(i>=5) TrackSigSyst12[i] = TrackSigSyst[i+1];
  cout << "TrackSigSyst : " << TrackSigSyst12[i] << endl;
  }

  double TrackBkgSyst[14] = {0.,};
  TrackBkgSyst[1] =  0.0580;
  TrackBkgSyst[2] =  0.0080;
  TrackBkgSyst[3] =  0.0759;
  TrackBkgSyst[4] =  0.0928;
  TrackBkgSyst[5] =  0.0581;
  TrackBkgSyst[6] =  0.0111;
  TrackBkgSyst[7] =  0.0181;
  TrackBkgSyst[8] =  0.0262;
  TrackBkgSyst[9] =  0.0205;
  TrackBkgSyst[10] = 0.0103;
  TrackBkgSyst[11] = 0.0002;
  TrackBkgSyst[12] = 0.0074;
  TrackBkgSyst[13] = 0.0115;
  
  double TrackBkgSystMerge[14] = {0};
  TrackBkgSystMerge[4]=(
      BinWidth[4]*TrackBkgSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*TrackBkgSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double TrackBkgSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) TrackBkgSyst12[i] = TrackBkgSyst[i];
    if(i==4) TrackBkgSyst12[4] = TrackBkgSystMerge[4];
    if(i>=5) TrackBkgSyst12[i] = TrackBkgSyst[i+1];
  cout << "TrackBkgSyst : " << TrackBkgSyst12[i] << endl;
  }

  double IDIsoSigSyst[14] = {0.,};
  IDIsoSigSyst[1] =  0.0676;
  IDIsoSigSyst[2] =  0.0137;
  IDIsoSigSyst[3] =  0.1018;
  IDIsoSigSyst[4] =  0.1294;
  IDIsoSigSyst[5] =  0.0874;
  IDIsoSigSyst[6] =  0.0171;
  IDIsoSigSyst[7] =  0.0455;
  IDIsoSigSyst[8] =  0.0909;
  IDIsoSigSyst[9] =  0.1220;
  IDIsoSigSyst[10] = 0.1445;
  IDIsoSigSyst[11] = 0.1606;
  IDIsoSigSyst[12] = 0.1712;
  IDIsoSigSyst[13] = 0.1764;
  
  double IDIsoSigSystMerge[14] = {0};
  IDIsoSigSystMerge[4]=(
      BinWidth[4]*IDIsoSigSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*IDIsoSigSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double IDIsoSigSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) IDIsoSigSyst12[i] = IDIsoSigSyst[i];
    if(i==4) IDIsoSigSyst12[4] = IDIsoSigSystMerge[4];
    if(i>=5) IDIsoSigSyst12[i] = IDIsoSigSyst[i+1];
  cout << "IDIsoSigSyst : " << IDIsoSigSyst12[i] << endl;
  }

  double IDIsoBkgSyst[14] = {0.,};
  IDIsoBkgSyst[1] =  0.0679;
  IDIsoBkgSyst[2] =  0.0040;
  IDIsoBkgSyst[3] =  0.0813;
  IDIsoBkgSyst[4] =  0.1067;
  IDIsoBkgSyst[5] =  0.0760;
  IDIsoBkgSyst[6] =  0.0277;
  IDIsoBkgSyst[7] =  0.0078;
  IDIsoBkgSyst[8] =  0.0259;
  IDIsoBkgSyst[9] =  0.0326;
  IDIsoBkgSyst[10] = 0.0350;
  IDIsoBkgSyst[11] = 0.0359;
  IDIsoBkgSyst[12] = 0.0363;
  IDIsoBkgSyst[13] = 0.0365;
  
  double IDIsoBkgSystMerge[14] = {0};
  //IDIsoBkgSystMerge[4]=(BinWidth[4]*IDIsoBkgSyst[4] + BinWidth[5]*IDIsoBkgSyst[5]) / (BinWidth[4] + BinWidth[5]);
  IDIsoBkgSystMerge[4]=(
      BinWidth[4]*IDIsoBkgSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*IDIsoBkgSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double IDIsoBkgSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) IDIsoBkgSyst12[i] = IDIsoBkgSyst[i];
    if(i==4) IDIsoBkgSyst12[4] = IDIsoBkgSystMerge[4];
    if(i>=5) IDIsoBkgSyst12[i] = IDIsoBkgSyst[i+1];
  cout << "IDIsoBkgSyst : " << IDIsoBkgSyst12[i] << endl;
  }

  double POGSyst[14] = {0.,};
  POGSyst[1] =  0.0754;
  POGSyst[2] =  0.1666;
  POGSyst[3] =  0.0602;
  POGSyst[4] =  0.2756;
  POGSyst[5] =  0.5310;
  POGSyst[6] =  0.4934;
  POGSyst[7] =  0.1987;
  POGSyst[8] =  0.1714;
  POGSyst[9] =  0.4491;
  POGSyst[10] = 0.6324;
  POGSyst[11] = 0.7448;
  POGSyst[12] = 0.8092;
  POGSyst[13] = 0.8389;
  
  double POGSystMerge[14] = {0};
  POGSystMerge[4]=(
      BinWidth[4]*POGSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*POGSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double POGSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) POGSyst12[i] = POGSyst[i];
    if(i==4) POGSyst12[4] = POGSystMerge[4];
    if(i>=5) POGSyst12[i] = POGSyst[i+1];
  cout << "POGSyst : " << POGSyst12[i] << endl;
  }
		
  double LepRecToySyst[14] = {0.,};
  LepRecToySyst[1]  =0.2024  ;
  LepRecToySyst[2]  =0.0975  ;
  LepRecToySyst[3]  =0.2092  ;
  LepRecToySyst[4]  =0.3244  ;
  LepRecToySyst[5]  =0.3395  ;
  LepRecToySyst[6]  =0.3180  ;
  LepRecToySyst[7]  =0.3084  ;
  LepRecToySyst[8]  =0.3233  ;
  LepRecToySyst[9]  =0.3573  ;
  LepRecToySyst[10] =0.4019  ;
  LepRecToySyst[11] =0.4451  ;
  LepRecToySyst[12] =0.4782  ;
  LepRecToySyst[13] =0.4958  ;
  
  double LepRecToySystMerge[14] = {0};
  LepRecToySystMerge[4]=(
      BinWidth[4]*LepRecToySyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*LepRecToySyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double LepRecToySyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) LepRecToySyst12[i] = LepRecToySyst[i];
    if(i==4) LepRecToySyst12[4] = LepRecToySystMerge[4];
    if(i>=5) LepRecToySyst12[i] = LepRecToySyst[i+1];
  cout << "LepRecToySyst : " << LepRecToySyst12[i] << endl;
  }


  double RecoilSyst[14] = {0.,};
  RecoilSyst[1]  =  0.0630;
  RecoilSyst[2]  =  0.0360;
  RecoilSyst[3]  =  0.0759;
  RecoilSyst[4]  =  0.1161;
  RecoilSyst[5]  =  0.1235;
  RecoilSyst[6]  =  0.1200;
  RecoilSyst[7]  =  0.1211;
  RecoilSyst[8]  =  0.1425;
  RecoilSyst[9]  =  0.1936;
  RecoilSyst[10] =  0.2608;
  RecoilSyst[11] =  0.3244;
  RecoilSyst[12] =  0.3720;
  RecoilSyst[13] =  0.3971;
  
  double RecoilSystMerge[14] = {0};
  RecoilSystMerge[4]=(
      BinWidth[4]*RecoilSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*RecoilSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double RecoilSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) RecoilSyst12[i] = RecoilSyst[i];
    if(i==4) RecoilSyst12[4] = RecoilSystMerge[4];
    if(i>=5) RecoilSyst12[i] = RecoilSyst[i+1];
  cout << "RecoilSyst : " << RecoilSyst12[i] << endl;
  }



  double ScaleSyst[14] = {0.,};
  ScaleSyst[1] =  0.4550;
  ScaleSyst[2] =  0.2724;
  ScaleSyst[3] =  0.0802;
  ScaleSyst[4] =  0.4758;
  ScaleSyst[5] =  0.7510;
  ScaleSyst[6] =  0.8612;
  ScaleSyst[7] =  0.8540;
  ScaleSyst[8] =  0.7908;
  ScaleSyst[9] =  0.7099;
  ScaleSyst[10] = 0.6312;
  ScaleSyst[11] = 0.5660;
  ScaleSyst[12] = 0.5202;
  ScaleSyst[13] = 0.4967;
  
  double ScaleSystMerge[14] = {0};
  //ScaleSystMerge[4]=(BinWidth[4]*ScaleSyst[4] + BinWidth[5]*ScaleSyst[5]) / (BinWidth[4] + BinWidth[5]);
  ScaleSystMerge[4]=(
      BinWidth[4]*ScaleSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + 
      BinWidth[5]*ScaleSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5))) 
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double ScaleSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) ScaleSyst12[i] = ScaleSyst[i];
    if(i==4) ScaleSyst12[4] = ScaleSystMerge[4];
    if(i>=5) ScaleSyst12[i] = ScaleSyst[i+1];
  cout << "ScaleSyst : " << ScaleSyst12[i] << endl;
  }

  double SmearSyst[14] = {0.,};
  SmearSyst[1] =  0.4146;
  SmearSyst[2] =  0.2165;
  SmearSyst[3] =  0.1626;
  SmearSyst[4] =  0.5080;
  SmearSyst[5] =  0.7051;
  SmearSyst[6] =  0.7199;
  SmearSyst[7] =  0.6292;
  SmearSyst[8] =  0.5122;
  SmearSyst[9] =  0.4098;
  SmearSyst[10] = 0.3293;
  SmearSyst[11] = 0.2695;
  SmearSyst[12] = 0.2297;
  SmearSyst[13] = 0.2099;
  
  double SmearSystMerge[14] = {0};
  SmearSystMerge[4]=(
      BinWidth[4]*SmearSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*SmearSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double SmearSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) SmearSyst12[i] = SmearSyst[i];
    if(i==4) SmearSyst12[4] = SmearSystMerge[4];
    if(i>=5) SmearSyst12[i] = SmearSyst[i+1];
  cout << "SmearSyst : " << SmearSyst12[i] << endl;
  }
 
  double FSRSyst[14] = {0.,};
  FSRSyst[1]  = 0.0037   ;
  FSRSyst[2]  = 0.0020   ;
  FSRSyst[3]  = 0.0011   ;
  FSRSyst[4]  = 0.0048   ;
  FSRSyst[5]  = 0.0074   ;
  FSRSyst[6]  = 0.0078   ;
  FSRSyst[7]  = 0.0062   ;
  FSRSyst[8]  = 0.0036   ;
  FSRSyst[9]  = 0.0010   ;
  FSRSyst[10] = 0.0012   ;
  FSRSyst[11] = 0.0029   ;
  FSRSyst[12] = 0.0040   ;
  FSRSyst[13] = 0.0046   ;
  
  double FSRSystMerge[14] = {0};
  FSRSystMerge[4]=(
      BinWidth[4]*FSRSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*FSRSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double FSRSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) FSRSyst12[i] = FSRSyst[i];
    if(i==4) FSRSyst12[4] = FSRSystMerge[4];
    if(i>=5) FSRSyst12[i] = FSRSyst[i+1];
  cout << "FSRSyst:   " <<i<<"\t"<< FSRSyst12[i] << endl;
  }
  

  double EWKSyst[14] = {0.,};
  EWKSyst[1]  =0.012739805   ;
  EWKSyst[2]  =0.008979764   ;
  EWKSyst[3]  =0.018967321   ;
  EWKSyst[4]  =0.031908746   ;
  EWKSyst[5]  =0.030499055   ;
  EWKSyst[6]  =0.016025577   ;
  EWKSyst[7]  =0.034669487   ;
  EWKSyst[8]  =0.055872249   ;
  EWKSyst[9]  =0.074163038   ;
  EWKSyst[10] =0.088341626   ;
  EWKSyst[11] =0.098635467   ;
  EWKSyst[12] =0.105338857   ;
  EWKSyst[13] =0.108641899   ;
  
  double EWKSystMerge[14] = {0};
  EWKSystMerge[4]=(
      BinWidth[4]*EWKSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) +
      BinWidth[5]*EWKSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)))
    / (BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5)));

  double EWKSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) EWKSyst12[i] = EWKSyst[i];
    if(i==4) EWKSyst12[4] = EWKSystMerge[4];
    if(i>=5) EWKSyst12[i] = EWKSyst[i+1];
  cout << "EWKSyst:   " <<i<<"\t"<< EWKSyst12[i] << endl;
  }
 

  double UnfBiasSyst[14] = {0.,};
  UnfBiasSyst[1]  = 0.8286 ;  // 0.69 ;
  UnfBiasSyst[2]  = 0.1332 ;  // 0.56 ;
  UnfBiasSyst[3]  = 1.2219 ;  // 0.62 ;
  UnfBiasSyst[4]  = 1.1937 ;  // 0.93 ;
  UnfBiasSyst[5]  = 2.0407 ;  // 1.24 ;
  UnfBiasSyst[6]  = 0.7443 ;  // 1.45 ;
  UnfBiasSyst[7]  = 0.3415 ;  // 1.81 ;
  UnfBiasSyst[8]  = 1.2224 ;  // 2.49 ;
  UnfBiasSyst[9]  = 1.9469 ;  // 2.77 ;
  UnfBiasSyst[10] = 1.7871 ;  // 2.49 ;
  UnfBiasSyst[11] = 8.7839 ;  // 2.47 ;
  UnfBiasSyst[12] = 3.3604 ;  // 2.86 ;
  UnfBiasSyst[13] = 4.4327 ;  // 3.07 ;
  
  double UnfBiasSystMerge[14] = {0};
  UnfBiasSystMerge[4]=(
      			BinWidth[4]*UnfBiasSyst[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) 
			+ 
			BinWidth[5]*UnfBiasSyst[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5))
			) 
    			/ 
			(
			 BinWidth[4]*(hWmToMu_RD->GetBinContent(4)/hWpToMu_RD->GetBinContent(4)) 
			 + BinWidth[5]*(hWmToMu_RD->GetBinContent(5)/hWpToMu_RD->GetBinContent(5))
			 );

  double UnfBiasSyst12[13] = {0};
  for (int i(1); i<13 ; i++)
  {
    if(i<4) UnfBiasSyst12[i] = UnfBiasSyst[i];
    if(i==4) UnfBiasSyst12[4] = UnfBiasSystMerge[4];
    if(i>=5) UnfBiasSyst12[i] = UnfBiasSyst[i+1];
  cout << "UnfBiasSyst12 : " << UnfBiasSyst12[i] << endl;
  }
  
  for( int ipt(1);ipt<=12;ipt++)
    {
      ///used for Normalized Xsec 12
      mTotalSyst12->SetBinContent(ipt, sqrt( 
           // + TMath::Power(mEffiToySyst12->GetBinContent(ipt),2)
           //+  TMath::Power(mMetResolSyst12->GetBinContent(ipt),2)
            TMath::Power(mQcdBckgrSyst12->GetBinContent(ipt),2)
           +  TMath::Power(mQcdShapeSyst12->GetBinContent(ipt),2)
           //// + TMath::Power(mEwkSyst12->GetBinContent(ipt),2) 
           //// + TMath::Power(mFsrSyst12->GetBinContent(ipt),2) 
           + TMath::Power(mSvdUnfSyst12->GetBinContent(ipt),2)
            ));

      //cout<<"mUnfoldBiasSyst12\t"<<UnfBiasSyst12[ipt]<<endl;
      cout<<"mQcdBckgrSyst12\t"<<mQcdBckgrSyst12->GetBinContent(ipt)<<endl;
  }
    for( int ipt(1);ipt<=12;ipt++)
    {
      ///used for Normalized Xsec 12
      mTotalUncer12->SetBinContent(ipt, sqrt( 
            TMath::Power(mTotalSyst12->GetBinContent(ipt),2)
	   + TMath::Power(mStat12->GetBinContent(ipt),2)
	    ));
	    }
    for( int ipt(1);ipt<=12;ipt++)
    {
      ///used for Normalized Xsec 12
      pTotalSyst12->SetBinContent(ipt, sqrt( 
            //// + TMath::Power(pEffiToySyst12->GetBinContent(ipt),2)
            ////+ TMath::Power(pMetResolSyst12->GetBinContent(ipt),2)
            TMath::Power(pQcdBckgrSyst12->GetBinContent(ipt),2)
           + TMath::Power(pQcdShapeSyst12->GetBinContent(ipt),2)
           //// + TMath::Power(pEwkSyst12->GetBinContent(ipt),2) 
           //// + TMath::Power(pFsrSyst12->GetBinContent(ipt),2) 
            + TMath::Power(pSvdUnfSyst12->GetBinContent(ipt),2)
	    ));
      cout<<"pQcdBckgrSyst12\t"<<pQcdBckgrSyst12->GetBinContent(ipt)<<endl;

  }
    for( int ipt(1);ipt<=12;ipt++)
    {
      ///used for Normalized Xsec 12
      pTotalUncer12->SetBinContent(ipt, sqrt( 
            TMath::Power(pTotalSyst12->GetBinContent(ipt),2)
	  + TMath::Power(pStat12->GetBinContent(ipt),2)
	    ));
	    }
    for(int ipt(1); ipt<=12;ipt++)
    {
      cout<<"mTotalSystErr12: "<<ipt<<"\t"<<mTotalSyst12->GetBinContent(ipt)<<endl;
    }
    for(int ipt(1); ipt<=12;ipt++)
    {
      cout<<"pTotalSystErr12: "<<ipt<<"\t"<<pTotalSyst12->GetBinContent(ipt)<<endl;
    }


  double WmToMu_RD[12]={0};
  double WpToMu_RD[12]={0};
  double WmToMuErr_RD[12]={0};
  double WpToMuErr_RD[12]={0};
 
  double WmToMu_Powheg[12]={0};
  double WpToMu_Powheg[12]={0};
  double WmToMuStatErr_Powheg[12]={0};
  double WpToMuStatErr_Powheg[12]={0};
 
  double WmToMu_Resbos[12]={0};
  double WpToMu_Resbos[12]={0};
  double WmToMuErr_Resbos[12]={0};
  double WpToMuErr_Resbos[12]={0};
 
  double WmToMu_FEWZ[12]={0};
  double WpToMu_FEWZ[12]={0};
  double WmToMuStatErr_FEWZ[12]={0};
  double WpToMuStatErr_FEWZ[12]={0};
 
  double WmWpratio_RD[12]={0};
  double WmWpratioErr_RD_NoCorr[12]={0};
  double WmWpratioErr_RD_FullCorr[12]={0};
  double WmWpratioErr_RD_Total[12]={0};

  double WmWpratio_Powheg[12]={0};
  double WmWpratioStatErr_Powheg[12]={0}; // Stat error by error propagation
  double WmWpratioPDFErr_Powheg[12]={0}; // W-/W+ PDF error
  double WmWpratioTotalErr_Powheg[12]={0}; // sqrt(StatErr^2 + PDFErr^2 + ScaleErr^2)
  
  double WmWpratio_FEWZ[12]={0};
  double WmWp_Up_ratio_FEWZ[12]={0};
  double WmWp_Down_ratio_FEWZ[12]={0};
  double WmWpratioStatErr_FEWZ[12]={0}; // Stat error by error propagation
  double WmWpratioPDFErr_FEWZ[12]={0}; // W-/W+ PDF error
  double WmWpratioScaleErr_FEWZ[12]={0}; // W-/W+ scale error
  double WmWpratioTotalErr_FEWZ[12]={0}; // sqrt(StatErr^2 + PDFErr^2 + ScaleErr^2)
 
  Double_t Resb_errMax[nBins-1];
  Double_t Resb_errMin[nBins-1];
  Double_t Resb_err[nBins-1];
 
  cout << fixed << setprecision(10) << endl;
  cout << " ===== Wpt minus and Wptplus Normalilzed Differntial cross-section and errors In Fiducial Volume ==="<< endl;
  cout << "Bin"<<
    "\t Wpt minus  " << " \t\tError " << 
    "\t\t\t Wpt plus  " << " \t\tError " << endl;
  for(int i=0; i<12; i++)
  {
    WmToMu_RD[i] = hWmToMu_RD->GetBinContent(i+1); 
    WpToMu_RD[i] = hWpToMu_RD->GetBinContent(i+1); 
    WmToMuErr_RD[i] = 0.01*hWmToMu_RD->GetBinContent(i+1)*mTotalUncer12->GetBinContent(i+1); 
    WpToMuErr_RD[i] = 0.01*hWpToMu_RD->GetBinContent(i+1)*pTotalUncer12->GetBinContent(i+1); 
    
    //cout <<i<<"\t"<< hWmToMu_RD->GetBinContent(i+1) << "\t"  << WmToMuErr_RD[i] << "\t\t" << hWpToMu_RD->GetBinContent(i+1) << "\t" << WpToMuErr_RD[i] <<endl;
   
    WmWpratio_RD[i] = WmToMu_RD[i] / WpToMu_RD[i] ; 
   // cout <<i<< 
   //   "\t" <<WmToMu_RD[i] << "\t" <<  WmToMuErr_RD[i] << 
   //   "\t\t" << WpToMu_RD[i] << "\t" << WpToMuErr_RD[i] << 
   //   "\t\t" << WmWpratio_RD[i] <<
   //   "\t\t" << WmToMuErr_RD[i]*WmToMuErr_RD[i]/WmToMu_RD[i]/WmToMu_RD[i] <<
   //   "\t\t" << WpToMuErr_RD[i]*WpToMuErr_RD[i]/WpToMu_RD[i]/WpToMu_RD[i] <<
   //   
   //   endl;
    
    WmWpratioErr_RD_NoCorr[i] = WmWpratio_RD[i]*sqrt(
	    WmToMuErr_RD[i]*WmToMuErr_RD[i]/WmToMu_RD[i]/WmToMu_RD[i] 
	  + WpToMuErr_RD[i]*WpToMuErr_RD[i]/WpToMu_RD[i]/WpToMu_RD[i]
	  ); 
    
    //cout << "WmWpratioErr_RD_NoCorr : " << WmWpratioErr_RD_NoCorr[i] << "\t in %\t"<< WmWpratioErr_RD_NoCorr[i]/WmWpratio_RD[i]*100 << endl; 

   // cout<< " & " << WmWpratioErr_RD_NoCorr[i]/WmWpratio_RD[i]*100;
    //cout<< " & " << WmWpratioErr_RD_NoCorr[i];
    
    WmWpratioErr_RD_FullCorr[i] = sqrt(
	  TrackSigSyst12[i+1]*TrackSigSyst12[i+1] 
	+ TrackBkgSyst12[i+1]*TrackBkgSyst12[i+1] 
	+ IDIsoSigSyst12[i+1]*IDIsoSigSyst12[i+1] 
	+ IDIsoBkgSyst12[i+1]*IDIsoBkgSyst12[i+1] 
	+ POGSyst12[i+1]*POGSyst12[i+1]
	+ LepRecToySyst12[i+1]*LepRecToySyst12[i+1] 
	 
	+ RecoilSyst12[i+1]*RecoilSyst12[i+1] 

	+ ScaleSyst12[i+1]*ScaleSyst12[i+1] 
	+ SmearSyst12[i+1]*SmearSyst12[i+1] 

        + FSRSyst12[i+1]*FSRSyst12[i+1] 
         
	+ EWKSyst12[i+1]*EWKSyst12[i+1] 
	 
	+ UnfBiasSyst12[i+1]*UnfBiasSyst12[i+1]	
	) * WmWpratio_RD[i]*0.01;
	

    //cout << "WmWpratioErr_RD_FullCorr LepRec syst:\t " <<i<<"\t"<< WmWpratioErr_RD_FullCorr[i] << endl; 
    //cout<< " & " << WmWpratioErr_RD_FullCorr[i];
    
    WmWpratioErr_RD_Total[i] = sqrt(WmWpratioErr_RD_NoCorr[i]*WmWpratioErr_RD_NoCorr[i] + WmWpratioErr_RD_FullCorr[i]*WmWpratioErr_RD_FullCorr[i]); 
    //cout<< " & " << WmWpratioErr_RD_Total[i];

    ///Powheg & Resbos ratio
    WmWpratio_Powheg[i] = hWmToMu_Powheg->GetBinContent(i+1) /hWpToMu_Powheg->GetBinContent(i+1) ; 
    //WmWpratio_Resbos[i] = hWmToMu_Resbos->GetBinContent(i+1) /hWpToMu_Resbos->GetBinContent(i+1) ; 
    WmWpratio_FEWZ[i] = hWmToMu_FEWZ->GetBinContent(i+1) /hWpToMu_FEWZ->GetBinContent(i+1) ; 
    WmWp_Up_ratio_FEWZ[i] = hWmToMu_Up_FEWZ->GetBinContent(i+1) /hWpToMu_Up_FEWZ->GetBinContent(i+1) ; 
    WmWp_Down_ratio_FEWZ[i] = hWmToMu_Down_FEWZ->GetBinContent(i+1) /hWpToMu_Down_FEWZ->GetBinContent(i+1) ; 
    WmWpratioScaleErr_FEWZ[i] = TMath::Max(fabs(WmWp_Up_ratio_FEWZ[i]-WmWpratio_FEWZ[i]),fabs(WmWp_Down_ratio_FEWZ[i]-WmWpratio_FEWZ[i]));

    // FEWZ ratio error propagation
    WmToMu_FEWZ[i] = hWmToMu_FEWZ->GetBinContent(i+1);
    WmToMuStatErr_FEWZ[i] = hWmToMu_FEWZ->GetBinError(i+1);
    
    WpToMu_FEWZ[i] = hWpToMu_FEWZ->GetBinContent(i+1);
    WpToMuStatErr_FEWZ[i] = hWpToMu_FEWZ->GetBinError(i+1);

    WmWpratioStatErr_FEWZ[i] = WmWpratio_FEWZ[i]* TMath::Sqrt(
	  (WmToMuStatErr_FEWZ[i]*WmToMuStatErr_FEWZ[i]/WmToMu_FEWZ[i]/WmToMu_FEWZ[i]) 
	+ (WpToMuStatErr_FEWZ[i]*WpToMuStatErr_FEWZ[i]/WpToMu_FEWZ[i]/WpToMu_FEWZ[i])
	); // FEWZ Stat Error propagation

    FEWZ_PDFUncer(WmWpratio_FEWZ,WmWpratioPDFErr_FEWZ);
   // cout << "WmWpratioPDFErr_FEWZ in main : " << WmWpratioPDFErr_FEWZ[i] << endl; // check PDF error number is correctly retured
 
    WmWpratioTotalErr_FEWZ[i] = sqrt(
	  WmWpratioStatErr_FEWZ[i]*WmWpratioStatErr_FEWZ[i] 
	+ WmWpratioPDFErr_FEWZ[i]*WmWpratioPDFErr_FEWZ[i] 
	+ WmWpratioScaleErr_FEWZ[i]*WmWpratioScaleErr_FEWZ[i]
	);

    // Powheg ratio error propagation
    WmToMu_Powheg[i] = hWmToMu_Powheg->GetBinContent(i+1);
    WmToMuStatErr_Powheg[i] = hWmToMu_Powheg->GetBinError(i+1);
    
    WpToMu_Powheg[i] = hWpToMu_Powheg->GetBinContent(i+1);
    WpToMuStatErr_Powheg[i] = hWpToMu_Powheg->GetBinError(i+1);

    WmWpratioStatErr_Powheg[i] =  WmWpratio_Powheg[i]* TMath::Sqrt(
	  (WmToMuStatErr_Powheg[i]*WmToMuStatErr_Powheg[i]/WmToMu_Powheg[i]/WmToMu_Powheg[i]) 
	+ (WpToMuStatErr_Powheg[i]*WpToMuStatErr_Powheg[i]/WpToMu_Powheg[i]/WpToMu_Powheg[i])
	);
 
    
    Powheg_PDFUncer(WmWpratio_Powheg,WmWpratioPDFErr_Powheg);
   // cout << "WmWpratioPDFErr_Powheg in main : " << WmWpratioPDFErr_Powheg[i] << "\t in %\t"<< 100*WmWpratioPDFErr_Powheg[i]/WmWpratio_Powheg[i]<<endl; // check PDF error number is correctly retured
    
    ///Total powheg error calculated here
    WmWpratioTotalErr_Powheg[i] = sqrt(WmWpratioStatErr_Powheg[i]*WmWpratioStatErr_Powheg[i] +WmWpratioPDFErr_Powheg[i]*WmWpratioPDFErr_Powheg[i]);

     // Resbos ratio error in normdiff stage
    double tmpVal,tmpDiff;
    double nomVal = hWmWpratio_Resbos[1]->GetBinContent(i+1);

    Resb_errMax[i] = -99999;
    Resb_errMin[i] = 990009;

    for (int j(0);j<7;j++)
    {
      tmpVal  = hWmWpratio_Resbos[j]->GetBinContent(i+1);
      tmpDiff = tmpVal - nomVal;
      if( tmpDiff > Resb_errMax[i]) Resb_errMax[i] = tmpDiff;
      if( tmpDiff < Resb_errMin[i]) Resb_errMin[i] = tmpDiff;
    }
//    cout << "Resbos center ratio : " << hWmWpratio_Resbos[1]->GetBinContent(i+1) << "\t error+ : " << Resb_errMax[i] << "\t error - : " << Resb_errMin[i] << endl;

    if (Resb_errMax[i] < 0) Resb_errMax[i] = 0.;
    if (Resb_errMin[i] > 0) Resb_errMin[i] = 0.;
    if (Resb_errMin[i] < 0) Resb_errMin[i] = -Resb_errMin[i];

    Resb_err[i] = TMath::Max(Resb_errMax[i],Resb_errMin[i]);
    cout<<i<<" Bin Resbos ratio : " << hWmWpratio_Resbos[1]->GetBinContent(i+1) << " \t Error : " <<Resb_errMin[i]<<"\t"<<Resb_errMax[i]<<"\t"<<Resb_err[i] <<endl; 
    

    //printf("WmWpratioErr_FEWZ : %.8f\n",WmWpratioErr_FEWZ[i]);
    //printf("WmWpratioErr_Powheg : %.8f\n",WmWpratioErr_Powheg[i]);
    //printf("WmWpratioErr_Resbos : %.8f\n",WmWpratioErr_Resbos[i]);
  
  }

  cout << " \\"<<"\\" << endl;
  cout << endl;
  
  
  //// Print Powheg errors
  cout <<"=====  Print Powheg errors "<< endl;
  cout << fixed << setprecision(2) << endl;
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWmWpratioStatErr_Powheg: " << WmWpratioStatErr_Powheg[i] <<"\t in %\t"<< 100*WmWpratioStatErr_Powheg[i]/WmWpratio_Powheg[i]<<endl; 
  }
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWmWpratioPDFErr_Powheg : " << WmWpratioPDFErr_Powheg[i] << "\t in %\t"<< 100*WmWpratioPDFErr_Powheg[i]/WmWpratio_Powheg[i]<<endl;
  }
 //// Print FEWZ errors
  cout <<"=====  Print FEWZ errors "<< endl;

  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWmWpratioStatErr_FEWZ : " << WmWpratioStatErr_FEWZ[i] << "\t in %\t"<< 100*WmWpratioStatErr_FEWZ[i]/WmWpratio_Powheg[i]<<endl;
  }

  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWmWpratioPDFErr_FEWZ : " << WmWpratioPDFErr_FEWZ[i] << "\t in %\t"<< 100*WmWpratioPDFErr_FEWZ[i]/WmWpratio_Powheg[i]<<endl;
  }
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWmWpratioScaleErr_FEWZ : " << WmWpratioScaleErr_FEWZ[i] << "\t in %\t"<< 100*WmWpratioScaleErr_FEWZ[i]/WmWpratio_Powheg[i]<<endl;
  }




  cout << fixed << setprecision(5) << endl;
  cout<<"Wpt minus FEWZ\t"<<"Wpt plus FEWZ"<<endl; 
  for(int i=0;i<12;i++)
  {
   cout<<hWmToMu_FEWZ->GetBinContent(i+1) <<"\t"<<hWpToMu_FEWZ->GetBinContent(i+1)<<endl ; 
  }
  
  cout << fixed << setprecision(3) << endl;
  cout << " ============ W-/W+ ratio and errors(Numbers) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error " << endl;
  for(int i=0;i<12;i++)
  {

    cout <<i+1<< "\t" <<  WmWpratio_RD[i] <<  " \\pm " << WmWpratioErr_RD_Total[i] <<  endl;
  }
  cout << " ============ W-/W+ ratio and errors(in %) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error %" << endl;
  for(int i=0;i<12;i++)
  {
    cout <<i+1<< "\t" <<  WmWpratio_RD[i] <<  " \\pm " << WmWpratioErr_RD_Total[i]/WmWpratio_RD[i]*100 <<  endl;
  }
  
  for(int i=0;i<12;i++)
  {
    hWmWpratio_RD->SetBinContent(i+1,WmWpratio_RD[i]);
    hWmWpratio_RD->SetBinError(i+1,WmWpratioErr_RD_Total[i]);
    
    hWmWpratio_Powheg->SetBinContent(i+1,WmWpratio_Powheg[i]);
    hWmWpratio_Powheg->SetBinError(i+1,WmWpratioTotalErr_Powheg[i]);
    
    hWmWpratio_Resbos_Central->SetBinContent(i+1,hWmWpratio_Resbos[1]->GetBinContent(i+1));
    hWmWpratio_Resbos_Central->SetBinError(i+1,Resb_err[i]);
    //hWmWpratio_Resbos[1]->SetBinError(i+1,Resb_err[i]);
    
    hWmWpratio_FEWZ->SetBinContent(i+1,WmWpratio_FEWZ[i]);
    hWmWpratio_FEWZ->SetBinError(i+1,WmWpratioTotalErr_FEWZ[i]);
  }

  // Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataTotalErr = new TH1D("hRatioDataTotalErr","hRatioDataTotalErr",nBins-1,WptBins);

  TH1D *hRatioResbosTotalErr = new TH1D("hRatioResbosTotalErr","hRatioResbosTotalErr",nBins-1,WptBins);

  TH1D *hRatioPowhegStatErr = new TH1D("hRatioPowhegStatErr","hRatioPowhegStatErr",nBins-1,WptBins);
  TH1D *hRatioPowhegPDFErr = new TH1D("hRatioPowhegPDFErr","hRatioPowhegPDFErr",nBins-1,WptBins);
  TH1D *hRatioPowhegTotalErr = new TH1D("hRatioPowhegTotalErr","hRatioPowhegTotalErr",nBins-1,WptBins);
  
  TH1D *hRatioFEWZStatErr = new TH1D("hRatioFEWZStatErr","hRatioFEWZStatErr",nBins-1,WptBins);
  TH1D *hRatioFEWZPDFErr = new TH1D("hRatioFEWZPDFErr","hRatioFEWZPDFErr",nBins-1,WptBins);
  TH1D *hRatioFEWZScaleErr = new TH1D("hRatioFEWZScaleErr","hRatioFEWZScaleErr",nBins-1,WptBins);
  TH1D *hRatioFEWZTotalErr = new TH1D("hRatioFEWZTotalErr","hRatioFEWZTotalErr",nBins-1,WptBins);

  for(int i(0); i<nBins-1; i++)
  {
    hRatioDataTotalErr->SetBinContent(i+1,1.);
    hRatioDataTotalErr->SetBinError(i+1,WmWpratioErr_RD_Total[i] / WmWpratio_RD[i]);
    
    //hRatioResbosTotalErr->SetBinContent(i+1,hWmWpratio_Resbos[1]->GetBinContent(i+1) / WmWpratio_RD[i]);
    //hRatioResbosTotalErr->SetBinError(i+1,hWmWpratio_Resbos[1]->GetBinError(i+1) / WmWpratio_RD[i]);
    hRatioResbosTotalErr->SetBinContent(i+1,hWmWpratio_Resbos_Central->GetBinContent(i+1) / WmWpratio_RD[i]);
    hRatioResbosTotalErr->SetBinError(i+1,hWmWpratio_Resbos_Central->GetBinError(i+1) / WmWpratio_RD[i]);
    
    hRatioPowhegStatErr->SetBinContent(i+1,WmWpratio_Powheg[i] / WmWpratio_RD[i]);
    hRatioPowhegStatErr->SetBinError(i+1,WmWpratioStatErr_Powheg[i] / WmWpratio_RD[i]);
    hRatioPowhegPDFErr->SetBinContent(i+1,WmWpratio_Powheg[i] / WmWpratio_RD[i]);
    hRatioPowhegPDFErr->SetBinError(i+1,(WmWpratioPDFErr_Powheg[i]+WmWpratioStatErr_Powheg[i]) / WmWpratio_RD[i]);
    
    hRatioFEWZStatErr->SetBinContent(i+1,WmWpratio_FEWZ[i] / WmWpratio_RD[i]);
    hRatioFEWZStatErr->SetBinError(i+1,WmWpratioStatErr_FEWZ[i] / WmWpratio_RD[i]);
    hRatioFEWZPDFErr->SetBinContent(i+1,WmWpratio_FEWZ[i] / WmWpratio_RD[i]);
    hRatioFEWZPDFErr->SetBinError(i+1,(WmWpratioStatErr_FEWZ[i]+WmWpratioPDFErr_FEWZ[i]) / WmWpratio_RD[i]);
    hRatioFEWZScaleErr->SetBinContent(i+1,WmWpratio_FEWZ[i] / WmWpratio_RD[i]);
    hRatioFEWZScaleErr->SetBinError(i+1,(WmWpratioStatErr_FEWZ[i]+WmWpratioPDFErr_FEWZ[i]+WmWpratioScaleErr_FEWZ[i]) / WmWpratio_RD[i]);
    hRatioFEWZTotalErr->SetBinContent(i+1,WmWpratio_FEWZ[i] / WmWpratio_RD[i]);
    hRatioFEWZTotalErr->SetBinError(i+1,WmWpratioTotalErr_FEWZ[i] / WmWpratio_RD[i]);
  }

  //Define TGraph
  TGraphErrors *tgData = new TGraphErrors(hWmWpratio_RD);
  TGraphErrors *tgDataStatSystRatioBand = new TGraphErrors(hRatioDataTotalErr);
  
  TGraphErrors *tgResbos = new TGraphErrors(hWmWpratio_Resbos_Central);
  //TGraphErrors *tgResbosPDFRatioBand = new TGraphErrors(hWmWpratio_Resbos_Central);
  TGraphErrors *tgResbosScaleRatioBand = new TGraphErrors(hRatioResbosTotalErr);
  
  TGraphErrors *tgPowheg = new TGraphErrors(hWmWpratio_Powheg);
  TGraphErrors *tgPowhegStatRatioBand = new TGraphErrors(hRatioPowhegStatErr);
  TGraphErrors *tgPowhegPDFRatioBand = new TGraphErrors(hRatioPowhegPDFErr);

  TGraphErrors *tgFEWZ = new TGraphErrors(hWmWpratio_FEWZ);
  TGraphErrors *tgFEWZStatRatioBand = new TGraphErrors(hRatioFEWZStatErr);
  TGraphErrors *tgFEWZPDFRatioBand = new TGraphErrors(hRatioFEWZPDFErr);
  TGraphErrors *tgFEWZScaleRatioBand = new TGraphErrors(hRatioFEWZScaleErr);

  //// Now design and Draw 
  gStyle->SetLineWidth(2.);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesSpacing(0.75);
  gStyle->SetHatchesLineWidth(2);

  // Marker Color
  Color_t MarkerColor_ResBos = kBlue+2;
  Color_t MarkerColor_POWHEG = kRed+3;
  Color_t MarkerColor_FEWZ = kGreen+4;

  // Marker Style (circle:20  square:21  triangle:22)
  Style_t MarkerStyle_ResBos = 21;
  Style_t MarkerStyle_POWHEG = 21;
  Style_t MarkerStyle_FEWZ = 21;

  Color_t BandColor_ResBos_scale = kBlue;
  Color_t BandColor_ResBos_PDF = kMagenta-9;
  Color_t BandColor_POWHEG_stat = kRed+2;
  Color_t BandColor_POWHEG_PDF = kYellow;
  Color_t BandColor_FEWZ_stat = kGreen+3;
  Color_t BandColor_FEWZ_PDF = kGreen+1;
  Color_t BandColor_FEWZ_scale = kCyan-9;
  
  // Band Color Transparent
  TColor *colResbos = gROOT->GetColor(kBlue);				// Used ResBos distribution and scale ratio band
  TColor *colResbosScale = gROOT->GetColor(BandColor_ResBos_scale);	// Used ResBos distribution and scale ratio band
  TColor *colResbosPDF = gROOT->GetColor(BandColor_ResBos_PDF);		// Used ResBos PDF ratio band
  TColor *colPowheg = gROOT->GetColor(kRed);				// Used Powheg distribution
  TColor *colPowhegStat = gROOT->GetColor(BandColor_POWHEG_stat);	// used Powheg ratio band
  TColor *colPowhegPDF = gROOT->GetColor(BandColor_POWHEG_PDF);		// used Powheg ratio band
  TColor *colFEWZ = gROOT->GetColor(kGreen);				// used FEWZ distribution
  TColor *colFEWZStat = gROOT->GetColor(BandColor_FEWZ_stat);		// used FEWZ stat ratio band
  TColor *colFEWZPDF = gROOT->GetColor(BandColor_FEWZ_PDF);		// used FEWZ PDF ratio band
  TColor *colFEWZScale = gROOT->GetColor(BandColor_FEWZ_scale);		// used FEWZ scale ratio band
  colResbosScale->SetAlpha(0.4);
  colResbosPDF->SetAlpha(0.4);
  colPowheg->SetAlpha(0.3);
  colPowhegStat->SetAlpha(0.6);
  colPowhegPDF->SetAlpha(0.6);
  colFEWZ->SetAlpha(0.5); 
  colFEWZStat->SetAlpha(0.6); 
  colFEWZPDF->SetAlpha(0.7); 
  colFEWZScale->SetAlpha(0.7); 


  TLegend *lL =new TLegend(0.25,0.17,0.65,0.45); lL->SetFillColor(0); lL->SetBorderSize(0);
  lL->AddEntry(tgData,"Data","PLE1");
  lL->AddEntry(tgResbos,"ResBos-P CT10 NNLL","f");
  lL->AddEntry(tgPowheg,"POWHEG CT10 NLO","f");
  lL->AddEntry(tgFEWZ,"FEWZ CT10 NNLO","f");

  TPaveText *channel = new TPaveText(0.27,0.48,0.67,0.54,"NDC");
  channel->SetBorderSize(0);
  channel->SetFillStyle(0);
  channel->AddText("W^{#font[122]{-}} #rightarrow #mu^{#font[122]{-}}#bar{#nu}_{#mu} / W^{+} #rightarrow #mu^{+}#nu_{#mu}");

  // Canvas for distribution
  TCanvas *lC1 = new TCanvas("Can","Can",50,50,W,H);
  lC1->SetFillColor(0);
  lC1->SetBorderMode(0);
  lC1->SetFrameFillStyle(0);
  lC1->SetFrameBorderMode(0);
  //lC1->SetLeftMargin( L/W );
  //lC1->SetLeftMargin( 0.15 );
  lC1->SetLeftMargin( 0.20 );
  lC1->SetRightMargin( R/W );
  lC1->SetTopMargin( T/H );
  lC1->SetBottomMargin( B/H );
  lC1->SetTickx(0);
  lC1->SetTicky(0);
  lC1->SetLogx();
  //lC1->SetLogy();

  // Frame setting 
  tgPowheg->GetYaxis()->SetRangeUser(0.,1.4);
  //tgPowheg->SetMinimum(1e-7);
  //tgPowheg->SetMaximum(5*TMath::MaxElement(n,tgPowheg->GetY()));
  tgPowheg->SetTitle("");
  tgPowheg->GetYaxis()->SetTitle("#left(#frac{1}{#lower[0.30]{#sigma}} #frac{d#sigma}{#lower[0.10]{p_{T}^{W#font[122]{\55}}}}#right) #lower[0.5]{#scale[2]{/}}#left(#frac{1}{#lower[0.30]{#sigma}} #frac{d#sigma}{#lower[0.07]{p_{T}^{W+}}}#right)");
  tgPowheg->GetYaxis()->SetTitleOffset(2.2);
  tgPowheg->GetYaxis()->SetTitleSize(0.04);
  tgPowheg->GetYaxis()->SetLabelSize(0.04);

  tgPowheg->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  tgPowheg->GetXaxis()->SetTitleSize(0.04);
  tgPowheg->GetXaxis()->SetTitleOffset(0.55);
  tgPowheg->GetXaxis()->SetLabelSize(0.04);

  //Data

  //Powheg
  tgPowheg->SetFillColor(kRed);
  tgFEWZ->SetFillColor(kGreen);
  tgResbos->SetFillColor(kBlue);

  //Draw Original Diff-Xsec Distribution
  tgPowheg->Draw("A2");
  tgResbos->Draw("2 same");
  tgFEWZ->Draw("2 same");
  tgData->Draw("P same");

  lL->Draw();
  channel->Draw();
  CMS_lumi(lC1,iPeriod,iPos);
  lC1->Update();
  lC1->RedrawAxis();
  lC1->GetFrame()->Draw();

  sprintf(tmpName,"WmMuWpMuNormFid12Bin.");
  lC1->SaveAs(tmpName+format);

  // Ratio plot style 
  tgDataStatSystRatioBand->SetMarkerStyle(20);
  tgDataStatSystRatioBand->SetMarkerColor(kBlack);
  tgDataStatSystRatioBand->SetMarkerSize(0.7);
  tgDataStatSystRatioBand->SetLineWidth(2.0);
  tgDataStatSystRatioBand->SetLineColor(kBlack);

  tgDataStatSystRatioBand->SetFillStyle(3354);
  tgDataStatSystRatioBand->SetFillColor(kGray+1);

  // Resbos Ratio plot style
  //tgResbosPDFRatioBand->SetMarkerColor(MarkerColor_ResBos);
  //tgResbosPDFRatioBand->SetMarkerStyle(MarkerStyle_ResBos);
  //tgResbosPDFRatioBand->SetFillColor(BandColor_ResBos_PDF);

  tgResbosScaleRatioBand->SetMarkerColor(MarkerColor_ResBos);
  tgResbosScaleRatioBand->SetMarkerStyle(MarkerStyle_ResBos);
  tgResbosScaleRatioBand->SetFillColor(BandColor_ResBos_scale);

  // Powheg Ratio plot style
  tgPowhegStatRatioBand->SetMarkerColor(MarkerColor_POWHEG);
  tgPowhegStatRatioBand->SetMarkerStyle(MarkerStyle_POWHEG);
  tgPowhegStatRatioBand->SetFillColor(BandColor_POWHEG_stat);

  tgPowhegPDFRatioBand->SetMarkerColor(MarkerColor_POWHEG);
  tgPowhegPDFRatioBand->SetMarkerStyle(MarkerStyle_POWHEG);
  tgPowhegPDFRatioBand->SetFillColor(BandColor_POWHEG_PDF);

  // FEWZ Ratio plot style
  tgFEWZStatRatioBand->SetMarkerColor(MarkerColor_FEWZ);
  tgFEWZStatRatioBand->SetMarkerStyle(MarkerStyle_FEWZ);
  tgFEWZStatRatioBand->SetFillColor(BandColor_FEWZ_stat);

  tgFEWZPDFRatioBand->SetMarkerColor(MarkerColor_FEWZ);
  tgFEWZPDFRatioBand->SetMarkerStyle(MarkerStyle_FEWZ);
  tgFEWZPDFRatioBand->SetFillColor(BandColor_FEWZ_PDF);

  tgFEWZScaleRatioBand->SetMarkerColor(MarkerColor_FEWZ);
  tgFEWZScaleRatioBand->SetMarkerStyle(MarkerStyle_FEWZ);
  tgFEWZScaleRatioBand->SetFillColor(BandColor_FEWZ_scale);

  // Canvas for Theory/data Ratio
  TCanvas *lC2 = new TCanvas("Can","Can",50,50,W,H); 
  lC2->SetFillColor(0);
  lC2->SetBorderMode(0);
  lC2->SetFrameFillStyle(0);
  lC2->SetFrameBorderMode(0);
  lC2->SetLeftMargin( L/W );
  lC2->SetRightMargin( R/W );
  lC2->SetTopMargin( T/H );
  lC2->SetBottomMargin( B/H );
  CMS_lumi(lC2,iPeriod,iPos);

  lC2->Divide(1,3,0,0);
  lC2->cd(1)->SetPad(0,0.66,0.96,0.945);
  //lC2->cd(1)->SetFillColor(2);
  lC2->cd(1)->SetTickx(1);
  lC2->cd(1)->SetTicky(1);
  lC2->cd(1)->SetLogx(1);
  
  TLegend *lResbos =new TLegend(0.18,0.15,0.60,0.30); lResbos->SetFillColor(0); lResbos->SetBorderSize(0);
  lResbos-> SetNColumns(2);
  lResbos->AddEntry(tgResbosScaleRatioBand,"ResBos scales","FP");
  lResbos->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  //lResbos->AddEntry(tgResbosPDFRatioBand,"ResBos PDF","F");
  lResbos->SetTextSize(0.08);

  TPaveText *tResBos = new TPaveText(0.18,0.82,0.34,0.92,"NDC");
  tResBos->SetBorderSize(0);
  tResBos->SetFillStyle(0);
  tResBos->SetTextSize(0.12);
  tResBos->AddText("#font[42]{ResBos-P}");
 
  TPaveText *tChannel = new TPaveText(0.40,0.82,0.72,0.93,"NDC");
  tChannel->SetBorderSize(0);
  tChannel->SetFillStyle(0);
  tChannel->SetTextSize(0.12);
  tChannel->AddText("W^{#font[122]{-}} #rightarrow #mu^{#font[122]{-}}#bar{#nu}_{#mu} / W^{+} #rightarrow #mu^{+}#nu_{#mu}");

  tgDataStatSystRatioBand->GetYaxis()->SetRangeUser(0.5,1.5);
  tgDataStatSystRatioBand->GetYaxis()->SetTitle("Theory/Data");
  tgDataStatSystRatioBand->GetYaxis()->CenterTitle();
  tgDataStatSystRatioBand->GetYaxis()->SetTitleOffset(1.3);
  tgDataStatSystRatioBand->GetYaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetYaxis()->SetTitleSize(33);
  tgDataStatSystRatioBand->GetYaxis()->SetLabelFont(43);
  tgDataStatSystRatioBand->GetYaxis()->SetLabelSize(29);
  tgDataStatSystRatioBand->GetYaxis()->SetNdivisions(405);
  tgDataStatSystRatioBand->GetXaxis()->SetRangeUser(0,600);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleOffset(0.6);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleSize(20);
  tgDataStatSystRatioBand->GetXaxis()->SetLabelSize(0.11);
  tgDataStatSystRatioBand->Draw("2 A");
  //tgResbosPDFRatioBand->Draw("2 P");
  tgResbosScaleRatioBand->Draw("2 P");
  lResbos->Draw();
  tResBos->Draw();
  tChannel->Draw();

  //Powheg Ratio plot
  lC2->cd(2)->SetPad(0,0.39,0.96,0.64);
  //lC2->cd(2)->SetFillColor(4);
  lC2->cd(2)->SetTickx(1);
  lC2->cd(2)->SetTicky(1);
  lC2->cd(2)->SetLogx(1);

  TLegend *lPowheg =new TLegend(0.18,0.04,0.60,0.34); lPowheg->SetFillColor(0); lPowheg->SetBorderSize(0);
  lPowheg-> SetNColumns(2);
  lPowheg->AddEntry(tgPowhegStatRatioBand,"POWHEG stat","FP");
  lPowheg->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  lPowheg->AddEntry(tgPowhegPDFRatioBand,"POWHEG PDF    ","F");
  //lPowheg->AddEntry(hRatioDataStatErr,"","");
  lPowheg->SetTextSize(0.09);

  TPaveText *tPowheg = new TPaveText(0.18,0.82,0.34,0.92,"NDC");
  tPowheg->SetBorderSize(0);
  tPowheg->SetFillStyle(0);
  tPowheg->SetTextSize(0.12);
  tPowheg->AddText("#font[42]{POWHEG}");
  
  tgDataStatSystRatioBand->Draw("2 A");
  tgPowhegPDFRatioBand->Draw("2 P");
  tgPowhegStatRatioBand->Draw("2 P");
  lPowheg->Draw();
  tPowheg->Draw();

  // FEWZ Ratio Plot
  lC2->cd(3)->SetPad(0,0,0.96,0.37);
  lC2->cd(3)->SetBottomMargin(0.26);
  //lC2->cd(3)->SetFillColor(3);
  lC2->cd(3)->SetTickx(1);
  lC2->cd(3)->SetTicky(1);
  lC2->cd(3)->SetLogx(1);

  TLegend *lFewz =new TLegend(0.18,0.28,0.60,0.52); lFewz->SetFillColor(0); lFewz->SetBorderSize(0);
  lFewz-> SetNColumns(2);
  lFewz->AddEntry(tgFEWZStatRatioBand,"FEWZ stat","FP");
  lFewz->AddEntry(tgDataStatSystRatioBand,"Data stat+syst","F");
  lFewz->AddEntry(tgFEWZPDFRatioBand,"FEWZ PDF","F");
  lFewz->AddEntry(tgFEWZScaleRatioBand,"FEWZ scales","F");
  lFewz->SetTextSize(0.06);

  TPaveText *tFewz = new TPaveText(0.18,0.82,0.34,0.92,"NDC");
  tFewz->SetBorderSize(0);
  tFewz->SetFillStyle(0);
  tFewz->SetTextSize(0.09);
  tFewz->AddText("#font[42]{FEWZ}");
 
  tgDataStatSystRatioBand->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  tgDataStatSystRatioBand->GetXaxis()->SetTitleOffset(1.5);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleFont(43);
  tgDataStatSystRatioBand->GetXaxis()->SetTitleSize(33);
  tgDataStatSystRatioBand->Draw("2 A");
  tgFEWZScaleRatioBand->Draw("2");
  tgFEWZPDFRatioBand->Draw("2");
  tgFEWZStatRatioBand->Draw("2P");
  lFewz->Draw();
  tFewz->Draw();

  sprintf(tmpName,"WmWpTheoryData.");
  lC2->SaveAs(tmpName+format);
  
  return 0;
}
 
void FEWZ_PDFUncer(double *WmWpRatio, double *Err)
{
  // FEWZ Normalized PDF error in % unit
  double PDFErr[12] = {0.,};

  PDFErr[0] = 1.14;  
  PDFErr[1] = 0.36;
  PDFErr[2] = 0.59;
  PDFErr[3] = 0.58;
  PDFErr[4] = 0.80;
  PDFErr[5] = 1.06;
  PDFErr[6] = 1.23;
  PDFErr[7] = 1.58;
  PDFErr[8] = 1.80;
  PDFErr[9] = 2.52;
  PDFErr[10] = 2.24;
  PDFErr[11] = 3.87;

  for(int i(0); i<12; i++)
  {
    Err[i] = WmWpRatio[i] * PDFErr[i] * 0.01;
    //cout << "PDFErr in function : " << Err[i] << endl;
  }
}
void Powheg_PDFUncer(double *WmWpRatio, double *Err)
{
  // Powheg Normalized PDF error in % unit
  double PDFErr[12] = {0.,};

  PDFErr[0] =0.34912  ;  
  PDFErr[1] =0.11815  ;
  PDFErr[2] =0.12652  ;
  PDFErr[3] =0.27171  ;
  PDFErr[4] =0.33522  ;
  PDFErr[5] =0.50769  ;
  PDFErr[6] =0.61623  ;
  PDFErr[7] =0.77008  ;
  PDFErr[8] =0.87044  ;
  PDFErr[9] =0.89737  ;
  PDFErr[10]=1.14706  ;
  PDFErr[11]=1.11763  ;

  for(int i(0); i<12; i++)
  {
    Err[i] = WmWpRatio[i] * PDFErr[i] * 0.01;
    //cout << "PDFErr in function : " << Err[i] << endl;
  }
}
