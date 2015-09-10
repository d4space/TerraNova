#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "../../Utils/const.h"

double BinWidth[14] ={0.0, 7.5-0, 12.5-7.5,17.5-12.5, 24.0-17.5, 30.0-24.0, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};
const int nBins = 14;
double WptLogBins[nBins] = {1.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};
double WptBins[nBins] = {0.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};

double ax[13]  = {4.25,10,15,20.75,27,35,45,60,90,130,170,220,425};
double aex[13] = {3.25,2.5,2.5,3.25,3,5,5,10,20,20,20,30,175};

const int n12Bins = 13;
double WptLog12Bins[n12Bins] = {1.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double Wpt12Bins[n12Bins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

double ax12[12]  = {4.25,10,15,23.75,35,45,60,90,130,170,220,425};
double aex12[12] = {3.25,2.5,2.5,6.25,5,5,10,20,20,20,30,175};

void drawDifference(TH1* iH0,TH1 *iH1,TH1* iH2, TGraphErrors* iH3, int chnl,TGraphErrors* iH4,TGraphAsymmErrors* iH5,TH1* StatErrBand,TGraphErrors* iH6){
  std::string lName = std::string(iH0->GetName());
  //TH1F *lHDiff  = new TH1F((lName+"Diff").c_str(),(lName+"Diff").c_str(),nBins-1,WptLogBins);// lHDiff->Sumw2();
  //12 Bin
  TH1F *lHDiff  = new TH1F((lName+"Diff").c_str(),(lName+"Diff").c_str(),n12Bins-1,WptLog12Bins);// lHDiff->Sumw2();
  
  TH1F *lXHDiff1 = new TH1F((lName+"XDiff1").c_str(),(lName+"XDiff1").c_str(),iH0->GetNbinsX(),iH0->GetXaxis()->GetXmin(),iH0->GetXaxis()->GetXmax());
  int i1 = 0;
  lXHDiff1->SetLineWidth(2); lXHDiff1->SetLineColor(kBlack); //lXHDiff1->SetLineStyle(2);
  
  StatErrBand->SetMarkerStyle(kFullCircle); StatErrBand->SetMarkerColor(kBlack);StatErrBand->SetMarkerSize(0.6);
  StatErrBand->SetLineWidth(2); StatErrBand->SetLineColor(kBlack);
  
  lHDiff->GetYaxis()->SetRangeUser(0.5,1.5); //Eleminus
  if(chnl==2)
    lHDiff->GetYaxis()->SetRangeUser(0.5,1.5);//Eleminus
  if(chnl==3)
    lHDiff->GetYaxis()->SetRangeUser(0.5,1.5);
  lHDiff->GetYaxis()->SetTitleOffset(0.4);
  lHDiff->GetYaxis()->SetTitleSize(0.12);
  lHDiff->GetYaxis()->SetLabelSize(0.12);
  lHDiff->GetYaxis()->CenterTitle();
  lHDiff->GetXaxis()->SetTitleOffset(1.0);
  lHDiff->GetXaxis()->SetTitleSize(0.12);
  lHDiff->GetXaxis()->SetLabelSize(0);
  if(chnl==3)
    lHDiff->GetXaxis()->SetLabelSize(0.12);
  lHDiff->GetYaxis()->SetNdivisions(405);
  if(chnl==3)
    lHDiff->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  lHDiff->GetYaxis()->SetTitle("Theory/Data");
  gStyle->SetOptStat(0);
  
  for(int i0 = 0; i0 < lHDiff->GetNbinsX()+1; i0++) {
    double lXCenter = lHDiff->GetBinCenter(i0);
    double lXVal     = iH0   ->GetBinContent(i0);
    lXHDiff1->SetBinContent(i0, 1.0);
    while(iH1->GetBinCenter(i1) < lXCenter) {i1++;}
    if(iH1->GetBinContent(i0) > 0) lHDiff->SetBinContent(i0,lXVal/(iH1->GetBinContent(i0)));
    if(iH1->GetBinContent(i0) > 0) lHDiff->SetBinError(i0,0.00001);
  }
  
  TGraphErrors* ErrBand = new TGraphErrors(iH2);
  ErrBand->SetFillColor(kBlack);
  ErrBand->SetFillStyle(3354);
  ErrBand->SetLineWidth(1);
  
  if (chnl == 1)
  {
    lHDiff->SetMarkerStyle(kOpenCircle);
    lHDiff->SetMarkerColor(kBlue);
    lHDiff->SetLineColor(kBlue);
  }
  if (chnl == 2)
  {
    lHDiff->SetMarkerStyle(kOpenTriangleUp);
    lHDiff->SetMarkerColor(kRed);
    lHDiff->SetLineColor(kRed);
  }
  if (chnl == 3)
  {
    lHDiff->SetMarkerStyle(kOpenSquare);
    lHDiff->SetMarkerColor(kGreen+3);
    lHDiff->SetLineColor(kGreen+3);
  }
  
  lHDiff->SetMarkerSize(0.8);


  lHDiff->SetTitle("");
  lHDiff->Draw("E");
  if (chnl == 1)
    iH5->Draw("2same");
  if (chnl == 2 || chnl == 3)
    iH4->Draw("2");
  if (chnl == 3)
    iH6->Draw("2");
  if (chnl == 2 || chnl == 3)
    iH3->Draw("2same");
  lXHDiff1->Draw("histsame");
  ErrBand->Draw("2same");
  lHDiff->Draw("Esame");
  StatErrBand->Draw("E1same");
}


int theoryStudy_separate(const TString BaseName)
{
  TString tmpTStr;
  char tmpName[30],tmpName_org[30];
  int Numb;

  TFile *f_Resbos;
  TFile *f_Resbos_12;
  
  TFile *f_Fewz;
  TFile *f_Fewz_12;
  
  TFile *f_Data;
  TFile *f_Data_PowhegErr;

  f_Resbos = new TFile("../../RstResbos/Resbos_"+BaseName+".root");
  f_Data = new TFile("../Result"+BaseName+"/Result_"+BaseName+".root");
  //f_Data_PowhegErr = new TFile("../../Systematics/WptXsecErrors/"+BaseName+"Errors.root");
  f_Data_PowhegErr = new TFile("../../Systematics_NormDiffXsec/WptXsecErrors/"+BaseName+"Errors.root");

  if (BaseName=="WpToMuNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WpToMuNu_13bin_dynamic_NNLO.root");
    f_Fewz_12 = new TFile("../../RstFEWZ_12Bin_Fiducial_dynamic/WpToMuNu_dynamic_NNLO.root");
    f_Resbos_12 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wplus.root");
  }
  if (BaseName=="WmToMuNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WmToMuNu_13bin_dynamic_NNLO.root");
    f_Fewz_12 = new TFile("../../RstFEWZ_12Bin_Fiducial_dynamic/WmToMuNu_dynamic_NNLO.root");
    f_Resbos_12 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wminus.root");
  }
  if (BaseName=="WpToEleNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WpToEleNu_13bin_dynamic_NNLO.root");
    f_Resbos_12 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wplus.root");
  }
  if (BaseName=="WmToEleNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WmToEleNu_13bin_dynamic_NNLO.root");
    f_Resbos_12 = new TFile("../../RstResbos_12BinFidVol/Resbos_Wminus.root");
  }

  TFile f_out("Wpt_"+BaseName+"NormDiffXsec_InFid_RDResBosPowhegFEWZ.root","recreate");

  TH1D *hWptBins_LinScale   = new TH1D("hWptBins_LinScale","hWptBins_LinScale",nBins-1,WptBins);hWptBins_LinScale->Sumw2();
  TH1D *hWpt12Bins_LinScale   = new TH1D("hWpt12Bins_LinScale","hWpt12Bins_LinScale",n12Bins-1,Wpt12Bins);hWpt12Bins_LinScale->Sumw2();
  
  
  
  /// Reading Yields and making diffXsec to merge the syst errors in %
  TH1D* hRD_Yield;
  hRD_Yield   = (TH1D*)f_Data->Get("BornEffCorr")->Clone("hData_Yield_BornLinScale");
  double RD_Yield[nBins]={0};
  for( int ipt(1);ipt<nBins;ipt++)
  {
    RD_Yield[ipt]=hRD_Yield->GetBinContent(ipt)/BinWidth[ipt]/18.4;
    cout<<ipt<<"\tRD_Yield \t"<<RD_Yield[ipt]<<"\t\t"<<BinWidth[ipt]<<endl;
  }
  TH1D* hPowheg_Yield;
  hPowheg_Yield   = (TH1D*)f_Data->Get("SVD_Born.Gen")->Clone("hPowheg_Xsec_BornLinScale");
  double Powheg_Yield[nBins]={0};
  for( int ipt(1);ipt<nBins;ipt++)
  {
    Powheg_Yield[ipt]=hPowheg_Yield->GetBinContent(ipt)/BinWidth[ipt]/18.4;
    cout<<ipt<<"\tPowheg_Yield \t"<<Powheg_Yield[ipt]<<"\t\t"<<BinWidth[ipt]<<endl;
  }


  ///=============Reading related stat, syst errors starts here==============
  //  
  TH1D* hData_StatErr		 = (TH1D*)f_Data_PowhegErr->Get("h_Stat")->Clone("hData_StatErr");
  TH1D *hData_StatErr12   = new TH1D("hData_StatErr12","hData_StatErr12",n12Bins-1,Wpt12Bins);hData_StatErr12->Sumw2();
  TH1D* hData_MetResolSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_met")->Clone("hData_MetResolSystErr");
  TH1D *hData_MetResolSystErr12   = new TH1D("hData_MetResolSystErr12","hData_MetResolSystErr12",n12Bins-1,Wpt12Bins);hData_MetResolSystErr12->Sumw2();
  TH1D* hData_EnMomScaleSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_scale")->Clone("hData_EnMomScaleSystErr");
  TH1D *hData_EnMomScaleSystErr12   = new TH1D("hData_EnMomScaleSystErr12","hData_EnMomScaleSystErr12",n12Bins-1,Wpt12Bins);hData_EnMomScaleSystErr12->Sumw2();
  TH1D* hData_EnMomSmearSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_smear")->Clone("hData_EnMomSmearSystErr");
  TH1D *hData_EnMomSmearSystErr12   = new TH1D("hData_EnMomSmearSystErr12","hData_EnMomSmearSystErr12",n12Bins-1,Wpt12Bins);hData_EnMomSmearSystErr12->Sumw2();
  TH1D* hData_QcdBckgrSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_qcdbckgr")->Clone("hData_QcdBckgrSystErr");
  TH1D *hData_QcdBckgrSystErr12   = new TH1D("hData_QcdBckgrSystErr12","hData_QcdBckgrSystErr12",n12Bins-1,Wpt12Bins);hData_QcdBckgrSystErr12->Sumw2();
  TH1D* hData_QcdShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_qcdshape")->Clone("hData_QcdShapeSystErr");
  TH1D *hData_QcdShapeSystErr12   = new TH1D("hData_QcdShapeSystErr12","hData_QcdShapeSystErr12",n12Bins-1,Wpt12Bins);hData_QcdShapeSystErr12->Sumw2();
  TH1D* hData_EwkSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_ewk")->Clone("hData_EwkSystErr");
  TH1D *hData_EwkSystErr12   = new TH1D("hData_EwkSystErr12","hData_EwkSystErr12",n12Bins-1,Wpt12Bins);hData_EwkSystErr12->Sumw2();
  TH1D* hData_FsrSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_fsr")->Clone("hData_FsrSystErr");
  TH1D *hData_FsrSystErr12   = new TH1D("hData_FsrSystErr12","hData_FsrSystErr12",n12Bins-1,Wpt12Bins);hData_FsrSystErr12->Sumw2();
  TH1D* hData_SvdUnfSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_SvdUnf")->Clone("hData_SvdUnfSystErr");
  TH1D *hData_SvdUnfSystErr12   = new TH1D("hData_SvdUnfSystErr12","hData_SvdUnfSystErr12",n12Bins-1,Wpt12Bins);hData_SvdUnfSystErr12->Sumw2();
  TH1D* hData_UnfoldBiasSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_UnfoldBias")->Clone("hData_UnfoldBiasSystErr");
  TH1D *hData_UnfoldBiasSystErr12   = new TH1D("hData_UnfoldBiasSystErr12","hData_UnfoldBiasSystErr12",n12Bins-1,Wpt12Bins);hData_UnfoldBiasSystErr12->Sumw2();
  
  ///Lepton Reconstruction related systematic errors
  //
  TH1D* hData_EffiToySystErr		 = (TH1D*)f_Data_PowhegErr->Get("h_toy")->Clone("hData_EffiToySystErr");
  TH1D *hData_EffiToySystErr12   = new TH1D("hData_EffiToySystErr12","hData_EffiToySystErr12",n12Bins-1,Wpt12Bins);hData_EffiToySystErr12->Sumw2();
  TH1D* hData_IDIsoSigShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_idisosig")->Clone("hData_IDIsoSigShapeSystErr");
  TH1D *hData_IDIsoSigShapeSystErr12   = new TH1D("hData_IDIsoSigShapeSystErr12","hData_IDIsoSigShapeSystErr12",n12Bins-1,Wpt12Bins);hData_IDIsoSigShapeSystErr12->Sumw2();
  TH1D* hData_IDIsoBkgrShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_idisobck")->Clone("hData_IDIsoBkgrShapeSystErr");
  TH1D *hData_IDIsoBkgrShapeSystErr12   = new TH1D("hData_IDIsoBkgrShapeSystErr12","hData_IDIsoBkgrShapeSystErr12",n12Bins-1,Wpt12Bins);hData_IDIsoBkgrShapeSystErr12->Sumw2();
  
  for( int ipt(1);ipt<n12Bins;ipt++)
  { 
    if(ipt<4)
    {
      hData_StatErr12->SetBinContent(ipt, hData_StatErr->GetBinContent(ipt) );
      hData_MetResolSystErr12->SetBinContent(ipt, hData_MetResolSystErr->GetBinContent(ipt) );
      hData_EnMomScaleSystErr12->SetBinContent(ipt, hData_EnMomScaleSystErr->GetBinContent(ipt) );
      hData_EnMomSmearSystErr12->SetBinContent(ipt, hData_EnMomSmearSystErr->GetBinContent(ipt) );
      hData_QcdBckgrSystErr12->SetBinContent(ipt, hData_QcdBckgrSystErr->GetBinContent(ipt) );
      hData_QcdShapeSystErr12->SetBinContent(ipt, hData_QcdShapeSystErr->GetBinContent(ipt) );
      hData_EwkSystErr12->SetBinContent(ipt, hData_EwkSystErr->GetBinContent(ipt) );
      hData_FsrSystErr12->SetBinContent(ipt, hData_FsrSystErr->GetBinContent(ipt) );
      hData_SvdUnfSystErr12->SetBinContent(ipt, hData_SvdUnfSystErr->GetBinContent(ipt) );
      hData_UnfoldBiasSystErr12->SetBinContent(ipt, hData_UnfoldBiasSystErr->GetBinContent(ipt) );
      
      ///12 bin Lepton Reconstruction related systematic errors
      //
      hData_EffiToySystErr12->SetBinContent(ipt, hData_EffiToySystErr->GetBinContent(ipt) );
      hData_IDIsoSigShapeSystErr12->SetBinContent(ipt, hData_IDIsoSigShapeSystErr->GetBinContent(ipt) );
      hData_IDIsoBkgrShapeSystErr12->SetBinContent(ipt, hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt) );
    }else if(ipt==4)
    {
      hData_StatErr12->SetBinContent(ipt, sqrt(hData_StatErr->GetBinContent(ipt)*hData_StatErr->GetBinContent(ipt) + hData_StatErr->GetBinContent(ipt+1)*hData_StatErr->GetBinContent(ipt+1) ) );
      //hData_MetResolSystErr12->SetBinContent(ipt, TMath::Max(hData_MetResolSystErr->GetBinContent(ipt),hData_MetResolSystErr->GetBinContent(ipt+1) ) );
      //hData_EnMomScaleSystErr12->SetBinContent(ipt, TMath::Max(hData_EnMomScaleSystErr->GetBinContent(ipt),hData_EnMomScaleSystErr->GetBinContent(ipt+1) ) );
      //hData_EnMomSmearSystErr12->SetBinContent(ipt, TMath::Max(hData_EnMomSmearSystErr->GetBinContent(ipt),hData_EnMomSmearSystErr->GetBinContent(ipt+1) ) );
      //hData_QcdBckgrSystErr12->SetBinContent(ipt, TMath::Max(hData_QcdBckgrSystErr->GetBinContent(ipt),hData_QcdBckgrSystErr->GetBinContent(ipt+1) ) );
      //hData_QcdShapeSystErr12->SetBinContent(ipt, TMath::Max(hData_QcdShapeSystErr->GetBinContent(ipt),hData_QcdShapeSystErr->GetBinContent(ipt+1) ) );
      //hData_EwkSystErr12->SetBinContent(ipt, TMath::Max(hData_EwkSystErr->GetBinContent(ipt),hData_EwkSystErr->GetBinContent(ipt+1) ) );
      //hData_FsrSystErr12->SetBinContent(ipt, TMath::Max(hData_FsrSystErr->GetBinContent(ipt),hData_FsrSystErr->GetBinContent(ipt+1) ) );
      //hData_SvdUnfSystErr12->SetBinContent(ipt, TMath::Max(hData_SvdUnfSystErr->GetBinContent(ipt),hData_SvdUnfSystErr->GetBinContent(ipt+1) ) );
      //hData_UnfoldBiasSystErr12->SetBinContent(ipt, TMath::Max(hData_UnfoldBiasSystErr->GetBinContent(ipt),hData_UnfoldBiasSystErr->GetBinContent(ipt+1) ) );
      //
      ///// 12 bin Lepton Reconstruction related systematic errors
      ////
      //hData_EffiToySystErr12->SetBinContent(ipt, TMath::Max(hData_EffiToySystErr->GetBinContent(ipt),hData_EffiToySystErr->GetBinContent(ipt+1) ) );
      //hData_IDIsoSigShapeSystErr12->SetBinContent(ipt, TMath::Max(hData_IDIsoSigShapeSystErr->GetBinContent(ipt),hData_IDIsoSigShapeSystErr->GetBinContent(ipt+1) ) );
      //hData_IDIsoBkgrShapeSystErr12->SetBinContent(ipt, TMath::Max(hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt),hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt+1) ) );
      
      hData_MetResolSystErr12->SetBinContent(ipt, (hData_MetResolSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_MetResolSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5])
	  				/(BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_EnMomScaleSystErr12->SetBinContent(ipt, (hData_EnMomScaleSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_EnMomScaleSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5])
	  				/(BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_EnMomSmearSystErr12->SetBinContent(ipt, (hData_EnMomSmearSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_EnMomSmearSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5])
	  				/(BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_QcdBckgrSystErr12->SetBinContent(ipt, (hData_QcdBckgrSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_QcdBckgrSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5])
	  				/(BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_QcdShapeSystErr12->SetBinContent(ipt, (hData_QcdShapeSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_QcdShapeSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_EwkSystErr12->SetBinContent(ipt, (hData_EwkSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_EwkSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_FsrSystErr12->SetBinContent(ipt, (hData_FsrSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_FsrSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_SvdUnfSystErr12->SetBinContent(ipt, (hData_SvdUnfSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_SvdUnfSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_UnfoldBiasSystErr12->SetBinContent(ipt, (hData_UnfoldBiasSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_UnfoldBiasSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      
      /// 12 bin Lepton Reconstruction related systematic errors
      //
      hData_EffiToySystErr12->SetBinContent(ipt, (hData_EffiToySystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_EffiToySystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_IDIsoSigShapeSystErr12->SetBinContent(ipt, (hData_IDIsoSigShapeSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_IDIsoSigShapeSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
      hData_IDIsoBkgrShapeSystErr12->SetBinContent(ipt, (hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	  				/ (BinWidth[4]*RD_Yield[4]+BinWidth[5]*RD_Yield[5]) );
    }
    else if(ipt>=5)
    {
      hData_StatErr12->SetBinContent(ipt, hData_StatErr->GetBinContent(ipt+1) );
      hData_MetResolSystErr12->SetBinContent(ipt, hData_MetResolSystErr->GetBinContent(ipt+1) );
      hData_EnMomScaleSystErr12->SetBinContent(ipt, hData_EnMomScaleSystErr->GetBinContent(ipt+1) );
      hData_EnMomSmearSystErr12->SetBinContent(ipt, hData_EnMomSmearSystErr->GetBinContent(ipt+1) );
      hData_QcdBckgrSystErr12->SetBinContent(ipt, hData_QcdBckgrSystErr->GetBinContent(ipt+1) );
      hData_QcdShapeSystErr12->SetBinContent(ipt, hData_QcdShapeSystErr->GetBinContent(ipt+1) );
      hData_EwkSystErr12->SetBinContent(ipt, hData_EwkSystErr->GetBinContent(ipt+1) );
      hData_FsrSystErr12->SetBinContent(ipt, hData_FsrSystErr->GetBinContent(ipt+1) );
      hData_SvdUnfSystErr12->SetBinContent(ipt, hData_SvdUnfSystErr->GetBinContent(ipt+1) );
      hData_UnfoldBiasSystErr12->SetBinContent(ipt, hData_UnfoldBiasSystErr->GetBinContent(ipt+1) );
      
      /// 12 bin Lepton Reconstruction related systematic errors
      //
      hData_EffiToySystErr12->SetBinContent(ipt, hData_EffiToySystErr->GetBinContent(ipt+1) );
      hData_IDIsoSigShapeSystErr12->SetBinContent(ipt, hData_IDIsoSigShapeSystErr->GetBinContent(ipt+1) );
      hData_IDIsoBkgrShapeSystErr12->SetBinContent(ipt, hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt+1) );
    
    }
  } 

  
 
  
  TH1D *hData_TotalSystErr   = new TH1D("hData_TotalSystErr","hData_TotalSystErr",nBins-1,WptBins);hData_TotalSystErr->Sumw2();
  TH1D *hData_TotalSystErr12   = new TH1D("hData_TotalSystErr12","hData_TotalSystErr12",n12Bins-1,Wpt12Bins);hData_TotalSystErr12->Sumw2();
  
  if(BaseName== "WpToMuNu" || BaseName== "WmToMuNu")
  {
    TH1D* hData_TrackSigShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_tracksig")->Clone("hData_TrackSigShapeSystErr");
    TH1D *hData_TrackSigShapeSystErr12   = new TH1D("hData_TrackSigShapeSystErr12","hData_TrackSigShapeSystErr12",n12Bins-1,Wpt12Bins);hData_TrackSigShapeSystErr12->Sumw2();
    TH1D* hData_TrackBkgrShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_trackbck")->Clone("hData_TrackBkgrShapeSystErr");
    TH1D *hData_TrackBkgrShapeSystErr12   = new TH1D("hData_TrackBkgrShapeSystErr12","hData_TrackBkgrShapeSystErr12",n12Bins-1,Wpt12Bins);hData_TrackBkgrShapeSystErr12->Sumw2();
    TH1D* hData_MuonPOGSystErr		 = (TH1D*)f_Data_PowhegErr->Get("h_POG")->Clone("hData_MuonPOGSystErr");
    TH1D *hData_MuonPOGSystErr12   = new TH1D("hData_MuonPOGSystErr12","hData_MuonPOGSystErr12",n12Bins-1,Wpt12Bins);hData_MuonPOGSystErr12->Sumw2();
    for( int ipt(1);ipt<n12Bins;ipt++)
    { 
      if(ipt<4)
      {
	hData_TrackSigShapeSystErr12->SetBinContent(ipt, hData_TrackSigShapeSystErr->GetBinContent(ipt) );
	hData_TrackBkgrShapeSystErr12->SetBinContent(ipt, hData_TrackBkgrShapeSystErr->GetBinContent(ipt) );
	hData_MuonPOGSystErr12->SetBinContent(ipt, hData_MuonPOGSystErr->GetBinContent(ipt) );
      }else if(ipt==4)
      {
	//hData_TrackSigShapeSystErr12->SetBinContent(ipt, TMath::Max(hData_TrackSigShapeSystErr->GetBinContent(ipt),hData_TrackSigShapeSystErr->GetBinContent(ipt+1) ) );
	//hData_TrackBkgrShapeSystErr12->SetBinContent(ipt, TMath::Max(hData_TrackBkgrShapeSystErr->GetBinContent(ipt),hData_TrackBkgrShapeSystErr->GetBinContent(ipt+1) ) );
	//hData_MuonPOGSystErr12->SetBinContent(ipt, TMath::Max(hData_MuonPOGSystErr->GetBinContent(ipt),hData_MuonPOGSystErr->GetBinContent(ipt+1) ) );
	hData_TrackSigShapeSystErr12->SetBinContent(ipt, (hData_TrackSigShapeSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4]+ hData_TrackSigShapeSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5])  
	    						/(BinWidth[4]*RD_Yield[4] + BinWidth[5]*RD_Yield[5]) );
	hData_TrackBkgrShapeSystErr12->SetBinContent(ipt, (hData_TrackBkgrShapeSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_TrackBkgrShapeSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	    					/ (BinWidth[4]*RD_Yield[4] + BinWidth[5]*RD_Yield[5]) );
	hData_MuonPOGSystErr12->SetBinContent(ipt, (hData_MuonPOGSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_MuonPOGSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	    					/ (BinWidth[4]*RD_Yield[4] + BinWidth[5]*RD_Yield[5]) );
      }else if(ipt>=5)
      {
	hData_TrackSigShapeSystErr12->SetBinContent(ipt, hData_TrackSigShapeSystErr->GetBinContent(ipt+1) );
	hData_TrackBkgrShapeSystErr12->SetBinContent(ipt, hData_TrackBkgrShapeSystErr->GetBinContent(ipt+1) );
	hData_MuonPOGSystErr12->SetBinContent(ipt, hData_MuonPOGSystErr->GetBinContent(ipt+1) );
      }
    } 
   
    //// Calculate total syst for muon
    for( int ipt(1);ipt<=nBins-1;ipt++)
    {

      ///used for Normalized Xsec
      hData_TotalSystErr->SetBinContent(ipt, sqrt( 
	    TMath::Power(hData_TrackSigShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_TrackBkgrShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MuonPOGSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EffiToySystErr->GetBinContent(ipt),2)
	    + TMath::Power( hData_IDIsoSigShapeSystErr->GetBinContent(ipt),2)
	    + TMath::Power( hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MetResolSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EnMomScaleSystErr->GetBinContent(ipt),2) 
	    + TMath::Power(hData_EnMomSmearSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdBckgrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EwkSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_FsrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_SvdUnfSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_UnfoldBiasSystErr->GetBinContent(ipt),2)   ));
            
            ///+ TMath::Power( hData_LumiSystErr->GetBinContent(ipt),2) )); ///do not use for Normalized Xsec
    }
     
    
    
    for( int ipt(1);ipt<=n12Bins-1;ipt++)
    {
    
      ///used for Normalized Xsec 12
      hData_TotalSystErr12->SetBinContent(ipt, sqrt( 
	    TMath::Power(hData_TrackSigShapeSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_TrackBkgrShapeSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MuonPOGSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EffiToySystErr12->GetBinContent(ipt),2)
	    + TMath::Power( hData_IDIsoSigShapeSystErr12->GetBinContent(ipt),2)
	    + TMath::Power( hData_IDIsoBkgrShapeSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MetResolSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EnMomScaleSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power(hData_EnMomSmearSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdBckgrSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdShapeSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EwkSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_FsrSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_SvdUnfSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_UnfoldBiasSystErr12->GetBinContent(ipt),2)   ));

  }
}
else if(BaseName== "WpToEleNu" || BaseName== "WmToEleNu")
  {
    TH1D* hData_BinningSystErr		 = (TH1D*)f_Data_PowhegErr->Get("h_bin")->Clone("hData_BinningSystErr");
    TH1D *hData_BinningSystErr12   = new TH1D("hData_BinningSystErr12","hData_BinningSystErr12",n12Bins-1,Wpt12Bins);hData_BinningSystErr12->Sumw2();
    for( int ipt(1);ipt<n12Bins;ipt++)
    { 
      if(ipt<4)
      {
	hData_BinningSystErr12->SetBinContent(ipt, hData_BinningSystErr->GetBinContent(ipt) );
      }else if(ipt==4)
      {
	//hData_BinningSystErr12->SetBinContent(ipt, TMath::Max(hData_BinningSystErr->GetBinContent(ipt),hData_BinningSystErr->GetBinContent(ipt+1) ) );
	hData_BinningSystErr12->SetBinContent(ipt, (hData_BinningSystErr->GetBinContent(ipt)*BinWidth[4]*RD_Yield[4] + hData_BinningSystErr->GetBinContent(ipt+1)*BinWidth[5]*RD_Yield[5]) 
	    				/ (BinWidth[4]*RD_Yield[4] + BinWidth[5]*RD_Yield[5]) );
      }else if(ipt>=5)
      {
	hData_BinningSystErr12->SetBinContent(ipt, hData_BinningSystErr->GetBinContent(ipt+1) );
      }
    } 
  
    //// Calculate total syst for electron
    for( int ipt(1);ipt<=nBins-1;ipt++)
    {

      ///used for Normalized Xsec
      hData_TotalSystErr->SetBinContent(ipt, sqrt( 
	    TMath::Power(hData_BinningSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EffiToySystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_IDIsoSigShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MetResolSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EnMomScaleSystErr->GetBinContent(ipt),2) 
	    + TMath::Power(hData_EnMomSmearSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdBckgrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EwkSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_FsrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_SvdUnfSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_UnfoldBiasSystErr->GetBinContent(ipt),2)  ));

             ///+ TMath::Power( hData_LumiSystErr->GetBinContent(ipt),2) )); ///do not use for Normalized Xsec         
    } 
    
    
    for( int ipt(1);ipt<=n12Bins-1;ipt++)
    {
      ///used for Normalized Xsec 12
      hData_TotalSystErr12->SetBinContent(ipt, sqrt( 
	    TMath::Power(hData_BinningSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EffiToySystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_IDIsoSigShapeSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_IDIsoBkgrShapeSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MetResolSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EnMomScaleSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power(hData_EnMomSmearSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdBckgrSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdShapeSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EwkSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_FsrSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_SvdUnfSystErr12->GetBinContent(ipt),2) 
	    + TMath::Power( hData_UnfoldBiasSystErr12->GetBinContent(ipt),2)  ));
    }
  }

  TH1D *hData_TotalUncer   = new TH1D("hData_TotalUncer","hData_TotalUncer",nBins-1,WptBins);hData_TotalUncer->Sumw2();
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hData_TotalUncer->SetBinContent(ipt, sqrt(TMath::Power(hData_StatErr->GetBinContent(ipt),2)+ TMath::Power(hData_TotalSystErr->GetBinContent(ipt),2)));
    cout<<"hData_TotalUncer: "<<ipt<<"\t"<<hData_TotalUncer->GetBinContent(ipt)<<endl;
  } 
  TH1D *hData_TotalUncer12   = new TH1D("hData_TotalUncer12","hData_TotalUncer12",n12Bins-1,Wpt12Bins);hData_TotalUncer12->Sumw2();
  for( int ipt(1);ipt<n12Bins;ipt++)
  {
    hData_TotalUncer12->SetBinContent(ipt, sqrt(TMath::Power(hData_StatErr12->GetBinContent(ipt),2)+ TMath::Power(hData_TotalSystErr12->GetBinContent(ipt),2)));
    cout<<"hData_TotalUncer 12 bin: "<<ipt<<"\t"<<hData_TotalUncer12->GetBinContent(ipt)<<endl;
  } 
 

  //for( int ipt(1);ipt<nBins;ipt++)
  for( int ipt(1);ipt<n12Bins;ipt++)
  {
    //cout<<"hData_TrackSigShapeSystErr: "<<ipt<<"\t"<<hData_TrackSigShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_TrackBkgrShapeSystErr: "<<ipt<<"\t"<<hData_TrackBkgrShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_IDIsoSigShapeSystErr: "<<ipt<<"\t"<<hData_IDIsoSigShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_IDIsoBkgrShapeSystErr: "<<ipt<<"\t"<<hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_MuonPOGSystErr: "<<ipt<<"\t"<<hData_MuonPOGSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_EffiToySystErr: "<<ipt<<"\t"<<hData_EffiToySystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_MetResolSystErr12: "<<ipt<<"\t"<<hData_MetResolSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_EnMomScaleSystErr12: "<<ipt<<"\t"<<hData_EnMomScaleSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_EnMomSmearSystErr: "<<ipt<<"\t"<<hData_EnMomSmearSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_QcdBckgrSystErr12: "<<ipt<<"\t"<<hData_QcdBckgrSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_QcdShapeSystErr12: "<<ipt<<"\t"<<hData_QcdShapeSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_EwkSystErr12: "<<ipt<<"\t"<<hData_EwkSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_FsrSystErr12: "<<ipt<<"\t"<<hData_FsrSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_SvdUnfSystErr12: "<<ipt<<"\t"<<hData_SvdUnfSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_UnfoldBiasSystErr12: "<<ipt<<"\t"<<hData_UnfoldBiasSystErr12->GetBinContent(ipt)<<endl;
    //cout<<"hData_StatErr12: "<<ipt<<"\t"<<hData_StatErr12->GetBinContent(ipt)<<endl;
    cout<<"hData_TotalSystErr12: "<<ipt<<"\t"<<hData_TotalSystErr12->GetBinContent(ipt)<<endl;

    //Total Lepton Effi syst
//      cout<< " & " << 
//      sqrt(  TMath::Power(hData_TrackSigShapeSystErr12->GetBinContent(ipt),2) 
//	  +  TMath::Power(hData_TrackBkgrShapeSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_IDIsoSigShapeSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_IDIsoBkgrShapeSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_MuonPOGSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_EffiToySystErr12->GetBinContent(ipt),2)
//	  ) ;
    
    //Scale Smear
//    cout<< " & " << 
//      sqrt(  TMath::Power(hData_EnMomScaleSystErr12->GetBinContent(ipt),2) 
//    	  +  TMath::Power(hData_EnMomSmearSystErr12->GetBinContent(ipt),2)
//    	  ) ;

  } 
  hData_TrackSigShapeSystErr12->Write();
  hData_TrackBkgrShapeSystErr12->Write();
  hData_IDIsoSigShapeSystErr12->Write();
  hData_IDIsoBkgrShapeSystErr12->Write();
  hData_MuonPOGSystErr12->Write();
  hData_EffiToySystErr12->Write();
  hData_StatErr12->Write();
  hData_MetResolSystErr12->Write(); 
  hData_EnMomScaleSystErr12->Write();
  hData_EnMomSmearSystErr12->Write();
  hData_QcdBckgrSystErr12->Write();
  hData_QcdShapeSystErr12->Write();
  hData_EwkSystErr12->Write();
  hData_FsrSystErr12->Write();
  hData_SvdUnfSystErr12->Write();
  hData_UnfoldBiasSystErr12->Write();
  
  cout << "======= Put in note=======" << endl;
  cout << fixed << setprecision(3);
  for(int ipt(1);ipt<n12Bins;ipt++)
  {
    //cout<< " & " << hData_TrackSigShapeSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_TrackBkgrShapeSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_IDIsoSigShapeSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_IDIsoBkgrShapeSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_MuonPOGSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_EffiToySystErr12->GetBinContent(ipt) ;
   
//    ///Total Lepton Effi syst
//    cout<< " & " << 
//      sqrt(  TMath::Power(hData_TrackSigShapeSystErr12->GetBinContent(ipt),2) 
//	  +  TMath::Power(hData_TrackBkgrShapeSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_IDIsoSigShapeSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_IDIsoBkgrShapeSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_MuonPOGSystErr12->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_EffiToySystErr12->GetBinContent(ipt),2)
//	  ) ;

    //cout<< " & " << hData_EnMomScaleSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_EnMomSmearSystErr12->GetBinContent(ipt) ;
    ///Scale Smear
    //cout<< " & " << 
    //  sqrt(  TMath::Power(hData_EnMomScaleSystErr12->GetBinContent(ipt),2) 
    //	  +  TMath::Power(hData_EnMomSmearSystErr12->GetBinContent(ipt),2)
    //	  ) ;

    //cout<< " & " << hData_MetResolSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_QcdBckgrSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_QcdShapeSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_EwkSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_SvdUnfSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_FsrSystErr12->GetBinContent(ipt) ;
    //cout<< " & " << hData_UnfoldBiasSystErr12->GetBinContent(ipt) ;
    
    //Total syst
    //cout<< " & " << hData_TotalSystErr12->GetBinContent(ipt) ;
    // Stat error
    //cout<< " & " << hData_StatErr12->GetBinContent(ipt) ;
    //Total Uncertainty
    cout<< " & " << hData_TotalUncer12->GetBinContent(ipt) ;



  }
    cout << " \\"<<"\\" << endl;
    cout << endl;



  TH1D* hPowheg_PDFErr		 = (TH1D*)f_Data_PowhegErr->Get("h_PowhegPDF")->Clone("hPowheg_PDFErr");
  TH1D *hPowheg_PDFErr12   = new TH1D("hPowheg_PDFErr12","hPowheg_PDFErr12",n12Bins-1,Wpt12Bins);hPowheg_PDFErr12->Sumw2();
  for( int ipt(1);ipt<n12Bins;ipt++)
  { 
    if(ipt<4)
    {
      hPowheg_PDFErr12->SetBinContent(ipt, hPowheg_PDFErr->GetBinContent(ipt) );
    }else if(ipt==4)
    {
      //hPowheg_PDFErr12->SetBinContent(ipt, TMath::Max(hPowheg_PDFErr->GetBinContent(ipt),hPowheg_PDFErr->GetBinContent(ipt+1) ) );
      hPowheg_PDFErr12->SetBinContent(ipt, (hPowheg_PDFErr->GetBinContent(ipt)*BinWidth[4]*Powheg_Yield[4] + hPowheg_PDFErr->GetBinContent(ipt+1)*BinWidth[5]*Powheg_Yield[5] )
	  				/(BinWidth[4]*Powheg_Yield[4] + BinWidth[5]*Powheg_Yield[5])	);
    }else if(ipt>=5)
    {
      hPowheg_PDFErr12->SetBinContent(ipt, hPowheg_PDFErr->GetBinContent(ipt+1) );
    }
  cout<<ipt<<"\t hPowheg_PDFErr12\t"<<hPowheg_PDFErr12->GetBinContent(ipt)<<endl;
  } 
  
  //TH1D* hData_TotalUncer = (TH1D*)f_Data_PowhegErr->Get("h_TotalUncer")->Clone("hData_TotalUncer");
  ///=============Reading related stat, syst errors Finished here==============
 




  /////==========Real Data Starts here=========================================
  //
  // Born Level RD Yield and Errors, and Converting to X sec yield and errors, making TGraph .
  
 

  TH1D* hData_Yield_BornLinScale;
  hData_Yield_BornLinScale   = (TH1D*)f_Data->Get("BornEffCorr")->Clone("hData_Yield_BornLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hData_Yield_BornLinScale->SetBinError(ipt, 0.01*hData_TotalUncer->GetBinContent(ipt)*hData_Yield_BornLinScale->GetBinContent(ipt) );
    //cout<<ipt<<"\t hData_Yield_BornLinScale\t"<< hData_Yield_BornLinScale->GetBinContent(ipt)<<endl;
  } 
  TH1D *hData_Yield_BornLinScale12   = new TH1D("hData_Yield_BornLinScale12","hData_Yield_BornLinScale12",n12Bins-1,Wpt12Bins);hData_Yield_BornLinScale12->Sumw2();
  for( int ipt(1);ipt<n12Bins;ipt++)
  { 
    if(ipt<4)
    {
      hData_Yield_BornLinScale12->SetBinContent(ipt, hData_Yield_BornLinScale->GetBinContent(ipt) );
    }else if(ipt==4)
    {
      hData_Yield_BornLinScale12->SetBinContent(ipt, hData_Yield_BornLinScale->GetBinContent(ipt) + hData_Yield_BornLinScale->GetBinContent(ipt+1) ) ;
    }else if(ipt>=5)
    {
      hData_Yield_BornLinScale12->SetBinContent(ipt, hData_Yield_BornLinScale->GetBinContent(ipt+1) );
    }
  } 
  for( int ipt(1);ipt<n12Bins;ipt++)
  {
    hData_Yield_BornLinScale12->SetBinError(ipt, 0.01*hData_TotalUncer12->GetBinContent(ipt)*hData_Yield_BornLinScale12->GetBinContent(ipt) );
  } 
  hData_Yield_BornLinScale->Scale(1./18.429);
  hData_Yield_BornLinScale12->Scale(1./18.429);

  TH1D* hData_Yield_ReconLinScale;
  hData_Yield_ReconLinScale = (TH1D*)f_Data->Get("data_Rec")->Clone("hData_Yield_ReconLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    double tmp = sqrt(hData_Yield_ReconLinScale->GetBinContent(ipt));
    hData_Yield_ReconLinScale->SetBinError(ipt,tmp);
  }
  TH1D *hData_Yield_ReconLinScale12   = new TH1D("hData_Yield_ReconLinScale12","hData_Yield_ReconLinScale12",n12Bins-1,Wpt12Bins);hData_Yield_ReconLinScale12->Sumw2();
  for( int ipt(1);ipt<n12Bins;ipt++)
  { 
    if(ipt<4)
    {
      hData_Yield_ReconLinScale12->SetBinContent(ipt, hData_Yield_ReconLinScale->GetBinContent(ipt) );
    }else if(ipt==4)
    {
      hData_Yield_ReconLinScale12->SetBinContent(ipt, hData_Yield_ReconLinScale->GetBinContent(ipt) + hData_Yield_ReconLinScale->GetBinContent(ipt+1) ) ;
    }else if(ipt>=5)
    {
      hData_Yield_ReconLinScale12->SetBinContent(ipt, hData_Yield_ReconLinScale->GetBinContent(ipt+1) );
    }
  } 
  for( int ipt(1);ipt<n12Bins;ipt++)
  {
    double tmp = sqrt(hData_Yield_ReconLinScale12->GetBinContent(ipt));
    hData_Yield_ReconLinScale12->SetBinError(ipt,tmp);
  }

  //Errors
  double Data_Xsec_BornStatErr[14];
  double Data_Xsec_BornTotalErr[14];
  double Data_Xsec_BornStatErr12[13];
  double Data_Xsec_BornTotalErr12[13];
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    Data_Xsec_BornStatErr[ipt] = hData_Yield_BornLinScale->GetBinContent(ipt) * ( hData_Yield_ReconLinScale->GetBinError(ipt)/hData_Yield_ReconLinScale->GetBinContent(ipt) );
    Data_Xsec_BornTotalErr[ipt] = hData_Yield_BornLinScale->GetBinError(ipt);
  } 
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    Data_Xsec_BornStatErr12[ipt] = hData_Yield_BornLinScale12->GetBinContent(ipt) * ( hData_Yield_ReconLinScale12->GetBinError(ipt)/hData_Yield_ReconLinScale12->GetBinContent(ipt) );
    Data_Xsec_BornTotalErr12[ipt] = hData_Yield_BornLinScale12->GetBinError(ipt);
  }

  // X-sec and error Log scale
  TH1D *hData_Xsec_BornLogScale     = new TH1D("hData_Xsec_BornLogScale","hData_Xsec_BornLogScale",nBins-1,WptLogBins);hData_Xsec_BornLogScale->Sumw2();
  TH1D *hData_Xsec_BornLogScale12     = new TH1D("hData_Xsec_BornLogScale12","hData_Xsec_BornLogScale12",n12Bins-1,WptLog12Bins);hData_Xsec_BornLogScale12->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hData_Xsec_BornLogScale->SetBinContent(ipt,hData_Yield_BornLinScale->GetBinContent(ipt)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt    ));
    hData_Xsec_BornLogScale->SetBinError(ipt,Data_Xsec_BornTotalErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  } 
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hData_Xsec_BornLogScale12->SetBinContent(ipt,hData_Yield_BornLinScale12->GetBinContent(ipt)/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt    ));
    hData_Xsec_BornLogScale12->SetBinError(ipt,Data_Xsec_BornTotalErr12[ipt]/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt));
  }
 
  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataErrBand = new TH1D("hRatioDataErrBand","hRatioDataErrBand",nBins-1,WptLogBins);hRatioDataErrBand->Sumw2();
  TH1D *hRatioDataStatErr = new TH1D("hRatioDataStatErr","hRatioDataStatErr",nBins-1,WptLogBins);hRatioDataStatErr->Sumw2();
  TH1D *hRatioDataErrBand12 = new TH1D("hRatioDataErrBand12","hRatioDataErrBand12",n12Bins-1,WptLog12Bins);hRatioDataErrBand12->Sumw2();
  TH1D *hRatioDataStatErr12 = new TH1D("hRatioDataStatErr12","hRatioDataStatErr12",n12Bins-1,WptLog12Bins);hRatioDataStatErr12->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioDataStatErr->SetBinContent(ipt,1.);
    hRatioDataStatErr->SetBinError(ipt,Data_Xsec_BornStatErr[ipt]/hData_Yield_BornLinScale->GetBinContent(ipt));
    
    hRatioDataErrBand->SetBinContent(ipt,1.);
    hRatioDataErrBand->SetBinError(ipt,hData_Yield_BornLinScale->GetBinError(ipt)/hData_Yield_BornLinScale->GetBinContent(ipt));
  }
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hRatioDataStatErr12->SetBinContent(ipt,1.);
    hRatioDataStatErr12->SetBinError(ipt,Data_Xsec_BornStatErr12[ipt]/hData_Yield_BornLinScale12->GetBinContent(ipt));
    
    hRatioDataErrBand12->SetBinContent(ipt,1.);
    hRatioDataErrBand12->SetBinError(ipt,hData_Yield_BornLinScale12->GetBinError(ipt)/hData_Yield_BornLinScale12->GetBinContent(ipt));
  }
  
  ////Normalize Data X-sec and its errors
  TH1D* hData_Xsec_BornLogScaleNorm = (TH1D*)hData_Xsec_BornLogScale->Clone("hData_Xsec_BornLogScaleNorm");
  //hData_Xsec_BornLogScaleNorm->Scale(1./hData_Xsec_BornLogScale->Integral());
  hData_Xsec_BornLogScaleNorm->Scale(1./hData_Yield_BornLinScale->Integral());
  TH1D* hData_Xsec_BornLogScaleNorm12 = (TH1D*)hData_Xsec_BornLogScale12->Clone("hData_Xsec_BornLogScaleNorm12");
  //hData_Xsec_BornLogScaleNorm12->Scale(1./hData_Xsec_BornLogScale12->Integral());
  hData_Xsec_BornLogScaleNorm12->Scale(1./hData_Yield_BornLinScale->Integral());
  
  hData_Xsec_BornLogScaleNorm12->Write();
  cout << fixed << setprecision(8);
  for(int i(1); i<=12; i++)
{
  cout << "hData_Xsec_BornLogScaleNorm12 Error : " << hData_Xsec_BornLogScaleNorm12->GetBinContent(i) << "\t" <<hData_Xsec_BornLogScaleNorm12->GetBinError(i) << endl;
}
  
  cout << fixed << setprecision(3);
  
  TH1D *hRatioDataStatErrNorm = (TH1D*)hRatioDataStatErr->Clone("hRatioDataStatErrNorm");
  //hRatioDataStatErrNorm->Scale(1./hData_Xsec_BornLogScale->Integral());
  hRatioDataStatErrNorm->Scale(1./hData_Yield_BornLinScale->Integral());
  TH1D *hRatioDataStatErrNorm12 = (TH1D*)hRatioDataStatErr12->Clone("hRatioDataStatErrNorm12");
  //hRatioDataStatErrNorm12->Scale(1./hData_Xsec_BornLogScale12->Integral());
  hRatioDataStatErrNorm12->Scale(1./hData_Yield_BornLinScale12->Integral());
  
  TH1D *hRatioDataErrBandNorm = (TH1D*)hRatioDataErrBand->Clone("hRatioDataErrBandNorm");
  //hRatioDataErrBandNorm->Scale(1./hData_Xsec_BornLogScale->Integral());
  hRatioDataErrBandNorm->Scale(1./hData_Yield_BornLinScale->Integral());
  TH1D *hRatioDataErrBandNorm12 = (TH1D*)hRatioDataErrBand12->Clone("hRatioDataErrBandNorm12");
  //hRatioDataErrBandNorm12->Scale(1./hData_Xsec_BornLogScale12->Integral());
  hRatioDataErrBandNorm12->Scale(1./hData_Yield_BornLinScale12->Integral());
  
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioDataStatErrNorm->SetBinContent(ipt,1.);
    hRatioDataErrBandNorm->SetBinContent(ipt,1.);
    //cout<<"Data Xsec After Norm: "<<hData_Xsec_BornLogScaleNorm->GetBinContent(ipt)<<"\t"<<hData_Xsec_BornLogScaleNorm->GetBinContent(ipt)<<endl;
  }
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hRatioDataStatErrNorm12->SetBinContent(ipt,1.);
    hRatioDataErrBandNorm12->SetBinContent(ipt,1.);
    //cout<<"Data Xsec After Norm: "<<hData_Xsec_BornLogScaleNorm->GetBinContent(ipt)<<"\t"<<hData_Xsec_BornLogScaleNorm->GetBinContent(ipt)<<endl;
  }
  

  hData_Xsec_BornLogScale->SetMarkerStyle(kFullCircle); 
  hData_Xsec_BornLogScale->SetMarkerColor(kBlack); 
  hData_Xsec_BornLogScale->SetMarkerSize(1);
  hData_Xsec_BornLogScale12->SetMarkerStyle(kFullCircle); 
  hData_Xsec_BornLogScale12->SetMarkerColor(kBlack); 
  hData_Xsec_BornLogScale12->SetMarkerSize(1);
  
  hData_Xsec_BornLogScaleNorm->SetMarkerStyle(kFullCircle); 
  hData_Xsec_BornLogScaleNorm->SetMarkerColor(kBlack); 
  hData_Xsec_BornLogScaleNorm->SetMarkerSize(1);
  hData_Xsec_BornLogScaleNorm12->SetMarkerStyle(kFullCircle); 
  hData_Xsec_BornLogScaleNorm12->SetMarkerColor(kBlack); 
  hData_Xsec_BornLogScaleNorm12->SetMarkerSize(1);
  
  hRatioDataStatErr->SetMarkerStyle(kFullCircle); 
  hRatioDataStatErr->SetMarkerColor(kBlack); 
  hRatioDataStatErr->SetMarkerSize(0.6);
  hRatioDataStatErr12->SetMarkerStyle(kFullCircle); 
  hRatioDataStatErr12->SetMarkerColor(kBlack); 
  hRatioDataStatErr12->SetMarkerSize(0.6);
  
  hRatioDataStatErrNorm->SetMarkerStyle(kFullCircle); 
  hRatioDataStatErrNorm->SetMarkerColor(kBlack); 
  hRatioDataStatErrNorm->SetMarkerSize(0.6);
  hRatioDataStatErrNorm12->SetMarkerStyle(kFullCircle); 
  hRatioDataStatErrNorm12->SetMarkerColor(kBlack); 
  hRatioDataStatErrNorm12->SetMarkerSize(0.6);
  
  /// TGraph
  TGraphErrors *Data_Xsec_Born = new TGraphErrors(hData_Xsec_BornLogScale);
  TGraphErrors *RatioDataErrBand = new TGraphErrors(hRatioDataErrBand);
  TGraphErrors *Data_Xsec_Born12 = new TGraphErrors(hData_Xsec_BornLogScale12);
  TGraphErrors *RatioDataErrBand12 = new TGraphErrors(hRatioDataErrBand12);
  //Normalized
  TGraphErrors *Data_Xsec_BornNorm = new TGraphErrors(hData_Xsec_BornLogScaleNorm);
  TGraphErrors *RatioDataErrBandNorm = new TGraphErrors(hRatioDataErrBandNorm);
  TGraphErrors *Data_Xsec_BornNorm12 = new TGraphErrors(hData_Xsec_BornLogScaleNorm12);
  TGraphErrors *RatioDataErrBandNorm12 = new TGraphErrors(hRatioDataErrBandNorm12);

  /////==========Real Data related stuff Finished here=========================================



  /////==========Powheg Starts here==============================================
  //
  // Born Level Powheg Yield and Errors, and Converting to X sec yield and errors, making TGraph .
  
  
  TH1D* hPowheg_Xsec_BornLinScale;
  hPowheg_Xsec_BornLinScale   = (TH1D*)f_Data->Get("SVD_Born.Gen")->Clone("hPowheg_Xsec_BornLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hPowheg_Xsec_BornLinScale->SetBinError(ipt, 0.01*hPowheg_PDFErr->GetBinContent(ipt)*hPowheg_Xsec_BornLinScale->GetBinContent(ipt) );
    
    cout<<"Powheg Yield:\t"<<hPowheg_Xsec_BornLinScale->GetBinContent(ipt)<<endl;
  } 


  TH1D *hPowheg_Xsec_BornLinScale12   = new TH1D("hPowheg_Xsec_BornLinScale12","hPowheg_Xsec_BornLinScale12",n12Bins-1,Wpt12Bins);hPowheg_Xsec_BornLinScale12->Sumw2();
  for( int ipt(1);ipt<n12Bins;ipt++)
  { 
    if(ipt<4)
    {
      hPowheg_Xsec_BornLinScale12->SetBinContent(ipt, hPowheg_Xsec_BornLinScale->GetBinContent(ipt) );
    }else if(ipt==4)
    {
      hPowheg_Xsec_BornLinScale12->SetBinContent(ipt, hPowheg_Xsec_BornLinScale->GetBinContent(ipt) + hPowheg_Xsec_BornLinScale->GetBinContent(ipt+1)  );
    }else if(ipt>=5)
    {
      hPowheg_Xsec_BornLinScale12->SetBinContent(ipt, hPowheg_Xsec_BornLinScale->GetBinContent(ipt+1) );
    }
  } 
  for( int ipt(1);ipt<n12Bins;ipt++)
  {
    hPowheg_Xsec_BornLinScale12->SetBinError(ipt, 0.01*hPowheg_PDFErr12->GetBinContent(ipt)*hPowheg_Xsec_BornLinScale12->GetBinContent(ipt) );
  } 
  hPowheg_Xsec_BornLinScale->Scale(1./18.429);
  hPowheg_Xsec_BornLinScale12->Scale(1./18.429);
    
  TH1D* hPowheg_Yield_BornAfterFidCut;
  hPowheg_Yield_BornAfterFidCut = (TH1D*)f_Data->Get("SVD_Born.Gen")->Clone("hPowheg_Yield_BornAfterFidCut");
  if(BaseName== "WpToMuNu")  { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Muon_WpToMuNu_S8); }
  if(BaseName== "WmToMuNu")  { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Muon_WmToMuNu_S8); }
  if(BaseName== "WpToEleNu") { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Ele_WpToEleNu_S8); }
  if(BaseName== "WmToEleNu") { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Ele_WmToEleNu_S8); }
  
  TH1D *hPowheg_Yield_BornAfterFidCut12   = new TH1D("hPowheg_Yield_BornAfterFidCut12","hPowheg_Yield_BornAfterFidCut12",n12Bins-1,Wpt12Bins);hPowheg_Yield_BornAfterFidCut12->Sumw2();
  for( int ipt(1);ipt<n12Bins;ipt++)
  { 
    if(ipt<4)
    {
      hPowheg_Yield_BornAfterFidCut12->SetBinContent(ipt, hPowheg_Yield_BornAfterFidCut->GetBinContent(ipt) );
    }else if(ipt==4)
    {
      hPowheg_Yield_BornAfterFidCut12->SetBinContent(ipt, hPowheg_Yield_BornAfterFidCut->GetBinContent(ipt) + hPowheg_Yield_BornAfterFidCut->GetBinContent(ipt+1)  );
    }else if(ipt>=5)
    {
      hPowheg_Yield_BornAfterFidCut12->SetBinContent(ipt, hPowheg_Yield_BornAfterFidCut->GetBinContent(ipt+1) );
    }
  } 

  //Errors
  double Powheg_Xsec_BornStatErr[14];
  double Powheg_Xsec_BornTotalErr[14];
  double Powheg_Xsec_BornStatErr12[13];
  double Powheg_Xsec_BornTotalErr12[13];
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    Powheg_Xsec_BornStatErr[ipt] = hPowheg_Xsec_BornLinScale->GetBinContent(ipt)*sqrt(hPowheg_Yield_BornAfterFidCut->GetBinContent(ipt))/hPowheg_Yield_BornAfterFidCut->GetBinContent(ipt);
    Powheg_Xsec_BornTotalErr[ipt] = sqrt(Powheg_Xsec_BornStatErr[ipt]*Powheg_Xsec_BornStatErr[ipt] + hPowheg_Xsec_BornLinScale->GetBinError(ipt)*hPowheg_Xsec_BornLinScale->GetBinError(ipt));
  }
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    Powheg_Xsec_BornStatErr12[ipt] = hPowheg_Xsec_BornLinScale12->GetBinContent(ipt)*sqrt(hPowheg_Yield_BornAfterFidCut12->GetBinContent(ipt))/hPowheg_Yield_BornAfterFidCut12->GetBinContent(ipt);
    Powheg_Xsec_BornTotalErr12[ipt] = sqrt(Powheg_Xsec_BornStatErr12[ipt]*Powheg_Xsec_BornStatErr12[ipt] + hPowheg_Xsec_BornLinScale12->GetBinError(ipt)*hPowheg_Xsec_BornLinScale12->GetBinError(ipt));
  }

  // X-sec and error Log scale
  TH1D *hPowheg_Xsec_BornLogScale   = new TH1D("hPowheg_Xsec_BornLogScale","hPowheg_Xsec_BornLogScale",nBins-1,WptLogBins);hPowheg_Xsec_BornLogScale->Sumw2();
  TH1D *hPowheg_Xsec_BornLogScale12   = new TH1D("hPowheg_Xsec_BornLogScale12","hPowheg_Xsec_BornLogScale12",n12Bins-1,WptLog12Bins);hPowheg_Xsec_BornLogScale12->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hPowheg_Xsec_BornLogScale->SetBinContent(ipt,hPowheg_Xsec_BornLinScale->GetBinContent(ipt)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
    hPowheg_Xsec_BornLogScale->SetBinError(ipt,Powheg_Xsec_BornTotalErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  }
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hPowheg_Xsec_BornLogScale12->SetBinContent(ipt,hPowheg_Xsec_BornLinScale12->GetBinContent(ipt)/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt));
    hPowheg_Xsec_BornLogScale12->SetBinError(ipt,Powheg_Xsec_BornTotalErr12[ipt]/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt));
  }

  
  ///   Theory/Data ratio plot errors related to Powheg
  TH1D *hRatioPowhegStatErrBand = new TH1D("hRatioPowhegStatErrBand","hRatioPowhegStatErrBand",nBins-1,WptLogBins);hRatioPowhegStatErrBand->Sumw2();
  TH1D *hRatioPowhegPDFErrBand = new TH1D("hRatioPowhegPDFErrBand","hRatioPowhegPDFErrBand",nBins-1,WptLogBins);hRatioPowhegPDFErrBand->Sumw2();
  TH1D *hRatioPowhegStatErrBand12 = new TH1D("hRatioPowhegStatErrBand12","hRatioPowhegStatErrBand12",n12Bins-1,WptLog12Bins);hRatioPowhegStatErrBand12->Sumw2();
  TH1D *hRatioPowhegPDFErrBand12 = new TH1D("hRatioPowhegPDFErrBand12","hRatioPowhegPDFErrBand12",n12Bins-1,WptLog12Bins);hRatioPowhegPDFErrBand12->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioPowhegStatErrBand->SetBinContent(ipt,hPowheg_Xsec_BornLogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioPowhegStatErrBand->SetBinError(ipt,Powheg_Xsec_BornStatErr[ipt]/hData_Yield_BornLinScale->GetBinContent(ipt));
    
    hRatioPowhegPDFErrBand->SetBinContent(ipt,hPowheg_Xsec_BornLogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioPowhegPDFErrBand->SetBinError(ipt,(Powheg_Xsec_BornStatErr[ipt]/hData_Yield_BornLinScale->GetBinContent(ipt))+(hPowheg_Xsec_BornLinScale->GetBinError(ipt)/hData_Yield_BornLinScale->GetBinContent(ipt)));
  }
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hRatioPowhegStatErrBand12->SetBinContent(ipt,hPowheg_Xsec_BornLogScale12->GetBinContent(ipt)/hData_Xsec_BornLogScale12->GetBinContent(ipt));
    hRatioPowhegStatErrBand12->SetBinError(ipt,Powheg_Xsec_BornStatErr12[ipt]/hData_Yield_BornLinScale12->GetBinContent(ipt));
    
    hRatioPowhegPDFErrBand12->SetBinContent(ipt,hPowheg_Xsec_BornLogScale12->GetBinContent(ipt)/hData_Xsec_BornLogScale12->GetBinContent(ipt));
    hRatioPowhegPDFErrBand12->SetBinError(ipt,(Powheg_Xsec_BornStatErr12[ipt]/hData_Yield_BornLinScale12->GetBinContent(ipt))+(hPowheg_Xsec_BornLinScale12->GetBinError(ipt)/hData_Yield_BornLinScale12->GetBinContent(ipt)));
  }
  
  
  ////Normalize Powheg X-sec and errors
  TH1D *hPowheg_Xsec_BornLogScaleNorm= (TH1D*)hPowheg_Xsec_BornLogScale->Clone("hPowheg_Xsec_BornLogScaleNorm"); 
  //hPowheg_Xsec_BornLogScaleNorm->Scale(1./hPowheg_Xsec_BornLogScale->Integral());
  hPowheg_Xsec_BornLogScaleNorm->Scale(1./hPowheg_Xsec_BornLinScale->Integral());
  TH1D *hPowheg_Xsec_BornLogScaleNorm12= (TH1D*)hPowheg_Xsec_BornLogScale12->Clone("hPowheg_Xsec_BornLogScaleNorm12"); 
  //hPowheg_Xsec_BornLogScaleNorm12->Scale(1./hPowheg_Xsec_BornLogScale12->Integral());
  hPowheg_Xsec_BornLogScaleNorm12->Scale(1./hPowheg_Xsec_BornLinScale12->Integral());
  hPowheg_Xsec_BornLogScaleNorm12->Write();
  
  TH1D *hRatioPowhegStatErrBandNorm= (TH1D*)hRatioPowhegStatErrBand->Clone("hRatioPowhegStatErrBandNorm"); 
  //hRatioPowhegStatErrBandNorm->Scale(1./hPowheg_Xsec_BornLogScale->Integral());
  TH1D *hRatioPowhegStatErrBandNorm12= (TH1D*)hRatioPowhegStatErrBand12->Clone("hRatioPowhegStatErrBandNorm12"); 
  //hRatioPowhegStatErrBandNorm->Scale(1./hPowheg_Xsec_BornLogScale->Integral());
  
  TH1D *hRatioPowhegPDFErrBandNorm= (TH1D*)hRatioPowhegPDFErrBand->Clone("hRatioPowhegPDFErrBandNorm"); 
  //hRatioPowhegPDFErrBandNorm->Scale(1./hPowheg_Xsec_BornLogScale->Integral());
  TH1D *hRatioPowhegPDFErrBandNorm12= (TH1D*)hRatioPowhegPDFErrBand12->Clone("hRatioPowhegPDFErrBandNorm12"); 
  //hRatioPowhegPDFErrBandNorm12->Scale(1./hPowheg_Xsec_BornLogScale12->Integral());
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioPowhegStatErrBandNorm->SetBinContent(ipt,hPowheg_Xsec_BornLogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
    hRatioPowhegPDFErrBandNorm->SetBinContent(ipt,hPowheg_Xsec_BornLogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
  }
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hRatioPowhegStatErrBandNorm12->SetBinContent(ipt,hPowheg_Xsec_BornLogScaleNorm12->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt));
    hRatioPowhegPDFErrBandNorm12->SetBinContent(ipt,hPowheg_Xsec_BornLogScaleNorm12->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt));
  }
  
  
  /// TGraph
  TGraphErrors *Powheg_Xsec_Born = new TGraphErrors(hPowheg_Xsec_BornLogScale);
  TGraphErrors *RatioPowhegStatErrBand = new TGraphErrors(hRatioPowhegStatErrBand);
  TGraphErrors *RatioPowhegPDFErrBand = new TGraphErrors(hRatioPowhegPDFErrBand);
  TGraphErrors *Powheg_Xsec_Born12 = new TGraphErrors(hPowheg_Xsec_BornLogScale12);
  TGraphErrors *RatioPowhegStatErrBand12 = new TGraphErrors(hRatioPowhegStatErrBand12);
  TGraphErrors *RatioPowhegPDFErrBand12 = new TGraphErrors(hRatioPowhegPDFErrBand12);
  /// Normalized
  TGraphErrors *Powheg_Xsec_BornNorm = new TGraphErrors(hPowheg_Xsec_BornLogScaleNorm);
  TGraphErrors *RatioPowhegStatErrBandNorm = new TGraphErrors(hRatioPowhegStatErrBandNorm);
  TGraphErrors *RatioPowhegPDFErrBandNorm = new TGraphErrors(hRatioPowhegPDFErrBandNorm);
  TGraphErrors *Powheg_Xsec_BornNorm12 = new TGraphErrors(hPowheg_Xsec_BornLogScaleNorm12);
  TGraphErrors *RatioPowhegStatErrBandNorm12 = new TGraphErrors(hRatioPowhegStatErrBandNorm12);
  TGraphErrors *RatioPowhegPDFErrBandNorm12 = new TGraphErrors(hRatioPowhegPDFErrBandNorm12);

  /////==========Powheg related stuff Finished here=========================================
  
  
  
  /////==========FEWZ Starts here==============================================
  //
  // FEWZ X sec value and errors, making TGraph .
  TH1D* hFEWZ_Xsec_LinScale;
  hFEWZ_Xsec_LinScale = (TH1D*)f_Fewz->Get("hxsec")->Clone("hFEWZ_Xsec_LinScale");
  
  TH1D* hFEWZ_Xsec_LinScale12;
  hFEWZ_Xsec_LinScale12 = (TH1D*)f_Fewz_12->Get("hxsec")->Clone("hFEWZ_Xsec_LinScale12");
  
  //Errors
  TH1D* hFEWZ_Xsec_ScaleError;
  hFEWZ_Xsec_ScaleError = (TH1D*)f_Fewz->Get("ScaleErr")->Clone("hFEWZ_Xsec_ScaleError");
  TH1D* hFEWZ_Xsec_ScaleError12;
  hFEWZ_Xsec_ScaleError12 = (TH1D*)f_Fewz_12->Get("Norm_ScaleErr")->Clone("hFEWZ_Xsec_ScaleError12");
  TH1D* hFEWZ_Xsec_PDFError;
  hFEWZ_Xsec_PDFError = (TH1D*)f_Fewz->Get("PDFErr")->Clone("hFEWZ_Xsec_PDFError");
  TH1D* hFEWZ_Xsec_PDFError12;
  hFEWZ_Xsec_PDFError12 = (TH1D*)f_Fewz_12->Get("Norm_PDFErr")->Clone("hFEWZ_Xsec_PDFError12");
  
  double FEWZ_Xsec_ScaleErr[14];
  double FEWZ_Xsec_PDFErr[14];
  double FEWZ_Xsec_TotalErr[14];
  double FEWZ_Xsec_ScaleErr12[13];
  double FEWZ_Xsec_PDFErr12[13];
  double FEWZ_Xsec_TotalErr12[13];
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    FEWZ_Xsec_ScaleErr[ipt] = 0.01*hFEWZ_Xsec_ScaleError->GetBinError(ipt)*hFEWZ_Xsec_LinScale->GetBinContent(ipt);
    FEWZ_Xsec_PDFErr[ipt]   = 0.01*hFEWZ_Xsec_PDFError->GetBinError(ipt)*hFEWZ_Xsec_LinScale->GetBinContent(ipt);
    
    FEWZ_Xsec_TotalErr[ipt] = sqrt( FEWZ_Xsec_ScaleErr[ipt]*FEWZ_Xsec_ScaleErr[ipt] + FEWZ_Xsec_PDFErr[ipt]*FEWZ_Xsec_PDFErr[ipt] );
  } 
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    FEWZ_Xsec_ScaleErr12[ipt] = 0.01*hFEWZ_Xsec_ScaleError12->GetBinError(ipt)*hFEWZ_Xsec_LinScale12->GetBinContent(ipt);
    FEWZ_Xsec_PDFErr12[ipt]   = 0.01*hFEWZ_Xsec_PDFError12->GetBinError(ipt)*hFEWZ_Xsec_LinScale12->GetBinContent(ipt);
    
    FEWZ_Xsec_TotalErr12[ipt] = sqrt( FEWZ_Xsec_ScaleErr12[ipt]*FEWZ_Xsec_ScaleErr12[ipt] + FEWZ_Xsec_PDFErr12[ipt]*FEWZ_Xsec_PDFErr12[ipt] );
  } 
  
  // X-sec and error Log scale
  TH1D *hFEWZ_Xsec_LogScale   = new TH1D("hFEWZ_Xsec_LogScale","hFEWZ_Xsec_LogScale",nBins-1,WptLogBins);hFEWZ_Xsec_LogScale->Sumw2();
  TH1D *hFEWZ_Xsec_LogScale12   = new TH1D("hFEWZ_Xsec_LogScale12","hFEWZ_Xsec_LogScale12",n12Bins-1,WptLog12Bins);hFEWZ_Xsec_LogScale12->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hFEWZ_Xsec_LogScale->SetBinContent(ipt,hFEWZ_Xsec_LinScale->GetBinContent(ipt)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
    hFEWZ_Xsec_LogScale->SetBinError(ipt,FEWZ_Xsec_TotalErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  } 
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hFEWZ_Xsec_LogScale12->SetBinContent(ipt,hFEWZ_Xsec_LinScale12->GetBinContent(ipt)/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt));
    hFEWZ_Xsec_LogScale12->SetBinError(ipt,FEWZ_Xsec_TotalErr12[ipt]/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt));
  } 
  
  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioFEWZStatErrBand = new TH1D("hRatioFEWZStatErrBand","hRatioFEWZStatErrBand",nBins-1,WptLogBins);hRatioFEWZStatErrBand->Sumw2();
  TH1D *hRatioFEWZScaleErr = new TH1D("hRatioFEWZScaleErr","hRatioFEWZScaleErr",nBins-1,WptLogBins);hRatioFEWZScaleErr->Sumw2();
  TH1D *hRatioFEWZQCDScaleErrBand = new TH1D("hRatioFEWZQCDScaleErrBand","hRatioFEWZQCDScaleErrBand",nBins-1,WptLogBins);hRatioFEWZQCDScaleErrBand->Sumw2();
  TH1D *hRatioFEWZScalePDFErrBand = new TH1D("hRatioFEWZScalePDFErrBand","hRatioFEWZScalePDFErrBand",nBins-1,WptLogBins);hRatioFEWZScalePDFErrBand->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioFEWZStatErrBand->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZStatErrBand->SetBinError(ipt,0.001);
    
    hRatioFEWZScaleErr->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZScaleErr->SetBinError(ipt,FEWZ_Xsec_ScaleErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  
    hRatioFEWZQCDScaleErrBand->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZQCDScaleErrBand->SetBinError(ipt,0.001 + hRatioFEWZScaleErr->GetBinError(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
  
    hRatioFEWZScalePDFErrBand->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZScalePDFErrBand->SetBinError(ipt,0.001+(hFEWZ_Xsec_LogScale->GetBinError(ipt)+hRatioFEWZScaleErr->GetBinError(ipt))/hData_Xsec_BornLogScale->GetBinContent(ipt));
  }
  TH1D *hRatioFEWZStatErrBand12 = new TH1D("hRatioFEWZStatErrBand12","hRatioFEWZStatErrBand12",n12Bins-1,WptLog12Bins);hRatioFEWZStatErrBand12->Sumw2();
  TH1D *hRatioFEWZScaleErr12 = new TH1D("hRatioFEWZScaleErr12","hRatioFEWZScaleErr12",n12Bins-1,WptLog12Bins);hRatioFEWZScaleErr12->Sumw2();
  TH1D *hRatioFEWZQCDScaleErrBand12 = new TH1D("hRatioFEWZQCDScaleErrBand12","hRatioFEWZQCDScaleErrBand12",n12Bins-1,WptLog12Bins);hRatioFEWZQCDScaleErrBand12->Sumw2();
  TH1D *hRatioFEWZScalePDFErrBand12 = new TH1D("hRatioFEWZScalePDFErrBand12","hRatioFEWZScalePDFErrBand12",n12Bins-1,WptLog12Bins);hRatioFEWZScalePDFErrBand12->Sumw2();
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hRatioFEWZStatErrBand12->SetBinContent(ipt,hFEWZ_Xsec_LogScale12->GetBinContent(ipt)/hData_Xsec_BornLogScale12->GetBinContent(ipt));
    hRatioFEWZStatErrBand12->SetBinError(ipt,0.001);
    
    hRatioFEWZScaleErr12->SetBinContent(ipt,hFEWZ_Xsec_LogScale12->GetBinContent(ipt)/hData_Xsec_BornLogScale12->GetBinContent(ipt));
    hRatioFEWZScaleErr12->SetBinError(ipt,FEWZ_Xsec_ScaleErr12[ipt]/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt));
  
    hRatioFEWZQCDScaleErrBand12->SetBinContent(ipt,hFEWZ_Xsec_LogScale12->GetBinContent(ipt)/hData_Xsec_BornLogScale12->GetBinContent(ipt));
    hRatioFEWZQCDScaleErrBand12->SetBinError(ipt,0.001 + hRatioFEWZScaleErr12->GetBinError(ipt)/hData_Xsec_BornLogScale12->GetBinContent(ipt));
  
    hRatioFEWZScalePDFErrBand12->SetBinContent(ipt,hFEWZ_Xsec_LogScale12->GetBinContent(ipt)/hData_Xsec_BornLogScale12->GetBinContent(ipt));
    hRatioFEWZScalePDFErrBand12->SetBinError(ipt,0.001+(hFEWZ_Xsec_LogScale12->GetBinError(ipt)+hRatioFEWZScaleErr12->GetBinError(ipt))/hData_Xsec_BornLogScale12->GetBinContent(ipt));
  }
 

  ////Normalize FEWZ X-sec and errors
  TH1D *hFEWZ_Xsec_LogScaleNorm = (TH1D*)hFEWZ_Xsec_LogScale->Clone("hFEWZ_Xsec_LogScaleNorm");
  //hFEWZ_Xsec_LogScaleNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());
  hFEWZ_Xsec_LogScaleNorm->Scale(1./hFEWZ_Xsec_LinScale->Integral());
  TH1D *hFEWZ_Xsec_LogScaleNorm12 = (TH1D*)hFEWZ_Xsec_LogScale12->Clone("hFEWZ_Xsec_LogScaleNorm12");
  //hFEWZ_Xsec_LogScaleNorm12->Scale(1./hFEWZ_Xsec_LogScale12->Integral());
  hFEWZ_Xsec_LogScaleNorm12->Scale(1./hFEWZ_Xsec_LinScale12->Integral());
  hFEWZ_Xsec_LogScaleNorm12->Write();
  

  TH1D *hRatioFEWZStatErrBandNorm = (TH1D*)hRatioFEWZStatErrBand->Clone("hRatioFEWZStatErrBandNorm");
  //hRatioFEWZStatErrBandNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());
  TH1D *hRatioFEWZStatErrBandNorm12 = (TH1D*)hRatioFEWZStatErrBand12->Clone("hRatioFEWZStatErrBandNorm12");
  //hRatioFEWZStatErrBandNorm12->Scale(1./hFEWZ_Xsec_LogScale12->Integral());
  
  TH1D *hRatioFEWZQCDScaleErrBandNorm = (TH1D*)hRatioFEWZQCDScaleErrBand->Clone("hRatioFEWZQCDScaleErrBandNorm");
  //hRatioFEWZQCDScaleErrBandNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());
  TH1D *hRatioFEWZQCDScaleErrBandNorm12 = (TH1D*)hRatioFEWZQCDScaleErrBand12->Clone("hRatioFEWZQCDScaleErrBandNorm12");
  //hRatioFEWZQCDScaleErrBandNorm12->Scale(1./hFEWZ_Xsec_LogScale12->Integral());
  
  TH1D *hRatioFEWZScalePDFErrBandNorm = (TH1D*)hRatioFEWZScalePDFErrBand->Clone("hRatioFEWZScalePDFErrBandNorm");
  //hRatioFEWZScalePDFErrBandNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());
  TH1D *hRatioFEWZScalePDFErrBandNorm12 = (TH1D*)hRatioFEWZScalePDFErrBand12->Clone("hRatioFEWZScalePDFErrBandNorm12");
  //hRatioFEWZScalePDFErrBandNorm12->Scale(1./hFEWZ_Xsec_LogScale12->Integral());
  
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioFEWZStatErrBandNorm->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
    hRatioFEWZQCDScaleErrBandNorm->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
    hRatioFEWZScalePDFErrBandNorm->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
  } 
  for( int ipt(1);ipt<=n12Bins-1;ipt++)
  {
    hRatioFEWZStatErrBandNorm12->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm12->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt));
    hRatioFEWZQCDScaleErrBandNorm12->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm12->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt));
    hRatioFEWZScalePDFErrBandNorm12->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm12->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt));
  } 

  /// TGraph
  TGraphErrors *FEWZ_Xsec = new TGraphErrors(hFEWZ_Xsec_LogScale);
  TGraphErrors* RatioFEWZStatErrBand = new TGraphErrors(hRatioFEWZStatErrBand);
  TGraphErrors* RatioFEWZQCDScaleErrBand = new TGraphErrors(hRatioFEWZQCDScaleErrBand);
  TGraphErrors* RatioFEWZScalePDFErrBand = new TGraphErrors(hRatioFEWZScalePDFErrBand);
  TGraphErrors *FEWZ_Xsec12 = new TGraphErrors(hFEWZ_Xsec_LogScale12);
  TGraphErrors* RatioFEWZStatErrBand12 = new TGraphErrors(hRatioFEWZStatErrBand12);
  TGraphErrors* RatioFEWZQCDScaleErrBand12 = new TGraphErrors(hRatioFEWZQCDScaleErrBand12);
  TGraphErrors* RatioFEWZScalePDFErrBand12 = new TGraphErrors(hRatioFEWZScalePDFErrBand12);
  /// Normalized
  TGraphErrors *FEWZ_XsecNorm = new TGraphErrors(hFEWZ_Xsec_LogScaleNorm);
  TGraphErrors* RatioFEWZStatErrBandNorm = new TGraphErrors(hRatioFEWZStatErrBandNorm);
  TGraphErrors* RatioFEWZQCDScaleErrBandNorm = new TGraphErrors(hRatioFEWZQCDScaleErrBandNorm);
  TGraphErrors* RatioFEWZScalePDFErrBandNorm = new TGraphErrors(hRatioFEWZScalePDFErrBandNorm);
  TGraphErrors *FEWZ_XsecNorm12 = new TGraphErrors(hFEWZ_Xsec_LogScaleNorm12);
  TGraphErrors* RatioFEWZStatErrBandNorm12 = new TGraphErrors(hRatioFEWZStatErrBandNorm12);
  TGraphErrors* RatioFEWZQCDScaleErrBandNorm12 = new TGraphErrors(hRatioFEWZQCDScaleErrBandNorm12);
  TGraphErrors* RatioFEWZScalePDFErrBandNorm12 = new TGraphErrors(hRatioFEWZScalePDFErrBandNorm12);


  /////==========FEWZ related stuff Finished here=========================================
  
  
  /////==========ResBos Starts here==============================================
  //
  // ResBos X sec value and errors, making TGraph .
 
  TH1D* hResBos30_CentralYield_LinScale;
  hResBos30_CentralYield_LinScale = (TH1D*)f_Resbos->Get("hResbos30")->Clone("hResBos30_CentralYield_LinScale");
  TH1D* hResBos30_CentralYield_LinScale12;
  hResBos30_CentralYield_LinScale12 = (TH1D*)f_Resbos_12->Get("hResbos30")->Clone("hResBos30_CentralYield_LinScale12");


  TH1D* AllResbos[7];
  for( int i(0);i<7;i++)
  {
    Numb = 29+i;
    sprintf(tmpName_org,"hResbos%d",Numb);
    sprintf(tmpName,"AllResbos_%d",i);
    AllResbos[i] = (TH1D*)f_Resbos->Get(tmpName_org)->Clone(tmpName);
  }
 
  //Errors
  Double_t Resb_errMax[nBins-1];
  Double_t Resb_errMin[nBins-1];
  double tmpVal,tmpDiff;

  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    double nomVal  = AllResbos[1]->GetBinContent(ipt+1);
    Resb_errMax[ipt] = -99999;
    Resb_errMin[ipt] = 990009;
    for (int i(0);i<7;i++)
    {
      tmpVal  = AllResbos[i]->GetBinContent(ipt+1);
      tmpDiff = tmpVal - nomVal;
      if( tmpDiff > Resb_errMax[ipt] ) Resb_errMax[ipt] = tmpDiff;
      if( tmpDiff < Resb_errMin[ipt] ) Resb_errMin[ipt] = tmpDiff;
    }

    if (Resb_errMax[ipt] < 0) Resb_errMax[ipt] = 0.;
    if (Resb_errMin[ipt] > 0) Resb_errMin[ipt] = 0.;
    if (Resb_errMin[ipt] < 0) Resb_errMin[ipt] = -Resb_errMin[ipt];
    Resb_errMax[ipt] = Resb_errMax[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
    Resb_errMin[ipt] = Resb_errMin[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
  }
  
  // X-sec and error Log scale
  TH1D *hResBos30_CentralXSec_LogScale   = new TH1D("hResBos30_CentralXSec_LogScale","hResBos30_CentralXSec_LogScale",nBins-1,WptLogBins);hResBos30_CentralXSec_LogScale->Sumw2();
  Double_t hResb30_CentralXsec[nBins-1];
  Double_t RatioResbVal[nBins-1],errResbosDataLo[nBins-1],errResbosDataHi[nBins-1];
  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    hResBos30_CentralXSec_LogScale->SetBinContent(ipt+1,hResBos30_CentralYield_LinScale->GetBinContent(ipt+1)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1) );
    hResb30_CentralXsec[ipt] = hResBos30_CentralYield_LinScale->GetBinContent(ipt+1)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
   
    RatioResbVal[ipt]=hResBos30_CentralXSec_LogScale->GetBinContent(ipt+1)/hData_Xsec_BornLogScale->GetBinContent(ipt+1);
    errResbosDataLo[ipt]=Resb_errMin[ipt]/hData_Xsec_BornLogScale->GetBinContent(ipt+1);
    errResbosDataHi[ipt]=Resb_errMax[ipt]/hData_Xsec_BornLogScale->GetBinContent(ipt+1);

  }
  ////12 bin
  TH1D* AllResbos12[7];
  for( int i(0);i<7;i++)
  {
    Numb = 29+i;
    sprintf(tmpName_org,"hResbos%d",Numb);
    sprintf(tmpName,"AllResbos12_%d",i);
    AllResbos12[i] = (TH1D*)f_Resbos_12->Get(tmpName_org)->Clone(tmpName);
  }
  Double_t Resb_errMax12[n12Bins-1];
  Double_t Resb_errMin12[n12Bins-1];
  double tmpVal12,tmpDiff12;

  for( int ipt(0);ipt<n12Bins-1;ipt++)
  {
    //double nomVal12  = AllResbos12[1]->GetBinContent(ipt+1);
    double nomVal12  = AllResbos12[1]->GetBinContent(ipt+1)/( AllResbos12[1]->Integral()*hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt+1));
    Resb_errMax12[ipt] = -99999;
    Resb_errMin12[ipt] = 990009;
    for (int i(0);i<7;i++)
    {
      //tmpVal12  = AllResbos12[i]->GetBinContent(ipt+1);
      tmpVal12  = AllResbos12[i]->GetBinContent(ipt+1)/( AllResbos12[i]->Integral()*hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt+1));
      tmpDiff12 = tmpVal12 - nomVal12;
      if( tmpDiff12 > Resb_errMax12[ipt] ) Resb_errMax12[ipt] = tmpDiff12;
      if( tmpDiff12 < Resb_errMin12[ipt] ) Resb_errMin12[ipt] = tmpDiff12;
    }

    if (Resb_errMax12[ipt] < 0) Resb_errMax12[ipt] = 0.;
    if (Resb_errMin12[ipt] > 0) Resb_errMin12[ipt] = 0.;
    if (Resb_errMin12[ipt] < 0) Resb_errMin12[ipt] = -Resb_errMin12[ipt];
    //Resb_errMax12[ipt] = Resb_errMax12[ipt]/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
    //Resb_errMin12[ipt] = Resb_errMin12[ipt]/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
    Resb_errMax12[ipt] = Resb_errMax12[ipt];
    Resb_errMin12[ipt] = Resb_errMin12[ipt];
  }
  
  //// 12 bin X-sec  and error Log scale
  TH1D *hResBos30_CentralXSec_LogScale12   = new TH1D("hResBos30_CentralXSec_LogScale12","hResBos30_CentralXSec_LogScale12",n12Bins-1,WptLog12Bins);hResBos30_CentralXSec_LogScale12->Sumw2();
  Double_t hResb30_CentralXsec12[n12Bins-1];
  Double_t RatioResbVal12[n12Bins-1],errResbosDataLo12[n12Bins-1],errResbosDataHi12[n12Bins-1];
  for( int ipt(0);ipt<n12Bins-1;ipt++)
  {
    hResBos30_CentralXSec_LogScale12->SetBinContent(ipt+1,hResBos30_CentralYield_LinScale12->GetBinContent(ipt+1)/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt+1) );
    hResb30_CentralXsec12[ipt] = hResBos30_CentralYield_LinScale12->GetBinContent(ipt+1)/hWpt12Bins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
   
    RatioResbVal12[ipt]=hResBos30_CentralXSec_LogScale12->GetBinContent(ipt+1)/hData_Xsec_BornLogScale12->GetBinContent(ipt+1);
    //errResbosDataLo12[ipt]=Resb_errMin12[ipt]/hData_Xsec_BornLogScale12->GetBinContent(ipt+1);
    //errResbosDataHi12[ipt]=Resb_errMax12[ipt]/hData_Xsec_BornLogScale12->GetBinContent(ipt+1);
    errResbosDataLo12[ipt]=Resb_errMin12[ipt];
    errResbosDataHi12[ipt]=Resb_errMax12[ipt];

  }
  
  ////Normalize ResBos X-sec and errors
  
  double resb30Total=0.0;
  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    //resb30Total=resb30Total + hResb30_CentralXsec[ipt];
    resb30Total=resb30Total + hResBos30_CentralYield_LinScale->GetBinContent(ipt+1);
  }
  double resb30Total12=0.0;
  for( int ipt(0);ipt<n12Bins-1;ipt++)
  {
    //resb30Total12=resb30Total12 + hResb30_CentralXsec12[ipt];
    resb30Total12=resb30Total12 + hResBos30_CentralYield_LinScale12->GetBinContent(ipt+1);
  }
  
  TH1D *hResBos30_CentralXSec_LogScaleNorm = (TH1D*)hResBos30_CentralXSec_LogScale->Clone("hResBos30_CentralXSec_LogScaleNorm");
  hResBos30_CentralXSec_LogScaleNorm->Scale(1./resb30Total);
  TH1D *hResBos30_CentralXSec_LogScaleNorm12 = (TH1D*)hResBos30_CentralXSec_LogScale12->Clone("hResBos30_CentralXSec_LogScaleNorm12");
  hResBos30_CentralXSec_LogScaleNorm12->Scale(1./resb30Total12);
  hResBos30_CentralXSec_LogScaleNorm12->Write();

  Double_t hResb30_CentralXsecNorm[nBins-1];
  Double_t Resb_errMaxNorm[nBins-1];
  Double_t Resb_errMinNorm[nBins-1];
  Double_t RatioResbValNorm[nBins-1],errResbosDataLoNorm[nBins-1],errResbosDataHiNorm[nBins-1];
  for( int ipt(0);ipt<nBins-1;ipt++)
  {

    hResb30_CentralXsecNorm[ipt] = hResb30_CentralXsec[ipt]/resb30Total;
    Resb_errMaxNorm[ipt] = Resb_errMax[ipt]/resb30Total;
    Resb_errMinNorm[ipt] = Resb_errMin[ipt]/resb30Total;
    
    RatioResbVal[ipt]=hResBos30_CentralXSec_LogScaleNorm->GetBinContent(ipt+1)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt+1);
    //RatioResbValNorm[ipt] = RatioResbVal[ipt]/resb30Total;
    errResbosDataLoNorm[ipt] =Resb_errMinNorm[ipt]/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt+1) ;
    errResbosDataHiNorm[ipt] =Resb_errMaxNorm[ipt]/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt+1) ;
  }
  Double_t hResb30_CentralXsecNorm12[n12Bins-1];
  Double_t Resb_errMaxNorm12[n12Bins-1];
  Double_t Resb_errMinNorm12[n12Bins-1];
  Double_t RatioResbValNorm12[n12Bins-1],errResbosDataLoNorm12[n12Bins-1],errResbosDataHiNorm12[n12Bins-1];
  for( int ipt(0);ipt<n12Bins-1;ipt++)
  {

    hResb30_CentralXsecNorm12[ipt] = hResb30_CentralXsec12[ipt]/resb30Total12;
    //Resb_errMaxNorm12[ipt] = Resb_errMax12[ipt]/resb30Total12;
    //Resb_errMinNorm12[ipt] = Resb_errMin12[ipt]/resb30Total12;
    Resb_errMaxNorm12[ipt] = Resb_errMax12[ipt];
    Resb_errMinNorm12[ipt] = Resb_errMin12[ipt];
    
    RatioResbVal12[ipt]=hResBos30_CentralXSec_LogScaleNorm12->GetBinContent(ipt+1)/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt+1);
    //RatioResbValNorm[ipt] = RatioResbVal[ipt]/resb30Total;
    
    RatioResbValNorm12[ipt]=hResBos30_CentralXSec_LogScaleNorm12->GetBinContent(ipt+1)/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt+1);
    errResbosDataLoNorm12[ipt] =Resb_errMinNorm12[ipt]/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt+1) ;
    errResbosDataHiNorm12[ipt] =Resb_errMaxNorm12[ipt]/hData_Xsec_BornLogScaleNorm12->GetBinContent(ipt+1) ;
  }

  ///   Theory/Data ratio plot errors related to Real Data
  TGraphAsymmErrors *Resb30_CentralXsec = new TGraphAsymmErrors(nBins-1, ax, hResb30_CentralXsec, aex, aex, Resb_errMin, Resb_errMax);
  TGraphAsymmErrors *RatioResbosErrBand = new TGraphAsymmErrors(nBins-1, ax, RatioResbVal, aex, aex, errResbosDataLo, errResbosDataHi);
  TGraphAsymmErrors *Resb30_CentralXsec12 = new TGraphAsymmErrors(n12Bins-1, ax12, hResb30_CentralXsec12, aex12, aex12, Resb_errMin12, Resb_errMax12);
  TGraphAsymmErrors *RatioResbosErrBand12 = new TGraphAsymmErrors(n12Bins-1, ax12, RatioResbVal12, aex12, aex12, errResbosDataLo12, errResbosDataHi12);
  ///   Normalized
  TGraphAsymmErrors *Resb30_CentralXsecNorm = new TGraphAsymmErrors(nBins-1, ax, hResb30_CentralXsecNorm, aex, aex, Resb_errMinNorm, Resb_errMaxNorm);
  TGraphAsymmErrors *RatioResbosErrBandNorm = new TGraphAsymmErrors(nBins-1, ax, RatioResbValNorm, aex, aex, errResbosDataLoNorm, errResbosDataHiNorm);
  TGraphAsymmErrors *Resb30_CentralXsecNorm12 = new TGraphAsymmErrors(n12Bins-1, ax12, hResb30_CentralXsecNorm12, aex12, aex12, Resb_errMinNorm12, Resb_errMaxNorm12);
  TGraphAsymmErrors *RatioResbosErrBandNorm12 = new TGraphAsymmErrors(n12Bins-1, ax12, RatioResbValNorm12, aex12, aex12, errResbosDataLoNorm12, errResbosDataHiNorm12);
  

  /////==========ResBos related stuff Finished here=========================================
  
  
  
  
  //// Now design and Draw 
  
  //Data

  RatioDataErrBand->SetFillColor(kBlack);
  RatioDataErrBand->SetFillStyle(3354);
  RatioDataErrBandNorm->SetFillColor(kBlack);
  RatioDataErrBandNorm->SetFillStyle(3354);
  RatioDataErrBand12->SetFillColor(kBlack);
  RatioDataErrBand12->SetFillStyle(3354);
  RatioDataErrBandNorm12->SetFillColor(kBlack);
  RatioDataErrBandNorm12->SetFillStyle(3354);
  
  //Powheg
  Powheg_Xsec_Born->SetFillColor(kRed);
  Powheg_Xsec_Born->SetFillStyle(3345);
  Powheg_Xsec_BornNorm->SetFillColor(kRed);
  Powheg_Xsec_BornNorm->SetFillStyle(3345);
  Powheg_Xsec_Born12->SetFillColor(kRed);
  Powheg_Xsec_Born12->SetFillStyle(3345);
  Powheg_Xsec_BornNorm12->SetFillColor(kRed);
  Powheg_Xsec_BornNorm12->SetFillStyle(3345);

  RatioPowhegStatErrBand->SetFillColor(kRed-10);
  RatioPowhegStatErrBand->SetFillStyle(3001);
  RatioPowhegStatErrBandNorm->SetFillColor(kRed-10);
  RatioPowhegStatErrBandNorm->SetFillStyle(3001);
  RatioPowhegStatErrBand12->SetFillColor(kRed-10);
  RatioPowhegStatErrBand12->SetFillStyle(3001);
  RatioPowhegStatErrBandNorm12->SetFillColor(kRed-10);
  RatioPowhegStatErrBandNorm12->SetFillStyle(3001);
  
  RatioPowhegPDFErrBand->SetFillColor(kRed+2);
  RatioPowhegPDFErrBand->SetFillStyle(3001);
  RatioPowhegPDFErrBandNorm->SetFillColor(kRed+2);
  RatioPowhegPDFErrBandNorm->SetFillStyle(3001);
  RatioPowhegPDFErrBand12->SetFillColor(kRed+2);
  RatioPowhegPDFErrBand12->SetFillStyle(3001);
  RatioPowhegPDFErrBandNorm12->SetFillColor(kRed+2);
  RatioPowhegPDFErrBandNorm12->SetFillStyle(3001);

  //FEWZ
  FEWZ_Xsec->SetFillColor(kGreen);
  FEWZ_Xsec->SetFillStyle(3444);
  FEWZ_XsecNorm->SetFillColor(kGreen);
  FEWZ_XsecNorm->SetFillStyle(3444);
  FEWZ_Xsec12->SetFillColor(kGreen);
  FEWZ_Xsec12->SetFillStyle(3444);
  FEWZ_XsecNorm12->SetFillColor(kGreen);
  FEWZ_XsecNorm12->SetFillStyle(3444);
  
  RatioFEWZStatErrBand->SetFillColor(kGreen);
  RatioFEWZStatErrBand->SetFillStyle(3001);
  RatioFEWZStatErrBandNorm->SetFillColor(kGreen);
  RatioFEWZStatErrBandNorm->SetFillStyle(3001);
  RatioFEWZStatErrBand12->SetFillColor(kGreen);
  RatioFEWZStatErrBand12->SetFillStyle(3001);
  RatioFEWZStatErrBandNorm12->SetFillColor(kGreen);
  RatioFEWZStatErrBandNorm12->SetFillStyle(3001);

  RatioFEWZQCDScaleErrBand->SetFillColor(kGreen+7);
  RatioFEWZQCDScaleErrBand->SetFillStyle(3001);
  RatioFEWZQCDScaleErrBandNorm->SetFillColor(kGreen+7);
  RatioFEWZQCDScaleErrBandNorm->SetFillStyle(3001);
  RatioFEWZQCDScaleErrBand12->SetFillColor(kGreen+7);
  RatioFEWZQCDScaleErrBand12->SetFillStyle(3001);
  RatioFEWZQCDScaleErrBandNorm12->SetFillColor(kGreen+7);
  RatioFEWZQCDScaleErrBandNorm12->SetFillStyle(3001);

  RatioFEWZScalePDFErrBand->SetFillColor(kGreen+3);
  RatioFEWZScalePDFErrBand->SetFillStyle(3001);
  RatioFEWZScalePDFErrBandNorm->SetFillColor(kGreen+3);
  RatioFEWZScalePDFErrBandNorm->SetFillStyle(3001);
  RatioFEWZScalePDFErrBand12->SetFillColor(kGreen+3);
  RatioFEWZScalePDFErrBand12->SetFillStyle(3001);
  RatioFEWZScalePDFErrBandNorm12->SetFillColor(kGreen+3);
  RatioFEWZScalePDFErrBandNorm12->SetFillStyle(3001);

  
  //ResBos
  Resb30_CentralXsec->SetFillColor(kBlue);
  Resb30_CentralXsec->SetFillStyle(3354);
  Resb30_CentralXsecNorm->SetFillColor(kBlue);
  Resb30_CentralXsecNorm->SetFillStyle(3354);
  Resb30_CentralXsec12->SetFillColor(kBlue);
  Resb30_CentralXsec12->SetFillStyle(3354);
  Resb30_CentralXsecNorm12->SetFillColor(kBlue);
  Resb30_CentralXsecNorm12->SetFillStyle(3354);
 
  RatioResbosErrBand->SetFillColor(kBlue-7);
  RatioResbosErrBand->SetFillStyle(3001);
  RatioResbosErrBandNorm->SetFillColor(kBlue-7);
  RatioResbosErrBandNorm->SetFillStyle(3001);
  RatioResbosErrBand12->SetFillColor(kBlue-7);
  RatioResbosErrBand12->SetFillStyle(3001);
  RatioResbosErrBandNorm12->SetFillColor(kBlue-7);
  RatioResbosErrBandNorm12->SetFillStyle(3001);
 

  TLegend *lL =new TLegend(0.2,0.2,0.52,0.38); lL->SetFillColor(0); lL->SetBorderSize(0);
  //lL->AddEntry(Data_Xsec_Born,"data","PL");
  //lL->AddEntry(Powheg_Xsec_Born,"POWHEG CT10 NLO","f");
  //lL->AddEntry(FEWZ_Xsec,"FEWZ CT10 NNLO","f");
  //lL->AddEntry(Resb30_CentralXsec,"ResBos CT10 NNLL","f");
  //lL->AddEntry(Data_Xsec_BornNorm,"data","PL");
  //lL->AddEntry(Powheg_Xsec_BornNorm,"POWHEG CT10 NLO","f");
  //lL->AddEntry(FEWZ_XsecNorm,"FEWZ CT10 NNLO","f");
  //lL->AddEntry(Resb30_CentralXsecNorm,"ResBos CT10 NNLL","f");
  lL->AddEntry(Data_Xsec_BornNorm12,"data","PL");
  lL->AddEntry(Powheg_Xsec_BornNorm12,"POWHEG CT10 NLO","f");
  lL->AddEntry(FEWZ_XsecNorm12,"FEWZ CT10 NNLO","f");
  lL->AddEntry(Resb30_CentralXsecNorm12,"ResBos CT10 NNLL","f");

  TPaveText *tb = new TPaveText(0.2,0.39,0.52,0.5,"NDC");
  tb->SetBorderSize(0);
  tb->SetFillStyle(0);
  tb->AddText("18.4 pb^{-1} at #sqrt{s} = 8 TeV");
  if (BaseName=="WpToMuNu")
    tb->AddText("W^{+} #rightarrow #mu^{+} #nu");
  if (BaseName=="WmToMuNu")
    tb->AddText("W^{-} #rightarrow #mu^{-} #bar{#nu}");
  if (BaseName=="WpToEleNu")
    tb->AddText("W^{+} #rightarrow e^{+} #nu");
  if (BaseName=="WmToEleNu")
    tb->AddText("W^{-} #rightarrow e^{-} #bar{#nu}");

  TCanvas *lC1 = new TCanvas("Can","Can",800,840); lC1->cd(); lC1->SetLogy();
  lC1->cd(1)->SetPad(0,0.05,0.95,1.0);
  lC1->cd(1)->SetTopMargin(0.05);
  lC1->cd(1)->SetBottomMargin(0.1);
  lC1->cd(1)->SetLeftMargin(0.15);
  lC1->cd(1)->SetRightMargin(0.07);
  lC1->cd(1)->SetTickx(1);
  lC1->cd(1)->SetTicky(1);
  gStyle->SetLineWidth(2.);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesSpacing(0.75);
  gStyle->SetHatchesLineWidth(2);
  gPad->SetLogx(1);
  gPad->SetLogy(1);

  ////Powheg_Xsec_Born->GetYaxis()->SetRangeUser(1e-3,5*hResb30_CentralXsec[0]);
  //Powheg_Xsec_Born->GetYaxis()->SetRangeUser(1e-6,5*hResb30_CentralXsec[0]);
  //Powheg_Xsec_Born->SetTitle("");
  //Powheg_Xsec_Born->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T}^{W} [pb/GeV]");
  //Powheg_Xsec_Born->GetYaxis()->SetTitleOffset(1.8);
  //Powheg_Xsec_Born->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  
  //Powheg_Xsec_Born->GetYaxis()->SetRangeUser(1e-3,5*hResb30_CentralXsec[0]);
  //Powheg_Xsec_BornNorm->GetYaxis()->SetRangeUser(1e-6,5*hResb30_CentralXsecNorm[0]);
  //Powheg_Xsec_BornNorm->SetTitle("");
  //Powheg_Xsec_BornNorm->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T}^{W} [pb/GeV]");
  //Powheg_Xsec_BornNorm->GetYaxis()->SetTitleOffset(1.8);
  //Powheg_Xsec_BornNorm->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  Powheg_Xsec_BornNorm12->GetYaxis()->SetRangeUser(1e-7,5*hResb30_CentralXsecNorm[0]);
  Powheg_Xsec_BornNorm12->SetTitle("");
  Powheg_Xsec_BornNorm12->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T}^{W} [pb/GeV]");
  Powheg_Xsec_BornNorm12->GetYaxis()->SetTitleOffset(1.8);
  Powheg_Xsec_BornNorm12->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  
  TPaveText *cmspre = new TPaveText(0.63,0.95,0.95,0.99,"NDC");
  cmspre->SetBorderSize(0);
  cmspre->SetFillStyle(0);
  cmspre->AddText("CMS Preliminary");
  //Powheg_Xsec_Born->Draw("A2");
  //FEWZ_Xsec->Draw("2");
  //Resb30_CentralXsec->Draw("2");
  //Data_Xsec_Born->Draw("p");
  //Powheg_Xsec_BornNorm->Draw("A2");
  //FEWZ_XsecNorm->Draw("2");
  //Resb30_CentralXsecNorm->Draw("2");
  //Data_Xsec_BornNorm->Draw("p");
  Powheg_Xsec_BornNorm12->Draw("A2");
  FEWZ_XsecNorm12->Draw("2");
  Resb30_CentralXsecNorm12->Draw("2");
  Data_Xsec_BornNorm12->Draw("p");

  lL->Draw();
  tb->Draw();
  cmspre->Draw();
  
  if(BaseName=="WpToMuNu")
    sprintf(tmpName,"wpmnCrossSecNormInFiducial.png");
  if(BaseName=="WmToMuNu")
    sprintf(tmpName,"wmmnCrossSecNormInFiducial.png");
  if(BaseName=="WpToEleNu")
    sprintf(tmpName,"wpenCrossSecNormInFiducial.png");
  if(BaseName=="WmToEleNu")
    sprintf(tmpName,"wmenCrossSecNormInFiducial.png");
  lC1->SaveAs(tmpName);
  
  TCanvas *lC2 = new TCanvas("Can","Can",800,800); lC2->cd(); lC2->SetLogy();
  lC2->Divide(1,3,0,0);
  lC2->cd(1)->SetPad(0,0.67,0.95,0.95);
  lC2->cd(1)->SetTopMargin(0.15);
  lC2->cd(1)->SetBottomMargin(0.01);
  lC2->cd(1)->SetLeftMargin(0.15);
  lC2->cd(1)->SetRightMargin(0.07);
  lC2->cd(1)->SetTickx(1);
  lC2->cd(1)->SetTicky(1);
  lC2->cd(1)->SetLogx(1);

  TPaveText *tb1 = new TPaveText(0.15,0.72,0.35,0.82,"NDC");
  tb1->SetBorderSize(0);
  tb1->SetFillStyle(0);
  tb1->SetTextSize(0.12);
  tb1->AddText("ResBos");
  TLegend *rL1 =new TLegend(0.2,0.1,0.58,0.30); rL1->SetFillColor(0); rL1->SetBorderSize(0);
  rL1-> SetNColumns(2);
  //rL1->AddEntry(RatioResbosErrBand,"Theory syst","F");
  //rL1->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  //hRatioDataStatErr->SetTitle("");
  //rL1->AddEntry(hRatioDataStatErr,"","");
  //rL1->AddEntry(RatioDataErrBand,"Data stat+syst","F");
  //rL1->AddEntry(RatioResbosErrBand12,"Theory syst","F");
  rL1->AddEntry(RatioResbosErrBandNorm12,"Theory syst","F");
  rL1->AddEntry(hRatioDataStatErr12,"Data stat","PLE1");
  hRatioDataStatErr12->SetTitle("");
  rL1->AddEntry(hRatioDataStatErr12,"","");
  rL1->AddEntry(RatioDataErrBand12,"Data stat+syst","F");
  ////Normalized RatioBand
  //rL1->AddEntry(RatioResbosErrBandNorm,"Theory syst","F");
  //rL1->AddEntry(hRatioDataStatErrNorm,"Data stat","PLE1");
  //hRatioDataStatErrNorm->SetTitle("");
  //rL1->AddEntry(hRatioDataStatErrNorm,"","");
  //rL1->AddEntry(RatioDataErrBandNorm,"Data stat+syst","F");
  rL1->SetTextSize(0.07);
  
  TLegend *tL1 =new TLegend(0.17,0.72,0.37,0.82); tL1->SetFillColor(0); tL1->SetBorderSize(0);
  //tL1->AddEntry(RatioResbosErrBand,"ResBos","F");
  //tL1->AddEntry(RatioResbosErrBand12,"ResBos","F");
  tL1->AddEntry(RatioResbosErrBandNorm12,"ResBos","F");
  ////Normalized RatioBand
  //tL1->AddEntry(RatioResbosErrBandNorm,"ResBos","F");
  tL1->SetTextSize(0.12);
  tL1->SetTextFont(2);

  TPaveText *cmspre2 = new TPaveText(0.05,0.92,0.85,0.99,"NDC");
  cmspre2->SetBorderSize(0);
  cmspre2->SetFillStyle(0);
  cmspre2->SetTextSize(0.10);
  cmspre2->AddText("CMS Preliminary, 18.4 pb^{-1} at #sqrt{s} = 8 TeV");

  TPaveText *tb4 = new TPaveText(0.35,0.72,0.67,0.83,"NDC");
  tb4->SetBorderSize(0);
  tb4->SetFillStyle(0);
  if (BaseName=="WpToMuNu")
    tb4->AddText("W^{+} #rightarrow #mu^{+} #nu");
  if (BaseName=="WmToMuNu")
    tb4->AddText("W^{-} #rightarrow #mu^{-} #bar{#nu}");
  if (BaseName=="WpToEleNu")
    tb4->AddText("W^{+} #rightarrow e^{+} #nu");
  if (BaseName=="WmToEleNu")
    tb4->AddText("W^{-} #rightarrow e^{-} #bar{#nu}");
  
  
  //// NOT Normalized RatioBand
  //drawDifference(hResBos30_CentralXSec_LogScale, hData_Xsec_BornLogScale, hRatioDataErrBand, RatioPowhegStatErrBand,1,RatioPowhegStatErrBand,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
  ////Normalized All RatioBand
  //drawDifference(hResBos30_CentralXSec_LogScaleNorm, hData_Xsec_BornLogScaleNorm, hRatioDataErrBandNorm, RatioPowhegStatErrBandNorm,1,RatioPowhegStatErrBandNorm,RatioResbosErrBandNorm,hRatioDataStatErrNorm,RatioFEWZQCDScaleErrBandNorm);
  ////Normalized only Xsec
  //drawDifference(hResBos30_CentralXSec_LogScaleNorm, hData_Xsec_BornLogScaleNorm, hRatioDataErrBand, RatioPowhegStatErrBand,1,RatioPowhegStatErrBand,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
  ////Normalized only Xsec  12 Bin
  //drawDifference(hResBos30_CentralXSec_LogScaleNorm12, hData_Xsec_BornLogScaleNorm12, hRatioDataErrBand12, RatioPowhegStatErrBand12,1,RatioPowhegStatErrBand12,RatioResbosErrBand12,hRatioDataStatErr12,RatioFEWZQCDScaleErrBand12);
  drawDifference(hResBos30_CentralXSec_LogScaleNorm12, hData_Xsec_BornLogScaleNorm12, hRatioDataErrBand12, RatioPowhegStatErrBand12,1,RatioPowhegStatErrBand12,RatioResbosErrBandNorm12,hRatioDataStatErr12,RatioFEWZQCDScaleErrBand12);
  rL1->Draw();
  //tb1->Draw();
  tL1->Draw();
  tb4->Draw();
  cmspre2->Draw();

  lC2->cd(2)->SetPad(0,0.37,0.95,0.65);
  lC2->cd(2)->SetTopMargin(0.01);
  lC2->cd(2)->SetBottomMargin(0.01);
  lC2->cd(2)->SetLeftMargin(0.15);
  lC2->cd(2)->SetRightMargin(0.07);
  lC2->cd(2)->SetTickx(1);
  lC2->cd(2)->SetTicky(1);
  lC2->cd(2)->SetLogx(1);

  TPaveText *tb2 = new TPaveText(0.15,0.82,0.35,0.92,"NDC");
  tb2->SetBorderSize(0);
  tb2->SetFillStyle(0);
  tb2->SetTextSize(0.12);
  tb2->AddText("Powheg");
  TLegend *rL2 =new TLegend(0.2,0.05,0.68,0.30); rL2->SetFillColor(0); rL2->SetBorderSize(0);
  rL2-> SetNColumns(2);
  //rL2->AddEntry(RatioPowhegPDFErrBand,"PDF    ","F");
  //rL2->AddEntry(RatioPowhegPDFErrBandNorm,"PDF    ","F");
  rL2->AddEntry(RatioPowhegPDFErrBandNorm12,"PDF    ","F");
  //rL2->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  rL2->AddEntry(hRatioDataStatErr12,"Data stat","PLE1");
  //rL2->AddEntry(hRatioDataStatErrNorm,"Data stat","PLE1");
  //rL2->AddEntry(RatioPowhegStatErrBand,"stat","F");
  //rL2->AddEntry(RatioPowhegStatErrBandNorm,"stat","F");
  rL2->AddEntry(RatioPowhegStatErrBandNorm12,"stat","F");
  //rL2->AddEntry(RatioDataErrBand,"Data stat+syst","F");
  rL2->AddEntry(RatioDataErrBand12,"Data stat+syst","F");
  //rL2->AddEntry(RatioDataErrBandNorm,"Data stat+syst","F");
  rL2->SetTextSize(0.07);

  TLegend *tL2 =new TLegend(0.17,0.85,0.37,0.95); tL2->SetFillColor(0); tL2->SetBorderSize(0);
  //tL2->AddEntry(RatioPowhegPDFErrBand,"POWHEG","F");
  ////Normalized RatioBand
  //tL2->AddEntry(RatioPowhegPDFErrBandNorm,"POWHEG","F");
  tL2->AddEntry(RatioPowhegPDFErrBandNorm12,"POWHEG","F");
  tL2->SetTextSize(0.12);
  tL2->SetTextFont(2);
  
  //// NOT Normalized RatioBand
  //drawDifference(hPowheg_Xsec_BornLogScale,hData_Xsec_BornLogScale,hRatioDataErrBand,RatioPowhegStatErrBand,2,RatioPowhegPDFErrBand,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
  ////Normalized All RatioBand
  //drawDifference(hPowheg_Xsec_BornLogScaleNorm,hData_Xsec_BornLogScaleNorm,hRatioDataErrBandNorm,RatioPowhegStatErrBandNorm,2,RatioPowhegPDFErrBandNorm,RatioResbosErrBandNorm,hRatioDataStatErrNorm,RatioFEWZQCDScaleErrBandNorm);
  ////Normalized Xsec, RatioPowhegStatErrBandNorm, RatioPowhegPDFErrBandNorm
  //drawDifference(hPowheg_Xsec_BornLogScaleNorm,hData_Xsec_BornLogScaleNorm,hRatioDataErrBand,RatioPowhegStatErrBandNorm,2,RatioPowhegPDFErrBandNorm,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
  ////   12Normalized Xsec, RatioPowhegStatErrBandNorm, RatioPowhegPDFErrBandNorm
  drawDifference(hPowheg_Xsec_BornLogScaleNorm12,hData_Xsec_BornLogScaleNorm12,hRatioDataErrBand12,RatioPowhegStatErrBandNorm12,2,RatioPowhegPDFErrBandNorm12,RatioResbosErrBand12,hRatioDataStatErr12,RatioFEWZQCDScaleErrBand12);
  rL2->Draw();
  tL2->Draw();

  lC2->cd(3)->SetPad(0,0.07,0.95,0.35);
  lC2->cd(3)->SetTopMargin(0.01);
  lC2->cd(3)->SetBottomMargin(0.05);
  lC2->cd(3)->SetLeftMargin(0.15);
  lC2->cd(3)->SetRightMargin(0.07);
  lC2->cd(3)->SetTickx(1);
  lC2->cd(3)->SetTicky(1);
  lC2->cd(3)->SetLogx(1);

  TPaveText *tb3 = new TPaveText(0.15,0.82,0.35,0.92,"NDC");
  tb3->SetBorderSize(0);
  tb3->SetFillStyle(0);
  tb3->SetTextSize(0.12);
  tb3->AddText("Fewz");
  TLegend *rL3 =new TLegend(0.2,0.1,0.58,0.30); rL3->SetFillColor(0); rL3->SetBorderSize(0);
  rL3-> SetNColumns(2);
  //rL3->AddEntry(RatioFEWZScalePDFErrBand,"PDF","F");
  //rL3->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  //rL3->AddEntry(RatioFEWZQCDScaleErrBand,"QCD scales","F");
  //rL3->AddEntry(RatioDataErrBand,"Data stat+syst","F");
  //rL3->AddEntry(RatioFEWZStatErrBand,"stat","F");
  //rL3->AddEntry(RatioFEWZScalePDFErrBandNorm,"PDF","F");
  //rL3->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  //rL3->AddEntry(RatioFEWZQCDScaleErrBandNorm,"QCD scales","F");
  //rL3->AddEntry(RatioDataErrBand,"Data stat+syst","F");
  //rL3->AddEntry(RatioFEWZStatErrBandNorm,"stat","F");
  rL3->AddEntry(RatioFEWZScalePDFErrBandNorm12,"PDF","F");
  rL3->AddEntry(hRatioDataStatErr12,"Data stat","PLE1");
  rL3->AddEntry(RatioFEWZQCDScaleErrBandNorm12,"QCD scales","F");
  rL3->AddEntry(RatioDataErrBand12,"Data stat+syst","F");
  rL3->AddEntry(RatioFEWZStatErrBandNorm12,"stat","F");
  rL3->SetTextSize(0.07);

  TLegend *tL3 =new TLegend(0.17,0.85,0.37,0.95); tL3->SetFillColor(0); tL3->SetBorderSize(0);
  //tL3->AddEntry(RatioFEWZScalePDFErrBand,"FEWZ","F");
  //tL3->AddEntry(RatioFEWZScalePDFErrBandNorm,"FEWZ","F");
  tL3->AddEntry(RatioFEWZScalePDFErrBandNorm12,"FEWZ","F");
  tL3->SetTextSize(0.12);
  tL3->SetTextFont(2);
  
  //drawDifference(hFEWZ_Xsec_LogScale,hData_Xsec_BornLogScale,hRatioDataErrBand,RatioFEWZStatErrBand,3,RatioFEWZScalePDFErrBand,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
  //Normalized Xsec, RatioFEWZStatErrBandNorm, RatioFEWZScalePDFErrBandNorm, RatioFEWZQCDScaleErrBandNorm
  //drawDifference(hFEWZ_Xsec_LogScaleNorm,hData_Xsec_BornLogScaleNorm,hRatioDataErrBand,RatioFEWZStatErrBandNorm,3,RatioFEWZScalePDFErrBandNorm,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBandNorm);
  ///12 bin
  drawDifference(hFEWZ_Xsec_LogScaleNorm12,hData_Xsec_BornLogScaleNorm12,hRatioDataErrBand12,RatioFEWZStatErrBandNorm12,3,RatioFEWZScalePDFErrBandNorm12,RatioResbosErrBand12,hRatioDataStatErr12,RatioFEWZQCDScaleErrBandNorm12);
  rL3->Draw();
  tL3->Draw();

  if(BaseName=="WpToMuNu")
    sprintf(tmpName,"wpmnRatioTheoDataNormInFiducial.png");
  if(BaseName=="WmToMuNu")
    sprintf(tmpName,"wmmnRatioTheoDataNormInFiducial.png");
  if(BaseName=="WpToEleNu")
    sprintf(tmpName,"wpenRatioTheoDataNormInFiducial.png");
  if(BaseName=="WmToEleNu")
    sprintf(tmpName,"wmenRatioTheoDataNormInFiducial.png");
  lC2->SaveAs(tmpName);
    
  return 0;
}
