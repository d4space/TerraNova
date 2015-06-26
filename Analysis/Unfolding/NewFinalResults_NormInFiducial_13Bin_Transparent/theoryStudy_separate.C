#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "../../Utils/const.h"

const int nBins = 14;
double WptLogBins[nBins] = {1.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};
double WptBins[nBins] = {0.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};
double BinWidth[14] ={0.0, 7.5-0, 12.5-7.5,17.5-12.5, 24.0-17.5, 30.0-24.0, 40.0-30.0, 50.0-40.0, 70.0-50.0, 110.0-70.0, 150.0-110.0, 190.0-150.0, 250.0-190.0, 600.0-250.0};

double ax[13]  = {4.25,10,15,20.75,27,35,45,60,90,130,170,220,425};
double aex[13] = {3.25,2.5,2.5,3.25,3,5,5,10,20,20,20,30,175};

void Powheg_NormPDFUncer(const TString BaseName, double *Xsec, double *Err);
void Powheg_NormDiffStatErr(const TString BaseName, double *Err);
void FEWZ_NormPDFUncer(const TString BaseName, double *Xsec, double *Err);
void NormDiffToyErr(double *_Xsec, double *_Err, double *_NormDiff_Err);

int theoryStudy_separate(const TString BaseName)
{
  TString tmpTStr;
  char tmpName[30],tmpName_org[30];
  int Numb;

  TFile *f_Resbos;
  TFile *f_Fewz;
  TFile *f_Data;
  TFile *f_Data_NormSyst;

  f_Resbos = new TFile("../../RstResbos/Resbos_"+BaseName+".root");
  f_Data = new TFile("../Result"+BaseName+"/Result_"+BaseName+".root");
  f_Data_NormSyst = new TFile("../../Systematics_NormDiffXsec/WptXsecErrors/"+BaseName+"Errors.root");
  
  if (BaseName=="WpToMuNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WpToMuNu_13bin_dynamic_NNLO.root");
  }
  if (BaseName=="WmToMuNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WmToMuNu_13bin_dynamic_NNLO.root");
  }
  if (BaseName=="WpToEleNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WpToEleNu_13bin_dynamic_NNLO.root");
  }
  if (BaseName=="WmToEleNu")
  {
    f_Fewz = new TFile("../../RstFEWZ/WmToEleNu_13bin_dynamic_NNLO.root");
  }


  TH1D *hWptBins_LinScale   = new TH1D("hWptBins_LinScale","hWptBins_LinScale",nBins-1,WptBins);
  TH1D *hWptBins_LogScale   = new TH1D("hWptBins_LogScale","hWptBins_LogScale",nBins-1,WptLogBins);
  

  ///=============Reading related stat, syst errors starts here==============
  //  
  TH1D* hData_StatErr	 = (TH1D*)f_Data_NormSyst->Get("h_Stat")->Clone("hData_StatErr");
  
  TH1D* hData_MetResolSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_met")->Clone("hData_MetResolSystErr");
  TH1D* hData_EnMomScaleSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_scale")->Clone("hData_EnMomScaleSystErr");
  TH1D* hData_EnMomSmearSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_smear")->Clone("hData_EnMomSmearSystErr");
  TH1D* hData_QcdBckgrSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_qcdbckgr")->Clone("hData_QcdBckgrSystErr");
  TH1D* hData_QcdShapeSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_qcdshape")->Clone("hData_QcdShapeSystErr");
  TH1D* hData_EwkSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_ewk")->Clone("hData_EwkSystErr");
  TH1D* hData_FsrSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_fsr")->Clone("hData_FsrSystErr");
  TH1D* hData_SvdUnfSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_SvdUnf")->Clone("hData_SvdUnfSystErr");
  TH1D* hData_UnfoldBiasSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_UnfoldBias")->Clone("hData_UnfoldBiasSystErr");
  TH1D* hData_LumiSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_LumiSyst")->Clone("hData_LumiSystErr");
  
  ///Lepton Reconstruction related systematic errors
  //
  TH1D* hData_EffiToySystErr		 = (TH1D*)f_Data_NormSyst->Get("h_toy")->Clone("hData_EffiToySystErr");
  TH1D* hData_IDIsoSigShapeSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_idisosig")->Clone("hData_IDIsoSigShapeSystErr");
  TH1D* hData_IDIsoBkgrShapeSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_idisobck")->Clone("hData_IDIsoBkgrShapeSystErr");
  
  TH1D *hData_TotalSystErr   = new TH1D("hData_TotalSystErr","hData_TotalSystErr",nBins-1,WptBins);
  
  if(BaseName== "WpToMuNu" || BaseName== "WmToMuNu")
  {
    TH1D* hData_TrackSigShapeSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_tracksig")->Clone("hData_TrackSigShapeSystErr");
    TH1D* hData_TrackBkgrShapeSystErr	 = (TH1D*)f_Data_NormSyst->Get("h_trackbck")->Clone("hData_TrackBkgrShapeSystErr");
    TH1D* hData_MuonPOGSystErr		 = (TH1D*)f_Data_NormSyst->Get("h_POG")->Clone("hData_MuonPOGSystErr");
   
    //// Calculate total syst for muon
    for( int ipt(1);ipt<=nBins-1;ipt++)
    {
      hData_TotalSystErr->SetBinContent(ipt, sqrt(
	    TMath::Power(hData_TrackSigShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_TrackBkgrShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_MuonPOGSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_EffiToySystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_IDIsoSigShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_MetResolSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_EnMomScaleSystErr->GetBinContent(ipt),2) +
	    TMath::Power(hData_EnMomSmearSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_QcdBckgrSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_QcdShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_EwkSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_FsrSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_SvdUnfSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_UnfoldBiasSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_LumiSystErr->GetBinContent(ipt),2) ));
    }
    
  }
  else if(BaseName== "WpToEleNu" || BaseName== "WmToEleNu")
  {
    TH1D* hData_BinningSystErr		 = (TH1D*)f_Data_NormSyst->Get("h_bin")->Clone("hData_BinningSystErr");
  
    //// Calculate total syst for electron
    for( int ipt(1);ipt<=nBins-1;ipt++)
    {
      hData_TotalSystErr->SetBinContent(ipt, sqrt(
	    TMath::Power(hData_BinningSystErr->GetBinContent(ipt),2) + 
	    TMath::Power( hData_EffiToySystErr->GetBinContent(ipt),2) + 
	    TMath::Power( hData_IDIsoSigShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_MetResolSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_EnMomScaleSystErr->GetBinContent(ipt),2) +
	    TMath::Power(hData_EnMomSmearSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_QcdBckgrSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_QcdShapeSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_EwkSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_FsrSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_SvdUnfSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_UnfoldBiasSystErr->GetBinContent(ipt),2) +
	    TMath::Power( hData_LumiSystErr->GetBinContent(ipt),2) ));
    }
  }

  TH1D *hData_TotalSyst   = new TH1D("hData_TotalSyst","hData_TotalSyst",nBins-1,WptBins);
  TH1D *hData_TotalUncer   = new TH1D("hData_TotalUncer","hData_TotalUncer",nBins-1,WptBins);
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hData_TotalSyst->SetBinContent(ipt, hData_TotalSystErr->GetBinContent(ipt));
    hData_TotalUncer->SetBinContent(ipt, sqrt(TMath::Power(hData_StatErr->GetBinContent(ipt),2)+ TMath::Power(hData_TotalSystErr->GetBinContent(ipt),2)));
    cout<<"hData_TotalUncer: "<<ipt<<"\t"<<hData_TotalUncer->GetBinContent(ipt)<<endl;
  } 
 

  for( int ipt(1);ipt<nBins;ipt++)
  {
    //cout<<"hData_TrackSigShapeSystErr: "<<ipt<<"\t"<<hData_TrackSigShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_TrackBkgrShapeSystErr: "<<ipt<<"\t"<<hData_TrackBkgrShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_IDIsoSigShapeSystErr: "<<ipt<<"\t"<<hData_IDIsoSigShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_IDIsoBkgrShapeSystErr: "<<ipt<<"\t"<<hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_MuonPOGSystErr: "<<ipt<<"\t"<<hData_MuonPOGSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_EffiToySystErr: "<<ipt<<"\t"<<hData_EffiToySystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_MetResolSystErr: "<<ipt<<"\t"<<hData_MetResolSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_EnMomScaleSystErr: "<<ipt<<"\t"<<hData_EnMomScaleSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_EnMomSmearSystErr: "<<ipt<<"\t"<<hData_EnMomSmearSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_QcdBckgrSystErr: "<<ipt<<"\t"<<hData_QcdBckgrSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_QcdShapeSystErr: "<<ipt<<"\t"<<hData_QcdShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_EwkSystErr: "<<ipt<<"\t"<<hData_EwkSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_FsrSystErr: "<<ipt<<"\t"<<hData_FsrSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_SvdUnfSystErr: "<<ipt<<"\t"<<hData_SvdUnfSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_UnfoldBiasSystErr: "<<ipt<<"\t"<<hData_UnfoldBiasSystErr->GetBinContent(ipt)<<endl;
  } 

  ///=============Reading related stat, syst errors Finished here==============
 




  /////==========Real Data Starts here=========================================
  //
  // Born Level RD Yield and Errors, and Converting to X sec yield and errors, making TGraph .

  TH1D* hData_Xsec;
  hData_Xsec   = (TH1D*)f_Data->Get("BornEffCorr")->Clone("hData_Xsec");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hData_Xsec->SetBinError(ipt, 0.01*hData_TotalUncer->GetBinContent(ipt)*hData_Xsec->GetBinContent(ipt) );
  } 
  hData_Xsec->Scale(1./18.429);

  // Data Total Xsec
  double Data_TotalXsec;
  Data_TotalXsec = hData_Xsec->Integral();

  //Errors
  double Data_Xsec_StatErr[14];
  double Data_Xsec_SystErr[14];
  double Data_Xsec_TotalErr[14];
  double Data_NormDiffXsec_StatErr[14];
  double Data_NormDiffXsec_SystErr[14];
  double Data_NormDiffXsec_TotalErr[14];
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    Data_Xsec_StatErr[ipt] = hData_Xsec->GetBinContent(ipt) * hData_StatErr->GetBinContent(ipt) * 0.01;
    Data_Xsec_SystErr[ipt] = hData_Xsec->GetBinContent(ipt) * hData_TotalSystErr->GetBinContent(ipt) * 0.01;
    Data_Xsec_TotalErr[ipt] = hData_Xsec->GetBinError(ipt);
    Data_NormDiffXsec_StatErr[ipt] = Data_Xsec_StatErr[ipt] / BinWidth[ipt] / Data_TotalXsec;
    Data_NormDiffXsec_SystErr[ipt] = Data_Xsec_SystErr[ipt] / BinWidth[ipt] / Data_TotalXsec;
    Data_NormDiffXsec_TotalErr[ipt] = Data_Xsec_TotalErr[ipt] / BinWidth[ipt] / Data_TotalXsec;
  }

  // Differential X-sec and error 
  TH1D *hData_DiffXsec     = new TH1D("hData_DiffXsec","hData_DiffXsec",nBins-1,WptLogBins);
  TH1D *hData_NormDiffXsec     = new TH1D("hData_NormDiffXsec","hData_NormDiffXsec",nBins-1,WptLogBins);
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hData_DiffXsec->SetBinContent(ipt,hData_Xsec->GetBinContent(ipt)/BinWidth[ipt]);
    hData_DiffXsec->SetBinError(ipt,Data_Xsec_TotalErr[ipt]/BinWidth[ipt]);
    hData_NormDiffXsec->SetBinContent(ipt,hData_Xsec->GetBinContent(ipt)/BinWidth[ipt]/Data_TotalXsec);
    hData_NormDiffXsec->SetBinError(ipt,Data_NormDiffXsec_TotalErr[ipt]);
    cout << scientific << "Data NormDiff Xsec : " << hData_NormDiffXsec->GetBinContent(ipt) << endl;
  }
  
  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataErrBand = new TH1D("hRatioDataErrBand","hRatioDataErrBand",nBins-1,WptLogBins);
  TH1D *hRatioDataStatErr = new TH1D("hRatioDataStatErr","hRatioDataStatErr",nBins-1,WptLogBins);
  TH1D *hRatioDataStatSystErr = new TH1D("hRatioDataStatSystErr","hRatioDataStatSystErr",nBins-1,WptLogBins);
  TH1D *hRatioDataTotalErr = new TH1D("hRatioDataTotalErr","hRatioDataTotalErr",nBins-1,WptLogBins);
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioDataErrBand->SetBinContent(ipt,1.);
    hRatioDataErrBand->SetBinError(ipt,hData_NormDiffXsec->GetBinError(ipt)/hData_NormDiffXsec->GetBinContent(ipt));
    
    hRatioDataStatErr->SetBinContent(ipt,1.);
    hRatioDataStatErr->SetBinError(ipt,Data_NormDiffXsec_StatErr[ipt]/hData_NormDiffXsec->GetBinContent(ipt));
    
    hRatioDataStatSystErr->SetBinContent(ipt,1.);
    hRatioDataStatSystErr->SetBinError(ipt,(Data_NormDiffXsec_StatErr[ipt]+Data_NormDiffXsec_SystErr[ipt])/hData_NormDiffXsec->GetBinContent(ipt));
    
    // New style
    hRatioDataTotalErr->SetBinContent(ipt,1.);
    hRatioDataTotalErr->SetBinError(ipt,Data_NormDiffXsec_TotalErr[ipt]/hData_NormDiffXsec->GetBinContent(ipt));
  }
  

  hData_NormDiffXsec->SetMarkerStyle(kFullCircle); 
  hData_NormDiffXsec->SetMarkerColor(kBlack); 
  hData_NormDiffXsec->SetMarkerSize(1);
  
  hRatioDataStatErr->SetMarkerStyle(kFullCircle); 
  hRatioDataStatErr->SetMarkerColor(kBlack); 
  hRatioDataStatErr->SetMarkerSize(0.6);
  
  /// TGraph
  TGraphErrors *Data_Xsec_Born = new TGraphErrors(hData_NormDiffXsec);
  TGraphErrors *RatioDataStatErrBand = new TGraphErrors(hRatioDataStatErr);
  TGraphErrors *RatioDataStatSystErrBand = new TGraphErrors(hRatioDataStatSystErr);
  TGraphErrors* DataRatio = new TGraphErrors(hRatioDataTotalErr);

  /////==========Real Data related stuff Finished here=========================================



  /////==========Powheg Starts here==============================================
  //
  // Born Level Powheg Yield and Errors, and Converting to X sec yield and errors, making TGraph .
  
  TH1D* hPowheg_Xsec;
  hPowheg_Xsec   = (TH1D*)f_Data->Get("SVD_Born.Gen")->Clone("hPowheg_Xsec");
  hPowheg_Xsec->Scale(1./18.429);
    
  TH1D* hPowheg_Yield_BornAfterFidCut;
  hPowheg_Yield_BornAfterFidCut = (TH1D*)f_Data->Get("SVD_Born.Gen")->Clone("hPowheg_Yield_BornAfterFidCut");
  if(BaseName== "WpToMuNu")  { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Muon_WpToMuNu_S8); }
  if(BaseName== "WmToMuNu")  { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Muon_WmToMuNu_S8); }
  if(BaseName== "WpToEleNu") { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Ele_WpToEleNu_S8); }
  if(BaseName== "WmToEleNu") { hPowheg_Yield_BornAfterFidCut->Scale(1./LumiWeight_Ele_WmToEleNu_S8); }

  // Calculate Total Xsec Powheg
  double Powheg_TotalXsec = hPowheg_Xsec->Integral();

  // Normalized Differential Xsec
  double Powheg_NormDiffXsec[14] = {0.,};
  for(int ipt(1); ipt<14; ipt++)
  {
    Powheg_NormDiffXsec[ipt] = hPowheg_Xsec->GetBinContent(ipt) / BinWidth[ipt] / Powheg_TotalXsec;
    cout << "Powheg NormDiff Xsec : " << Powheg_NormDiffXsec[ipt] << endl;
  }

  // Normalized Differential Errors
  double Powheg_NormDiff_StatErr[14]; // Toy
  double Powheg_NormDiff_PDFErr[14]; 
  double Powheg_NormDiff_TotalErr[14];
   
  Powheg_NormDiffStatErr(BaseName,Powheg_NormDiff_StatErr);
  Powheg_NormPDFUncer(BaseName,Powheg_NormDiffXsec,Powheg_NormDiff_PDFErr);
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    Powheg_NormDiff_TotalErr[ipt] = sqrt(Powheg_NormDiff_StatErr[ipt]*Powheg_NormDiff_StatErr[ipt] + Powheg_NormDiff_PDFErr[ipt]*Powheg_NormDiff_PDFErr[ipt]);
    cout << "Powheg NormDiff Stat Err : " << Powheg_NormDiff_StatErr[ipt] <<"\t PDF Err : " << Powheg_NormDiff_PDFErr[ipt] << "\t TotalErr : " << Powheg_NormDiff_TotalErr[ipt] << endl;
  }

  // Normalized Differential X-sec and error hist
  TH1D *hPowheg_NormDiffXsec   = new TH1D("hPowheg_NormDiffXsec","",nBins-1,WptLogBins);
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hPowheg_NormDiffXsec->SetBinContent(ipt,Powheg_NormDiffXsec[ipt]);
    hPowheg_NormDiffXsec->SetBinError(ipt,Powheg_NormDiff_TotalErr[ipt]);
  }

  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioPowhegStatErrBand = new TH1D("hRatioPowhegStatErrBand","hRatioPowhegStatErrBand",nBins-1,WptLogBins);
  TH1D *hRatioPowhegStatPDFErrBand = new TH1D("hRatioPowhegStatPDFErrBand","hRatioPowhegStatPDFErrBand",nBins-1,WptLogBins);
  TH1D *hRatioPowhegTotalErrBand = new TH1D("hRatioPowhegTotalErrBand","hRatioPowhegTotalErrBand",nBins-1,WptLogBins);
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioPowhegStatErrBand->SetBinContent(ipt,hPowheg_NormDiffXsec->GetBinContent(ipt)/hData_NormDiffXsec->GetBinContent(ipt));
    hRatioPowhegStatErrBand->SetBinError(ipt,(Powheg_NormDiff_StatErr[ipt]/hData_NormDiffXsec->GetBinContent(ipt)));

    hRatioPowhegStatPDFErrBand->SetBinContent(ipt,hPowheg_NormDiffXsec->GetBinContent(ipt)/hData_NormDiffXsec->GetBinContent(ipt));
    hRatioPowhegStatPDFErrBand->SetBinError(ipt,(Powheg_NormDiff_StatErr[ipt]+Powheg_NormDiff_PDFErr[ipt])/hData_NormDiffXsec->GetBinContent(ipt));

    hRatioPowhegTotalErrBand->SetBinContent(ipt,hPowheg_NormDiffXsec->GetBinContent(ipt)/hData_NormDiffXsec->GetBinContent(ipt)); 
    hRatioPowhegTotalErrBand->SetBinError(ipt,(Powheg_NormDiff_TotalErr[ipt]/hData_NormDiffXsec->GetBinContent(ipt))); // Total : sqrt(stat^2 + pdf^2) 

    //cout <<"PowhegRatio and error : " <<  hRatioPowhegTotalErrBand->GetBinContent(ipt) << " +- " << hRatioPowhegTotalErrBand->GetBinError(ipt) << endl;
  }
  
  /// TGraph
  TGraphErrors *Powheg_Xsec_Born = new TGraphErrors(hPowheg_NormDiffXsec);
  TGraphErrors *RatioPowhegStatErrBand = new TGraphErrors(hRatioPowhegStatErrBand);
  TGraphErrors *RatioPowhegStatPDFErrBand = new TGraphErrors(hRatioPowhegStatPDFErrBand);
  TGraphErrors* PowhegRatio = new TGraphErrors(hRatioPowhegTotalErrBand);

  /////==========Powheg related stuff Finished here=========================================
  
  
  
  /////==========FEWZ Starts here==============================================
  //
  // FEWZ X sec value and errors, making TGraph .
  TH1D* hFEWZ_Xsec;
  TH1D* hFEWZ_Xsec_up;
  TH1D* hFEWZ_Xsec_down;
  hFEWZ_Xsec = (TH1D*)f_Fewz->Get("hxsec")->Clone("hFEWZ_Xsec");
  hFEWZ_Xsec_up = (TH1D*)f_Fewz->Get("hxsec_up")->Clone("hFEWZ_Xsec_up");
  hFEWZ_Xsec_down = (TH1D*)f_Fewz->Get("hxsec_down")->Clone("hFEWZ_Xsec_down");

  // Calculate Total Xsec FEWZ
  double FEWZ_TotalXsec = hFEWZ_Xsec->Integral();
  double FEWZ_TotalXsec_up = hFEWZ_Xsec_up->Integral();
  double FEWZ_TotalXsec_down = hFEWZ_Xsec_down->Integral();
  cout << Form("FEWZ Total Xsec : %.2f, up : %.2f, down : %.2f",FEWZ_TotalXsec,FEWZ_TotalXsec_up,FEWZ_TotalXsec_down) << endl;
 
  // Normalized Differential Xsec
  double FEWZ_NormDiffXsec[14] ={0.,};
  double FEWZ_NormDiffXsec_up[14] ={0.,};
  double FEWZ_NormDiffXsec_down[14] ={0.,};
  for(int ipt(1); ipt<=nBins-1;ipt++)
  {
    FEWZ_NormDiffXsec[ipt] = hFEWZ_Xsec->GetBinContent(ipt) / BinWidth[ipt] / FEWZ_TotalXsec;
    FEWZ_NormDiffXsec_up[ipt] = hFEWZ_Xsec_up->GetBinContent(ipt) / BinWidth[ipt] / FEWZ_TotalXsec_up;
    FEWZ_NormDiffXsec_down[ipt] = hFEWZ_Xsec_up->GetBinContent(ipt) / BinWidth[ipt] / FEWZ_TotalXsec_down;
  }

  //Normalized Differential Errors
  double FEWZXsec[14]={0.,}; // for Toy
  double FEWZ_StatErr[14]={0.,}; // for Toy
  double FEWZ_NormDiff_StatErr[14]={0.,}; // Toy
  double FEWZ_NormDiff_ScaleErr[14]={0.,}; // Max(Abs(Normalized up-nominal),Abs(normalized down- nominal))
  double FEWZ_NormDiff_PDFErr[14]={0.,}; // Normalized eigenvector
  double FEWZ_NormDiff_TotalErr[14]={0.,}; // sqrt(stat^2+scale^2+pdf^2)

  // Stat Error
  for(int i(1); i<14; i++)
  {
    FEWZXsec[i] = hFEWZ_Xsec->GetBinContent(i);
    FEWZ_StatErr[i] = hFEWZ_Xsec->GetBinError(i);
  }
  NormDiffToyErr(FEWZXsec,FEWZ_StatErr,FEWZ_NormDiff_StatErr);
  
  // PDF error
  FEWZ_NormPDFUncer(BaseName, FEWZ_NormDiffXsec,FEWZ_NormDiff_PDFErr);

  // Scale Error
  for(int i(1); i<14; i++)
  {
    FEWZ_NormDiff_ScaleErr[i] = TMath::Max(fabs(FEWZ_NormDiffXsec_up[i]-FEWZ_NormDiffXsec[i]) , fabs(FEWZ_NormDiffXsec_down[i]-FEWZ_NormDiffXsec[i]));
  }

  // Total Error
  for(int i(1); i<14; i++)
  {
    FEWZ_NormDiff_TotalErr[i] = sqrt(FEWZ_NormDiff_StatErr[i]*FEWZ_NormDiff_StatErr[i] + FEWZ_NormDiff_PDFErr[i]*FEWZ_NormDiff_PDFErr[i] + FEWZ_NormDiff_ScaleErr[i]*FEWZ_NormDiff_ScaleErr[i]);
  }

  // Normalized Errors check
  for(int i(1); i<14; i++)
  {
    cout<< Form("FEWZ NormDiff Xsec : %f\t Stat : %f\t PDF : %f\t Scale : %f\t TotalUncer : %f",FEWZ_NormDiffXsec[i],FEWZ_NormDiff_StatErr[i],FEWZ_NormDiff_PDFErr[i],FEWZ_NormDiff_ScaleErr[i],FEWZ_NormDiff_TotalErr[i]) << endl;
  }

  // Normalized Differential X-sec and error in hist
  TH1D *hFEWZ_NormDiffXsec   = new TH1D("hFEWZ_NormDiffXsec","hFEWZ_NormDiffXsec",nBins-1,WptLogBins);
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hFEWZ_NormDiffXsec->SetBinContent(ipt,FEWZ_NormDiffXsec[ipt]);
    hFEWZ_NormDiffXsec->SetBinError(ipt,FEWZ_NormDiff_TotalErr[ipt]);
  } 
 
  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioFEWZStatErrBand = new TH1D("hRatioFEWZStatErrBand","hRatioFEWZStatErrBand",nBins-1,WptLogBins);
  TH1D *hRatioFEWZStatPDFErrBand = new TH1D("hRatioFEWZStatPDFErrBand","hRatioFEWZStatPDFErrBand",nBins-1,WptLogBins);
  TH1D *hRatioFEWZStatPDFScaleErrBand = new TH1D("hRatioFEWZStatPDFScaleErrBand","hRatioFEWZStatPDFScaleErrBand",nBins-1,WptLogBins);
  TH1D *hRatioFEWZTotalErrBand = new TH1D("hRatioFEWZTotalErrBand","hRatioFEWZTotalErrBand",nBins-1,WptLogBins);
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioFEWZStatErrBand->SetBinContent(ipt,hFEWZ_NormDiffXsec->GetBinContent(ipt)/hData_NormDiffXsec->GetBinContent(ipt));
    hRatioFEWZStatErrBand->SetBinError(ipt,FEWZ_NormDiff_StatErr[ipt]/hData_NormDiffXsec->GetBinContent(ipt));

    hRatioFEWZStatPDFErrBand->SetBinContent(ipt,hFEWZ_NormDiffXsec->GetBinContent(ipt)/hData_NormDiffXsec->GetBinContent(ipt));
    hRatioFEWZStatPDFErrBand->SetBinError(ipt,(FEWZ_NormDiff_StatErr[ipt]+FEWZ_NormDiff_PDFErr[ipt])/hData_NormDiffXsec->GetBinContent(ipt));

    hRatioFEWZStatPDFScaleErrBand->SetBinContent(ipt,hFEWZ_NormDiffXsec->GetBinContent(ipt)/hData_NormDiffXsec->GetBinContent(ipt));
    hRatioFEWZStatPDFScaleErrBand->SetBinError(ipt,(FEWZ_NormDiff_StatErr[ipt]+FEWZ_NormDiff_PDFErr[ipt]+FEWZ_NormDiff_ScaleErr[ipt])/hData_NormDiffXsec->GetBinContent(ipt));

    hRatioFEWZTotalErrBand->SetBinContent(ipt,hFEWZ_NormDiffXsec->GetBinContent(ipt)/hData_NormDiffXsec->GetBinContent(ipt));
    hRatioFEWZTotalErrBand->SetBinError(ipt,FEWZ_NormDiff_TotalErr[ipt]/hData_NormDiffXsec->GetBinContent(ipt));
  }
  
  /// TGraph
  TGraphErrors *FEWZ_Xsec = new TGraphErrors(hFEWZ_NormDiffXsec);
  TGraphErrors* RatioFEWZStatErrBand = new TGraphErrors(hRatioFEWZStatErrBand);
  TGraphErrors* RatioFEWZStatPDFErrBand = new TGraphErrors(hRatioFEWZStatPDFErrBand);
  TGraphErrors* RatioFEWZStatPDFScaleErrBand = new TGraphErrors(hRatioFEWZStatPDFScaleErrBand);
  TGraphErrors* FEWZRatio = new TGraphErrors(hRatioFEWZTotalErrBand);

  /////==========FEWZ related stuff Finished here=========================================
  
  
  /////==========ResBos Starts here==============================================
  //
  // ResBos X sec value and errors, making TGraph .
 
  TH1D* hResBos30_Xsec;
  hResBos30_Xsec = (TH1D*)f_Resbos->Get("hResbos30")->Clone("hResBos30_Xsec");
  double ResBos_TotalXsec = hResBos30_Xsec->Integral(); //Resbos Total Xsec center grid

  // Calculate NormDiff Xsec Resbos
  double ResBos_NormDiffXsec[14];
  for(int i(1); i<14; i++)
  {
    ResBos_NormDiffXsec[i] = hResBos30_Xsec->GetBinContent(i) / BinWidth[i] / ResBos_TotalXsec;
    cout << "ResBos NormDiffXsec : " << ResBos_NormDiffXsec[i] << endl;
  }

  // Calculate NormDiff Xsec Resbos other grids
  TH1D* AllResbos[7];
  double AllResbos_TotalXsec[7];
  double AllResbos_NormDiffXsec[7][14];
  for( int i(0);i<7;i++)
  {
    Numb = 29+i;
    sprintf(tmpName_org,"hResbos%d",Numb);
    sprintf(tmpName,"AllResbos_%d",i);
    AllResbos[i] = (TH1D*)f_Resbos->Get(tmpName_org)->Clone(tmpName);
    AllResbos_TotalXsec[i] = AllResbos[i]->Integral(); // Total Xsec for Resbos grid
    cout << "AllResbos TotalXsec : " << AllResbos_TotalXsec[i] << endl;
    for(int j(1); j<14; j++)
    {
      AllResbos_NormDiffXsec[i][j]=AllResbos[i]->GetBinContent(j) / BinWidth[j]/ AllResbos_TotalXsec[i];
      cout << "Resbos NormDiff : " << AllResbos[i]->GetBinContent(j) / BinWidth[j]/ AllResbos_TotalXsec[i]  << endl;
    }
  }


  //Errors
  Double_t Resb_errMax[nBins-1];
  Double_t Resb_errMin[nBins-1];
  Double_t NormResb_errMax[nBins-1];
  Double_t NormResb_errMin[nBins-1];
  double tmpVal,tmpDiff;

  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    double nomVal  = ResBos_NormDiffXsec[ipt+1];
    Resb_errMax[ipt] = -99999;
    Resb_errMin[ipt] = 990009;
    for (int i(0);i<7;i++)
    {
      tmpVal  = AllResbos_NormDiffXsec[i][ipt+1];
      tmpDiff = tmpVal - nomVal;
      if( tmpDiff > Resb_errMax[ipt] ) Resb_errMax[ipt] = tmpDiff;
      if( tmpDiff < Resb_errMin[ipt] ) Resb_errMin[ipt] = tmpDiff;
    }
  cout << "Resbos Xsec : " << hResBos30_Xsec->GetBinContent(ipt+1) << "\t error+ : " << Resb_errMax[ipt] << "\t error - : " << Resb_errMin[ipt] << endl; 

    if (Resb_errMax[ipt] < 0) Resb_errMax[ipt] = 0.;
    if (Resb_errMin[ipt] > 0) Resb_errMin[ipt] = 0.;
    if (Resb_errMin[ipt] < 0) Resb_errMin[ipt] = -Resb_errMin[ipt];
    Resb_errMax[ipt] = Resb_errMax[ipt];
    Resb_errMin[ipt] = Resb_errMin[ipt];
  }

  // Normalized Differential X-sec and error 
  TH1D *hResBos30_NormDiffXsec   = new TH1D("hResBos30_NormDiffXsec","hResBos30_NormDiffXsec",nBins-1,WptLogBins);
  TH1D *hResBosRatio   = new TH1D("hResBosRatio","hResBosRatio",nBins-1,WptLogBins);
  Double_t hResb30_NormDiffXsec[nBins-1];
  Double_t hResb30_Error[nBins-1];
  Double_t RatioResbVal[nBins-1],errResbosDataLo[nBins-1],errResbosDataHi[nBins-1];
  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    hResBos30_NormDiffXsec->SetBinContent(ipt+1,ResBos_NormDiffXsec[ipt+1]);

    hResb30_NormDiffXsec[ipt] = hResBos30_NormDiffXsec->GetBinContent(ipt+1);
   
    RatioResbVal[ipt]=hResBos30_NormDiffXsec->GetBinContent(ipt+1)/hData_NormDiffXsec->GetBinContent(ipt+1);
    errResbosDataLo[ipt]=Resb_errMin[ipt]/hData_NormDiffXsec->GetBinContent(ipt+1);
    errResbosDataHi[ipt]=Resb_errMax[ipt]/hData_NormDiffXsec->GetBinContent(ipt+1);
    
    hResb30_Error[ipt] = (errResbosDataLo[ipt] > errResbosDataHi[ipt]) ? errResbosDataLo[ipt] : errResbosDataHi[ipt];
    hResBosRatio->SetBinContent(ipt+1,RatioResbVal[ipt]);
    hResBosRatio->SetBinError(ipt+1,hResb30_Error[ipt]);

    cout << "Resbos NormDiffXsec : " << hResb30_NormDiffXsec[ipt] << "\t errorHi : " << Resb_errMax[ipt] / hResb30_NormDiffXsec[ipt] *100 << "\t errorLo : " << Resb_errMin[ipt] / hResb30_NormDiffXsec[ipt] *100  << endl;
    //cout << "resbos ratio and error : " << RatioResbVal[ipt] << " +- " << hResb30_Error[ipt] << endl; 
  }
  ///   Theory/Data ratio plot errors related to Real Data
  TGraphAsymmErrors *Resb30_NormDiffXsec = new TGraphAsymmErrors(nBins-1, ax, hResb30_NormDiffXsec, aex, aex, Resb_errMin, Resb_errMax);
  TGraphAsymmErrors *RatioResbosErrBand = new TGraphAsymmErrors(nBins-1, ax, RatioResbVal, aex, aex, errResbosDataLo, errResbosDataHi);
  TGraphErrors* ResBosRatio = new TGraphErrors(hResBosRatio);
  

  /////==========ResBos related stuff Finished here=========================================
 
  
  //// Now design and Draw 
  gStyle->SetLineWidth(2.);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesSpacing(0.75);
  gStyle->SetHatchesLineWidth(2);

  //Color Transparent
  TColor *colRed = gROOT->GetColor(kRed);
  TColor *colBlue = gROOT->GetColor(kBlue);
  TColor *colGreen = gROOT->GetColor(kGreen);
  colRed->SetAlpha(0.1);
  colBlue->SetAlpha(0.2);
  colGreen->SetAlpha(0.2);

  TPaveText *cmspre = new TPaveText(0.6,0.90,0.95,0.95,"NDC");
  cmspre->SetBorderSize(0);
  cmspre->SetTextSize(0.04);
  cmspre->SetFillStyle(0);
  cmspre->AddText("CMS Preliminary");

  TLegend *lL =new TLegend(0.2,0.17,0.52,0.40); lL->SetFillColor(0); lL->SetBorderSize(0);
  lL->AddEntry(Data_Xsec_Born,"data","PL");
  lL->AddEntry(Powheg_Xsec_Born,"POWHEG CT10 NLO","f");
  lL->AddEntry(FEWZ_Xsec,"FEWZ CT10 NNLO","f");
  lL->AddEntry(Resb30_NormDiffXsec,"ResBos CT10 NNLL","f");

  TPaveText *tb = new TPaveText(0.2,0.42,0.52,0.58,"NDC");
  tb->SetBorderSize(0);
  tb->SetFillStyle(0);
  tb->AddText("L = 18.4 pb^{-1}, #sqrt{s} = 8 TeV");
  if (BaseName=="WpToMuNu")
    tb->AddText("W^{+} #rightarrow #mu^{+} #nu");
  if (BaseName=="WmToMuNu")
    tb->AddText("W^{-} #rightarrow #mu^{-} #bar{#nu}");
  if (BaseName=="WpToEleNu")
    tb->AddText("W^{+} #rightarrow e^{+} #nu");
  if (BaseName=="WmToEleNu")
    tb->AddText("W^{-} #rightarrow e^{-} #bar{#nu}");

  // Canvas for distribution
  TCanvas *lC1 = new TCanvas("Can","Can",800,840); lC1->cd(); lC1->SetLogy();
  lC1->cd(1)->SetPad(0.05,0.05,0.95,0.95);
  lC1->cd(1)->SetTopMargin(0.1);
  lC1->cd(1)->SetBottomMargin(0.1);
  lC1->cd(1)->SetLeftMargin(0.17);
  lC1->cd(1)->SetRightMargin(0.07);
  lC1->cd(1)->SetTickx(1);
  lC1->cd(1)->SetTicky(1);
  gPad->SetLogx(1);
  gPad->SetLogy(1);

  // Frame setting 
  Powheg_Xsec_Born->GetYaxis()->SetRangeUser(1e-7,5*hResb30_NormDiffXsec[1]);
  Powheg_Xsec_Born->SetTitle("");
  Powheg_Xsec_Born->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T}^{W} [pb/GeV]");
  Powheg_Xsec_Born->GetYaxis()->SetTitleOffset(1.45);
  Powheg_Xsec_Born->GetYaxis()->SetTitleSize(0.05);

  Powheg_Xsec_Born->GetXaxis()->SetLabelSize(0.03);
  Powheg_Xsec_Born->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  Powheg_Xsec_Born->GetXaxis()->SetTitleSize(0.04);
  Powheg_Xsec_Born->GetXaxis()->SetTitleOffset(0.55);

  //Data

  //Powheg
  Powheg_Xsec_Born->SetFillColor(kRed);
  FEWZ_Xsec->SetFillColor(kGreen);
  Resb30_NormDiffXsec->SetFillColor(kBlue);

  //Draw Original Diff-Xsec Distribution
  Powheg_Xsec_Born->Draw("A2");
  FEWZ_Xsec->Draw("2");
  Resb30_NormDiffXsec->Draw("2");
  Data_Xsec_Born->Draw("P");
  
  lL->Draw();
  tb->Draw();
  cmspre->Draw();

  if(BaseName=="WpToMuNu")
    sprintf(tmpName,"wpmnCrossSec.png");
  if(BaseName=="WmToMuNu")
    sprintf(tmpName,"wmmnCrossSec.png");
  if(BaseName=="WpToEleNu")
    sprintf(tmpName,"wpenCrossSec.png");
  if(BaseName=="WmToEleNu")
    sprintf(tmpName,"wmenCrossSec.png");
  lC1->SaveAs(tmpName);

  // Ratio plot style 
  RatioDataStatErrBand->SetMarkerStyle(20);
  RatioDataStatErrBand->SetMarkerColor(kBlack);
  RatioDataStatErrBand->SetMarkerSize(0.7);
  RatioDataStatErrBand->SetLineWidth(2.0);
  RatioDataStatErrBand->SetLineColor(kBlack);
   
  DataRatio->SetFillStyle(3354);
  DataRatio->SetFillColor(kGray+1);

// Resbos Ratio plot style
  RatioResbosErrBand->SetMarkerStyle(20);
  RatioResbosErrBand->SetMarkerColor(kBlue+2);

  RatioResbosErrBand->SetFillColor(kBlue-9);
  //RatioResbosErrBand->SetFillStyle(3001);

// Powheg Ratio plot style
  RatioPowhegStatErrBand->SetMarkerStyle(22);
  RatioPowhegStatErrBand->SetMarkerColor(kRed+2);

  RatioPowhegStatErrBand->SetFillColor(kRed-10);
  //RatioPowhegStatErrBand->SetFillStyle(3001);

  RatioPowhegStatPDFErrBand->SetFillColor(kRed-2);
  //RatioPowhegStatPDFErrBand->SetFillStyle(3001);

// FEWZ Ratio plot style
  RatioFEWZStatErrBand->SetFillColor(kGreen-10);
  //RatioFEWZStatErrBand->SetFillStyle(3001);
  
  RatioFEWZStatPDFErrBand->SetFillColor(kGreen-6);
  //RatioFEWZStatPDFErrBand->SetFillStyle(3001);

  RatioFEWZStatPDFScaleErrBand->SetFillColor(kCyan-7);
  //RatioFEWZStatPDFScaleErrBand->SetFillStyle(3001);

  RatioFEWZStatErrBand->SetMarkerStyle(21);
  RatioFEWZStatErrBand->SetMarkerColor(kGreen+2);

  
  // Canvas for Theory/data Ratio
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
  rL1->AddEntry(RatioResbosErrBand,"Theory syst","F");
  rL1->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  hRatioDataStatErr->SetTitle("");
  rL1->AddEntry(hRatioDataStatErr,"","");
  rL1->AddEntry(DataRatio,"Data stat+syst","F");
  rL1->SetTextSize(0.07);

  TLegend *tL1 =new TLegend(0.17,0.72,0.37,0.82); tL1->SetFillColor(0); tL1->SetBorderSize(0);
  tL1->AddEntry(RatioResbosErrBand,"ResBos","F");
  tL1->SetTextSize(0.12);
  tL1->SetTextFont(2);

  TPaveText *cmspre2 = new TPaveText(0.6,0.90,0.97,0.95,"NDC");
  cmspre2->SetBorderSize(0);
  cmspre2->SetFillStyle(0);
  cmspre2->SetTextSize(0.13);
  cmspre2->AddText("CMS Preliminary");

  TPaveText *lumi = new TPaveText(0.10,0.57,0.60,0.67,"NDC");
  lumi->SetBorderSize(0);
  lumi->SetFillStyle(0);
  lumi->SetTextSize(0.10);
  lumi->AddText("L = 18.4 pb^{-1}, #sqrt{s} = 8 TeV");

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

  TH1D* hRatioResbosDummy = new TH1D("","",nBins-1,WptLogBins);
  hRatioResbosDummy->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatioResbosDummy->GetYaxis()->SetTitle("Theory/Data");
  hRatioResbosDummy->GetYaxis()->CenterTitle();
  hRatioResbosDummy->GetYaxis()->SetTitleOffset(0.4);
  hRatioResbosDummy->GetYaxis()->SetTitleSize(0.12);
  hRatioResbosDummy->GetYaxis()->SetLabelSize(0.10);
  hRatioResbosDummy->GetYaxis()->SetNdivisions(405);
  hRatioResbosDummy->GetXaxis()->SetTitleOffset(0.6);
  hRatioResbosDummy->GetXaxis()->SetTitleSize(0.08);
  hRatioResbosDummy->GetXaxis()->SetLabelSize(0);
  hRatioResbosDummy->Draw();
  RatioResbosErrBand->Draw("2 P");
  DataRatio->Draw("2");
  RatioDataStatErrBand->Draw("P E");
  rL1->Draw();
  tL1->Draw();
  tb4->Draw();
  cmspre2->Draw();
  lumi->Draw();

  //Powheg Ratio plot
  lC2->cd(2)->SetPad(0,0.39,0.95,0.65);
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
  rL2->AddEntry(RatioPowhegStatErrBand,"stat","F");
  rL2->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  rL2->AddEntry(RatioPowhegStatPDFErrBand,"PDF    ","F");
  rL2->AddEntry(DataRatio,"Data stat+syst","F");
  rL2->SetTextSize(0.07);

  TLegend *tL2 =new TLegend(0.17,0.85,0.37,0.95); tL2->SetFillColor(0); tL2->SetBorderSize(0);
  tL2->AddEntry(RatioPowhegStatErrBand,"POWHEG","F");
  tL2->SetTextSize(0.12);
  tL2->SetTextFont(2);
  
  TH1D* hRatioPowhegDummy = new TH1D("","",nBins-1,WptLogBins);
  hRatioPowhegDummy->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatioPowhegDummy->GetYaxis()->SetTitle("Theory/Data");
  hRatioPowhegDummy->GetYaxis()->CenterTitle();
  hRatioPowhegDummy->GetYaxis()->SetTitleOffset(0.4);
  hRatioPowhegDummy->GetYaxis()->SetTitleSize(0.12);
  hRatioPowhegDummy->GetYaxis()->SetLabelSize(0.10);
  hRatioPowhegDummy->GetYaxis()->SetNdivisions(405);
  hRatioPowhegDummy->GetXaxis()->SetTitleOffset(0.6);
  hRatioPowhegDummy->GetXaxis()->SetTitleSize(0.08);
  hRatioPowhegDummy->GetXaxis()->SetLabelSize(0);
  hRatioPowhegDummy->Draw();
  RatioPowhegStatPDFErrBand->Draw("2");
  RatioPowhegStatErrBand->Draw("2 P");
  DataRatio->Draw("2");
  RatioDataStatErrBand->Draw("P E");
  rL2->Draw();
  tL2->Draw();
  

  // FEWZ Ratio Plot
  lC2->cd(3)->SetPad(0,0.07,0.95,0.37);
  lC2->cd(3)->SetTopMargin(0.01);
  lC2->cd(3)->SetBottomMargin(0.15);
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
  TLegend *rL3 =new TLegend(0.2,0.20,0.58,0.40); rL3->SetFillColor(0); rL3->SetBorderSize(0);
  rL3-> SetNColumns(2);
  rL3->AddEntry(RatioFEWZStatErrBand,"stat","F");
  rL3->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  rL3->AddEntry(RatioFEWZStatPDFErrBand,"PDF","F");
  rL3->AddEntry(DataRatio,"Data stat+syst","F");
  rL3->AddEntry(RatioFEWZStatPDFScaleErrBand,"QCD scales","F");
  rL3->SetTextSize(0.07);

  TLegend *tL3 =new TLegend(0.17,0.85,0.37,0.95); tL3->SetFillColor(0); tL3->SetBorderSize(0);
  tL3->AddEntry(RatioFEWZStatErrBand,"FEWZ","F");
  tL3->SetTextSize(0.12);
  tL3->SetTextFont(2);

  TH1D* hRatioFEWZDummy = new TH1D("","",nBins-1,WptLogBins);
  hRatioFEWZDummy->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatioFEWZDummy->GetYaxis()->SetTitle("Theory/Data");
  hRatioFEWZDummy->GetYaxis()->CenterTitle();
  hRatioFEWZDummy->GetYaxis()->SetTitleOffset(0.4);
  hRatioFEWZDummy->GetYaxis()->SetTitleSize(0.12);
  hRatioFEWZDummy->GetYaxis()->SetLabelSize(0.10);
  hRatioFEWZDummy->GetYaxis()->SetNdivisions(405);
  hRatioFEWZDummy->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  hRatioFEWZDummy->GetXaxis()->SetTitleOffset(0.6);
  hRatioFEWZDummy->GetXaxis()->SetTitleSize(0.08);
  hRatioFEWZDummy->GetXaxis()->SetLabelSize(0.08);
  hRatioFEWZDummy->Draw();
  RatioFEWZStatPDFScaleErrBand->Draw("2");
  RatioFEWZStatPDFErrBand->Draw("2");
  RatioFEWZStatErrBand->Draw("2 P");
  DataRatio->Draw("2");
  RatioDataStatErrBand->Draw("P E");
  rL3->Draw();
  tL3->Draw();

  if(BaseName=="WpToMuNu")
    sprintf(tmpName,"wpmnRatioTheoData.png");
  if(BaseName=="WmToMuNu")
    sprintf(tmpName,"wmmnRatioTheoData.png");
  if(BaseName=="WpToEleNu")
    sprintf(tmpName,"wpenRatioTheoData.png");
  if(BaseName=="WmToEleNu")
    sprintf(tmpName,"wmenRatioTheoData.png");
  lC2->SaveAs(tmpName);
 
 /* 
  // Start New Style 
  TCanvas *Cnew = new TCanvas("NewCanvas","NewCanvas",800,900);
  //Cnew->Divide(1,3,0,0);
  Cnew->Divide(1,2,0,0);
  Cnew->cd(1)->SetPad(0,0.50,0.95,0.95);
  Cnew->cd(1)->SetTopMargin(0.05);
  Cnew->cd(1)->SetLeftMargin(0.15);
  Cnew->cd(1)->SetRightMargin(0.07);
  Cnew->cd(1)->SetTickx(1);
  Cnew->cd(1)->SetTicky(1);
  Cnew->cd(1)->SetLogx();
  Cnew->cd(1)->SetLogy();
  gStyle->SetLineWidth(2.);

  //Draw Original Diff-Xsec Distribution
  Powheg_Xsec_Born->Draw("A2");
  FEWZ_Xsec->Draw("2");
  Resb30_NormDiffXsec->Draw("2");
  Data_Xsec_Born->Draw("P");

  // Draw new style ratio pad2 (Theory / Data)
  Cnew->cd(2)->SetPad(0,0.05,0.95,0.45);
  Cnew->cd(2)->SetBottomMargin(0.1);
  Cnew->cd(2)->SetLeftMargin(0.15);
  Cnew->cd(2)->SetRightMargin(0.07);
  Cnew->cd(2)->SetTickx(1);
  Cnew->cd(2)->SetTicky(1);
  Cnew->cd(2)->SetLogx();
  gStyle->SetLineWidth(2.);

  // set canvas
  TH1D *hRatioDummy2 = new TH1D("hRatioDummy2","",nBins-1,WptLogBins);
  hRatioDummy2->GetYaxis()->SetRangeUser(0.6,1.5);
  hRatioDummy2->GetYaxis()->SetTitle("Theory / Data");
  hRatioDummy2->GetYaxis()->CenterTitle();
  hRatioDummy2->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  hRatioDummy2->GetYaxis()->SetTitleSize(0.08);
  hRatioDummy2->GetYaxis()->SetLabelSize(0.05);
  hRatioDummy2->GetYaxis()->SetTitleOffset(0.56);
  hRatioDummy2->GetYaxis()->SetNdivisions(405);
  hRatioDummy2->GetXaxis()->SetTitleSize(0.05);
  hRatioDummy2->GetXaxis()->SetLabelSize(0.07);

 
  DataRatio->SetFillColor(kBlack);
  //DataRatio->SetFillStyle(1001);
  DataRatio->SetFillStyle(3001);

  PowhegRatio->SetFillColor(kRed);
  //PowhegRatio->SetFillStyle(3002);
  PowhegRatio->SetFillStyle(3004);
  PowhegRatio->SetMarkerStyle(22); // triangle style
  PowhegRatio->SetMarkerColor(kRed+2);
  PowhegRatio->SetMarkerSize(1.0);
  PowhegRatio->SetLineColor(kRed);

  FEWZRatio->SetFillColor(kGreen);
  //FEWZRatio->SetFillStyle(3002);
  FEWZRatio->SetFillStyle(3005);
  FEWZRatio->SetMarkerStyle(21); // square
  FEWZRatio->SetMarkerColor(kGreen+2);
  FEWZRatio->SetMarkerSize(1.0);
  FEWZRatio->SetLineColor(kGreen);

  ResBosRatio->SetFillColor(kBlue);
  ResBosRatio->SetFillStyle(3013);
  ResBosRatio->SetMarkerStyle(34); // cross
  ResBosRatio->SetMarkerColor(kBlue+2);
  ResBosRatio->SetMarkerSize(1.0);
  ResBosRatio->SetLineColor(kBlue);

  // Draw line
  TLine* line = new TLine(0,1,600,1);
  line -> SetLineColor(kBlack);
  line -> SetLineWidth(1);
  line->SetLineStyle(1);
  
  // Draw
  hRatioDummy2->Draw();
  DataRatio->Draw("2");
  line->Draw("L");
  PowhegRatio->Draw("5 P");
  ResBosRatio->Draw("5 P");
  FEWZRatio->Draw("5 P");

  //gPad->RedrawAxis();

  Cnew->SaveAs(BaseName+"_NewStyle.png");
  */

  return 0;

}

void Powheg_NormPDFUncer(const TString BaseName, double *Xsec, double *Err)
{
  // Powheg Norm PDF error in % unit
  double PDFErr[14] = {0.,};

  if (BaseName=="WpToMuNu")
  {
    PDFErr[1] = 0.547639;  
    PDFErr[2] = 0.11084 ;
    PDFErr[3] = 0.423104;
    PDFErr[4] = 0.485296;
    PDFErr[5] = 0.593284;
    PDFErr[6] = 0.523514;
    PDFErr[7] = 0.456424;
    PDFErr[8] = 0.390481;
    PDFErr[9] = 0.423705;
    PDFErr[10] =1.07401 ;
    PDFErr[11] =1.92522 ;
    PDFErr[12] =2.57608 ;
    PDFErr[13] =3.85306 ;
  }
  if (BaseName=="WmToMuNu")
  {
    PDFErr[1] = 0.805695;  
    PDFErr[2] = 0.193853;
    PDFErr[3] = 0.529425;
    PDFErr[4] = 0.667961;
    PDFErr[5] = 0.773128;
    PDFErr[6] = 0.69784 ;
    PDFErr[7] = 0.707969;
    PDFErr[8] = 0.529504;
    PDFErr[9] = 0.780084;
    PDFErr[10] =1.43238 ;
    PDFErr[11] =2.21681 ;
    PDFErr[12] =3.18584 ;
    PDFErr[13] =3.50682 ;
  }
  if (BaseName=="WpToEleNu")
  {
    PDFErr[1] = 0.561873;  
    PDFErr[2] = 0.103198;
    PDFErr[3] = 0.396217;
    PDFErr[4] = 0.520109;
    PDFErr[5] = 0.586403;
    PDFErr[6] = 0.563382;
    PDFErr[7] = 0.456273;
    PDFErr[8] = 0.36933 ;
    PDFErr[9] = 0.444609;
    PDFErr[10] =1.07102 ;
    PDFErr[11] =1.46382 ;
    PDFErr[12] =2.20702 ;
    PDFErr[13] =2.55988 ;
  }
  if (BaseName=="WmToEleNu")
  {
    PDFErr[1] = 0.746642;  
    PDFErr[2] = 0.203343;
    PDFErr[3] = 0.514733;
    PDFErr[4] = 0.645094;
    PDFErr[5] = 0.721366;
    PDFErr[6] = 0.575011;
    PDFErr[7] = 0.526136;
    PDFErr[8] = 0.492826;
    PDFErr[9] = 0.717755;
    PDFErr[10] =1.33446 ;
    PDFErr[11] =2.2485  ;
    PDFErr[12] =2.47387 ;
    PDFErr[13] =3.70401 ;
  }

  for(int i(1); i<14; i++)
  {
    Err[i] = Xsec[i] * PDFErr[i] * 0.01;
    cout << "PDFErr in function : " << Err[i] << endl;
  }
}

void FEWZ_NormPDFUncer(const TString BaseName, double *Xsec, double *Err)
{
  // FEWZ Norm PDF error in % unit
  double PDFErr[14] = {0.,};

  if (BaseName=="WpToMuNu")
  {
    PDFErr[1] = 2.3193;  
    PDFErr[2] = 0.9559;
    PDFErr[3] = 1.2003;
    PDFErr[4] = 1.4857;
    PDFErr[5] = 1.8138;
    PDFErr[6] = 2.1539;
    PDFErr[7] = 2.6047;
    PDFErr[8] = 3.1196;
    PDFErr[9] = 4.0807;
    PDFErr[10] = 5.1796;
    PDFErr[11] =6.1866;
    PDFErr[12] =8.4432;
    PDFErr[13] =9.9301;
  }
  if (BaseName=="WmToMuNu")
  {
    PDFErr[1] = 3.2060 ;  
    PDFErr[2] = 1.1572 ;
    PDFErr[3] = 1.4176 ;
    PDFErr[4] = 1.8612 ;
    PDFErr[5] = 2.2587 ;
    PDFErr[6] = 2.6728 ;
    PDFErr[7] = 3.2931 ;
    PDFErr[8] = 3.8921 ;
    PDFErr[9] = 5.0757 ;
    PDFErr[10] =6.4401 ;
    PDFErr[11] =8.0128 ;
    PDFErr[12] =9.7073 ;
    PDFErr[13] =11.4727;
  }
  if (BaseName=="WpToEleNu")
  {
    PDFErr[1] = 2.1399;  
    PDFErr[2] = 0.9745;
    PDFErr[3] = 1.1343;
    PDFErr[4] = 1.3962;
    PDFErr[5] = 1.7307;
    PDFErr[6] = 2.0348;
    PDFErr[7] = 2.4141;
    PDFErr[8] = 2.9025;
    PDFErr[9] = 3.8044;
    PDFErr[10] =5.0325;
    PDFErr[11] =6.2783;
    PDFErr[12] =7.4780;
    PDFErr[13] =9.4248;
  }
  if (BaseName=="WmToEleNu")
  {
    PDFErr[1] = 2.9607 ;  
    PDFErr[2] = 1.0371 ;
    PDFErr[3] = 1.4474 ;
    PDFErr[4] = 1.7560 ;
    PDFErr[5] = 2.0934 ;
    PDFErr[6] = 2.5630 ;
    PDFErr[7] = 3.0683 ;
    PDFErr[8] = 3.6932 ;
    PDFErr[9] = 4.5908 ;
    PDFErr[10] =6.1002 ;
    PDFErr[11] =7.6093 ;
    PDFErr[12] =8.9306 ;
    PDFErr[13] =11.1906;
  }

  for(int i(1); i<14; i++)
  {
    Err[i] = Xsec[i] * PDFErr[i] * 0.01;
    cout << "PDFErr in function : " << Err[i] << endl;
  }
}

void Powheg_NormDiffStatErr(const TString BaseName, double *Err)
{
  // Powheg Norm Stat error
  double StatErr[14] = {0.,};

  if (BaseName=="WpToMuNu")
  {
    StatErr[1] = 0.0000903743 ;  
    StatErr[2] = 0.0001154903 ;
    StatErr[3] = 0.0000936253 ;
    StatErr[4] = 0.0000658495 ;
    StatErr[5] = 0.0000553982 ;
    StatErr[6] = 0.0000337565 ;
    StatErr[7] = 0.0000256861 ;
    StatErr[8] = 0.0000128208 ;
    StatErr[9] = 0.0000049989 ;
    StatErr[10] =0.0000024574 ;
    StatErr[11] =0.0000013903 ;
    StatErr[12] =0.0000006696 ;
    StatErr[13] =0.0000000808 ;
  }
  if (BaseName=="WmToMuNu")
  {
    StatErr[1] = 0.0000631993 ;  
    StatErr[2] = 0.0000813040 ;
    StatErr[3] = 0.0000664781 ;
    StatErr[4] = 0.0000466908 ;
    StatErr[5] = 0.0000400361 ;
    StatErr[6] = 0.0000244169 ;
    StatErr[7] = 0.0000186636 ;
    StatErr[8] = 0.0000092752 ;
    StatErr[9] = 0.0000036048 ;
    StatErr[10] =0.0000017660 ;
    StatErr[11] =0.0000010185 ;
    StatErr[12] =0.0000004666 ;
    StatErr[13] =0.0000000537 ;
  }
  if (BaseName=="WpToEleNu")
  {
    StatErr[1] = 0.0000563160 ;  
    StatErr[2] = 0.0000718086 ;
    StatErr[3] = 0.0000581501 ;
    StatErr[4] = 0.0000407893 ;
    StatErr[5] = 0.0000344829 ;
    StatErr[6] = 0.0000209113 ;
    StatErr[7] = 0.0000157895 ;
    StatErr[8] = 0.0000078463 ;
    StatErr[9] = 0.0000031049 ;
    StatErr[10] =0.0000015731 ;
    StatErr[11] =0.0000008835 ;
    StatErr[12] =0.0000004167 ;
    StatErr[13] =0.0000000501 ;
  }
  if (BaseName=="WmToEleNu")
  {
    StatErr[1] = 0.0000602392;  
    StatErr[2] = 0.0000772968;
    StatErr[3] = 0.0000631681;
    StatErr[4] = 0.0000444215;
    StatErr[5] = 0.0000378046;
    StatErr[6] = 0.0000230441;
    StatErr[7] = 0.0000174330;
    StatErr[8] = 0.0000087733;
    StatErr[9] = 0.0000034816;
    StatErr[10] =0.0000017404;
    StatErr[11] =0.0000009611;
    StatErr[12] =0.0000004441;
    StatErr[13] =0.0000000489;
  }

  for(int i(1); i<14; i++)
  {
    Err[i] = StatErr[i];
    cout << "StatErr in function : " << Err[i] << endl;
  }
}

void NormDiffToyErr(double *_Xsec, double *_Err, double *_NormDiff_Err)
{
  int Ntoy = 100000;

  double mx2[14] = {0.,};
  double m2x[14] = {0.,};
  double temp[14] = {0.,};

  // Check _Xsec and _StatErr argument is correct 
//  for(int j=0; j<12; ++j)
//  {
    //cout << Form("_Xsec : %.4f \t _StatErr : %.4f",_Xsec[j],_StatErr[j]) << endl;
//  }

  // Here Calculate Normalized Stat Error by using toy
  for(int e=0; e<Ntoy; ++e)
  {
    for(int i=1; i<14; ++i)
    {   
        temp[i]=gRandom->Gaus(_Xsec[i], _Err[i]);
    }   
    double sum=0;
    for(int i=1; i<14; ++i) sum+=temp[i];
    for(int i=1; i<14; ++i)
    {   
      temp[i]/=sum;
      m2x[i]+=temp[i];
      mx2[i]+=temp[i]*temp[i];
    }   
  }

  double rms[14];
  for(int i=1; i<14; ++i)
  {
    m2x[i]/=Ntoy;
    mx2[i]/=Ntoy;
    rms[i] = sqrt(mx2[i]-m2x[i]*m2x[i]);
  }
 
  // return value
  for(int i=1; i<14; ++i)
  {
    _NormDiff_Err[i] = rms[i]/ BinWidth[i]; // return rms to Norm_Err
  }
  for(int i=1; i<14; ++i)
  {
    cout << Form("rms/BinWidth = %.8f",rms[i]/BinWidth[i]) << endl;
  }
}

