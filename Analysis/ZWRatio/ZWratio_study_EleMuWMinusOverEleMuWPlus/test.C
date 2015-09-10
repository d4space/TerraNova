#include <iostream>

void test(TString BaseName){

  // READ root files
  TFile *f_in;
  if (BaseName=="WpToMuNu")
  {
    f_in = new TFile("Wpt_WpToMuNuNormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  }
  if (BaseName=="WmToMuNu")
  {
    f_in = new TFile("Wpt_WmToMuNuNormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  }
  if (BaseName=="WpToEleNu")
  {
    f_in = new TFile("Wpt_WpToEleNuNormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  }
  if (BaseName=="WmToEleNu")
  {
    f_in = new TFile("Wpt_WmToEleNuNormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  }

  // Get NormDiff xsec histogram
  TH1D* h_orig = (TH1D*)f_in->Get("hData_Xsec_BornLogScaleNorm12");
  
  for(int i(1); i<13; i++) 
  {
    cout <<BaseName << " NormDiffXsec : " << h_orig->GetBinContent(i) << endl;
  }


  // Common Fiducial SF root file for Ele
  if (BaseName == "WpToEleNu" || BaseName == "WmToEleNu")
  {
    TFile *f_RatioSF = new TFile("CommonFid_EleRatio_fromPowhegMC.root");

    TH1D* h_RatioSF; 
    if (BaseName=="WpToEleNu")
    {
      h_RatioSF = (TH1D*)f_RatioSF->Get("h1_ElePRatioPreFSR_norm");
    }
    if (BaseName=="WmToEleNu")
    {
      h_RatioSF = (TH1D*)f_RatioSF->Get("h1_EleMRatioPreFSR_norm");
    }

    for(int i(1); i<13; i++) 
    {
      cout << "h_RatioSF : " << h_RatioSF->GetBinContent(i) << endl;
    }

    TH1D* h_AfterSF = new TH1D("h_AfterSF","h_AfterSF",12,0,12);
    for(int i(1); i<13; i++) 
    {
      h_AfterSF->SetBinContent(i,h_orig->GetBinContent(i)*h_RatioSF->GetBinContent(i));
      cout <<BaseName << " NormDiffXsec after SF applied : " << h_AfterSF->GetBinContent(i) << endl;
    }
  }

  // Save histogram to same format
  TH1D* h_NormDiffXsec = new TH1D("h_NormDiffXsec_CommonFid","h_NormDiffXsec_CommonFid",12,0,12);
  for(int i(1); i<13; i++) 
  {
    h_NormDiffXsec->SetBinContent(i,h_orig->GetBinContent(i));
    if (BaseName == "WpToEleNu" || BaseName == "WmToEleNu") h_NormDiffXsec->SetBinContent(i,h_AfterSF->GetBinContent(i));
    cout << "h_NormDiffXsec : " << h_NormDiffXsec->GetBinContent(i) << endl;
  }

  // Total Uncertainty Calculation
  TH1D* h_Stat = (TH1D*)f_in->Get("hData_StatErr12")->Clone("h_Stat");
  
  // systematic errors
  TH1D* h_Toy = (TH1D*)f_in->Get("hData_EffiToySystErr12")->Clone("h_Toy");
  TH1D* h_IDIsoSig = (TH1D*)f_in->Get("hData_IDIsoSigShapeSystErr12")->Clone("h_IDIsoSig");
  TH1D* h_IDIsoBkg = (TH1D*)f_in->Get("hData_IDIsoBkgrShapeSystErr12")->Clone("h_IDIsoBkg");
  if (BaseName == "WpToMuNu" || BaseName == "WmToMuNu")
  {
    TH1D* h_TrackSig = (TH1D*)f_in->Get("hData_TrackSigShapeSystErr12")->Clone("h_TrackSig");
    TH1D* h_TrackBkg = (TH1D*)f_in->Get("hData_TrackBkgrShapeSystErr12")->Clone("h_TrackBkg");
    TH1D* h_MuonPOG = (TH1D*)f_in->Get("hData_MuonPOGSystErr12")->Clone("h_MuonPOG");
  }
  if (BaseName == "WpToEleNu" || BaseName == "WmToEleNu")
  {
    h_Binning = (TH1D*)f_in->Get("hData_BinningSystErr12")->Clone("h_Binning");
  }
  TH1D* h_MetResol = (TH1D*)f_in->Get("hData_MetResolSystErr12")->Clone("h_MetResol");
  TH1D* h_EnMomScale = (TH1D*)f_in->Get("hData_EnMomScaleSystErr12")->Clone("h_EnMomScale");
  TH1D* h_EnMomSmear = (TH1D*)f_in->Get("hData_EnMomSmearSystErr12")->Clone("h_EnMomSmear");
  TH1D* h_QCDBkg = (TH1D*)f_in->Get("hData_QcdBckgrSystErr12")->Clone("h_QCDBkg");
  TH1D* h_QCDShape = (TH1D*)f_in->Get("hData_QcdShapeSystErr12")->Clone("h_QCDShape");
  TH1D* h_EWK = (TH1D*)f_in->Get("hData_EwkSystErr12")->Clone("h_EWK");
  TH1D* h_FSR = (TH1D*)f_in->Get("hData_FsrSystErr12")->Clone("h_FSR");
  TH1D* h_SVDUnf = (TH1D*)f_in->Get("hData_SvdUnfSystErr12")->Clone("h_SVDUnf");
  TH1D* h_UnfBias = (TH1D*)f_in->Get("hData_UnfoldBiasSystErr12")->Clone("h_UnfBias");

  TH1D* h_TotalSyst = new TH1D("h_TotalSyst","h_TotalSyst",12,0,12);
  TH1D* h_TotalUncer = new TH1D("h_TotalUncer","h_TotalUncer",12,0,12);

  for(int i(1); i<13; i++)
  {
    if (BaseName == "WpToMuNu" || BaseName == "WmToMuNu")
    {
    h_TotalSyst->SetBinContent(i, sqrt(
	  TMath::Power(h_Toy->GetBinContent(i),2)
	  + TMath::Power(h_IDIsoSig->GetBinContent(i),2)
	  + TMath::Power(h_IDIsoBkg->GetBinContent(i),2)
	  + TMath::Power(h_TrackSig->GetBinContent(i),2)
	  + TMath::Power(h_TrackBkg->GetBinContent(i),2)
	  + TMath::Power(h_MuonPOG->GetBinContent(i),2)
	  + TMath::Power(h_MetResol->GetBinContent(i),2)
	  + TMath::Power(h_EnMomScale->GetBinContent(i),2)
	  + TMath::Power(h_EnMomSmear->GetBinContent(i),2)
	  + TMath::Power(h_QCDBkg->GetBinContent(i),2)
	  + TMath::Power(h_QCDShape->GetBinContent(i),2)
	  + TMath::Power(h_EWK->GetBinContent(i),2)
	  + TMath::Power(h_FSR->GetBinContent(i),2)
	  + TMath::Power(h_SVDUnf->GetBinContent(i),2)
	  + TMath::Power(h_UnfBias->GetBinContent(i),2)
	  ));
    }
    if (BaseName == "WpToEleNu" || BaseName == "WmToEleNu")
    {
    h_TotalSyst->SetBinContent(i, sqrt(
	  TMath::Power(h_Toy->GetBinContent(i),2)
	  + TMath::Power(h_IDIsoSig->GetBinContent(i),2)
	  + TMath::Power(h_IDIsoBkg->GetBinContent(i),2)
	  + TMath::Power(h_Binning->GetBinContent(i),2)
	  + TMath::Power(h_MetResol->GetBinContent(i),2)
	  + TMath::Power(h_EnMomScale->GetBinContent(i),2)
	  + TMath::Power(h_EnMomSmear->GetBinContent(i),2)
	  + TMath::Power(h_QCDBkg->GetBinContent(i),2)
	  + TMath::Power(h_QCDShape->GetBinContent(i),2)
	  + TMath::Power(h_EWK->GetBinContent(i),2)
	  + TMath::Power(h_FSR->GetBinContent(i),2)
	  + TMath::Power(h_SVDUnf->GetBinContent(i),2)
	  + TMath::Power(h_UnfBias->GetBinContent(i),2)
	  ));
    }
    cout << "hTotalSyst : " << h_TotalSyst->GetBinContent(i) << endl;
  }

  // Total Uncertainty
  for(int i(1); i<13; i++)
  {
    h_TotalUncer->SetBinContent(i,sqrt(TMath::Power(h_TotalSyst->GetBinContent(i),2) + TMath::Power(h_Stat->GetBinContent(i),2) ) );
    cout << "hTotalUncer : " << h_TotalUncer->GetBinContent(i) << endl;
  }

  // Save TotalUncer to h_NormDiffXsec
  for(int i(1); i<13; i++)
  {
    h_NormDiffXsec->SetBinError(i,h_NormDiffXsec->GetBinContent(i) * h_TotalUncer->GetBinContent(i) * 0.01); // % error converting to number
  }

  // Output file create and histogram save
  TFile* f_out = new TFile(BaseName+"_CommonFid.root","recreate");
  h_NormDiffXsec->Write();

}

