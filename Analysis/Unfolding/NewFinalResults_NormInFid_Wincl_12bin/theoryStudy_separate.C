#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "../../Utils/const.h"

//const int nBins = 14;
//double WptLogBins[nBins] = {1.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};
//double WptBins[nBins] = {0.0,7.5,12.5,17.5,24,30,40,50,70,110,150,190,250,600};
//
//double ax[13]  = {4.25,10,15,20.75,27,35,45,60,90,130,170,220,425};
//double aex[13] = {3.25,2.5,2.5,3.25,3,5,5,10,20,20,20,30,175};

const int nBins = 13;
double WptLogBins[nBins] = {1.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
double WptBins[nBins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};

double ax[12]  = {4.25,10,15,23.75,35,45,60,90,130,170,220,425};
double aex[12] = {3.25,2.5,2.5,6.25,5,5,10,20,20,20,30,175};


void drawDifference(TH1* iH0,TH1 *iH1,TH1* iH2, TGraphErrors* iH3, int chnl,TGraphErrors* iH4,TGraphAsymmErrors* iH5,TH1* StatErrBand,TGraphErrors* iH6){
  std::string lName = std::string(iH0->GetName());
  TH1F *lHDiff  = new TH1F((lName+"Diff").c_str(),(lName+"Diff").c_str(),nBins-1,WptLogBins);// lHDiff->Sumw2();
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
  TFile *f_Fewz;
  TFile *f_Data;

  f_Resbos = new TFile("Result"+BaseName+"/Result_"+BaseName+"_Theory.root");
  f_Data = new TFile("Result"+BaseName+"/Result_"+BaseName+"_12Bin.root");
  f_Data_PowhegErr = new TFile("Result"+BaseName+"/Result_"+BaseName+"_12Bin.root");
  f_Fewz = new TFile("Result"+BaseName+"/Result_"+BaseName+"_Theory.root");



  TH1D *hWptBins_LinScale   = new TH1D("hWptBins_LinScale","hWptBins_LinScale",nBins-1,WptBins);hWptBins_LinScale->Sumw2();
  TH1D *hWptBins_LogScale   = new TH1D("hWptBins_LogScale","hWptBins_LogScale",nBins-1,WptLogBins);hWptBins_LogScale->Sumw2();
  

  ///=============Reading related stat, syst errors starts here==============
  //  
  TH1D* hData_StatErr		 = (TH1D*)f_Data_PowhegErr->Get("h_Stat")->Clone("hData_StatErr");
  
  //TH1D* hData_TotalReconEff	 = (TH1D*)f_Data_PowhegErr->Get("h_TotalEff")->Clone("hData_TotalReconEff");
  TH1D* hData_EffiToySystErr		 = (TH1D*)f_Data_PowhegErr->Get("h_toy")->Clone("hData_EffiToySystErr");
  TH1D* hData_IDIsoSigShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_idisosig")->Clone("hData_IDIsoSigShapeSystErr");
  TH1D* hData_IDIsoBkgrShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_idisobck")->Clone("hData_IDIsoBkgrShapeSystErr");

  TH1D* hData_MetResolSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_met")->Clone("hData_MetResolSystErr");
  TH1D* hData_EnMomScaleSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_scale")->Clone("hData_EnMomScaleSystErr");
  TH1D* hData_EnMomSmearSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_smear")->Clone("hData_EnMomSmearSystErr");
  TH1D* hData_QcdBckgrSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_qcdbckgr")->Clone("hData_QcdBckgrSystErr");
  TH1D* hData_QcdShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_qcdshape")->Clone("hData_QcdShapeSystErr");
  TH1D* hData_EwkSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_ewk")->Clone("hData_EwkSystErr");
  TH1D* hData_FsrSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_fsr")->Clone("hData_FsrSystErr");
  TH1D* hData_SvdUnfSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_SvdUnf")->Clone("hData_SvdUnfSystErr");
  TH1D* hData_UnfoldBiasSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_UnfoldBias")->Clone("hData_UnfoldBiasSystErr");
  
  ///Lepton Reconstruction related systematic errors
  //
  
  TH1D *hData_TotalSystErr   = new TH1D("hData_TotalSystErr","hData_TotalSystErr",nBins-1,WptBins);hData_TotalSystErr->Sumw2();
  
  //if(BaseName== "WpToMuNu" || BaseName== "WmToMuNu")
  if(BaseName== "WInclToMuNu")
  {
    TH1D* hData_TrackSigShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_tracksig")->Clone("hData_TrackSigShapeSystErr");
    TH1D* hData_TrackBkgrShapeSystErr	 = (TH1D*)f_Data_PowhegErr->Get("h_trackbck")->Clone("hData_TrackBkgrShapeSystErr");
    TH1D* hData_MuonPOGSystErr		 = (TH1D*)f_Data_PowhegErr->Get("h_POG")->Clone("hData_MuonPOGSystErr");
    
    //// Calculate total syst for muon
    for( int ipt(1);ipt<=nBins-1;ipt++)
    {
      ///used for Normalized Xsec
      hData_TotalSystErr->SetBinContent(ipt, sqrt( 
	      TMath::Power( hData_TrackSigShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_TrackBkgrShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MuonPOGSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EffiToySystErr->GetBinContent(ipt),2)
	    + TMath::Power( hData_IDIsoSigShapeSystErr->GetBinContent(ipt),2)
	    + TMath::Power( hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_MetResolSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EnMomScaleSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EnMomSmearSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdBckgrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EwkSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_FsrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_SvdUnfSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_UnfoldBiasSystErr->GetBinContent(ipt),2)  ));
	    
    }
    
  }
  //else if(BaseName== "WpToEleNu" || BaseName== "WmToEleNu")
  else if(BaseName== "WInclToEleNu" )
  {
    TH1D* hData_BinningSystErr		 = (TH1D*)f_Data_PowhegErr->Get("h_bin")->Clone("hData_BinningSystErr");
  
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
	    + TMath::Power( hData_EnMomSmearSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdBckgrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_QcdShapeSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_EwkSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_FsrSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_SvdUnfSystErr->GetBinContent(ipt),2) 
	    + TMath::Power( hData_UnfoldBiasSystErr->GetBinContent(ipt),2) ));
            

    }
  }

  TH1D *hData_TotalUncer   = new TH1D("hData_TotalUncer","hData_TotalUncer",nBins-1,WptBins);hData_TotalUncer->Sumw2();
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hData_TotalUncer->SetBinContent(ipt, sqrt(TMath::Power(hData_StatErr->GetBinContent(ipt),2)+ TMath::Power(hData_TotalSystErr->GetBinContent(ipt),2)));
    cout<<"hData_TotalUncer: "<<ipt<<"\t"<<hData_TotalUncer->GetBinContent(ipt)<<endl;
  } 
 

  for( int ipt(1);ipt<nBins;ipt++)
  {
    //cout<<"hData_TrackSigShapeSystErr: "<<ipt<<"\t"<<hData_TrackSigShapeSystErr->GetBinContent(ipt)<<endl;
    // cout<<"hData_TrackBkgrShapeSystErr: "<<ipt<<"\t"<<hData_TrackBkgrShapeSystErr->GetBinContent(ipt)<<endl;
    // cout<<"hData_IDIsoSigShapeSystErr: "<<ipt<<"\t"<<hData_IDIsoSigShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_IDIsoBkgrShapeSystErr: "<<ipt<<"\t"<<hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_MuonPOGSystErr: "<<ipt<<"\t"<<hData_MuonPOGSystErr->GetBinContent(ipt)<<endl;
   // cout<<"hData_EffiToySystErr: "<<ipt<<"\t"<<hData_EffiToySystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_MetResolSystErr: "<<ipt<<"\t"<<hData_MetResolSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_EnMomScaleSystErr: "<<ipt<<"\t"<<hData_EnMomScaleSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_EnMomSmearSystErr: "<<ipt<<"\t"<<hData_EnMomSmearSystErr->GetBinContent(ipt)<<endl;
   // cout<<"hData_QcdBckgrSystErr: "<<ipt<<"\t"<<hData_QcdBckgrSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_QcdShapeSystErr: "<<ipt<<"\t"<<hData_QcdShapeSystErr->GetBinContent(ipt)<<endl;
   // cout<<"hData_EwkSystErr: "<<ipt<<"\t"<<hData_EwkSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_FsrSystErr: "<<ipt<<"\t"<<hData_FsrSystErr->GetBinContent(ipt)<<endl;
    // cout<<"hData_SvdUnfSystErr: "<<ipt<<"\t"<<hData_SvdUnfSystErr->GetBinContent(ipt)<<endl;
    //cout<<"hData_UnfoldBiasSystErr: "<<ipt<<"\t"<<hData_UnfoldBiasSystErr->GetBinContent(ipt)<<endl;
    //cout<<"Stat error: "<<ipt<<"\t"<<hData_StatErr->GetBinContent(ipt)<<endl;
    //cout<<"Total uncertainty error: "<<ipt<<"\t"<<hData_TotalUncer->GetBinContent(ipt)<<endl;
  } 
  cout << "======= Put in note=======" << endl;
  cout << fixed << setprecision(4);
  for(int ipt(1);ipt<nBins;ipt++)
  {
    //cout<< " & " << hData_TrackSigShapeSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_TrackBkgrShapeSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_IDIsoSigShapeSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_MuonPOGSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_EffiToySystErr->GetBinContent(ipt) ;
   
//    ///Total Lepton Effi syst
//    cout<< " & " << 
//      sqrt(  TMath::Power(hData_TrackSigShapeSystErr->GetBinContent(ipt),2) 
//	  +  TMath::Power(hData_TrackBkgrShapeSystErr->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_IDIsoSigShapeSystErr->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_IDIsoBkgrShapeSystErr->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_MuonPOGSystErr->GetBinContent(ipt),2)
//	  +  TMath::Power(hData_EffiToySystErr->GetBinContent(ipt),2)
//	  ) ;

    //cout<< " & " << hData_EnMomScaleSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_EnMomSmearSystErr->GetBinContent(ipt) ;
    ///Scale Smear
//    cout<< " & " << 
//      sqrt(  TMath::Power(hData_EnMomScaleSystErr->GetBinContent(ipt),2) 
//    	  +  TMath::Power(hData_EnMomSmearSystErr->GetBinContent(ipt),2)
//    	  ) ;

    //cout<< " & " << hData_MetResolSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_QcdBckgrSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_QcdShapeSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_EwkSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_SvdUnfSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_FsrSystErr->GetBinContent(ipt) ;
    //cout<< " & " << hData_UnfoldBiasSystErr->GetBinContent(ipt) ;
    
    //Total syst
    //cout<< " & " << hData_TotalSystErr->GetBinContent(ipt) ;
    // Stat error
    //cout<< " & " << hData_StatErr->GetBinContent(ipt) ;
    //Total Uncertainty
    cout<< " & " << hData_TotalUncer->GetBinContent(ipt) ;



  }
    cout << " \\"<<"\\" << endl;
    cout << endl;


  TH1D* hPowheg_PDFErr		 = (TH1D*)f_Data->Get("h_PowhegPDF")->Clone("hPowheg_PDFErr");
  TH1D* hPowheg_StatErr		 = (TH1D*)f_Data->Get("h_PowhegStat")->Clone("hPowheg_StatErr");
  
  //TH1D* hData_TotalUncer = (TH1D*)f_Data_PowhegErr->Get("h_TotalUncer")->Clone("hData_TotalUncer");
  ///=============Reading related stat, syst errors Finished here==============
 




  /////==========Real Data Starts here=========================================
  //
  // Born Level RD Yield and Errors, and Converting to X sec yield and errors, making TGraph .
  
 

  TH1D* hData_Xsec_BornLinScale;
  hData_Xsec_BornLinScale   = (TH1D*)f_Data->Get("BornEffCorr")->Clone("hData_Xsec_BornLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hData_Xsec_BornLinScale->SetBinError(ipt, 0.01*hData_TotalUncer->GetBinContent(ipt)*hData_Xsec_BornLinScale->GetBinContent(ipt) );
    //hData_Xsec_BornLinScale->SetBinError(ipt, hData_TotalUncer->GetBinContent(ipt) );
  } 
  hData_Xsec_BornLinScale->Scale(1./18.429);

  TH1D* hData_Yield_ReconLinScale;
  hData_Yield_ReconLinScale = (TH1D*)f_Data->Get("data_Rec")->Clone("hData_Yield_ReconLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    double tmp = sqrt(hData_Yield_ReconLinScale->GetBinContent(ipt));
    hData_Yield_ReconLinScale->SetBinError(ipt,tmp);



  }
    cout << " \\"<<"\\" << endl;
    cout << endl;


  TH1D* hPowheg_PDFErr		 = (TH1D*)f_Data->Get("h_PowhegPDF")->Clone("hPowheg_PDFErr");
  TH1D* hPowheg_StatErr		 = (TH1D*)f_Data->Get("h_PowhegStat")->Clone("hPowheg_StatErr");
  
  //TH1D* hData_TotalUncer = (TH1D*)f_Data_PowhegErr->Get("h_TotalUncer")->Clone("hData_TotalUncer");
  ///=============Reading related stat, syst errors Finished here==============
 




  /////==========Real Data Starts here=========================================
  //
  // Born Level RD Yield and Errors, and Converting to X sec yield and errors, making TGraph .
  
 

  TH1D* hData_Xsec_BornLinScale;
  hData_Xsec_BornLinScale   = (TH1D*)f_Data->Get("BornEffCorr")->Clone("hData_Xsec_BornLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hData_Xsec_BornLinScale->SetBinError(ipt, 0.01*hData_TotalUncer->GetBinContent(ipt)*hData_Xsec_BornLinScale->GetBinContent(ipt) );
    //hData_Xsec_BornLinScale->SetBinError(ipt, hData_TotalUncer->GetBinContent(ipt) );
  } 
  hData_Xsec_BornLinScale->Scale(1./18.429);

  TH1D* hData_Yield_ReconLinScale;
  hData_Yield_ReconLinScale = (TH1D*)f_Data->Get("data_Rec")->Clone("hData_Yield_ReconLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    double tmp = sqrt(hData_Yield_ReconLinScale->GetBinContent(ipt));
    hData_Yield_ReconLinScale->SetBinError(ipt,tmp);
  }

  //Errors
  double Data_Xsec_BornStatErr[nBins-1];
  double Data_Xsec_BornTotalErr[nBins-1];
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    Data_Xsec_BornStatErr[ipt] = hData_Xsec_BornLinScale->GetBinContent(ipt) * ( hData_Yield_ReconLinScale->GetBinError(ipt)/hData_Yield_ReconLinScale->GetBinContent(ipt) );
    Data_Xsec_BornTotalErr[ipt] = hData_Xsec_BornLinScale->GetBinError(ipt);
  }

  // X-sec and error Log scale
  TH1D *hData_Xsec_BornLogScale     = new TH1D("hData_Xsec_BornLogScale","hData_Xsec_BornLogScale",nBins-1,WptLogBins);hData_Xsec_BornLogScale->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hData_Xsec_BornLogScale->SetBinContent(ipt,hData_Xsec_BornLinScale->GetBinContent(ipt)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt    ));
    hData_Xsec_BornLogScale->SetBinError(ipt,Data_Xsec_BornTotalErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  }
  
  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataErrBand = new TH1D("hRatioDataErrBand","hRatioDataErrBand",nBins-1,WptLogBins);hRatioDataErrBand->Sumw2();
  TH1D *hRatioDataStatErr = new TH1D("hRatioDataStatErr","hRatioDataStatErr",nBins-1,WptLogBins);hRatioDataStatErr->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioDataErrBand->SetBinContent(ipt,1.);
    hRatioDataErrBand->SetBinError(ipt,hData_Xsec_BornLinScale->GetBinError(ipt)/hData_Xsec_BornLinScale->GetBinContent(ipt));
    
    hRatioDataStatErr->SetBinContent(ipt,1.);
    hRatioDataStatErr->SetBinError(ipt,Data_Xsec_BornStatErr[ipt]/hData_Xsec_BornLinScale->GetBinContent(ipt));
  }
  

  ////Normalize Data X-sec and its errors
  TH1D* hData_Xsec_BornLogScaleNorm = (TH1D*)hData_Xsec_BornLogScale->Clone("hData_Xsec_BornLogScaleNorm");
  //hData_Xsec_BornLogScaleNorm->Scale(1./hData_Xsec_BornLogScale->Integral());
  hData_Xsec_BornLogScaleNorm->Scale(1./hData_Xsec_BornLinScale->Integral());
  
  TFile f_out("Wpt_NormDiffXsec_InFid_RDResBosPowhegFEWZ.root","recreate");
  hData_Xsec_BornLogScaleNorm->Write();

  TH1D *hRatioDataStatErrNorm = (TH1D*)hRatioDataStatErr->Clone("hRatioDataStatErrNorm");
  //hRatioDataStatErrNorm->Scale(1./hData_Xsec_BornLogScale->Integral());
  hRatioDataStatErrNorm->Scale(1./hData_Xsec_BornLinScale->Integral());
  
  TH1D *hRatioDataErrBandNorm = (TH1D*)hRatioDataErrBand->Clone("hRatioDataErrBandNorm");
  //hRatioDataErrBandNorm->Scale(1./hData_Xsec_BornLogScale->Integral());
  hRatioDataErrBandNorm->Scale(1./hData_Xsec_BornLinScale->Integral());

  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioDataStatErrNorm->SetBinContent(ipt,1.);
    //hRatioDataStatErrNorm->SetBinError(ipt,Data_Xsec_BornStatErr[ipt]/hData_Xsec_BornLinScale->GetBinContent(ipt)/hData_Xsec_BornLinScale->Integral());
    
    hRatioDataErrBandNorm->SetBinContent(ipt,1.);
    //hRatioDataErrBandNorm->SetBinError(ipt,hData_Xsec_BornLinScale->GetBinError(ipt)/hData_Xsec_BornLinScale->GetBinContent(ipt)/hData_Xsec_BornLinScale->Integral());
    //cout<<"Data Xsec After Norm: "<<hData_Xsec_BornLogScaleNorm->GetBinContent(ipt)<<"\t"<<hData_Xsec_BornLogScaleNorm->GetBinContent(ipt)<<endl;
  }


  hData_Xsec_BornLogScale->SetMarkerStyle(kFullCircle); 
  hData_Xsec_BornLogScale->SetMarkerColor(kBlack); 
  hData_Xsec_BornLogScale->SetMarkerSize(1);
  
  hData_Xsec_BornLogScaleNorm->SetMarkerStyle(kFullCircle);
  hData_Xsec_BornLogScaleNorm->SetMarkerColor(kBlack);
  hData_Xsec_BornLogScaleNorm->SetMarkerSize(1);


  hRatioDataStatErr->SetMarkerStyle(kFullCircle); 
  hRatioDataStatErr->SetMarkerColor(kBlack); 
  hRatioDataStatErr->SetMarkerSize(0.6);
 
  hRatioDataStatErrNorm->SetMarkerStyle(kFullCircle);
  hRatioDataStatErrNorm->SetMarkerColor(kBlack);
  hRatioDataStatErrNorm->SetMarkerSize(0.6);

  
  /// TGraph
  TGraphErrors *Data_Xsec_Born = new TGraphErrors(hData_Xsec_BornLogScale);
  TGraphErrors *RatioDataErrBand = new TGraphErrors(hRatioDataErrBand);

  //Normalized
  TGraphErrors *Data_Xsec_BornNorm = new TGraphErrors(hData_Xsec_BornLogScaleNorm);
  TGraphErrors *RatioDataErrBandNorm = new TGraphErrors(hRatioDataErrBandNorm);
  
  /////==========Real Data related stuff Finished here=========================================



  /////==========Powheg Starts here==============================================
  //
  // Born Level Powheg Yield and Errors, and Converting to X sec yield and errors, making TGraph .
  
  
  TH1D* hPowheg_Xsec_BornLinScale;
  hPowheg_Xsec_BornLinScale   = (TH1D*)f_Data->Get("SVD_Born.Gen")->Clone("hPowheg_Xsec_BornLinScale");
  for( int ipt(1);ipt<nBins;ipt++)
  {
    hPowheg_Xsec_BornLinScale->SetBinError(ipt, hPowheg_PDFErr->GetBinContent(ipt) );
  } 
  //hPowheg_Xsec_BornLinScale->Scale(1./18.429);
    
  TH1D* hPowheg_Yield_BornAfterFidCut;
  hPowheg_Yield_BornAfterFidCut = (TH1D*)f_Data->Get("hPowheg_Yield_BornAfterFidCut")->Clone("hPowheg_Yield_BornAfterFidCut");

  //Errors
  double Powheg_Xsec_BornStatErr[nBins-1];
  double Powheg_Xsec_BornTotalErr[nBins-1];
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    Powheg_Xsec_BornStatErr[ipt] = hPowheg_StatErr->GetBinContent(ipt);
    Powheg_Xsec_BornTotalErr[ipt] = sqrt(Powheg_Xsec_BornStatErr[ipt]*Powheg_Xsec_BornStatErr[ipt] + hPowheg_Xsec_BornLinScale->GetBinError(ipt)*hPowheg_Xsec_BornLinScale->GetBinError(ipt));
    cout<<"hPowheg_StatErr \t"<<hPowheg_StatErr->GetBinContent(ipt)<<endl;
  }

  // X-sec and error Log scale
  TH1D *hPowheg_Xsec_BornLogScale   = new TH1D("hPowheg_Xsec_BornLogScale","hPowheg_Xsec_BornLogScale",nBins-1,WptLogBins);hPowheg_Xsec_BornLogScale->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hPowheg_Xsec_BornLogScale->SetBinContent(ipt,hPowheg_Xsec_BornLinScale->GetBinContent(ipt)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
    hPowheg_Xsec_BornLogScale->SetBinError(ipt,Powheg_Xsec_BornTotalErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  }

  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioPowhegStatErrBand = new TH1D("hRatioPowhegStatErrBand","hRatioPowhegStatErrBand",nBins-1,WptLogBins);hRatioPowhegStatErrBand->Sumw2();
  TH1D *hRatioPowhegPDFErrBand = new TH1D("hRatioPowhegPDFErrBand","hRatioPowhegPDFErrBand",nBins-1,WptLogBins);hRatioPowhegPDFErrBand->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioPowhegStatErrBand->SetBinContent(ipt,hPowheg_Xsec_BornLogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioPowhegStatErrBand->SetBinError(ipt,Powheg_Xsec_BornStatErr[ipt]/hData_Xsec_BornLinScale->GetBinContent(ipt));
    
    hRatioPowhegPDFErrBand->SetBinContent(ipt,hPowheg_Xsec_BornLogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioPowhegPDFErrBand->SetBinError(ipt,(Powheg_Xsec_BornStatErr[ipt]/hData_Xsec_BornLinScale->GetBinContent(ipt))+(hPowheg_Xsec_BornLinScale->GetBinError(ipt)/hData_Xsec_BornLinScale->GetBinContent(ipt)));
  }
 
  ////Normalize Powheg X-sec and errors
  TH1D *hPowheg_Xsec_BornLogScaleNorm= (TH1D*)hPowheg_Xsec_BornLogScale->Clone("hPowheg_Xsec_BornLogScaleNorm");
  //hPowheg_Xsec_BornLogScaleNorm->Scale(1./hPowheg_Xsec_BornLogScale->Integral());
  hPowheg_Xsec_BornLogScaleNorm->Scale(1./hPowheg_Xsec_BornLinScale->Integral());
  hPowheg_Xsec_BornLogScaleNorm->Write();

  TH1D *hRatioPowhegStatErrBandNorm= (TH1D*)hRatioPowhegStatErrBand->Clone("hRatioPowhegStatErrBandNorm");
  //hRatioPowhegStatErrBandNorm->Scale(1./hPowheg_Xsec_BornLogScale->Integral());

  TH1D *hRatioPowhegPDFErrBandNorm= (TH1D*)hRatioPowhegPDFErrBand->Clone("hRatioPowhegPDFErrBandNorm");
  //hRatioPowhegPDFErrBandNorm->Scale(1./hPowheg_Xsec_BornLogScale->Integral());

  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioPowhegStatErrBandNorm->SetBinContent(ipt,hPowheg_Xsec_BornLogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
    hRatioPowhegPDFErrBandNorm->SetBinContent(ipt,hPowheg_Xsec_BornLogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
  }



  /// TGraph
  TGraphErrors *Powheg_Xsec_Born = new TGraphErrors(hPowheg_Xsec_BornLogScale);
  TGraphErrors *RatioPowhegStatErrBand = new TGraphErrors(hRatioPowhegStatErrBand);
  TGraphErrors *RatioPowhegPDFErrBand = new TGraphErrors(hRatioPowhegPDFErrBand);

  /// Normalized
  TGraphErrors *Powheg_Xsec_BornNorm = new TGraphErrors(hPowheg_Xsec_BornLogScaleNorm);
  TGraphErrors *RatioPowhegStatErrBandNorm = new TGraphErrors(hRatioPowhegStatErrBandNorm);
  TGraphErrors *RatioPowhegPDFErrBandNorm = new TGraphErrors(hRatioPowhegPDFErrBandNorm);
  

  /////==========Powheg related stuff Finished here=========================================
  
  
  
  /////==========FEWZ Starts here==============================================
  //
  // FEWZ X sec value and errors, making TGraph .
  TH1D* hFEWZ_Xsec_LinScale;
  hFEWZ_Xsec_LinScale = (TH1D*)f_Fewz->Get("hxsec")->Clone("hFEWZ_Xsec_LinScale");
  
  //Errors
  TH1D* hFEWZ_Xsec_ScaleError;
  hFEWZ_Xsec_ScaleError = (TH1D*)f_Fewz->Get("ScaleErr")->Clone("hFEWZ_Xsec_ScaleError");
  TH1D* hFEWZ_Xsec_PDFError;
  hFEWZ_Xsec_PDFError = (TH1D*)f_Fewz->Get("PDFErr")->Clone("hFEWZ_Xsec_PDFError");
  
  double FEWZ_Xsec_ScaleErr[nBins-1];
  double FEWZ_Xsec_PDFErr[nBins-1];
  double FEWZ_Xsec_TotalErr[nBins-1];
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    FEWZ_Xsec_ScaleErr[ipt] = hFEWZ_Xsec_ScaleError->GetBinError(ipt);
    FEWZ_Xsec_PDFErr[ipt]   = hFEWZ_Xsec_PDFError->GetBinError(ipt);
    
    FEWZ_Xsec_TotalErr[ipt] = sqrt( FEWZ_Xsec_ScaleErr[ipt]*FEWZ_Xsec_ScaleErr[ipt] + FEWZ_Xsec_PDFErr[ipt]*FEWZ_Xsec_PDFErr[ipt] );
  }
  
  // X-sec and error Log scale
  TH1D *hFEWZ_Xsec_LogScale   = new TH1D("hFEWZ_Xsec_LogScale","hFEWZ_Xsec_LogScale",nBins-1,WptLogBins);hFEWZ_Xsec_LogScale->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hFEWZ_Xsec_LogScale->SetBinContent(ipt,hFEWZ_Xsec_LinScale->GetBinContent(ipt)/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
    hFEWZ_Xsec_LogScale->SetBinError(ipt,FEWZ_Xsec_TotalErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  } 
  ///   Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioFEWZStatErrBand = new TH1D("hRatioFEWZStatErrBand","hRatioFEWZStatErrBand",nBins-1,WptLogBins);hRatioFEWZStatErrBand->Sumw2();
  TH1D *hRatioFEWZScaleErr = new TH1D("hRatioFEWZScaleErr","hRatioFEWZScaleErr",nBins-1,WptLogBins);hRatioFEWZScaleErr->Sumw2();
  TH1D *hRatioFEWZQCDScaleErrBand = new TH1D("hRatioFEWZQCDScalePDFErrBand","hRatioFEWZQCDScalePDFErrBand",nBins-1,WptLogBins);hRatioFEWZQCDScalePDFErrBand->Sumw2();
  TH1D *hRatioFEWZScalePDFErrBand = new TH1D("hRatioFEWZScalePDFErrBand","hRatioFEWZScalePDFErrBand",nBins-1,WptLogBins);hRatioFEWZScalePDFErrBand->Sumw2();
  for( int ipt(1);ipt<=nBins-1;ipt++)
  {
    hRatioFEWZStatErrBand->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZStatErrBand->SetBinError(ipt,0.01);
    
    hRatioFEWZScaleErr->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZScaleErr->SetBinError(ipt,FEWZ_Xsec_ScaleErr[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt));
  
    hRatioFEWZQCDScaleErrBand->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZQCDScaleErrBand->SetBinError(ipt,0.01 + hRatioFEWZScaleErr->GetBinError(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
  
    hRatioFEWZScalePDFErrBand->SetBinContent(ipt,hFEWZ_Xsec_LogScale->GetBinContent(ipt)/hData_Xsec_BornLogScale->GetBinContent(ipt));
    hRatioFEWZScalePDFErrBand->SetBinError(ipt,0.01+(hFEWZ_Xsec_LogScale->GetBinError(ipt)+hRatioFEWZScaleErr->GetBinError(ipt))/hData_Xsec_BornLogScale->GetBinContent(ipt));
  }
 
  ////Normalize FEWZ X-sec and errors
  TH1D *hFEWZ_Xsec_LogScaleNorm = (TH1D*)hFEWZ_Xsec_LogScale->Clone("hFEWZ_Xsec_LogScaleNorm");
  //hFEWZ_Xsec_LogScaleNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());
  hFEWZ_Xsec_LogScaleNorm->Scale(1./hFEWZ_Xsec_LinScale->Integral());
  hFEWZ_Xsec_LogScaleNorm->Write();

  TH1D *hRatioFEWZStatErrBandNorm = (TH1D*)hRatioFEWZStatErrBand->Clone("hRatioFEWZStatErrBandNorm");
  //hRatioFEWZStatErrBandNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());

  TH1D *hRatioFEWZQCDScaleErrBandNorm = (TH1D*)hRatioFEWZQCDScaleErrBand->Clone("hRatioFEWZQCDScaleErrBandNorm");
  //hRatioFEWZQCDScaleErrBandNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());

  TH1D *hRatioFEWZScalePDFErrBandNorm = (TH1D*)hRatioFEWZScalePDFErrBand->Clone("hRatioFEWZScalePDFErrBandNorm");
  //hRatioFEWZScalePDFErrBandNorm->Scale(1./hFEWZ_Xsec_LogScale->Integral());

   for( int ipt(1);ipt<=nBins-1;ipt++)
   {
     hRatioFEWZStatErrBandNorm->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
     hRatioFEWZQCDScaleErrBandNorm->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
     hRatioFEWZScalePDFErrBandNorm->SetBinContent(ipt,hFEWZ_Xsec_LogScaleNorm->GetBinContent(ipt)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt));
   }


  /// TGraph
  TGraphErrors *FEWZ_Xsec = new TGraphErrors(hFEWZ_Xsec_LogScale);
  TGraphErrors* RatioFEWZStatErrBand = new TGraphErrors(hRatioFEWZStatErrBand);
  TGraphErrors* RatioFEWZQCDScaleErrBand = new TGraphErrors(hRatioFEWZQCDScaleErrBand);
  TGraphErrors* RatioFEWZScalePDFErrBand = new TGraphErrors(hRatioFEWZScalePDFErrBand);

  /// Normalized
  TGraphErrors *FEWZ_XsecNorm = new TGraphErrors(hFEWZ_Xsec_LogScaleNorm);
  TGraphErrors* RatioFEWZStatErrBandNorm = new TGraphErrors(hRatioFEWZStatErrBandNorm);
  TGraphErrors* RatioFEWZQCDScaleErrBandNorm = new TGraphErrors(hRatioFEWZQCDScaleErrBandNorm);
  TGraphErrors* RatioFEWZScalePDFErrBandNorm = new TGraphErrors(hRatioFEWZScalePDFErrBandNorm);


  /////==========FEWZ related stuff Finished here=========================================
  
  
  /////==========ResBos Starts here==============================================
  //
  // ResBos X sec value and errors, making TGraph .
 
  TH1D* hResBos30_CentralYield_LinScale;
  hResBos30_CentralYield_LinScale = (TH1D*)f_Resbos->Get("hResbos30")->Clone("hResBos30_CentralYield_LinScale");

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
    //double nomVal  = AllResbos[1]->GetBinContent(ipt+1);
    double nomVal  = AllResbos[1]->GetBinContent(ipt+1)/(AllResbos[1]->Integral()*hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1));
    Resb_errMax[ipt] = -99999;
    Resb_errMin[ipt] = 990009;
    for (int i(0);i<7;i++)
    {
      //tmpVal  = AllResbos[i]->GetBinContent(ipt+1);
      tmpVal  = AllResbos[i]->GetBinContent(ipt+1)/(AllResbos[i]->Integral()*hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1));
      tmpDiff = tmpVal - nomVal;
      if( tmpDiff > Resb_errMax[ipt] ) Resb_errMax[ipt] = tmpDiff;
      if( tmpDiff < Resb_errMin[ipt] ) Resb_errMin[ipt] = tmpDiff;
    }

    if (Resb_errMax[ipt] < 0) Resb_errMax[ipt] = 0.;
    if (Resb_errMin[ipt] > 0) Resb_errMin[ipt] = 0.;
    if (Resb_errMin[ipt] < 0) Resb_errMin[ipt] = -Resb_errMin[ipt];
    //Resb_errMax[ipt] = Resb_errMax[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
    //Resb_errMin[ipt] = Resb_errMin[ipt]/hWptBins_LinScale->GetXaxis()->GetBinWidth(ipt+1);
    Resb_errMax[ipt] = Resb_errMax[ipt];
    Resb_errMin[ipt] = Resb_errMin[ipt];
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
    //errResbosDataLo[ipt]=Resb_errMin[ipt];
    //errResbosDataHi[ipt]=Resb_errMax[ipt];

  }
  
  
  ////Normalize ResBos X-sec and errors
  double resb30Total=0.0;
  for( int ipt(0);ipt<nBins-1;ipt++)
  {
    //resb30Total=resb30Total + hResb30_CentralXsec[ipt];
    resb30Total=resb30Total + hResBos30_CentralYield_LinScale->GetBinContent(ipt+1);
  }
  TH1D *hResBos30_CentralXSec_LogScaleNorm = (TH1D*)hResBos30_CentralXSec_LogScale->Clone("hResBos30_CentralXSec_LogScaleNorm");
  hResBos30_CentralXSec_LogScaleNorm->Scale(1./resb30Total);
  hResBos30_CentralXSec_LogScaleNorm->Write();

  Double_t hResb30_CentralXsecNorm[nBins-1];
  Double_t Resb_errMaxNorm[nBins-1];
  Double_t Resb_errMinNorm[nBins-1];
  Double_t RatioResbValNorm[nBins-1],errResbosDataLoNorm[nBins-1],errResbosDataHiNorm[nBins-1];
  for( int ipt(0);ipt<nBins-1;ipt++)
  {

    hResb30_CentralXsecNorm[ipt] = hResb30_CentralXsec[ipt]/resb30Total;
    //Resb_errMaxNorm[ipt] = Resb_errMax[ipt]/resb30Total;
    //Resb_errMinNorm[ipt] = Resb_errMin[ipt]/resb30Total;
    Resb_errMaxNorm[ipt] = Resb_errMax[ipt];
    Resb_errMinNorm[ipt] = Resb_errMin[ipt];

    RatioResbVal[ipt]=hResBos30_CentralXSec_LogScaleNorm->GetBinContent(ipt+1)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt+1);
    //RatioResbValNorm[ipt] = RatioResbVal[ipt]/resb30Total;
    
    RatioResbValNorm[ipt]=hResBos30_CentralXSec_LogScaleNorm->GetBinContent(ipt+1)/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt+1);
    errResbosDataLoNorm[ipt] =Resb_errMinNorm[ipt]/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt+1) ;
    errResbosDataHiNorm[ipt] =Resb_errMaxNorm[ipt]/hData_Xsec_BornLogScaleNorm->GetBinContent(ipt+1) ;
  }
 
   
  
  ///   Theory/Data ratio plot errors related to Real Data
  TGraphAsymmErrors *Resb30_CentralXsec = new TGraphAsymmErrors(nBins-1, ax, hResb30_CentralXsec, aex, aex, Resb_errMin, Resb_errMax);
  TGraphAsymmErrors *RatioResbosErrBand = new TGraphAsymmErrors(nBins-1, ax, RatioResbVal, aex, aex, errResbosDataLo, errResbosDataHi);
 ///   Normalized
  TGraphAsymmErrors *Resb30_CentralXsecNorm = new TGraphAsymmErrors(nBins-1, ax, hResb30_CentralXsecNorm, aex, aex, Resb_errMinNorm, Resb_errMaxNorm);
  TGraphAsymmErrors *RatioResbosErrBandNorm = new TGraphAsymmErrors(nBins-1, ax, RatioResbValNorm, aex, aex, errResbosDataLoNorm, errResbosDataHiNorm);
 

  /////==========ResBos related stuff Finished here=========================================
  
  
  
  
  //// Now design and Draw 
  
  //Data

  RatioDataErrBand->SetFillColor(kBlack);
  RatioDataErrBand->SetFillStyle(3354);
  
  RatioDataErrBandNorm->SetFillColor(kBlack);
  RatioDataErrBandNorm->SetFillStyle(3354);

  //Powheg
  Powheg_Xsec_Born->SetFillColor(kRed);
  Powheg_Xsec_Born->SetFillStyle(3345);

  Powheg_Xsec_BornNorm->SetFillColor(kRed);
  Powheg_Xsec_BornNorm->SetFillStyle(3345);

  RatioPowhegStatErrBand->SetFillColor(kRed-10);
  RatioPowhegStatErrBand->SetFillStyle(3001);

  RatioPowhegStatErrBandNorm->SetFillColor(kRed-10);
  RatioPowhegStatErrBandNorm->SetFillStyle(3001);

  RatioPowhegPDFErrBand->SetFillColor(kRed+2);
  RatioPowhegPDFErrBand->SetFillStyle(3001);

  RatioPowhegPDFErrBandNorm->SetFillColor(kRed+2);
  RatioPowhegPDFErrBandNorm->SetFillStyle(3001);


  //FEWZ
  FEWZ_Xsec->SetFillColor(kGreen);
  FEWZ_Xsec->SetFillStyle(3444);
 
  FEWZ_XsecNorm->SetFillColor(kGreen);
  FEWZ_XsecNorm->SetFillStyle(3444);

  RatioFEWZStatErrBand->SetFillColor(kGreen);
  RatioFEWZStatErrBand->SetFillStyle(3001);

  RatioFEWZStatErrBandNorm->SetFillColor(kGreen);
  RatioFEWZStatErrBandNorm->SetFillStyle(3001);

  RatioFEWZQCDScaleErrBand->SetFillColor(kGreen+7);
  RatioFEWZQCDScaleErrBand->SetFillStyle(3001);

  RatioFEWZQCDScaleErrBandNorm->SetFillColor(kGreen+7);
  RatioFEWZQCDScaleErrBandNorm->SetFillStyle(3001);

  RatioFEWZScalePDFErrBand->SetFillColor(kGreen+3);
  RatioFEWZScalePDFErrBand->SetFillStyle(3001);

  RatioFEWZScalePDFErrBandNorm->SetFillColor(kGreen+3);
  RatioFEWZScalePDFErrBandNorm->SetFillStyle(3001);

  
  //ResBos
  Resb30_CentralXsec->SetFillColor(kBlue);
  Resb30_CentralXsec->SetFillStyle(3354);

  Resb30_CentralXsecNorm->SetFillColor(kBlue);
  Resb30_CentralXsecNorm->SetFillStyle(3354);

  RatioResbosErrBand->SetFillColor(kBlue-7);
  RatioResbosErrBand->SetFillStyle(3001);
 
  RatioResbosErrBandNorm->SetFillColor(kBlue-7);
  RatioResbosErrBandNorm->SetFillStyle(3001);


  TLegend *lL =new TLegend(0.2,0.2,0.52,0.38); lL->SetFillColor(0); lL->SetBorderSize(0);
  //lL->AddEntry(Data_Xsec_Born,"data","PL");
  //lL->AddEntry(Powheg_Xsec_Born,"POWHEG CT10 NLO","f");
  //lL->AddEntry(FEWZ_Xsec,"FEWZ CT10 NNLO","f");
  //lL->AddEntry(Resb30_CentralXsec,"ResBos CT10 NNLL","f");
  lL->AddEntry(Data_Xsec_BornNorm,"data","PL");
  lL->AddEntry(Powheg_Xsec_BornNorm,"POWHEG CT10 NLO","f");
  lL->AddEntry(FEWZ_XsecNorm,"FEWZ CT10 NNLO","f");
  lL->AddEntry(Resb30_CentralXsecNorm,"ResBos CT10 NNLL","f");

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

  if (BaseName=="WInclToMuNu")
     tb->AddText("W #rightarrow #mu #nu");
  if (BaseName=="WInclToEleNu")
     tb->AddText("W #rightarrow e #nu");

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

  //Powheg_Xsec_Born->GetYaxis()->SetRangeUser(1e-3,5*hResb30_CentralXsec[0]);
  //Powheg_Xsec_Born->SetTitle("");
  //Powheg_Xsec_Born->GetYaxis()->SetTitle("d#sigma/dp_{T}^{W} [pb/GeV]");
  //Powheg_Xsec_Born->GetYaxis()->SetTitleOffset(1.8);
  //Powheg_Xsec_Born->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  Powheg_Xsec_BornNorm->GetYaxis()->SetRangeUser(1e-7,5*hResb30_CentralXsecNorm[0]);
  Powheg_Xsec_BornNorm->SetTitle("");
  Powheg_Xsec_BornNorm->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T}^{W} [pb/GeV]");
  Powheg_Xsec_BornNorm->GetYaxis()->SetTitleOffset(1.8);
  Powheg_Xsec_BornNorm->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  
  TPaveText *cmspre = new TPaveText(0.63,0.95,0.95,0.99,"NDC");
  cmspre->SetBorderSize(0);
  cmspre->SetFillStyle(0);
  cmspre->AddText("CMS Preliminary");
  //Powheg_Xsec_Born->Draw("A2");
  //FEWZ_Xsec->Draw("2");
  //Resb30_CentralXsec->Draw("2");
  //Data_Xsec_Born->Draw("p");
  Powheg_Xsec_BornNorm->Draw("A2");
  FEWZ_XsecNorm->Draw("2");
  Resb30_CentralXsecNorm->Draw("2");
  Data_Xsec_BornNorm->Draw("p");

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
  
  if(BaseName=="WInclToMuNu")
    sprintf(tmpName,"winclmnCrossSecNormInFiducial.png");
  if(BaseName=="WInclToEleNu")
    sprintf(tmpName,"winclenCrossSecNormInFiducial.png");
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
  rL1->AddEntry(RatioResbosErrBandNorm,"Theory syst","F");
  rL1->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  hRatioDataStatErr->SetTitle("");
  rL1->AddEntry(hRatioDataStatErr,"","");
  rL1->AddEntry(RatioDataErrBand,"Data stat+syst","F");
  rL1->SetTextSize(0.07);
  
  TLegend *tL1 =new TLegend(0.17,0.72,0.37,0.82); tL1->SetFillColor(0); tL1->SetBorderSize(0);
  //tL1->AddEntry(RatioResbosErrBand,"ResBos","F");
  tL1->AddEntry(RatioResbosErrBandNorm,"ResBos","F");
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
 
  if (BaseName=="WInclToMuNu")
      tb4->AddText("W #rightarrow #mu #nu");
  if (BaseName=="WInclToEleNu")
      tb4->AddText("W #rightarrow e #nu");

  
  //drawDifference(hResBos30_CentralXSec_LogScaleNorm, hData_Xsec_BornLogScaleNorm, hRatioDataErrBand, RatioPowhegStatErrBand,1,RatioPowhegStatErrBand,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
  drawDifference(hResBos30_CentralXSec_LogScaleNorm, hData_Xsec_BornLogScaleNorm, hRatioDataErrBand, RatioPowhegStatErrBand,1,RatioPowhegStatErrBand,RatioResbosErrBandNorm,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
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
  rL2->AddEntry(RatioPowhegPDFErrBandNorm,"PDF    ","F");
  rL2->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  //rL2->AddEntry(RatioPowhegStatErrBand,"stat","F");
  rL2->AddEntry(RatioPowhegStatErrBandNorm,"stat","F");
  rL2->AddEntry(RatioDataErrBand,"Data stat+syst","F");
  rL2->SetTextSize(0.07);

  TLegend *tL2 =new TLegend(0.17,0.85,0.37,0.95); tL2->SetFillColor(0); tL2->SetBorderSize(0);
  //tL2->AddEntry(RatioPowhegPDFErrBand,"POWHEG","F");
  tL2->AddEntry(RatioPowhegPDFErrBandNorm,"POWHEG","F");
  tL2->SetTextSize(0.12);
  tL2->SetTextFont(2);
  
  //drawDifference(hPowhegLog,hDataLog,hDataErrBand,pRatio,2,pRatioPDF,ResbosErrBand,hStatErr,fScaleRatio);
  //drawDifference(Powheg_Xsec_Born,hData_Xsec_BornLogScale,hRatioDataErrBand,RatioPowhegStatErrBand,2,RatioPowhegPDFErrBand,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
  drawDifference(hPowheg_Xsec_BornLogScaleNorm,hData_Xsec_BornLogScaleNorm,hRatioDataErrBand,RatioPowhegStatErrBandNorm,2,RatioPowhegPDFErrBandNorm,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBand);
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
  rL3->AddEntry(RatioFEWZScalePDFErrBandNorm,"PDF","F");
  rL3->AddEntry(hRatioDataStatErr,"Data stat","PLE1");
  //rL3->AddEntry(RatioFEWZQCDScaleErrBand,"QCD scales","F");
  rL3->AddEntry(RatioFEWZQCDScaleErrBandNorm,"QCD scales","F");
  rL3->AddEntry(RatioDataErrBand,"Data stat+syst","F");
  //rL3->AddEntry(RatioFEWZStatErrBand,"stat","F");
  rL3->AddEntry(RatioFEWZStatErrBandNorm,"stat","F");
  rL3->SetTextSize(0.07);

  TLegend *tL3 =new TLegend(0.17,0.85,0.37,0.95); tL3->SetFillColor(0); tL3->SetBorderSize(0);
  //tL3->AddEntry(RatioFEWZScalePDFErrBand,"FEWZ","F");
  tL3->AddEntry(RatioFEWZScalePDFErrBandNorm,"FEWZ","F");
  tL3->SetTextSize(0.12);
  tL3->SetTextFont(2);
  
  drawDifference(hFEWZ_Xsec_LogScaleNorm,hData_Xsec_BornLogScaleNorm,hRatioDataErrBand,RatioFEWZStatErrBandNorm,3,RatioFEWZScalePDFErrBandNorm,RatioResbosErrBand,hRatioDataStatErr,RatioFEWZQCDScaleErrBandNorm);
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

  if(BaseName=="WInclToMuNu")
     sprintf(tmpName,"winclmnRatioTheoDataNormInFiducial.png");
  if(BaseName=="WInclToEleNu")
     sprintf(tmpName,"winclenRatioTheoDataNormInFiducial.png");
  lC2->SaveAs(tmpName);

  return 0;
}
