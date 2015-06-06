#include "../Utils/MyTools.hh"
#include "../Utils/CPlot.hh"
#include "../Utils/MitStyleRemix.hh"
#include "../Utils/const.h"
#include <TFile.h>

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooChi2Var.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"

//Make MET difference plots
void WZ_ZmassComp_Gaus(const TString Mode,//Channel - Muon or Electron
    const TString corrName,
    const TString outputDir 
    )
{
  TH1D* makeDiffHist(TH1D* h1, TH1D* h2, const TString name);
  const TString format("png"); 

  TFile *fnameW_MC;
  TFile *fnameW_RD;
  TFile *fnameZ;

  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;

  fnameW_MC = new TFile("Muon2012LoPU/Muon_DYToMuMu_S8.root");
  fnameW_RD = new TFile("Muon2012LoPU/Muon_RD_LowPU.root");
  fnameZ = new TFile("histo_Zmass_Corr.root");

  CPlot *plotMC;
  CPlot *plotRD;
  CPlot *plotMeanMC;
  CPlot *plotMeanRD;
  CPlot *plotWidthMC;
  CPlot *plotWidthRD;
  
  CPlot *plotMeanMC_Wpt;
  CPlot *plotMeanMC_Zpt;
  CPlot *plotWidthMC_Zpt;
  CPlot *plotWidthMC_Wpt;

  TH1D *hMC_Wpt[4];
  TH1D *hRD_Wpt[4];
  TH1D *hMC_Zpt[4];
  TH1D *hRD_Zpt[4];

  RooDataHist *h1_MC_Wpt;
  RooDataHist *h1_RD_Wpt;
  RooDataHist *h1_MC_Zpt;
  RooDataHist *h1_RD_Zpt;

  TH1D *h_meanMC_Wpt = new TH1D("h_meanMC_Wpt","h_meanMC_Wpt",4,0,4); h_meanMC_Wpt->Sumw2();
  TH1D *h_meanMC_Zpt = new TH1D("h_meanMC_Zpt","h_meanMC_Zpt",4,0,4); h_meanMC_Zpt->Sumw2();
  TH1D *h_meanRD_Wpt = new TH1D("h_meanRD_Wpt","h_meanRD_Wpt",4,0,4); h_meanRD_Wpt->Sumw2();
  TH1D *h_meanRD_Zpt = new TH1D("h_meanRD_Zpt","h_meanRD_Zpt",4,0,4); h_meanRD_Zpt->Sumw2();

  TH1D *h_widthMC_Wpt = new TH1D("h_widthMC_Wpt","h_widthMC_Wpt",4,0,4); h_widthMC_Wpt->Sumw2();
  TH1D *h_widthMC_Zpt = new TH1D("h_widthMC_Zpt","h_widthMC_Zpt",4,0,4); h_widthMC_Zpt->Sumw2();
  TH1D *h_widthRD_Wpt = new TH1D("h_widthRD_Wpt","h_widthRD_Wpt",4,0,4); h_widthRD_Wpt->Sumw2();
  TH1D *h_widthRD_Zpt = new TH1D("h_widthRD_Zpt","h_widthRD_Zpt",4,0,4); h_widthRD_Zpt->Sumw2();

  char histName[30];
  char tmpName[30];

  TString plotTitle = "MC Total, before correction";

  TCanvas *myCan;
  myCan = MakeCanvas("myCan","myCan",800,600);

  //Loop start
  //i=0-Total, 1-BB, 2-BE, 3-EE
  for(int i(0);i<4;i++){
    sprintf(histName,"hMC_Wpt_%d",i);
    sprintf(tmpName,"h1_Zmass");
    if(i==1) sprintf(tmpName,"h1_Zmass_BB");
    if(i==2) sprintf(tmpName,"h1_Zmass_BE");
    if(i==3) sprintf(tmpName,"h1_Zmass_EE");

    if(corrName=="CorrTotalRegion"){
      sprintf(tmpName,"h1_ZmassCorr");
      if(i==1) sprintf(tmpName,"h1_ZmassCorr_BB");
      if(i==2) sprintf(tmpName,"h1_ZmassCorr_BE");
      if(i==3) sprintf(tmpName,"h1_ZmassCorr_EE");
    }

    hMC_Wpt[i]= (TH1D*)fnameW_MC->Get(tmpName)->Clone(histName); hMC_Wpt[i]->Sumw2();

    sprintf(histName,"hRD_Wpt_%d",i);
    hRD_Wpt[i]= (TH1D*)fnameW_RD->Get(tmpName)->Clone(histName); hRD_Wpt[i]->Sumw2();

    sprintf(histName,"hMC_Zpt_%d",i);
    sprintf(tmpName,"h1_ZmassMC");
    if(i==1) sprintf(tmpName,"h1_ZmassMC_BB");
    if(i==2) sprintf(tmpName,"h1_ZmassMC_BE");
    if(i==3) sprintf(tmpName,"h1_ZmassMC_EE");
    
    if(corrName=="CorrTotalRegion"){
      sprintf(tmpName,"h1_ZmassMCcorr");
      if(i==1) sprintf(tmpName,"h1_ZmassMCcorr_BB");
      if(i==2) sprintf(tmpName,"h1_ZmassMCcorr_BE");
      if(i==3) sprintf(tmpName,"h1_ZmassMCcorr_EE");
    }

    hMC_Zpt[i]= (TH1D*)fnameZ->Get(tmpName)->Clone(histName); hMC_Zpt[i]->Sumw2();

    sprintf(histName,"hRD_Zpt_%d",i);
    sprintf(tmpName,"h1_ZmassData");
    if(i==1) sprintf(tmpName,"h1_ZmassData_BB");
    if(i==2) sprintf(tmpName,"h1_ZmassData_BE");
    if(i==3) sprintf(tmpName,"h1_ZmassData_EE");
    
    if(corrName=="CorrTotalRegion"){
      sprintf(tmpName,"h1_ZmassData");
      if(i==1) sprintf(tmpName,"h1_ZmassData_BB");
      if(i==2) sprintf(tmpName,"h1_ZmassData_BE");
      if(i==3) sprintf(tmpName,"h1_ZmassData_EE");
    }

    hRD_Zpt[i]= (TH1D*)fnameZ->Get(tmpName)->Clone(histName); hRD_Zpt[i]->Sumw2();

    hMC_Wpt[i] -> SetMarkerSize(0.9);
    hMC_Wpt[i] -> SetMarkerColor(kRed);
    hMC_Wpt[i] -> SetLineColor(kRed);
    hRD_Wpt[i] -> SetMarkerSize(0.9);
    hRD_Wpt[i] -> SetMarkerColor(kRed);
    hRD_Wpt[i] -> SetLineColor(kRed);
    
    hMC_Zpt[i] -> SetMarkerSize(0.9);
    hMC_Zpt[i] -> SetMarkerColor(kBlue);
    hMC_Zpt[i] -> SetLineColor(kBlue);
    hRD_Zpt[i] -> SetMarkerSize(0.9);
    hRD_Zpt[i] -> SetMarkerColor(kBlue);
    hRD_Zpt[i] -> SetLineColor(kBlue);
    
    RooRealVar x("x", "x",80,100);
    x.setBins(20);
    x.setRange("R0",86,96);

    h1_MC_Wpt = new RooDataHist("h1_MC_Wpt","h1_MC_Wpt",RooArgSet(x),hMC_Wpt[i]);
    h1_RD_Wpt = new RooDataHist("h1_RD_Wpt","h1_RD_Wpt",RooArgSet(x),hRD_Wpt[i]);
    h1_MC_Zpt = new RooDataHist("h1_MC_Zpt","h1_MC_Zpt",RooArgSet(x),hMC_Zpt[i]);
    h1_RD_Zpt = new RooDataHist("h1_RD_Zpt","h1_RD_Zpt",RooArgSet(x),hRD_Zpt[i]);

    //==================================
    //Gauss function
    //==================================
    RooRealVar meanMC_Wpt("meanMC_Wpt","",91.2,80,100);
    RooRealVar meanRD_Wpt("meanRD_Wpt","",91.2,80,100);
    RooRealVar meanMC_Zpt("meanMC_Zpt","",91.2,80,100);
    RooRealVar meanRD_Zpt("meanRD_Zpt","",91.2,80,100);
    RooRealVar sigmaMC_Wpt("sigmaMC_Wpt","",2,0.,5.);
    RooRealVar sigmaRD_Wpt("sigmaRD_Wpt","",2,0.,5.);
    RooRealVar sigmaMC_Zpt("sigmaMC_Wpt","",2,0.,5.);
    RooRealVar sigmaRD_Zpt("sigmaRD_Wpt","",2,0.,5.);

    RooGaussian mcModel_Wpt("mcModel_Wpt", "",x,meanMC_Wpt,sigmaMC_Wpt);
    RooGaussian mcModel_Zpt("mcModel_Zpt", "",x,meanMC_Zpt,sigmaMC_Zpt);
    RooGaussian dataModel_Wpt("dataModel_Wpt", "",x,meanRD_Wpt,sigmaRD_Wpt);
    RooGaussian dataModel_Zpt("dataModel_Zpt", "",x,meanRD_Zpt,sigmaRD_Zpt);
	         
    //==================================
    //Fit Zmass distribution with Gaussian function
    //==================================
    mcModel_Wpt.fitTo(*h1_MC_Wpt,Range("R0"));
    mcModel_Zpt.fitTo(*h1_MC_Zpt,Range("R0"));
    dataModel_Wpt.fitTo(*h1_RD_Wpt,Range("R0"));
    dataModel_Zpt.fitTo(*h1_RD_Zpt,Range("R0"));

    //==================================
    //Fill histograms: mean and width values
    //==================================
    h_meanMC_Wpt->SetBinContent(i+1,meanMC_Wpt.getVal());
    h_meanRD_Wpt->SetBinContent(i+1,meanRD_Wpt.getVal());
    h_meanMC_Zpt->SetBinContent(i+1,meanMC_Zpt.getVal());
    h_meanRD_Zpt->SetBinContent(i+1,meanRD_Zpt.getVal());

    h_meanMC_Wpt->SetBinError(i+1,meanMC_Wpt.getError());
    h_meanRD_Wpt->SetBinError(i+1,meanRD_Wpt.getError());
    h_meanMC_Zpt->SetBinError(i+1,meanMC_Zpt.getError());
    h_meanRD_Zpt->SetBinError(i+1,meanRD_Zpt.getError());

    h_widthMC_Wpt->SetBinContent(i+1,sigmaMC_Wpt.getVal());
    h_widthRD_Wpt->SetBinContent(i+1,sigmaRD_Wpt.getVal());
    h_widthMC_Zpt->SetBinContent(i+1,sigmaMC_Zpt.getVal());
    h_widthRD_Zpt->SetBinContent(i+1,sigmaRD_Zpt.getVal());

    h_widthMC_Wpt->SetBinError(i+1,sigmaMC_Wpt.getError());
    h_widthRD_Wpt->SetBinError(i+1,sigmaRD_Wpt.getError());
    h_widthMC_Zpt->SetBinError(i+1,sigmaMC_Zpt.getError());
    h_widthRD_Zpt->SetBinError(i+1,sigmaRD_Zpt.getError());

    RooPlot* mcframe = x.frame(Bins(20));
    RooPlot* dataframe = x.frame(Bins(20));

    //==================================
    //MC: Draw Wpt and Zpt fit results
    //==================================
    sprintf(tmpName,"MC_Total_noCorr");
    if(i==1){sprintf(tmpName,"MC_BB_noCorr"); plotTitle = "MC BB, before correction";}
    if(i==2){sprintf(tmpName,"MC_BE_noCorr"); plotTitle = "MC BE, before correction";}
    if(i==3){sprintf(tmpName,"MC_EE_noCorr"); plotTitle = "MC EE, before correction";}
    if(corrName=="CorrTotalRegion")
    {
      sprintf(tmpName,"MC_Total_Corr"); plotTitle = "MC Total, after correction";
      if(i==1){sprintf(tmpName,"MC_BB_Corr"); plotTitle = "MC BB, after correction";}
      if(i==2){sprintf(tmpName,"MC_BE_Corr"); plotTitle = "MC BE, after correction";}
      if(i==3){sprintf(tmpName,"MC_EE_Corr"); plotTitle = "MC EE, after correction";}
    }
    
    sprintf(histName,"Events / %.1f",hMC_Wpt[i]->GetBinWidth(1));
    plotMC = new CPlot(tmpName,mcframe,plotTitle,"M_{#mu#mu} [GeV]",histName);
    plotMC->setOutDir(CPlot::sOutDir);
    plotMC->SetLegend(0.65,0.7,0.88,0.82);
    h1_MC_Wpt->plotOn(mcframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    mcModel_Wpt.plotOn(mcframe,LineColor(kRed));
    h1_MC_Zpt->plotOn(mcframe,LineColor(kBlue),MarkerColor(kBlue),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    mcModel_Zpt.plotOn(mcframe,LineColor(kBlue));
    plotMC->GetLegend()->AddEntry(hMC_Wpt[i],"Wpt, Z #rightarrow #mu#mu","pl");
    plotMC->GetLegend()->AddEntry(hMC_Zpt[i],"Zpt, Z #rightarrow #mu#mu","pl");

    //plotMC->AddTextBox("MC fit",0.20,0.88,0.35,0.91,0);
    sprintf(tmpName,"Wpt: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanMC_Wpt.getVal(),meanMC_Wpt.getError(),sigmaMC_Wpt.getVal(),sigmaMC_Wpt.getError());
    plotMC->AddTextBox(tmpName,0.20,0.85,0.60,0.88,0);
    sprintf(tmpName,"Zpt: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanMC_Zpt.getVal(),meanMC_Zpt.getError(),sigmaMC_Zpt.getVal(),sigmaMC_Zpt.getError());
    plotMC->AddTextBox(tmpName,0.20,0.82,0.60,0.85,0);
    plotMC->SetYRange(0.,1.2*(hMC_Zpt[i]->GetMaximum()));
    plotMC->Draw(myCan,kTRUE,"png");

    //==================================
    //Data: Draw Wpt and Zpt fit results
    //==================================
    sprintf(tmpName,"RD_Total_noCorr"); plotTitle = "Data Total";
    if(i==1){sprintf(tmpName,"RD_BB_noCorr"); plotTitle = "Data BB";}
    if(i==2){sprintf(tmpName,"RD_BE_noCorr"); plotTitle = "Data BE";}
    if(i==3){sprintf(tmpName,"RD_EE_noCorr"); plotTitle = "Data EE";}
    if(corrName=="CorrTotalRegion")
    {
      sprintf(tmpName,"RD_Total_Corr"); plotTitle = "Data Total";
      if(i==1){sprintf(tmpName,"RD_BB_Corr"); plotTitle = "Data BB";}
      if(i==2){sprintf(tmpName,"RD_BE_Corr"); plotTitle = "Data BE";}
      if(i==3){sprintf(tmpName,"RD_EE_Corr"); plotTitle = "Data EE";}
    }
    
    sprintf(histName,"Events / %.1f",hRD_Wpt[i]->GetBinWidth(1));
    plotRD = new CPlot(tmpName,dataframe,plotTitle,"M_{#mu#mu} [GeV]",histName);
    plotRD->setOutDir(CPlot::sOutDir);
    plotRD->SetLegend(0.65,0.7,0.88,0.82);
    h1_RD_Wpt->plotOn(dataframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    dataModel_Wpt.plotOn(dataframe,LineColor(kRed));
    h1_RD_Zpt->plotOn(dataframe,LineColor(kBlue),MarkerColor(kBlue),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    dataModel_Zpt.plotOn(dataframe,LineColor(kBlue));
    plotRD->GetLegend()->AddEntry(hRD_Wpt[i],"Wpt, Data","pl");
    plotRD->GetLegend()->AddEntry(hRD_Zpt[i],"Zpt, Data","pl");

    //plotRD->AddTextBox("Data fit",0.20,0.88,0.35,0.91,0);
    sprintf(tmpName,"Wpt: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanRD_Wpt.getVal(),meanRD_Wpt.getError(),sigmaRD_Wpt.getVal(),sigmaRD_Wpt.getError());
    plotRD->AddTextBox(tmpName,0.20,0.85,0.60,0.88,0);
    sprintf(tmpName,"Zpt: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanRD_Zpt.getVal(),meanRD_Zpt.getError(),sigmaRD_Zpt.getVal(),sigmaRD_Zpt.getError());
    plotRD->AddTextBox(tmpName,0.20,0.82,0.60,0.85,0);
    plotRD->SetYRange(0.,1.2*(hRD_Zpt[i]->GetMaximum()));
    plotRD->Draw(myCan,kTRUE,"png");
  }

  //==================================
  //Wpt histograms: Set muon regions to the x-axis
  //==================================
  h_meanMC_Wpt->GetXaxis()->SetBinLabel(1,"Total");
  h_meanMC_Wpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_meanMC_Wpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_meanMC_Wpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");
  h_widthMC_Wpt->GetXaxis()->SetBinLabel(1,"Total");
  h_widthMC_Wpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_widthMC_Wpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_widthMC_Wpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");
  
  h_meanRD_Wpt->GetXaxis()->SetBinLabel(1,"Total");
  h_meanRD_Wpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_meanRD_Wpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_meanRD_Wpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");
  h_widthRD_Wpt->GetXaxis()->SetBinLabel(1,"Total");
  h_widthRD_Wpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_widthRD_Wpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_widthRD_Wpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");

  h_meanMC_Zpt->GetXaxis()->SetBinLabel(1,"Total");
  h_meanMC_Zpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_meanMC_Zpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_meanMC_Zpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");
  h_widthMC_Zpt->GetXaxis()->SetBinLabel(1,"Total");
  h_widthMC_Zpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_widthMC_Zpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_widthMC_Zpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");
  
  h_meanRD_Zpt->GetXaxis()->SetBinLabel(1,"Total");
  h_meanRD_Zpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_meanRD_Zpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_meanRD_Zpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");
  h_widthRD_Zpt->GetXaxis()->SetBinLabel(1,"Total");
  h_widthRD_Zpt->GetXaxis()->SetBinLabel(2,"Barrel-Barrel");
  h_widthRD_Zpt->GetXaxis()->SetBinLabel(3,"Barrel-Endcap");
  h_widthRD_Zpt->GetXaxis()->SetBinLabel(4,"Endcap-Endcap");
/*
  h_meanMC_Wpt->SetLabelSize(.06);
  h_meanMC_Wpt->SetLabelSize(.06);
  h_meanMC_Wpt->SetLabelSize(.06);
  h_meanMC_Wpt->SetLabelSize(.06);
  h_widthMC_Wpt->SetLabelSize(.06);
  h_widthMC_Wpt->SetLabelSize(.06);
  h_widthMC_Wpt->SetLabelSize(.06);
  h_widthMC_Wpt->SetLabelSize(.06);
  
  h_meanRD_Wpt->SetLabelSize(.06);
  h_meanRD_Wpt->SetLabelSize(.06);
  h_meanRD_Wpt->SetLabelSize(.06);
  h_meanRD_Wpt->SetLabelSize(.06);
  h_widthRD_Wpt->SetLabelSize(.06);
  h_widthRD_Wpt->SetLabelSize(.06);
  h_widthRD_Wpt->SetLabelSize(.06);
  h_widthRD_Wpt->SetLabelSize(.06);

  h_meanMC_Zpt->SetLabelSize(.06);
  h_meanMC_Zpt->SetLabelSize(.06);
  h_meanMC_Zpt->SetLabelSize(.06);
  h_meanMC_Zpt->SetLabelSize(.06);
  h_widthMC_Zpt->SetLabelSize(.06);
  h_widthMC_Zpt->SetLabelSize(.06);
  h_widthMC_Zpt->SetLabelSize(.06);
  h_widthMC_Zpt->SetLabelSize(.06);
  
  h_meanRD_Zpt->SetLabelSize(.06);
  h_meanRD_Zpt->SetLabelSize(.06);
  h_meanRD_Zpt->SetLabelSize(.06);
  h_meanRD_Zpt->SetLabelSize(.06);
  h_widthRD_Zpt->SetLabelSize(.06);
  h_widthRD_Zpt->SetLabelSize(.06);
  h_widthRD_Zpt->SetLabelSize(.06);
  h_widthRD_Zpt->SetLabelSize(.06);
*/
  //==================================
  //Wpt: Draw Mean plots
  //==================================
  plotTitle = "Wpt: mean of M(#mu^{+}#mu^{-}), before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: mean of M(#mu^{+}#mu^{-}), after correction";
  TString histoName = "Wpt_mean_" + corrName;
  plotMeanMC_Wpt= new CPlot(histoName,plotTitle,"#mu regions","Mean of M(#mu^{+}#mu^{-})");
  plotMeanMC_Wpt->setOutDir(CPlot::sOutDir);
  plotMeanMC_Wpt->AddHist1D(h_meanRD_Wpt,"Data","E1",kBlack);
  plotMeanMC_Wpt->AddHist1D(h_meanMC_Wpt,"Z#rightarrow #mu#mu","E1",kRed);
  
  plotMeanMC_Wpt->SetLegend(0.7,0.75,0.85,0.88);
  plotMeanMC_Wpt->SetYRange(90,92);
  plotMeanMC_Wpt->Draw(myCan,kTRUE,format);

  //==================================
  //Zpt: Draw Mean plots
  //==================================
  plotTitle = "Zpt: mean of M(#mu^{+}#mu^{-}), before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Zpt: mean of M(#mu^{+}#mu^{-}), after correction";
  histoName = "Zpt_mean_" + corrName;
  plotMeanMC_Zpt= new CPlot(histoName,plotTitle,"#mu regions","Mean of M(#mu^{+}#mu^{-})");
  plotMeanMC_Zpt->setOutDir(CPlot::sOutDir);
  plotMeanMC_Zpt->AddHist1D(h_meanRD_Zpt,"Data","E1",kBlack);
  plotMeanMC_Zpt->AddHist1D(h_meanMC_Zpt,"Z#rightarrow #mu#mu","E1",kRed);
  
  plotMeanMC_Zpt->SetLegend(0.7,0.75,0.85,0.88);
  plotMeanMC_Zpt->SetYRange(90,92);
  plotMeanMC_Zpt->Draw(myCan,kTRUE,format);

  //==================================
  //Wpt: Draw Width plots
  //==================================
  plotTitle = "Wpt: width of M(#mu^{+}#mu^{-}), before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: width of M(#mu^{+}#mu^{-}), after correction";
  histoName = "Wpt_width_" + corrName;
  plotWidthMC_Wpt= new CPlot(histoName,plotTitle,"#mu regions","Width of M(#mu^{+}#mu^{-})");
  plotWidthMC_Wpt->setOutDir(CPlot::sOutDir);
  plotWidthMC_Wpt->AddHist1D(h_widthRD_Wpt,"Data","E1",kBlack);
  plotWidthMC_Wpt->AddHist1D(h_widthMC_Wpt,"Z#rightarrow #mu#mu","E1",kRed);
  
  plotWidthMC_Wpt->SetLegend(0.7,0.75,0.85,0.88);
  plotWidthMC_Wpt->SetYRange(1.6,3.4);
  plotWidthMC_Wpt->Draw(myCan,kTRUE,format);

  //==================================
  //Zpt: Draw Width plots
  //==================================
  plotTitle = "Zpt: width of M(#mu^{+}#mu^{-}), before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Zpt: width of M(#mu^{+}#mu^{-}), after correction";
  histoName = "Zpt_width_" + corrName;
  plotWidthMC_Zpt= new CPlot(histoName,plotTitle,"#mu regions","Width of M(#mu^{+}#mu^{-})");
  plotWidthMC_Zpt->setOutDir(CPlot::sOutDir);
  plotWidthMC_Zpt->AddHist1D(h_widthRD_Zpt,"Data","E1",kBlack);
  plotWidthMC_Zpt->AddHist1D(h_widthMC_Zpt,"Z#rightarrow #mu#mu","E1",kRed);
  
  plotWidthMC_Zpt->SetLegend(0.7,0.75,0.85,0.88);
  plotWidthMC_Zpt->SetYRange(1.6,3.4);
  plotWidthMC_Zpt->Draw(myCan,kTRUE,format);

  //==================================
  //MC: Draw Mean plots
  //==================================
  plotTitle = "MC: mean of M(#mu^{+}#mu^{-}), before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "MC: mean of M(#mu^{+}#mu^{-}), after correction";
  histoName = "MC_mean_" + corrName;
  plotMeanMC= new CPlot(histoName,plotTitle,"#mu regions","Mean of M(#mu^{+}#mu^{-})");
  plotMeanMC->setOutDir(CPlot::sOutDir);
  plotMeanMC->AddHist1D(h_meanMC_Wpt,"Wpt, Z#rightarrow #mu#mu","E1",kRed);
  plotMeanMC->AddHist1D(h_meanMC_Zpt,"Zpt, Z#rightarrow #mu#mu","E1",kBlue);
  
  plotMeanMC->SetLegend(0.7,0.75,0.85,0.88);
  plotMeanMC->SetYRange(90,92);
  plotMeanMC->Draw(myCan,kTRUE,format);

  //==================================
  //Data: Draw Mean plots
  //==================================
  plotTitle = "Data: mean of M(#mu^{+}#mu^{-})";
  histoName = "RD_mean_" + corrName;
  plotMeanRD= new CPlot(histoName,plotTitle,"#mu regions","Mean of M(#mu^{+}#mu^{-})");
  plotMeanRD->setOutDir(CPlot::sOutDir);
  plotMeanRD->AddHist1D(h_meanRD_Wpt,"Wpt, Data","E1",kRed);
  plotMeanRD->AddHist1D(h_meanRD_Zpt,"Zpt, Data","E1",kBlue);
  
  plotMeanRD->SetLegend(0.7,0.75,0.85,0.88);
  plotMeanRD->SetYRange(90,92);
  plotMeanRD->Draw(myCan,kTRUE,format);

  //==================================
  //MC: Draw Width plots
  //==================================
  plotTitle = "MC: width of M(#mu^{+}#mu^{-}), before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "MC: width of M(#mu^{+}#mu^{-}), after correction";
  histoName = "MC_width_" + corrName;
  plotWidthMC= new CPlot(histoName,plotTitle,"#mu regions","Width of M(#mu^{+}#mu^{-})");
  plotWidthMC->setOutDir(CPlot::sOutDir);
  plotWidthMC->AddHist1D(h_widthMC_Wpt,"Wpt, Z#rightarrow #mu#mu","E1",kRed);
  plotWidthMC->AddHist1D(h_widthMC_Zpt,"Zpt, Z#rightarrow #mu#mu","E1",kBlue);
  
  plotWidthMC->SetLegend(0.7,0.75,0.85,0.88);
  plotWidthMC->SetYRange(1.6,3.4);
  plotWidthMC->Draw(myCan,kTRUE,format);

  //==================================
  //Data: Draw Width plots
  //==================================
  plotTitle = "Data: width of M(#mu^{+}#mu^{-})";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Data: width of M(#mu^{+}#mu^{-})";
  histoName = "RD_width_" + corrName;
  plotWidthRD= new CPlot(histoName,plotTitle,"#mu regions","Width of M(#mu^{+}#mu^{-})");
  plotWidthRD->setOutDir(CPlot::sOutDir);
  plotWidthRD->AddHist1D(h_widthRD_Wpt,"Wpt, Data","E1",kRed);
  plotWidthRD->AddHist1D(h_widthRD_Zpt,"Zpt, Data","E1",kBlue);
  
  plotWidthRD->SetLegend(0.7,0.75,0.85,0.88);
  plotWidthRD->SetYRange(1.6,3.4);
  plotWidthRD->Draw(myCan,kTRUE,format);
}

TH1D *makeDiffHist(TH1D* h1, TH1D* h2, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=h1->GetNbinsX(); ibin++) {

    Double_t diff = (h1->GetBinContent(ibin)-h2->GetBinContent(ibin));

    Double_t err = sqrt(h1->GetBinContent(ibin));
    if(err==0) err= sqrt(h2->GetBinContent(ibin));

    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
  }
  hDiff->GetYaxis()->SetTitleOffset(0.42);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  return hDiff;
}
