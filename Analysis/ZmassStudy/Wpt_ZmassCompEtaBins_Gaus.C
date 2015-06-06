#include "../Utils/MyTools.hh"
#include "../Utils/CPlot.hh"
#include "../Utils/MitStyleRemix.hh"
#include "../Utils/const.h"
#include <TFile.h>

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooChi2Var.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"

void Wpt_ZmassCompEtaBins_Gaus(const TString Mode,//Channel - Muon or Electron
    const TString corrName,
    const TString outputDir 
    )
{
  TString plotTitle;
  TString mu_etaRange[6];
  mu_etaRange[0] = "-2.1 #geq eta < -1.4";
  mu_etaRange[1] = "-1.4 #geq eta < -0.7";
  mu_etaRange[2] = "-0.7 #geq eta < 0";
  mu_etaRange[3] = "0 #geq eta < 0.7";
  mu_etaRange[4] = "0.7 #geq eta < 1.4";
  mu_etaRange[5] = "1.4 #geq eta < 2.1";

  TH1D* makeDiffHist(TH1D* h1, TH1D* h2, const TString name);
  const TString format("png"); 
  Int_t ratioColor = kGray+2;

  TFile *fname_MC;
  TFile *fname_RD;

  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;

  fname_MC = new TFile("Muon2012LoPU/Muon_DYToMuMu_S8.root");
  fname_RD = new TFile("Muon2012LoPU/Muon_RD_LowPU.root");

  if(Mode=="Electron")
  {
    fname_MC = new TFile("Electron2012LoPU/Ele_DYToEE_S8.root");
    fname_RD = new TFile("Electron2012LoPU/Ele_RD_LowPU.root");
  }

  CPlot *plotMllEtaBinP;
  CPlot *plotMllEtaBinM;
  CPlot *plotMllEtameanP;
  CPlot *plotMllEtameanM;
  CPlot *plotMllEtawidthP;
  CPlot *plotMllEtawidthM;

  TH1D *hMCetaBinP[ScaleBins];
  TH1D *hRDetaBinP[ScaleBins];
  TH1D *hMCetaBinM[ScaleBins];
  TH1D *hRDetaBinM[ScaleBins];

  RooDataHist *h1_MCetaBinP;
  RooDataHist *h1_RDetaBinP;
  RooDataHist *h1_MCetaBinM;
  RooDataHist *h1_RDetaBinM;

  TH1D *hMCmeanp = new TH1D("hMCmeanp","hMCmeanp",ScaleBins,-2.1,2.1);hMCmeanp->Sumw2();
  TH1D *hRDmeanp = new TH1D("hRDmeanp","hRDmeanp",ScaleBins,-2.1,2.1);hRDmeanp->Sumw2();
  TH1D *hMCmeanm = new TH1D("hMCmeanm","hMCmeanm",ScaleBins,-2.1,2.1);hMCmeanm->Sumw2();
  TH1D *hRDmeanm = new TH1D("hRDmeanm","hRDmeanm",ScaleBins,-2.1,2.1);hRDmeanm->Sumw2();

  TH1D *hMCwidthp = new TH1D("hMCwidthp","hMCwidthp",ScaleBins,-2.1,2.1);hMCwidthp->Sumw2();
  TH1D *hRDwidthp = new TH1D("hRDwidthp","hRDwidthp",ScaleBins,-2.1,2.1);hRDwidthp->Sumw2();
  TH1D *hMCwidthm = new TH1D("hMCwidthm","hMCwidthm",ScaleBins,-2.1,2.1);hMCwidthm->Sumw2();
  TH1D *hRDwidthm = new TH1D("hRDwidthm","hRDwidthm",ScaleBins,-2.1,2.1);hRDwidthm->Sumw2();

  char histName[50];
  char tmpName[50];

  TCanvas *myCan;
  myCan = MakeCanvas("myCan","myCan",800,600);

  myCan->SetPad(0,0,1.0,1.0);
  myCan->SetTopMargin(0.11);
  myCan->SetBottomMargin(0.15);
  myCan->SetLeftMargin(0.15);  
  myCan->SetRightMargin(0.05);  
  myCan->SetTickx(1);
  myCan->SetTicky(1);  

  //=============================
  //Read Zmass histograms for each pt bin
  //=============================
  for(int i(0);i<ScaleBins;i++){
    sprintf(tmpName,"h1_Zmass_muEtaP_%d",i);
    if(corrName=="CorrTotalRegion")
      sprintf(tmpName,"h1_ZmassCorr_muEtaP_%d",i);

    if(outputDir=="Wpt_ZmassPlotsEtaBins_noOverLap_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_noOverLap_muEtaP_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_noOverLap_muEtaP_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_LeadingLept_noOverLap_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_LeadingLept_noOverLap_muEtaP_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_LeadingLept_noOverLap_muEtaP_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_LeadingLept_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_LeadingLept_muEtaP_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_LeadingLept_muEtaP_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_TrailingLept_noOverLap_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_TrailingLept_noOverLap_muEtaP_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_TrailingLept_noOverLap_muEtaP_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_TrailingLept_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_TrailingLept_muEtaP_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_TrailingLept_muEtaP_%d",i);
    }

    sprintf(histName,"hMCetaBinP_%d",i);
    hMCetaBinP[i]= (TH1D*)fname_MC->Get(tmpName)->Clone(histName); hMCetaBinP[i]->Sumw2();
    sprintf(histName,"hRDetaBinP_%d",i);
    hRDetaBinP[i]= (TH1D*)fname_RD->Get(tmpName)->Clone(histName); hRDetaBinP[i]->Sumw2();

    sprintf(tmpName,"h1_Zmass_muEtaM_%d",i);
    if(corrName=="CorrTotalRegion")
      sprintf(tmpName,"h1_ZmassCorr_muEtaM_%d",i);

    if(outputDir=="Wpt_ZmassPlotsEtaBins_noOverLap_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_noOverLap_muEtaM_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_noOverLap_muEtaM_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_LeadingLept_noOverLap_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_LeadingLept_noOverLap_muEtaM_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_LeadingLept_noOverLap_muEtaM_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_LeadingLept_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_LeadingLept_muEtaM_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_LeadingLept_muEtaM_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_TrailingLept_noOverLap_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_TrailingLept_noOverLap_muEtaM_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_TrailingLept_noOverLap_muEtaM_%d",i);
    }

    if(outputDir=="Wpt_ZmassPlotsEtaBins_TrailingLept_Gaus")
    {
      sprintf(tmpName,"h1_Zmass_TrailingLept_muEtaM_%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassCorr_TrailingLept_muEtaM_%d",i);
    }

    sprintf(histName,"hMCetaBinM_%d",i);
    hMCetaBinM[i]= (TH1D*)fname_MC->Get(tmpName)->Clone(histName); hMCetaBinM[i]->Sumw2();
    sprintf(histName,"hRDetaBinM_%d",i);
    hRDetaBinM[i]= (TH1D*)fname_RD->Get(tmpName)->Clone(histName); hRDetaBinM[i]->Sumw2();
    
    hMCetaBinP[i] -> SetMarkerSize(0.9);
    hMCetaBinP[i] -> SetMarkerColor(kRed);
    hMCetaBinP[i] -> SetLineColor(kRed);
    hRDetaBinP[i] -> SetMarkerSize(0.9);
    hRDetaBinP[i] -> SetMarkerColor(kBlack);
    hRDetaBinP[i] -> SetLineColor(kBlack);
    hMCetaBinM[i] -> SetMarkerSize(0.9);
    hMCetaBinM[i] -> SetMarkerColor(kRed);
    hMCetaBinM[i] -> SetLineColor(kRed);
    hRDetaBinM[i] -> SetMarkerSize(0.9);
    hRDetaBinM[i] -> SetMarkerColor(kBlack);
    hRDetaBinM[i] -> SetLineColor(kBlack);

    RooRealVar x("x", "x",80,100);
    x.setBins(40);
    x.setRange("R0",86,96);

    h1_MCetaBinP = new RooDataHist("h1_MCetaBinP","h1_MCetaBinP",RooArgSet(x),hMCetaBinP[i]);
    h1_MCetaBinM = new RooDataHist("h1_MCetaBinM","h1_MCetaBinM",RooArgSet(x),hMCetaBinM[i]);
    h1_RDetaBinP = new RooDataHist("h1_RDetaBinP","h1_RDetaBinP",RooArgSet(x),hRDetaBinP[i]);
    h1_RDetaBinM = new RooDataHist("h1_RDetaBinM","h1_RDetaBinM",RooArgSet(x),hRDetaBinM[i]);

    //=============================
    //Gauss function
    //=============================
    RooRealVar meanMCp("meanMCp","",91.2,80,100);
    RooRealVar meanRDp("meanRDp","",91.2,80,100);
    RooRealVar meanMCm("meanMCm","",91.2,80,100);
    RooRealVar meanRDm("meanRDm","",91.2,80,100);
    RooRealVar sigmaMCp("sigmaMCp","",5,-50.,50.);
    RooRealVar sigmaRDp("sigmaRDp","",5,-50.,50.);
    RooRealVar sigmaMCm("sigmaMCm","",5,-50.,50.);
    RooRealVar sigmaRDm("sigmaRDm","",5,-50.,50.);

    RooGaussian mcModelp("mcModelp", "",x,meanMCp,sigmaMCp);
    RooGaussian mcModelm("mcModelm", "",x,meanMCm,sigmaMCm);
    RooGaussian dataModelp("dataModelp", "",x,meanRDp,sigmaRDp);
    RooGaussian dataModelm("dataModelm", "",x,meanRDm,sigmaRDm);

    //==================================
    //Fit Zmass distribution with Gaussian function
    //=============================
    mcModelp.fitTo(*h1_MCetaBinP,Range("R0"));
    mcModelm.fitTo(*h1_MCetaBinM,Range("R0"));
    dataModelp.fitTo(*h1_RDetaBinP,Range("R0"));
    dataModelm.fitTo(*h1_RDetaBinM,Range("R0"));

    //==================================
    //Fill 2D plots: etaBins and mean values
    //=============================
    hMCmeanp->SetBinContent(i+1,meanMCp.getVal());
    hRDmeanp->SetBinContent(i+1,meanRDp.getVal());
    hMCmeanm->SetBinContent(i+1,meanMCm.getVal());
    hRDmeanm->SetBinContent(i+1,meanRDm.getVal());

    hMCmeanp->SetBinError(i+1,meanMCp.getError());
    hRDmeanp->SetBinError(i+1,meanRDp.getError());
    hMCmeanm->SetBinError(i+1,meanMCm.getError());
    hRDmeanm->SetBinError(i+1,meanRDm.getError());

    hMCwidthp->SetBinContent(i+1,sigmaMCp.getVal());
    hRDwidthp->SetBinContent(i+1,sigmaRDp.getVal());
    hMCwidthm->SetBinContent(i+1,sigmaMCm.getVal());
    hRDwidthm->SetBinContent(i+1,sigmaRDm.getVal());

    hMCwidthp->SetBinError(i+1,sigmaMCp.getError());
    hRDwidthp->SetBinError(i+1,sigmaRDp.getError());
    hMCwidthm->SetBinError(i+1,sigmaMCm.getError());
    hRDwidthm->SetBinError(i+1,sigmaRDm.getError());

    cout<<meanMCp.getVal()<<"\t"<<meanMCm.getVal()<<"\t"<<meanRDp.getVal()<<"\t"<<meanRDm.getVal()<<endl;
    RooPlot* pframe = x.frame(Bins(40));
    RooPlot* mframe = x.frame(Bins(40));

    //==================================
    //Draw muon plus
    //==================================
    sprintf(tmpName,"pFit_etaBin%d_noCorr",i);
    plotTitle = "Wpt: " + mu_etaRange[i] + " of #mu^{+}, before correction";
    if(corrName=="CorrTotalRegion")
    {
      sprintf(tmpName,"pFit_etaBin%d_Corr",i);
      plotTitle = "Wpt: " + mu_etaRange[i] + " of #mu^{+}, after correction";
    }
    sprintf(histName,"Events / %.1f",hMCetaBinP[i]->GetBinWidth(1));
    plotMllEtaBinP = new CPlot(tmpName,pframe,plotTitle,"M_{#mu#mu} [GeV]",histName);
    plotMllEtaBinP->setOutDir(CPlot::sOutDir);

    plotMllEtaBinP->SetLegend(0.7,0.7,0.88,0.82);
    h1_MCetaBinP->plotOn(pframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    mcModelp.plotOn(pframe,LineColor(kRed));
    h1_RDetaBinP->plotOn(pframe,LineColor(kBlack),MarkerColor(kBlack),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    dataModelp.plotOn(pframe,LineColor(kBlack));
    plotMllEtaBinP->GetLegend()->AddEntry(hMCetaBinP[i],"Z #rightarrow #mu#mu","pl");
    plotMllEtaBinP->GetLegend()->AddEntry(hRDetaBinP[i],"Data","pl");

    sprintf(tmpName,"MC: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanMCp.getVal(),meanMCp.getError(),sigmaMCp.getVal(),sigmaMCp.getError());
    plotMllEtaBinP->AddTextBox(tmpName,0.20,0.83,0.6,0.88,0);
    sprintf(tmpName,"Data: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanRDp.getVal(),meanRDp.getError(),sigmaRDp.getVal(),sigmaRDp.getError());
    plotMllEtaBinP->AddTextBox(tmpName,0.20,0.78,0.6,0.83,0);

    plotMllEtaBinP->SetYRange(0.,1.4*TMath::Max(hRDetaBinP[i]->GetMaximum(),hMCetaBinP[i]->GetMaximum()));

    plotMllEtaBinP->Draw(myCan,kTRUE,"png");

    //==================================
    //Draw muon minus
    //=============================
    sprintf(tmpName,"mFit_etaBin%d_noCorr",i);
    plotTitle = "Wpt: " + mu_etaRange[i] + " of #mu^{-}, before correction";
    if(corrName=="CorrTotalRegion")
    {
      sprintf(tmpName,"mFit_etaBin%d_Corr",i);
      plotTitle = "Wpt: " + mu_etaRange[i] + " of #mu^{-}, after correction";
    }
    sprintf(histName,"Events / %.1f",hMCetaBinM[i]->GetBinWidth(1));
    plotMllEtaBinM = new CPlot(tmpName,mframe,plotTitle,"M_{#mu#mu} [GeV]",histName);
    plotMllEtaBinM->setOutDir(CPlot::sOutDir);

    plotMllEtaBinM->SetLegend(0.7,0.7,0.88,0.82);
    h1_MCetaBinM->plotOn(mframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    mcModelm.plotOn(mframe,LineColor(kRed));
    h1_RDetaBinM->plotOn(mframe,LineColor(kBlack),MarkerColor(kBlack),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    dataModelm.plotOn(mframe,LineColor(kBlack));
    plotMllEtaBinM->GetLegend()->AddEntry(hMCetaBinM[i],"Z #rightarrow #mu#mu","pl");
    plotMllEtaBinM->GetLegend()->AddEntry(hRDetaBinM[i],"Data","pl");

    sprintf(tmpName,"MC: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanMCm.getVal(),meanMCm.getError(),sigmaMCm.getVal(),sigmaMCm.getError());
    plotMllEtaBinM->AddTextBox(tmpName,0.20,0.83,0.6,0.88,0);
    sprintf(tmpName,"Data: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanRDm.getVal(),meanRDm.getError(),sigmaRDm.getVal(),sigmaRDm.getError());
    plotMllEtaBinM->AddTextBox(tmpName,0.20,0.78,0.6,0.83,0);

    plotMllEtaBinM->SetYRange(0.,1.4*TMath::Max(hRDetaBinM[i]->GetMaximum(),hMCetaBinM[i]->GetMaximum()));

    plotMllEtaBinM->Draw(myCan,kTRUE,"png");
  }

  //==================================
  //Save Mean and Width histograms to root file
  //==================================
  TString filename = outputDir + "/Wpt_MeanWidth_Gaus_" + corrName + ".root";
  TFile *outfile = new TFile(filename,"RECREATE");
  hMCmeanp ->Write();
  hMCmeanm ->Write();
  hMCwidthp->Write();
  hMCwidthm->Write();
  hRDmeanp ->Write();
  hRDmeanm ->Write();
  hRDwidthp->Write();
  hRDwidthm->Write();
  outfile->Close();

  //==================================
  //Draw 2D plots
  //==================================
  TString histoName = "plusMean_" + corrName;
  plotTitle = "Wpt: mean at #eta of #mu^{+}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: mean at #eta of #mu^{+}, after correction";
  plotMllEtameanP= new CPlot(histoName,plotTitle,"#eta of #mu^{+}","Mean of M(#mu^{+}#mu^{-})");
  plotMllEtameanP->setOutDir(CPlot::sOutDir);
  plotMllEtameanP->AddHist1D(hMCmeanp,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllEtameanP->AddHist1D(hRDmeanp,"Data","E1",kBlack);
                
  plotMllEtameanP->SetLegend(0.7,0.7,0.88,0.82);
  plotMllEtameanP->SetYRange(90,92);
  plotMllEtameanP->Draw(myCan,kTRUE,format);

  histoName = "minusMean_" + corrName;
  plotTitle = "Wpt: mean at #eta of #mu^{-}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: mean at #eta of #mu^{-}, after correction";
  plotMllEtameanM= new CPlot(histoName,plotTitle,"#eta of #mu^{-}","Mean of M(#mu^{+}#mu^{-})");
  plotMllEtameanM->setOutDir(CPlot::sOutDir);
  plotMllEtameanM->AddHist1D(hMCmeanm,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllEtameanM->AddHist1D(hRDmeanm,"Data","E1",kBlack);
                
  plotMllEtameanM->SetLegend(0.7,0.7,0.88,0.82);
  plotMllEtameanM->SetYRange(90,92);
  plotMllEtameanM->Draw(myCan,kTRUE,format);

  histoName = "plusWidth_" + corrName;
  plotTitle = "Wpt: width at #eta of #mu^{+}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: width at #eta of #mu^{+}, after correction";
  plotMllEtawidthP= new CPlot(histoName,plotTitle,"#eta of #mu^{+}","Width of M(#mu^{+}#mu^{-})");
  plotMllEtawidthP->setOutDir(CPlot::sOutDir);
  plotMllEtawidthP->AddHist1D(hMCwidthp,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllEtawidthP->AddHist1D(hRDwidthp,"Data","E1",kBlack);
                 
  plotMllEtawidthP->SetLegend(0.7,0.7,0.88,0.82);
  plotMllEtawidthP->SetYRange(1.5,3);
  plotMllEtawidthP->Draw(myCan,kTRUE,format);

  histoName = "minusWidth_" + corrName;
  plotTitle = "Wpt: width at #eta of #mu^{-}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: width at #eta of #mu^{-}, after correction";
  plotMllEtawidthM= new CPlot(histoName,plotTitle,"#eta of #mu^{-}","Width of M(#mu^{+}#mu^{-})");
  plotMllEtawidthM->setOutDir(CPlot::sOutDir);
  plotMllEtawidthM->AddHist1D(hMCwidthm,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllEtawidthM->AddHist1D(hRDwidthm,"Data","E1",kBlack);
                 
  plotMllEtawidthM->SetLegend(0.7,0.7,0.88,0.82);
  plotMllEtawidthM->SetYRange(1.5,3);
  plotMllEtawidthM->Draw(myCan,kTRUE,format);
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
