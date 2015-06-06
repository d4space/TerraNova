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

void Wpt_ZmassCompPtBins_Gaus(const TString Mode,//Channel - Muon or Electron
    const TString corrName,
    const TString outputDir 
    )
{
  TString plotTitle;
  TString mup_ptRange[4];
  TString mum_ptRange[4];
  mup_ptRange[0] = "20 #geq p_{T}^{#mu^{+}} < 30";
  mup_ptRange[1] = "30 #geq p_{T}^{#mu^{+}} < 40";
  mup_ptRange[2] = "40 #geq p_{T}^{#mu^{+}} < 45";
  mup_ptRange[3] = "p_{T}^{#mu^{+}} > 45";
  mum_ptRange[0] = "20 #geq p_{T}^{#mu^{-}} < 30";
  mum_ptRange[1] = "30 #geq p_{T}^{#mu^{-}} < 40";
  mum_ptRange[2] = "40 #geq p_{T}^{#mu^{-}} < 45";
  mum_ptRange[3] = "p_{T}^{#mu^{-}} > 45";

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

  CPlot *plotMllPtBinP;
  CPlot *plotMllPtBinM;
  CPlot *plotMllPtmeanP;
  CPlot *plotMllPtmeanM;
  CPlot *plotMllPtwidthP;
  CPlot *plotMllPtwidthM;

  TH1D *hMCptBinP[4];
  TH1D *hRDptBinP[4];
  TH1D *hMCptBinM[4];
  TH1D *hRDptBinM[4];

  RooDataHist *h1_MCptBinP;
  RooDataHist *h1_RDptBinP;
  RooDataHist *h1_MCptBinM;
  RooDataHist *h1_RDptBinM;

  TH1D *hMCmeanp = new TH1D("hMCmeanp","hMCmeanp",4,0,4);hMCmeanp->Sumw2();
  TH1D *hRDmeanp = new TH1D("hRDmeanp","hRDmeanp",4,0,4);hRDmeanp->Sumw2();
  TH1D *hMCmeanm = new TH1D("hMCmeanm","hMCmeanm",4,0,4);hMCmeanm->Sumw2();
  TH1D *hRDmeanm = new TH1D("hRDmeanm","hRDmeanm",4,0,4);hRDmeanm->Sumw2();

  TH1D *hMCwidthp = new TH1D("hMCwidthp","hMCwidthp",4,0,4);hMCwidthp->Sumw2();
  TH1D *hRDwidthp = new TH1D("hRDwidthp","hRDwidthp",4,0,4);hRDwidthp->Sumw2();
  TH1D *hMCwidthm = new TH1D("hMCwidthm","hMCwidthm",4,0,4);hMCwidthm->Sumw2();
  TH1D *hRDwidthm = new TH1D("hRDwidthm","hRDwidthm",4,0,4);hRDwidthm->Sumw2();

  char histName[30];
  char tmpName[30];

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
  for(int i(0);i<4;i++){
    sprintf(tmpName,"h1_Zmass_muPtP_%d",i);
    if(corrName=="CorrTotalRegion")
      sprintf(tmpName,"h1_ZmassCorr_muPtP_%d",i);
    sprintf(histName,"hMCptBinP_%d",i);
    hMCptBinP[i]= (TH1D*)fname_MC->Get(tmpName)->Clone(histName); hMCptBinP[i]->Sumw2();
    sprintf(histName,"hRDptBinP_%d",i);
    hRDptBinP[i]= (TH1D*)fname_RD->Get(tmpName)->Clone(histName); hRDptBinP[i]->Sumw2();

    sprintf(tmpName,"h1_Zmass_muPtM_%d",i);
    if(corrName=="CorrTotalRegion")
      sprintf(tmpName,"h1_ZmassCorr_muPtM_%d",i);
    sprintf(histName,"hMCptBinM_%d",i);
    hMCptBinM[i]= (TH1D*)fname_MC->Get(tmpName)->Clone(histName); hMCptBinM[i]->Sumw2();
    sprintf(histName,"hRDptBinM_%d",i);
    hRDptBinM[i]= (TH1D*)fname_RD->Get(tmpName)->Clone(histName); hRDptBinM[i]->Sumw2();

    hMCptBinP[i] -> SetMarkerSize(0.9);
    hMCptBinP[i] -> SetMarkerColor(kRed);
    hMCptBinP[i] -> SetLineColor(kRed);
    hRDptBinP[i] -> SetMarkerSize(0.9);
    hRDptBinP[i] -> SetMarkerColor(kBlack);
    hRDptBinP[i] -> SetLineColor(kBlack);
    hMCptBinM[i] -> SetMarkerSize(0.9);
    hMCptBinM[i] -> SetMarkerColor(kRed);
    hMCptBinM[i] -> SetLineColor(kRed);
    hRDptBinM[i] -> SetMarkerSize(0.9);
    hRDptBinM[i] -> SetMarkerColor(kBlack);
    hRDptBinM[i] -> SetLineColor(kBlack);

    RooRealVar x("x", "x",80,100);
    x.setBins(40);
    x.setRange("R0",86,96);

    h1_MCptBinP = new RooDataHist("h1_MCptBinP","h1_MCptBinP",RooArgSet(x),hMCptBinP[i]);
    h1_MCptBinM = new RooDataHist("h1_MCptBinM","h1_MCptBinM",RooArgSet(x),hMCptBinM[i]);
    h1_RDptBinP = new RooDataHist("h1_RDptBinP","h1_RDptBinP",RooArgSet(x),hRDptBinP[i]);
    h1_RDptBinM = new RooDataHist("h1_RDptBinM","h1_RDptBinM",RooArgSet(x),hRDptBinM[i]);

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
    mcModelp.fitTo(*h1_MCptBinP,Range("R0"));
    mcModelm.fitTo(*h1_MCptBinM,Range("R0"));
    dataModelp.fitTo(*h1_RDptBinP,Range("R0"));
    dataModelm.fitTo(*h1_RDptBinM,Range("R0"));

    //==================================
    //Fill 2D plots: ptBins and mean values
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

    RooPlot* pframe = x.frame(Bins(40));
    RooPlot* mframe = x.frame(Bins(40));

    //==================================
    //Draw muon plus
    //==================================
    sprintf(tmpName,"pFit_ptBin%d_noCorr",i);
    plotTitle = "Wpt: " + mup_ptRange[i] + ", before correction";
    if(corrName=="CorrTotalRegion")
    {
      sprintf(tmpName,"pFit_ptBin%d_Corr",i);
      plotTitle = "Wpt: " + mup_ptRange[i] + ", after correction";
    }
    sprintf(histName,"Events / %.1f",hMCptBinP[i]->GetBinWidth(1));
    plotMllPtBinP = new CPlot(tmpName,pframe,plotTitle,"M_{#mu#mu} [GeV]",histName);
    plotMllPtBinP->setOutDir(CPlot::sOutDir);

    plotMllPtBinP->SetLegend(0.7,0.7,0.88,0.82);
    h1_MCptBinP->plotOn(pframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    mcModelp.plotOn(pframe,LineColor(kRed));
    h1_RDptBinP->plotOn(pframe,LineColor(kBlack),MarkerColor(kBlack),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    dataModelp.plotOn(pframe,LineColor(kBlack));
    plotMllPtBinP->GetLegend()->AddEntry(hMCptBinP[i],"Z #rightarrow #mu#mu","pl");
    plotMllPtBinP->GetLegend()->AddEntry(hRDptBinP[i],"Data","pl");

    sprintf(tmpName,"MC: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanMCp.getVal(),meanMCp.getError(),sigmaMCp.getVal(),sigmaMCp.getError());
    plotMllPtBinP->AddTextBox(tmpName,0.20,0.83,0.6,0.88,0);
    sprintf(tmpName,"Data: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanRDp.getVal(),meanRDp.getError(),sigmaRDp.getVal(),sigmaRDp.getError());
    plotMllPtBinP->AddTextBox(tmpName,0.20,0.78,0.6,0.83,0);

    plotMllPtBinP->SetYRange(0.,1.4*(hRDptBinP[i]->GetMaximum()));

    plotMllPtBinP->Draw(myCan,kTRUE,"png");

    //==================================
    //Draw muon minus
    //=============================
    sprintf(tmpName,"mFit_ptBin%d_noCorr",i);
    plotTitle = "Wpt: " + mum_ptRange[i] + ", before correction";
    if(corrName=="CorrTotalRegion")
    {
      sprintf(tmpName,"mFit_ptBin%d_Corr",i);
      plotTitle = "Wpt: " + mum_ptRange[i] + ", after correction";
    }
    sprintf(histName,"Events / %.1f",hMCptBinM[i]->GetBinWidth(1));
    plotMllPtBinM = new CPlot(tmpName,mframe,plotTitle,"M_{#mu#mu} [GeV]",histName);
    plotMllPtBinM->setOutDir(CPlot::sOutDir);

    plotMllPtBinM->SetLegend(0.7,0.7,0.88,0.82);
    h1_MCptBinM->plotOn(mframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    mcModelm.plotOn(mframe,LineColor(kRed));
    h1_RDptBinM->plotOn(mframe,LineColor(kBlack),MarkerColor(kBlack),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
    dataModelm.plotOn(mframe,LineColor(kBlack));
    plotMllPtBinM->GetLegend()->AddEntry(hMCptBinM[i],"Z #rightarrow #mu#mu","pl");
    plotMllPtBinM->GetLegend()->AddEntry(hRDptBinM[i],"Data","pl");

    sprintf(tmpName,"MC: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanMCm.getVal(),meanMCm.getError(),sigmaMCm.getVal(),sigmaMCm.getError());
    plotMllPtBinM->AddTextBox(tmpName,0.20,0.83,0.6,0.88,0);
    sprintf(tmpName,"Data: #mu=%.2f #pm %.2f, #sigma=%.2f #pm %.2f",meanRDm.getVal(),meanRDm.getError(),sigmaRDm.getVal(),sigmaRDm.getError());
    plotMllPtBinM->AddTextBox(tmpName,0.20,0.78,0.6,0.83,0);

    plotMllPtBinM->SetYRange(0.,1.4*(hRDptBinM[i]->GetMaximum()));

    plotMllPtBinM->Draw(myCan,kTRUE,"png");
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
  //Wpt histograms: Set muon pt regions to the x-axis
  //==================================
  for(int ipt(1);ipt<=4;ipt++){
    hMCmeanp  -> GetXaxis()->SetBinLabel(ipt,mup_ptRange[ipt-1]);
    hMCmeanm  -> GetXaxis()->SetBinLabel(ipt,mum_ptRange[ipt-1]);
    hMCwidthp -> GetXaxis()->SetBinLabel(ipt,mup_ptRange[ipt-1]);
    hMCwidthm -> GetXaxis()->SetBinLabel(ipt,mum_ptRange[ipt-1]);
  }

  //==================================
  //Set x axis label size
  //==================================
  hMCmeanp ->SetLabelSize(.06);
  hMCmeanm ->SetLabelSize(.06);
  hMCwidthp->SetLabelSize(.06);
  hMCwidthm->SetLabelSize(.06);

  //==================================
  //Draw 2D plots
  //==================================
  TString histoName = "plusMean_" + corrName;
  plotTitle = "Wpt: mean at p_{T}^{#mu^{+}}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: mean at p_{T}^{#mu^{+}}, after correction";
  plotMllPtmeanP= new CPlot(histoName,plotTitle,"","Mean of M(#mu^{+}#mu^{-})");
  plotMllPtmeanP->setOutDir(CPlot::sOutDir);
  plotMllPtmeanP->AddHist1D(hMCmeanp,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllPtmeanP->AddHist1D(hRDmeanp,"Data","E1",kBlack);
                
  plotMllPtmeanP->SetLegend(0.7,0.7,0.88,0.82);
  plotMllPtmeanP->SetYRange(90,92);
  plotMllPtmeanP->Draw(myCan,kTRUE,format);

  histoName = "minusMean_" + corrName;
  plotTitle = "Wpt: mean at p_{T}^{#mu^{-}}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: mean at p_{T}^{#mu^{-}}, after correction";
  plotMllPtmeanM= new CPlot(histoName,plotTitle,"","Mean of M(#mu^{+}#mu^{-})");
  plotMllPtmeanM->setOutDir(CPlot::sOutDir);
  plotMllPtmeanM->AddHist1D(hMCmeanm,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllPtmeanM->AddHist1D(hRDmeanm,"Data","E1",kBlack);
                
  plotMllPtmeanM->SetLegend(0.7,0.7,0.88,0.82);
  plotMllPtmeanM->SetYRange(90,92);
  plotMllPtmeanM->Draw(myCan,kTRUE,format);

  histoName = "plusWidth_" + corrName;
  plotTitle = "Wpt: width at p_{T}^{#mu^{+}}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: width at p_{T}^{#mu^{+}}, after correction";
  plotMllPtwidthP= new CPlot(histoName,plotTitle,"","Width of M(#mu^{+}#mu^{-})");
  plotMllPtwidthP->setOutDir(CPlot::sOutDir);
  plotMllPtwidthP->AddHist1D(hMCwidthp,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllPtwidthP->AddHist1D(hRDwidthp,"Data","E1",kBlack);
                 
  plotMllPtwidthP->SetLegend(0.7,0.7,0.88,0.82);
  plotMllPtwidthP->SetYRange(1.5,3);
  plotMllPtwidthP->Draw(myCan,kTRUE,format);

  histoName = "minusWidth_" + corrName;
  plotTitle = "Wpt: width at p_{T}^{#mu^{-}}, before correction";
  if(corrName=="CorrTotalRegion")
    plotTitle = "Wpt: width at p_{T}^{#mu^{-}}, after correction";
  plotMllPtwidthM= new CPlot(histoName,plotTitle,"","Width of M(#mu^{+}#mu^{-})");
  plotMllPtwidthM->setOutDir(CPlot::sOutDir);
  plotMllPtwidthM->AddHist1D(hMCwidthm,"Z#rightarrow #mu#mu","E1",kRed);
  plotMllPtwidthM->AddHist1D(hRDwidthm,"Data","E1",kBlack);
                 
  plotMllPtwidthM->SetLegend(0.7,0.7,0.88,0.82);
  plotMllPtwidthM->SetYRange(1.5,3);
  plotMllPtwidthM->Draw(myCan,kTRUE,format);
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
