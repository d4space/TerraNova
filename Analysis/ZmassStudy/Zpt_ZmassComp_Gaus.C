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

//Make MET difference plots
void Zpt_ZmassComp_Gaus(const TString Mode,//Channel - Muon or Electron
    const TString corrName,
    const TString outputDir 
    )
{
  TH1D* makeDiffHist(TH1D* h1, TH1D* h2, const TString name);
  const TString format("png"); 
  Int_t ratioColor = kGray+2;

  TFile *fname;

  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;

  fname = new TFile("./histo_Zmass_Corr.root");

  CPlot *plotDiff;
  CPlot *plotComp;
  CPlot *pplotEtaBin;
  CPlot *mplotEtaBin;
  CPlot *plotMllEtaP;
  CPlot *plotMllEtaM;

  TH1D *hChi;
  TH1D *hMC;
  TH1D *hRD;

  TH1D *hMCetaBinP[6];
  TH1D *hRDetaBinP[6];
  TH1D *hMCetaBinM[6];
  TH1D *hRDetaBinM[6];

  RooDataHist *h1_MCetaBinP;
  RooDataHist *h1_RDetaBinP;
  RooDataHist *h1_MCetaBinM;
  RooDataHist *h1_RDetaBinM;

  TH1D *hMCp = new TH1D("hMCp","hMCp",6,-2.1,2.1);hMCp->Sumw2();
  TH1D *hRDp = new TH1D("hRDp","hRDp",6,-2.1,2.1);hRDp->Sumw2();
  TH1D *hMCm = new TH1D("hMCm","hMCm",6,-2.1,2.1);hMCm->Sumw2();
  TH1D *hRDm = new TH1D("hRDm","hRDm",6,-2.1,2.1);hRDm->Sumw2();

  char histName[30];
  char tmpName[30];

  TCanvas *myCan;
  myCan = MakeCanvas("myCan","myCan",800,600);

  TCanvas *c;
  c = MakeCanvas("c","c",800,600);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.15);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);
  
  sprintf(histName,"h1_ZmassMC");
  if(corrName=="noCorrBB")
    sprintf(histName,"h1_ZmassMC_BB");
  if(corrName=="noCorrBE")
    sprintf(histName,"h1_ZmassMC_BE");
  if(corrName=="noCorrEE")
    sprintf(histName,"h1_ZmassMC_EE");

  if(corrName=="CorrTotalRegion")
    sprintf(histName,"h1_ZmassMCcorr");
  if(corrName=="CorrBB")
    sprintf(histName,"h1_ZmassMCcorr_BB");
  if(corrName=="CorrBE")
    sprintf(histName,"h1_ZmassMCcorr_BE");
  if(corrName=="CorrEE")
    sprintf(histName,"h1_ZmassMCcorr_EE");

  hMC=(TH1D*)fname->Get(histName)->Clone("hMC");
  hMC->Sumw2();
  cout<<"MC Bin Contents: "<<hMC->Integral()<<endl;
  hMC->SetLineStyle(1);
  hMC->SetMarkerSize(0.1);

  sprintf(histName,"h1_ZmassData");
  if(corrName=="noCorrBB")
    sprintf(histName,"h1_ZmassData_BB");
  if(corrName=="noCorrBE")
    sprintf(histName,"h1_ZmassData_BE");
  if(corrName=="noCorrEE")
    sprintf(histName,"h1_ZmassData_EE");

  if(corrName=="CorrTotalRegion")
    sprintf(histName,"h1_ZmassData");
  if(corrName=="CorrBB")
    sprintf(histName,"h1_ZmassData_BB");
  if(corrName=="CorrBE")
    sprintf(histName,"h1_ZmassData_BE");
  if(corrName=="CorrEE")
    sprintf(histName,"h1_ZmassData_EE");

  hRD=(TH1D*)fname->Get(histName)->Clone("hRD");
  hRD->Sumw2();
  cout<<"RD Bin Contents: "<<hRD->Integral()<<endl;
  hRD->SetLineStyle(1);
  hRD->SetMarkerSize(0.1);
  cout<<corrName<<"\t MC: "<<hMC->Integral()<<"\t RD: "<<hRD->Integral()<<endl;

  TString histoName = "Zmass_" + corrName;

  plotComp= new CPlot(histoName,"","","Events");
  plotComp->setOutDir(CPlot::sOutDir);

  sprintf(tmpName,"Z#rightarrow #mu#mu");
  if(Mode=="Electron")
    sprintf(tmpName,"Z#rightarrow ee");

  plotComp->AddToStack(hMC,tmpName,kOrange-2,kOrange-3);
  plotComp->AddHist1D(hRD,"Data","",kBlack);
  plotComp->SetYRange(0.,1.2*(hRD->GetMaximum()));

  plotComp->SetLegend(0.68,0.57,0.9,0.8);
  plotComp->Draw(c,kFALSE,format,1);

  hChi = makeDiffHist(hMC,hRD,"hChi");
  hChi->SetMarkerStyle(kFullCircle);
  hChi->SetMarkerSize(0.1);

  sprintf(tmpName,"Z mass, #eta-whole range");
  if (corrName=="CorrBB" || corrName=="noCorrBB")
    sprintf(tmpName,"Z mass, Barrel-Barrel");
  if (corrName=="CorrBE" || corrName=="noCorrBE")
    sprintf(tmpName,"Z mass, Barrel-Endcap");
  if (corrName=="CorrEE" || corrName=="noCorrEE")
    sprintf(tmpName,"Z mass, Endcap-Endcap");

  plotDiff=new CPlot(histoName,"",tmpName,"#chi");
  plotDiff->setOutDir(CPlot::sOutDir);
  plotDiff->AddHist1D(hChi,"EX0",ratioColor);
  plotDiff->SetYRange(-8,8);
  plotDiff->AddLine(80, 0,100, 0,kBlack,1);
  plotDiff->AddLine(80, 5,100, 5,kBlack,3);
  plotDiff->AddLine(80,-5,100,-5,kBlack,3);
  plotDiff->Draw(c,kTRUE,format,2);

  //=============================
  //Read Zmass histograms for each eta bin
  if(corrName=="noCorrTotalRegion" || corrName=="CorrTotalRegion"){
    for(int i(0);i<6;i++){
      sprintf(tmpName,"h1_ZmassMC_muP%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassMCcorr_muP%d",i);
      sprintf(histName,"hMCetaBinP_%d",i);
      hMCetaBinP[i]= (TH1D*)fname->Get(tmpName)->Clone(histName); hMCetaBinP[i]->Sumw2();
      
      sprintf(tmpName,"h1_ZmassData_muP%d",i);
      sprintf(histName,"hRDetaBinP_%d",i);
      hRDetaBinP[i]= (TH1D*)fname->Get(tmpName)->Clone(histName); hRDetaBinP[i]->Sumw2();

      sprintf(tmpName,"h1_ZmassMC_muM%d",i);
      if(corrName=="CorrTotalRegion")
	sprintf(tmpName,"h1_ZmassMCcorr_muM%d",i);
      sprintf(histName,"hMCetaBinM_%d",i);
      hMCetaBinM[i]= (TH1D*)fname->Get(tmpName)->Clone(histName); hMCetaBinM[i]->Sumw2();
      
      sprintf(tmpName,"h1_ZmassData_muM%d",i);
      sprintf(histName,"hRDetaBinM_%d",i);
      hRDetaBinM[i]= (TH1D*)fname->Get(tmpName)->Clone(histName); hRDetaBinM[i]->Sumw2();

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

      //Gauss function
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
      mcModelp.fitTo(*h1_MCetaBinP,Range("R0"));
      mcModelm.fitTo(*h1_MCetaBinM,Range("R0"));
      dataModelp.fitTo(*h1_RDetaBinP,Range("R0"));
      dataModelm.fitTo(*h1_RDetaBinM,Range("R0"));

      //==================================
      //Fill 2D plots: etaBins and mean values
      hMCp->SetBinContent(i+1,meanMCp.getVal());
      hRDp->SetBinContent(i+1,meanRDp.getVal());
      hMCm->SetBinContent(i+1,meanMCm.getVal());
      hRDm->SetBinContent(i+1,meanRDm.getVal());

      hMCp->SetBinError(i+1,meanMCp.getError());
      hRDp->SetBinError(i+1,meanRDp.getError());
      hMCm->SetBinError(i+1,meanMCm.getError());
      hRDm->SetBinError(i+1,meanRDm.getError());

      RooPlot* pframe = x.frame(Bins(40));
      RooPlot* mframe = x.frame(Bins(40));

      //==================================
      //Draw muon plus
      TCanvas *etaCanp = new TCanvas("etaCanp","Can",800,800);
      sprintf(tmpName,"pFit_etaBin%d_noCorr",i);
      if(corrName=="CorrTotalRegion")
        sprintf(tmpName,"pFit_etaBin%d_Corr",i);
      sprintf(histName,"Events / %.1f",hMCetaBinP[i]->GetBinWidth(1));
      pplotEtaBin = new CPlot(tmpName,pframe,"","M_{#mu#mu}",histName);
      pplotEtaBin->setOutDir(CPlot::sOutDir);

      pplotEtaBin->SetLegend(0.7,0.7,0.88,0.82);
      h1_MCetaBinP->plotOn(pframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
      mcModelp.plotOn(pframe,LineColor(kRed));
      h1_RDetaBinP->plotOn(pframe,LineColor(kBlack),MarkerColor(kBlack),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
      dataModelp.plotOn(pframe,LineColor(kBlack));
      pplotEtaBin->GetLegend()->AddEntry(hMCetaBinP[i],"Z #rightarrow #mu#mu","pl");
      pplotEtaBin->GetLegend()->AddEntry(hRDetaBinP[i],"Data","pl");

      pplotEtaBin->AddTextBox("#mu^{+} channel",0.18,0.88,0.35,0.91,0);
      sprintf(tmpName,"MC: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanMCp.getVal(),meanMCp.getError(),sigmaMCp.getVal(),sigmaMCp.getError());
      pplotEtaBin->AddTextBox(tmpName,0.18,0.85,0.55,0.88,0);
      sprintf(tmpName,"Data: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanRDp.getVal(),meanRDp.getError(),sigmaRDp.getVal(),sigmaRDp.getError());
      pplotEtaBin->AddTextBox(tmpName,0.18,0.82,0.55,0.85,0);

      pplotEtaBin->SetYRange(0.,1.2*(hRDetaBinP[i]->GetMaximum()));

      pplotEtaBin->Draw(etaCanp,kTRUE,"png");

      //==================================
      //Draw muon minus
      TCanvas *etaCanm = new TCanvas("etaCanm","Can",800,800);
      sprintf(tmpName,"mFit_etaBin%d_noCorr",i);
      if(corrName=="CorrTotalRegion")
        sprintf(tmpName,"mFit_etaBin%d_Corr",i);
      sprintf(histName,"Events / %.1f",hMCetaBinM[i]->GetBinWidth(1));
      mplotEtaBin = new CPlot(tmpName,mframe,"","M_{#mu#mu}",histName);
      mplotEtaBin->setOutDir(CPlot::sOutDir);

      mplotEtaBin->SetLegend(0.7,0.7,0.88,0.82);
      h1_MCetaBinM->plotOn(mframe,LineColor(kRed),MarkerColor(kRed),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
      mcModelm.plotOn(mframe,LineColor(kRed));
      h1_RDetaBinM->plotOn(mframe,LineColor(kBlack),MarkerColor(kBlack),MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"),DataError(RooAbsData::SumW2));
      dataModelm.plotOn(mframe,LineColor(kBlack));
      mplotEtaBin->GetLegend()->AddEntry(hMCetaBinM[i],"Z #rightarrow #mu#mu","pl");
      mplotEtaBin->GetLegend()->AddEntry(hRDetaBinM[i],"Data","pl");

      mplotEtaBin->AddTextBox("#mu^{-} channel",0.18,0.88,0.35,0.91,0);
      sprintf(tmpName,"MC: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanMCm.getVal(),meanMCm.getError(),sigmaMCm.getVal(),sigmaMCm.getError());
      mplotEtaBin->AddTextBox(tmpName,0.18,0.85,0.55,0.88,0);
      sprintf(tmpName,"Data: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f",meanRDm.getVal(),meanRDm.getError(),sigmaRDm.getVal(),sigmaRDm.getError());
      mplotEtaBin->AddTextBox(tmpName,0.18,0.82,0.55,0.85,0);

      mplotEtaBin->SetYRange(0.,1.2*(hRDetaBinM[i]->GetMaximum()));

      mplotEtaBin->Draw(etaCanm,kTRUE,"png");
    }

    //Draw 2D plots
    histoName = "pMassEta_" + corrName;
    plotMllEtaP= new CPlot(histoName,"","#eta of #mu^{+}","Mean of M(#mu^{+}#mu^{-})");
    plotMllEtaP->setOutDir(CPlot::sOutDir);
    plotMllEtaP->AddHist1D(hMCp,"Z#rightarrow #mu#mu","E1",kRed);
    plotMllEtaP->AddHist1D(hRDp,"Data","E1",kBlack);

    plotMllEtaP->SetLegend(0.7,0.75,0.85,0.88);
    plotMllEtaP->SetYRange(90,92);
    plotMllEtaP->Draw(myCan,kTRUE,format);

    histoName = "mMassEta_" + corrName;
    plotMllEtaM= new CPlot(histoName,"","#eta of #mu^{-}","Mean of M(#mu^{+}#mu^{-})");
    plotMllEtaM->setOutDir(CPlot::sOutDir);
    plotMllEtaM->AddHist1D(hMCm,"Z#rightarrow #mu#mu","E1",kRed);
    plotMllEtaM->AddHist1D(hRDm,"Data","E1",kBlack);

    plotMllEtaM->SetLegend(0.7,0.75,0.85,0.88);
    plotMllEtaM->SetYRange(90,92);
    plotMllEtaM->Draw(myCan,kTRUE,format);
  }
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
