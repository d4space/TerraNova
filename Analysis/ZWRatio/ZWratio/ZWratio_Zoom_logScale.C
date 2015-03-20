#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
#include "./Utils/const.h"
using namespace std;

int ZWratio_Zoom_logScale()
{
  const int nBins = 13;
  //double WptBins[nBins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
  double WptBins[nBins] = {1.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};


  ///Data
  TFile *f_WinclMu_RD = new TFile("../WptIncl_NormDiffXsec_InFid/Wpt_NormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  TFile *f_ZinclMu_RD = new TFile("../Zpt_RealData_12Bin/ZptXsecErrors/ZptXsecErrors_FidVolume.root");
  TH1D* hWpt_RD = f_WinclMu_RD->Get("hData_Xsec_BornLogScaleNorm");
  TH1D* hZpt_RD = f_ZinclMu_RD->Get("hZptDiffXsec12BinInFidNorm");
  TH1D* hWZratio_RD = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  ///Powheg
  TFile *f_WinclMu_Powheg = new TFile("../Wpt_PowhegPreFSR_12Bin/root/InclWToMuNu_PowhegWpt.root");
  TFile *f_ZinclMu_Powheg = new TFile("../Zpt_PowhegPreFSR_12Bin/root/ZToMuMu_Powheg.root");
  TH1D* hWpt_Powheg = f_WinclMu_Powheg->Get("hxsec_NormDiff");
  TH1D* hZpt_Powheg = f_ZinclMu_Powheg->Get("hxsec_NormDiff");
  TH1D* hWZratio_Powheg = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  ///ResBos
  //TFile *f_WinclMu_Resbos = new TFile("/u/user/sangilpark/TheoryNormalization/root/WpTincl_Resbos.root");
  //TFile *f_ZinclMu_Resbos = new TFile("/u/user/sangilpark/TheoryNormalization/root/ZpT_Resbos.root");
  //TH1D* hWpt_Resbos = f_WinclMu_Resbos->Get("NormDiffXsec_Resbos_12bin");
  //TH1D* hZpt_Resbos = f_ZinclMu_Resbos->Get("NormDiffXsec_Resbos_12bin");
  //TH1D* hWZratio_Resbos = new TH1D("W/Z ratio", "", nBins-1,WptBins);

  ///FEWZ
  TFile *f_ZinclMu_FEWZ = new TFile("../ZpT_FEWZ_12Bin/root/ZToMuMu_FEWZ.root");
  TFile *f_WinclMu_FEWZ = new TFile("../WpT_FEWZ_12Bin/root/WinclToMuNu_FEWZ.root");
  TH1D* hZpt_FEWZ = f_ZinclMu_FEWZ->Get("hxsec_NormDiff");
  TH1D* hWpt_FEWZ = f_WinclMu_FEWZ->Get("hxsec_NormDiff");
  TH1D* hZpt_Up_FEWZ = f_ZinclMu_FEWZ->Get("hxsec_NormDiff_up");
  TH1D* hZpt_Down_FEWZ = f_ZinclMu_FEWZ->Get("hxsec_NormDiff_down");
  TH1D* hWpt_Up_FEWZ = f_WinclMu_FEWZ->Get("hxsec_NormDiff_up");
  TH1D* hWpt_Down_FEWZ = f_WinclMu_FEWZ->Get("hxsec_NormDiff_down");
  TH1D* hWZratio_FEWZ = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  double Wpt_RD[12]={0};
  double Zpt_RD[12]={0};
  double WptErr_RD[12]={0};
  double ZptErr_RD[12]={0};
 
  double Wpt_Powheg[12]={0};
  double Zpt_Powheg[12]={0};
  double WptErr_Powheg[12]={0};
  double ZptErr_Powheg[12]={0};
 
  double Wpt_Resbos[12]={0};
  double Zpt_Resbos[12]={0};
  double WptErr_Resbos[12]={0};
  double ZptErr_Resbos[12]={0};
 
  double Wpt_FEWZ[12]={0};
  double Zpt_FEWZ[12]={0};
  double WptErr_FEWZ[12]={0};
  double ZptErr_FEWZ[12]={0};
 
  double WZratio_RD[12]={0};
  double WZratioErr_RD[12]={0};

  double WZratio_Powheg[12]={0};
  double WZratioErr_Powheg[12]={0};
  
  double WZratio_Resbos[12]={0};
  double WZratioErr_Resbos[12]={0};
  
  double WZratio_FEWZ[12]={0};
  double WZratio_Up_FEWZ[12]={0};
  double WZratio_Down_FEWZ[12]={0};
  double WZratioErr_FEWZ[12]={0};
  double WZratioScaleErr_FEWZ[12]={0};
  double WZratioTotalErr_FEWZ[12]={0};
  
  cout << fixed << setprecision(8) << endl;
  cout << " ===== Wpt and Zpt Normalilzed Differntial cross-section and errors In Fiducial Volume ==="<< endl;
  cout << "Bin"<<
    "\t Wpt  " << " \t\tError " << 
    "\t\t\t Zpt  " << " \t\tError " << endl;
  for(int i=0;i<12;i++)
  {
    Wpt_RD[i] = hWpt_RD->GetBinContent(i+1); 
    Zpt_RD[i] = hZpt_RD->GetBinContent(i+1); 
    WptErr_RD[i] = hWpt_RD->GetBinError(i+1); 
    ZptErr_RD[i] = hZpt_RD->GetBinError(i+1); 
    
//    cout <<i<< "\t" <<Wpt_RD[i] << "\t" <<  WptErr_RD[i] << "\t\t" << Zpt_RD[i] << "\t" << ZptErr_RD[i] << endl;
  
    // Data ratio and error propagation
    //WZratio_RD[i] = Wpt_RD[i] / Zpt_RD[i] ; 
    WZratio_RD[i] = Zpt_RD[i] / Wpt_RD[i] ; 
    WZratioErr_RD[i] = WZratio_RD[i] * sqrt(WptErr_RD[i]*WptErr_RD[i]/Wpt_RD[i]/Wpt_RD[i] + ZptErr_RD[i]*ZptErr_RD[i]/Zpt_RD[i]/Zpt_RD[i]); 
   
    ///Powheg,FEWZ,Resbos ratio
    WZratio_Powheg[i] = hZpt_Powheg->GetBinContent(i+1) / hWpt_Powheg->GetBinContent(i+1) ; 
    //WZratio_Resbos[i] = hZpt_Resbos->GetBinContent(i+1) / hWpt_Resbos->GetBinContent(i+1) ; 
    WZratio_FEWZ[i]   = hZpt_FEWZ->GetBinContent(i+1)  / hWpt_FEWZ->GetBinContent(i+1) ; 
    WZratio_Up_FEWZ[i]   = hZpt_Up_FEWZ->GetBinContent(i+1)  / hWpt_Up_FEWZ->GetBinContent(i+1) ; 
    WZratio_Down_FEWZ[i]   = hZpt_Down_FEWZ->GetBinContent(i+1)  / hWpt_Down_FEWZ->GetBinContent(i+1) ; 

    WZratioScaleErr_FEWZ[i] = TMath::Max(fabs(WZratio_Up_FEWZ[i]-WZratio_FEWZ[i]),fabs(WZratio_Down_FEWZ[i]-WZratio_FEWZ[i]));

    // Powheg ratio error propagation
    Wpt_Powheg[i] = hWpt_Powheg->GetBinContent(i+1);
    WptErr_Powheg[i] = hWpt_Powheg->GetBinError(i+1);

    Zpt_Powheg[i] = hZpt_Powheg->GetBinContent(i+1);
    ZptErr_Powheg[i] = hZpt_Powheg->GetBinError(i+1);

    WZratioErr_Powheg[i] = WZratio_Powheg[i] * TMath::Sqrt((ZptErr_Powheg[i]*ZptErr_Powheg[i]/Zpt_Powheg[i]/Zpt_Powheg[i] + WptErr_Powheg[i]*WptErr_Powheg[i]/Wpt_Powheg[i]/Wpt_Powheg[i]));

    /*
    // Resbos ratio error propagation
    Wpt_Resbos[i] = hWpt_Resbos->GetBinContent(i+1);
    WptErr_Resbos[i] = hWpt_Resbos->GetBinError(i+1);

    Zpt_Resbos[i] = hZpt_Resbos->GetBinContent(i+1);
    ZptErr_Resbos[i] = hZpt_Resbos->GetBinError(i+1);

    WZratioErr_Resbos[i] = WZratio_Resbos[i] * TMath::Sqrt((ZptErr_Resbos[i]*ZptErr_Resbos[i]/Zpt_Resbos[i]/Zpt_Resbos[i] + WptErr_Resbos[i]*WptErr_Resbos[i]/Wpt_Resbos[i]/Wpt_Resbos[i]));
*/
    // FEWZ ratio error propagation
    Wpt_FEWZ[i] = hWpt_FEWZ->GetBinContent(i+1);
    WptErr_FEWZ[i] = hWpt_FEWZ->GetBinError(i+1);

    Zpt_FEWZ[i] = hZpt_FEWZ->GetBinContent(i+1);
    ZptErr_FEWZ[i] = hZpt_FEWZ->GetBinError(i+1);

    WZratioErr_FEWZ[i] = WZratio_FEWZ[i] * TMath::Sqrt((ZptErr_FEWZ[i]*ZptErr_FEWZ[i]/Zpt_FEWZ[i]/Zpt_FEWZ[i] + WptErr_FEWZ[i]*WptErr_FEWZ[i]/Wpt_FEWZ[i]/Wpt_FEWZ[i]));

    WZratioTotalErr_FEWZ[i] = sqrt(WZratioErr_FEWZ[i]*WZratioErr_FEWZ[i] + WZratioScaleErr_FEWZ[i]*WZratioScaleErr_FEWZ[i]); 
    
    cout << Form("Wpt_FEWZ : %.8f +- %.8f \t Zpt_FEWZ : %.8f +- %.8f",Wpt_FEWZ[i],WptErr_FEWZ[i],Zpt_FEWZ[i],ZptErr_FEWZ[i]) << endl;
  }

  cout<<"Wpt FEWZ\t"<<"Zpt FEWZ"<<endl; 
  for(int i=0;i<12;i++)
  {
   cout<<hWpt_FEWZ->GetBinContent(i+1) <<"\t"<<hZpt_FEWZ->GetBinContent(i+1)<<endl ; 
  }
  
  cout << fixed << setprecision(3) << endl;
  cout << " ============ Wpt/Zpt ratio and errors(Number) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error " << endl;
  for(int i=0;i<12;i++)
  {

    cout <<i+1<< "\t" <<  WZratio_RD[i] <<  " \\pm " << WZratioErr_RD[i] <<  endl;
  } 
  
  cout << " ============ Wpt/Zpt ratio and errors(in %) " << endl;
  cout << "Bin"<<"\tRatio " << " \t    Error % " << endl;
  for(int i=0;i<12;i++)
  {

    cout <<i+1<< "\t" <<  WZratio_RD[i] <<  " \\pm " << WZratioErr_RD[i]/WZratio_RD[i]*100 <<  endl;
  } 
  
  for(int i=0;i<12;i++)
  {

    hWZratio_RD->SetBinContent(i+1,WZratio_RD[i]);
    hWZratio_RD->SetBinError(i+1,WZratioErr_RD[i]);
    
    hWZratio_Powheg->SetBinContent(i+1,WZratio_Powheg[i]);
    hWZratio_Powheg->SetBinError(i+1,WZratioErr_Powheg[i]);
  
    //hWZratio_Resbos->SetBinContent(i+1,WZratio_Resbos[i]);
    //hWZratio_Resbos->SetBinError(i+1,WZratioErr_Resbos[i]);
    
    hWZratio_FEWZ->SetBinContent(i+1,WZratio_FEWZ[i]);
    hWZratio_FEWZ->SetBinError(i+1,WZratioTotalErr_FEWZ[i]);
  }
  
// Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataTotalErr = new TH1D("hRatioDataTotalErr","hRatioDataTotalErr",nBins-1,WptBins);
  TH1D *hRatioPowhegTotalErr = new TH1D("hRatioPowhegTotalErr","hRatioPowhegTotalErr",nBins-1,WptBins);
  //TH1D *hRatioResbosTotalErr = new TH1D("hRatioResbosTotalErr","hRatioResbosTotalErr",nBins-1,WptBins);
  TH1D *hRatioFEWZTotalErr = new TH1D("hRatioFEWZTotalErr","hRatioFEWZTotalErr",nBins-1,WptBins);

  for(int i(0); i<nBins-1; i++)
  {
    hRatioDataTotalErr->SetBinContent(i+1,1.);
    hRatioDataTotalErr->SetBinError(i+1,WZratioErr_RD[i] / WZratio_RD[i]);

    hRatioPowhegTotalErr->SetBinContent(i+1,WZratio_Powheg[i] / WZratio_RD[i]);
    hRatioPowhegTotalErr->SetBinError(i+1,WZratioErr_Powheg[i] /WZratio_RD[i]);

    //hRatioResbosTotalErr->SetBinContent(i+1,WZratio_Resbos[i] / WZratio_RD[i]);
    //hRatioResbosTotalErr->SetBinError(i+1,WZratioErr_Resbos[i] / WZratio_RD[i]);

    hRatioFEWZTotalErr->SetBinContent(i+1,WZratio_FEWZ[i] / WZratio_RD[i]);
    hRatioFEWZTotalErr->SetBinError(i+1,WZratioTotalErr_FEWZ[i] / WZratio_RD[i]);
  }

  TFile *f_out = new TFile("./RatioZW.root","recreate");
  f_out->cd();
  hWZratio_RD->Write();
 
  // Draw
  TGraphErrors *tgWZratio_Powheg = new TGraphErrors(hWZratio_Powheg);
  //TGraphErrors *tgWZratio_Resbos = new TGraphErrors(hWZratio_Resbos);
  TGraphErrors *tgWZratio_FEWZ = new TGraphErrors(hWZratio_FEWZ);
 
  TGraphErrors* tgRatioData = new TGraphErrors(hRatioDataTotalErr);
  TGraphErrors* tgRatioPowheg = new TGraphErrors(hRatioPowhegTotalErr);
  //TGraphErrors* tgRatioResbos = new TGraphErrors(hRatioResbosTotalErr);
  TGraphErrors* tgRatioFEWZ = new TGraphErrors(hRatioFEWZTotalErr);
 
  TLegend *L1 = new TLegend(0.2,0.85,0.4,0.7);L1->SetFillColor(0); L1->SetBorderSize(0);
  L1->AddEntry(hWZratio_RD,"Data","PL");
  L1->AddEntry(tgWZratio_Powheg,"Powheg","f");
  //L1->AddEntry(tgWZratio_Resbos,"ResBos","f");
  L1->AddEntry(tgWZratio_FEWZ,"FEWZ","f");
  
  TCanvas *C1 = new TCanvas("can1","can1",800,900);
  C1->Divide(1,4,0,0);
  
  C1->cd(1)->SetPad(0,0.50,1.0,1.0);
  C1->cd(1)->SetTopMargin(0.1);
  C1->cd(1)->SetBottomMargin(0.01);
  C1->cd(1)->SetLeftMargin(0.15);
  C1->cd(1)->SetRightMargin(0.07);
  C1->cd(1)->SetTickx(1);
  C1->cd(1)->SetTicky(1);
  C1->cd(1)->SetLogx(1);


  TPaveText *cmspre = new TPaveText(0.35,0.96,0.8,0.96,"NDC");
  cmspre->SetBorderSize(0);
  cmspre->SetFillStyle(0);
  cmspre->SetTextSize(0.04);
  cmspre->AddText("CMS Preliminary, 18.4 pb^{-1} at #sqrt{s} = 8 TeV");
  
  TPaveText *Wchannel = new TPaveText(0.2,0.65,0.4,0.65,"NDC");
  Wchannel->SetBorderSize(0);
  Wchannel->SetFillStyle(0);
  Wchannel->SetTextSize(0.04);
  Wchannel->SetTextColor(kBlue);
  Wchannel->AddText("W #rightarrow #mu #nu");
  
  TPaveText *Zchannel = new TPaveText(0.2,0.60,0.42,0.60,"NDC");
  Zchannel->SetBorderSize(0);
  Zchannel->SetFillStyle(0);
  Zchannel->SetTextSize(0.04);
  Zchannel->SetTextColor(kBlue);
  Zchannel->AddText("Z #rightarrow #mu^{+} #mu^{-}");
  
  TPaveText *LeptonCut = new TPaveText(0.45,0.62,0.65,0.62,"NDC");
  LeptonCut->SetBorderSize(0);
  LeptonCut->SetFillStyle(0);
  LeptonCut->SetTextSize(0.04);
  LeptonCut->SetTextColor(kBlue);
  LeptonCut->AddText("p_{T}>20 GeV, |#eta |<2.1");

  hWZratio_RD->GetYaxis()->SetRangeUser(-1.2,6.);
  hWZratio_RD->GetYaxis()->SetTitleOffset(1.2);
  hWZratio_RD->GetYaxis()->SetLabelSize(0.035);
  hWZratio_RD->GetYaxis()->SetNdivisions(405);
  hWZratio_RD->GetYaxis()->SetTitle("(#frac{1}{#sigma^{Z}} #frac{d#sigma^{Z}}{p_{T}^{Z}})/(#frac{1}{#sigma^{W}} #frac{d#sigma^{W}}{p_{T}^{W}})");
  hWZratio_RD->GetYaxis()->SetTitleSize(0.04);
  hWZratio_RD->GetXaxis()->SetTitleOffset(0.6);
  hWZratio_RD->GetXaxis()->SetTitleSize(0.1);
  hWZratio_RD->GetXaxis()->SetLabelSize(0);
  hWZratio_RD->GetXaxis()->SetTitle("");
  hWZratio_RD->SetStats(0);
  hWZratio_RD->SetMarkerStyle(20);
  hWZratio_RD->SetMarkerColor(kBlack);
  hWZratio_RD->Draw("E1");

 // Draw TGraph style
  tgWZratio_Powheg->SetFillStyle(3004);
  tgWZratio_Powheg->SetFillColor(kRed);
  tgWZratio_Powheg->SetLineColor(kRed+2);

  //tgWZratio_Resbos->SetFillStyle(3013);
  //tgWZratio_Resbos->SetFillColor(kBlue);
  //tgWZratio_Resbos->SetLineColor(kBlue+2);

  tgWZratio_FEWZ->SetFillStyle(3005);
  tgWZratio_FEWZ->SetFillColor(kGreen);
  tgWZratio_FEWZ->SetLineColor(kGreen+2);
  
  tgWZratio_FEWZ->Draw("5");
  tgWZratio_Powheg->Draw("5");
  //tgWZratio_Resbos->Draw("5");
  
  L1->Draw();
  cmspre->Draw();
  //Wchannel->Draw();
  //Zchannel->Draw();
  //LeptonCut->Draw();
  
  C1->cd(2)->SetPad(0.028,0.505,0.76,0.805);
  
  C1->cd(2)->SetLogx(1);

  TH1D* hWZratio_RD_Zoom = (TH1D*)hWZratio_RD->Clone("hWZratio_RD_Zoom");
  
  hWZratio_RD_Zoom->GetYaxis()->SetTitle("");
  hWZratio_RD_Zoom->GetXaxis()->SetRangeUser(1,150);
  hWZratio_RD_Zoom->GetYaxis()->SetRangeUser(0.5,1.7);
  hWZratio_RD_Zoom->GetYaxis()->SetLabelSize(0.05);
  hWZratio_RD_Zoom->GetXaxis()->SetTitle("");
  hWZratio_RD_Zoom->Draw(" E1");
  
  tgWZratio_Powheg->Draw("5");
  //tgWZratio_Resbos->Draw("5");
  tgWZratio_FEWZ->Draw("5");
  gPad->RedrawAxis();
  
  C1->cd(3)->SetPad(0,0.1,1,0.5);
  C1->cd(3)->SetTopMargin(0.1);
  C1->cd(3)->SetBottomMargin(0.01);
  C1->cd(3)->SetLeftMargin(0.15);
  C1->cd(3)->SetRightMargin(0.07);
  C1->cd(3)->SetTickx(1);
  C1->cd(3)->SetTicky(1);
  C1->cd(3)->SetLogx(1);

  TH1D *hRatioDummy = new TH1D("hRatioDummy","",nBins-1,WptBins);
  
  // set canvas range and XY axis title
  hRatioDummy->GetYaxis()->SetRangeUser(0.1,3);
  hRatioDummy->GetYaxis()->SetTitle("Theory / Data");
  hRatioDummy->GetXaxis()->SetTitle("p_{T}^{V} [GeV]");
  hRatioDummy->GetYaxis()->SetTitleSize(0.07);
  hRatioDummy->GetYaxis()->SetTitleOffset(0.58);
  hRatioDummy->GetYaxis()->SetLabelSize(0.05);
  hRatioDummy->GetYaxis()->SetNdivisions(605);
  hRatioDummy->GetYaxis()->CenterTitle();
  hRatioDummy->GetXaxis()->SetTitleSize(0.07);
  hRatioDummy->GetXaxis()->SetLabelSize(0.05);

  // FEWZ Ratio plot setting
  tgRatioData->SetFillColor(kBlack);
  tgRatioData->SetFillStyle(3001);

  tgRatioPowheg->SetFillColor(kRed);
  tgRatioPowheg->SetLineColor(kRed+2);
  tgRatioPowheg->SetFillStyle(3004);
  tgRatioPowheg->SetMarkerStyle(22);
  tgRatioPowheg->SetMarkerColor(kRed+2);

  //tgRatioResbos->SetFillColor(kBlue);
  //tgRatioResbos->SetFillStyle(3013);

  tgRatioFEWZ->SetFillColor(kGreen);
  tgRatioFEWZ->SetLineColor(kGreen+2);
  tgRatioFEWZ->SetFillStyle(3005);
  tgRatioFEWZ->SetMarkerStyle(21);
  tgRatioFEWZ->SetMarkerColor(kGreen+2);
 

   // Draw Ratio plot 
  hRatioDummy->Draw();
  tgRatioData->Draw("2");
  tgRatioFEWZ->Draw("5 P");  
  tgRatioPowheg->Draw("5 P");
  //tgRatioResbos->Draw("5 P");
  gPad->RedrawAxis();
 
  // Set Zoomed Ratio plot
  C1->cd(4)->SetPad(0.028,0.056,0.76,0.405);
  
  C1->cd(4)->SetLogx(1);
  
  TH1D *hRatioDummy_Zoom = (TH1D*)hRatioDummy->Clone("hRatioDummy_Zoom");
  hRatioDummy_Zoom->GetYaxis()->SetRangeUser(0.5,1.5); // Y axis range
  hRatioDummy_Zoom->GetXaxis()->SetRangeUser(1,150); // X axis range
  hRatioDummy_Zoom->GetYaxis()->SetTitle("Theory / Data");
  hRatioDummy_Zoom->GetYaxis()->SetTitleOffset(0.68);
  hRatioDummy_Zoom->GetXaxis()->SetLabelSize(0.05);
  hRatioDummy_Zoom->GetXaxis()->SetTitle("");
  
  hRatioDummy_Zoom->Draw();
  tgRatioData->Draw("2");
  tgRatioFEWZ->Draw("5 P");  
  tgRatioPowheg->Draw("5 P");
  gPad->RedrawAxis();

  C1->SaveAs("RatioNormZW_Fid_Zoom_logScale.png");

  return 0;
}
