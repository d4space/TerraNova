#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TH1D.h"
#include <TGraphErrors.h>             // graph class
#include <TGraphAsymmErrors.h>        // graph class
#include <TLatex.h>
using namespace std;

void FEWZ_PDFUncer(double *ZWRatio, double *Err);
void Powheg_PDFUncer(double *ZWRatio, double *Err);

int ZWratio_Zoom_logScale()
{
  const int nBins = 13;
  //double WptBins[nBins] = {0.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};
  double WptBins[nBins] = {1.0,7.5,12.5,17.5,30,40,50,70,110,150,190,250,600};


  ///Data
  TFile *f_WinclMu_RD = new TFile("../WptIncl_NormDiffXsec_InFid/Wpt_NormDiffXsec_InFid_RDResBosPowhegFEWZ.root");
  TFile *f_ZinclMu_RD = new TFile("../Zpt_RealData_12Bin/ZptXsecErrors/ZptXsecErrors_FidVolume.root");
  TH1D* hWpt_RD = (TH1D*)f_WinclMu_RD->Get("hData_Xsec_BornLogScaleNorm");
  TH1D* hZpt_RD = (TH1D*)f_ZinclMu_RD->Get("hZptDiffXsec12BinInFidNorm");
  TH1D* hWZratio_RD = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  ///Powheg
  TFile *f_WinclMu_Powheg = new TFile("../Wpt_PowhegPreFSR_12Bin/root/InclWToMuNu_PowhegWpt.root");
  TFile *f_ZinclMu_Powheg = new TFile("../Zpt_PowhegPreFSR_12Bin/root/ZToMuMu_Powheg.root");
  TH1D* hWpt_Powheg = (TH1D*)f_WinclMu_Powheg->Get("hxsec_NormDiff");
  TH1D* hZpt_Powheg = (TH1D*)f_ZinclMu_Powheg->Get("hxsec_NormDiff");
  TH1D* hWZratio_Powheg = new TH1D("W/Z ratio", "", nBins-1,WptBins);
 
  ///ResBos
  TFile *f_WinclMu_Resbos = new TFile("../WpT_Resbos_12Bin/root/WpTincl_Resbos.root");
  TFile *f_ZinclMu_Resbos = new TFile("../ZpT_Resbos_12Bin/root/ZToMuMu_Resbos.root");
  TH1D* hWpt_Resbos = f_WinclMu_Resbos->Get("NormDiffXsec_Resbos_12bin");
  TH1D* hZpt_Resbos = f_ZinclMu_Resbos->Get("NormDiffXsec");
  TH1D* hWZratio_Resbos = new TH1D("W/Z ratio", "", nBins-1,WptBins); 

  ///FEWZ
  TFile *f_ZinclMu_FEWZ = new TFile("../ZpT_FEWZ_12Bin/root/ZToMuMu_FEWZ.root");
  TFile *f_WinclMu_FEWZ = new TFile("../WpT_FEWZ_12Bin/root/WinclToMuNu_FEWZ.root");
  TH1D* hZpt_FEWZ = (TH1D*)f_ZinclMu_FEWZ->Get("hxsec_NormDiff");
  TH1D* hWpt_FEWZ = (TH1D*)f_WinclMu_FEWZ->Get("hxsec_NormDiff");
  TH1D* hZpt_Up_FEWZ = (TH1D*)f_ZinclMu_FEWZ->Get("hxsec_NormDiff_up");
  TH1D* hZpt_Down_FEWZ = (TH1D*)f_ZinclMu_FEWZ->Get("hxsec_NormDiff_down");
  TH1D* hWpt_Up_FEWZ = (TH1D*)f_WinclMu_FEWZ->Get("hxsec_NormDiff_up");
  TH1D* hWpt_Down_FEWZ = (TH1D*)f_WinclMu_FEWZ->Get("hxsec_NormDiff_down");
  TH1D* hWZratio_FEWZ = new TH1D("W/Z ratio", "", nBins-1,WptBins);
  
  double Wpt_RD[12]={0};
  double Zpt_RD[12]={0};
  double WptErr_RD[12]={0};
  double ZptErr_RD[12]={0};
 
  double Wpt_Powheg[12]={0};
  double Zpt_Powheg[12]={0};
  double WptStatErr_Powheg[12]={0};
  double ZptStatErr_Powheg[12]={0};
 
  double Wpt_Resbos[12]={0};
  double Zpt_Resbos[12]={0};
  double WptErr_Resbos[12]={0};
  double ZptErr_Resbos[12]={0};
 
  double Wpt_FEWZ[12]={0};
  double Zpt_FEWZ[12]={0};
  double WptStatErr_FEWZ[12]={0};
  double ZptStatErr_FEWZ[12]={0};
 
  double WZratio_RD[12]={0};
  double WZratioErr_RD[12]={0};

  double WZratio_Powheg[12]={0};
  double WZratioStatErr_Powheg[12]={0};
  double WZratioPDFErr_Powheg[12]={0};
  double WZratioTotalErr_Powheg[12]={0}; // sqrt(Stat^2 + PDF^2 + Scale^2)
  
  double WZratio_Resbos[12]={0};
  double WZratioErr_Resbos[12]={0};
  
  double WZratio_FEWZ[12]={0};
  double WZratio_Up_FEWZ[12]={0};
  double WZratio_Down_FEWZ[12]={0};
  double WZratioStatErr_FEWZ[12]={0}; // Error propagation
  double WZratioPDFErr_FEWZ[12]={0}; // Norm PDF %
  double WZratioScaleErr_FEWZ[12]={0}; // Norm Scale subtraction
  double WZratioTotalErr_FEWZ[12]={0}; // sqrt(Stat^2 + PDF^2 + Scale^2)
  
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
    WZratio_Resbos[i] = hZpt_Resbos->GetBinContent(i+1) / hWpt_Resbos->GetBinContent(i+1) ; 
    WZratio_FEWZ[i]   = hZpt_FEWZ->GetBinContent(i+1)  / hWpt_FEWZ->GetBinContent(i+1) ; 
    WZratio_Up_FEWZ[i]   = hZpt_Up_FEWZ->GetBinContent(i+1)  / hWpt_Up_FEWZ->GetBinContent(i+1) ; 
    WZratio_Down_FEWZ[i]   = hZpt_Down_FEWZ->GetBinContent(i+1)  / hWpt_Down_FEWZ->GetBinContent(i+1) ; 

    WZratioScaleErr_FEWZ[i] = TMath::Max(fabs(WZratio_Up_FEWZ[i]-WZratio_FEWZ[i]),fabs(WZratio_Down_FEWZ[i]-WZratio_FEWZ[i]));

    // Powheg ratio error propagation
    Wpt_Powheg[i] = hWpt_Powheg->GetBinContent(i+1);
    WptStatErr_Powheg[i] = hWpt_Powheg->GetBinError(i+1);

    Zpt_Powheg[i] = hZpt_Powheg->GetBinContent(i+1);
    ZptStatErr_Powheg[i] = hZpt_Powheg->GetBinError(i+1);

    WZratioStatErr_Powheg[i] = WZratio_Powheg[i] * TMath::Sqrt((ZptStatErr_Powheg[i]*ZptStatErr_Powheg[i]/Zpt_Powheg[i]/Zpt_Powheg[i] + WptStatErr_Powheg[i]*WptStatErr_Powheg[i]/Wpt_Powheg[i]/Wpt_Powheg[i]));

    Powheg_PDFUncer(WZratio_Powheg,WZratioPDFErr_Powheg);
    
    WZratioTotalErr_Powheg[i] = sqrt(WZratioStatErr_Powheg[i]*WZratioStatErr_Powheg[i] + WZratioPDFErr_Powheg[i]*WZratioPDFErr_Powheg[i]); //  
    cout << Form("ZW Powheg Ratio : %.4f \t Stat : %.4f \t PDF : %.4f ",WZratio_Powheg[i],WZratioStatErr_Powheg[i],WZratioPDFErr_Powheg[i]) << endl;
    
    //*
    // Resbos ratio error propagation
    Wpt_Resbos[i] = hWpt_Resbos->GetBinContent(i+1);
    WptErr_Resbos[i] = hWpt_Resbos->GetBinError(i+1);

    Zpt_Resbos[i] = hZpt_Resbos->GetBinContent(i+1);
    ZptErr_Resbos[i] = hZpt_Resbos->GetBinError(i+1);

    WZratioErr_Resbos[i] = WZratio_Resbos[i] * TMath::Sqrt((ZptErr_Resbos[i]*ZptErr_Resbos[i]/Zpt_Resbos[i]/Zpt_Resbos[i] + WptErr_Resbos[i]*WptErr_Resbos[i]/Wpt_Resbos[i]/Wpt_Resbos[i]));
//*/
    // FEWZ ratio error propagation
    Wpt_FEWZ[i] = hWpt_FEWZ->GetBinContent(i+1);
    WptStatErr_FEWZ[i] = hWpt_FEWZ->GetBinError(i+1);

    Zpt_FEWZ[i] = hZpt_FEWZ->GetBinContent(i+1);
    ZptStatErr_FEWZ[i] = hZpt_FEWZ->GetBinError(i+1);

    WZratioStatErr_FEWZ[i] = WZratio_FEWZ[i] * TMath::Sqrt((ZptStatErr_FEWZ[i]*ZptStatErr_FEWZ[i]/Zpt_FEWZ[i]/Zpt_FEWZ[i] + WptStatErr_FEWZ[i]*WptStatErr_FEWZ[i]/Wpt_FEWZ[i]/Wpt_FEWZ[i]));

    FEWZ_PDFUncer(WZratio_FEWZ,WZratioPDFErr_FEWZ);
    //cout << "WZratioPDFErr_FEWZ in main : " << WZratioPDFErr_FEWZ[i] << endl; // check PDF error number is correctly returned

    WZratioTotalErr_FEWZ[i] = sqrt(WZratioStatErr_FEWZ[i]*WZratioStatErr_FEWZ[i] + WZratioPDFErr_FEWZ[i]*WZratioPDFErr_FEWZ[i] + WZratioScaleErr_FEWZ[i]*WZratioScaleErr_FEWZ[i]); //  
    
    //cout << Form("Wpt_FEWZ : %.8f +- %.8f \t Zpt_FEWZ : %.8f +- %.8f",Wpt_FEWZ[i],WptStatErr_FEWZ[i],Zpt_FEWZ[i],ZptStatErr_FEWZ[i]) << endl;
  }


  //// Print Powheg errors
  cout <<"=====  Print Powheg errors "<< endl;
  cout << fixed << setprecision(2) << endl;
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioStatErr_Powheg: " << WZratioStatErr_Powheg[i] <<"\t in %\t"<< 100*WZratioStatErr_Powheg[i]/WZratio_Powheg[i]<<endl; 
  }
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioPDFErr_Powheg : " << WZratioPDFErr_Powheg[i] << "\t in %\t"<< 100*WZratioPDFErr_Powheg[i]/WZratio_Powheg[i]<<endl;
  }
 //// Print FEWZ errors
  cout <<"=====  Print FEWZ errors "<< endl;

  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioStatErr_FEWZ : " << WZratioStatErr_FEWZ[i] << "\t in %\t"<< 100*WZratioStatErr_FEWZ[i]/WZratio_FEWZ[i]<<endl;
  }

  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioPDFErr_FEWZ : " << WZratioPDFErr_FEWZ[i] << "\t in %\t"<< 100*WZratioPDFErr_FEWZ[i]/WZratio_FEWZ[i]<<endl;
  }
  for(int i=0;i<12;i++)
  {
    cout << i<<"\tWZratioScaleErr_FEWZ : " << WZratioScaleErr_FEWZ[i] << "\t in %\t"<< 100*WZratioScaleErr_FEWZ[i]/WZratio_FEWZ[i]<<endl;
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
    hWZratio_Powheg->SetBinError(i+1,WZratioTotalErr_Powheg[i]);
  
    hWZratio_Resbos->SetBinContent(i+1,WZratio_Resbos[i]);
    //hWZratio_Resbos->SetBinError(i+1,WZratioErr_Resbos[i]);
    hWZratio_Resbos->SetBinError(i+1,0);
    
    hWZratio_FEWZ->SetBinContent(i+1,WZratio_FEWZ[i]);
    hWZratio_FEWZ->SetBinError(i+1,WZratioTotalErr_FEWZ[i]);
  }
  
// Theory/Data ratio plot errors related to Real Data
  TH1D *hRatioDataTotalErr = new TH1D("hRatioDataTotalErr","hRatioDataTotalErr",nBins-1,WptBins);
  TH1D *hRatioPowhegTotalErr = new TH1D("hRatioPowhegTotalErr","hRatioPowhegTotalErr",nBins-1,WptBins);
  TH1D *hRatioResbosTotalErr = new TH1D("hRatioResbosTotalErr","hRatioResbosTotalErr",nBins-1,WptBins);
  TH1D *hRatioFEWZTotalErr = new TH1D("hRatioFEWZTotalErr","hRatioFEWZTotalErr",nBins-1,WptBins);

  for(int i(0); i<nBins-1; i++)
  {
    hRatioDataTotalErr->SetBinContent(i+1,1.);
    hRatioDataTotalErr->SetBinError(i+1,WZratioErr_RD[i] / WZratio_RD[i]);

    hRatioPowhegTotalErr->SetBinContent(i+1,WZratio_Powheg[i] / WZratio_RD[i]);
    hRatioPowhegTotalErr->SetBinError(i+1,WZratioTotalErr_Powheg[i] /WZratio_RD[i]);

    hRatioResbosTotalErr->SetBinContent(i+1,WZratio_Resbos[i] / WZratio_RD[i]);
    //hRatioResbosTotalErr->SetBinError(i+1,WZratioErr_Resbos[i] / WZratio_RD[i]);
    hRatioResbosTotalErr->SetBinError(i+1,0);

    hRatioFEWZTotalErr->SetBinContent(i+1,WZratio_FEWZ[i] / WZratio_RD[i]);
    hRatioFEWZTotalErr->SetBinError(i+1,WZratioTotalErr_FEWZ[i] / WZratio_RD[i]);
  }

  TFile *f_out = new TFile("./RatioZW.root","recreate");
  f_out->cd();
  hWZratio_RD->Write();


  //Color Transparent
  TColor *colRed = gROOT->GetColor(kRed);
  TColor *colBlue = gROOT->GetColor(kBlue);
  TColor *colGreen = gROOT->GetColor(kGreen);
  colRed->SetAlpha(0.2);
  colBlue->SetAlpha(0.2);
  colGreen->SetAlpha(0.2);
  gStyle->SetOptStat(0); 

  // Draw
  TGraphErrors *tgWZratio_Powheg = new TGraphErrors(hWZratio_Powheg);
  TGraphErrors *tgWZratio_Resbos = new TGraphErrors(hWZratio_Resbos);
  TGraphErrors *tgWZratio_FEWZ = new TGraphErrors(hWZratio_FEWZ);
 
  TGraphErrors* tgRatioData = new TGraphErrors(hRatioDataTotalErr);
  TGraphErrors* tgRatioPowheg = new TGraphErrors(hRatioPowhegTotalErr);
  TGraphErrors* tgRatioResbos = new TGraphErrors(hRatioResbosTotalErr);
  TGraphErrors* tgRatioFEWZ = new TGraphErrors(hRatioFEWZTotalErr);
 
  TLegend *L1 = new TLegend(0.2,0.65,0.5,0.85);L1->SetFillColor(0); L1->SetBorderSize(0);
  L1->AddEntry(hWZratio_RD,"data","PL");
  L1->AddEntry(tgWZratio_Powheg,"Powheg","f");
  L1->AddEntry(tgWZratio_Resbos,"ResBos","f");
  L1->AddEntry(tgWZratio_FEWZ,"FEWZ","f");
  
  TCanvas *C1 = new TCanvas("can1","can1",800,900);
  C1->Divide(1,4,0,0);
  
  C1->cd(1)->SetPad(0,0.50,0.95,0.95);
  C1->cd(1)->SetTopMargin(0.1);
  C1->cd(1)->SetBottomMargin(0.01);
  C1->cd(1)->SetLeftMargin(0.15);
  C1->cd(1)->SetRightMargin(0.07);
  C1->cd(1)->SetTickx(1);
  //C1->cd(1)->SetTicky(1);
  C1->cd(1)->SetTicky(2);
  C1->cd(1)->SetLogx(1);


  //TPaveText *cmspre = new TPaveText(0.35,0.96,0.8,0.96,"NDC");
  TPaveText *cmspre = new TPaveText(0.70,0.93,0.95,0.95,"NDC");
  cmspre->SetBorderSize(0);
  cmspre->SetFillStyle(0);
  cmspre->SetTextSize(0.06);
  cmspre->AddText("CMS Preliminary");
  
  TPaveText *lumi = new TPaveText(0.35,0.83,0.8,0.85,"NDC");
  lumi->SetBorderSize(0);
  lumi->SetFillStyle(0);
  lumi->SetTextSize(0.05);
  lumi->AddText("L = 18.4 pb^{-1}, #sqrt{s} = 8 TeV");
  
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
  //hWZratio_RD->GetYaxis()->SetTitleOffset(1.2);
  hWZratio_RD->GetYaxis()->SetTitleOffset(0.8);
  //hWZratio_RD->GetYaxis()->SetLabelSize(0.035);
  hWZratio_RD->GetYaxis()->SetLabelSize(0.033);
  hWZratio_RD->GetYaxis()->SetTitleSize(0.053);
  //hWZratio_RD->GetYaxis()->SetNdivisions(405);
  hWZratio_RD->GetYaxis()->SetNdivisions(410);
  hWZratio_RD->GetYaxis()->SetTitle("(#frac{1}{#sigma^{Z}} #frac{d#sigma^{Z}}{p_{T}^{Z}})/(#frac{1}{#sigma^{W}} #frac{d#sigma^{W}}{p_{T}^{W}})");
  hWZratio_RD->GetXaxis()->SetTitleOffset(0.6);
  hWZratio_RD->GetXaxis()->SetTitleSize(0.1);
  hWZratio_RD->GetXaxis()->SetLabelSize(0);
  hWZratio_RD->GetXaxis()->SetTitle("");
  hWZratio_RD->SetStats(0);
  hWZratio_RD->SetMarkerStyle(20);
  hWZratio_RD->SetMarkerColor(kBlack);
  hWZratio_RD->Draw("E1");

 // Draw TGraph style
  //tgWZratio_Powheg->SetFillStyle(3004);
  tgWZratio_Powheg->SetFillColor(kRed);
  tgWZratio_Powheg->SetLineColor(kRed+2);

  //tgWZratio_Resbos->SetFillStyle(3013);
  tgWZratio_Resbos->SetFillColor(kBlue);
  tgWZratio_Resbos->SetLineColor(kBlue+2);

  //tgWZratio_FEWZ->SetFillStyle(3005);
  tgWZratio_FEWZ->SetFillColor(kGreen);
  tgWZratio_FEWZ->SetLineColor(kGreen+2);
  
  tgWZratio_FEWZ->Draw("5");
  tgWZratio_Powheg->Draw("5");
  tgWZratio_Resbos->Draw("5");
  
  L1->Draw();
  cmspre->Draw();
  lumi->Draw();
  //Wchannel->Draw();
  //Zchannel->Draw();
  //LeptonCut->Draw();
  
  //C1->cd(2)->SetPad(0.028,0.505,0.76,0.805);
  //C1->cd(2)->SetPad(0.028,0.505,0.76,0.770);
  //C1->cd(2)->SetPad(0.028,0.505,0.76,0.760);
  C1->cd(2)->SetPad(0.078,0.504,0.722,0.730);
  
  C1->cd(2)->SetLogx(1);

  TH1D* hWZratio_RD_Zoom = (TH1D*)hWZratio_RD->Clone("hWZratio_RD_Zoom");
  
  hWZratio_RD_Zoom->GetYaxis()->SetTitle("");
  hWZratio_RD_Zoom->GetXaxis()->SetRangeUser(1,150);
  hWZratio_RD_Zoom->GetYaxis()->SetRangeUser(0.4,1.7);
  //hWZratio_RD_Zoom->GetYaxis()->SetLabelSize(0.05);
  hWZratio_RD_Zoom->GetYaxis()->SetLabelSize(0.065);
  hWZratio_RD_Zoom->GetYaxis()->SetNdivisions(405);
  hWZratio_RD_Zoom->GetXaxis()->SetTitle("");
  hWZratio_RD_Zoom->Draw(" E1");
  
  tgWZratio_Powheg->Draw("5");
  //tgWZratio_Resbos->Draw("5");
  tgWZratio_FEWZ->Draw("5");
  gPad->RedrawAxis();
  
  C1->cd(3)->SetPad(0.,0.1,0.95,0.5);
  C1->cd(3)->SetTopMargin(0.1);
  //C1->cd(3)->SetBottomMargin(0.01);
  C1->cd(3)->SetBottomMargin(0.1);
  C1->cd(3)->SetLeftMargin(0.15);
  C1->cd(3)->SetRightMargin(0.07);
  C1->cd(3)->SetTickx(1);
  C1->cd(3)->SetTicky(1);
  C1->cd(3)->SetLogx(1);

  TH1D *hRatioDummy = new TH1D("hRatioDummy","",nBins-1,WptBins);
  
  // set canvas range and XY axis title
  hRatioDummy->GetYaxis()->SetRangeUser(0.1,3);
 // hRatioDummy->GetYaxis()->SetTitle("Theory / Data");
  hRatioDummy->GetYaxis()->SetTitleSize(0.07);
  hRatioDummy->GetYaxis()->SetTitleOffset(0.58);
  //hRatioDummy->GetYaxis()->SetLabelSize(0.05);
  hRatioDummy->GetYaxis()->SetLabelSize(0.045);
  hRatioDummy->GetYaxis()->SetNdivisions(605);
  hRatioDummy->GetYaxis()->CenterTitle();
  
  hRatioDummy->GetXaxis()->SetTitleSize(0.07);
  hRatioDummy->GetXaxis()->SetTitleOffset(0.5);
  hRatioDummy->GetXaxis()->SetLabelSize(0);
  hRatioDummy->GetXaxis()->SetTitle("p_{T}^{V} [GeV]");

  // FEWZ Ratio plot setting
  tgRatioData->SetFillColor(kGray+2);
  tgRatioData->SetFillStyle(3354);

  tgRatioPowheg->SetFillColor(kRed);
  tgRatioPowheg->SetLineColor(kRed+2);
  tgRatioPowheg->SetMarkerStyle(22);
  tgRatioPowheg->SetMarkerColor(kRed+2);

  tgRatioResbos->SetFillColor(kBlue);
  tgRatioResbos->SetLineColor(kBlue+2);
  tgRatioResbos->SetMarkerStyle(20);
  tgRatioResbos->SetMarkerColor(kBlue+2);

  tgRatioFEWZ->SetFillColor(kGreen);
  tgRatioFEWZ->SetLineColor(kGreen+2);
  tgRatioFEWZ->SetMarkerStyle(21);
  tgRatioFEWZ->SetMarkerColor(kGreen+2);
 

   // Draw Ratio plot 
  hRatioDummy->Draw();
  tgRatioData->Draw("2");
  tgRatioFEWZ->Draw("5 P");  
  tgRatioPowheg->Draw("5 P");
  tgRatioResbos->Draw("5 P");
  gPad->RedrawAxis();
 
  // Set Zoomed Ratio plot
  //C1->cd(4)->SetPad(0.028,0.056,0.76,0.405);
  C1->cd(4)->SetPad(0.078,0.115,0.723,0.364);
  
  C1->cd(4)->SetLogx(1);
  
  TH1D *hRatioDummy_Zoom = (TH1D*)hRatioDummy->Clone("hRatioDummy_Zoom");
  hRatioDummy_Zoom->GetYaxis()->SetRangeUser(0.5,1.5); // Y axis range
  hRatioDummy_Zoom->GetXaxis()->SetRangeUser(1,150); // X axis range
  hRatioDummy_Zoom->GetYaxis()->SetTitle("Theory / Data");
  hRatioDummy_Zoom->GetYaxis()->SetTitleOffset(0.43);
  hRatioDummy_Zoom->GetYaxis()->SetTitleSize(0.11);
  hRatioDummy_Zoom->GetYaxis()->SetLabelSize(0.06);
  //hRatioDummy_Zoom->GetXaxis()->SetLabelSize(0.05);
  
  hRatioDummy_Zoom->GetXaxis()->SetLabelSize(0.08);
  hRatioDummy_Zoom->GetXaxis()->SetTitle("");
  
  hRatioDummy_Zoom->Draw();
  tgRatioData->Draw("2");
  tgRatioFEWZ->Draw("5 P");  
  tgRatioPowheg->Draw("5 P");
  tgRatioResbos->Draw("5 P");
  gPad->RedrawAxis();

  //C1->SaveAs("RatioNormZW_Fid_Zoom_logScale.png");
  C1->SaveAs("RatioNormZW_Fid_Zoom_logScale.pdf");

  return 0;
}

void FEWZ_PDFUncer(double *ZWRatio, double *Err)
{
  // FEWZ Normalized PDF error in % unit
  double PDFErr[12] = {0.,};

  PDFErr[0] = 2.51;  
  PDFErr[1] = 0.61;
  PDFErr[2] = 0.96;
  PDFErr[3] = 1.13;
  PDFErr[4] = 1.44;
  PDFErr[5] = 1.61;
  PDFErr[6] = 1.82;
  PDFErr[7] = 2.07;
  PDFErr[8] = 2.42;
  PDFErr[9] = 2.34;
  PDFErr[10] =2.30;
  PDFErr[11] =2.47;

  for(int i(0); i<12; i++)
  {
    Err[i] = ZWRatio[i] * PDFErr[i] * 0.01;
    //cout << "PDFErr in function : " << Err[i] << endl;
  }
}
void Powheg_PDFUncer(double *ZWRatio, double *Err)
{
  // Powheg Normalized PDF error in % unit
  double PDFErr[12] = {0.,};

  PDFErr[0] = 0.70904 ;  
  PDFErr[1] = 0.08829 ;
  PDFErr[2] = 0.26802 ;
  PDFErr[3] = 0.48621 ;
  PDFErr[4] = 0.74545 ;
  PDFErr[5] = 0.78758 ;
  PDFErr[6] = 0.89615 ;
  PDFErr[7] = 1.05564 ;
  PDFErr[8] = 1.23073 ;
  PDFErr[9] = 1.09974 ;
  PDFErr[10]= 1.36105 ;
  PDFErr[11]= 2.64723 ;

  for(int i(0); i<12; i++)
  {
    Err[i] = ZWRatio[i] * PDFErr[i] * 0.01;
    //cout << "PDFErr in function : " << Err[i] << endl;
  }
}
