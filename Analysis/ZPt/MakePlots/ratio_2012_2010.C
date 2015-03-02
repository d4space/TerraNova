#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "readStandardFile.C"

#include <algorithm>


int colorForId(int i) {
  const int colors[] = {kGreen+3, kRed, kAzure+5, kBlue,
			kBlack, kGray
  };
  return colors[i];
}

int fillStyleForId(int i) {
  const int fill[] = {1001,3001,3002
                      //const int fill[] = {1001,3011,3012
  };
  return fill[i];
}


double computeErrRatio(double A,
                       double errA,
                       double B,
                       double errB) ;

//------------------------------------------------------------------------------
// Print the title of a histogram
//------------------------------------------------------------------------------
void PrintIt(TPad *pad, TString title)//char *title)
{
  TLatex *latex = new TLatex();
  latex->SetTextFont(  42);
  latex->SetTextSize(0.05);
  
  // Get the most recent changes
  pad->Update();
  

  double xmin = pad->GetUxmin();
  double xmax = pad->GetUxmax();
  double ymin = pad->GetUymin();
  double ymax = pad->GetUymax();
 
  double xpos = xmin + 0.50*(xmax - xmin);
  double ypos = ymax + 0.05*(ymax - ymin);

  if (pad->GetLogy())  ypos = pow(10,ypos);
  if (pad->GetLogx())  xpos = pow(10,xpos);

  latex->SetTextAlign(22);
  latex->DrawLatex(xpos,ypos,title);
}

//
void lumiInfo(double px, double py, const char* lept="l",double sizef=1.0) {
  char eta[40], pt[40], pdfSet[40], factScale[40], renormScale[40];
  
  TLatex *plabel1 = new TLatex(px,py,"#int L = 36 pb^{-1} at #sqrt{s}=7 TeV");
  plabel1 -> SetNDC();
  plabel1 -> SetTextFont(42);
  plabel1 -> SetTextColor(kBlack);
  plabel1 -> SetTextSize(0.05*sizef);
  plabel1 -> SetTextAlign(22);
  plabel1 -> SetTextAngle(0);

  TLatex *plabel2 = new TLatex(px+0.01,py-0.1,"#int L = 18.4 pb^{-1} at #sqrt{s}=8 TeV");
  plabel2 -> SetNDC();
  plabel2 -> SetTextFont(42);
  plabel2 -> SetTextColor(kBlack);
  plabel2 -> SetTextSize(0.05*sizef);
  plabel2 -> SetTextAlign(22);

  //Then for each plot, pick a nice spot and draw
  plabel1 -> Draw();
  plabel2 -> Draw();

}


//void ratio_fewz2() {
void ratio_2012_2010() {
  
  //TString pathData   = "/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/results/data/";
  //TString pathTheory = "/data/uftrig01b/nwickram/root/CMSSW_5_2_4_patch3/results/thoery/";
  TString pathData   = "./data/";
  TString pathTheory = "./thoery/";
  
  std::cout << " reading file " << pathData + "data_2010.txt" << std::endl;
  TGraphAsymmErrors* gData7TeV = readStandardFilePtTGE(pathData + "data_2010.txt", 0);

  std::cout << " reading file " << pathData + "data_2012.txt" << std::endl;
  TGraphAsymmErrors* gData8TeV = readStandardFilePtTGE(pathData + "data_2012.txt", 0); 


  std::cout << " reading file " << pathTheory + "powheg_2010.txt" << std::endl;
  TGraphAsymmErrors* gPowheg7TeV = readStandardFilePtTGE(pathTheory + "powheg_2010.txt", 0);

  std::cout << " reading file " << pathTheory + "powheg_2012.txt" << std::endl;
  TGraphAsymmErrors* gPowheg8TeV = readStandardFilePtTGE(pathTheory + "powheg_2012.txt", 0); 
  

  std::cout << " reading file " << pathTheory + "Zpt_7TeV_fewz_CTEQ12nnlo.txt" << std::endl;
  TGraphAsymmErrors* gFEWZ_NNLO_7TeV = readStandardFilePtTGE(pathTheory + "Zpt_7TeV_fewz_CTEQ12nnlo.txt", 8);

  std::cout << " reading file " << pathTheory + "Zpt_8TeV_fewz_CTEQ12nnlo.txt" << std::endl;
  TGraphAsymmErrors* gFEWZ_NNLO_8TeV = readStandardFilePtTGE(pathTheory + "Zpt_8TeV_fewz_CTEQ12nnlo.txt", 8); 
  
  std::cout << " reading file " << pathTheory + "Zpt_8over7TeV_fewz_CTEQ12nnlo.txt" << std::endl;
  TGraphAsymmErrors* gFEWZ_NNLO_8over7TeV = readStandardFilePtTGE(pathTheory + "Zpt_8over7TeV_fewz2_CTEQ12nnlo.txt", 8); 
  
  TGraphAsymmErrors* ratioLowPt  = new TGraphAsymmErrors( 8);
  TGraphAsymmErrors* ratioHighPt = new TGraphAsymmErrors(10);

  TGraphAsymmErrors* ratioLowPt_Powheg= new TGraphAsymmErrors( 8);
  TGraphAsymmErrors* ratioHighPt_Fewz = new TGraphAsymmErrors(10);

  double px7tev,py7tev,erry7tev;
  double px8tev,py8tev,erry8tev;

  double px7tev_pow,py7tev_pow,erry7tev_pow;
  double px8tev_pow,py8tev_pow,erry8tev_pow;

  double px7tev_fewz,py7tev_fewz,erry7tev_fewz;
  double px8tev_fewz,py8tev_fewz,erry8tev_fewz;
  double erry8over7tev_fewz_l, erry8over7tev_fewz_h;

  for (int j=0; j<18; j++) {

    gData7TeV->GetPoint(j,px7tev,py7tev);
    if (py7tev==0) {
      std::cout << "Error: no x-sec value for 7 TeV at bin " << j << "\n";
      continue;
    }

    gData8TeV->GetPoint(j,px8tev,py8tev);
    if (py8tev==0)  {
      std::cout << "Error: no x-sec value for 8 TeV at bin " << j << "\n";
      continue;
    }
    if (px7tev != px8tev) {
      std::cout << "Error: different x values for bin " << j << "\n";
      continue;
    }

    erry8tev = gData8TeV->GetErrorYhigh(j);
    erry7tev = gData7TeV->GetErrorYhigh(j);
    printf("%d) %e | %e+/-%e - %e+/-%e\n",j, px7tev, py7tev,erry7tev, py8tev,erry8tev);
    double error = computeErrRatio(py8tev,erry8tev,py7tev,erry7tev);

    if ( j<8 ) {
      ratioLowPt -> SetPoint(j,px7tev,py8tev/py7tev);
      ratioLowPt -> SetPointError(j,0,0,error,error);
    }
    else {
      ratioHighPt -> SetPoint(j-8,px7tev,py8tev/py7tev);
      ratioHighPt -> SetPointError(j-8,0,0,error,error);
    }

    // powheg
    gPowheg7TeV->GetPoint(j,px7tev_pow,py7tev_pow);
    gPowheg8TeV->GetPoint(j,px8tev_pow,py8tev_pow);

    erry8tev_pow = gPowheg8TeV->GetErrorYhigh(j);
    erry7tev_pow = gPowheg7TeV->GetErrorYhigh(j);

    double error_pow = computeErrRatio(py8tev_pow,erry8tev_pow,py7tev_pow,erry7tev_pow);
    if ( j<8 ) {
      ratioLowPt_Powheg -> SetPoint(j,px7tev,py8tev_pow/py7tev_pow);
      ratioLowPt_Powheg -> SetPointError(j,0,0,error_pow,error_pow);
    }
  }

  
  for (int j=0; j<=9; j++) {
    // fewz
    //  gFEWZ_NNLO_7TeV->GetPoint(j,px7tev_fewz,py7tev_fewz);
    gFEWZ_NNLO_8over7TeV->GetPoint(j,px8tev_fewz,py8tev_fewz);

    erry8over7tev_fewz_h = gFEWZ_NNLO_8over7TeV -> GetErrorYhigh(j);
    erry8over7tev_fewz_l = gFEWZ_NNLO_8over7TeV -> GetErrorYlow(j);
    printf("%d) %e | %e | %e+/-%e\n",j,px8tev_fewz,py8tev_fewz, erry8over7tev_fewz_l,erry8over7tev_fewz_h );

    //   double error_fewz = computeErrRatio(py8tev_fewz,erry8tev_fewz,py7tev_fewz,erry7tev_fewz);
    //  ratioHighPt_Fewz -> SetPoint(j,px7tev_fewz,py8tev_fewz/py7tev_fewz);
    //ratioHighPt_Fewz -> SetPointError(j,0,0,error_fewz,error_fewz);
      ratioHighPt_Fewz -> SetPoint(j,px8tev_fewz,py8tev_fewz);
      ratioHighPt_Fewz -> SetPointError(j,0,0,erry8over7tev_fewz_l,erry8over7tev_fewz_h);
    
  }

    

  // low pT
  TCanvas* c1=new TCanvas("c1","c1",750,600);
  //c1->SetLogx();
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.20);

  ratioLowPt->SetTitle("");
  //ratioLowPt->GetXaxis()->SetTitle("q_{T} [GeV/c]");
  ratioLowPt->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  ratioLowPt->GetYaxis()->SetTitle("Ratio 8 TeV/7 TeV");
  ratioLowPt->GetYaxis()->CenterTitle();
  ratioLowPt->GetYaxis()->SetRangeUser(0,3.5);

  ratioLowPt->SetLineWidth(1);
  ratioLowPt->SetFillColor(kBlack);
  ratioLowPt->SetLineColor(kBlack);
  ratioLowPt->SetFillStyle(kBlack);
  ratioLowPt->SetMarkerStyle(20);
  ratioLowPt->Draw("AP");

  //ratioLowPt_Powheg->SetLineWidth(4);
  ratioLowPt_Powheg->SetFillColor(kGreen);
  ratioLowPt_Powheg->SetLineColor(kGreen);
  ratioLowPt_Powheg->SetFillStyle(3001);
  ratioLowPt_Powheg->Draw("3");

  ratioLowPt->Draw("P");

  TLine* lineLowPt = new TLine(0,1,20,1);
  lineLowPt -> SetLineStyle(2);
  lineLowPt -> SetLineWidth(2);
  lineLowPt -> Draw ("same");

  TLegend* tlLowPt = new TLegend(0.28, 0.64, 0.47, 0.83);
  tlLowPt->SetFillStyle (0);            
  tlLowPt->SetBorderSize(0);
  tlLowPt->SetTextFont(42);
  tlLowPt->SetTextSize(0.035);
  tlLowPt -> AddEntry(ratioLowPt, "data","lp");
  tlLowPt -> AddEntry(ratioLowPt_Powheg, "Powheg","f");
  tlLowPt -> Draw("same");

  lumiInfo(0.69,0.80,"",0.75);
  PrintIt(c1,"CMS Preliminary");
  //c1->SaveAs("diffXSec_muons_2012_vs_muons_2010_lowPt.png");


  // high pT
  TCanvas* c2=new TCanvas("c2","c2",750,600);
  c2->SetLogx();
  c2->SetRightMargin(0.02);
  c2->SetLeftMargin(0.20);

  ratioHighPt->SetTitle("");
  //ratioHighPt->GetXaxis()->SetTitle("q_{T} [GeV/c]");
  ratioHighPt->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  ratioHighPt->GetYaxis()->SetTitle("Ratio 8 TeV/7 TeV");
  ratioHighPt->GetYaxis()->CenterTitle();
  ratioHighPt->GetYaxis()->SetRangeUser(0,3.5);
  ratioHighPt->GetXaxis()->SetRangeUser(20,600);


  ratioHighPt->SetLineWidth(1);
  ratioHighPt->SetFillColor(kBlack);
  ratioHighPt->SetLineColor(kBlack);
  ratioHighPt->SetFillStyle(kBlack);
  ratioHighPt->SetMarkerStyle(20);
  ratioHighPt->Draw("AP");

  //ratioHighPt_Fewz->SetLineWidth(4);
  ratioHighPt_Fewz->SetFillColor(kRed);
  ratioHighPt_Fewz->SetLineColor(kRed);
  ratioHighPt_Fewz->SetFillStyle(3001);
  ratioHighPt_Fewz->Draw("3");

  ratioHighPt->Draw("P");


  TLine* lineHighPt = new TLine(25,1,460,1);
  lineHighPt -> SetLineStyle(2);
  lineHighPt -> SetLineWidth(2);
  lineHighPt -> Draw ("same");

  TLegend* tlHighPt = new TLegend(0.28, 0.71, 0.47, 0.90);
  tlHighPt->SetFillStyle (0);            
  tlHighPt->SetBorderSize(0);
  tlHighPt->SetTextFont(42);
  tlHighPt->SetTextSize(0.035);
  tlHighPt -> AddEntry(ratioHighPt, "data","lp");
  //  tlHighPt -> AddEntry(ratioHighPt_Fewz, "FEWZ NNLO+MSTW2008NNLO","f");
  tlHighPt -> AddEntry(ratioHighPt_Fewz, "FEWZ NNLO+CTEQ12NNLO","f");
  tlHighPt -> Draw("same");


  lumiInfo(0.42,0.62,"",0.75);
  PrintIt(c2,"CMS Preliminary");

  //c2->SaveAs("diffXSec_muons_2012_vs_muons_2010_highPt.png");
}

double computeErrRatio(double A,
                       double errA,
                       double B,
                       double errB) {

  double dRdA = 1/B;
  double dRdB = A/(B*B);
  
  double errorSquared = (dRdA*dRdA*errA*errA) + (dRdB*dRdB*errB*errB);
  
  return TMath::Sqrt( errorSquared );

}
